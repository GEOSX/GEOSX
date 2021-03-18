/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file Mesh.cpp
 */

#include "MultiscaleMeshLevel.hpp"

#include "common/TypeDispatch.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshData.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshUtils.hpp"
#include "linearAlgebra/multiscale/PartitionerBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/utilities/MeshMapUtilities.hpp"

namespace geosx
{

using namespace dataRepository;

namespace multiscale
{

MultiscaleMeshLevel::MultiscaleMeshLevel( string name )
  : m_name( std::move( name ) ),
  m_root( m_name, m_rootNode ),
  m_cellManager( m_name + "_cellManager", &m_root ),
  m_nodeManager( m_name + "_nodeManager", &m_root )
{
  m_cellManager.registerExtrinsicData< meshData::OrigElementRegion >( m_name );
  m_cellManager.registerExtrinsicData< meshData::OrigElementSubRegion >( m_name );
  m_cellManager.registerExtrinsicData< meshData::OrigElementIndex >( m_name );

  m_nodeManager.registerExtrinsicData< meshData::OrigNodeIndex >( m_name );
}

namespace fineGridConstruction
{

localIndex assignNodeLocalIndices( NodeManager & nodeManager,
                                   string const & localIndexKey,
                                   string const & globalIndexKey,
                                   globalIndex const rankOffset,
                                   localIndex const numLocalNodes )
{
  arrayView1d< globalIndex const > const nodeGlobalIndex = nodeManager.getReference< array1d< globalIndex > >( globalIndexKey ).toViewConst();
  arrayView1d< localIndex > const nodeLocalIndex = nodeManager.getReference< array1d< localIndex > >( localIndexKey ).toView();
  localIndex numPresentNodes = numLocalNodes;

  // Assign new local indices s.t. locally owned objects come first, followed by ghosted ones.
  // We do this by looping and filtering over all nodes of the original mesh, as this is the primary use case.
  // It might pessimize the case of a small flow region in a larger mechanics mesh, but it's less critical.
  // Note: there doesn't seem to be a good way to do this in parallel ... or does it?
  forAll< serialPolicy >( nodeManager.size(), [&]( localIndex const k )
  {
    globalIndex const newGlobalIndex = nodeGlobalIndex[k];
    if( newGlobalIndex >= 0 )
    {
      localIndex const localIndexTmp = static_cast< localIndex >( newGlobalIndex - rankOffset );
      localIndex const newLocalIndex = ( 0 <= localIndexTmp && localIndexTmp < numLocalNodes ) ? localIndexTmp : numPresentNodes++;
      nodeLocalIndex[k] = newLocalIndex;
    }
  } );

  return numPresentNodes;
}

void copyNodeData( NodeManager const & nodeManager,
                   string const & localIndexKey,
                   string const & globalIndexKey,
                   MultiscaleMeshObjectManager & msNodeManager )
{
  arrayView1d< localIndex const > const nodeLocalIndex = nodeManager.getReference< array1d< localIndex > >( localIndexKey );
  arrayView1d< globalIndex const > const nodeGlobalIndex = nodeManager.getReference< array1d< globalIndex > >( globalIndexKey );

  arrayView1d< integer const > const origGhostRank  = nodeManager.ghostRank().toViewConst();
  arrayView1d< integer const > const origIsExternal = nodeManager.isExternal().toViewConst();
  arrayView1d< integer const > const origBoundary   = nodeManager.getDomainBoundaryIndicator().toViewConst();

  arrayView1d< integer >     const msGhostRank     = msNodeManager.ghostRank().toView();
  arrayView1d< integer >     const msIsExternal    = msNodeManager.isExternal().toView();
  arrayView1d< integer >     const msBoundary      = msNodeManager.getDomainBoundaryIndicator().toView();
  arrayView1d< globalIndex > const msLocalToGlobal = msNodeManager.localToGlobalMap().toView();
  arrayView1d< localIndex >  const msOrigIndex     = msNodeManager.getExtrinsicData< meshData::OrigNodeIndex >().toView();

  forAll< parallelHostPolicy >( nodeManager.size(), [=]( localIndex const k )
  {
    localIndex const nodeIdx = nodeLocalIndex[k];
    if( nodeIdx >= 0 )
    {
      msLocalToGlobal[nodeIdx] = nodeGlobalIndex[k];
      msGhostRank[nodeIdx] = origGhostRank[k];
      msIsExternal[nodeIdx] = origIsExternal[k];
      msBoundary[nodeIdx] = origBoundary[k];
      msOrigIndex[nodeIdx] = k;
    }
  } );
  msNodeManager.constructGlobalToLocalMap();
  msNodeManager.setMaxGlobalIndex();
}

void populateNodeManager( string const & localIndexKey,
                          std::vector< int > const & neighborRanks,
                          MeshLevel & mesh,
                          std::vector< string > const & regions,
                          MultiscaleMeshObjectManager & msNodeManager )
{
  // Use a temporary DofManager to generate new contiguous global index numbering
  string const dofName = "globalIndex";
  DofManager dofManager( "nodeIndexMgr" );
  dofManager.setMesh( mesh );
  dofManager.addField( dofName, DofManager::Location::Node, 1, regions );
  dofManager.reorderByRank();

  string const globalIndexKey = dofManager.getKey( dofName );
  NodeManager & nodeManager = mesh.getNodeManager();

  localIndex const numLocalNodes = dofManager.numLocalDofs( dofName );
  localIndex const numPresentNodes = assignNodeLocalIndices( nodeManager,
                                                             localIndexKey,
                                                             globalIndexKey,
                                                             dofManager.rankOffset(),
                                                             numLocalNodes );

  // Setup neighbor data
  for( int const rank : neighborRanks )
  {
    msNodeManager.addNeighbor( rank );
  }

  // Populate the new node manager, now that we know its total size
  msNodeManager.resize( numPresentNodes );
  msNodeManager.setNumOwnedObjects( numLocalNodes );
  copyNodeData( nodeManager, localIndexKey, globalIndexKey, msNodeManager );
  meshUtils::copySets( nodeManager, localIndexKey, msNodeManager );
  meshUtils::copyNeighborData( nodeManager, localIndexKey, neighborRanks, msNodeManager, meshUtils::filterArray< localIndex > );

  // Remove index fields on the GEOSX mesh
  dofManager.clear();
}

localIndex assignCellLocalIndices( ElementRegionManager & elemManager,
                                   std::vector< string > const & regions,
                                   string const & localIndexKey,
                                   string const & globalIndexKey,
                                   globalIndex const rankOffset,
                                   localIndex const numLocalCells )
{
  localIndex numPresentCells = numLocalCells;

  // Assign new local indices s.t. locally owned objects come first, followed by ghosted ones.
  elemManager.forElementSubRegions( regions, [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    // Allocate a element-based field on the GEOSX mesh that will keep a mapping to new cell local indices
    arrayView1d< localIndex > const cellLocalIndex = subRegion.getReference< array1d< localIndex > >( localIndexKey ).toView();
    arrayView1d< globalIndex const > const cellGlobalIndex = subRegion.getReference< array1d< globalIndex > >( globalIndexKey ).toViewConst();

    forAll< serialPolicy >( subRegion.size(), [&]( localIndex const ei )
    {
      localIndex const localIndexTmp = static_cast< localIndex >( cellGlobalIndex[ei] - rankOffset );
      cellLocalIndex[ei] = ( 0 <= localIndexTmp && localIndexTmp < numLocalCells ) ? localIndexTmp : numPresentCells++;
    } );
  } );

  return numPresentCells;
}

void copyCellData( ElementRegionManager const & elemManager,
                   std::vector< string > const & regions,
                   string const & localIndexKey,
                   string const & globalIndexKey,
                   MultiscaleMeshObjectManager & msCellManager )
{
  arrayView1d< integer >     const msGhostRank        = msCellManager.ghostRank().toView();
  arrayView1d< integer >     const msIsExternal       = msCellManager.isExternal().toView();
  arrayView1d< integer >     const msBoundary         = msCellManager.getDomainBoundaryIndicator().toView();
  arrayView1d< globalIndex > const msLocalToGlobal    = msCellManager.localToGlobalMap().toView();
  arrayView1d< localIndex >  const msOrigRegIndex     = msCellManager.getExtrinsicData< meshData::OrigElementRegion >().toView();
  arrayView1d< localIndex >  const msOrigSubRegIndex  = msCellManager.getExtrinsicData< meshData::OrigElementSubRegion >().toView();
  arrayView1d< localIndex >  const msOrigElemIndex    = msCellManager.getExtrinsicData< meshData::OrigElementIndex >().toView();

  elemManager.forElementSubRegionsComplete( regions, [&]( localIndex,
                                                          localIndex const er,
                                                          localIndex const esr,
                                                          ElementRegionBase const &,
                                                          ElementSubRegionBase const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex = subRegion.getReference< array1d< localIndex > >( localIndexKey ).toViewConst();
    arrayView1d< globalIndex const > const cellGlobalIndex = subRegion.getReference< array1d< globalIndex > >( globalIndexKey ).toViewConst();

    arrayView1d< integer const > const origGhostRank  = subRegion.ghostRank().toViewConst();
    arrayView1d< integer const > const origIsExternal = subRegion.isExternal().toViewConst();
    arrayView1d< integer const > const origBoundary   = subRegion.getDomainBoundaryIndicator().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei )
    {
      localIndex const cellIdx = cellLocalIndex[ei];
      msLocalToGlobal[cellIdx] = cellGlobalIndex[ei];
      msGhostRank[cellIdx] = origGhostRank[ei];
      msIsExternal[cellIdx] = origIsExternal[ei];
      msBoundary[cellIdx] = origBoundary[ei];
      msOrigRegIndex[cellIdx] = er;
      msOrigSubRegIndex[cellIdx] = esr;
      msOrigElemIndex[cellIdx] = ei;
    } );
  } );
  msCellManager.constructGlobalToLocalMap();
  msCellManager.setMaxGlobalIndex();
}

void populateCellManager( string const & localIndexKey,
                          std::vector< int > const & neighborRanks,
                          MeshLevel & mesh,
                          std::vector< string > const & regions,
                          MultiscaleMeshObjectManager & msCellManager )
{
  // Use a temporary DofManager to generate new contiguous global index numbering
  string const dofName = "globalIndex";
  DofManager dofManager( "cellIndexMgr" );
  dofManager.setMesh( mesh );
  dofManager.addField( dofName, DofManager::Location::Elem, 1, regions );
  dofManager.reorderByRank();

  string const globalIndexKey = dofManager.getKey( dofName );
  ElementRegionManager & elemManager = mesh.getElemManager();

  localIndex const numLocalCells = dofManager.numLocalDofs( dofName );
  localIndex const numPresentCells = assignCellLocalIndices( elemManager,
                                                             regions,
                                                             localIndexKey,
                                                             globalIndexKey,
                                                             dofManager.rankOffset(),
                                                             numLocalCells );

  // Setup neighbor data
  for( int const rank : neighborRanks )
  {
    msCellManager.addNeighbor( rank );
  }

  // Populate the new cell manager, now that we know its total size
  msCellManager.resize( numPresentCells );
  msCellManager.setNumOwnedObjects( numLocalCells );
  copyCellData( elemManager, regions, localIndexKey, globalIndexKey, msCellManager );
  elemManager.forElementSubRegions( regions, [&]( localIndex, ElementSubRegionBase const & subRegion )
  {
    meshUtils::copySets( subRegion, localIndexKey, msCellManager );
    meshUtils::copyNeighborData( subRegion, localIndexKey, neighborRanks, msCellManager, meshUtils::filterArray< localIndex > );
  } );

  // Make sure to remove temporary global index field on the GEOSX mesh
  dofManager.clear();
}

void buildNodeToCellMap( string const & cellLocalIndexKey,
                         MeshLevel const & mesh,
                         MultiscaleMeshObjectManager & msNodeManager )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const origNodeIndex = msNodeManager.getExtrinsicData< meshData::OrigNodeIndex >().toViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex const > > cellLocalIndex =
    elemManager.constructArrayViewAccessor< localIndex, 1 >( cellLocalIndexKey );

  ArrayOfArraysView< localIndex const > const elemRegion = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const elemSubRegion = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const elemIndex = nodeManager.elementList().toViewConst();

  // Count the length of each sub-array
  array1d< localIndex > rowCounts( msNodeManager.size() );
  forAll< parallelHostPolicy >( msNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const k )
  {
    rowCounts[k] = elemIndex.sizeOfArray( origNodeIndex[k] );
  } );

  // Resize
  MultiscaleMeshObjectManager::MapType & nodeToCell = msNodeManager.toDualRelation();
  nodeToCell.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  forAll< parallelHostPolicy >( msNodeManager.size(), [&]( localIndex const k )
  {
    localIndex const nodeIdx = origNodeIndex[k];
    arraySlice1d< localIndex const > const er = elemRegion[nodeIdx];
    arraySlice1d< localIndex const > const es = elemSubRegion[nodeIdx];
    arraySlice1d< localIndex const > const ei = elemIndex[nodeIdx];
    for( localIndex a = 0; a < ei.size(); ++a )
    {
      nodeToCell.insertIntoSet( k, cellLocalIndex[er[a]][es[a]][ei[a]] );
    }
  } );
}

void buildCellToNodeMap( string const & cellLocalIndexKey,
                         string const & nodeLocalIndexKey,
                         MeshLevel const & mesh,
                         std::vector< string > const & regions,
                         MultiscaleMeshObjectManager & msCellManager )
{
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const nodeLocalIndex =
    nodeManager.getReference< array1d< localIndex > >( nodeLocalIndexKey ).toViewConst();

  // Count the length of each row
  array1d< localIndex > rowCounts( msCellManager.size() );
  elemManager.forElementSubRegions( regions, [&]( localIndex, auto const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex =
      subRegion.template getReference< array1d< localIndex > >( cellLocalIndexKey ).toViewConst();
    auto const nodes = subRegion.nodeList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [=, rowCounts = rowCounts.toView()]( localIndex const ei )
    {
      rowCounts[cellLocalIndex[ei]] = meshMapUtilities::size1( nodes, ei );
    } );
  } );

  // Resize the map
  MultiscaleMeshObjectManager::MapType & cellToNode = msCellManager.toDualRelation();
  cellToNode.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  elemManager.forElementSubRegions( regions, [&]( localIndex, auto const & subRegion )
  {
    arrayView1d< localIndex const > const cellLocalIndex =
      subRegion.template getReference< array1d< localIndex > >( cellLocalIndexKey ).toViewConst();
    auto const nodes = subRegion.nodeList().toViewConst();

    forAll< parallelHostPolicy >( subRegion.size(), [=, cellToNode = cellToNode.toView()]( localIndex const ei )
    {
      localIndex const cellIdx = cellLocalIndex[ei];
      for( localIndex a = 0; a < meshMapUtilities::size1( nodes, ei ); ++a )
      {
        localIndex const k = meshMapUtilities::value( nodes, ei, a );
        cellToNode.insertIntoSet( cellIdx, nodeLocalIndex[k] );
      }
    } );
  } );
}

} // namespace fineGridConstruction

void MultiscaleMeshLevel::buildFineMesh( MeshLevel & mesh,
                                         std::vector< string > const & regions )
{
  GEOSX_MARK_FUNCTION;

  // TODO: get rid of getParent() use
  m_sourceMesh = &mesh;
  m_regions = regions;
  m_domain = dynamicCast< DomainPartition * >( &mesh.getParent().getParent().getParent() );
  std::vector< int > neighborRanks = m_domain->getNeighborRanks();

  string const cellLocalIndexKey = m_name + "_cell_localIndex";
  string const nodeLocalIndexKey = m_name + "_node_localIndex";

  // Allocate fields on the GEOSX mesh that will keep mappings to new local indices
  mesh.getNodeManager().registerWrapper< array1d< localIndex > >( nodeLocalIndexKey ).setApplyDefaultValue( -1 );
  mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    subRegion.registerWrapper< array1d< localIndex > >( cellLocalIndexKey ).setApplyDefaultValue( -1 );
  } );

  // Generate new contiguous local/global numberings for the target subset of nodes/cells
  fineGridConstruction::populateNodeManager( nodeLocalIndexKey, neighborRanks, mesh, regions, m_nodeManager );
  fineGridConstruction::populateCellManager( cellLocalIndexKey, neighborRanks, mesh, regions, m_cellManager );

  // Extract appropriately renumbered cell-node maps into new managers
  fineGridConstruction::buildNodeToCellMap( cellLocalIndexKey, mesh, m_nodeManager );
  fineGridConstruction::buildCellToNodeMap( cellLocalIndexKey, nodeLocalIndexKey, mesh, regions, m_cellManager );

  // Remove local index fields used during construction from GEOSX mesh
  mesh.getNodeManager().deregisterWrapper( nodeLocalIndexKey );
  mesh.getElemManager().forElementSubRegions( regions, [&]( localIndex, ElementSubRegionBase & subRegion )
  {
    subRegion.deregisterWrapper( cellLocalIndexKey );
  } );
}

namespace coarseGridConstruction
{

void buildFineCellLists( MultiscaleMeshObjectManager const & fineCellManager,
                         MultiscaleMeshObjectManager & coarseCellManager )
{
  ArrayOfArrays< localIndex > & fineCellLists = coarseCellManager.getExtrinsicData< meshData::FineCellLocalIndices >();
  arrayView1d< localIndex const > const coarseCellLocalIndex = fineCellManager.getExtrinsicData< meshData::CoarseCellLocalIndex >();

  // Calculate the size of each list
  array1d< localIndex > fineCellListSizes( coarseCellManager.size() );
  forAll< parallelHostPolicy >( fineCellManager.size(), [=, sizes = fineCellListSizes.toView()] ( localIndex const ic )
  {
    RAJA::atomicInc( parallelHostAtomic{}, &sizes[coarseCellLocalIndex[ic]] );
  } );

  // Allocate space for each list
  fineCellLists.resizeFromCapacities< parallelHostPolicy >( fineCellListSizes.size(), fineCellListSizes.data() );

  // Populate the lists
  forAll< parallelHostPolicy >( fineCellManager.size(), [=, lists = fineCellLists.toView()] ( localIndex const ic )
  {
    lists.emplaceBackAtomic< parallelHostAtomic >( coarseCellLocalIndex[ic], ic );
  } );

  // Sort the lists for potentially better performance when using to sub-index larger arrays
  forAll< parallelHostPolicy >( coarseCellManager.size(), [lists = fineCellLists.toView()]( localIndex const icc )
  {
    arraySlice1d< localIndex > const list = lists[icc];
    LvArray::sortedArrayManipulation::makeSorted( list.begin(), list.end() );
  } );
}

void fillBasicCellData( MultiscaleMeshObjectManager const & fineCellManager,
                        MultiscaleMeshObjectManager & coarseCellManager )
{
  arrayView1d< localIndex const > const coarseCellIndex = fineCellManager.getExtrinsicData< meshData::CoarseCellLocalIndex >().toViewConst();
  arrayView1d< integer const > const fineGhostRank = fineCellManager.ghostRank();
  arrayView1d< integer const > const fineIsExternal = fineCellManager.isExternal();
  arrayView1d< integer const > const fineDomainBoundary = fineCellManager.getDomainBoundaryIndicator();

  arrayView1d< integer > const coarseGhostRank = coarseCellManager.ghostRank();
  arrayView1d< integer > const coarseIsExternal = coarseCellManager.isExternal();
  arrayView1d< integer > const coarseDomainBoundary = coarseCellManager.getDomainBoundaryIndicator();

  forAll< parallelHostPolicy >( fineCellManager.size(), [=]( localIndex const icf )
  {
    localIndex const icc = coarseCellIndex[icf];
    RAJA::atomicMax( parallelHostAtomic{}, &coarseGhostRank[icc], fineGhostRank[icf] );
    RAJA::atomicMax( parallelHostAtomic{}, &coarseIsExternal[icc], fineIsExternal[icf] );
    RAJA::atomicMax( parallelHostAtomic{}, &coarseDomainBoundary[icc], fineDomainBoundary[icf] );
  } );
}

void buildCellLocalToGlobalMaps( std::set< globalIndex > const & ghostGlobalIndices,
                                 globalIndex const rankOffset,
                                 MultiscaleMeshObjectManager & coarseCellManager )
{
  arrayView1d< globalIndex > const coarseLocalToGlobal = coarseCellManager.localToGlobalMap();
  {
    localIndex icc = 0;
    for(; icc < coarseCellManager.numOwnedObjects(); ++icc )
    {
      coarseLocalToGlobal[icc] = rankOffset + icc;
    }
    for( globalIndex const coarseGlobalIndex : ghostGlobalIndices )
    {
      coarseLocalToGlobal[icc++] = coarseGlobalIndex;
    }
  }
  coarseCellManager.constructGlobalToLocalMap();
  coarseCellManager.setMaxGlobalIndex();
}

template< typename T >
struct SetCompare
{
  ArrayOfSetsView< T const > const & sets;
  bool operator()( localIndex const i, localIndex const j ) const
  {
    arraySlice1d< T const > const si = sets[i];
    arraySlice1d< T const > const sj = sets[j];
    return std::lexicographical_compare( si.begin(), si.end(), sj.begin(), sj.end() );
  }
};

array1d< localIndex > findCoarseNodes( MultiscaleMeshObjectManager const & fineNodeManager,
                                       MultiscaleMeshObjectManager const & fineCellManager,
                                       ArrayOfSetsView< globalIndex const > const & nodeToSubdomain )
{
  // Construct a list of "skeleton" nodes (those with 3 or more adjacent subdomains)
  array1d< localIndex > skelNodes;
  for( localIndex inf = 0; inf < fineNodeManager.size(); ++inf )
  {
    if( nodeToSubdomain.sizeOfSet( inf ) >= 3 )
    {
      skelNodes.emplace_back( inf );
    }
  }

  // Sort skeleton nodes according to subdomain lists so as to locate nodes of identical adjacencies
  SetCompare< globalIndex > const adjacencyComp{ nodeToSubdomain.toViewConst() };
  std::sort( skelNodes.begin(), skelNodes.end(), adjacencyComp );

  // Identify "features" (groups of skeleton nodes with the same subdomain adjacency)
  array1d< localIndex > const featureIndex( fineNodeManager.size() );
  featureIndex.setValues< serialPolicy >( -1 );

  ArrayOfArrays< localIndex > featureNodes;
  featureNodes.reserve( skelNodes.size() ); // overallocate to avoid reallocation
  featureNodes.reserveValues( skelNodes.size() ); // precise allocation

  localIndex numFeatures = 0;
  featureNodes.appendArray( 0 );
  featureNodes.emplaceBack( numFeatures, skelNodes[0] );
  for( localIndex i = 1; i < skelNodes.size(); ++i )
  {
    if( adjacencyComp( skelNodes[i-1], skelNodes[i] ) )
    {
      ++numFeatures;
      featureNodes.appendArray( 0 );
    }
    featureNodes.emplaceBack( numFeatures, skelNodes[i] );
    featureIndex[skelNodes[i]] = numFeatures;
  }
  ++numFeatures;

  // Construct feature-to-feature adjacency
  ArrayOfSets< localIndex > const featureAdjacency( numFeatures, 64 );
  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = fineNodeManager.toDualRelation().toViewConst();
  MultiscaleMeshObjectManager::MapViewConst const cellToNode = fineCellManager.toDualRelation().toViewConst();

  forAll< parallelHostPolicy >( numFeatures, [nodeToCell, cellToNode,
                                              featureNodes = featureNodes.toViewConst(),
                                              featureIndex = featureIndex.toViewConst(),
                                              featureAdjacency = featureAdjacency.toView()]( localIndex const f )
  {
    for( localIndex const inf : featureNodes[f] )
    {
      meshUtils::forUniqueNeighbors< 128 >( inf, nodeToCell, cellToNode, [&]( localIndex const nbrIdx )
      {
        if( nbrIdx != inf && featureIndex[nbrIdx] >= 0 )
        {
          featureAdjacency.insertIntoSet( f, featureIndex[nbrIdx] );
        }
      } );
    }
  } );

  // Choose features that represent coarse nodes (highest adjacency among neighbors)
  array1d< integer > const isCoarseNode( numFeatures );
  forAll< parallelHostPolicy >( numFeatures, [isCoarseNode = isCoarseNode.toView(),
                                              featureNodes = featureNodes.toViewConst(),
                                              featureAdjacency = featureAdjacency.toViewConst(),
                                              nodeToSubdomain = nodeToSubdomain.toViewConst()]( localIndex const f )
  {
    arraySlice1d< globalIndex const > const subs = nodeToSubdomain[ featureNodes( f, 0 ) ];
    for( localIndex const f_nbr : featureAdjacency[f] )
    {
      arraySlice1d< globalIndex const > const subs_nbr = nodeToSubdomain[ featureNodes( f_nbr, 0 ) ];
      if( std::includes( subs_nbr.begin(), subs_nbr.end(), subs.begin(), subs.end() ) )
      {
        // discard feature if its subdomain adjacency is fully included in any of its direct neighbors
        return;
      }
    }
    // if not discarded, it is a coarse node
    isCoarseNode[f] = 1;
  } );

  // Make a list of fine-scale indices of coarse nodes that are locally owned
  array1d< localIndex > coarseNodes;
  for( localIndex f = 0; f < numFeatures; ++f )
  {
    if( isCoarseNode[f] == 1 )
    {
      arraySlice1d< localIndex const > const nodes = featureNodes[f];
      for( localIndex inf : nodes )
      {
        coarseNodes.emplace_back( inf );
      }
    }
  }

  return coarseNodes;
}

void fillBasicNodeData( MultiscaleMeshObjectManager const & fineNodeManager,
                        MultiscaleMeshObjectManager & coarseNodeManager )
{
  arrayView1d< localIndex const > const fineNodeIndex = coarseNodeManager.getExtrinsicData< meshData::FineNodeLocalIndex >().toViewConst();
  arrayView1d< integer const > const fineGhostRank = fineNodeManager.ghostRank();
  arrayView1d< integer const > const fineIsExternal = fineNodeManager.isExternal();
  arrayView1d< integer const > const fineDomainBoundary = fineNodeManager.getDomainBoundaryIndicator();

  arrayView1d< integer > const coarseGhostRank = coarseNodeManager.ghostRank();
  arrayView1d< integer > const coarseIsExternal = coarseNodeManager.isExternal();
  arrayView1d< integer > const coarseDomainBoundary = coarseNodeManager.getDomainBoundaryIndicator();

  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=]( localIndex const inc )
  {
    localIndex const inf = fineNodeIndex[inc];
    coarseGhostRank[inc] = fineGhostRank[inf];
    coarseIsExternal[inc] = fineIsExternal[inf];
    coarseDomainBoundary[inc] = fineDomainBoundary[inf];
  } );
}

ArrayOfSets< globalIndex >
buildFineNodeToGlobalSubdomainMap( MultiscaleMeshObjectManager const & fineNodeManager,
                                   MultiscaleMeshObjectManager const & fineCellManager,
                                   MultiscaleMeshObjectManager const & coarseCellManager,
                                   arrayView1d< string const > const & boundaryNodeSets )
{
  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = fineNodeManager.toDualRelation().toViewConst();

  // count the row lengths
  array1d< localIndex > rowCounts( fineNodeManager.size() );
  forAll< parallelHostPolicy >( fineNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inf )
  {
    rowCounts[inf] = nodeToCell.sizeOfSet( inf );
  } );
  for( string const & setName : boundaryNodeSets )
  {
    SortedArrayView< localIndex const > const set = fineNodeManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, rowCounts = rowCounts.toView()]( localIndex const i )
    {
      ++rowCounts[set[i]];
    } );
  }

  // Resize from row lengths
  ArrayOfSets< globalIndex > nodeToSubdomain;
  nodeToSubdomain.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  arrayView1d< globalIndex const > const coarseCellGlobalIndex = fineCellManager.getExtrinsicData< meshData::CoarseCellGlobalIndex >();
  forAll< parallelHostPolicy >( fineNodeManager.size(), [=, nodeToSub = nodeToSubdomain.toView()]( localIndex const inf )
  {
    for( localIndex const icf : nodeToCell[inf] )
    {
      nodeToSub.insertIntoSet( inf, coarseCellGlobalIndex[icf] );
    }
  } );
  globalIndex numSubdomains = coarseCellManager.maxGlobalIndex() + 1;
  for( string const & setName : boundaryNodeSets )
  {
    SortedArrayView< localIndex const > const set = fineNodeManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, nodeToSub = nodeToSubdomain.toView()]( localIndex const inf )
    {
      nodeToSub.insertIntoSet( set[inf], numSubdomains );
    } );
    ++numSubdomains;
  }

  return nodeToSubdomain;
}

void buildNodeToCellMap( MultiscaleMeshObjectManager const & fineCellManager,
                         MultiscaleMeshObjectManager const & fineNodeManager,
                         MultiscaleMeshObjectManager & coarseNodeManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< localIndex const > const coarseCellLocalIndex = fineCellManager.getExtrinsicData< meshData::CoarseCellLocalIndex >();
  arrayView1d< localIndex const > const fineNodeLocalIndex = coarseNodeManager.getExtrinsicData< meshData::FineNodeLocalIndex >();
  MultiscaleMeshObjectManager::MapViewConst const fineNodeToCell = fineNodeManager.toDualRelation().toViewConst();

  // First pass: count length of each sub-array in order to do exact allocation
  array1d< localIndex > rowCounts( coarseNodeManager.size() );
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inc )
  {
    localIndex coarseCellIndices[32];
    localIndex count = 0;
    for( localIndex const icf : fineNodeToCell[fineNodeLocalIndex[inc]] )
    {
      coarseCellIndices[count++] = coarseCellLocalIndex[icf];
    }
    rowCounts[inc] = LvArray::sortedArrayManipulation::makeSortedUnique( coarseCellIndices, coarseCellIndices + count );
  } );

  // Resize the map
  MultiscaleMeshObjectManager::MapType & coarseNodeToCell = coarseNodeManager.toDualRelation();
  coarseNodeToCell.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Second pass: fill the map
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, coarseNodeToCell = coarseNodeToCell.toView()]( localIndex const inc )
  {
    for( localIndex const icf : fineNodeToCell[fineNodeLocalIndex[inc]] )
    {
      coarseNodeToCell.insertIntoSet( inc, coarseCellLocalIndex[icf] );
    }
  } );
}

void buildCellToNodeMap( MultiscaleMeshObjectManager const & coarseNodeManager,
                         MultiscaleMeshObjectManager & coarseCellManager )
{
  GEOSX_MARK_FUNCTION;

  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = coarseNodeManager.toDualRelation().toViewConst();
  MultiscaleMeshObjectManager::MapType & cellToNode = coarseCellManager.toDualRelation();

  // First pass: count the row lengths in transpose map
  array1d< localIndex > rowCounts( cellToNode.size() );
  forAll< parallelHostPolicy >( nodeToCell.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inc )
  {
    for( localIndex const icc : nodeToCell[inc] )
    {
      RAJA::atomicInc< parallelHostAtomic >( &rowCounts[icc] );
    }
  } );

  // Create and resize the temporary map
  ArrayOfArrays< localIndex > cellToNodeTemp;
  cellToNodeTemp.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Second pass: fill the map
  forAll< parallelHostPolicy >( nodeToCell.size(), [=, cellToNode = cellToNodeTemp.toView()]( localIndex const inc )
  {
    for( localIndex const icc : nodeToCell[inc] )
    {
      cellToNode.emplaceBackAtomic< parallelHostAtomic >( icc, inc );
    }
  } );

  // Move the temp map and sort the entries
  cellToNode.assimilate( std::move( cellToNodeTemp ), LvArray::sortedArrayManipulation::UNSORTED_NO_DUPLICATES );
}

} // namespace coarseGridConstruction

void MultiscaleMeshLevel::buildCoarseCells( MultiscaleMeshLevel & fineMesh,
                                            LinearSolverParameters::Multiscale::Coarsening const & coarse_params )
{
  GEOSX_MARK_FUNCTION;

  MultiscaleMeshObjectManager & fineCellManager = fineMesh.cellManager();

  // Allocate arrays to hold partitioning info (local and global)
  arrayView1d< localIndex > const & coarseLocalIndex =
    fineCellManager.registerExtrinsicData< meshData::CoarseCellLocalIndex >( m_name ).referenceAsView();
  arrayView1d< globalIndex > const & coarseGlobalIndex =
    fineCellManager.registerExtrinsicData< meshData::CoarseCellGlobalIndex >( m_name ).referenceAsView();

  // Generate the partitioning locally
  std::unique_ptr< PartitionerBase > partitioner = PartitionerBase::create( coarse_params );
  localIndex const numLocalCoarseCells = partitioner->generate( fineMesh, coarseLocalIndex );
  GEOSX_LOG_RANK( "Multiscale " << m_name << ": generated " << numLocalCoarseCells << " local coarse cells" );

  // Compute global number of partitions
  globalIndex const rankOffset = MpiWrapper::prefixSum< globalIndex >( numLocalCoarseCells );

  // Fill in partition global index for locally owned cells
  for( localIndex icf = 0; icf < fineCellManager.numOwnedObjects(); ++icf )
  {
    coarseGlobalIndex[icf] = coarseLocalIndex[icf] + rankOffset;
  }

  // Synchronize partition global index across ranks
  string_array fieldNames;
  fieldNames.emplace_back( meshData::CoarseCellGlobalIndex::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames, fineCellManager, m_domain->getNeighbors(), false );

  // Scan ghosted cells and collect all new partition global indices
  std::set< globalIndex > ghostGlobalIndices;
  for( localIndex icf = fineCellManager.numOwnedObjects(); icf < fineCellManager.size(); ++icf )
  {
    ghostGlobalIndices.insert( coarseGlobalIndex[icf] );
  }
  localIndex const numPresentCoarseCells = numLocalCoarseCells + LvArray::integerConversion< localIndex >( ghostGlobalIndices.size() );

  // Resize and start populating coarse cell manager
  m_cellManager.resize( numPresentCoarseCells );
  m_cellManager.setNumOwnedObjects( numLocalCoarseCells );

  // Populate coarse local-global maps
  coarseGridConstruction::buildCellLocalToGlobalMaps( ghostGlobalIndices, rankOffset, m_cellManager );
  GEOSX_LOG_RANK_0( "Multiscale " << m_name << ": generated " << m_cellManager.maxGlobalIndex() + 1 << " global coarse cells" );

  // finish filling the partition array for ghosted fine cells
  for( localIndex ic = fineCellManager.numOwnedObjects(); ic < fineCellManager.size(); ++ic )
  {
    coarseLocalIndex[ic] = m_cellManager.globalToLocalMap( coarseGlobalIndex[ic] );
  }

  m_cellManager.registerWrapper< meshData::FineCellLocalIndices::type >( meshData::FineCellLocalIndices::key() );
  coarseGridConstruction::buildFineCellLists( fineCellManager, m_cellManager );
  coarseGridConstruction::fillBasicCellData( fineCellManager, m_cellManager );

  // Populate neighbor data and sets
  std::vector< int > neighborRanks = m_domain->getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    m_cellManager.addNeighbor( rank );
  }
  meshUtils::copySets( fineMesh.cellManager(),
                       meshData::CoarseCellLocalIndex::key(),
                       m_cellManager );
  meshUtils::copyNeighborData( fineMesh.cellManager(),
                               meshData::CoarseCellLocalIndex::key(),
                               neighborRanks,
                               m_cellManager,
                               meshUtils::filterArrayUnique< localIndex > );
}

void MultiscaleMeshLevel::buildCoarseNodes( MultiscaleMeshLevel & fineMesh,
                                            arrayView1d< string const > const & boundaryNodeSets )
{
  GEOSX_MARK_FUNCTION;

  MultiscaleMeshObjectManager & fineNodeManager = fineMesh.nodeManager();
  MultiscaleMeshObjectManager & fineCellManager = fineMesh.cellManager();

  // Create the array manually in order to specify default capacity (no suitable post-resize facility exists)
  ArrayOfSets< globalIndex > & nodeToSubdomain =
    fineNodeManager.registerWrapper< meshData::NodeToCoarseSubdomain::type >( meshData::NodeToCoarseSubdomain::key() ).reference();

  // Build and sync an adjacency map of fine nodes to coarse subdomains (including global boundaries)
  nodeToSubdomain = coarseGridConstruction::buildFineNodeToGlobalSubdomainMap( fineNodeManager, fineCellManager, m_cellManager, boundaryNodeSets );
  array1d< string > fields;
  fields.emplace_back( meshData::NodeToCoarseSubdomain::key() );
  CommunicationTools::getInstance().synchronizeFields( fields, fineNodeManager, m_domain->getNeighbors(), false );

  // Find all locally present coarse nodes
  array1d< localIndex > const coarseNodes =
    coarseGridConstruction::findCoarseNodes( fineNodeManager, fineCellManager, nodeToSubdomain.toViewConst() );

  // Reorder them so as to have all local nodes precede ghosted
  arrayView1d< integer const > const nodeGhostRank = fineNodeManager.ghostRank();
  auto localEnd = std::partition( coarseNodes.begin(), coarseNodes.end(),
                                  [=]( localIndex const inf ){ return nodeGhostRank[inf] < 0; } );

  localIndex const numLocalCoarseNodes = std::distance( coarseNodes.begin(), localEnd );
  localIndex const numPresentCoarseNodes = coarseNodes.size();
  globalIndex const rankOffset = MpiWrapper::prefixSum< globalIndex >( numLocalCoarseNodes );

  GEOSX_LOG_RANK( "Multiscale " << m_name << ": generated " << numLocalCoarseNodes << " local coarse nodes" );

  m_nodeManager.resize( numPresentCoarseNodes );
  m_nodeManager.setNumOwnedObjects( numLocalCoarseNodes );

  // Finally, build coarse node maps
  arrayView1d< localIndex > const coarseNodeLocalIndex =
    fineNodeManager.registerExtrinsicData< meshData::CoarseNodeLocalIndex >( m_name ).reference().toView();
  arrayView1d< globalIndex > const coarseNodeGlobalIndex =
    fineNodeManager.registerExtrinsicData< meshData::CoarseNodeGlobalIndex >( m_name ).reference().toView();
  arrayView1d< localIndex > const fineNodeLocalIndex =
    m_nodeManager.registerExtrinsicData< meshData::FineNodeLocalIndex >( m_name ).reference().toView();
  arrayView1d< globalIndex > const coarseNodeLocalToGlobal =
    m_nodeManager.localToGlobalMap();

  // Fill the local part
  forAll< parallelHostPolicy >( numLocalCoarseNodes, [=, coarseNodes = coarseNodes.toViewConst()]( localIndex const i )
  {
    localIndex const inf = coarseNodes[i];
    coarseNodeLocalIndex[inf] = i;
    coarseNodeGlobalIndex[inf] = rankOffset + i;
    coarseNodeLocalToGlobal[i] = rankOffset + i;
    fineNodeLocalIndex[i] = inf;
  } );

  // Sync across ranks
  string_array fieldNames;
  fieldNames.emplace_back( meshData::CoarseNodeGlobalIndex::key() );
  CommunicationTools::getInstance().synchronizeFields( fieldNames, fineNodeManager, m_domain->getNeighbors(), false );

  // Fill the ghosted part
  forAll< parallelHostPolicy >( numLocalCoarseNodes, numPresentCoarseNodes,
                                [=, coarseNodes = coarseNodes.toViewConst()]( localIndex const i )
  {
    localIndex const inf = coarseNodes[i];
    coarseNodeLocalIndex[inf] = i;
    coarseNodeLocalToGlobal[i] = coarseNodeGlobalIndex[inf];
    fineNodeLocalIndex[i] = inf;
  } );

  m_nodeManager.constructGlobalToLocalMap();
  m_nodeManager.setMaxGlobalIndex();
  coarseGridConstruction::fillBasicNodeData( fineNodeManager, m_nodeManager );

  GEOSX_LOG_RANK_0( "Multiscale " << m_name << ": generated " << m_nodeManager.maxGlobalIndex() + 1 << " global coarse nodes" );

  // Populate neighbor data and sets
  std::vector< int > neighborRanks = m_domain->getNeighborRanks();
  for( int const rank : neighborRanks )
  {
    m_nodeManager.addNeighbor( rank );
  }
  meshUtils::copySets( fineNodeManager,
                       meshData::CoarseNodeLocalIndex::key(),
                       m_nodeManager );
  meshUtils::copyNeighborData( fineNodeManager,
                               meshData::CoarseNodeLocalIndex::key(),
                               neighborRanks,
                               m_nodeManager,
                               meshUtils::filterArray< localIndex > );
}

void MultiscaleMeshLevel::buildCoarseMesh( MultiscaleMeshLevel & fineMesh,
                                           LinearSolverParameters::Multiscale::Coarsening const & coarse_params,
                                           array1d< string > const & boundaryNodeSets )
{
  GEOSX_MARK_FUNCTION;
  m_fineMesh = &fineMesh;
  m_domain = fineMesh.m_domain;
  buildCoarseCells( fineMesh, coarse_params );
  buildCoarseNodes( fineMesh, boundaryNodeSets );
  coarseGridConstruction::buildNodeToCellMap( fineMesh.cellManager(), fineMesh.nodeManager(), m_nodeManager );
  coarseGridConstruction::buildCellToNodeMap( m_nodeManager, m_cellManager );
}

void MultiscaleMeshLevel::writeCellData( std::vector< string > const & fieldNames ) const
{
  if( m_fineMesh )
  {
    writeCellDataCoarse( fieldNames );
  }
  else
  {
    writeCellDataFine( fieldNames );
  }
}

void MultiscaleMeshLevel::writeNodeData( std::vector< string > const & fieldNames ) const
{
  if( m_fineMesh )
  {
    writeNodeDataCoarse( fieldNames );
  }
  else
  {
    writeNodeDataFine( fieldNames );
  }
}

void MultiscaleMeshLevel::writeCellDataFine( std::vector< string > const & fieldNames ) const
{
  GEOSX_ASSERT( m_sourceMesh != nullptr );
  arrayView1d< localIndex const > const origRegion    = m_cellManager.getExtrinsicData< meshData::OrigElementRegion >();
  arrayView1d< localIndex const > const origSubRegion = m_cellManager.getExtrinsicData< meshData::OrigElementSubRegion >();
  arrayView1d< localIndex const > const origIndex     = m_cellManager.getExtrinsicData< meshData::OrigElementIndex >();

  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_cellManager.getWrapperBase( fieldName );
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), false, [&]( auto array )
    {
      using ArrayType = decltype( array );
      using ArrayViewType = typename ArrayType::ParentClass;

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const sourceFieldName = m_name + '_' + wrapper.getName();
      m_sourceMesh->getElemManager().forElementSubRegions( m_regions, [&]( localIndex, ElementSubRegionBase & subRegion )
      {
        subRegion.registerWrapper< ArrayType >( sourceFieldName ).copyWrapperAttributes( wrapper );
        auto & dstField = subRegion.getReference< ArrayType >( sourceFieldName );

        // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
        auto dims = srcField.dimsArray();
        dims[0] = dstField.size( 0 );
        dstField.resize( ArrayType::NDIM, dims.data );
      } );

      auto accessor = m_sourceMesh->getElemManager().constructViewAccessor< ArrayType, ArrayViewType >( sourceFieldName );
      forAll< parallelHostPolicy >( m_cellManager.size(), [=, dstField = accessor.toNestedView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField[origRegion[ic]][origSubRegion[ic]]( origIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    } );
  }
}

void MultiscaleMeshLevel::writeCellDataCoarse( std::vector< string > const & fieldNames ) const
{
  GEOSX_ASSERT( m_fineMesh != nullptr );
  arrayView1d< localIndex const > const coarseCellIndex = m_fineMesh->cellManager().getExtrinsicData< meshData::CoarseCellLocalIndex >();

  std::vector< string > fineFieldNames;
  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_cellManager.getWrapperBase( fieldName );
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), false, [&]( auto array )
    {
      using ArrayType = decltype( array );
      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const fineFieldName = m_name + '_' + wrapper.getName();
      fineFieldNames.push_back( fineFieldName );
      m_fineMesh->cellManager().registerWrapper< ArrayType >( fineFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = m_fineMesh->cellManager().getReference< ArrayType >( fineFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( dstField.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( dstField[ ic ],
                                              [&]( auto & dstVal, auto const ... indices )
        {
          dstVal = srcField( coarseCellIndex[ic], indices ... );
        } );
      } );
    } );
  }

  // Recursively call on finer levels
  m_fineMesh->writeCellData( fineFieldNames );
}

void MultiscaleMeshLevel::writeNodeDataFine( std::vector< string > const & fieldNames ) const
{
  GEOSX_ASSERT( m_sourceMesh != nullptr );
  arrayView1d< localIndex const > const origNodeIndex = m_nodeManager.getExtrinsicData< meshData::OrigNodeIndex >();

  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_nodeManager.getWrapperBase( fieldName );
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), false, [&]( auto array )
    {
      using ArrayType = decltype( array );

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const sourceFieldName = m_name + '_' + wrapper.getName();
      m_sourceMesh->getNodeManager().registerWrapper< ArrayType >( sourceFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = m_sourceMesh->getNodeManager().getReference< ArrayType >( sourceFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( m_nodeManager.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField( origNodeIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    } );
  }
}

void MultiscaleMeshLevel::writeNodeDataCoarse( std::vector< string > const & fieldNames ) const
{
  GEOSX_ASSERT( m_fineMesh != nullptr );
  arrayView1d< localIndex const > const fineNodeIndex = m_nodeManager.getExtrinsicData< meshData::FineNodeLocalIndex >();

  std::vector< string > fineFieldNames;
  for( string const & fieldName : fieldNames )
  {
    WrapperBase const & wrapper = m_nodeManager.getWrapperBase( fieldName );
    types::dispatch( types::StandardArrays{}, wrapper.getTypeId(), false, [&]( auto array )
    {
      using ArrayType = decltype( array );

      auto const & typedWrapper = dynamicCast< Wrapper< ArrayType > const & >( wrapper );
      auto const & srcField = typedWrapper.reference().toViewConst();

      string const fineFieldName = m_name + '_' + wrapper.getName();
      fineFieldNames.push_back( fineFieldName );
      m_fineMesh->nodeManager().registerWrapper< ArrayType >( fineFieldName ).copyWrapperAttributes( wrapper );
      auto & dstField = m_fineMesh->nodeManager().getReference< ArrayType >( fineFieldName );

      // Hack to resize all dimensions of dstField correctly, depends on impl details of Array (dimsArray)
      auto dims = srcField.dimsArray();
      dims[0] = dstField.size( 0 );
      dstField.resize( ArrayType::NDIM, dims.data );

      forAll< parallelHostPolicy >( m_nodeManager.size(), [=, dstField = dstField.toView()]( localIndex const ic )
      {
        LvArray::forValuesInSliceWithIndices( srcField[ ic ],
                                              [&]( auto const & sourceVal, auto const ... indices )
        {
          dstField( fineNodeIndex[ic], indices ... ) = sourceVal;
        } );
      } );
    } );
  }

  // Recursively call on finer levels
  m_fineMesh->writeNodeData( fineFieldNames );
}

} // namespace multiscale
} // namespace geosx
