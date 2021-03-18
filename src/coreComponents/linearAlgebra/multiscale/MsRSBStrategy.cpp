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
 * @file MsRSBStrategy.cpp
 */

#include "MsRSBStrategy.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshData.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshUtils.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/utilities/TransposeOperator.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
MsRSBStrategy< LAI >::MsRSBStrategy( string name,
                                     LinearSolverParameters::Multiscale params )
  : MultiscaleStrategy< LAI >( std::move( name ), std::move( params ) ),
  m_mesh( m_name )
{}

template< typename LAI >
std::unique_ptr< LinearOperator< typename LAI::ParallelVector > >
makeRestriction( LinearSolverParameters::Multiscale const & params,
                 typename LAI::ParallelMatrix const & prolongation )
{
  std::unique_ptr< LinearOperator< typename LAI::ParallelVector > > restriction;
  if( params.galerkin )
  {
    // Make a transpose operator with a reference to P, which will be computed later
    restriction = std::make_unique< TransposeOperator< LAI > >( prolongation );
  }
  else
  {
    // Make an explicit transpose of tentative (initial) P
    typename LAI::ParallelMatrix R;
    prolongation.transpose( R );
    restriction = std::make_unique< typename LAI::ParallelMatrix >( std::move( R ) );
  }
  return restriction;
}

template< typename Matrix >
Matrix makeJacobiIterationMatrix( Matrix const & fine_mat,
                                  localIndex const numComp,
                                  real64 const omega )
{
  GEOSX_MARK_FUNCTION;

  // 1. Apply SC approximation
  Matrix iterMat = LAIHelperFunctions::separateComponentFilter( fine_mat, numComp );

  // 2. Compute -w * Dinv * Asc;
  typename Matrix::Vector diag;
  diag.createWithLocalSize( fine_mat.numLocalRows(), fine_mat.getComm() );
  iterMat.extractDiagonal( diag );
  diag.reciprocal();
  diag.scale( -omega );
  iterMat.leftScale( diag );

  // 3. Compute I - w * Dinv * Asc
  diag.set( 1.0 );
  iterMat.addDiagonal( diag );
  return iterMat;
}

ArrayOfSets< localIndex >
buildFineNodeToLocalSubdomainMap( MultiscaleMeshObjectManager const & fineNodeManager,
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
  ArrayOfSets< localIndex > nodeToSubdomain;
  nodeToSubdomain.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  arrayView1d< localIndex const > const coarseCellLocalIndex = fineCellManager.getExtrinsicData< meshData::CoarseCellLocalIndex >();
  forAll< parallelHostPolicy >( fineNodeManager.size(), [=, nodeToSub = nodeToSubdomain.toView()]( localIndex const inf )
  {
    for( localIndex const icf : nodeToCell[inf] )
    {
      nodeToSub.insertIntoSet( inf, coarseCellLocalIndex[icf] );
    }
  } );
  localIndex numSubdomains = coarseCellManager.size();
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

ArrayOfSets< localIndex >
buildCoarseNodeToLocalSubdomainMap( MultiscaleMeshObjectManager const & coarseNodeManager,
                                    MultiscaleMeshObjectManager const & coarseCellManager,
                                    arrayView1d< string const > const & boundaryNodeSets )
{
  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = coarseNodeManager.toDualRelation().toViewConst();

  // count the row lengths
  array1d< localIndex > rowCounts( coarseNodeManager.size() );
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, rowCounts = rowCounts.toView()]( localIndex const inc )
  {
    rowCounts[inc] = nodeToCell.sizeOfSet( inc );
  } );
  for( string const & setName : boundaryNodeSets )
  {
    SortedArrayView< localIndex const > const set = coarseNodeManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, rowCounts = rowCounts.toView()]( localIndex const i )
    {
      ++rowCounts[set[i]];
    } );
  }

  // Resize from row lengths
  ArrayOfSets< localIndex > nodeToSubdomain;
  nodeToSubdomain.resizeFromCapacities< parallelHostPolicy >( rowCounts.size(), rowCounts.data() );

  // Fill the map
  forAll< parallelHostPolicy >( coarseNodeManager.size(), [=, nodeToSub = nodeToSubdomain.toView()]( localIndex const inc )
  {
    nodeToSub.insertIntoSet( inc, nodeToCell[inc].begin(), nodeToCell[inc].end() );
  } );
  localIndex numSubdomains = coarseCellManager.size();
  for( string const & setName : boundaryNodeSets )
  {
    SortedArrayView< localIndex const > const set = coarseNodeManager.getSet( setName ).toViewConst();
    forAll< parallelHostPolicy >( set.size(), [=, nodeToSub = nodeToSubdomain.toView()]( localIndex const inc )
    {
      nodeToSub.insertIntoSet( set[inc], numSubdomains );
    } );
    ++numSubdomains;
  }

  return nodeToSubdomain;
}

/**
 * @brief Build the basic sparsity pattern for prolongation.
 *
 * Support of a coarse nodal basis function is defined as the set of fine-scale nodes
 * that are adjacent exclusively to subdomains (coarse cells or boundaries) that are
 * also adjacent to that coarse node.
 */
ArrayOfSets< localIndex >
buildNodalSupports( MultiscaleMeshLevel const & fine,
                    MultiscaleMeshLevel const & coarse,
                    arrayView1d< string const > const & boundaryNodeSets )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfSets< localIndex > const fineNodeToCoarseCell =
    buildFineNodeToLocalSubdomainMap( fine.nodeManager(), fine.cellManager(), coarse.cellManager(), {} );
  ArrayOfSets< localIndex > const fineNodeToSubdomain =
    buildFineNodeToLocalSubdomainMap( fine.nodeManager(), fine.cellManager(), coarse.cellManager(), boundaryNodeSets );
  ArrayOfSets< localIndex > const coarseNodeToSubdomain =
    buildCoarseNodeToLocalSubdomainMap( coarse.nodeManager(), coarse.cellManager(), boundaryNodeSets );

  ArrayOfSetsView< localIndex const > const coarseCellToNode = coarse.cellManager().toDualRelation().toViewConst();
  arrayView1d< localIndex const > const coarseNodeIndex = fine.nodeManager().getExtrinsicData< meshData::CoarseNodeLocalIndex >();

  // Algorithm:
  // Loop over all fine nodes.
  // If node is a coarse node, immediately assign to its own support.
  // Otherwise, get a list of adjacent coarse cells.
  // If list is length 1, assign the node to supports of all coarse nodes adjacent to that coarse cell.
  // Otherwise, collect a unique list of candidate coarse nodes by visiting them through coarse cells.
  // For each candidate, check that fine node's subdomain list is included in the candidates subdomain list.
  // Otherwise, discard the candidate.
  //
  // All above is done twice: once to count (or get upper bound on) row lengths, once to actually build supports.
  // For the last case, don't need to check inclusion when counting, just use number of candidates as upper bound.

  // Count row lengths
  array1d< localIndex > rowLengths( fine.nodeManager().size() );
  forAll< parallelHostPolicy >( fine.nodeManager().size(), [coarseNodeIndex, coarseCellToNode,
                                                            rowLengths = rowLengths.toView(),
                                                            fineNodeToSubdomain = fineNodeToSubdomain.toViewConst(),
                                                            fineNodeToCoarseCell = fineNodeToCoarseCell.toViewConst()]( localIndex const inf )
  {
    if( coarseNodeIndex[inf] >= 0 )
    {
      rowLengths[inf] = 1;
    }
    else if( fineNodeToSubdomain.sizeOfSet( inf ) == 1 )
    {
      rowLengths[inf] = coarseCellToNode.sizeOfSet( fineNodeToSubdomain( inf, 0 ) );
    }
    else
    {
      localIndex numCoarseNodes = 0;
      meshUtils::forUniqueNeighbors< 128 >( inf, fineNodeToCoarseCell, coarseCellToNode, [&]( localIndex )
      {
        ++numCoarseNodes;
      } );
      rowLengths[inf] = numCoarseNodes;
    }
  } );

  // Resize
  ArrayOfSets< localIndex > supports;
  supports.resizeFromCapacities< parallelHostPolicy >( rowLengths.size(), rowLengths.data() );

  // Fill the map
  forAll< parallelHostPolicy >( fine.nodeManager().size(), [coarseNodeIndex, coarseCellToNode,
                                                            supports = supports.toView(),
                                                            fineNodeToSubdomain = fineNodeToSubdomain.toViewConst(),
                                                            fineNodeToCoarseCell = fineNodeToCoarseCell.toViewConst(),
                                                            coarseNodeToSubdomain = coarseNodeToSubdomain.toViewConst()]
                                  ( localIndex const inf )
  {
    if( coarseNodeIndex[inf] >= 0 )
    {
      supports.insertIntoSet( inf, coarseNodeIndex[inf] );
    }
    else if( fineNodeToSubdomain.sizeOfSet( inf ) == 1 )
    {
      arraySlice1d< localIndex const > const coarseNodes = coarseCellToNode[fineNodeToSubdomain( inf, 0 )];
      supports.insertIntoSet( inf, coarseNodes.begin(), coarseNodes.end() );
    }
    else
    {
      arraySlice1d< localIndex const > const fsubs = fineNodeToSubdomain[inf];
      meshUtils::forUniqueNeighbors< 128 >( inf, fineNodeToCoarseCell, coarseCellToNode, [&]( localIndex const inc )
      {
        arraySlice1d< localIndex const > const csubs = coarseNodeToSubdomain[inc];
        if( std::includes( csubs.begin(), csubs.end(), fsubs.begin(), fsubs.end() ) )
        {
          supports.insertIntoSet( inf, inc );
        }
      } );
    }
  } );

  return supports;
}

SparsityPattern< globalIndex >
buildProlongationSparsity( MultiscaleMeshObjectManager const & fineManager,
                           MultiscaleMeshObjectManager const & coarseManager,
                           ArrayOfSetsView< localIndex const > const & supports,
                           integer const numDofComp )
{
  GEOSX_MARK_FUNCTION;

  // This assumes an SDC-type pattern (i.e. no coupling between dof components on the same node)
  array1d< localIndex > rowLengths( numDofComp * fineManager.numOwnedObjects() );
  forAll< parallelHostPolicy >( fineManager.numOwnedObjects(), [=, rowLengths = rowLengths.toView()]( localIndex const k )
  {
    for( localIndex i = 0; i < numDofComp; ++i )
    {
      rowLengths[k * numDofComp + i] = supports.sizeOfSet( k );
    }
  } );

  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( rowLengths.size(),
                                                         numDofComp * ( coarseManager.maxGlobalIndex() + 1 ),
                                                         rowLengths.data() );

  arrayView1d< globalIndex const > const coarseLocalToGlobal = coarseManager.localToGlobalMap();
  forAll< parallelHostPolicy >( fineManager.numOwnedObjects(), [=, pattern = pattern.toView()]( localIndex const k )
  {
    localIndex const fineOffset = k * numDofComp;
    for( localIndex const inc : supports[k] )
    {
      globalIndex const coarseOffset = coarseLocalToGlobal[inc] * numDofComp;
      for( localIndex i = 0; i < numDofComp; ++i )
      {
        pattern.insertNonZero( fineOffset + i, coarseOffset + i );
      }
    }
  } );

  return pattern;
}

ArrayOfSets< localIndex >
buildLocalNodalConnectivity( MultiscaleMeshObjectManager const & nodeManager,
                             MultiscaleMeshObjectManager const & cellManager )
{
  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = nodeManager.toDualRelation().toViewConst();
  MultiscaleMeshObjectManager::MapViewConst const cellToNode = cellManager.toDualRelation().toViewConst();
  arrayView1d< integer const > const cellGhostRank = cellManager.ghostRank();

  // Collect row sizes
  array1d< localIndex > rowLength( nodeManager.size() );
  forAll< parallelHostPolicy >( nodeManager.size(), [=, rowLength = rowLength.toView()]( localIndex const k )
  {
    meshUtils::forUniqueNeighbors< 128 >( k, nodeToCell, cellToNode, cellGhostRank, [&]( localIndex const )
    {
      ++rowLength[k];
    } );
  } );

  // Resize
  ArrayOfSets< localIndex > conn( nodeManager.size() );
  conn.resizeFromCapacities< parallelHostPolicy >( rowLength.size(), rowLength.data() );

  // Fill the map
  forAll< parallelHostPolicy >( nodeManager.size(), [=, conn = conn.toView()]( localIndex const k )
  {
    meshUtils::forUniqueNeighbors< 128 >( k, nodeToCell, cellToNode, cellGhostRank, [&]( localIndex const n )
    {
      conn.insertIntoSet( k, n );
    } );
  } );

  return conn;
}

array1d< localIndex >
makeSeededPartition( ArrayOfSetsView< localIndex const > const & connectivity,
                     arrayView1d< localIndex const > const & seeds,
                     ArrayOfSetsView< localIndex const > const & supports )
{
  localIndex const numParts = seeds.size();
  localIndex const numNodes = connectivity.size();

  array1d< localIndex > part( numNodes );
  part.setValues< parallelHostPolicy >( -1 );

  // Initialize the partitions and expansion front
  array1d< localIndex > front;
  front.reserve( numNodes );
  forAll< serialPolicy >( numParts, [&]( localIndex const ip )
  {
    part[seeds[ip]] = ip;
    for( localIndex const k : connectivity[seeds[ip]] )
    {
      if( part[k] < 0 && supports.contains( k, ip ) )
      {
        front.emplace_back( k );
      }
    }
  } );

  // Use AoA with 1 array for its atomic emplace capability
  // TODO: this might not be efficient due to thread contention;
  //       may need to use a serial policy instead (benchmark me)
  ArrayOfArrays< localIndex > unassigned( 1, 8 * numNodes );
  array1d< localIndex > newPart;

  integer numIter = 0;
  while( !front.empty() )
  {
    // Make the list unique
    localIndex const numFrontNodes = LvArray::sortedArrayManipulation::makeSortedUnique( front.begin(), front.end() );
    front.resize( numFrontNodes );
    newPart.resize( numFrontNodes );
    newPart.setValues< parallelHostPolicy >( -1 );

    // Pass 1: assign partitions to the front nodes based on majority among neighbors
    RAJA::ReduceSum< parallelHostReduce, localIndex > numAssn = 0;
    forAll< parallelHostPolicy >( front.size(), [=, front = front.toViewConst(),
                                                 part = part.toViewConst(),
                                                 newPart = newPart.toView()]( localIndex const i )
    {
      localIndex const k = front[i];
      localIndex maxCount = 0;
      meshUtils::forUniqueNeighborValues< 128 >( k, connectivity, part,
                                                 []( localIndex const p ){ return p >= 0; }, // only assigned nodes
                                                 [&]( localIndex const p, localIndex const count )
      {
        if( p >= 0 && count > maxCount && supports.contains( k, p ) )
        {
          newPart[i] = p;
          maxCount = count;
        }
      } );

      if( maxCount > 0 )
      {
        numAssn += 1;
      }
    } );

    // Terminate the loop as soon as no new assignments are made
    if( numAssn.get() == 0 )
    {
      break;
    }

    // Pass 2: copy new front partition assignments back into full partition array
    forAll< parallelHostPolicy >( front.size(), [front = front.toViewConst(),
                                                 newPart = newPart.toViewConst(),
                                                 part = part.toView()] ( localIndex const i )
    {
      part[front[i]] = newPart[i];
    } );

    // Pass 3: build a list of unassigned neighbor indices to become the new front
    forAll< parallelHostPolicy >( front.size(), [=, front = front.toViewConst(),
                                                 part = part.toViewConst(),
                                                 unassigned = unassigned.toView()]( localIndex const i )
    {
      localIndex const k = front[i];
      meshUtils::forUniqueNeighborValues< 128 >( k, connectivity,
                                                 []( localIndex const _ ){ return _; }, // just unique neighbor indices
                                                 [&]( localIndex const n ){ return part[n] < 0; }, // only unassigned nodes
                                                 [&]( localIndex const n )
      {
        unassigned.emplaceBackAtomic< parallelHostAtomic >( 0, n );
      } );
    } );

    // Make the new front expansion list
    front.clear();
    front.insert( 0, unassigned[0].begin(), unassigned[0].end() );
    unassigned.clearArray( 0 );

    ++numIter;
  }

  if( !front.empty() )
  {
    for( localIndex const k : front )
    {
      std::ostringstream os;
      os << "Unassigned node : " << k << "\n  Belongs to supports of : " << supports[k] << "\n  Neighbor nodes : ";
      for( localIndex n : connectivity[k] )
      {
        os << n << '[' << part[n] << "], ";
      }
      GEOSX_LOG_RANK( oss.str() );
    }
  }

  return part;
}

template< typename Matrix >
Matrix buildTentativeProlongation( MultiscaleMeshLevel & fineMesh,
                                   MultiscaleMeshLevel const & coarseMesh,
                                   integer const numDofComp,
                                   arrayView1d< string const > const & boundarySets,
                                   MPI_Comm const & comm )
{
  GEOSX_MARK_FUNCTION;

  // Build support regions and tentative prolongation
  ArrayOfSets< localIndex > const supports = buildNodalSupports( fineMesh, coarseMesh, boundarySets );
  ArrayOfSets< localIndex > const nodalConn = buildLocalNodalConnectivity( fineMesh.nodeManager(), fineMesh.cellManager() );
  arrayView1d< localIndex const > const coarseNodes = coarseMesh.nodeManager().getExtrinsicData< meshData::FineNodeLocalIndex >().toViewConst();
  array1d< localIndex > const initPart = makeSeededPartition( nodalConn.toViewConst(), coarseNodes, supports.toViewConst() );

  SparsityPattern< globalIndex > localPattern =
    buildProlongationSparsity( fineMesh.nodeManager(), coarseMesh.nodeManager(), supports.toViewConst(), numDofComp );

  // Construct the tentative prolongation, consuming the sparsity pattern
  CRSMatrix< real64, globalIndex > localMatrix;
  localMatrix.assimilate< parallelHostPolicy >( std::move( localPattern ) );

  // Add initial unity values
  arrayView1d< globalIndex const > const coarseLocalToGlobal = coarseMesh.nodeManager().localToGlobalMap();
  forAll< parallelHostPolicy >( fineMesh.nodeManager().numOwnedObjects(), [=, localMatrix = localMatrix.toViewConstSizes()]( localIndex const inf )
  {
    if( initPart[inf] >= 0 )
    {
      real64 const value = 1.0;
      for( integer c = 0; c < numDofComp; ++c )
      {
        globalIndex const col = coarseLocalToGlobal[initPart[inf]] * numDofComp + c;
        localMatrix.addToRow< serialAtomic >( inf * numDofComp + c, &col, &value, 1 );
      }
    }
  } );

  Matrix mat;
  mat.create( localMatrix.toViewConst(), coarseMesh.nodeManager().numOwnedObjects() * numDofComp, comm );
  return mat;
}

template< typename Matrix >
void plotProlongation( Matrix const & prolongation,
                       string const & prefix,
                       integer const numDofComp,
                       MultiscaleMeshLevel & fineMesh )
{
  std::vector< string > bNames{ "X ", "Y ", "Z " };
  std::vector< string > cNames{ " x", " y", " z" };

  globalIndex const numNodes = prolongation.numGlobalCols() / numDofComp;
  int const paddedSize = LvArray::integerConversion< int >( std::to_string( numNodes ).size() );

  std::vector< arrayView3d< real64 > > views;
  std::vector< string > names;

  for( globalIndex icn = 0; icn < numNodes; ++icn )
  {
    string const name = prefix + "_P_" + stringutilities::padValue( icn, paddedSize );
    auto arr = std::make_unique< array3d< real64 > >( fineMesh.nodeManager().size(), numDofComp, numDofComp );
    auto & wrapper = fineMesh.nodeManager().registerWrapper( name, std::move( arr ) ).
      setPlotLevel( dataRepository::PlotLevel::LEVEL_0 ).
      setDimLabels( 1, { bNames.begin(), bNames.begin() + numDofComp } ).
      setDimLabels( 2, { cNames.begin(), cNames.begin() + numDofComp } );
    views.push_back( wrapper.referenceAsView() );
    names.push_back( name );
  }

  array1d< globalIndex > colIndices;
  array1d< real64 > values;

  for( localIndex localRow = 0; localRow < prolongation.numLocalRows(); ++localRow )
  {
    globalIndex const globalRow = prolongation.ilower() + localRow;
    localIndex const rowNode = localRow / numDofComp;
    integer const rowComp = static_cast< integer >( localRow % numDofComp );

    localIndex const numValues = prolongation.rowLength( globalRow );
    colIndices.resize( numValues );
    values.resize( numValues );
    prolongation.getRowCopy( globalRow, colIndices, values );

    for( localIndex i = 0; i < numValues; ++i )
    {
      globalIndex const colNode = colIndices[i] / numDofComp;
      integer const colComp = static_cast< integer >( colIndices[i] % numDofComp );
      views[colNode]( rowNode, colComp, rowComp ) = values[i];
    }
  }

  fineMesh.writeNodeData( names );
  for( string const & name : names )
  {
    fineMesh.nodeManager().deregisterWrapper( name );
  }
}

template< typename LAI >
void MsRSBStrategy< LAI >::initializeCoarseLevel( MultiscaleStrategy< LAI > & fine_level )
{
  GEOSX_MARK_FUNCTION;

  MsRSBStrategy< LAI > & fine = dynamicCast< MsRSBStrategy< LAI > & >( fine_level );
  m_numComp = fine.m_numComp;
  m_location = fine.m_location;

  // Coarsen mesh
  m_mesh.buildCoarseMesh( fine.mesh(), m_params.coarsening, m_params.boundarySets );

  // Write data back to GEOSX for visualization and debug
  if( m_params.debugLevel >= 1 )
  {
    m_mesh.writeCellData( { ObjectManagerBase::viewKeyStruct::ghostRankString() } );
    m_mesh.writeNodeData( { meshData::FineNodeLocalIndex::key() } );
    fine.mesh().writeCellData( { meshData::CoarseCellLocalIndex::key(),
                                 meshData::CoarseCellGlobalIndex::key(),
                                 ObjectManagerBase::viewKeyStruct::ghostRankString() } );
    fine.mesh().writeNodeData( { meshData::CoarseNodeLocalIndex::key(),
                                 meshData::CoarseNodeGlobalIndex::key(),
                                 ObjectManagerBase::viewKeyStruct::ghostRankString() } );
  }

  // Create a "fake" coarse matrix (no data, just correct sizes/comms)
  MultiscaleMeshObjectManager const & mgr = m_location == DofManager::Location::Node ? m_mesh.nodeManager() : m_mesh.cellManager();
  localIndex const numLocalDof = mgr.numOwnedObjects() * m_numComp;
  m_matrix.createWithLocalSize( numLocalDof, numLocalDof, 0, fine.matrix().getComm() );

  // Build initial (tentative) prolongation operator
  m_prolongation = buildTentativeProlongation< Matrix >( fine.mesh(), m_mesh, m_numComp, m_params.boundarySets, fine.matrix().getComm() );
  if( m_params.debugLevel >= 2 )
  {
    m_prolongation.write( m_name + "_P_init.mat", LAIOutputFormat::MATRIX_MARKET );
    plotProlongation( m_prolongation, m_name, m_numComp, fine.mesh() );
  }

  m_restriction = makeRestriction< LAI >( m_params, m_prolongation );
}

template< typename LAI >
void MsRSBStrategy< LAI >::initializeFineLevel( MeshLevel & mesh,
                                                DofManager const & dofManager,
                                                string const & fieldName,
                                                MPI_Comm const & comm )
{
  GEOSX_MARK_FUNCTION;

  m_numComp = dofManager.numComponents( fieldName );
  m_location = dofManager.location( fieldName );
  m_mesh.buildFineMesh( mesh, dofManager.regions( fieldName ) );

  // Create a "fake" fine matrix (no data, just correct sizes/comms for use at coarse level init)
  localIndex const numLocalDof = dofManager.numLocalDofs( fieldName );
  m_matrix.createWithLocalSize( numLocalDof, numLocalDof, 0, comm );
}

template< typename LAI >
void MsRSBStrategy< LAI >::compute( Matrix const & fine_mat )
{
  GEOSX_MARK_FUNCTION;

  Matrix const iterMat = makeJacobiIterationMatrix( fine_mat, m_numComp, m_params.msrsb.relaxation );

  Matrix Pnew( m_prolongation );
  Vector rescale;
  rescale.createWithLocalSize( m_prolongation.numLocalRows(), m_prolongation.getComm() );
//  for( integer iter = 0; iter < m_params.msrsb.maxSmoothingIter; ++iter )
//  {
//    // TODO: multiply while preserving sparsity pattern of output matrix
//    iterMat.multiply( m_prolongation, Pnew );
//    // TODO: compute update norm, check convergence
//  }

  /// Compute coarse operator
  if( m_params.galerkin )
  {
    fine_mat.multiplyPtAP( m_prolongation, m_matrix );
  }
  else
  {
    Matrix const & restriction = dynamicCast< Matrix const & >( *m_restriction );
    fine_mat.multiplyRAP( m_prolongation, restriction, m_matrix );
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MsRSBStrategy< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MsRSBStrategy< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MsRSBStrategy< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geosx
