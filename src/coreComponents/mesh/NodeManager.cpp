/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file NodeManager.cpp
 */

#include "NodeManager.hpp"
#include "FaceManager.hpp"
#include "EdgeManager.hpp"
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{

using namespace dataRepository;

// *********************************************************************************************************************
/**
 * @return
 */
//START_SPHINX_REFPOS_REG
NodeManager::NodeManager( string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  //END_SPHINX_REFPOS_REG
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::faceListString(), &m_toFacesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString(), &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString(), &elementList() );

}


NodeManager::~NodeManager()
{}


void NodeManager::resize( localIndex const newSize )
{
  m_toFacesRelation.resize( newSize, 2 * getFaceMapOverallocation() );
  m_toEdgesRelation.resize( newSize, 2 * getEdgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementSubRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * getElemMapOverAllocation() );
  ObjectManagerBase::resize( newSize );
}


void NodeManager::setNodesInformation( CellBlockManagerABC const & cellBlockManager )
{
  this->resize( cellBlockManager.numNodes() );

  m_referencePosition = cellBlockManager.getNodesPositions();

  m_localToGlobalMap = cellBlockManager.getNodeLocalToGlobal();

  for( const auto & nameArray: cellBlockManager.getNodeSets() )
  {
    auto & array = m_sets.registerWrapper< SortedArray< localIndex > >( nameArray.first ).reference();
    array = nameArray.second;
  }
}

void NodeManager::setEdgeMaps( CellBlockManagerABC const & cellBlockManager, EdgeManager const & edgeManager )
{
  GEOSX_MARK_FUNCTION;

  m_toEdgesRelation.base() = cellBlockManager.getNodeToEdges();

  m_toEdgesRelation.setRelatedObject( edgeManager );
}


void NodeManager::setFaceMaps( CellBlockManagerABC const & cellBlockManager, FaceManager const & faceManager )
{
  GEOSX_MARK_FUNCTION;

  m_toFacesRelation = cellBlockManager.getNodeToFaces();

  m_toFacesRelation.setRelatedObject( faceManager );
}

void NodeManager::setElementMaps( CellBlockManagerABC const & cellBlockManager, ElementRegionManager const & elementRegionManager )
{
  const localIndex totalNodeElems = cellBlockManager.numNodes();

  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;
  localIndex const numNodes = size();

  toElementList = cellBlockManager.getNodeToElements();

  // Resize the node to elem map.
  toElementRegionList.resize( 0 );
  toElementSubRegionList.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  toElementRegionList.reserve( entriesToReserve );
  toElementSubRegionList.reserve( entriesToReserve );

  // FIXME I'm not sure that there is a strong needs that `toElementRegionList` and `toElementSubRegionList`
  //       get the same capacities as `toElementList`.
  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeElems + numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );

    const localIndex numElementsPerNode = toElementList[ nodeID ].size() + getElemMapOverAllocation();
    toElementRegionList.setCapacityOfArray( nodeID, numElementsPerNode );
    toElementSubRegionList.setCapacityOfArray( nodeID, numElementsPerNode );
  }

  elementRegionManager.forElementSubRegionsComplete< CellElementSubRegion >(
    [&toElementRegionList, &toElementSubRegionList] ( localIndex const er, localIndex const esr, ElementRegionBase const &, CellElementSubRegion const & subRegion )
  {
    arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemToNodeMap = subRegion.nodeList();
    for( localIndex elemIdx = 0; elemIdx < subRegion.size(); ++elemIdx )
    {
      for( localIndex localNodeIdx = 0; localNodeIdx < subRegion.numNodesPerElement(); ++localNodeIdx )
      {
        // localNodeIdx is the node index in the referential of each cell (0 to 7 for a cube, e.g.).
        // While nodeIdx is the global index of the node.
        localIndex const nodeIdx = elemToNodeMap( elemIdx, localNodeIdx );
        toElementRegionList.emplaceBack( nodeIdx, er );
        toElementSubRegionList.emplaceBack( nodeIdx, esr );
      }
    }
  } );

  this->m_toElements.setElementRegionManager( elementRegionManager );
}

void NodeManager::compressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toFacesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}


void NodeManager::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementSubRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString() ));

  if( this->hasWrapper( "usedFaces" ) )
  {
    exclusionList.insert( this->getWrapperIndex( "usedFaces" ));
  }
}


localIndex NodeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}


localIndex NodeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}


template< bool DOPACK >
localIndex NodeManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toEdgesRelation.toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToEdges,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_toEdgesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toFacesRelation.toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToFaces,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_toFacesRelation.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );
  return packedSize;
}


localIndex NodeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const )
{
  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::edgeListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toEdgesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_toEdgesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::faceListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}


void NodeManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toEdgesRelation,
                                    m_toEdgesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation,
                                    m_toFacesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );

}


void NodeManager::depopulateUpMaps( std::set< localIndex > const & receivedNodes,
                                    array2d< localIndex > const & edgesToNodes,
                                    ArrayOfArraysView< localIndex const > const & facesToNodes,
                                    ElementRegionManager const & elemRegionManager )
{

  ObjectManagerBase::cleanUpMap( receivedNodes, m_toEdgesRelation.toView(), edgesToNodes );
  ObjectManagerBase::cleanUpMap( receivedNodes, m_toFacesRelation.toView(), facesToNodes );

  for( auto const & targetIndex : receivedNodes )
  {
    std::set< std::tuple< localIndex, localIndex, localIndex > > eraseList;
    for( localIndex k=0; k<m_toElements.m_toElementRegion.sizeOfArray( targetIndex ); ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      CellElementSubRegion const & subRegion = elemRegionManager.getRegion( elemRegionIndex ).
                                                 getSubRegion< CellElementSubRegion >( elemSubRegionIndex );
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const downmap = subRegion.nodeList();
      bool hasTargetIndex = false;

      for( localIndex a=0; a<downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap( elemIndex, a );
        if( compositeLocalIndex==targetIndex )
        {
          hasTargetIndex=true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert( std::make_tuple( elemRegionIndex, elemSubRegionIndex, elemIndex ) );
      }
    }
    for( auto const & val : eraseList )
    {
      erase( m_toElements, targetIndex, std::get< 0 >( val ), std::get< 1 >( val ), std::get< 2 >( val ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, NodeManager, string const &, Group * const )

}
