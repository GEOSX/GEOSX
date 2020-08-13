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
NodeManager::NodeManager( std::string const & name, Group * const parent ) :
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString, &m_referencePosition );
  //END_SPHINX_REFPOS_REG
  this->registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::faceListString, &m_toFacesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString,
                         &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString,
                         &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString, &elementList() );
}

NodeManager::~NodeManager()
{}

void
NodeManager::resize( localIndex const newSize )
{
  m_toFacesRelation.resize( newSize, 2 * getFaceMapOverallocation() );
  m_toEdgesRelation.resize( newSize, 2 * getEdgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementSubRegion.resize( newSize,
                                            2 * getElemMapOverAllocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * getElemMapOverAllocation() );
  ObjectManagerBase::resize( newSize );
}

void
NodeManager::SetEdgeMaps( EdgeManager const * const edgeManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const & edgeToNodeMap = edgeManager->nodeList();
  localIndex const numEdges = edgeToNodeMap.size( 0 );
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toEdgesTemp( numNodes, edgeManager->maxEdgesPerNode() );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeEdges = 0;

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID ) {
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 0 ),
                                                         edgeID );
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 1 ),
                                                         edgeID );
    totalNodeEdges += 2;
  } );

  // Resize the node to edge map.
  m_toEdgesRelation.resize( 0 );

  // Reserve space for the number of current nodes plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  m_toEdgesRelation.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() +
    numNodes * getEdgeMapOverallocation() * ( 1 + 2 * overAllocationFactor );
  m_toEdgesRelation.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toEdgesRelation.appendSet( toEdgesTemp.sizeOfArray( nodeID ) +
                                 getEdgeMapOverallocation() );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_toEdgesRelation.toView();
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    localIndex * const edges = toEdgesTemp[nodeID];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueEdges =
      LvArray::sortedArrayManipulation::makeSortedUnique( edges,
                                                          edges + numNodeEdges );
    toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges );
  } );

  m_toEdgesRelation.SetRelatedObject( edgeManager );
}

void
NodeManager::SetFaceMaps( FaceManager const * const faceManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArraysView< localIndex const > const & faceToNodes =
    faceManager->nodeList().toViewConst();
  localIndex const numFaces = faceToNodes.size();
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toFacesTemp( numNodes, faceManager->maxFacesPerNode() );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeFaces = 0;

  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID ) {
    localIndex const numFaceNodes = faceToNodes.sizeOfArray( faceID );
    totalNodeFaces += numFaceNodes;
    for( localIndex a = 0; a < numFaceNodes; ++a )
    {
      toFacesTemp.emplaceBackAtomic< parallelHostAtomic >( faceToNodes( faceID, a ),
                                                           faceID );
    }
  } );

  // Resize the node to face map.
  m_toFacesRelation.resize( 0 );

  // Reserve space for the number of nodes faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  m_toFacesRelation.reserve( entriesToReserve );

  // Reserve space for the total number of node faces + extra space for existing nodes + even more space for new nodes.
  localIndex const valuesToReserve = totalNodeFaces.get() +
    numNodes * FaceManager::nodeMapExtraSpacePerFace() *
      ( 1 + 2 * overAllocationFactor );
  m_toFacesRelation.reserveValues( valuesToReserve );

  // Append the individual arrays.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toFacesRelation.appendSet( toFacesTemp.sizeOfArray( nodeID ) +
                                 getFaceMapOverallocation() );
  }

  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    localIndex * const faces = toFacesTemp[nodeID];
    localIndex const numNodeFaces = toFacesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueFaces =
      LvArray::sortedArrayManipulation::makeSortedUnique( faces,
                                                          faces + numNodeFaces );
    m_toFacesRelation.insertIntoSet( nodeID, faces, faces + numUniqueFaces );
  } );

  m_toFacesRelation.SetRelatedObject( faceManager );
}

void
NodeManager::SetElementMaps( ElementRegionManager const * const elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList =
    m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;
  localIndex const numNodes = size();

  // The number of elements attached to the each node.
  array1d< localIndex > elemsPerNode( numNodes );

  // The total number of elements, the sum of elemsPerNode.
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeElems = 0;

  elementRegionManager->forElementSubRegions< CellElementSubRegion >(
    [&elemsPerNode, &totalNodeElems]( CellElementSubRegion const & subRegion ) {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemToNodeMap =
        subRegion.nodeList();
      forAll< parallelHostPolicy >(
        subRegion.size(),
        [&elemsPerNode, totalNodeElems, &elemToNodeMap, &subRegion](
          localIndex const k ) {
          localIndex const numIndependedNodes =
            subRegion.numIndependentNodesPerElement();
          totalNodeElems += numIndependedNodes;
          for( localIndex a = 0; a < numIndependedNodes; ++a )
          {
            localIndex const nodeIndex = elemToNodeMap( k, a );
            RAJA::atomicInc< parallelHostAtomic >( &elemsPerNode[nodeIndex] );
          }
        } );
    } );

  // Resize the node to elem map.
  toElementRegionList.resize( 0 );
  toElementSubRegionList.resize( 0 );
  toElementList.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  toElementRegionList.reserve( entriesToReserve );
  toElementSubRegionList.reserve( entriesToReserve );
  toElementList.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeElems.get() +
    numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );
  toElementList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );
    toElementList.appendArray( 0 );

    toElementRegionList.setCapacityOfArray(
      nodeID,
      elemsPerNode[nodeID] + getElemMapOverAllocation() );
    toElementSubRegionList.setCapacityOfArray(
      nodeID,
      elemsPerNode[nodeID] + getElemMapOverAllocation() );
    toElementList.setCapacityOfArray(
      nodeID,
      elemsPerNode[nodeID] + getElemMapOverAllocation() );
  }

  // Populate the element maps.
  // Note that this can't be done in parallel because the three element lists must be in the same order.
  // If this becomes a bottleneck create a temporary ArrayOfArrays of tuples and insert into that first then copy over.
  elementRegionManager->forElementSubRegionsComplete< CellElementSubRegion >(
    [&toElementRegionList, &toElementSubRegionList, &toElementList](
      localIndex const er,
      localIndex const esr,
      ElementRegionBase const &,
      CellElementSubRegion const & subRegion ) {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemToNodeMap =
        subRegion.nodeList();
      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        for( localIndex a = 0; a < subRegion.numIndependentNodesPerElement(); ++a )
        {
          localIndex const nodeIndex = elemToNodeMap( k, a );
          toElementRegionList.emplaceBack( nodeIndex, er );
          toElementSubRegionList.emplaceBack( nodeIndex, esr );
          toElementList.emplaceBack( nodeIndex, k );
        }
      }
    } );

  this->m_toElements.setElementRegionManager( elementRegionManager );
}

void
NodeManager::CompressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toFacesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}

void
NodeManager::ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::elementRegionListString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::elementSubRegionListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString ) );

  if( this->hasWrapper( "usedFaces" ) )
  {
    exclusionList.insert( this->getWrapperIndex( "usedFaces" ) );
  }
}

localIndex
NodeManager::PackUpDownMapsSize(
  arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

localIndex
NodeManager::PackUpDownMaps(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex
NodeManager::PackUpDownMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString ) );
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               m_toEdgesRelation.toArrayOfArraysView(),
                               m_unmappedGlobalIndicesInToEdges,
                               packList,
                               this->localToGlobalMap(),
                               m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString ) );
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               m_toFacesRelation.toArrayOfArraysView(),
                               m_unmappedGlobalIndicesInToFaces,
                               packList,
                               this->localToGlobalMap(),
                               m_toFacesRelation.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );
  return packedSize;
}

localIndex
NodeManager::UnpackUpDownMaps( buffer_unit_type const *& buffer,
                               localIndex_array & packList,
                               bool const overwriteUpMaps,
                               bool const )
{
  localIndex unPackedSize = 0;

  string temp;
  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::edgeListString, "" );
  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_toEdgesRelation,
                       packList,
                       m_unmappedGlobalIndicesInToEdges,
                       this->globalToLocalMap(),
                       m_toEdgesRelation.RelatedObjectGlobalToLocal(),
                       overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::faceListString, "" );
  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_toFacesRelation,
                       packList,
                       m_unmappedGlobalIndicesInToFaces,
                       this->globalToLocalMap(),
                       m_toFacesRelation.RelatedObjectGlobalToLocal(),
                       overwriteUpMaps );

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString, "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void
NodeManager::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_toEdgesRelation,
                                    m_toEdgesRelation.RelatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_toFacesRelation.RelatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void
NodeManager::depopulateUpMaps(
  std::set< localIndex > const & receivedNodes,
  array2d< localIndex > const & edgesToNodes,
  ArrayOfArraysView< localIndex const > const & facesToNodes,
  ElementRegionManager const & elemRegionManager )
{
  ObjectManagerBase::CleanUpMap( receivedNodes,
                                 m_toEdgesRelation.toView(),
                                 edgesToNodes );
  ObjectManagerBase::CleanUpMap( receivedNodes,
                                 m_toFacesRelation.toView(),
                                 facesToNodes );

  for( auto const & targetIndex : receivedNodes )
  {
    std::set< std::tuple< localIndex, localIndex, localIndex > > eraseList;
    for( localIndex k = 0;
         k < m_toElements.m_toElementRegion.sizeOfArray( targetIndex );
         ++k )
    {
      localIndex const elemRegionIndex =
        m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex =
        m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex = m_toElements.m_toElementIndex[targetIndex][k];

      CellElementSubRegion const * subRegion =
        elemRegionManager.GetRegion( elemRegionIndex )
          ->GetSubRegion< CellElementSubRegion >( elemSubRegionIndex );
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & downmap =
        subRegion->nodeList();
      bool hasTargetIndex = false;

      for( localIndex a = 0; a < downmap.size( 1 ); ++a )
      {
        localIndex const compositeLocalIndex = downmap( elemIndex, a );
        if( compositeLocalIndex == targetIndex )
        {
          hasTargetIndex = true;
        }
      }
      if( !hasTargetIndex )
      {
        eraseList.insert(
          std::make_tuple( elemRegionIndex, elemSubRegionIndex, elemIndex ) );
      }
    }
    for( auto const & val : eraseList )
    {
      erase( m_toElements,
             targetIndex,
             std::get< 0 >( val ),
             std::get< 1 >( val ),
             std::get< 2 >( val ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        NodeManager,
                        std::string const &,
                        Group * const )

}  // namespace geosx
