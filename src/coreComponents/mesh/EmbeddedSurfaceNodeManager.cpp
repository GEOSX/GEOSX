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
 * @file EmbeddedSurfaceNodeManager.cpp
 */

#include "EmbeddedSurfaceNodeManager.hpp"
#include "FaceManager.hpp"
#include "EdgeManager.hpp"
#include "ToElementRelation.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"

namespace geosx
{

using namespace dataRepository;

EmbeddedSurfaceNodeManager::EmbeddedSurfaceNodeManager( string const & name,
                                                        Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_referencePosition( 0, 3 )
{
  registerWrapper( viewKeyStruct::referencePositionString(), &m_referencePosition );
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_toEdgesRelation );
  this->registerWrapper( viewKeyStruct::elementRegionListString(), &elementRegionList() );
  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &elementSubRegionList() );
  this->registerWrapper( viewKeyStruct::elementListString(), &elementList() );
  this->registerWrapper( viewKeyStruct::parentEdgeGlobalIndexString(), &m_parentEdgeGlobalIndex );
}


EmbeddedSurfaceNodeManager::~EmbeddedSurfaceNodeManager()
{}


void EmbeddedSurfaceNodeManager::resize( localIndex const newSize )
{
  m_toEdgesRelation.resize( newSize, 2 * getEdgeMapOverallocation() );
  m_toElements.m_toElementRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementSubRegion.resize( newSize, 2 * getElemMapOverAllocation() );
  m_toElements.m_toElementIndex.resize( newSize, 2 * getElemMapOverAllocation() );
  ObjectManagerBase::resize( newSize );
}


void EmbeddedSurfaceNodeManager::setEdgeMaps( EdgeManager const & embSurfEdgeManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const edgeToNodeMap = embSurfEdgeManager.nodeList();
  localIndex const numEdges = edgeToNodeMap.size( 0 );
  localIndex const numNodes = size();

  ArrayOfArrays< localIndex > toEdgesTemp( numNodes, EdgeManager::maxEdgesPerNode() );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeEdges = 0;

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID )
  {
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 0 ), edgeID );
    toEdgesTemp.emplaceBackAtomic< parallelHostAtomic >( edgeToNodeMap( edgeID, 1 ), edgeID );
    totalNodeEdges += 2;
  } );

  // Resize the node to edge map.
  m_toEdgesRelation.resize( 0 );

  // Reserve space for the number of current nodes plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numNodes;
  m_toEdgesRelation.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalNodeEdges.get() + numNodes * getEdgeMapOverallocation() * ( 1 + 2 * overAllocationFactor );
  m_toEdgesRelation.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    m_toEdgesRelation.appendSet( toEdgesTemp.sizeOfArray( nodeID ) + getEdgeMapOverallocation() );
  }

  ArrayOfSetsView< localIndex > const & toEdgesView = m_toEdgesRelation.toView();
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex * const edges = toEdgesTemp[ nodeID ];
    localIndex const numNodeEdges = toEdgesTemp.sizeOfArray( nodeID );
    localIndex const numUniqueEdges = LvArray::sortedArrayManipulation::makeSortedUnique( edges, edges + numNodeEdges );
    toEdgesView.insertIntoSet( nodeID, edges, edges + numUniqueEdges );
  } );

  m_toEdgesRelation.setRelatedObject( embSurfEdgeManager );
}

void EmbeddedSurfaceNodeManager::setElementMaps( ElementRegionManager const & elementRegionManager )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< localIndex > & toElementRegionList = m_toElements.m_toElementRegion;
  ArrayOfArrays< localIndex > & toElementSubRegionList = m_toElements.m_toElementSubRegion;
  ArrayOfArrays< localIndex > & toElementList = m_toElements.m_toElementIndex;
  localIndex const numNodes = size();

  // The number of elements attached to the each node.
  array1d< localIndex > elemsPerNode( numNodes );

  // The total number of elements, the sum of elemsPerNode.
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalNodeElems = 0;

  elementRegionManager.
    forElementSubRegions< EmbeddedSurfaceSubRegion >( [&elemsPerNode, &totalNodeElems]( EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    forAll< parallelHostPolicy >( subRegion.size(), [&elemsPerNode, totalNodeElems, &elemToNodeMap ] ( localIndex const k )
    {
      localIndex const numElementNodes = elemToNodeMap.sizeOfArray( k );
      totalNodeElems += numElementNodes;
      for( localIndex a = 0; a < numElementNodes; ++a )
      {
        localIndex const nodeIndex = elemToNodeMap( k, a );
        RAJA::atomicInc< parallelHostAtomic >( &elemsPerNode[ nodeIndex ] );
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
  localIndex const valuesToReserve = totalNodeElems.get() + numNodes * getElemMapOverAllocation() * ( 1 + 2 * overAllocationFactor );
  toElementRegionList.reserveValues( valuesToReserve );
  toElementSubRegionList.reserveValues( valuesToReserve );
  toElementList.reserveValues( valuesToReserve );

  // Append an array for each node with capacity to hold the appropriate number of elements plus some wiggle room.
  for( localIndex nodeID = 0; nodeID < numNodes; ++nodeID )
  {
    toElementRegionList.appendArray( 0 );
    toElementSubRegionList.appendArray( 0 );
    toElementList.appendArray( 0 );

    toElementRegionList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
    toElementSubRegionList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
    toElementList.setCapacityOfArray( nodeID, elemsPerNode[ nodeID ] + getElemMapOverAllocation() );
  }

  // Populate the element maps.
  // Note that this can't be done in parallel because the three element lists must be in the same order.
  // If this becomes a bottleneck create a temporary ArrayOfArrays of tuples and insert into that first then copy over.
  elementRegionManager.
    forElementSubRegionsComplete< EmbeddedSurfaceSubRegion >( [&toElementRegionList, &toElementSubRegionList, &toElementList]
                                                                ( localIndex const er, localIndex const esr, ElementRegionBase const &,
                                                                EmbeddedSurfaceSubRegion const & subRegion )
  {
    EmbeddedSurfaceSubRegion::NodeMapType const & elemToNodeMap = subRegion.nodeList();
    for( localIndex k = 0; k < subRegion.size(); ++k )
    {
      for( localIndex a=0; a<elemToNodeMap.sizeOfArray( k ); ++a )
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


void EmbeddedSurfaceNodeManager::compressRelationMaps()
{
  m_toEdgesRelation.compress();
  m_toElements.m_toElementRegion.compress();
  m_toElements.m_toElementSubRegion.compress();
  m_toElements.m_toElementIndex.compress();
}


void EmbeddedSurfaceNodeManager::appendNode( arraySlice1d< real64 const > const & pointCoord,
                                             integer const & pointGhostRank )
{
  if( pointGhostRank < 0 )
  {
    localIndex nodeIndex =  this->size();
    this->resize( nodeIndex + 1 );
    LvArray::tensorOps::copy< 3 >( m_referencePosition[nodeIndex], pointCoord );
    m_ghostRank[ nodeIndex ] = pointGhostRank;
  }
}

void EmbeddedSurfaceNodeManager::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementSubRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString() ));
}

localIndex EmbeddedSurfaceNodeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}


localIndex EmbeddedSurfaceNodeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                                       arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}


template< bool DOPACK >
localIndex EmbeddedSurfaceNodeManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                                              arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );
  return packedSize;
}


localIndex EmbeddedSurfaceNodeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                                         localIndex_array & packList,
                                                         bool const overwriteUpMaps,
                                                         bool const )
{
  localIndex unPackedSize = 0;

  string temp;

  unPackedSize += bufferOps::Unpack( buffer, temp );
  GEOSX_ERROR_IF( temp != viewKeyStruct::elementListString(), "" );
  unPackedSize += bufferOps::Unpack( buffer,
                                     this->m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}


REGISTER_CATALOG_ENTRY( ObjectManagerBase, EmbeddedSurfaceNodeManager, string const &, Group * const )

}
