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
 * @file NeighborCommunicator.cpp
 *
 */

#include "NeighborCommunicator.hpp"
#include "MPI_iCommData.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "mesh/MeshLevel.hpp"

namespace geosx
{

using namespace dataRepository;

NeighborCommunicator::NeighborCommunicator( int rank ):
  m_neighborRank( rank ),
  m_sendBufferSize(),
  m_receiveBufferSize(),
  m_sendBuffer{ maxComm },
  m_receiveBuffer{ maxComm }
{ }

void NeighborCommunicator::mpiISendReceive( buffer_unit_type const * const sendBuffer,
                                            int const sendSize,
                                            MPI_Request & sendRequest,
                                            buffer_unit_type * const receiveBuffer,
                                            int const receiveSize,
                                            MPI_Request & receiveRequest,
                                            int const commID,
                                            MPI_Comm mpiComm )
{
  int const sendTag = CommTag( MpiWrapper::commRank(), m_neighborRank, commID );
  //m_rank * m_size + m_neighborRank + m_size*m_size*commID;
  MpiWrapper::iSend( const_cast< buffer_unit_type * >(sendBuffer),
                     sendSize,
                     m_neighborRank,
                     sendTag,
                     mpiComm,
                     &sendRequest );

  int const receiveTag = CommTag( m_neighborRank, MpiWrapper::commRank(), commID );
  //m_neighborRank * m_size + m_rank + m_size*m_size*commID;
  MpiWrapper::iRecv( receiveBuffer,
                     receiveSize,
                     m_neighborRank,
                     receiveTag,
                     mpiComm,
                     &receiveRequest );
}


void NeighborCommunicator::mpiISendReceiveBufferSizes( int const commID,
                                                       MPI_Request & mpiSendRequest,
                                                       MPI_Request & mpiRecvRequest,
                                                       MPI_Comm mpiComm )
{
//  m_sendBufferSize[commID] = LvArray::integerConversion<int>( m_sendBuffer[commID].size());
  mpiISendReceive( &m_sendBufferSize[commID], 1, mpiSendRequest,
                   &m_receiveBufferSize[commID],
                   1, mpiRecvRequest,
                   commID,
                   mpiComm );
}

void NeighborCommunicator::mpiISendReceiveBuffers( int const commID,
                                                   MPI_Request & mpiSendRequest,
                                                   MPI_Request & mpiRecvRequest,
                                                   MPI_Comm mpiComm )
{
  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  mpiISendReceive( m_sendBuffer[commID].data(),
                   LvArray::integerConversion< int >( m_sendBuffer[commID].size()),
                   mpiSendRequest,
                   m_receiveBuffer[commID].data(),
                   LvArray::integerConversion< int >( m_receiveBuffer[commID].size()),
                   mpiRecvRequest,
                   commID,
                   mpiComm );

}


void NeighborCommunicator::mpiWaitAll( int const GEOSX_UNUSED_PARAM( commID ),
                                       MPI_Request & mpiSendRequest,
                                       MPI_Status & mpiSendStatus,
                                       MPI_Request & mpiRecvRequest,
                                       MPI_Status & mpiReceiveStatus )

{
  MpiWrapper::waitAll( 1, &mpiRecvRequest, &mpiReceiveStatus );
  MpiWrapper::waitAll( 1, &mpiSendRequest, &mpiSendStatus );
}

void NeighborCommunicator::clear()
{
  for( int i = 0; i < maxComm; ++i )
  {
    m_sendBuffer[i].clear();
    m_receiveBuffer[i].clear();
  }
}

void NeighborCommunicator::addNeighborGroupToMesh( MeshLevel & mesh ) const
{
  mesh.getNodeManager().addNeighbor( m_neighborRank );
  mesh.getEdgeManager().addNeighbor( m_neighborRank );
  mesh.getFaceManager().addNeighbor( m_neighborRank );

  mesh.getElemManager().forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & elementSubRegion )
  {
    elementSubRegion.addNeighbor( m_neighborRank );
  } );
}

int NeighborCommunicator::postSizeRecv( int const commID,
                                        MPI_Request & mpiRecvSizeRequest )
{
  int const recvTag = 101; //CommTag( m_neighborRank, MpiWrapper::commRank(), commID );
  return MpiWrapper::iRecv( &m_receiveBufferSize[commID],
                            1,
                            m_neighborRank,
                            recvTag,
                            MPI_COMM_GEOSX,
                            &mpiRecvSizeRequest );
}

int NeighborCommunicator::postSizeSend( int const commID,
                                        MPI_Request & mpiSendSizeRequest )
{
  int const sendTag = 101; //CommTag( m_neighborRank, MpiWrapper::commRank(), commID );
  return MpiWrapper::iSend( &m_sendBufferSize[commID],
                            1,
                            m_neighborRank,
                            sendTag,
                            MPI_COMM_GEOSX,
                            &mpiSendSizeRequest );
}

int NeighborCommunicator::postRecv( int const commID,
                                    MPI_Request & mpRecvRequest )
{
  int const recvTag = 102; //CommTag( m_neighborRank, MpiWrapper::commRank(), commID );
  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );
  return MpiWrapper::iRecv( m_receiveBuffer[commID].data(),
                            m_receiveBufferSize[commID],
                            m_neighborRank,
                            recvTag,
                            MPI_COMM_GEOSX,
                            &mpRecvRequest );
}

int NeighborCommunicator::postSend( int const commID,
                                    MPI_Request & mpiSendRequest )
{
  int const sendTag = 102; //CommTag( m_neighborRank, MpiWrapper::commRank(), commID );
  return MpiWrapper::iSend( m_sendBuffer[commID].data(),
                            m_sendBufferSize[commID],
                            m_neighborRank,
                            sendTag,
                            MPI_COMM_GEOSX,
                            &mpiSendRequest );
}

using ElemAdjListViewType = ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > >;
using ElemAdjListRefWrapType = ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > >;
using ElemAdjListRefType = ElementRegionManager::ElementReferenceAccessor< localIndex_array >;

inline int GhostSize( NodeManager & nodeManager, arrayView1d< localIndex const > const nodeAdjacencyList,
                      EdgeManager & edgeManager, arrayView1d< localIndex const > const edgeAdjacencyList,
                      FaceManager & faceManager, arrayView1d< localIndex const > const faceAdjacencyList,
                      ElementRegionManager & elemManager, ElemAdjListViewType const & elementAdjacencyList )
{
  int bufferSize = 0;
  bufferSize += nodeManager.packGlobalMapsSize( nodeAdjacencyList, 0 );
  bufferSize += edgeManager.packGlobalMapsSize( edgeAdjacencyList, 0 );
  bufferSize += faceManager.packGlobalMapsSize( faceAdjacencyList, 0 );
  bufferSize += elemManager.PackGlobalMapsSize( elementAdjacencyList );
  bufferSize += nodeManager.packUpDownMapsSize( nodeAdjacencyList );
  bufferSize += edgeManager.packUpDownMapsSize( edgeAdjacencyList );
  bufferSize += faceManager.packUpDownMapsSize( faceAdjacencyList );
  bufferSize += elemManager.PackUpDownMapsSize( elementAdjacencyList );
  parallelDeviceEvents events;
  bufferSize += nodeManager.packSize( {}, nodeAdjacencyList, 0, false, events );
  bufferSize += edgeManager.packSize( {}, edgeAdjacencyList, 0, false, events );
  bufferSize += faceManager.packSize( {}, faceAdjacencyList, 0, false, events );
  bufferSize += elemManager.PackSize( {}, elementAdjacencyList );
  waitAllDeviceEvents( events );
  return bufferSize;
}

inline int PackGhosts( buffer_unit_type * sendBufferPtr,
                       NodeManager & nodeManager, arrayView1d< localIndex const > const nodeAdjacencyList,
                       EdgeManager & edgeManager, arrayView1d< localIndex const > const edgeAdjacencyList,
                       FaceManager & faceManager, arrayView1d< localIndex const > const faceAdjacencyList,
                       ElementRegionManager & elemManager, ElemAdjListViewType const & elementAdjacencyList )
{
  int packedSize = 0;
  packedSize += nodeManager.packGlobalMaps( sendBufferPtr, nodeAdjacencyList, 0 );
  packedSize += edgeManager.packGlobalMaps( sendBufferPtr, edgeAdjacencyList, 0 );
  packedSize += faceManager.packGlobalMaps( sendBufferPtr, faceAdjacencyList, 0 );
  packedSize += elemManager.PackGlobalMaps( sendBufferPtr, elementAdjacencyList );
  packedSize += nodeManager.packUpDownMaps( sendBufferPtr, nodeAdjacencyList );
  packedSize += edgeManager.packUpDownMaps( sendBufferPtr, edgeAdjacencyList );
  packedSize += faceManager.packUpDownMaps( sendBufferPtr, faceAdjacencyList );
  packedSize += elemManager.PackUpDownMaps( sendBufferPtr, elementAdjacencyList );
  parallelDeviceEvents events;
  packedSize += nodeManager.pack( sendBufferPtr, {}, nodeAdjacencyList, 0, false, events );
  packedSize += edgeManager.pack( sendBufferPtr, {}, edgeAdjacencyList, 0, false, events );
  packedSize += faceManager.pack( sendBufferPtr, {}, faceAdjacencyList, 0, false, events );
  packedSize += elemManager.Pack( sendBufferPtr, {}, elementAdjacencyList );
  waitAllDeviceEvents( events );
  return packedSize;
}

void NeighborCommunicator::prepareAndSendGhosts( bool const GEOSX_UNUSED_PARAM( contactActive ),
                                                 integer const depth,
                                                 MeshLevel & mesh,
                                                 int const commID,
                                                 MPI_Request & mpiRecvSizeRequest,
                                                 MPI_Request & mpiSendSizeRequest,
                                                 MPI_Request & mpiSendRequest )
{
  GEOSX_MARK_FUNCTION;

  this->postSizeRecv( commID,
                      mpiRecvSizeRequest ); // post recv for buffer size from neighbor.

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  localIndex_array & nodeAdjacencyList = nodeManager.getNeighborData( m_neighborRank ).adjacencyList();
  localIndex_array & edgeAdjacencyList = edgeManager.getNeighborData( m_neighborRank ).adjacencyList();
  localIndex_array & faceAdjacencyList = faceManager.getNeighborData( m_neighborRank ).adjacencyList();

  {
    ElemAdjListRefWrapType elementAdjacencyList =
      elemManager.constructReferenceAccessor< localIndex_array >( ObjectManagerBase::viewKeyStruct::adjacencyListString(),
                                                                  std::to_string( this->m_neighborRank ) );

    mesh.generateAdjacencyLists( nodeManager.getNeighborData( m_neighborRank ).matchedPartitionBoundary(),
                                 nodeAdjacencyList,
                                 edgeAdjacencyList,
                                 faceAdjacencyList,
                                 elementAdjacencyList,
                                 depth );
  }

  ElemAdjListViewType const elemAdjacencyList =
    elemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( ObjectManagerBase::viewKeyStruct::adjacencyListString(),
                                                                                           std::to_string( this->m_neighborRank ) );

  int const bufferSize = GhostSize( nodeManager, nodeAdjacencyList,
                                    edgeManager, edgeAdjacencyList,
                                    faceManager, faceAdjacencyList,
                                    elemManager, elemAdjacencyList );

  this->resizeSendBuffer( commID, bufferSize );
  this->postSizeSend( commID,
                      mpiSendSizeRequest );

  buffer_type & sendBuff = sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuff.data();

  int const packedSize = PackGhosts( sendBufferPtr,
                                     nodeManager, nodeAdjacencyList,
                                     edgeManager, edgeAdjacencyList,
                                     faceManager, faceAdjacencyList,
                                     elemManager, elemAdjacencyList );

  GEOSX_ERROR_IF_NE( bufferSize, packedSize );

  this->postSend( commID, mpiSendRequest );
}

void NeighborCommunicator::unpackGhosts( MeshLevel & mesh,
                                         int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  buffer_type const & receiveBuff = receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuff.data();

  int unpackedSize = 0;

  localIndex_array nodeUnpackList;
  unpackedSize += nodeManager.unpackGlobalMaps( receiveBufferPtr, nodeUnpackList, 0 );

  localIndex_array edgeUnpackList;
  unpackedSize += edgeManager.unpackGlobalMaps( receiveBufferPtr, edgeUnpackList, 0 );

  localIndex_array faceUnpackList;
  unpackedSize += faceManager.unpackGlobalMaps( receiveBufferPtr, faceUnpackList, 0 );

  ElemAdjListRefType elementAdjacencyReceiveListArray =
    elemManager.constructReferenceAccessor< localIndex_array >( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString(),
                                                                std::to_string( this->m_neighborRank ) );
  unpackedSize += elemManager.UnpackGlobalMaps( receiveBufferPtr,
                                                elementAdjacencyReceiveListArray );

  ElemAdjListViewType elementAdjacencyReceiveList =
    elemManager.constructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >( ObjectManagerBase::viewKeyStruct::ghostsToReceiveString(),
                                                                                           std::to_string( this->m_neighborRank ) );

  unpackedSize += nodeManager.unpackUpDownMaps( receiveBufferPtr, nodeUnpackList, false, false );
  unpackedSize += edgeManager.unpackUpDownMaps( receiveBufferPtr, edgeUnpackList, false, false );
  unpackedSize += faceManager.unpackUpDownMaps( receiveBufferPtr, faceUnpackList, false, false );
  unpackedSize += elemManager.UnpackUpDownMaps( receiveBufferPtr, elementAdjacencyReceiveListArray, false );

  parallelDeviceEvents events;
  unpackedSize += nodeManager.unpack( receiveBufferPtr, nodeUnpackList, 0, false, events );
  unpackedSize += edgeManager.unpack( receiveBufferPtr, edgeUnpackList, 0, false, events );
  unpackedSize += faceManager.unpack( receiveBufferPtr, faceUnpackList, 0, false, events );
  unpackedSize += elemManager.Unpack( receiveBufferPtr, elementAdjacencyReceiveList );
  waitAllDeviceEvents( events );
}

void NeighborCommunicator::prepareAndSendSyncLists( MeshLevel const & mesh,
                                                    int const commID,
                                                    MPI_Request & mpiRecvSizeRequest,
                                                    MPI_Request & mpiSendSizeRequest,
                                                    MPI_Request & mpiSendRequest
                                                    )
{
  GEOSX_MARK_FUNCTION;

  this->postSizeRecv( commID,
                      mpiRecvSizeRequest );

  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const nodeGhostsToReceive = nodeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  arrayView1d< localIndex const > const edgeGhostsToReceive = edgeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  arrayView1d< localIndex const > const faceGhostsToReceive = faceManager.getNeighborData( m_neighborRank ).ghostsToReceive();

  arrayView1d< globalIndex const > const nodeLocalToGlobal = nodeManager.localToGlobalMap();
  arrayView1d< globalIndex const > const edgeLocalToGlobal = edgeManager.localToGlobalMap();
  arrayView1d< globalIndex const > const faceLocalToGlobal = faceManager.localToGlobalMap();

  buffer_type & sendBuff = sendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuff.data();

  int bufferSize = 0;
  bufferSize += bufferOps::Pack< false >( sendBufferPtr,
                                          nodeGhostsToReceive.toSliceConst(),
                                          nullptr,
                                          nodeGhostsToReceive.size(),
                                          nodeLocalToGlobal.toSliceConst() );

  bufferSize += bufferOps::Pack< false >( sendBufferPtr,
                                          edgeGhostsToReceive.toSliceConst(),
                                          nullptr,
                                          edgeGhostsToReceive.size(),
                                          edgeLocalToGlobal.toSliceConst() );

  bufferSize += bufferOps::Pack< false >( sendBufferPtr,
                                          faceGhostsToReceive.toSliceConst(),
                                          nullptr,
                                          faceGhostsToReceive.size(),
                                          faceLocalToGlobal.toSliceConst() );

  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
  {
    arrayView1d< localIndex const > const ghostsToReceive = subRegion.getNeighborData( m_neighborRank ).ghostsToReceive();
    arrayView1d< globalIndex const > const localToGlobal = subRegion.localToGlobalMap();
    bufferSize += bufferOps::Pack< false >( sendBufferPtr,
                                            ghostsToReceive.toSliceConst(),
                                            nullptr,
                                            ghostsToReceive.size(),
                                            localToGlobal.toSliceConst() );
  } );

  this->resizeSendBuffer( commID, bufferSize );
  this->postSizeSend( commID, mpiSendSizeRequest );

  int packedSize = 0;
  packedSize += bufferOps::Pack< true >( sendBufferPtr,
                                         nodeGhostsToReceive.toSliceConst(),
                                         nullptr,
                                         nodeGhostsToReceive.size(),
                                         nodeLocalToGlobal.toSliceConst() );

  packedSize += bufferOps::Pack< true >( sendBufferPtr,
                                         edgeGhostsToReceive.toSliceConst(),
                                         nullptr,
                                         edgeGhostsToReceive.size(),
                                         edgeLocalToGlobal.toSliceConst() );

  packedSize += bufferOps::Pack< true >( sendBufferPtr,
                                         faceGhostsToReceive.toSliceConst(),
                                         nullptr,
                                         faceGhostsToReceive.size(),
                                         faceLocalToGlobal.toSliceConst() );

  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
  {
    arrayView1d< localIndex const > const ghostsToReceive = subRegion.getNeighborData( m_neighborRank ).ghostsToReceive();
    arrayView1d< globalIndex const > const localToGlobal = subRegion.localToGlobalMap();
    packedSize += bufferOps::Pack< true >( sendBufferPtr,
                                           ghostsToReceive.toSliceConst(),
                                           nullptr,
                                           ghostsToReceive.size(),
                                           localToGlobal.toSliceConst() );
  } );

  GEOSX_ERROR_IF( bufferSize != packedSize, "Allocated Buffer Size is not equal to packed buffer size" );

  this->postSend( commID, mpiSendRequest );
}

void NeighborCommunicator::unpackAndRebuildSyncLists( MeshLevel & mesh,
                                                      int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  localIndex_array & nodeGhostsToSend = nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  localIndex_array & edgeGhostsToSend = edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  localIndex_array & faceGhostsToSend = faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  buffer_type const & receiveBuff = receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuff.data();

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             nodeGhostsToSend,
                             nodeManager.globalToLocalMap() );

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             edgeGhostsToSend,
                             edgeManager.globalToLocalMap() );

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             faceGhostsToSend,
                             faceManager.globalToLocalMap() );

  nodeManager.setGhostRankForSenders( m_neighborRank );
  edgeManager.setGhostRankForSenders( m_neighborRank );
  faceManager.setGhostRankForSenders( m_neighborRank );

  elemManager.forElementSubRegions< ElementSubRegionBase >( [&] ( ElementSubRegionBase & subRegion )
  {
    bufferOps::UnpackSyncList( receiveBufferPtr,
                               subRegion.getNeighborData( m_neighborRank ).ghostsToSend(),
                               subRegion.globalToLocalMap() );

    subRegion.setGhostRankForSenders( m_neighborRank );
  } );
}


int NeighborCommunicator::packCommSizeForSync( string_array const & fieldNames,
                                               ObjectManagerBase const & manager,
                                               int const commID,
                                               bool onDevice,
                                               parallelDeviceEvents & events )
{
  arrayView1d< localIndex const > const & ghostsToSend = manager.getNeighborData( m_neighborRank ).ghostsToSend();
  localIndex const bufferSize = manager.packSize( fieldNames, ghostsToSend, 0, onDevice, events );
  m_sendBufferSize[commID] = LvArray::integerConversion< int >( bufferSize );
  return m_sendBufferSize[commID];
}


int NeighborCommunicator::packCommSizeForSync( std::map< string, string_array > const & fieldNames,
                                               MeshLevel const & mesh,
                                               int const commID,
                                               bool onDevice,
                                               parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const & nodeGhostsToSend = nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & edgeGhostsToSend = edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & faceGhostsToSend = faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  int bufferSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    bufferSize += nodeManager.packSize( fieldNames.at( "node" ), nodeGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    bufferSize += edgeManager.packSize( fieldNames.at( "edge" ), edgeGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    bufferSize += faceManager.packSize( fieldNames.at( "face" ), faceGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
    {
      bufferSize += subRegion.packSize( fieldNames.at( "elems" ), subRegion.getNeighborData( m_neighborRank ).ghostsToSend(), 0, onDevice, events );
    } );
  }

  this->m_sendBufferSize[commID] = bufferSize;
  return bufferSize;
}


void NeighborCommunicator::packCommBufferForSync( string_array const & fieldNames,
                                                  ObjectManagerBase const & manager,
                                                  int const commID,
                                                  bool onDevice,
                                                  parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< localIndex const > const & ghostsToSend = manager.getNeighborData( m_neighborRank ).ghostsToSend();

  buffer_type & sendBuff = sendBuffer( commID );
  localIndex const bufferSize =  LvArray::integerConversion< localIndex >( sendBuff.size() );
  buffer_unit_type * sendBufferPtr = sendBuff.data();

  localIndex const packedSize = manager.pack( sendBufferPtr, fieldNames, ghostsToSend, 0, onDevice, events );
  GEOSX_ERROR_IF_NE( bufferSize, packedSize );
}


void NeighborCommunicator::packCommBufferForSync( std::map< string, string_array > const & fieldNames,
                                                  MeshLevel const & mesh,
                                                  int const commID,
                                                  bool onDevice,
                                                  parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = mesh.getNodeManager();
  EdgeManager const & edgeManager = mesh.getEdgeManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView1d< localIndex const > const & nodeGhostsToSend = nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & edgeGhostsToSend = edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & faceGhostsToSend = faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  buffer_type & sendBuff = sendBuffer( commID );
  localIndex const bufferSize =  LvArray::integerConversion< localIndex >( sendBuff.size() );
  buffer_unit_type * sendBufferPtr = sendBuff.data();

  localIndex packedSize = 0;
  if( fieldNames.count( "node" ) > 0 )
  {
    packedSize += nodeManager.pack( sendBufferPtr, fieldNames.at( "node" ), nodeGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    packedSize += edgeManager.pack( sendBufferPtr, fieldNames.at( "edge" ), edgeGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    packedSize += faceManager.pack( sendBufferPtr, fieldNames.at( "face" ), faceGhostsToSend, 0, onDevice, events );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion )
    {
      packedSize += subRegion.pack( sendBufferPtr, fieldNames.at( "elems" ), subRegion.getNeighborData( m_neighborRank ).ghostsToSend(), 0, onDevice, events );
    } );
  }

  GEOSX_ERROR_IF_NE( bufferSize, packedSize );
}


void NeighborCommunicator::unpackBufferForSync( ObjectManagerBase & manager,
                                                int const commID,
                                                bool onDevice,
                                                parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  buffer_type const & receiveBuff = receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuff.data();

  array1d< localIndex > & ghostsToReceive = manager.getNeighborData( m_neighborRank ).ghostsToReceive();
  manager.unpack( receiveBufferPtr, ghostsToReceive, 0, onDevice, events );
}


void NeighborCommunicator::unpackBufferForSync( std::map< string, string_array > const & fieldNames,
                                                MeshLevel & mesh,
                                                int const commID,
                                                bool onDevice,
                                                parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  buffer_type const & receiveBuff = receiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuff.data();

  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  array1d< localIndex > & nodeGhostsToReceive = nodeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  array1d< localIndex > & edgeGhostsToReceive = edgeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  array1d< localIndex > & faceGhostsToReceive = faceManager.getNeighborData( m_neighborRank ).ghostsToReceive();

  int unpackedSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    unpackedSize += nodeManager.unpack( receiveBufferPtr, nodeGhostsToReceive, 0, onDevice, events );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    unpackedSize += edgeManager.unpack( receiveBufferPtr, edgeGhostsToReceive, 0, onDevice, events );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    unpackedSize += faceManager.unpack( receiveBufferPtr, faceGhostsToReceive, 0, onDevice, events );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >( [&] ( ElementSubRegionBase & subRegion )
    {
      unpackedSize += subRegion.unpack( receiveBufferPtr, subRegion.getNeighborData( m_neighborRank ).ghostsToReceive(), 0, onDevice, events );
    } );
  }
}


} /* namespace geosx */
