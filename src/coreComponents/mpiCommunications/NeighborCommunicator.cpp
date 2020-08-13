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
 * @file NeighborCommunicator.cpp
 *
 */

#include "mpiCommunications/NeighborCommunicator.hpp"

#include "common/TimingMacros.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "mesh/MeshLevel.hpp"
#include <sys/time.h>

namespace geosx
{
using namespace dataRepository;

NeighborCommunicator::NeighborCommunicator() :
  m_neighborRank( -1 ),
  m_sendBufferSize(),
  m_receiveBufferSize(),
  m_sendBuffer { maxComm },
  m_receiveBuffer { maxComm },
  m_mpiSendBufferRequest(),
  m_mpiRecvBufferRequest(),
  m_mpiSendBufferStatus(),
  m_mpiRecvBufferStatus()
{}

void
NeighborCommunicator::MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                                        int const sendSize,
                                        MPI_Request & sendRequest,
                                        buffer_unit_type * const receiveBuffer,
                                        int const receiveSize,
                                        MPI_Request & receiveRequest,
                                        int const commID,
                                        MPI_Comm mpiComm )
{
  int const sendTag = CommTag( MpiWrapper::Comm_rank(), m_neighborRank, commID );
  //m_rank * m_size + m_neighborRank + m_size*m_size*commID;
  MpiWrapper::iSend( const_cast< buffer_unit_type * >( sendBuffer ),
                     sendSize,
                     m_neighborRank,
                     sendTag,
                     mpiComm,
                     &sendRequest );

  int const receiveTag = CommTag( m_neighborRank, MpiWrapper::Comm_rank(), commID );
  //m_neighborRank * m_size + m_rank + m_size*m_size*commID;
  MpiWrapper::iRecv( receiveBuffer,
                     receiveSize,
                     m_neighborRank,
                     receiveTag,
                     mpiComm,
                     &receiveRequest );
}

void
NeighborCommunicator::MPI_iSendReceiveBufferSizes( int const commID,
                                                   MPI_Comm mpiComm )
{
  MPI_iSendReceiveBufferSizes( commID,
                               m_mpiSendBufferRequest[commID],
                               m_mpiRecvBufferRequest[commID],
                               mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceiveBufferSizes( int const commID,
                                                   MPI_Request & mpiSendRequest,
                                                   MPI_Request & mpiRecvRequest,
                                                   MPI_Comm mpiComm )
{
  //  m_sendBufferSize[commID] = LvArray::integerConversion<int>( m_sendBuffer[commID].size());
  MPI_iSendReceive( &m_sendBufferSize[commID],
                    1,
                    mpiSendRequest,
                    &m_receiveBufferSize[commID],
                    1,
                    mpiRecvRequest,
                    commID,
                    mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceiveBuffers( int const commID,
                                               MPI_Comm mpiComm )
{
  MPI_iSendReceiveBuffers( commID,
                           m_mpiSendBufferRequest[commID],
                           m_mpiRecvBufferRequest[commID],
                           mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceiveBuffers( int const commID,
                                               MPI_Request & mpiSendRequest,
                                               MPI_Request & mpiRecvRequest,
                                               MPI_Comm mpiComm )
{
  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive(
    m_sendBuffer[commID].data(),
    LvArray::integerConversion< int >( m_sendBuffer[commID].size() ),
    mpiSendRequest,
    m_receiveBuffer[commID].data(),
    LvArray::integerConversion< int >( m_receiveBuffer[commID].size() ),
    mpiRecvRequest,
    commID,
    mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceive( int const commID, MPI_Comm mpiComm )
{
  MPI_iSendReceive( commID,
                    m_mpiSendBufferRequest[commID],
                    m_mpiRecvBufferRequest[commID],
                    mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceive( int const commID,
                                        MPI_Request & mpiSendRequest,
                                        MPI_Request & mpiRecvRequest,
                                        MPI_Comm mpiComm )
{
  MPI_iSendReceiveBufferSizes( commID, mpiComm );

  MpiWrapper::Waitall( 1,
                       &( m_mpiRecvBufferRequest[commID] ),
                       &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1,
                       &( m_mpiSendBufferRequest[commID] ),
                       &( m_mpiSendBufferStatus[commID] ) );

  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( m_sendBuffer[commID].data(),
                    m_sendBufferSize[commID],
                    mpiSendRequest,
                    m_receiveBuffer[commID].data(),
                    m_receiveBufferSize[commID],
                    mpiRecvRequest,
                    commID,
                    mpiComm );
}

void
NeighborCommunicator::MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                                        int const sendSize,
                                        int const commID,
                                        MPI_Comm mpiComm )
{
  MPI_iSendReceive( &sendSize,
                    1,
                    m_mpiSendBufferRequest[commID],
                    &m_receiveBufferSize[commID],
                    1,
                    m_mpiRecvBufferRequest[commID],
                    commID,
                    mpiComm );

  MpiWrapper::Waitall( 1,
                       &( m_mpiRecvBufferRequest[commID] ),
                       &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1,
                       &( m_mpiSendBufferRequest[commID] ),
                       &( m_mpiSendBufferStatus[commID] ) );

  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( sendBuffer,
                    sendSize,
                    m_mpiSendBufferRequest[commID],
                    m_receiveBuffer[commID].data(),
                    m_receiveBufferSize[commID],
                    m_mpiRecvBufferRequest[commID],
                    commID,
                    mpiComm );
}

void
NeighborCommunicator::MPI_WaitAll( int const GEOSX_UNUSED_PARAM( commID ),
                                   MPI_Request & mpiSendRequest,
                                   MPI_Status & mpiSendStatus,
                                   MPI_Request & mpiRecvRequest,
                                   MPI_Status & mpiReceiveStatus )

{
  MpiWrapper::Waitall( 1, &mpiRecvRequest, &mpiReceiveStatus );
  MpiWrapper::Waitall( 1, &mpiSendRequest, &mpiSendStatus );
}

void
NeighborCommunicator::MPI_WaitAll( int const commID )
{
  MpiWrapper::Waitall( 1,
                       &( m_mpiRecvBufferRequest[commID] ),
                       &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Waitall( 1,
                       &( m_mpiSendBufferRequest[commID] ),
                       &( m_mpiSendBufferStatus[commID] ) );
}

void
NeighborCommunicator::Clear()
{
  for( int i = 0; i < maxComm; ++i )
  {
    m_sendBuffer[i].clear();
    m_receiveBuffer[i].clear();
  }
}

void
NeighborCommunicator::AddNeighborGroupToMesh( MeshLevel & mesh ) const
{
  mesh.getNodeManager()->addNeighbor( m_neighborRank );
  mesh.getEdgeManager()->addNeighbor( m_neighborRank );
  mesh.getFaceManager()->addNeighbor( m_neighborRank );

  mesh.getElemManager()->forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase & elementSubRegion ) {
      elementSubRegion.addNeighbor( m_neighborRank );
    } );
}

int
NeighborCommunicator::PostSizeRecv( int const commID )
{
  int const recvTag =
    101;  //CommTag( m_neighborRank, MpiWrapper::Comm_rank(), commID );
  return MpiWrapper::iRecv( &m_receiveBufferSize[commID],
                            1,
                            m_neighborRank,
                            recvTag,
                            MPI_COMM_GEOSX,
                            &m_mpiRecvSizeRequest[commID] );
}

MPI_Request
NeighborCommunicator::GetSizeRecvRequest( int const commID )
{
  return m_mpiRecvSizeRequest[commID];
}

int
NeighborCommunicator::PostSizeSend( int const commID )
{
  int const sendTag =
    101;  //CommTag( m_neighborRank, MpiWrapper::Comm_rank(), commID );
  return MpiWrapper::iSend( &m_sendBufferSize[commID],
                            1,
                            m_neighborRank,
                            sendTag,
                            MPI_COMM_GEOSX,
                            &m_mpiSendSizeRequest[commID] );
}

int
NeighborCommunicator::PostRecv( int const commID )
{
  int const recvTag =
    102;  //CommTag( m_neighborRank, MpiWrapper::Comm_rank(), commID );
  m_receiveBuffer[commID].resize( m_receiveBufferSize[commID] );
  return MpiWrapper::iRecv( m_receiveBuffer[commID].data(),
                            m_receiveBufferSize[commID],
                            m_neighborRank,
                            recvTag,
                            MPI_COMM_GEOSX,
                            &m_mpiRecvBufferRequest[commID] );
}

MPI_Request
NeighborCommunicator::GetRecvRequest( int const commID )
{
  return m_mpiRecvBufferRequest[commID];
}

int
NeighborCommunicator::PostSend( int const commID )
{
  int const sendTag =
    102;  //CommTag( m_neighborRank, MpiWrapper::Comm_rank(), commID );
  return MpiWrapper::iSend( m_sendBuffer[commID].data(),
                            m_sendBufferSize[commID],
                            m_neighborRank,
                            sendTag,
                            MPI_COMM_GEOSX,
                            &m_mpiSendBufferRequest[commID] );
}

using ElemAdjListViewType =
  ElementRegionManager::ElementViewAccessor< arrayView1d< localIndex > >;
using ElemAdjListRefWrapType =
  ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > >;
using ElemAdjListRefType =
  ElementRegionManager::ElementReferenceAccessor< localIndex_array >;

inline int
GhostSize( NodeManager & nodeManager,
           localIndex_array & nodeAdjacencyList,
           EdgeManager & edgeManager,
           localIndex_array & edgeAdjacencyList,
           FaceManager & faceManager,
           localIndex_array & faceAdjacencyList,
           ElementRegionManager & elemManager,
           ElemAdjListViewType const & elementAdjacencyList )
{
  int bufferSize = 0;
  bufferSize += nodeManager.PackGlobalMapsSize( nodeAdjacencyList, 0 );
  bufferSize += edgeManager.PackGlobalMapsSize( edgeAdjacencyList, 0 );
  bufferSize += faceManager.PackGlobalMapsSize( faceAdjacencyList, 0 );
  bufferSize += elemManager.PackGlobalMapsSize( elementAdjacencyList );
  bufferSize += nodeManager.PackUpDownMapsSize( nodeAdjacencyList );
  bufferSize += edgeManager.PackUpDownMapsSize( edgeAdjacencyList );
  bufferSize += faceManager.PackUpDownMapsSize( faceAdjacencyList );
  bufferSize += elemManager.PackUpDownMapsSize( elementAdjacencyList );
  bufferSize += nodeManager.PackSize( {}, nodeAdjacencyList, 0 );
  bufferSize += edgeManager.PackSize( {}, edgeAdjacencyList, 0 );
  bufferSize += faceManager.PackSize( {}, faceAdjacencyList, 0 );
  bufferSize += elemManager.PackSize( {}, elementAdjacencyList );
  return bufferSize;
}

inline int
PackGhosts( buffer_unit_type * sendBufferPtr,
            NodeManager & nodeManager,
            localIndex_array & nodeAdjacencyList,
            EdgeManager & edgeManager,
            localIndex_array & edgeAdjacencyList,
            FaceManager & faceManager,
            localIndex_array & faceAdjacencyList,
            ElementRegionManager & elemManager,
            ElemAdjListViewType const & elementAdjacencyList )
{
  int packedSize = 0;
  packedSize += nodeManager.PackGlobalMaps( sendBufferPtr, nodeAdjacencyList, 0 );
  packedSize += edgeManager.PackGlobalMaps( sendBufferPtr, edgeAdjacencyList, 0 );
  packedSize += faceManager.PackGlobalMaps( sendBufferPtr, faceAdjacencyList, 0 );
  packedSize += elemManager.PackGlobalMaps( sendBufferPtr, elementAdjacencyList );
  packedSize += nodeManager.PackUpDownMaps( sendBufferPtr, nodeAdjacencyList );
  packedSize += edgeManager.PackUpDownMaps( sendBufferPtr, edgeAdjacencyList );
  packedSize += faceManager.PackUpDownMaps( sendBufferPtr, faceAdjacencyList );
  packedSize += elemManager.PackUpDownMaps( sendBufferPtr, elementAdjacencyList );
  packedSize += nodeManager.Pack( sendBufferPtr, {}, nodeAdjacencyList, 0 );
  packedSize += edgeManager.Pack( sendBufferPtr, {}, edgeAdjacencyList, 0 );
  packedSize += faceManager.Pack( sendBufferPtr, {}, faceAdjacencyList, 0 );
  packedSize += elemManager.Pack( sendBufferPtr, {}, elementAdjacencyList );
  return packedSize;
}

void
NeighborCommunicator::PrepareAndSendGhosts(
  bool const GEOSX_UNUSED_PARAM( contactActive ),
  integer const depth,
  MeshLevel & mesh,
  int const commID )
{
  GEOSX_MARK_FUNCTION;

  this->PostSizeRecv( commID );  // post recv for buffer size from neighbor.

  NodeManager & nodeManager = *( mesh.getNodeManager() );
  EdgeManager & edgeManager = *( mesh.getEdgeManager() );
  FaceManager & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager & elemManager = *( mesh.getElemManager() );

  localIndex_array & nodeAdjacencyList =
    nodeManager.getNeighborData( m_neighborRank ).adjacencyList();
  localIndex_array & edgeAdjacencyList =
    edgeManager.getNeighborData( m_neighborRank ).adjacencyList();
  localIndex_array & faceAdjacencyList =
    faceManager.getNeighborData( m_neighborRank ).adjacencyList();

  {
    ElemAdjListRefWrapType elementAdjacencyList =
      elemManager.ConstructReferenceAccessor< localIndex_array >(
        ObjectManagerBase::viewKeyStruct::adjacencyListString,
        std::to_string( this->m_neighborRank ) );

    mesh.GenerateAdjacencyLists(
      nodeManager.getNeighborData( m_neighborRank ).matchedPartitionBoundary(),
      nodeAdjacencyList,
      edgeAdjacencyList,
      faceAdjacencyList,
      elementAdjacencyList,
      depth );
  }

  ElemAdjListViewType const elemAdjacencyList =
    elemManager.ConstructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >(
      ObjectManagerBase::viewKeyStruct::adjacencyListString,
      std::to_string( this->m_neighborRank ) );

  int bufferSize = GhostSize( nodeManager,
                              nodeAdjacencyList,
                              edgeManager,
                              edgeAdjacencyList,
                              faceManager,
                              faceAdjacencyList,
                              elemManager,
                              elemAdjacencyList );

  this->resizeSendBuffer( commID, bufferSize );
  this->PostSizeSend( commID );

  buffer_type & sendBuffer = SendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = PackGhosts( sendBufferPtr,
                               nodeManager,
                               nodeAdjacencyList,
                               edgeManager,
                               edgeAdjacencyList,
                               faceManager,
                               faceAdjacencyList,
                               elemManager,
                               elemAdjacencyList );

  GEOSX_ERROR_IF_NE( bufferSize, packedSize );

  this->PostSend( commID );
}

void
NeighborCommunicator::UnpackGhosts( MeshLevel & mesh, int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = *( mesh.getNodeManager() );
  EdgeManager & edgeManager = *( mesh.getEdgeManager() );
  FaceManager & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager & elemManager = *( mesh.getElemManager() );

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  int unpackedSize = 0;

  localIndex_array nodeUnpackList;
  unpackedSize +=
    nodeManager.UnpackGlobalMaps( receiveBufferPtr, nodeUnpackList, 0 );

  localIndex_array edgeUnpackList;
  unpackedSize +=
    edgeManager.UnpackGlobalMaps( receiveBufferPtr, edgeUnpackList, 0 );

  localIndex_array faceUnpackList;
  unpackedSize +=
    faceManager.UnpackGlobalMaps( receiveBufferPtr, faceUnpackList, 0 );

  ElemAdjListRefType elementAdjacencyReceiveListArray =
    elemManager.ConstructReferenceAccessor< localIndex_array >(
      ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
      std::to_string( this->m_neighborRank ) );
  unpackedSize += elemManager.UnpackGlobalMaps( receiveBufferPtr,
                                                elementAdjacencyReceiveListArray );

  ElemAdjListViewType elementAdjacencyReceiveList =
    elemManager.ConstructViewAccessor< array1d< localIndex >, arrayView1d< localIndex > >(
      ObjectManagerBase::viewKeyStruct::ghostsToReceiveString,
      std::to_string( this->m_neighborRank ) );

  unpackedSize +=
    nodeManager.UnpackUpDownMaps( receiveBufferPtr, nodeUnpackList, false, false );
  unpackedSize +=
    edgeManager.UnpackUpDownMaps( receiveBufferPtr, edgeUnpackList, false, false );
  unpackedSize +=
    faceManager.UnpackUpDownMaps( receiveBufferPtr, faceUnpackList, false, false );
  unpackedSize += elemManager.UnpackUpDownMaps( receiveBufferPtr,
                                                elementAdjacencyReceiveListArray,
                                                false );

  unpackedSize += nodeManager.Unpack( receiveBufferPtr, nodeUnpackList, 0 );
  unpackedSize += edgeManager.Unpack( receiveBufferPtr, edgeUnpackList, 0 );
  unpackedSize += faceManager.Unpack( receiveBufferPtr, faceUnpackList, 0 );
  unpackedSize +=
    elemManager.Unpack( receiveBufferPtr, elementAdjacencyReceiveList );
}

void
NeighborCommunicator::PrepareAndSendSyncLists( MeshLevel const & mesh,
                                               int const commID )
{
  GEOSX_MARK_FUNCTION;

  this->PostSizeRecv( commID );

  NodeManager const & nodeManager = *( mesh.getNodeManager() );
  EdgeManager const & edgeManager = *( mesh.getEdgeManager() );
  FaceManager const & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager const & elemManager = *( mesh.getElemManager() );

  arraySlice1d< localIndex const > const nodeGhostsToReceive =
    nodeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  arraySlice1d< localIndex const > const edgeGhostsToReceive =
    edgeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  arraySlice1d< localIndex const > const faceGhostsToReceive =
    faceManager.getNeighborData( m_neighborRank ).ghostsToReceive();

  buffer_type & sendBuffer = SendBuffer( commID );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int bufferSize = 0;
  bufferSize +=
    bufferOps::Pack< false >( sendBufferPtr,
                              nodeGhostsToReceive,
                              nullptr,
                              nodeGhostsToReceive.size(),
                              nodeManager.localToGlobalMap().toSliceConst() );

  bufferSize +=
    bufferOps::Pack< false >( sendBufferPtr,
                              edgeGhostsToReceive,
                              nullptr,
                              edgeGhostsToReceive.size(),
                              edgeManager.localToGlobalMap().toSliceConst() );

  bufferSize +=
    bufferOps::Pack< false >( sendBufferPtr,
                              faceGhostsToReceive,
                              nullptr,
                              faceGhostsToReceive.size(),
                              faceManager.localToGlobalMap().toSliceConst() );

  elemManager.forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase const & subRegion ) {
      arraySlice1d< localIndex const > const ghostsToReceive =
        subRegion.getNeighborData( m_neighborRank ).ghostsToReceive();
      bufferSize +=
        bufferOps::Pack< false >( sendBufferPtr,
                                  ghostsToReceive,
                                  nullptr,
                                  ghostsToReceive.size(),
                                  subRegion.localToGlobalMap().toSliceConst() );
    } );

  this->resizeSendBuffer( commID, bufferSize );
  this->PostSizeSend( commID );

  int packedSize = 0;
  packedSize +=
    bufferOps::Pack< true >( sendBufferPtr,
                             nodeGhostsToReceive,
                             nullptr,
                             nodeGhostsToReceive.size(),
                             nodeManager.localToGlobalMap().toSliceConst() );

  packedSize +=
    bufferOps::Pack< true >( sendBufferPtr,
                             edgeGhostsToReceive,
                             nullptr,
                             edgeGhostsToReceive.size(),
                             edgeManager.localToGlobalMap().toSliceConst() );

  packedSize +=
    bufferOps::Pack< true >( sendBufferPtr,
                             faceGhostsToReceive,
                             nullptr,
                             faceGhostsToReceive.size(),
                             faceManager.localToGlobalMap().toSliceConst() );

  elemManager.forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase const & subRegion ) {
      arraySlice1d< localIndex const > const ghostsToReceive =
        subRegion.getNeighborData( m_neighborRank ).ghostsToReceive();
      packedSize +=
        bufferOps::Pack< true >( sendBufferPtr,
                                 ghostsToReceive,
                                 nullptr,
                                 ghostsToReceive.size(),
                                 subRegion.localToGlobalMap().toSliceConst() );
    } );

  GEOSX_ERROR_IF( bufferSize != packedSize,
                  "Allocated Buffer Size is not equal to packed buffer size" );

  this->PostSend( commID );
}

void
NeighborCommunicator::UnpackAndRebuildSyncLists( MeshLevel & mesh,
                                                 int const commID )
{
  GEOSX_MARK_FUNCTION;

  NodeManager & nodeManager = *( mesh.getNodeManager() );
  EdgeManager & edgeManager = *( mesh.getEdgeManager() );
  FaceManager & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager & elemManager = *( mesh.getElemManager() );

  localIndex_array & nodeGhostsToSend =
    nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  localIndex_array & edgeGhostsToSend =
    edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  localIndex_array & faceGhostsToSend =
    faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             nodeGhostsToSend,
                             nodeManager.globalToLocalMap() );

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             edgeGhostsToSend,
                             edgeManager.globalToLocalMap() );

  bufferOps::UnpackSyncList( receiveBufferPtr,
                             faceGhostsToSend,
                             faceManager.globalToLocalMap() );

  nodeManager.SetGhostRankForSenders( m_neighborRank );
  edgeManager.SetGhostRankForSenders( m_neighborRank );
  faceManager.SetGhostRankForSenders( m_neighborRank );

  elemManager.forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase & subRegion ) {
      bufferOps::UnpackSyncList(
        receiveBufferPtr,
        subRegion.getNeighborData( m_neighborRank ).ghostsToSend(),
        subRegion.globalToLocalMap() );

      subRegion.SetGhostRankForSenders( m_neighborRank );
    } );
}

int
NeighborCommunicator::PackCommSizeForSync(
  std::map< string, string_array > const & fieldNames,
  MeshLevel const & mesh,
  int const commID,
  bool on_device )
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = *( mesh.getNodeManager() );
  EdgeManager const & edgeManager = *( mesh.getEdgeManager() );
  FaceManager const & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager const & elemManager = *( mesh.getElemManager() );

  arrayView1d< localIndex const > const & nodeGhostsToSend =
    nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & edgeGhostsToSend =
    edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & faceGhostsToSend =
    faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  int bufferSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    bufferSize +=
      nodeManager.PackSize( fieldNames.at( "node" ), nodeGhostsToSend, 0, on_device );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    bufferSize +=
      edgeManager.PackSize( fieldNames.at( "edge" ), edgeGhostsToSend, 0, on_device );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    bufferSize +=
      faceManager.PackSize( fieldNames.at( "face" ), faceGhostsToSend, 0, on_device );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >(
      [&]( ElementSubRegionBase const & subRegion ) {
        bufferSize += subRegion.PackSize(
          fieldNames.at( "elems" ),
          subRegion.getNeighborData( m_neighborRank ).ghostsToSend(),
          0,
          on_device );
      } );
  }

  this->m_sendBufferSize[commID] = bufferSize;
  return bufferSize;
}

void
NeighborCommunicator::PackCommBufferForSync(
  std::map< string, string_array > const & fieldNames,
  MeshLevel const & mesh,
  int const commID,
  bool on_device )
{
  GEOSX_MARK_FUNCTION;

  NodeManager const & nodeManager = *( mesh.getNodeManager() );
  EdgeManager const & edgeManager = *( mesh.getEdgeManager() );
  FaceManager const & faceManager = *( mesh.getFaceManager() );
  ElementRegionManager const & elemManager = *( mesh.getElemManager() );

  arrayView1d< localIndex const > const & nodeGhostsToSend =
    nodeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & edgeGhostsToSend =
    edgeManager.getNeighborData( m_neighborRank ).ghostsToSend();
  arrayView1d< localIndex const > const & faceGhostsToSend =
    faceManager.getNeighborData( m_neighborRank ).ghostsToSend();

  buffer_type & sendBuffer = SendBuffer( commID );
  int const bufferSize = LvArray::integerConversion< int >( sendBuffer.size() );
  buffer_unit_type * sendBufferPtr = sendBuffer.data();

  int packedSize = 0;
  if( fieldNames.count( "node" ) > 0 )
  {
    packedSize += nodeManager.Pack( sendBufferPtr,
                                    fieldNames.at( "node" ),
                                    nodeGhostsToSend,
                                    0,
                                    on_device );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    packedSize += edgeManager.Pack( sendBufferPtr,
                                    fieldNames.at( "edge" ),
                                    edgeGhostsToSend,
                                    0,
                                    on_device );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    packedSize += faceManager.Pack( sendBufferPtr,
                                    fieldNames.at( "face" ),
                                    faceGhostsToSend,
                                    0,
                                    on_device );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions( [&]( ElementSubRegionBase const & subRegion ) {
      packedSize +=
        subRegion.Pack( sendBufferPtr,
                        fieldNames.at( "elems" ),
                        subRegion.getNeighborData( m_neighborRank ).ghostsToSend(),
                        0,
                        on_device );
    } );
  }

  GEOSX_ERROR_IF_NE( bufferSize, packedSize );
}

void
NeighborCommunicator::SendRecvBuffers( int const commID )
{
  this->MPI_iSendReceive( commID, MPI_COMM_GEOSX );
}

void
NeighborCommunicator::UnpackBufferForSync(
  std::map< string, string_array > const & fieldNames,
  MeshLevel * const mesh,
  int const commID,
  bool on_device )
{
  GEOSX_MARK_FUNCTION;

  buffer_type const & receiveBuffer = ReceiveBuffer( commID );
  buffer_unit_type const * receiveBufferPtr = receiveBuffer.data();

  NodeManager & nodeManager = *( mesh->getNodeManager() );
  EdgeManager & edgeManager = *( mesh->getEdgeManager() );
  FaceManager & faceManager = *( mesh->getFaceManager() );
  ElementRegionManager & elemManager = *( mesh->getElemManager() );

  array1d< localIndex > & nodeGhostsToReceive =
    nodeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  array1d< localIndex > & edgeGhostsToReceive =
    edgeManager.getNeighborData( m_neighborRank ).ghostsToReceive();
  array1d< localIndex > & faceGhostsToReceive =
    faceManager.getNeighborData( m_neighborRank ).ghostsToReceive();

  int unpackedSize = 0;

  if( fieldNames.count( "node" ) > 0 )
  {
    unpackedSize +=
      nodeManager.Unpack( receiveBufferPtr, nodeGhostsToReceive, 0, on_device );
  }

  if( fieldNames.count( "edge" ) > 0 )
  {
    unpackedSize +=
      edgeManager.Unpack( receiveBufferPtr, edgeGhostsToReceive, 0, on_device );
  }

  if( fieldNames.count( "face" ) > 0 )
  {
    unpackedSize +=
      faceManager.Unpack( receiveBufferPtr, faceGhostsToReceive, 0, on_device );
  }

  if( fieldNames.count( "elems" ) > 0 )
  {
    elemManager.forElementSubRegions< ElementSubRegionBase >(
      [&]( ElementSubRegionBase & subRegion ) {
        unpackedSize += subRegion.Unpack(
          receiveBufferPtr,
          subRegion.getNeighborData( m_neighborRank ).ghostsToReceive(),
          0,
          on_device );
      } );
  }
}

} /* namespace geosx */
