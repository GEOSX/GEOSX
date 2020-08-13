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
 * @file NeighborCommunicator.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_
#define GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_

#include "MpiWrapper.hpp"

#include "common/DataTypes.hpp"
#include "dataRepository/ReferenceWrapper.hpp"
#include "LvArray/src/limits.hpp"

namespace geosx
{
inline int
CommTag( int const GEOSX_UNUSED_PARAM( senderRank ),
         int const GEOSX_UNUSED_PARAM( receiverRank ),
         int const comm )
{
  //  int m_size;
  //  MPI_Comm_size( MPI_COMM_GEOSX, &m_size );
  //  return senderRank * m_size + receiverRank + m_size * m_size * comm;
  return comm;
}

class MeshLevel;
class ObjectManagerBase;
class NeighborCommunicator
{
public:
  NeighborCommunicator();

  NeighborCommunicator( int rank ) :
    NeighborCommunicator()
  {
    SetNeighborRank( rank );
  }

  void
  MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                    int const sendSize,
                    MPI_Request & sendRequest,
                    buffer_unit_type * const receiveBuffer,
                    int const receiveSize,
                    MPI_Request & receiveRequest,
                    int const commID,
                    MPI_Comm mpiComm );

  void
  MPI_iSendReceiveBufferSizes( int const commID, MPI_Comm mpiComm );

  void
  MPI_iSendReceiveBufferSizes( int const commID,
                               MPI_Request & mpiSendRequest,
                               MPI_Request & mpiRecvRequest,
                               MPI_Comm mpiComm );

  void
  MPI_iSendReceiveBuffers( int const commID, MPI_Comm mpiComm );

  void
  MPI_iSendReceiveBuffers( int const commID,
                           MPI_Request & mpiSendRequest,
                           MPI_Request & mpiRecvRequest,
                           MPI_Comm mpiComm );

  void
  MPI_iSendReceive( int const commID, MPI_Comm mpiComm );

  void
  MPI_iSendReceive( int const commID,
                    MPI_Request & mpiSendRequest,
                    MPI_Request & mpiRecvRequest,
                    MPI_Comm mpiComm );

  template< typename T >
  void
  MPI_iSendReceive( T const * const sendBuffer,
                    int const sendSize,
                    MPI_Request & sendReq,
                    T * const recvBuffer,
                    int const recvSize,
                    MPI_Request & recvReq,
                    int const commID,
                    MPI_Comm mpiComm )
  {
    MPI_iSendReceive( reinterpret_cast< buffer_unit_type const * >( sendBuffer ),
                      sendSize * sizeof( T ),
                      sendReq,
                      reinterpret_cast< buffer_unit_type * >( recvBuffer ),
                      recvSize * sizeof( T ),
                      recvReq,
                      commID,
                      mpiComm );
  }

  template< typename T >
  void
  MPI_iSendReceive( array1d< T > const & sendBuffer,
                    MPI_Request & sendReq,
                    array1d< T > & recvBuffer,
                    MPI_Request & recvReq,
                    int const commID,
                    MPI_Comm mpiComm );

  template< typename T >
  void
  MPI_iSendReceive( array1d< T > const & sendBuffer,
                    array1d< T > & recvBuffer,
                    int const commID,
                    MPI_Comm mpiComm )
  {
    MPI_iSendReceive( sendBuffer,
                      m_mpiSendBufferRequest[commID],
                      recvBuffer,
                      m_mpiRecvBufferRequest[commID],
                      commID,
                      mpiComm );
  }

  template< typename T >
  void
  MPI_iSendReceiveWait( array1d< T > const & sendBuffer,
                        array1d< T > & recvBuffer,
                        int const commID,
                        MPI_Comm mpiComm )
  {
    MPI_iSendReceive( sendBuffer, recvBuffer, commID, mpiComm );
    MPI_WaitAll( commID );
  }

  void
  MPI_iSendReceive( buffer_unit_type const * const sendBuffer,
                    int const sendSize,
                    int const commID,
                    MPI_Comm mpiComm );

  void
  MPI_WaitAll( int const commID,
               MPI_Request & mpiSendRequest,
               MPI_Status & mpiSendStatus,
               MPI_Request & mpiRecvRequest,
               MPI_Status & mpiReceiveStatus );

  void
  MPI_WaitAll( int const commID );

  /**
   * @brief Post a non-blocking recieve for the buffer size
   *        information for a later call to PostRecv.
   * @note Get the receive MPI_Request by calling
   *       GetSizeRecvRequest with the same CommID
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in
   * @return The return code from the interal call to
   *         MPI_iRecv.
   */
  int
  PostSizeRecv( int const commID );

  /**
   * @brief Get the MPI_Request associated with the last
   *        size-recieve request.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The MPI_Request for the size receive
   */
  MPI_Request
  GetSizeRecvRequest( int const commID );

  /**
   * @brief Post a non-blocking send for the buffer
   *        size information so the neighbor knows
   *        how much data it will recieve from a
   *        subsequent call to PostSend.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iSend.
   */
  int
  PostSizeSend( int const commID );

  /**
   * @brief Post a non-blocking receive for the actual
   *        information being recieved from the neighbor.
   * @note Get the receive request by calling GetRecvRequest
   *       with the same commID.
   * @note Make sure the request for any required size recvs
   *       has completed prior to calling this function, as
   *       it uses the size recv information to resize the
   *       internal recv buffer.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iRecv.
   */
  int
  PostRecv( int const commID );

  /**
   * @brief Get the MPI_Request associated with the last
   *        recieve request posted.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The MPI_Request for the receive.
   */
  MPI_Request
  GetRecvRequest( int const commID );

  /**
   * @brief Post a non-blocking send for the actual information
   *        being sent to the neighboring rank.
   * @param commID The identifier for the pseudo-comm
   *               the communication is taking place in.
   * @return The return code from the internal call to
   *         MPI_iSend.
   */
  int
  PostSend( int const commID );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange ghost information.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void
  PrepareAndSendGhosts( bool const contactActive,
                        int const depth,
                        MeshLevel & meshLevel,
                        int const commID );

  /**
   * Unpack the receive buffer and process ghosting
   *  information recieved from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void
  UnpackGhosts( MeshLevel & meshLevel, int const commID );

  /**
   * Posts non-blocking sends to m_neighborRank for
   *  both the size and regular communication buffers
   *  to exchange synchronization lists.
   *  Additionally posts a non-blocking recv for size
   *  information from m_neighborRank, this size recv
   *  must be completed before PostRecv is called in order
   *  to correctly resize the receive buffer.
   */
  void
  PrepareAndSendSyncLists( MeshLevel const & meshLevel, int const commID );

  /**
   * Unpack the receive buffer and process synchronization
   *  list information recieved from m_neighborRank.
   *  This must be called after PostRecv is called, and
   *  the request associated with that recv has
   *  completed (retrieve the request using GetRecvRequest)
   */
  void
  UnpackAndRebuildSyncLists( MeshLevel & meshLevel, int const CommID );

  void
  PackCommBufferForSync( std::map< string, string_array > const & fieldNames,
                         MeshLevel const & meshLevel,
                         int const commID,
                         bool on_device = false );

  int
  PackCommSizeForSync( std::map< string, string_array > const & fieldNames,
                       MeshLevel const & meshLevel,
                       int const commID,
                       bool on_device = false );

  void
  SendRecvBuffers( int const commID );

  void
  UnpackBufferForSync( std::map< string, string_array > const & fieldNames,
                       MeshLevel * const meshLevel,
                       int const commID,
                       bool on_device = false );

  void
  SetNeighborRank( int const rank )
  {
    m_neighborRank = rank;
  }
  int
  NeighborRank() const
  {
    return m_neighborRank;
  }

  void
  Clear();

  static int constexpr maxComm = 100;

  buffer_type const &
  ReceiveBuffer( int commID ) const
  {
    return m_receiveBuffer[commID];
  }
  buffer_type &
  ReceiveBuffer( int commID )
  {
    return m_receiveBuffer[commID];
  }

  int const &
  ReceiveBufferSize( int commID ) const
  {
    return m_receiveBufferSize[commID];
  }
  int &
  ReceiveBufferSize( int commID )
  {
    return m_receiveBufferSize[commID];
  }

  buffer_type const &
  SendBuffer( int commID ) const
  {
    return m_sendBuffer[commID];
  }
  buffer_type &
  SendBuffer( int commID )
  {
    return m_sendBuffer[commID];
  }

  void
  resizeSendBuffer( int const commID, int const newSize )
  {
    m_sendBufferSize[commID] = newSize;
    m_sendBuffer[commID].resize( newSize );
  }

  void
  resizeRecvBuffer( int const commID, int const newSize )
  {
    m_receiveBufferSize[commID] = newSize;
    m_receiveBuffer[commID].resize( newSize );
  }

  void
  AddNeighborGroupToMesh( MeshLevel & mesh ) const;

private:
  int m_neighborRank;

  int m_sendBufferSize[maxComm];
  int m_receiveBufferSize[maxComm];

  std::vector< buffer_type > m_sendBuffer;
  std::vector< buffer_type > m_receiveBuffer;

  MPI_Request m_mpiSendBufferRequest[maxComm];
  MPI_Request m_mpiRecvBufferRequest[maxComm];

  MPI_Request m_mpiSendSizeRequest[maxComm];
  MPI_Request m_mpiRecvSizeRequest[maxComm];

  MPI_Status m_mpiSendBufferStatus[maxComm];
  MPI_Status m_mpiRecvBufferStatus[maxComm];
};

template< typename T >
void
NeighborCommunicator::MPI_iSendReceive( array1d< T > const & sendBuffer,
                                        MPI_Request & sendReq,
                                        array1d< T > & recvBuffer,
                                        MPI_Request & recvReq,
                                        int const commID,
                                        MPI_Comm mpiComm )
{
  m_sendBufferSize[commID] = LvArray::integerConversion< int >( sendBuffer.size() );

  MPI_iSendReceive( &m_sendBufferSize[commID],
                    1,
                    sendReq,
                    &m_receiveBufferSize[commID],
                    1,
                    recvReq,
                    commID,
                    mpiComm );

  MpiWrapper::Wait( &( recvReq ), &( m_mpiRecvBufferStatus[commID] ) );
  MpiWrapper::Wait( &( sendReq ), &( m_mpiSendBufferStatus[commID] ) );

  recvBuffer.resize( m_receiveBufferSize[commID] );

  MPI_iSendReceive( sendBuffer.data(),
                    m_sendBufferSize[commID],
                    sendReq,
                    recvBuffer.data(),
                    m_receiveBufferSize[commID],
                    recvReq,
                    commID,
                    mpiComm );
}

} /* namespace geosx */

#endif /* GEOSX_MPICOMMUNICATIONS_NEIGHBORCOMMUNICATOR_HPP_ */
