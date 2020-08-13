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
 * @file CommunicationTools.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_
#define GEOSX_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_

#include "MpiWrapper.hpp"

#include "common/DataTypes.hpp"

#include <set>

namespace geosx
{
class ObjectManagerBase;
class NeighborCommunicator;
class MeshLevel;
class ElementRegionManager;

class MPI_iCommData;

class CommunicationTools
{
public:
  CommunicationTools();
  ~CommunicationTools();

  static void
  AssignGlobalIndices( ObjectManagerBase & object,
                       ObjectManagerBase const & compositionObject,
                       std::vector< NeighborCommunicator > & neighbors );

  static void
  AssignNewGlobalIndices( ObjectManagerBase & object,
                          std::set< localIndex > const & indexList );

  static void
  AssignNewGlobalIndices(
    ElementRegionManager & elementManager,
    std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems );

  static void
  FindGhosts( MeshLevel & meshLevel,
              std::vector< NeighborCommunicator > & neighbors,
              bool use_nonblocking );

  static std::set< int > &
  getFreeCommIDs();
  static int
  reserveCommID();
  static void
  releaseCommID( int & ID );

  static void
  FindMatchedPartitionBoundaryObjects(
    ObjectManagerBase * const group,
    std::vector< NeighborCommunicator > & allNeighbors );

  static void
  SynchronizeFields( const std::map< string, string_array > & fieldNames,
                     MeshLevel * const mesh,
                     std::vector< NeighborCommunicator > & allNeighbors,
                     bool on_device = false );

  static void
  SynchronizePackSendRecvSizes(
    const std::map< string, string_array > & fieldNames,
    MeshLevel * const mesh,
    std::vector< NeighborCommunicator > & neighbors,
    MPI_iCommData & icomm,
    bool on_device = false );

  static void
  SynchronizePackSendRecv(
    const std::map< string, string_array > & fieldNames,
    MeshLevel * const mesh,
    std::vector< NeighborCommunicator > & allNeighbors,
    MPI_iCommData & icomm,
    bool on_device = false );

  static void
  SynchronizeUnpack( MeshLevel * const mesh,
                     std::vector< NeighborCommunicator > & neighbors,
                     MPI_iCommData & icomm,
                     bool on_device = false );
};

class MPI_iCommData
{
public:
  MPI_iCommData() :
    size( 0 ),
    commID( -1 ),
    sizeCommID( -1 ),
    fieldNames(),
    mpiSendBufferRequest(),
    mpiRecvBufferRequest(),
    mpiSendBufferStatus(),
    mpiRecvBufferStatus()
  {
    commID = CommunicationTools::reserveCommID();
    sizeCommID = CommunicationTools::reserveCommID();
  }

  ~MPI_iCommData()
  {
    if( commID >= 0 )
    {
      CommunicationTools::releaseCommID( commID );
    }

    if( sizeCommID >= 0 )
    {
      CommunicationTools::releaseCommID( sizeCommID );
    }
  }

  void
  resize( localIndex numMessages )
  {
    mpiSendBufferRequest.resize( numMessages );
    mpiRecvBufferRequest.resize( numMessages );
    mpiSendBufferStatus.resize( numMessages );
    mpiRecvBufferStatus.resize( numMessages );
    mpiSizeSendBufferRequest.resize( numMessages );
    mpiSizeRecvBufferRequest.resize( numMessages );
    mpiSizeSendBufferStatus.resize( numMessages );
    mpiSizeRecvBufferStatus.resize( numMessages );
    size = static_cast< int >( numMessages );
  }

  int size;
  int commID;
  int sizeCommID;
  std::map< string, string_array > fieldNames;

  array1d< MPI_Request > mpiSendBufferRequest;
  array1d< MPI_Request > mpiRecvBufferRequest;
  array1d< MPI_Status > mpiSendBufferStatus;
  array1d< MPI_Status > mpiRecvBufferStatus;

  array1d< MPI_Request > mpiSizeSendBufferRequest;
  array1d< MPI_Request > mpiSizeRecvBufferRequest;
  array1d< MPI_Status > mpiSizeSendBufferStatus;
  array1d< MPI_Status > mpiSizeRecvBufferStatus;
};

} /* namespace geosx */

#endif /* GEOSX_MPICOMMUNICATIONS_COMMUNICATIONTOOLS_HPP_ */
