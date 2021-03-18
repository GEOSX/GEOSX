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
 * @file CommunicationTools.cpp
 *
 */

#include "mesh/mpiCommunications/CommunicationTools.hpp"

#include "common/TimingMacros.hpp"
#include "mesh/mpiCommunications/MPI_iCommData.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/MeshLevel.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include <algorithm>

namespace geosx
{


CommunicationTools * CommunicationTools::m_instance = nullptr;

using namespace dataRepository;

CommunicationTools::CommunicationTools()
{
  for( int i = 0; i < NeighborCommunicator::maxComm; ++i )
  {
    m_freeCommIDs.insert( i );
  }

  GEOSX_ERROR_IF( m_instance != nullptr, "Only one CommunicationTools can exist at a time." );
  m_instance = this;
}

CommunicationTools::~CommunicationTools()
{
  GEOSX_ERROR_IF( m_instance != this, "m_instance != this should not be possible." );
  m_instance = nullptr;
}

CommunicationTools & CommunicationTools::getInstance()
{
  GEOSX_ERROR_IF( m_instance == nullptr,
                  "CommunicationTools has not been constructed, or is already been destructed." );
  return *m_instance;
}


void CommunicationTools::assignGlobalIndices( ObjectManagerBase & object,
                                              NodeManager const & compositionObject,
                                              std::vector< NeighborCommunicator > & neighbors )
{
  GEOSX_MARK_FUNCTION;
  arrayView1d< integer > const & ghostRank = object.ghostRank();
  ghostRank.setValues< serialPolicy >( -2 );

  int const commRank = MpiWrapper::commRank();

  localIndex const numberOfObjectsHere = object.size();
  globalIndex const offset = MpiWrapper::prefixSum< globalIndex >( numberOfObjectsHere );

  arrayView1d< globalIndex > const & localToGlobal = object.localToGlobalMap();

  // set the global indices as if they were all local to this process
  for( localIndex a = 0; a < object.size(); ++a )
  {
    localToGlobal[a] = offset + a;
  }

  // get the relation to the composition object used that will be used to identify the main object. For example,
  // a face can be identified by its nodes.
  std::vector< std::vector< globalIndex > > objectToCompositionObject;
  object.extractMapFromObjectForAssignGlobalIndexNumbers( compositionObject, objectToCompositionObject );

  // now arrange the data from objectToCompositionObject into a map "indexByFirstCompositionIndex", such that the key
  // is the lowest global index of the composition object that make up this object. The value of the map is a pair, with
  // the array being the remaining composition object global indices, and the second being the global index of the
  // object
  // itself.
  map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, localIndex > > > indexByFirstCompositionIndex;

  localIndex bufferSize = 0;
  for( std::size_t a = 0; a < objectToCompositionObject.size(); ++a )
  {
    if( objectToCompositionObject[a].size() > 0 )
    {
      // set nodelist array
      std::vector< globalIndex > const & nodeList = objectToCompositionObject[a];

      // grab the first global index of the composition objects
      const globalIndex firstCompositionIndex = nodeList[0];

      // create a temporary to hold the pair
      std::pair< std::vector< globalIndex >, globalIndex > tempComp;

      // fill the array with the remaining composition object global indices
      tempComp.first.insert( tempComp.first.begin(), nodeList.begin() + 1, nodeList.end() );

      // set the second value of the pair to the localIndex of the object.
      tempComp.second = a;

      // push the tempComp onto the map.
      indexByFirstCompositionIndex[firstCompositionIndex].emplace_back( std::move( tempComp ) );
      bufferSize += 2 + nodeList.size();
    }
  }

  globalIndex_array objectToCompositionObjectSendBuffer;
  objectToCompositionObjectSendBuffer.reserve( bufferSize );

  // put the map into a buffer
  for( std::size_t a = 0; a < objectToCompositionObject.size(); ++a )
  {
    if( objectToCompositionObject[a].size() > 0 )
    {
      std::vector< globalIndex > const & nodeList = objectToCompositionObject[a];
      objectToCompositionObjectSendBuffer.emplace_back( nodeList.size() );
      objectToCompositionObjectSendBuffer.emplace_back( localToGlobal[a] );
      for( std::size_t b = 0; b < nodeList.size(); ++b )
      {
        objectToCompositionObjectSendBuffer.emplace_back( nodeList[b] );
      }
    }
  }

  MPI_iCommData commData( getCommID() );
  commData.resize( neighbors.size() );

  array1d< int >  receiveBufferSizes( neighbors.size());
  array1d< globalIndex_array > receiveBuffers( neighbors.size());

  int const sendSize = LvArray::integerConversion< int const >( objectToCompositionObjectSendBuffer.size() );

  for( std::size_t neighborIndex = 0; neighborIndex < neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    neighbor.mpiISendReceive( &sendSize,
                              1,
                              commData.mpiSendBufferSizeRequest( neighborIndex ),
                              &(receiveBufferSizes[neighborIndex]),
                              1,
                              commData.mpiRecvBufferSizeRequest( neighborIndex ),
                              commData.commID(),
                              MPI_COMM_GEOSX );
  }


  for( std::size_t count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    receiveBuffers[neighborIndex].resize( receiveBufferSizes[neighborIndex] );
    neighbor.mpiISendReceive( objectToCompositionObjectSendBuffer.data(),
                              sendSize,
                              commData.mpiSendBufferRequest( neighborIndex ),
                              receiveBuffers[neighborIndex].data(),
                              receiveBufferSizes[neighborIndex],
                              commData.mpiRecvBufferRequest( neighborIndex ),
                              commData.commID(),
                              MPI_COMM_GEOSX );

  }

  // unpack the data from neighbor->tempNeighborData.neighborNumbers[DomainPartition::FiniteElementNodeManager] to
  // the local arrays

  // object to receive the neighbor data
  // this baby is an Array (for each neighbor) of maps, with the key of lowest composition index, and a value
  // containing an array containing the std::pairs of the remaining composition indices, and the globalIndex of the
  // object.
  std::vector< map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, globalIndex > > > >
  neighborCompositionObjects( neighbors.size() );

  for( std::size_t count=0; count<neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( commData.size(),
                         commData.mpiRecvBufferRequest(),
                         &neighborIndex,
                         commData.mpiRecvBufferStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    globalIndex const * recBuffer = receiveBuffers[neighborIndex].data();
    localIndex recBufferSize = receiveBufferSizes[neighborIndex];
    globalIndex const * endBuffer = recBuffer + recBufferSize;
    // iterate over data that was just received
    while( recBuffer < endBuffer )
    {
      // the first thing packed was the data size for a given object
      localIndex dataSize = LvArray::integerConversion< localIndex >( *recBuffer++ );

      // the second thing packed was the globalIndex of that object
      const globalIndex neighborGlobalIndex = *( recBuffer++ );

      // the global indices of the composition objects were next. they are ordered, so the lowest one is first.
      const globalIndex firstCompositionIndex = *( recBuffer++ );

      // the remaining composition object indices.
      std::vector< globalIndex > temp;
      for( localIndex b = 1; b < dataSize; ++b )
      {
        temp.emplace_back( *( recBuffer++ ) );
      }

      // fill neighborCompositionObjects
      std::pair< std::vector< globalIndex >, globalIndex >
      tempComp( std::make_pair( std::move( temp ), std::move( neighborGlobalIndex ) ) );

      neighborCompositionObjects[neighborIndex][firstCompositionIndex].emplace_back( tempComp );
    }

    // Set iterators to the beginning of each indexByFirstCompositionIndex,
    // and neighborCompositionObjects[neighborNum].
    map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, localIndex > > >::const_iterator
      iter_local = indexByFirstCompositionIndex.begin();
    map< globalIndex, std::vector< std::pair< std::vector< globalIndex >, globalIndex > > >::const_iterator
      iter_neighbor = neighborCompositionObjects[neighborIndex].begin();

    // now we continue the while loop as long as both of our iterators are in range.
    while( iter_local != indexByFirstCompositionIndex.end() &&
           iter_neighbor != neighborCompositionObjects[neighborIndex].end() )
    {
      // check to see if the map keys (first composition index) are the same.
      if( iter_local->first == iter_neighbor->first )
      {
        // first we loop over all local composition arrays (objects with the matched key)
        for( std::vector< std::pair< std::vector< globalIndex >, localIndex > >::const_iterator
             iter_local2 = iter_local->second.begin();
             iter_local2 != iter_local->second.end(); ++iter_local2 )
        {
          // and loop over all of the neighbor composition arrays (objects with the matched key)
          for( std::vector< std::pair< std::vector< globalIndex >, globalIndex > >::const_iterator
               iter_neighbor2 = iter_neighbor->second.begin();
               iter_neighbor2 != iter_neighbor->second.end();
               ++iter_neighbor2 )
          {
            // now compare the composition arrays
            if( iter_local2->first.size() == iter_neighbor2->first.size() &&
                std::equal( iter_local2->first.begin(),
                            iter_local2->first.end(),
                            iter_neighbor2->first.begin() ) )
            {
              // they are equal, so we need to overwrite the global index for the object
              if( iter_neighbor2->second < localToGlobal[iter_local2->second] )
              {
                if( neighbor.neighborRank() < commRank )
                {
                  localToGlobal[iter_local2->second] = iter_neighbor2->second;
                  ghostRank[iter_local2->second] = neighbor.neighborRank();
                }
                else
                {
                  ghostRank[iter_local2->second] = -1;
                }

              }

              // we should break out of the iter_local2 loop since we aren't going to find another match.
              break;
            }
          }
        }
        ++iter_local;
        ++iter_neighbor;
      }
      else if( iter_local->first < iter_neighbor->first )
      {
        ++iter_local;
      }
      else if( iter_local->first > iter_neighbor->first )
      {
        ++iter_neighbor;
      }
    }
  }

  object.constructGlobalToLocalMap();

  object.setMaxGlobalIndex();
}

void CommunicationTools::assignNewGlobalIndices( ObjectManagerBase & object,
                                                 std::set< localIndex > const & indexList )
{
  // TODO: This should be done with a prefix sum!
  int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOSX );
  localIndex numberOfNewObjectsHere = indexList.size();
  localIndex_array numberOfNewObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MpiWrapper::allGather( numberOfNewObjectsHere, numberOfNewObjects );


  glocalIndexOffset[0] = 0;
  for( int rank = 1; rank < commSize; ++rank )
  {
    glocalIndexOffset[rank] = glocalIndexOffset[rank - 1] + numberOfNewObjects[rank - 1];
  }

  arrayView1d< globalIndex > const & localToGlobal = object.localToGlobalMap();

  localIndex nIndicesAssigned = 0;
  for( localIndex const newLocalIndex : indexList )
  {
    GEOSX_ERROR_IF( localToGlobal[newLocalIndex] != -1,
                    "Local object " << newLocalIndex << " should be new but already has a global index "
                                    << localToGlobal[newLocalIndex] );

    localToGlobal[newLocalIndex] = object.maxGlobalIndex() + glocalIndexOffset[thisRank] + nIndicesAssigned + 1;
    object.updateGlobalToLocalMap( newLocalIndex );

    nIndicesAssigned += 1;
  }

  object.setMaxGlobalIndex();
}

void
CommunicationTools::
  assignNewGlobalIndices( ElementRegionManager & elementManager,
                          std::map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems )
{
  // TODO: This should be done with a prefix sum!
  int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  localIndex numberOfNewObjectsHere = 0;

  for( auto const & iter : newElems )
  {
    std::set< localIndex > const & indexList = iter.second;
    numberOfNewObjectsHere += indexList.size();
  }

  localIndex_array numberOfNewObjects( commSize );
  localIndex_array glocalIndexOffset( commSize );
  MpiWrapper::allGather( numberOfNewObjectsHere, numberOfNewObjects );


  glocalIndexOffset[0] = 0;
  for( int rank = 1; rank < commSize; ++rank )
  {
    glocalIndexOffset[rank] = glocalIndexOffset[rank - 1] + numberOfNewObjects[rank - 1];
  }

  localIndex nIndicesAssigned = 0;
  for( auto const & iter : newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set< localIndex > const & indexList = iter.second;

    ElementSubRegionBase & subRegion = elementManager.getRegion( er ).getSubRegion( esr );
    arrayView1d< globalIndex > const & localToGlobal = subRegion.localToGlobalMap();

    for( localIndex const newLocalIndex : indexList )
    {
      GEOSX_ERROR_IF( localToGlobal[newLocalIndex] != -1,
                      "Local object " << newLocalIndex << " should be new but already has a global index "
                                      << localToGlobal[newLocalIndex] );

      localToGlobal[newLocalIndex] = elementManager.maxGlobalIndex() + glocalIndexOffset[thisRank] + nIndicesAssigned + 1;
      subRegion.updateGlobalToLocalMap( newLocalIndex );

      nIndicesAssigned += 1;
    }
  }

  elementManager.setMaxGlobalIndex();
}

void
CommunicationTools::
  findMatchedPartitionBoundaryObjects( ObjectManagerBase & objectManager,
                                       std::vector< NeighborCommunicator > & allNeighbors )
{
  GEOSX_MARK_FUNCTION;
  arrayView1d< integer > const & domainBoundaryIndicator = objectManager.getDomainBoundaryIndicator();

  array1d< globalIndex > globalPartitionBoundaryObjectsIndices;
  objectManager.constructGlobalListOfBoundaryObjects( globalPartitionBoundaryObjectsIndices );


  // send the size of the partitionBoundaryObjects to neighbors
  {
    array1d< array1d< globalIndex > > neighborPartitionBoundaryObjects( allNeighbors.size() );

    MPI_iCommData commData( getCommID() );
    int const commID = commData.commID();
    std::size_t const numNeighbors = allNeighbors.size();
    commData.resize( numNeighbors );
    for( std::size_t i=0; i<numNeighbors; ++i )
    {
      allNeighbors[i].mpiISendReceiveSizes( globalPartitionBoundaryObjectsIndices,
                                            commData.mpiSendBufferSizeRequest( i ),
                                            commData.mpiRecvBufferSizeRequest( i ),
                                            commID,
                                            MPI_COMM_GEOSX );
    }

    MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferSizeRequest(), commData.mpiSendBufferSizeStatus() );
    MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferSizeRequest(), commData.mpiRecvBufferSizeStatus() );

    for( std::size_t i=0; i<numNeighbors; ++i )
    {
      allNeighbors[i].mpiISendReceiveData( globalPartitionBoundaryObjectsIndices,
                                           commData.mpiSendBufferRequest( i ),
                                           neighborPartitionBoundaryObjects[i],
                                           commData.mpiRecvBufferRequest( i ),
                                           commID,
                                           MPI_COMM_GEOSX );
    }
    MpiWrapper::waitAll( numNeighbors, commData.mpiSendBufferRequest(), commData.mpiSendBufferStatus() );
    MpiWrapper::waitAll( numNeighbors, commData.mpiRecvBufferRequest(), commData.mpiRecvBufferStatus() );

    for( std::size_t i=0; i<allNeighbors.size(); ++i )
    {
      NeighborCommunicator & neighbor = allNeighbors[i];
      localIndex_array & matchedPartitionBoundaryObjects = objectManager.getNeighborData( neighbor.neighborRank() ).matchedPartitionBoundary();

      localIndex localCounter = 0;
      localIndex neighborCounter = 0;
      while( localCounter < globalPartitionBoundaryObjectsIndices.size() &&
             neighborCounter < neighborPartitionBoundaryObjects[i].size() )
      {
        if( globalPartitionBoundaryObjectsIndices[localCounter] == neighborPartitionBoundaryObjects[i][neighborCounter] )
        {
          localIndex const localMatchedIndex = objectManager.globalToLocalMap( globalPartitionBoundaryObjectsIndices[localCounter] );
          matchedPartitionBoundaryObjects.emplace_back( localMatchedIndex );
          domainBoundaryIndicator[ localMatchedIndex ] = 2;
          ++localCounter;
          ++neighborCounter;
        }
        else if( globalPartitionBoundaryObjectsIndices[localCounter] > neighborPartitionBoundaryObjects[i][neighborCounter] )
        {
          ++neighborCounter;
        }
        else
        {
          ++localCounter;
        }
      }
    }
  }
}

/**
 * @brief Check that the provided object's ghosts are consistent with the neighbors.
 * @param objectManager the owner of the objects to check.
 * @param neighbors list of all the neighbors.
 */
void verifyGhostingConsistency( ObjectManagerBase const & objectManager,
                                std::vector< NeighborCommunicator > const & neighbors )
{
  arrayView1d< integer const > const & ghostRank = objectManager.ghostRank();

  /// Variable to track if an error has occurred.
  bool error = false;

  /// For each neighbor make sure that the ghost rank is consistent with the send and receive lists.
  for( NeighborCommunicator const & neighbor : neighbors )
  {
    int const neighborRank = neighbor.neighborRank();
    NeighborData const & neighborData = objectManager.getNeighborData( neighborRank );

    arrayView1d< localIndex const > const & recvList = neighborData.ghostsToReceive();
    for( localIndex const recvIdx : recvList )
    {
      if( ghostRank[ recvIdx ] != neighborRank )
      {
        error = true;
        GEOSX_LOG_RANK( "Receiving " << recvIdx << " from " << neighborRank <<
                        " but ghostRank[ " << recvIdx << " ] is " << ghostRank[ recvIdx ] );
      }
    }

    arrayView1d< localIndex const > const & sendList = neighborData.ghostsToSend();
    for( localIndex const sendIdx : sendList )
    {
      if( ghostRank[ sendIdx ] != -1 )
      {
        error = true;
        GEOSX_LOG_RANK( "Sending " << sendIdx << " to " << neighborRank <<
                        " but ghostRank[ " << sendIdx << " ] is " << ghostRank[ sendIdx ] );
      }
    }

    arrayView1d< std::pair< globalIndex, int > const > const & nonLocalGhosts = neighborData.nonLocalGhosts();
    if( nonLocalGhosts.size() != 0 )
    {
      error = true;
      GEOSX_LOG_RANK( "Expected to send 0 non local ghosts to rank " << neighborRank <<
                      " but sending " << nonLocalGhosts.size() );
    }
  }

  GEOSX_ERROR_IF( error, "Encountered a ghosting inconsistency in " << objectManager.getName() );
}

/**
 * @brief Remove the given indices from the communication list.
 * @param indicesToAdd the local indices of objects to be removed.
 * @param commIndices the local indices of the existing objects to be communicated.
 */
void removeFromCommList( std::vector< localIndex > const & indicesToRemove, array1d< localIndex > & commIndices )
{
  localIndex * const itr = std::remove_if( commIndices.begin(), commIndices.end(), [&indicesToRemove]( localIndex const idx )
  {
    return std::find( indicesToRemove.begin(), indicesToRemove.end(), idx ) != indicesToRemove.end();
  } );

  localIndex const nRemoved = commIndices.end() - itr;
  GEOSX_ERROR_IF_NE( nRemoved, localIndex( indicesToRemove.size() ) );
  commIndices.resize( commIndices.size() - nRemoved );
}

/**
 * @brief Fix up second neighbor ghosting issues by modifying the receive lists and ghost rank.
 * @param objectManager the owner of the objects to fix up.
 * @param neighbors array of neighbors.
 */
void fixReceiveLists( ObjectManagerBase & objectManager,
                      std::vector< NeighborCommunicator > const & neighbors )
{
  int nonLocalGhostsTag = 45;

  std::vector< MPI_Request > nonLocalGhostsRequests( neighbors.size() );

  /// For each neighbor send them the indices of their ghosts that they mistakenly believe are owned by this rank.
  for( std::size_t i = 0; i < neighbors.size(); ++i )
  {
    int const neighborRank = neighbors[ i ].neighborRank();

    MpiWrapper::iSend( objectManager.getNeighborData( neighborRank ).nonLocalGhosts().toViewConst(),
                       neighborRank,
                       nonLocalGhostsTag,
                       MPI_COMM_GEOSX,
                       &nonLocalGhostsRequests[ i ] );
  }

  for( NeighborCommunicator const & neighbor : neighbors )
  {
    int const neighborRank = neighbor.neighborRank();

    /// Receive the lists of ghosts we mistakenly though were owned by this neighbor.
    array1d< std::pair< globalIndex, int > > ghostsFromSecondNeighbor;
    MpiWrapper::recv( ghostsFromSecondNeighbor,
                      neighborRank,
                      nonLocalGhostsTag,
                      MPI_COMM_GEOSX,
                      MPI_STATUS_IGNORE );

    /// Array of ghosts to fix.
    std::vector< localIndex > ghostsToFix;

    /// Map from owning MPI rank to an array of local objects we need to fix.
    std::unordered_map< int, std::vector< localIndex > > ghostsBySecondNeighbor;

    arrayView1d< integer > const & ghostRank = objectManager.ghostRank();

    /// Populate ghostsToFix and ghostsBySecondNeighbor while also updating ghostRank.
    for( std::pair< globalIndex, int > const & pair : ghostsFromSecondNeighbor )
    {
      localIndex const lid = objectManager.globalToLocalMap( pair.first );
      ghostsBySecondNeighbor[ pair.second ].emplace_back( lid );
      ghostsToFix.emplace_back( lid );
      ghostRank[ lid ] = pair.second;
    }

    /// Remove the ghosts to fix from the neighbor's receive list.
    removeFromCommList( ghostsToFix, objectManager.getNeighborData( neighborRank ).ghostsToReceive() );

    /// Iterate over the ranks that own the objects. For each rank add the new objects to the receive list.
    for( std::pair< int const, std::vector< localIndex > > const & pair : ghostsBySecondNeighbor )
    {
      array1d< localIndex > & trueOwnerRecvList = objectManager.getNeighborData( pair.first ).ghostsToReceive();
      trueOwnerRecvList.insert( trueOwnerRecvList.size(), pair.second.begin(), pair.second.end() );
    }
  }

  /// Wait on the initial send requests.
  MpiWrapper::waitAll( nonLocalGhostsRequests.size(), nonLocalGhostsRequests.data(), MPI_STATUSES_IGNORE );
}

/**
 * @brief Remove neighbors if there is no communication between them
 * @param nodeManager the NodeManager.
 * @param edgeManager the EdgeManager.
 * @param faceManager the FaceManager.
 * @param elemManager the ElementRegionManager.
 * @param neighbors the list of NeighborCommunicators, may be modified.
 */
void removeUnusedNeighbors( NodeManager & nodeManager,
                            EdgeManager & edgeManager,
                            FaceManager & faceManager,
                            ElementRegionManager & elemManager,
                            std::vector< NeighborCommunicator > & neighbors )
{
  for( std::size_t i = 0; i < neighbors.size(); )
  {
    int const neighborRank = neighbors[ i ].neighborRank();

    bool used = false;

    used = used || nodeManager.getNeighborData( neighborRank ).communicationExists();

    used = used || edgeManager.getNeighborData( neighborRank ).communicationExists();

    used = used || faceManager.getNeighborData( neighborRank ).communicationExists();

    elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
    {
      used = used || subRegion.getNeighborData( neighborRank ).communicationExists();
    } );

    if( used )
    {
      ++i;
    }
    else
    {
      nodeManager.removeNeighbor( neighborRank );
      edgeManager.removeNeighbor( neighborRank );
      faceManager.removeNeighbor( neighborRank );

      elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
      {
        subRegion.removeNeighbor( neighborRank );
      } );

      neighbors.erase( neighbors.begin() + i );
    }
  }
}

/**
 * @brief Complete each communication phase in order or with a wait any.
 * @param n the number of requests in each phase.
 * @param phases list of phases.
 * @param unorderedComms if true complete the communications of each phase in the order they are received.
 */
void waitOrderedOrWaitAll( int const n, std::vector< std::function< MPI_Request ( int ) > > const & phases, bool const unorderedComms )
{
  if( unorderedComms )
  {
    MpiWrapper::activeWaitSomeCompletePhase( n, phases );
  }
  else
  {
    MpiWrapper::activeWaitOrderedCompletePhase( n, phases );
  }
}

void CommunicationTools::setupGhosts( MeshLevel & meshLevel,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      bool const unorderedComms )
{
  GEOSX_MARK_FUNCTION;
  MPI_iCommData commData( getCommID() );
  commData.resize( neighbors.size() );

  NodeManager & nodeManager = meshLevel.getNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  auto sendGhosts = [&] ( int idx )
  {
    neighbors[idx].prepareAndSendGhosts( false,
                                         1,
                                         meshLevel,
                                         commData.commID(),
                                         commData.mpiRecvBufferSizeRequest( idx ),
                                         commData.mpiSendBufferSizeRequest( idx ),
                                         commData.mpiSendBufferRequest( idx ) );

    return commData.mpiRecvBufferSizeRequest( idx );
  };
  auto postRecv = [&] ( int idx )
  {
    neighbors[idx].postRecv( commData.commID(),
                             commData.mpiRecvBufferRequest( idx ) );
    return commData.mpiRecvBufferRequest( idx );
  };
  auto unpackGhosts = [&] ( int idx )
  {
    neighbors[idx].unpackGhosts( meshLevel, commData.commID() );
    return MPI_REQUEST_NULL;
  };

  waitOrderedOrWaitAll( neighbors.size(), { sendGhosts, postRecv, unpackGhosts }, unorderedComms );

  nodeManager.setReceiveLists();
  edgeManager.setReceiveLists();
  faceManager.setReceiveLists();

  // at present removing this barrier allows a nondeterministic mpi error to happen on lassen
  //   it occurs less than 5% of the time and happens when a process enters the recv phase in
  //   the sync list exchange while another process still has not unpacked the ghosts received from
  //   the first process. Depending on the mpi implementation the sync send from the first process
  //   can be recv'd by the second process instead of the ghost send which has already been sent but
  //   not necessarily received.
  // Some restructuring to ensure this can't happen ( can also probably just change the send/recv tagging )
  //   can eliminate this. But at present runtimes are the same in either case, as time is mostly just
  //   shifted from the waitall in UnpackAndRebuildSyncLists since the processes are more 'in-sync' when
  //   hitting that point after introducing this barrier.
  MpiWrapper::barrier();

  auto sendSyncLists = [&] ( int idx )
  {
    MPI_Request & mpiSizeRecvRequest = commData.mpiRecvBufferSizeRequest( idx );
    MPI_Request & mpiSizeSendRequest = commData.mpiSendBufferSizeRequest( idx );
    MPI_Request & mpiSendRequest = commData.mpiSendBufferRequest( idx );
    neighbors[idx].prepareAndSendSyncLists( meshLevel,
                                            commData.commID(),
                                            mpiSizeRecvRequest,
                                            mpiSizeSendRequest,
                                            mpiSendRequest );
    return mpiSizeRecvRequest;
  };
  auto rebuildSyncLists = [&] ( int idx )
  {
    neighbors[idx].unpackAndRebuildSyncLists( meshLevel,
                                              commData.commID() );
    return MPI_REQUEST_NULL;
  };

  waitOrderedOrWaitAll( neighbors.size(), { sendSyncLists, postRecv, rebuildSyncLists }, unorderedComms );

  fixReceiveLists( nodeManager, neighbors );
  fixReceiveLists( edgeManager, neighbors );
  fixReceiveLists( faceManager, neighbors );

  waitOrderedOrWaitAll( neighbors.size(), { sendSyncLists, postRecv, rebuildSyncLists }, unorderedComms );

  nodeManager.fixUpDownMaps( false );
  verifyGhostingConsistency( nodeManager, neighbors );
  edgeManager.fixUpDownMaps( false );
  verifyGhostingConsistency( edgeManager, neighbors );
  faceManager.fixUpDownMaps( false );
  verifyGhostingConsistency( faceManager, neighbors );
  elemManager.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase & subRegion )
  {
    subRegion.fixUpDownMaps( false );
    verifyGhostingConsistency( subRegion, neighbors );
  } );

  removeUnusedNeighbors( nodeManager, edgeManager, faceManager, elemManager, neighbors );

  nodeManager.compressRelationMaps();
  edgeManager.compressRelationMaps();
  faceManager.compressRelationMaps();
}

void CommunicationTools::synchronizePackSendRecvSizes( string_array const & fieldNames,
                                                       ObjectManagerBase & manager,
                                                       std::vector< NeighborCommunicator > & neighbors,
                                                       MPI_iCommData & icomm,
                                                       bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  icomm.fieldNames().insert( { manager.getName(), fieldNames } );
  icomm.resize( neighbors.size() );

  parallelDeviceEvents events;
  for( std::size_t neighborIndex = 0; neighborIndex < neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    int const bufferSize = neighbor.packCommSizeForSync( fieldNames, manager, icomm.commID(), onDevice, events );

    neighbor.mpiISendReceiveBufferSizes( icomm.commID(),
                                         icomm.mpiSendBufferSizeRequest( neighborIndex ),
                                         icomm.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

    neighbor.resizeSendBuffer( icomm.commID(), bufferSize );
  }
  waitAllDeviceEvents( events );
}

void CommunicationTools::synchronizePackSendRecvSizes( const std::map< string, string_array > & fieldNames,
                                                       MeshLevel & mesh,
                                                       std::vector< NeighborCommunicator > & neighbors,
                                                       MPI_iCommData & icomm,
                                                       bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  icomm.fieldNames().insert( fieldNames.begin(), fieldNames.end() );
  icomm.resize( neighbors.size() );

  parallelDeviceEvents events;
  for( std::size_t neighborIndex = 0; neighborIndex < neighbors.size(); ++neighborIndex )
  {
    NeighborCommunicator & neighbor = neighbors[neighborIndex];
    int const bufferSize = neighbor.packCommSizeForSync( fieldNames, mesh, icomm.commID(), onDevice, events );

    neighbor.mpiISendReceiveBufferSizes( icomm.commID(),
                                         icomm.mpiSendBufferSizeRequest( neighborIndex ),
                                         icomm.mpiRecvBufferSizeRequest( neighborIndex ),
                                         MPI_COMM_GEOSX );

    neighbor.resizeSendBuffer( icomm.commID(), bufferSize );
  }
  waitAllDeviceEvents( events );
}

void CommunicationTools::asyncPack( string_array const & fieldNames,
                                    ObjectManagerBase & manager,
                                    std::vector< NeighborCommunicator > & neighbors,
                                    MPI_iCommData & icomm,
                                    bool onDevice,
                                    parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;
  for( NeighborCommunicator & neighbor : neighbors )
  {
    neighbor.packCommBufferForSync( fieldNames, manager, icomm.commID(), onDevice, events );
  }
}

void CommunicationTools::asyncPack( const std::map< string, string_array > & fieldNames,
                                    MeshLevel & mesh,
                                    std::vector< NeighborCommunicator > & neighbors,
                                    MPI_iCommData & icomm,
                                    bool onDevice,
                                    parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;
  for( NeighborCommunicator & neighbor : neighbors )
  {
    neighbor.packCommBufferForSync( fieldNames, mesh, icomm.commID(), onDevice, events );
  }
}

void CommunicationTools::asyncSendRecv( std::vector< NeighborCommunicator > & neighbors,
                                        MPI_iCommData & icomm,
                                        bool onDevice,
                                        parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;
  if( onDevice )
  {
    waitAllDeviceEvents( events );
  }

  // could swap this to test and make this function call async as well, only launch the sends/recvs for
  // those we've already recv'd sizing for, go back to some usefule compute / launch some other compute, then
  // check this again
  for( std::size_t count = 0; count < neighbors.size(); ++count )
  {
    int neighborIndex;
    MpiWrapper::waitAny( icomm.size(),
                         icomm.mpiRecvBufferSizeRequest(),
                         &neighborIndex,
                         icomm.mpiRecvBufferSizeStatus() );

    NeighborCommunicator & neighbor = neighbors[neighborIndex];

    neighbor.mpiISendReceiveBuffers( icomm.commID(),
                                     icomm.mpiSendBufferRequest( neighborIndex ),
                                     icomm.mpiRecvBufferRequest( neighborIndex ),
                                     MPI_COMM_GEOSX );
  }
}

void CommunicationTools::synchronizePackSendRecv( string_array const & fieldNames,
                                                  ObjectManagerBase & manager,
                                                  std::vector< NeighborCommunicator > & neighbors,
                                                  MPI_iCommData & icomm,
                                                  bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  parallelDeviceEvents events;
  asyncPack( fieldNames, manager, neighbors, icomm, onDevice, events );
  asyncSendRecv( neighbors, icomm, onDevice, events );
}

void CommunicationTools::synchronizePackSendRecv( const std::map< string, string_array > & fieldNames,
                                                  MeshLevel & mesh,
                                                  std::vector< NeighborCommunicator > & neighbors,
                                                  MPI_iCommData & icomm,
                                                  bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  parallelDeviceEvents events;
  asyncPack( fieldNames, mesh, neighbors, icomm, onDevice, events );
  asyncSendRecv( neighbors, icomm, onDevice, events );
}

bool CommunicationTools::asyncUnpack( ObjectManagerBase & manager,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      MPI_iCommData & icomm,
                                      bool onDevice,
                                      parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  int recvCount = 0;
  std::vector< int > neighborIndices;
  neighborIndices.reserve( icomm.size() );
  MpiWrapper::testSome( icomm.size(),
                        icomm.mpiRecvBufferRequest(),
                        &recvCount,
                        &neighborIndices[0],
                        icomm.mpiRecvBufferStatus() );

  for( int recvIdx = 0; recvIdx < recvCount; ++recvIdx )
  {
    NeighborCommunicator & neighbor = neighbors[ neighborIndices[ recvIdx ] ];
    neighbor.unpackBufferForSync( manager, icomm.commID(), onDevice, events );
  }

  // we don't want to check if the request has completed,
  //  we want to check that we've processed the resulting buffer
  //  which means that we've tested the request and it has been
  //  deallocated and set to MPI_REQUEST_NULL
  int allDone = true;
  const MPI_Request * reqs = icomm.mpiRecvBufferRequest( );
  for( int idx = 0; idx < icomm.size(); ++idx )
  {
    if( reqs[ idx ] != MPI_REQUEST_NULL )
    {
      allDone = false;
      break;
    }
  }

  return allDone;
}

bool CommunicationTools::asyncUnpack( MeshLevel & mesh,
                                      std::vector< NeighborCommunicator > & neighbors,
                                      MPI_iCommData & icomm,
                                      bool onDevice,
                                      parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  int recvCount = 0;
  std::vector< int > neighborIndices;
  neighborIndices.reserve( icomm.size() );
  MpiWrapper::testSome( icomm.size(),
                        icomm.mpiRecvBufferRequest(),
                        &recvCount,
                        &neighborIndices[0],
                        icomm.mpiRecvBufferStatus() );

  for( int recvIdx = 0; recvIdx < recvCount; ++recvIdx )
  {
    NeighborCommunicator & neighbor = neighbors[ neighborIndices[ recvIdx ] ];
    neighbor.unpackBufferForSync( icomm.getFieldNames(), mesh, icomm.commID(), onDevice, events );
  }

  // we don't want to check if the request has completed,
  //  we want to check that we've processed the resulting buffer
  //  which means that we've tested the request and it has been
  //  deallocated and set to MPI_REQUEST_NULL
  int allDone = true;
  const MPI_Request * reqs = icomm.mpiRecvBufferRequest( );
  for( int idx = 0; idx < icomm.size(); ++idx )
  {
    if( reqs[ idx ] != MPI_REQUEST_NULL )
    {
      allDone = false;
      break;
    }
  }

  return allDone;
}

void CommunicationTools::finalizeUnpack( ObjectManagerBase & manager,
                                         std::vector< NeighborCommunicator > & neighbors,
                                         MPI_iCommData & icomm,
                                         bool onDevice,
                                         parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  // poll mpi for completion then wait 10 nanoseconds 6,000,000,000 times (60 sec timeout)
  GEOSX_ASYNC_WAIT( 6000000000, 10, asyncUnpack( manager, neighbors, icomm, onDevice, events ) );
  if( onDevice )
  {
    waitAllDeviceEvents( events );
  }

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferSizeRequest(),
                       icomm.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferRequest(),
                       icomm.mpiSendBufferStatus() );
}

void CommunicationTools::finalizeUnpack( MeshLevel & mesh,
                                         std::vector< NeighborCommunicator > & neighbors,
                                         MPI_iCommData & icomm,
                                         bool onDevice,
                                         parallelDeviceEvents & events )
{
  GEOSX_MARK_FUNCTION;

  // poll mpi for completion then wait 10 nanoseconds 6,000,000,000 times (60 sec timeout)
  GEOSX_ASYNC_WAIT( 6000000000, 10, asyncUnpack( mesh, neighbors, icomm, onDevice, events ) );
  if( onDevice )
  {
    waitAllDeviceEvents( events );
  }

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferSizeRequest(),
                       icomm.mpiSendBufferSizeStatus() );

  MpiWrapper::waitAll( icomm.size(),
                       icomm.mpiSendBufferRequest(),
                       icomm.mpiSendBufferStatus() );

}

void CommunicationTools::synchronizeUnpack( ObjectManagerBase & manager,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            MPI_iCommData & icomm,
                                            bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  parallelDeviceEvents events;
  finalizeUnpack( manager, neighbors, icomm, onDevice, events );
}

void CommunicationTools::synchronizeUnpack( MeshLevel & mesh,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            MPI_iCommData & icomm,
                                            bool onDevice )
{
  GEOSX_MARK_FUNCTION;
  parallelDeviceEvents events;
  finalizeUnpack( mesh, neighbors, icomm, onDevice, events );
}

void CommunicationTools::synchronizeFields( string_array const & fieldNames,
                                            ObjectManagerBase & manager,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            bool onDevice )
{
  MPI_iCommData icomm( getCommID() );
  icomm.resize( neighbors.size() );
  synchronizePackSendRecvSizes( fieldNames, manager, neighbors, icomm, onDevice );
  synchronizePackSendRecv( fieldNames, manager, neighbors, icomm, onDevice );
  synchronizeUnpack( manager, neighbors, icomm, onDevice );
}

void CommunicationTools::synchronizeFields( const std::map< string, string_array > & fieldNames,
                                            MeshLevel & mesh,
                                            std::vector< NeighborCommunicator > & neighbors,
                                            bool onDevice )
{
  MPI_iCommData icomm( getCommID() );
  icomm.resize( neighbors.size() );
  synchronizePackSendRecvSizes( fieldNames, mesh, neighbors, icomm, onDevice );
  synchronizePackSendRecv( fieldNames, mesh, neighbors, icomm, onDevice );
  synchronizeUnpack( mesh, neighbors, icomm, onDevice );
}


} /* namespace geosx */
