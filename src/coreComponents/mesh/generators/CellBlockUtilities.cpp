/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "CellBlockUtilities.hpp"

#include "CellBlockManagerABC.hpp"

#include "common/GEOS_RAJA_Interface.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

void getFaceNodes( string const & elementType,
                   localIndex const iElement,
                   localIndex const iFace,
                   array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elementToNodes,
                   array1d< localIndex > & nodeIndices )
{
  if( elementType == "C3D8" )
  {
    nodeIndices.resize( 4 );
    if( iFace == 0 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
    else if( iFace == 1 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
      nodeIndices[3] = elementToNodes[iElement][1];
    }
    else if( iFace == 2 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][4];
      nodeIndices[2] = elementToNodes[iElement][6];
      nodeIndices[3] = elementToNodes[iElement][2];

    }
    else if( iFace == 3 )
    {
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][7];
      nodeIndices[3] = elementToNodes[iElement][5];
    }
    else if( iFace == 4 )
    {
      nodeIndices[0] = elementToNodes[iElement][3];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][6];
      nodeIndices[3] = elementToNodes[iElement][7];
    }
    else if( iFace == 5 )
    {
      nodeIndices[0] = elementToNodes[iElement][4];
      nodeIndices[1] = elementToNodes[iElement][5];
      nodeIndices[2] = elementToNodes[iElement][7];
      nodeIndices[3] = elementToNodes[iElement][6];
    }
  }
  else if( elementType == "C3D6" )
  {
    if( iFace == 0 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
    else if( iFace == 1 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
      nodeIndices[3] = elementToNodes[iElement][1];
    }
    else if( iFace == 2 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 3 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][5];
    }
    else if( iFace == 4 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][2];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][5];
      nodeIndices[3] = elementToNodes[iElement][4];
    }
  }
  else if( elementType == "C3D4" )
  {
    nodeIndices.resize( 3 );
    if( iFace == 0 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][1];
    }
    else if( iFace == 1 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][3];
    }
    else if( iFace == 2 )
    {
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][2];
    }
    else if( iFace == 3 )
    {
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][3];
    }
  }
  else if( elementType == "C3D5" )
  {
    if( iFace == 0 )
    {
      nodeIndices.resize( 4 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][2];
      nodeIndices[3] = elementToNodes[iElement][3];
    }
    else if( iFace == 1 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][0];
      nodeIndices[1] = elementToNodes[iElement][1];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 2 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][1];
      nodeIndices[1] = elementToNodes[iElement][2];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 3 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][2];
      nodeIndices[1] = elementToNodes[iElement][3];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
    else if( iFace == 4 )
    {
      nodeIndices.resize( 3 );
      nodeIndices[0] = elementToNodes[iElement][3];
      nodeIndices[1] = elementToNodes[iElement][0];
      nodeIndices[2] = elementToNodes[iElement][4];
    }
  }
  else
  {
    GEOSX_ERROR( "Error. Don't know what kind of element this is and cannot build faces." );
  }
}

/**
 * @class EdgeBuilder
 * @brief This class stores the data necessary to construct the various edge maps.
 */
struct EdgeBuilder
{
  /**
   * @brief Constructor.
   * @param [in] n1_ the greater of the two node indices that comprise the edge.
   * @param [in] faceID_ the ID of the face this edge came from.
   * @param [in] faceLocalEdgeIndex_ the face local index of this edge.
   */
  EdgeBuilder( localIndex const n1_,
               localIndex const faceID_,
               localIndex const faceLocalEdgeIndex_ ):
    n1( int32_t( n1_ ) ),
    faceID( int32_t( faceID_ ) ),
    faceLocalEdgeIndex( int32_t( faceLocalEdgeIndex_ ) )
  {}

  /**
   * @brief Return true if the two EdgeBuilders share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator==( EdgeBuilder const & rhs ) const
  { return n1 == rhs.n1; }

  /**
   * @brief Return true if the two EdgeBuilders don't share the same greatest node index.
   * @param [in] rhs the EdgeBuilder to compare against.
   */
  bool operator!=( EdgeBuilder const & rhs ) const
  { return !this->operator==( rhs ); }

  /// The larger of the two node indices that comprise the edge.
  int32_t n1;
  /// The face the edge came from.
  int32_t faceID;
  /// The face local index of the edge.
  int32_t faceLocalEdgeIndex;
};

/**
 * @brief Add an edge to the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 * @param [in] edgeID the ID of the edge to add.
 * @param [in] firstNodeID the ID of the first node of the edge.
 * @param [in] firstMatch the index of the first EdgeBuilder that describes this edge in edgesByLowestNode[ firstNodeID ].
 * @param [in] numMatches the number of EdgeBuilders that describe this edge in edgesByLowestNode[ firstNodeID ].
 */
void addEdge( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
              ArrayOfArraysView< localIndex > const & faceToEdgeMap,
              ArrayOfSetsView< localIndex > const & edgeToFaceMap,
              arrayView2d< localIndex > const & edgeToNodeMap,
              localIndex const edgeID,
              localIndex const firstNodeID,
              localIndex const firstMatch,
              localIndex const numMatches )
{
  GEOSX_ASSERT_GE( edgeToFaceMap.capacityOfSet( edgeID ), numMatches );

  // Populate the edge to node map.
  edgeToNodeMap( edgeID, 0 ) = firstNodeID;
  edgeToNodeMap( edgeID, 1 ) = edgesByLowestNode( firstNodeID, firstMatch ).n1;

  // Loop through all the matches and fill in the face to edge and edge to face maps.
  for( localIndex i = 0; i < numMatches; ++i )
  {
    localIndex const faceID = edgesByLowestNode( firstNodeID, firstMatch + i ).faceID;
    localIndex const faceLocalEdgeIndex = edgesByLowestNode( firstNodeID, firstMatch + i ).faceLocalEdgeIndex;

    faceToEdgeMap( faceID, faceLocalEdgeIndex ) = edgeID;
    edgeToFaceMap.insertIntoSet( edgeID, faceID );
  }
}

/**
 * @brief Populate the edgesByLowestNode map.
 * @param [in] numNodes The number of nodes.
 * @param [in] faceToNodeMap a map that associates an ordered list of nodes with each face.
 * @param [in/out] edgesByLowestNode of size numNodes, where each sub array has been preallocated to hold
 *        *enough* space.
 * For each edge of each face, this function gets the lowest node in the edge n0, creates an EdgeBuilder
 * associated with the edge and then appends the EdgeBuilder to edgesByLowestNode[ n0 ]. Finally it sorts
 * the contents of each sub-array of edgesByLowestNode from least to greatest.
 */
ArrayOfArrays< EdgeBuilder > createEdgesByLowestNode( localIndex numNodes,
                                                      ArrayOfArraysView< localIndex const > const & faceToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  ArrayOfArrays< EdgeBuilder > edgesByLowestNode( numNodes, 2 * CellBlockManagerABC::maxEdgesPerNode() );

  localIndex const numFaces = faceToNodeMap.size();

  // loop over all the faces.
  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const faceID )
  {
    localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceID );

    // loop over all the nodes in the face. there will be an edge for each node.
    for( localIndex a=0; a< numNodesInFace; ++a )
    {
      // sort the nodes in order of index value.
      localIndex node0 = faceToNodeMap( faceID, a );
      localIndex node1 = faceToNodeMap( faceID, ( a + 1 ) % numNodesInFace );
      if( node0 > node1 )
        std::swap( node0, node1 );

      // And append the edge to edgesByLowestNode.
      edgesByLowestNode.emplaceBackAtomic< parallelHostAtomic >( node0, node1, faceID, a );
    }
  } );

  // This comparator is not attached to EdgeBuilder because of potential inconsistencies with operator==
  auto const comp = []( EdgeBuilder const & e0, EdgeBuilder const & e1 ) -> bool
  {
    if( e0.n1 < e1.n1 )
    {
      return true;
    }
    if( e0.n1 > e1.n1 )
    {
      return false;
    }
    return e0.faceID < e1.faceID;
  };

  // Loop over all the nodes and sort the associated edges.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    EdgeBuilder * const edges = edgesByLowestNode[nodeID];
    std::sort( edges, edges + edgesByLowestNode.sizeOfArray( nodeID ), comp );
  } );

  return edgesByLowestNode;
}

/**
 * @brief Return the total number of unique edges and fill in the uniqueEdgeOffsets array.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [out] uniqueEdgeOffsets an array of size numNodes + 1. After this function returns node i contains
 * edges with IDs ranging from uniqueEdgeOffsets[ i ] to uniqueEdgeOffsets[ i + 1 ] - 1.
 */
localIndex calculateTotalNumberOfEdges( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                                        arrayView1d< localIndex > const & uniqueEdgeOffsets )
{
  localIndex const numNodes = edgesByLowestNode.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );

  uniqueEdgeOffsets[0] = 0;

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // If there are no edges associated with this node we can skip it.
    if( numEdges == 0 )
      return;

    localIndex & numUniqueEdges = uniqueEdgeOffsets[ nodeID + 1 ];
    numUniqueEdges = 0;

    // Otherwise since edgesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries.
    localIndex j = 0;
    for(; j < numEdges - 1; ++j )
    {
      numUniqueEdges += edgesByLowestNode( nodeID, j ) != edgesByLowestNode( nodeID, j + 1 );
    }

    numUniqueEdges += j == numEdges - 1;
  } );

  // At this point uniqueEdgeOffsets[ i ] holds the number of unique edges associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique edge offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueEdgeOffsets.begin(), uniqueEdgeOffsets.end() );

  return uniqueEdgeOffsets.back();
}


/**
 * @brief Populate the face to edge map, edge to face map, and edge to node map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an array containing the unique ID of the first edge associated with each node.
 * @param [in] faceToNodeMap the map from faces to nodes.
 * @param [in/out] faceToEdgeMap the map from face IDs to edge IDs.
 * @param [in/out] edgeToFacemap the map from edgeIDs to faceIDs.
 * @param [in/out] edgeToNodeMap the map from edgeIDs to nodeIDs.
 */
void populateEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                       arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                       ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                       ArrayOfArrays< localIndex > & faceToEdgeMap,
                       ArrayOfSets< localIndex > & edgeToFaceMap,
                       arrayView2d< localIndex > const & edgeToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numFaces = faceToNodeMap.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueEdgeOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numFaces, faceToEdgeMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToFaceMap.size() );
  GEOSX_ERROR_IF_NE( numUniqueEdges, edgeToNodeMap.size( 0 ) );

  // The face to edge map has the same shape as the face to node map, so we can resize appropriately.
  localIndex totalSize = 0;
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    totalSize += faceToNodeMap.sizeOfArray( faceID );
  }

  // Resize the face to edge map
  faceToEdgeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numFaces;
  faceToEdgeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face edges + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalSize + numFaces * CellBlockManagerABC::edgeMapExtraSpacePerFace() * ( 1 + 2 * overAllocationFactor );
  faceToEdgeMap.reserveValues( valuesToReserve );
  for( localIndex faceID = 0; faceID < numFaces; ++faceID )
  {
    faceToEdgeMap.appendArray( faceToNodeMap.sizeOfArray( faceID ) );
    faceToEdgeMap.setCapacityOfArray( faceToEdgeMap.size() - 1,
                                      faceToNodeMap.sizeOfArray( faceID ) + CellBlockManagerABC::edgeMapExtraSpacePerFace() );
  }

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }
      // Then add the edge.
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, numMatches );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      addEdge( edgesByLowestNode, faceToEdgeMap.toView(), edgeToFaceMap.toView(), edgeToNodeMap, curEdgeID, nodeID, j, 1 );
    }
  } );
}

/**
 * @brief Resize the edge to face map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * param [out] edgeToFaceMap the map from edges to faces. This function resizes the array appropriately.
 */
void resizeEdgeToFaceMap( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                          arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                          ArrayOfSets< localIndex > & edgeToFaceMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = edgesByLowestNode.size();
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  array1d< localIndex > numFacesPerEdge( numUniqueEdges );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalEdgeFaces( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID )
  {
    localIndex curEdgeID = uniqueEdgeOffsets[ nodeID ];
    localIndex const numEdges = edgesByLowestNode.sizeOfArray( nodeID );

    // loop over all the EdgeBuilders associated with the node
    localIndex j = 0;
    while( j < numEdges - 1 )
    {
      // Find the number of EdgeBuilders that describe the same edge
      localIndex numMatches = 1;
      while( edgesByLowestNode( nodeID, j ) == edgesByLowestNode( nodeID, j + numMatches ) )
      {
        ++numMatches;
        if( j + numMatches == numEdges )
          break;
      }

      // The number of matches is the number of faces associated with this edge.
      numFacesPerEdge( curEdgeID ) = numMatches;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
      ++curEdgeID;
      j += numMatches;
    }

    if( j == numEdges - 1 )
    {
      numFacesPerEdge( curEdgeID ) = 1;
      totalEdgeFaces += numFacesPerEdge( curEdgeID );
    }
  } );

  // Resize the edge to face map
  edgeToFaceMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueEdges;
  edgeToFaceMap.reserve( entriesToReserve );

  // Reserve space for the total number of edge faces + extra space for existing edges + even more space for new edges.
  localIndex const valuesToReserve = totalEdgeFaces.get() + numUniqueEdges * CellBlockManagerABC::faceMapExtraSpacePerEdge() * ( 1 + 2 * overAllocationFactor );
  edgeToFaceMap.reserveValues( valuesToReserve );

  // Append the individual sets.
  for( localIndex faceID = 0; faceID < numUniqueEdges; ++faceID )
  {
    edgeToFaceMap.appendSet( numFacesPerEdge[ faceID ] + CellBlockManagerABC::faceMapExtraSpacePerEdge() );
  }
}

/**
 * @brief Resize the edge to face map.
 * @param [in] edgesByLowestNode and array of size numNodes of arrays of EdgeBuilders associated with each node.
 * @param [in] uniqueEdgeOffsets an containing the unique edge IDs for each node in edgesByLowestNode.
 * param [out] edgeToNodesMap the map from edges to nodes. This function resizes the array appropriately.
 * param [out] edgeToFaceMap the map from edges to faces. This function resizes the array appropriately.
 */
void resizeEdgeMaps( ArrayOfArraysView< EdgeBuilder const > const & edgesByLowestNode,
                     arrayView1d< localIndex const > const & uniqueEdgeOffsets,
                     array2d< localIndex > & edgeToNodesMap,
                     ArrayOfSets< localIndex > & edgeToFaceMap )
{
  GEOSX_MARK_FUNCTION;

  resizeEdgeToFaceMap( edgesByLowestNode, uniqueEdgeOffsets, edgeToFaceMap );

  // Each face may belong to _maximum_ 2 elements.
  // If it belongs to only one, we put `-1` for the undefined value.
  localIndex const numUniqueEdges = uniqueEdgeOffsets.back();
  edgeToNodesMap.resize( numUniqueEdges, 2 );
}

localIndex buildEdgeMaps( localIndex numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfSets< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap )
{
  ArrayOfArrays< EdgeBuilder > edgesByLowestNode = createEdgesByLowestNode( numNodes, faceToNodeMap );

  array1d< localIndex > uniqueEdgeOffsets( numNodes + 1 );
  localIndex const numEdges = calculateTotalNumberOfEdges( edgesByLowestNode.toViewConst(), uniqueEdgeOffsets );

  resizeEdgeMaps( edgesByLowestNode.toViewConst(),
                  uniqueEdgeOffsets,
                  edgeToNodeMap,
                  edgeToFaceMap );

  populateEdgeMaps( edgesByLowestNode.toViewConst(),
                    uniqueEdgeOffsets,
                    faceToNodeMap,
                    faceToEdgeMap,
                    edgeToFaceMap,
                    edgeToNodeMap );

  return numEdges;
}

}
