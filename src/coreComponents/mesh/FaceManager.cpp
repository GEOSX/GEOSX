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
 * @file FaceManager.cpp
 */

#include "mesh/ExtrinsicMeshData.hpp"
#include "FaceManager.hpp"
#include "NodeManager.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "common/Logger.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

FaceManager::FaceManager( string const &, Group * const parent ) :
  ObjectManagerBase( "FaceManager", parent )
{
  this->registerWrapper( viewKeyStruct::nodeListString, &m_nodeList );
  this->registerWrapper( viewKeyStruct::edgeListString, &m_edgeList );
  //  m_nodeList.SetRelatedObject( parent->getGroup<NodeManager>(MeshLevel::groupStructKeys::nodeManagerString));

  this
    ->registerWrapper( viewKeyStruct::elementRegionListString,
                       &m_toElements.m_toElementRegion )
    ->setApplyDefaultValue( -1 );

  this
    ->registerWrapper( viewKeyStruct::elementSubRegionListString,
                       &m_toElements.m_toElementSubRegion )
    ->setApplyDefaultValue( -1 );

  this
    ->registerWrapper( viewKeyStruct::elementListString,
                       &m_toElements.m_toElementIndex )
    ->setApplyDefaultValue( -1 );

  this->registerWrapper( viewKeyStruct::faceAreaString, &m_faceArea );

  this->registerWrapper( viewKeyStruct::faceCenterString, &m_faceCenter )
    ->reference()
    .resizeDimension< 1 >( 3 );

  this->registerWrapper( viewKeyStruct::faceNormalString, &m_faceNormal )
    ->reference()
    .resizeDimension< 1 >( 3 );

  this
    ->registerWrapper( viewKeyStruct::faceRotationMatrixString,
                       &m_faceRotationMatrix )
    ->reference()
    .resizeDimension< 1, 2 >( 3, 3 );

  m_toElements.resize( 0, 2 );

  //0-based; note that the following field is ALSO 0
  //for faces that are not external faces, so check isExternal before using
  //  this->AddKeylessDataField<localIndex>("externalFaceIndex", true, true);
  //
  //  this->AddKeylessDataField<R1Tensor>("FaceCenter",true,true);
}

FaceManager::~FaceManager()
{}

void
FaceManager::resize( localIndex const newSize )
{
  m_nodeList.resize( newSize, 2 * nodeMapExtraSpacePerFace() );
  m_edgeList.resize( newSize, 2 * edgeMapExtraSpacePerFace() );
  ObjectManagerBase::resize( newSize );
}

/**
 * @class FaceBuilder
 * @brief This class stores the data necessary to construct the various face maps.
 */
struct FaceBuilder
{
  /**
   * @brief Constructor.
   * @param [in] n1_ the second smallest node index that comprise the face.
   * @param [in] n2_ the third smallest node index that comprise the face.
   * @param [in] er_ the element region this face came from.
   * @param [in] esr_ the element sub region this face came from.
   * @param [in] k_ the element this face came from.
   * @param [in] elementLocalFaceIndex_ the element local index of this face.
   */
  FaceBuilder( localIndex const n1_,
               localIndex const n2_,
               localIndex const er_,
               localIndex const esr_,
               localIndex const k_,
               localIndex const elementLocalFaceIndex_ ) :
    n1( int32_t( n1_ ) ),
    n2( int32_t( n2_ ) ),
    er( int32_t( er_ ) ),
    esr( int32_t( esr_ ) ),
    k( int32_t( k_ ) ),
    elementLocalFaceIndex( int32_t( elementLocalFaceIndex_ ) )
  {}

  /**
   * @brief Imposes an ordering on FaceBuilders. First compares n1, then n2, then er, esr, and k.
   * @param [in] rhs the FaceBuilder to compare against.
   * @return true if argument faceBuilder has lower second or third nodes or smaller region, subregion or element
   */
  bool
  operator<( FaceBuilder const & rhs ) const
  {
    if( n1 < rhs.n1 ) return true;
    if( n1 > rhs.n1 ) return false;
    if( n2 < rhs.n2 ) return true;
    if( n2 > rhs.n2 ) return false;
    if( er < rhs.er ) return true;
    if( er > rhs.er ) return false;
    if( esr < rhs.esr ) return true;
    if( esr > rhs.esr ) return false;
    return k < rhs.k;
  }

  /**
   * @brief Return true if the two FaceBuilders share the same second and third smallest nodes.
   * @param [in] rhs the FaceBuilder to compare against.
   * @return true if face second and third smallest node that comprise the face are the same
   */
  bool
  operator==( FaceBuilder const & rhs ) const
  {
    return n1 == rhs.n1 && n2 == rhs.n2;
  }

  /// index of the second smallest node of the face
  int32_t n1;
  /// index of the thirs smallest node of the face
  int32_t n2;
  /// element region label
  int32_t er;
  /// element subregion label
  int32_t esr;
  /// element label the face belongs to index
  int32_t k;
  /// the element local index of this face.
  int32_t elementLocalFaceIndex;
};

/**
 * @brief Get the three smallest values.
 * @param [in] values the array to search.
 * @param [out] minValues an array that will hold the three smallest values, from least to greatest.
 */
void
findSmallestThreeValues( arrayView1d< localIndex const > const & values,
                         localIndex ( &minValues )[3] )
{
  localIndex const n = values.size();
  GEOSX_ASSERT_GE( n, 3 );

  // Pick out the first three values
  minValues[0] = values[0];
  minValues[1] = values[1];
  minValues[2] = values[2];

  // and sort them
  if( minValues[0] > minValues[1] ) std::swap( minValues[0], minValues[1] );
  if( minValues[1] > minValues[2] ) std::swap( minValues[1], minValues[2] );
  if( minValues[0] > minValues[1] ) std::swap( minValues[0], minValues[1] );

  for( localIndex i = 3; i < n; ++i )
  {
    // If the current value is greater than the last minimum values skip it.
    if( values[i] > minValues[2] ) continue;

    minValues[2] = values[i];
    if( minValues[1] > minValues[2] ) std::swap( minValues[1], minValues[2] );
    if( minValues[0] > minValues[1] ) std::swap( minValues[0], minValues[1] );
  }
}

/**
 * @brief Populate the facesByLowestNode map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [inout] facesByLowestNode of size numNodes, where each sub array has been preallocated to hold
 *        *enough* space.
 * For each face of each element, this function gets the three lowest nodes in the face {n0, n1, n2}, creates
 * an EdgeBuilder associated with the face from n1 and n2 and then appends the EdgeBuilder to facesByLowestNode[ n0 ].
 * Finally it sorts the contents of each sub-array of facesByLowestNode from least to greatest.
 */
void
createFacesByLowestNode( ElementRegionManager const & elementManager,
                         ArrayOfArraysView< FaceBuilder > const & facesByLowestNode )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();

  // loop over all the regions
  for( typename dataRepository::indexType er = 0; er < elementManager.numRegions();
       ++er )
  {
    ElementRegionBase const & elemRegion = *elementManager.GetRegion( er );

    // loop over all the subregions
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >(
      [&]( localIndex const esr, CellElementSubRegion const & subRegion ) {
        localIndex const numFacesPerElement = subRegion.numFacesPerElement();
        localIndex const numElements = subRegion.size();

        // Begin the parallel region so that tempNodeList and lowestNodes are thread private.
        PRAGMA_OMP( "omp parallel" )
        {
          localIndex_array tempNodeList;
          localIndex lowestNodes[3];

          // Loop over all the elements.
          PRAGMA_OMP( "omp for" )
          for( localIndex k = 0; k < numElements; ++k )
          {
            for( localIndex elementLocalFaceIndex = 0;
                 elementLocalFaceIndex < numFacesPerElement;
                 ++elementLocalFaceIndex )
            {
              subRegion.GetFaceNodes( k, elementLocalFaceIndex, tempNodeList );
              findSmallestThreeValues( tempNodeList, lowestNodes );

              facesByLowestNode.emplaceBackAtomic< parallelHostAtomic >(
                lowestNodes[0],
                lowestNodes[1],
                lowestNodes[2],
                er,
                esr,
                k,
                elementLocalFaceIndex );
            }
          }
        }
      } );
  }

  // Loop over all the nodes and sort the associated faces.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    FaceBuilder * const faces = facesByLowestNode[nodeID];
    std::sort( faces, faces + facesByLowestNode.sizeOfArray( nodeID ) );
  } );
}

/**
 * @brief Return the total number of unique faces and fill in the uniqueFaceOffsets array.
 * @param [in] facesByLowestNode and array of size numNodes of arrays of FaceBuilders associated with each node.
 * @param [out] uniqueFaceOffsets an array of size numNodes + 1. After this function returns node i contains
 *              faces with IDs ranging from uniqueFaceOffsets[ i ] to uniqueFaceOffsets[ i + 1 ] - 1.
 * @return return total number of faces
 */
localIndex
calculateTotalNumberOfFaces(
  ArrayOfArraysView< FaceBuilder const > const & facesByLowestNode,
  arrayView1d< localIndex > const & uniqueFaceOffsets )
{
  localIndex const numNodes = facesByLowestNode.size();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );

  uniqueFaceOffsets[0] = 0;

  // Loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );

    // If there are no faces associated with this node we can skip it.
    if( numFaces == 0 ) return;

    localIndex & numUniqueFaces = uniqueFaceOffsets[nodeID + 1];
    numUniqueFaces = 0;

    // Otherwise since facesByLowestNode[ nodeID ] is sorted we can compare subsequent entries
    // count up the unique entries. Since each face can appear at most twice if we find a match we
    // can skip the next entry as well.
    localIndex j = 0;
    for( ; j < numFaces - 1; ++j )
    {
      ++numUniqueFaces;
      j += facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 );
    }

    numUniqueFaces += j == numFaces - 1;
  } );

  // At this point uniqueFaceOffsets[ i ] holds the number of unique face associated with node i - 1.
  // Perform an inplace prefix-sum to get the unique face offset.
  RAJA::inclusive_scan_inplace< parallelHostPolicy >( uniqueFaceOffsets.begin(),
                                                      uniqueFaceOffsets.end() );
  return uniqueFaceOffsets.back();
}

/**
 * @brief Resize the face to node map.
 * @param [in] elementManager the ElementRegionManager.
 * @param [in] facesByLowestNode and array of size numNodes of arrays of FaceBuilders associated with each node.
 * @param [in] uniqueFaceOffsets an containing the unique face IDs for each node in facesByLowestNode.
 * @param [out] faceToNodeMap the map from faces to nodes. This function resizes the array appropriately.
 */
void
resizeFaceToNodeMap( ElementRegionManager const & elementManager,
                     ArrayOfArraysView< FaceBuilder const > const & facesByLowestNode,
                     arrayView1d< localIndex const > const & uniqueFaceOffsets,
                     ArrayOfArrays< localIndex > & faceToNodeMap )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  array1d< localIndex > numNodesPerFace( numUniqueFaces );
  RAJA::ReduceSum< parallelHostReduce, localIndex > totalFaceNodes( 0.0 );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    localIndex curFaceID = uniqueFaceOffsets[nodeID];
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );

    // If there are no faces associated with this node we can skip it.
    if( numFaces == 0 ) return;

    // loop over all the FaceBuilders associated with the node
    localIndex j = 0;
    for( ; j < numFaces - 1; ++j )
    {
      localIndex const er = facesByLowestNode( nodeID, j ).er;
      localIndex const esr = facesByLowestNode( nodeID, j ).esr;
      localIndex const k = facesByLowestNode( nodeID, j ).k;
      localIndex const elementLocalFaceIndex =
        facesByLowestNode( nodeID, j ).elementLocalFaceIndex;

      // Get the number of face nodes from the subregion.
      CellElementSubRegion const & subRegion =
        *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
      numNodesPerFace[curFaceID] =
        subRegion.GetNumFaceNodes( k, elementLocalFaceIndex );
      totalFaceNodes += numNodesPerFace[curFaceID];

      j += facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 );
      ++curFaceID;
    }

    if( j == numFaces - 1 )
    {
      localIndex const er = facesByLowestNode( nodeID, j ).er;
      localIndex const esr = facesByLowestNode( nodeID, j ).esr;
      localIndex const k = facesByLowestNode( nodeID, j ).k;
      localIndex const elementLocalFaceIndex =
        facesByLowestNode( nodeID, j ).elementLocalFaceIndex;

      // Get the number of face nodes from the subregion.
      CellElementSubRegion const & subRegion =
        *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
      numNodesPerFace[curFaceID] =
        subRegion.GetNumFaceNodes( k, elementLocalFaceIndex );
      totalFaceNodes += numNodesPerFace[curFaceID];
    }
  } );

  // Resize the face to node map.
  faceToNodeMap.resize( 0 );

  // Reserve space for the number of current faces plus some extra.
  double const overAllocationFactor = 0.3;
  localIndex const entriesToReserve = ( 1 + overAllocationFactor ) * numUniqueFaces;
  faceToNodeMap.reserve( entriesToReserve );

  // Reserve space for the total number of face nodes + extra space for existing faces + even more space for new faces.
  localIndex const valuesToReserve = totalFaceNodes.get() +
    numUniqueFaces * FaceManager::nodeMapExtraSpacePerFace() *
      ( 1 + 2 * overAllocationFactor );
  faceToNodeMap.reserveValues( valuesToReserve );

  // Append the individual arrays.
  for( localIndex faceID = 0; faceID < numUniqueFaces; ++faceID )
  {
    faceToNodeMap.appendArray( numNodesPerFace[faceID] );
    faceToNodeMap.setCapacityOfArray(
      faceToNodeMap.size() - 1,
      numNodesPerFace[faceID] + FaceManager::nodeMapExtraSpacePerFace() );
  }
}

/**
 * @brief Add an interior face to the element lists, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] faceID the ID of the face to add.
 * @param [in] fb0 the FaceBuilder associated with the first element of the current face.
 * @param [in] fb1 the FaceBuilder associated with the second element of the current face.
 * @param [inout] elemRegionList the face to element region map.
 * @param [inout] elemSubRegionList the face to element subregion map.
 * @param [inout] elemList the face to element map.
 * @param [inout] nodeList the face to node map.
 */
void
addInteriorFace( ElementRegionManager & elementManager,
                 localIndex const faceID,
                 FaceBuilder const & fb0,
                 FaceBuilder const & fb1,
                 arrayView2d< localIndex > const & elemRegionList,
                 arrayView2d< localIndex > const & elemSubRegionList,
                 arrayView2d< localIndex > const & elemList,
                 ArrayOfArrays< localIndex > & nodeList )
{
  // Handle the first element.
  {
    localIndex const er = fb0.er;
    localIndex const esr = fb0.esr;
    localIndex const k = fb0.k;
    localIndex const elementLocalFaceIndex = fb0.elementLocalFaceIndex;

    // Get the subRegion associated with the element.
    CellElementSubRegion & subRegion =
      *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );

    // The first element defines the node ordering for the face.
    localIndex const numFaceNodes =
      subRegion.GetFaceNodes( k, elementLocalFaceIndex, nodeList[faceID] );
    GEOSX_ASSERT_EQ( numFaceNodes, nodeList.sizeOfArray( faceID ) );
    GEOSX_DEBUG_VAR( numFaceNodes );

    // Add the face to the element to face map.
    subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

    // Populate the face to element maps.
    elemRegionList( faceID, 0 ) = er;
    elemSubRegionList( faceID, 0 ) = esr;
    elemList( faceID, 0 ) = k;
  }

  // Handle the second element.
  {
    localIndex const er = fb1.er;
    localIndex const esr = fb1.esr;
    localIndex const k = fb1.k;
    localIndex const elementLocalFaceIndex = fb1.elementLocalFaceIndex;

    CellElementSubRegion & subRegion =
      *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
    subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

    elemRegionList( faceID, 1 ) = er;
    elemSubRegionList( faceID, 1 ) = esr;
    elemList( faceID, 1 ) = k;
  }
}

/**
 * @brief Add a boundary face to the element lists, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] faceID the ID of the face to add.
 * @param [in] fb the FaceBuilder associated with the first element of the current face.
 * @param [inout] elemRegionList the face to element region map.
 * @param [inout] elemSubRegionList the face to element subregion map.
 * @param [inout] elemList the face to element map.
 * @param [inout] nodeList the face to node map.
 */
void
addBoundaryFace( ElementRegionManager & elementManager,
                 localIndex const faceID,
                 FaceBuilder const & fb,
                 arrayView2d< localIndex > const & elemRegionList,
                 arrayView2d< localIndex > const & elemSubRegionList,
                 arrayView2d< localIndex > const & elemList,
                 ArrayOfArrays< localIndex > & nodeList )
{
  localIndex const er = fb.er;
  localIndex const esr = fb.esr;
  localIndex const k = fb.k;
  localIndex const elementLocalFaceIndex = fb.elementLocalFaceIndex;

  // Get the subRegion associated with the element.
  CellElementSubRegion & subRegion =
    *elementManager.GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );

  // Get the nodes associated with the face.
  localIndex const numFaceNodes =
    subRegion.GetFaceNodes( k, elementLocalFaceIndex, nodeList[faceID] );
  GEOSX_ASSERT_EQ( numFaceNodes, nodeList.sizeOfArray( faceID ) );
  GEOSX_DEBUG_VAR( numFaceNodes );

  // Add the face to the element to face map.
  subRegion.faceList()( k, elementLocalFaceIndex ) = faceID;

  // Populate the face to element maps.
  elemRegionList( faceID, 0 ) = er;
  elemSubRegionList( faceID, 0 ) = esr;
  elemList( faceID, 0 ) = k;

  elemRegionList( faceID, 1 ) = -1;
  elemSubRegionList( faceID, 1 ) = -1;
  elemList( faceID, 1 ) = -1;
}

/**
 * @brief Populate the face to element maps, the face to node map, and the element to face map.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level.
 * @param [in] facesByLowestNode and array of size numNodes of arrays of FaceBuilders associated with each node.
 * @param [in] uniqueFaceOffsets an array containing the unique ID of the first face associated with each node.
 * @param [inout] elemRegionList the face to element region map.
 * @param [inout] elemSubRegionList the face to element subregion map.
 * @param [inout] elemList the face to element map.
 * @param [inout] nodeList the face to node map.
 */
void
populateMaps( ElementRegionManager & elementManager,
              ArrayOfArraysView< FaceBuilder const > const & facesByLowestNode,
              arrayView1d< localIndex const > const & uniqueFaceOffsets,
              arrayView2d< localIndex > const & elemRegionList,
              arrayView2d< localIndex > const & elemSubRegionList,
              arrayView2d< localIndex > const & elemList,
              ArrayOfArrays< localIndex > & nodeList )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numNodes = facesByLowestNode.size();
  localIndex const numUniqueFaces = uniqueFaceOffsets.back();
  GEOSX_ERROR_IF_NE( numNodes, uniqueFaceOffsets.size() - 1 );
  GEOSX_ERROR_IF_NE( numUniqueFaces, elemRegionList.size( 0 ) );
  GEOSX_ERROR_IF_NE( numUniqueFaces, elemSubRegionList.size( 0 ) );
  GEOSX_ERROR_IF_NE( numUniqueFaces, elemList.size( 0 ) );
  GEOSX_ERROR_IF_NE( numUniqueFaces, nodeList.size() );

  // loop over all the nodes.
  forAll< parallelHostPolicy >( numNodes, [&]( localIndex const nodeID ) {
    localIndex curFaceID = uniqueFaceOffsets[nodeID];
    localIndex const numFaces = facesByLowestNode.sizeOfArray( nodeID );

    // loop over all the FaceBuilders associated with the node
    localIndex j = 0;
    for( ; j < numFaces - 1; ++j )
    {
      // If two subsequent FaceBuilders compare equal then they describe an interior face.
      if( facesByLowestNode( nodeID, j ) == facesByLowestNode( nodeID, j + 1 ) )
      {
        addInteriorFace( elementManager,
                         curFaceID,
                         facesByLowestNode( nodeID, j ),
                         facesByLowestNode( nodeID, j + 1 ),
                         elemRegionList,
                         elemSubRegionList,
                         elemList,
                         nodeList );
        ++j;
      }
      // Otherwise it's a boundary face.
      else
      {
        addBoundaryFace( elementManager,
                         curFaceID,
                         facesByLowestNode( nodeID, j ),
                         elemRegionList,
                         elemSubRegionList,
                         elemList,
                         nodeList );
      }

      ++curFaceID;
    }

    if( j == numFaces - 1 )
    {
      addBoundaryFace( elementManager,
                       curFaceID,
                       facesByLowestNode( nodeID, j ),
                       elemRegionList,
                       elemSubRegionList,
                       elemList,
                       nodeList );
    }
  } );
}

void
FaceManager::BuildFaces( NodeManager * const nodeManager,
                         ElementRegionManager * const elementManager )
{
  GEOSX_MARK_FUNCTION;

  m_toElements.setElementRegionManager( elementManager );

  localIndex const numNodes = nodeManager->size();

  ArrayOfArrays< FaceBuilder > facesByLowestNode( numNodes, 2 * maxFacesPerNode() );
  createFacesByLowestNode( *elementManager, facesByLowestNode.toView() );

  array1d< localIndex > uniqueFaceOffsets( numNodes + 1 );
  localIndex const numFaces =
    calculateTotalNumberOfFaces( facesByLowestNode.toViewConst(),
                                 uniqueFaceOffsets );

  resizeFaceToNodeMap( *elementManager,
                       facesByLowestNode.toViewConst(),
                       uniqueFaceOffsets,
                       nodeList() );

  resize( numFaces );

  populateMaps( *elementManager,
                facesByLowestNode.toViewConst(),
                uniqueFaceOffsets,
                m_toElements.m_toElementRegion,
                m_toElements.m_toElementSubRegion,
                m_toElements.m_toElementIndex,
                nodeList() );

  // First create the sets
  auto const & nodeSets = nodeManager->sets().wrappers();
  for( localIndex i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    CreateSet( setName );
  }

  // Then loop over them in parallel and fill them in.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void {
    auto const & setWrapper = nodeSets[i];
    std::string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const & targetSet =
      nodeManager->sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    ConstructSetFromSetAndMap( targetSet, m_nodeList.toViewConst(), setName );
  } );

  SetDomainBoundaryObjects( nodeManager );

  computeGeometry( nodeManager );
}

void
FaceManager::computeGeometry( NodeManager const * const nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  forAll< parallelHostPolicy >( this->size(), [&]( localIndex const faceID ) {
    m_faceArea[faceID] =
      computationalGeometry::Centroid_3DPolygon( m_nodeList[faceID],
                                                 X,
                                                 m_faceCenter[faceID],
                                                 m_faceNormal[faceID] );

    // This needs to be done somewhere else, also we probably shouldn't be orienting the normals like this.
    // Set normal orientation according to a global criterion
    computationalGeometry::FixNormalOrientation_3D( m_faceNormal[faceID] );

    // Compute the local rotation matrix according to the normal vector
    computationalGeometry::RotationMatrix_3D( m_faceNormal[faceID],
                                              m_faceRotationMatrix[faceID] );
  } );
}

void
FaceManager::SetDomainBoundaryObjects( NodeManager * const nodeManager )
{
  // Set value of domainBounaryIndicator to one if it is found to have only one elements that it
  // is connected to.
  integer_array & faceDomainBoundaryIndicator =
    this->getReference< integer_array >( viewKeys.domainBoundaryIndicator );
  faceDomainBoundaryIndicator.setValues< serialPolicy >( 0 );

  arrayView2d< localIndex const > const & elemRegionList = this->elementRegionList();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const kf ) {
    if( elemRegionList[kf][1] == -1 )
    {
      faceDomainBoundaryIndicator( kf ) = 1;
    }
  } );

  integer_array & nodeDomainBoundaryIndicator =
    nodeManager->getReference< integer_array >(
      nodeManager->viewKeys.domainBoundaryIndicator );
  nodeDomainBoundaryIndicator.setValues< serialPolicy >( 0 );

  ArrayOfArraysView< localIndex const > const & faceToNodesMap =
    this->nodeList().toViewConst();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const k ) {
    if( faceDomainBoundaryIndicator[k] == 1 )
    {
      localIndex const numNodes = faceToNodesMap.sizeOfArray( k );
      for( localIndex a = 0; a < numNodes; ++a )
      {
        nodeDomainBoundaryIndicator[faceToNodesMap( k, a )] = 1;
      }
    }
  } );
}

void
FaceManager::SetIsExternal()
{
  integer_array const & isDomainBoundary =
    this->getReference< integer_array >( viewKeys.domainBoundaryIndicator );

  m_isExternal.setValues< serialPolicy >( 0 );
  for( localIndex k = 0; k < size(); ++k )
  {
    if( isDomainBoundary[k] == 1 )
    {
      m_isExternal[k] = 1;
    }
  }
}

localIndex
FaceManager::getMaxFaceNodes() const
{
  localIndex maxSize = 0;
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    nodeList().toViewConst();
  for( localIndex kf = 0; kf < size(); ++kf )
  {
    maxSize = std::max( maxSize, faceToNodeMap.sizeOfArray( kf ) );
  }

  return maxSize;
}

void
FaceManager::SortAllFaceNodes( NodeManager const * const nodeManager,
                               ElementRegionManager const * const elemManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const & elemRegionList = elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();

  const indexType max_face_nodes = getMaxFaceNodes();
  GEOSX_ERROR_IF( max_face_nodes >= MAX_FACE_NODES,
                  "More nodes on a face than expected!" );

  elemManager->forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion const & subRegion ) {
      subRegion.calculateElementCenters( X );
    } );

  ArrayOfArraysView< localIndex > const & faceToNodeMap = nodeList().toView();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const kf ) {
    ElementRegionBase const * const elemRegion =
      elemManager->GetRegion( elemRegionList[kf][0] );
    CellElementSubRegion const * const subRegion =
      elemRegion->GetSubRegion< CellElementSubRegion >( elemSubRegionList[kf][0] );
    const localIndex numFaceNodes = faceToNodeMap.sizeOfArray( kf );
    SortFaceNodes( X,
                   subRegion->getElementCenter()[elemList( kf, 0 )],
                   faceToNodeMap[kf],
                   numFaceNodes );
  } );
}

void
FaceManager::SortFaceNodes(
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
  arraySlice1d< real64 const > const elementCenter,
  localIndex * const faceNodes,
  localIndex const numFaceNodes )
{
  localIndex const firstNodeIndex = faceNodes[0];

  // get face center (average vertex location)
  R1Tensor fc( 0 );
  for( localIndex n = 0; n < numFaceNodes; ++n )
  {
    fc += X[faceNodes[n]];
  }
  fc /= realT( numFaceNodes );

  R1Tensor ex, ey, ez;
  // Approximate face normal direction (unscaled)

  if( numFaceNodes == 2 )  //2D only.
  {
    ex = X[faceNodes[1]];
    ex -= X[faceNodes[0]];
    ey = elementCenter;
    ey -= fc;

    ez.Cross( ex, ey );
    // The element should be on the right hand side of the vector from node 0 to
    // node 1.
    // This ensure that the normal vector of an external face points to outside
    // the element.
    if( ez[2] > 0 )
    {
      localIndex itemp = faceNodes[0];
      faceNodes[0] = faceNodes[1];
      faceNodes[1] = itemp;
    }
  }
  else
  {
    ez = fc;
    ez -= elementCenter;

    /// Approximate in-plane axis
    ex = X[faceNodes[0]];
    ex -= fc;

    ex /= ex.L2_Norm();
    ey.Cross( ez, ex );
    ey /= ey.L2_Norm();

    std::pair< realT, localIndex > thetaOrder[MAX_FACE_NODES];

    /// Sort nodes counterclockwise around face center
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      R1Tensor v = X[faceNodes[n]];
      v -= fc;
      thetaOrder[n] = std::pair< realT, localIndex >( atan2( Dot( v, ey ), Dot( v, ex ) ),
                                                      faceNodes[n] );
    }

    std::sort( thetaOrder, thetaOrder + numFaceNodes );

    // Reorder nodes on face
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex tempFaceNodes[MAX_FACE_NODES];

    localIndex firstIndexIndex = 0;
    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }

    for( localIndex n = 0; n < numFaceNodes; ++n )
    {
      const localIndex index = ( firstIndexIndex + n ) % numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }
  }
}

void
FaceManager::ExtractMapFromObjectForAssignGlobalIndexNumbers(
  ObjectManagerBase const * const nodeManager,
  std::vector< std::vector< globalIndex > > & globalFaceNodes )
{
  GEOSX_MARK_FUNCTION;
  nodeManager->CheckTypeID( typeid( NodeManager ) );

  localIndex const numFaces = size();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    this->nodeList().toViewConst();
  arrayView1d< integer const > const & isDomainBoundary =
    this->getReference< integer_array >( viewKeys.domainBoundaryIndicator );

  globalFaceNodes.resize( numFaces );

  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const & faceID ) {
    std::vector< globalIndex > & curFaceGlobalNodes = globalFaceNodes[faceID];

    if( isDomainBoundary( faceID ) )
    {
      localIndex const numNodes = faceToNodeMap.sizeOfArray( faceID );
      curFaceGlobalNodes.resize( numNodes );

      for( localIndex a = 0; a < numNodes; ++a )
      {
        curFaceGlobalNodes[a] =
          nodeManager->localToGlobalMap()( faceToNodeMap( faceID, a ) );
      }

      std::sort( curFaceGlobalNodes.begin(), curFaceGlobalNodes.end() );
    }
    else
    {
      curFaceGlobalNodes.resize( 0 );
    }
  } );
}

void
FaceManager::ViewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::elementRegionListString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::elementSubRegionListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString ) );
}

localIndex
FaceManager::PackUpDownMapsSize(
  arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

localIndex
FaceManager::PackUpDownMaps(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex
FaceManager::PackUpDownMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString ) );

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_nodeList.Base().toViewConst(),
                                           m_unmappedGlobalIndicesInToNodes,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_nodeList.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_edgeList.Base().toViewConst(),
                                           m_unmappedGlobalIndicesInToEdges,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_edgeList.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );

  return packedSize;
}

localIndex
FaceManager::UnpackUpDownMaps(
  buffer_unit_type const *& buffer,
  localIndex_array & packList,
  bool const overwriteUpMaps,
  bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  // GEOSX_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_nodeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_nodeList.RelatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOSX_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_edgeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_edgeList.RelatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOSX_ERROR_IF_NE( elementListString, viewKeyStruct::elementListString );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void
FaceManager::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_nodeList,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_edgeList,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );
}

void
FaceManager::compressRelationMaps()
{
  m_nodeList.compress();
  m_edgeList.compress();
}

void
FaceManager::enforceStateFieldConsistencyPostTopologyChange(
  std::set< localIndex > const & targetIndices )
{
  arrayView1d< localIndex const > const & childFaceIndices =
    getExtrinsicData< extrinsicMeshData::ChildIndex >();

  ObjectManagerBase::enforceStateFieldConsistencyPostTopologyChange( targetIndices );

  for( localIndex const targetIndex : targetIndices )
  {
    localIndex const childIndex = childFaceIndices[targetIndex];
    if( childIndex != -1 )
    {
      LvArray::tensorOps::scaledCopy< 3 >( m_faceNormal[targetIndex],
                                           m_faceNormal[childIndex],
                                           -1 );
    }
  }
}

void
FaceManager::depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                               ElementRegionManager const & elemRegionManager )
{
  for( auto const & targetIndex : receivedFaces )
  {
    for( localIndex k = 0; k < m_toElements.m_toElementRegion.size( 1 ); ++k )
    {
      localIndex const elemRegionIndex =
        m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex =
        m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex = m_toElements.m_toElementIndex[targetIndex][k];

      if( elemRegionIndex != -1 && elemSubRegionIndex != -1 && elemIndex != -1 )
      {
        CellElementSubRegion const * subRegion =
          elemRegionManager.GetRegion( elemRegionIndex )
            ->GetSubRegion< CellElementSubRegion >( elemSubRegionIndex );
        array2d< localIndex > const & downmap = subRegion->faceList();
        bool hasTargetIndex = false;

        for( localIndex a = 0; a < downmap.size( 1 ); ++a )
        {
          localIndex const compositeLocalIndex = downmap[elemIndex][a];
          if( compositeLocalIndex == targetIndex )
          {
            hasTargetIndex = true;
          }
        }
        if( !hasTargetIndex )
        {
          m_toElements.m_toElementRegion[targetIndex][k] = -1;
          m_toElements.m_toElementSubRegion[targetIndex][k] = -1;
          m_toElements.m_toElementIndex[targetIndex][k] = -1;
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        FaceManager,
                        std::string const &,
                        Group * const )

}  // namespace geosx
