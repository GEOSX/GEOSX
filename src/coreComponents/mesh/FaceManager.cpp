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
 * @file FaceManager.cpp
 */

#include "mesh/ExtrinsicMeshData.hpp"
#include "FaceManager.hpp"
#include "NodeManager.hpp"
#include "BufferOps.hpp"
#include "common/TimingMacros.hpp"
#include "ElementRegionManager.hpp"
#include "utilities/ComputationalGeometry.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/Logger.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

FaceManager::FaceManager( string const &, Group * const parent ):
  ObjectManagerBase( "FaceManager", parent )
{
  this->registerWrapper( viewKeyStruct::nodeListString(), &m_nodeList );
  this->registerWrapper( viewKeyStruct::edgeListString(), &m_edgeList );

  this->registerWrapper( viewKeyStruct::elementRegionListString(), &m_toElements.m_toElementRegion ).
    setApplyDefaultValue( -1 );

  this->registerWrapper( viewKeyStruct::elementSubRegionListString(), &m_toElements.m_toElementSubRegion ).
    setApplyDefaultValue( -1 );

  // Do we really want this to be resized and accessed by anyone?
  this->registerWrapper( viewKeyStruct::elementListString(), &m_toElements.m_toElementIndex ).
    setApplyDefaultValue( -1 );

  this->registerWrapper( viewKeyStruct::faceAreaString(), &m_faceArea );

  this->registerWrapper( viewKeyStruct::faceCenterString(), &m_faceCenter ).
    reference().resizeDimension< 1 >( 3 );

  this->registerWrapper( viewKeyStruct::faceNormalString(), &m_faceNormal ).
    reference().resizeDimension< 1 >( 3 );

  m_toElements.resize( 0, 2 );

}

FaceManager::~FaceManager()
{}

void FaceManager::resize( localIndex const newSize )
{
  m_nodeList.resize( newSize, 2 * nodeMapExtraSpacePerFace() );
  m_edgeList.resize( newSize, 2 * edgeMapExtraSpacePerFace() );
  ObjectManagerBase::resize( newSize );
}

/**
 * @brief Populate the mappings from face to element region, and from face to element subregion.
 * @param [in] elementManager the ElementRegionManager associated with this mesh level. Regions and subregions come from it.
 * @param [out] elemRegionList the face to element region map.
 * @param [out] elemSubRegionList the face to element subregion map.
 *
 * @note @p elemRegionList and @p elemSubRegionList needs to have the correct dimensions. Values are all overwritten.
 */
void populateRegions( ElementRegionManager const & elementManager,
                      arrayView2d< localIndex > const & elemRegionList,
                      arrayView2d< localIndex > const & elemSubRegionList )
{
  GEOSX_ERROR_IF_NE( 2, elemRegionList.size( 1 ) );
  GEOSX_ERROR_IF_NE( 2, elemSubRegionList.size( 1 ) );

  elemRegionList.setValues< serialPolicy >( -1 );
  elemSubRegionList.setValues< serialPolicy >( -1 );

  for( typename dataRepository::indexType er = 0; er < elementManager.numRegions(); ++er )
  {
    ElementRegionBase const & elemRegion = elementManager.getRegion( er );

    // loop over all the subregions
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >(
      [&]( localIndex const esr, CellElementSubRegion const & subRegion )
    {
      FixedOneToManyRelation const & elementToFaces = subRegion.faceList();
      for( localIndex iElement = 0; iElement < subRegion.size(); ++iElement )
      {
        for( localIndex iFace = 0; iFace < subRegion.numFacesPerElement(); ++iFace )
        {
          localIndex const & faceId = elementToFaces( iElement, iFace );
          // Storing in the first free slot (i.e. not -1)
          if( elemSubRegionList( faceId, 0 ) == -1 )
          {
            elemSubRegionList( faceId, 0 ) = esr;
          }
          else
          {
            elemSubRegionList( faceId, 1 ) = esr;
          }
          if( elemRegionList( faceId, 0 ) == -1 )
          {
            elemRegionList( faceId, 0 ) = er;
          }
          else
          {
            elemRegionList( faceId, 1 ) = er;
          }
        }
      }
    } );
  }
}

void FaceManager::buildFaces( CellBlockManagerABC const & cellBlockManager,
                              NodeManager & nodeManager,
                              ElementRegionManager & elementManager )
{
  GEOSX_MARK_FUNCTION;

  m_toElements.setElementRegionManager( elementManager );

  resize( cellBlockManager.numFaces() );

  populateRegions( elementManager,
                   m_toElements.m_toElementRegion,
                   m_toElements.m_toElementSubRegion );

  m_toElements.m_toElementIndex = cellBlockManager.getFaceToElements(); // TODO do better
  nodeList().base() = cellBlockManager.getFaceToNodes(); // TODO

  // First create the sets
  auto const & nodeSets = nodeManager.sets().wrappers();
  for( localIndex i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    createSet( setName );
  }

  // Then loop over them in parallel and fill them in.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const & targetSet = nodeManager.sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    constructSetFromSetAndMap( targetSet, m_nodeList.toViewConst(), setName );
  } );

  setDomainBoundaryObjects( nodeManager );

  computeGeometry( nodeManager );
}


void FaceManager::computeGeometry( NodeManager const & nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  // loop over faces and calculate faceArea, faceNormal and faceCenter
  forAll< parallelHostPolicy >( this->size(), [&]( localIndex const faceID )
  {
    m_faceArea[ faceID ] = computationalGeometry::Centroid_3DPolygon( m_nodeList[ faceID ],
                                                                      X,
                                                                      m_faceCenter[ faceID ],
                                                                      m_faceNormal[ faceID ] );

  } );
}

void FaceManager::setDomainBoundaryObjects( NodeManager & nodeManager )
{
  // Set value of domainBoundaryIndicator to one if it is found to have only one elements that it
  // is connected to.
  arrayView1d< integer > const & faceDomainBoundaryIndicator = this->getDomainBoundaryIndicator();
  faceDomainBoundaryIndicator.zero();

  arrayView2d< localIndex const > const elemRegionList = this->elementRegionList();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const kf )
  {
    if( elemRegionList[kf][1] == -1 )
    {
      faceDomainBoundaryIndicator( kf ) = 1;
    }
  } );

  arrayView1d< integer > const & nodeDomainBoundaryIndicator = nodeManager.getDomainBoundaryIndicator();
  nodeDomainBoundaryIndicator.zero();

  ArrayOfArraysView< localIndex const > const faceToNodesMap = this->nodeList().toViewConst();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const k )
  {
    if( faceDomainBoundaryIndicator[k] == 1 )
    {
      localIndex const numNodes = faceToNodesMap.sizeOfArray( k );
      for( localIndex a=0; a< numNodes; ++a )
      {
        nodeDomainBoundaryIndicator[faceToNodesMap( k, a )] = 1;
      }
    }
  } );
}

void FaceManager::setIsExternal()
{
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();

  m_isExternal.zero();
  for( localIndex k=0; k<size(); ++k )
  {
    if( isDomainBoundary[k]==1 )
    {
      m_isExternal[k] = 1;
    }
  }
}

localIndex FaceManager::getMaxFaceNodes() const
{
  localIndex maxSize = 0;
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = nodeList().toViewConst();
  for( localIndex kf =0; kf < size(); ++kf )
  {
    maxSize = std::max( maxSize, faceToNodeMap.sizeOfArray( kf ) );
  }

  return maxSize;
}

void FaceManager::sortAllFaceNodes( NodeManager const & nodeManager,
                                    ElementRegionManager const & elemManager )
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< localIndex const > const elemRegionList = elementRegionList();
  arrayView2d< localIndex const > const elemSubRegionList = elementSubRegionList();
  arrayView2d< localIndex const > const elemList = elementList();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  const indexType max_face_nodes = getMaxFaceNodes();
  GEOSX_ERROR_IF( max_face_nodes >= MAX_FACE_NODES, "More nodes on a face than expected!" );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion )
  { subRegion.calculateElementCenters( X ); } );

  ArrayOfArraysView< localIndex > const & faceToNodeMap = nodeList().toView();

  forAll< parallelHostPolicy >( size(), [&]( localIndex const kf )
  {
    ElementRegionBase const & elemRegion = elemManager.getRegion( elemRegionList[kf][0] );
    CellElementSubRegion const & subRegion = elemRegion.getSubRegion< CellElementSubRegion >( elemSubRegionList[kf][0] );
    localIndex const numFaceNodes = faceToNodeMap.sizeOfArray( kf );
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();
    sortFaceNodes( X, elemCenter[ elemList( kf, 0 ) ], faceToNodeMap[ kf ], numFaceNodes );
  } );
}

void FaceManager::sortFaceNodes( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                 arraySlice1d< real64 const > const elementCenter,
                                 localIndex * const faceNodes,
                                 localIndex const numFaceNodes )
{
  localIndex const firstNodeIndex = faceNodes[0];

  // get face center (average vertex location)
  real64 fc[3] = { 0 };
  for( localIndex n =0; n < numFaceNodes; ++n )
  {
    LvArray::tensorOps::add< 3 >( fc, X[faceNodes[n]] );
  }
  LvArray::tensorOps::scale< 3 >( fc, 1.0 / numFaceNodes );

  //real64 ex[3], ey[3], ez[3];
  // Approximate face normal direction (unscaled)

  if( numFaceNodes == 2 )  //2D only.
  {
    real64 ex[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[1]] );
    LvArray::tensorOps::subtract< 3 >( ex, X[faceNodes[0]] );

    real64 ey[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( elementCenter );
    LvArray::tensorOps::subtract< 3 >( ey, fc );

    real64 ez[3];
    LvArray::tensorOps::crossProduct( ez, ex, ey );

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
    real64 ez[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
    LvArray::tensorOps::subtract< 3 >( ez, elementCenter );

    /// Approximate in-plane axis
    real64 ex[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[0]] );
    LvArray::tensorOps::subtract< 3 >( ex, fc );
    LvArray::tensorOps::normalize< 3 >( ex );

    real64 ey[3];
    LvArray::tensorOps::crossProduct( ey, ez, ex );
    LvArray::tensorOps::normalize< 3 >( ey );

    std::pair< real64, localIndex > thetaOrder[MAX_FACE_NODES];

    /// Sort nodes counterclockwise around face center
    for( localIndex n =0; n < numFaceNodes; ++n )
    {
      real64 v[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[faceNodes[n]] );
      LvArray::tensorOps::subtract< 3 >( v, fc );
      thetaOrder[n] = std::make_pair( atan2( LvArray::tensorOps::AiBi< 3 >( v, ey ), LvArray::tensorOps::AiBi< 3 >( v, ex ) ), faceNodes[n] );
    }

    std::sort( thetaOrder, thetaOrder + numFaceNodes );

    // Reorder nodes on face
    for( localIndex n =0; n < numFaceNodes; ++n )
    {
      faceNodes[n] = thetaOrder[n].second;
    }

    localIndex tempFaceNodes[MAX_FACE_NODES];

    localIndex firstIndexIndex = 0;
    for( localIndex n =0; n < numFaceNodes; ++n )
    {
      tempFaceNodes[n] = thetaOrder[n].second;
      if( tempFaceNodes[n] == firstNodeIndex )
      {
        firstIndexIndex = n;
      }
    }

    for( localIndex n=0; n < numFaceNodes; ++n )
    {
      const localIndex index = (firstIndexIndex + n) % numFaceNodes;
      faceNodes[n] = tempFaceNodes[index];
    }
  }
}

void FaceManager::extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                                   std::vector< std::vector< globalIndex > > & globalFaceNodes )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numFaces = size();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = this->nodeList().toViewConst();
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();

  globalFaceNodes.resize( numFaces );

  forAll< parallelHostPolicy >( numFaces, [&]( localIndex const & faceID )
  {
    std::vector< globalIndex > & curFaceGlobalNodes = globalFaceNodes[ faceID ];

    if( isDomainBoundary( faceID ) )
    {
      localIndex const numNodes = faceToNodeMap.sizeOfArray( faceID );
      curFaceGlobalNodes.resize( numNodes );

      for( localIndex a = 0; a < numNodes; ++a )
      {
        curFaceGlobalNodes[ a ]= nodeManager.localToGlobalMap()( faceToNodeMap( faceID, a ) );
      }

      std::sort( curFaceGlobalNodes.begin(), curFaceGlobalNodes.end() );
    }
    else
    {
      curFaceGlobalNodes.resize( 0 );
    }
  } );
}

void FaceManager::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::viewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementSubRegionListString() ));
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::elementListString() ));
}

localIndex FaceManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}

localIndex FaceManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex FaceManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString() ) );

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_nodeList.base().toViewConst(),
                                           m_unmappedGlobalIndicesInToNodes,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_nodeList.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_edgeList.base().toViewConst(),
                                           m_unmappedGlobalIndicesInToEdges,
                                           packList,
                                           this->localToGlobalMap(),
                                           m_edgeList.relatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::elementListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           this->m_toElements,
                                           packList,
                                           m_toElements.getElementRegionManager() );

  return packedSize;
}

localIndex FaceManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  // GEOSX_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_nodeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_nodeList.relatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOSX_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_edgeList,
                                     packList,
                                     m_unmappedGlobalIndicesInToEdges,
                                     this->globalToLocalMap(),
                                     m_edgeList.relatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOSX_ERROR_IF_NE( elementListString, viewKeyStruct::elementListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toElements,
                                     packList,
                                     m_toElements.getElementRegionManager(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void FaceManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_nodeList,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_edgeList,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );
}

void FaceManager::compressRelationMaps()
{
  m_nodeList.compress();
  m_edgeList.compress();
}

void FaceManager::enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices )
{
  arrayView1d< localIndex const > const childFaceIndices = getExtrinsicData< extrinsicMeshData::ChildIndex >();

  ObjectManagerBase::enforceStateFieldConsistencyPostTopologyChange ( targetIndices );

  for( localIndex const targetIndex : targetIndices )
  {
    localIndex const childIndex = childFaceIndices[targetIndex];
    if( childIndex != -1 )
    {
      LvArray::tensorOps::scaledCopy< 3 >( m_faceNormal[ targetIndex ], m_faceNormal[ childIndex ], -1 );
    }
  }
}

void FaceManager::depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                                    ElementRegionManager const & elemRegionManager )
{
  for( auto const & targetIndex : receivedFaces )
  {
    for( localIndex k=0; k<m_toElements.m_toElementRegion.size( 1 ); ++k )
    {
      localIndex const elemRegionIndex    = m_toElements.m_toElementRegion[targetIndex][k];
      localIndex const elemSubRegionIndex = m_toElements.m_toElementSubRegion[targetIndex][k];
      localIndex const elemIndex          = m_toElements.m_toElementIndex[targetIndex][k];

      if( elemRegionIndex!=-1 && elemSubRegionIndex!=-1 && elemIndex!=-1 )
      {
        CellElementSubRegion const & subRegion = elemRegionManager.getRegion( elemRegionIndex ).
                                                   getSubRegion< CellElementSubRegion >( elemSubRegionIndex );
        array2d< localIndex > const & downmap = subRegion.faceList();
        bool hasTargetIndex = false;

        for( localIndex a=0; a<downmap.size( 1 ); ++a )
        {
          localIndex const compositeLocalIndex = downmap[elemIndex][a];
          if( compositeLocalIndex==targetIndex )
          {
            hasTargetIndex=true;
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

REGISTER_CATALOG_ENTRY( ObjectManagerBase, FaceManager, string const &, Group * const )

}
