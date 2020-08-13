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
 * @file FaceElementSubRegion.cpp
 */

#include "FaceElementSubRegion.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

FaceElementSubRegion::FaceElementSubRegion( string const & name,
                                            dataRepository::Group * const parent ) :
  ElementSubRegionBase( name, parent ),
  m_unmappedGlobalIndicesInToNodes(),
  m_unmappedGlobalIndicesInToEdges(),
  m_unmappedGlobalIndicesInToFaces(),
  m_faceElementsToCells(),
  m_newFaceElements(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_toFacesRelation(),
  m_elementAperture(),
  m_elementArea(),
  m_elementRotationMatrix()
{
  SetElementType( "C3D8" );

  registerWrapper( viewKeyStruct::dNdXString, &m_dNdX )
    ->setSizedFromParent( 1 )
    ->reference()
    .resizeDimension< 3 >( 3 );

  registerWrapper( viewKeyStruct::detJString, &m_detJ )
    ->setSizedFromParent( 1 )
    ->reference();

  registerWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation )
    ->setDescription( "Map to the nodes attached to each FaceElement." );

  registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation )
    ->setDescription( "Map to the edges attached to each FaceElement." );

  registerWrapper( viewKeyStruct::faceListString, &m_toFacesRelation )
    ->setDescription( "Map to the faces attached to each FaceElement." )
    ->reference()
    .resize( 0, 2 );

  registerWrapper( viewKeyStruct::elementApertureString, &m_elementAperture )
    ->setApplyDefaultValue( -1.0 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_0 )
    ->setDescription( "The aperture of each FaceElement." );

  registerWrapper( viewKeyStruct::elementAreaString, &m_elementArea )
    ->setApplyDefaultValue( -1.0 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_2 )
    ->setDescription( "The area of each FaceElement." );

  registerWrapper( viewKeyStruct::elementRotationMatrixString,
                   &m_elementRotationMatrix )
    ->setApplyDefaultValue( 0.0 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_2 )
    ->setDescription( "The rotation matrix of each FaceElement." )
    ->reference()
    .resizeDimension< 1, 2 >( 3, 3 );

  registerWrapper( viewKeyStruct::faceElementsToCellRegionsString,
                   &m_faceElementsToCells.m_toElementRegion )
    ->setApplyDefaultValue( -1 )
    ->setPlotLevel( PlotLevel::NOPLOT )
    ->setDescription(
      "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::faceElementsToCellSubRegionsString,
                   &m_faceElementsToCells.m_toElementSubRegion )
    ->setApplyDefaultValue( -1 )
    ->setPlotLevel( PlotLevel::NOPLOT )
    ->setDescription(
      "A map of face element local indices to the cell local indices" );

  registerWrapper( viewKeyStruct::faceElementsToCellIndexString,
                   &m_faceElementsToCells.m_toElementIndex )
    ->setApplyDefaultValue( -1 )
    ->setPlotLevel( PlotLevel::NOPLOT )
    ->setDescription(
      "A map of face element local indices to the cell local indices" );

  registerWrapper< real64_array >( viewKeyStruct::creationMassString )
    ->setApplyDefaultValue( 0.0 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )
    ->setDescription(
      "The amount of remaining mass that was introduced when the FaceElement "
      "was created." );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  registerWrapper( viewKeyStruct::separationCoeffString, &m_separationCoefficient )
    ->setApplyDefaultValue( 0.0 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_1 )
    ->setDescription(
      "Scalar indicator of level of separation for a fracturing face." );
#endif

  m_faceElementsToCells.resize( 0, 2 );
  m_faceElementsToCells.setElementRegionManager(
    getParent()
      ->getParent()
      ->getParent()
      ->getParent()
      ->group_cast< ElementRegionManager * >() );

  m_numNodesPerElement = 8;
}

FaceElementSubRegion::~FaceElementSubRegion()
{}

void
FaceElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
  this->m_toEdgesRelation.SetRelatedObject( mesh->getEdgeManager() );
  this->m_toFacesRelation.SetRelatedObject( mesh->getFaceManager() );
}

void
FaceElementSubRegion::CalculateElementGeometricQuantities(
  localIndex const k,
  arrayView1d< real64 const > const & faceArea )
{
  m_elementArea[k] = faceArea[m_toFacesRelation[k][0]];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[m_toFacesRelation[k][0]];
}

void
FaceElementSubRegion::CalculateElementGeometricQuantities(
  localIndex const k,
  arrayView1d< real64 const > const & faceArea,
  arrayView3d< real64 const > const & faceRotationMatrix )
{
  localIndex const faceID = m_toFacesRelation( k, 0 );
  m_elementArea[k] = faceArea[faceID];
  m_elementVolume[k] = m_elementAperture[k] * faceArea[faceID];
  LvArray::tensorOps::copy< 3, 3 >( m_elementRotationMatrix[k],
                                    faceRotationMatrix[faceID] );
}

void
FaceElementSubRegion::CalculateElementGeometricQuantities(
  NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
  FaceManager const & faceManager )
{
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView3d< real64 const > const & faceRotationMatrix =
    faceManager.faceRotationMatrix();

  forAll< serialPolicy >( this->size(), [=]( localIndex const k ) {
    CalculateElementGeometricQuantities( k, faceArea, faceRotationMatrix );
  } );
}

localIndex
FaceElementSubRegion::PackUpDownMapsSize(
  arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

localIndex
FaceElementSubRegion::PackUpDownMaps(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex
FaceElementSubRegion::PackUpDownMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString ) );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               m_toNodesRelation.Base(),
                               m_unmappedGlobalIndicesInToNodes,
                               packList,
                               this->localToGlobalMap(),
                               m_toNodesRelation.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::edgeListString ) );
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               m_toEdgesRelation.Base(),
                               m_unmappedGlobalIndicesInToEdges,
                               packList,
                               this->localToGlobalMap(),
                               m_toEdgesRelation.RelatedObjectLocalToGlobal() );

  packedSize +=
    bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString ) );
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               m_toFacesRelation.Base().toViewConst(),
                               m_unmappedGlobalIndicesInToFaces,
                               packList,
                               this->localToGlobalMap(),
                               m_toFacesRelation.RelatedObjectLocalToGlobal() );

  packedSize += bufferOps::Pack< DOPACK >(
    buffer,
    string( viewKeyStruct::faceElementsToCellRegionsString ) );
  packedSize +=
    bufferOps::Pack< DOPACK >( buffer,
                               this->m_faceElementsToCells,
                               packList,
                               m_faceElementsToCells.getElementRegionManager() );

  return packedSize;
}

localIndex
FaceElementSubRegion::UnpackUpDownMaps(
  buffer_unit_type const *& buffer,
  localIndex_array & packList,
  bool const overwriteUpMaps,
  bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString );

  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_toNodesRelation,
                       packList,
                       m_unmappedGlobalIndicesInToNodes,
                       this->globalToLocalMap(),
                       m_toNodesRelation.RelatedObjectGlobalToLocal() );

  string edgeListString;
  unPackedSize += bufferOps::Unpack( buffer, edgeListString );
  GEOSX_ERROR_IF_NE( edgeListString, viewKeyStruct::edgeListString );

  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_toEdgesRelation,
                       packList,
                       m_unmappedGlobalIndicesInToEdges,
                       this->globalToLocalMap(),
                       m_toEdgesRelation.RelatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOSX_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString );

  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_toFacesRelation.Base(),
                       packList,
                       m_unmappedGlobalIndicesInToFaces,
                       this->globalToLocalMap(),
                       m_toFacesRelation.RelatedObjectGlobalToLocal() );

  string elementListString;
  unPackedSize += bufferOps::Unpack( buffer, elementListString );
  GEOSX_ERROR_IF_NE( elementListString,
                     viewKeyStruct::faceElementsToCellRegionsString );

  unPackedSize +=
    bufferOps::Unpack( buffer,
                       m_faceElementsToCells,
                       packList,
                       m_faceElementsToCells.getElementRegionManager(),
                       overwriteUpMaps );

  return unPackedSize;
}

void
FaceElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toEdgesRelation,
                                    m_unmappedGlobalIndicesInToEdges,
                                    clearIfUnmapped );

  ObjectManagerBase::FixUpDownMaps( m_toFacesRelation,
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void
FaceElementSubRegion::inheritGhostRankFromParentFace(
  FaceManager const * const faceManager,
  std::set< localIndex > const & indices )
{
  arrayView1d< integer const > const & faceGhostRank = faceManager->ghostRank();
  for( localIndex const & index : indices )
  {
    m_ghostRank[index] = faceGhostRank[m_toFacesRelation[index][0]];
  }
}

void
FaceElementSubRegion::ViewPackingExclusionList(
  SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::edgeListString ) );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::faceListString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::faceElementsToCellRegionsString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::faceElementsToCellSubRegionsString ) );
  exclusionList.insert(
    this->getWrapperIndex( viewKeyStruct::faceElementsToCellIndexString ) );
}

} /* namespace geosx */
