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
 * @file EmbeddedSurfaceSubRegion.cpp
 */

#include "EmbeddedSurfaceSubRegion.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "NodeManager.hpp"
#include "MeshLevel.hpp"
#include "BufferOps.hpp"

namespace geosx
{
using namespace dataRepository;

EmbeddedSurfaceSubRegion::EmbeddedSurfaceSubRegion( string const & name,
                                                    dataRepository::Group * const parent ) :
  ElementSubRegionBase( name, parent ),
  m_normalVector(),
  m_tangentVector1(),
  m_tangentVector2(),
  m_embeddedSurfaceToRegion(),
  m_embeddedSurfaceToSubRegion(),
  m_embeddedSurfaceToCell(),
  m_toNodesRelation(),
  m_toEdgesRelation(),
  m_elementAperture(),
  m_elementArea(),
  m_numOfJumpEnrichments( 3 )
{
  registerWrapper( viewKeyStruct::regionListString, &m_embeddedSurfaceToRegion )
    ->setDescription( "Map to the region cut by each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::subregionListString,
                   &m_embeddedSurfaceToSubRegion )
    ->setDescription( "Map to the subregion cut by each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::nodeListString, &m_toNodesRelation )
    ->setDescription( "Map to the nodes attached to each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::edgeListString, &m_toEdgesRelation )
    ->setDescription( "Map to the edges." );

  registerWrapper( viewKeyStruct::cellListString, &m_embeddedSurfaceToCell )
    ->setDescription( "Map to the cells." );

  registerWrapper( viewKeyStruct::normalVectorString, &m_normalVector )
    ->setDescription( "Unit normal vector to the embedded surface." );

  registerWrapper( viewKeyStruct::elementApertureString, &m_elementAperture )
    ->setApplyDefaultValue( 1.0e-5 )
    ->setPlotLevel( dataRepository::PlotLevel::LEVEL_0 )
    ->setDescription( "The aperture of each EmbeddedSurface." );

  registerWrapper( viewKeyStruct::elementAreaString, &m_elementArea )
    ->setApplyDefaultValue( -1.0 )
    ->setDescription( "The area of each EmbeddedSurface element." );
}

EmbeddedSurfaceSubRegion::~EmbeddedSurfaceSubRegion()
{}

void
EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities( localIndex const k )
{
  m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
}

void
EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities(
  NodeManager const & GEOSX_UNUSED_PARAM( nodeManager ),
  FaceManager const & GEOSX_UNUSED_PARAM( facemanager ) )
{
  // loop over the elements
  forAll< serialPolicy >( this->size(), [=]( localIndex const k ) {
    m_elementVolume[k] = m_elementAperture[k] * m_elementArea[k];
  } );
}

void
EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface( localIndex const cellIndex,
                                                 R1Tensor normalVector )
{
  m_embeddedSurfaceToCell.emplace_back( cellIndex );
  m_normalVector.emplace_back( normalVector );

  // resize
  this->resize( this->size() + 1 );
}

bool
EmbeddedSurfaceSubRegion::AddNewEmbeddedSurface(
  localIndex const cellIndex,
  localIndex const regionIndex,
  localIndex const subRegionIndex,
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager,
  FixedOneToManyRelation const & cellToEdges,
  BoundedPlane const * fracture )
{
  /* The goal is to add an embeddedSurfaceElem if it is contained within the BoundedPlane
   *
   * A. Identify whether the cell falls within the bounded plane or not
   *
   * we loop over each edge:
   *
   *   1. check if it is cut by the plane using the Dot product between the distance of each node
   *   from the origin and the normal vector.
   *   2. If an edge is cut by the plane we look for the intersection between a line and a plane.
   *   3. Once we have the intersection we check if it falls inside the bounded plane.
   *
   * Only elements for which all intersection points fall within the fracture plane limits will be added.
   * If the fracture does not cut through the entire element we will just chop it (it's a discretization error).
   *
   * B. Once we know the element has to be added we compute it's geometric properties:
   * - Surface Area: this is trivial given the intersection points as long as they are ordered either CW or CCW.
   * - centre:
   * - Volume:
   * - Heaviside:
   */

  bool addEmbeddedElem = true;
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord =
    nodeManager.referencePosition();
  EdgeManager::NodeMapType::ViewTypeConst const & edgeToNodes =
    edgeManager.nodeList();
  R1Tensor origin = fracture->getCenter();
  R1Tensor normalVector = fracture->getNormal();
  localIndex edgeIndex;
  R1Tensor lineDir, dist, point;
  real64 prodScalarProd;

  array1d< R1Tensor > intersectionPoints;
  for( localIndex ke = 0; ke < cellToEdges.size( 1 ); ke++ )
  {
    edgeIndex = cellToEdges[cellIndex][ke];
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][0]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd = LvArray::tensorOps::AiBi< 3 >( dist, normalVector );
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][1]] );
    LvArray::tensorOps::subtract< 3 >( dist, origin );
    prodScalarProd *= LvArray::tensorOps::AiBi< 3 >( dist, normalVector );

    // check if the plane intersects the edge
    if( prodScalarProd < 0 )
    {
      lineDir =
        LVARRAY_TENSOROPS_INIT_LOCAL_3( nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( lineDir,
                                         nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::normalize< 3 >( lineDir );
      //find the intersection point
      point = computationalGeometry::LinePlaneIntersection(
        lineDir,
        nodesCoord[edgeToNodes[edgeIndex][0]],
        normalVector,
        origin );

      // Check if the point is inside the fracture (bounded plane)
      if( !fracture->IsCoordInObject( point ) )
      {
        addEmbeddedElem = false;
      }
      intersectionPoints.emplace_back( point );
    }
  }  //end of edge loop

  if( addEmbeddedElem )
  {
    m_embeddedSurfaceToCell.emplace_back( cellIndex );
    m_embeddedSurfaceToRegion.emplace_back( regionIndex );
    m_embeddedSurfaceToSubRegion.emplace_back( subRegionIndex );
    m_normalVector.emplace_back( normalVector );
    m_tangentVector1.emplace_back( fracture->getWidthVector() );
    m_tangentVector2.emplace_back( fracture->getLengthVector() );
    // resize
    this->resize( this->size() + 1 );
    this->CalculateElementGeometricQuantities( intersectionPoints,
                                               this->size() - 1 );
  }
  return addEmbeddedElem;
}

void
EmbeddedSurfaceSubRegion::inheritGhostRank(
  array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank )
{
  arrayView1d< integer > const & ghostRank = this->ghostRank();
  for( localIndex k = 0; k < size(); ++k )
  {
    ghostRank[k] =
      cellGhostRank[m_embeddedSurfaceToRegion[k]]
                   [m_embeddedSurfaceToSubRegion[k]][m_embeddedSurfaceToCell[k]];
  }
}

void
EmbeddedSurfaceSubRegion::CalculateElementGeometricQuantities(
  array1d< R1Tensor > const intersectionPoints,
  localIndex const k )
{
  for( localIndex p = 0; p < intersectionPoints.size(); p++ )
  {
    // TODO change to LvArray::tensorOps::add
    m_elementCenter( k, 0 ) += intersectionPoints[p][0];
    m_elementCenter( k, 1 ) += intersectionPoints[p][1];
    m_elementCenter( k, 2 ) += intersectionPoints[p][2];
  }

  m_elementArea[k] =
    computationalGeometry::ComputeSurfaceArea( intersectionPoints,
                                               intersectionPoints.size(),
                                               m_normalVector[k] );

  LvArray::tensorOps::scale< 3 >( m_elementCenter[k],
                                  1.0 / intersectionPoints.size() );
  this->CalculateElementGeometricQuantities( k );
}

void
EmbeddedSurfaceSubRegion::setupRelatedObjectsInRelations(
  MeshLevel const * const mesh )
{
  this->m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

real64
EmbeddedSurfaceSubRegion::ComputeHeavisideFunction(
  ArraySlice< real64 const, 1, nodes::REFERENCE_POSITION_USD - 1 > const nodeCoord,
  localIndex const k ) const
{
  real64 heaviside;
  R1Tensor distanceVector;
  LvArray::tensorOps::copy< 3 >( distanceVector, nodeCoord );
  LvArray::tensorOps::subtract< 3 >( distanceVector, m_elementCenter[k] );

  heaviside =
    LvArray::tensorOps::AiBi< 3 >( distanceVector, m_normalVector[k] ) > 0 ? 1 : 0;

  return heaviside;
}

void
EmbeddedSurfaceSubRegion::getIntersectionPoints(
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager,
  ElementRegionManager const & elemManager,
  array1d< R1Tensor > & intersectionPoints,
  array1d< localIndex > & connectivityList,
  array1d< int > & offSet ) const
{
  offSet.resize( size() + 1 );
  offSet.setValues< serialPolicy >( 0 );
  for( localIndex k = 0; k < size(); k++ )
  {
    ComputeIntersectionPoints( nodeManager,
                               edgeManager,
                               elemManager,
                               intersectionPoints,
                               connectivityList,
                               offSet,
                               k );
  }
}

void
EmbeddedSurfaceSubRegion::ComputeIntersectionPoints(
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager,
  ElementRegionManager const & elemManager,
  array1d< R1Tensor > & intersectionPoints,
  array1d< localIndex > & connectivityList,
  array1d< int > & offSet,
  localIndex const k ) const
{
  // I ll use this for plotting
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord =
    nodeManager.referencePosition();
  EdgeManager::NodeMapType::ViewTypeConst const & edgeToNodes =
    edgeManager.nodeList();

  FixedOneToManyRelation const & cellToEdges =
    elemManager.GetRegion( m_embeddedSurfaceToRegion[k] )
      ->GetSubRegion< CellElementSubRegion >( m_embeddedSurfaceToSubRegion[k] )
      ->edgeList();

  // Initialize variables
  localIndex edgeIndex;
  R1Tensor lineDir, dist, point;
  real64 prodScalarProd;
  bool isNew;
  R1Tensor distance, planeCentre;
  array1d< R1Tensor > localPoints;

  int count = 0;
  if( k > 0 )
  {
    count = offSet[k];
  }

  for( localIndex ke = 0; ke < cellToEdges.size( 1 ); ke++ )
  {
    edgeIndex = cellToEdges[m_embeddedSurfaceToCell[k]][ke];
    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][0]] );
    LvArray::tensorOps::subtract< 3 >( dist, m_elementCenter[k] );
    prodScalarProd = LvArray::tensorOps::AiBi< 3 >( dist, m_normalVector[k] );

    LvArray::tensorOps::copy< 3 >( dist, nodesCoord[edgeToNodes[edgeIndex][1]] );
    LvArray::tensorOps::subtract< 3 >( dist, m_elementCenter[k] );
    prodScalarProd *= LvArray::tensorOps::AiBi< 3 >( dist, m_normalVector[k] );

    if( prodScalarProd < 0 )
    {
      count += 1;
      LvArray::tensorOps::copy< 3 >( lineDir, nodesCoord[edgeToNodes[edgeIndex][0]] );
      LvArray::tensorOps::subtract< 3 >( lineDir,
                                         nodesCoord[edgeToNodes[edgeIndex][1]] );
      LvArray::tensorOps::normalize< 3 >( lineDir );
      LvArray::tensorOps::copy< 3 >( planeCentre, m_elementCenter[k] );
      point = computationalGeometry::LinePlaneIntersection(
        lineDir,
        nodesCoord[edgeToNodes[edgeIndex][0]],
        m_normalVector[k],
        planeCentre );

      localPoints.emplace_back( point );

      isNew = true;
      for( int i = 0; i < intersectionPoints.size(); i++ )
      {
        LvArray::tensorOps::copy< 3 >( distance, point );
        LvArray::tensorOps::subtract< 3 >( distance, intersectionPoints[i] );
        if( LvArray::tensorOps::l2Norm< 3 >( distance ) < 1e-9 )
        {
          isNew = false;
          //pointIndex = i;
          break;
        }
      }

      if( isNew )
      {
        intersectionPoints.emplace_back( point );
        //pointIndex = intersectionPoints.size() - 1;
      }
    }
  }  //end of edge loop

  // Reorder the points CCW and then add the correct index to the connectivity list
  localPoints = computationalGeometry::orderPointsCCW( localPoints,
                                                       localPoints.size(),
                                                       m_normalVector[k] );
  for( localIndex j = 0; j < localPoints.size(); j++ )
  {
    for( localIndex h = 0; h < intersectionPoints.size(); h++ )
    {
      LvArray::tensorOps::copy< 3 >( distance, localPoints[j] );
      LvArray::tensorOps::subtract< 3 >( distance, intersectionPoints[h] );
      if( LvArray::tensorOps::l2Norm< 3 >( distance ) < 1e-9 )
      {
        connectivityList.emplace_back( h );
      }
    }
  }
  offSet[k + 1] = count;
}

} /* namespace geosx */
