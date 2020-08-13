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
 * @file BoundedPlane.cpp
 */

#include "BoundedPlane.hpp"

namespace geosx
{
using namespace dataRepository;

BoundedPlane::BoundedPlane( const std::string & name, Group * const parent ) :
  SimpleGeometricObjectBase( name, parent ),
  m_origin { 0.0, 0.0, 0.0 },
  m_normal { 0.0, 0.0, 1.0 },
  m_lengthVector { 0.0, 0.0, 0.0 },
  m_widthVector { 0.0, 0.0, 0.0 }
{
  registerWrapper( viewKeyStruct::originString, &m_origin )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Origin point (x,y,z) of the plane (basically, any point on the plane)" );

  registerWrapper( viewKeyStruct::normalString, &m_normal )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)" );

  registerWrapper( viewKeyStruct::mLengthVectorString, &m_lengthVector )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Tangent vector defining the orthonormal basis along with the normal." );

  registerWrapper( viewKeyStruct::mWidthVectorString, &m_widthVector )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Tangent vector defining the orthonormal basis along with the normal." );

  registerWrapper( viewKeyStruct::dimensionsString, &m_dimensions )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Length and width of the bounded plane" );

  m_points.resize( 4 );
}

BoundedPlane::~BoundedPlane()
{}

void
BoundedPlane::PostProcessInput()
{
  // Make sure that you have an orthonormal basis.
  m_normal.Normalize();
  m_lengthVector.Normalize();
  m_widthVector.Normalize();

  //Check if they are all orthogonal
  R1Tensor vector = Cross( m_lengthVector, m_widthVector );
  GEOSX_ERROR_IF( std::fabs( std::fabs( Dot( m_normal, vector ) ) - 1 ) > 1e-15 ||
                    std::fabs( Dot( m_widthVector, m_lengthVector ) ) > 1e-15,
                  "Error: the 3 vectors provided in the BoundedPlane do not "
                  "form an orthonormal basis!" );
  GEOSX_ERROR_IF( m_dimensions.size() != 2,
                  "Error: Need to provide both length and width!" );

  findRectangleLimits();
}

void
BoundedPlane::findRectangleLimits()
{
  R1Tensor lengthVec = m_lengthVector;
  R1Tensor widthVec = m_widthVector;
  lengthVec *= 0.5 * m_dimensions[0];
  widthVec *= 0.5 * m_dimensions[1];

  for( int i = 0; i < 4; i++ ) m_points[i] = m_origin;

  m_points[0] -= lengthVec;
  m_points[0] -= widthVec;

  m_points[1] += lengthVec;
  m_points[1] -= widthVec;

  m_points[2] += lengthVec;
  m_points[2] += widthVec;

  m_points[3] -= lengthVec;
  m_points[3] += widthVec;

  if( getLogLevel() > 1 )
  {
    GEOSX_LOG_RANK_0( "Point A: " << m_points[0] );
    GEOSX_LOG_RANK_0( "Point B: " << m_points[1] );
    GEOSX_LOG_RANK_0( "Point C: " << m_points[2] );
    GEOSX_LOG_RANK_0( "Point D: " << m_points[3] );
  }
}

bool
BoundedPlane::IsCoordInObject( const R1Tensor & coord ) const
{
  bool isInside = true;

  R1Tensor dummy = coord;
  dummy -= m_origin;

  // 1. Check if point is on the plane
  if( std::abs( Dot( dummy, m_normal ) ) < 1e-15 )
  {
    R1Tensor vec = coord;
    R1Tensor abVec = m_points[1];
    R1Tensor adVec = m_points[3];

    vec -= m_points[0];
    abVec -= m_points[0];
    adVec -= m_points[0];

    // 2. Check if it is inside the rectangle
    if( Dot( vec, abVec ) < 0 || Dot( vec, abVec ) > Dot( abVec, abVec ) )
      isInside = false;

    if( Dot( vec, adVec ) < 0 || Dot( vec, adVec ) > Dot( adVec, adVec ) )
      isInside = false;
  }
  else
  {
    isInside = false;
  }

  return isInside;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase,
                        BoundedPlane,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
