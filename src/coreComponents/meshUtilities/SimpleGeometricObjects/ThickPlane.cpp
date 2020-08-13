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
 * @file ThickPlane.cpp
 */

#include "ThickPlane.hpp"

namespace geosx
{
using namespace dataRepository;

ThickPlane::ThickPlane( const std::string & name, Group * const parent ) :
  SimpleGeometricObjectBase( name, parent ),
  m_origin { 0.0, 0.0, 0.0 },
  m_normal { 0.0, 0.0, 1.0 },
  m_thickness { 0.0 }
{
  registerWrapper( viewKeyStruct::originString, &m_origin )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Origin point (x,y,z) of the plane (basically, any point on the plane)" );

  registerWrapper( viewKeyStruct::normalString, &m_normal )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)" );

  registerWrapper( viewKeyStruct::thicknessString, &m_thickness )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "The total thickness of the plane (with half to each side)" );
}

ThickPlane::~ThickPlane()
{}

void
ThickPlane::PostProcessInput()
{
  m_thickness *= 0.5;  // actually store the half-thickness
  GEOSX_ERROR_IF( m_thickness <= 0,
                  "Error: the plane appears to have zero or negative thickness" );

  m_normal.Normalize();
  GEOSX_ERROR_IF( std::fabs( m_normal.L2_Norm() - 1.0 ) > 1e-15,
                  "Error: could not properly normalize input normal." );
}

bool
ThickPlane::IsCoordInObject( const R1Tensor & coord ) const
{
  real64 normalDistance = 0.0;
  for( int i = 0; i < 3; ++i )
  {
    normalDistance += m_normal[i] * ( coord[i] - m_origin[i] );
  }

  return std::fabs( normalDistance ) <= m_thickness;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase,
                        ThickPlane,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
