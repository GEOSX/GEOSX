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
 * @file Cylinder.cpp
 * @brief Generate a cylindrical geometry.
 * @param Center point of one (upper or lower) face of the cylinder
 * @param Center point of the other face of the cylinder
 * @param Radius of the cylinder
 */

#include "Cylinder.hpp"

namespace geosx
{
using namespace dataRepository;

Cylinder::Cylinder( const std::string & name, Group * const parent ) :
  SimpleGeometricObjectBase( name, parent ),
  m_point1 { 0.0, 0.0, 0.0 },
  m_point2 { 0.0, 0.0, 0.0 },
  m_radius { 0.0 }
{
  registerWrapper( viewKeyStruct::point1String, &m_point1 )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Center point of one (upper or lower) face of the cylinder" );

  registerWrapper( viewKeyStruct::point2String, &m_point2 )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Center point of the other face of the cylinder" );

  registerWrapper( viewKeyStruct::radiusString, &m_radius )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Radius of the cylinder" );
}

Cylinder::~Cylinder()
{}

bool
Cylinder::IsCoordInObject( const R1Tensor & coord ) const
{
  bool rval = false;

  R1Tensor axisVector, coord_minus_point1, projection, distance;
  realT height;

  axisVector = m_point2;
  axisVector -= m_point1;
  height = axisVector.Normalize();

  coord_minus_point1 = coord;
  coord_minus_point1 -= m_point1;

  projection = axisVector;
  projection *= Dot( axisVector, coord_minus_point1 );

  distance = coord_minus_point1;
  distance -= projection;

  if( distance.L2_Norm() < m_radius && projection.L2_Norm() < height )
  {
    rval = true;
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( SimpleGeometricObjectBase,
                        Cylinder,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
