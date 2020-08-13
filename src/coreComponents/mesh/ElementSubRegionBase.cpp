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
 * @file CellBase.cpp
 */

#include "ElementSubRegionBase.hpp"

namespace geosx
{
using namespace dataRepository;

ElementSubRegionBase::ElementSubRegionBase( string const & name, Group * const parent ) :
  ObjectManagerBase( name, parent ),
  m_constitutiveModels( groupKeyStruct::constitutiveModelsString, this ),
  m_numNodesPerElement(),
  m_numEdgesPerElement(),
  m_numFacesPerElement(),
  m_elementCenter(),
  m_elementVolume()
{
  RegisterGroup( groupKeyStruct::constitutiveModelsString, &m_constitutiveModels )
    ->setSizedFromParent( 1 );

  registerWrapper( viewKeyStruct::numNodesPerElementString, &m_numNodesPerElement );

  registerWrapper( viewKeyStruct::numEdgesPerElementString, &m_numEdgesPerElement );

  registerWrapper( viewKeyStruct::numFacesPerElementString, &m_numFacesPerElement );

  registerWrapper( viewKeyStruct::elementCenterString, &m_elementCenter )
    ->setPlotLevel( PlotLevel::LEVEL_1 )
    ->reference()
    .resizeDimension< 1 >( 3 );

  registerWrapper( viewKeyStruct::elementVolumeString, &m_elementVolume )
    ->setPlotLevel( PlotLevel::LEVEL_1 );
}

ElementSubRegionBase::~ElementSubRegionBase()
{}

void
ElementSubRegionBase::SetElementType( string const & elementType )
{
  m_elementTypeString = elementType;
}

std::vector< int >
ElementSubRegionBase::getVTKNodeOrdering() const
{
  if( !m_elementTypeString.compare( 0, 4, "C3D4" ) ) return { 1, 0, 2, 3 };
  if( !m_elementTypeString.compare( 0, 4, "C3D8" ) )
    return { 0, 1, 3, 2, 4, 5, 7, 6 };
  if( !m_elementTypeString.compare( 0, 4, "C3D6" ) )
    return { 0, 3, 4, 1, 2, 5, 0, 0 };
  if( !m_elementTypeString.compare( 0, 4, "C3D5" ) )
    return { 0, 3, 2, 1, 4, 0, 0, 0 };
  if( !m_elementTypeString.compare( 0, 4, "BEAM" ) ) return { 0, 1 };

  GEOSX_ERROR( "Unrecognized elementType: " << m_elementTypeString );
  return {};
}

} /* namespace geosx */
