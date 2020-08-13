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
 * @file EmbeddedSurfaceRegion.cpp
 */

#include "EmbeddedSurfaceRegion.hpp"

#include "EdgeManager.hpp"

namespace geosx
{
using namespace dataRepository;

EmbeddedSurfaceRegion::EmbeddedSurfaceRegion( string const & name,
                                              Group * const parent ) :
  ElementRegionBase( name, parent )
{
  this->GetGroup( viewKeyStruct::elementSubRegions )
    ->RegisterGroup< EmbeddedSurfaceSubRegion >( "default" );

  registerWrapper( viewKeyStruct::defaultApertureString, &m_defaultAperture )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "The default aperture of for new embedded surface Elements." );
}

EmbeddedSurfaceRegion::~EmbeddedSurfaceRegion()
{}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        EmbeddedSurfaceRegion,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
