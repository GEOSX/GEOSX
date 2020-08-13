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

/*
 * @file FiniteVolumeManager.cpp
 *
 */

#include "FiniteVolumeManager.hpp"

#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;

FiniteVolumeManager::FiniteVolumeManager( string const & name, Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FiniteVolumeManager::~FiniteVolumeManager()
{}

Group *
FiniteVolumeManager::CreateChild( string const & childKey,
                                  string const & childName )
{
  std::unique_ptr< FluxApproximationBase > approx =
    FluxApproximationBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup< FluxApproximationBase >( childName, std::move( approx ) );
}

void
FiniteVolumeManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FluxApproximationBase here
  for( auto & catalogIter : FluxApproximationBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

FluxApproximationBase const &
FiniteVolumeManager::getFluxApproximation(
  std::string const & name ) const
{
  return getGroupReference< FluxApproximationBase >( name );
}

FluxApproximationBase &
FiniteVolumeManager::getFluxApproximation(
  std::string const & name )
{
  return getGroupReference< FluxApproximationBase >( name );
}

}  // namespace geosx
