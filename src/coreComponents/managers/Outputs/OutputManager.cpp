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
 * @file OutputManager.cpp
 */

#include "OutputManager.hpp"
#include "SiloOutput.hpp"

namespace geosx
{
using namespace dataRepository;

OutputManager::OutputManager( std::string const & name, Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );
}

OutputManager::~OutputManager()
{}

Group *
OutputManager::CreateChild( string const & childKey, string const & childName )
{
  GEOSX_LOG_RANK_0( "Adding Output: " << childKey << ", " << childName );
  std::unique_ptr< OutputBase > output =
    OutputBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup< OutputBase >( childName, std::move( output ) );
}

void
OutputManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from OutputBase here
  for( auto & catalogIter : OutputBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace geosx */
