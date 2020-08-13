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
 * @file FunctionManager.cpp
 */

#include "FunctionManager.hpp"

#include "dataRepository/Group.hpp"
#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
using namespace dataRepository;

FunctionManager::FunctionManager( const std::string & name, Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FunctionManager::~FunctionManager()
{
  // TODO Auto-generated destructor stub
}

Group *
FunctionManager::CreateChild( string const & functionCatalogKey,
                              string const & functionName )
{
  GEOSX_LOG_RANK_0( "   " << functionCatalogKey << ": " << functionName );
  std::unique_ptr< FunctionBase > function =
    FunctionBase::CatalogInterface::Factory( functionCatalogKey, functionName, this );
  return this->RegisterGroup< FunctionBase >( functionName, std::move( function ) );
}

void
FunctionManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from FunctionBase here
  for( auto & catalogIter : FunctionBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

}  // namespace geosx
