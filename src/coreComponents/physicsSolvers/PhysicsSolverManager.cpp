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
 * @file PhysicsSolverManager.cpp
 */

#include "PhysicsSolverManager.hpp"

#include "SolverBase.hpp"

namespace geosx
{
using namespace dataRepository;

PhysicsSolverManager::PhysicsSolverManager( std::string const & name,
                                            Group * const parent ) :
  Group( name, parent ),
  m_gravityVector( R1Tensor( 0.0 ) )
{
  setInputFlags( InputFlags::REQUIRED );

  this->registerWrapper( viewKeyStruct::gravityVectorString, &m_gravityVector )
    ->setApplyDefaultValue( { 0.0, 0.0, -9.81 } )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Gravity vector used in the physics solvers" );
}

PhysicsSolverManager::~PhysicsSolverManager()
{}

//START_SPHINX_INCLUDE_00
Group *
PhysicsSolverManager::CreateChild( string const & childKey,
                                   string const & childName )
{
  Group * rval = nullptr;
  if( SolverBase::CatalogInterface::hasKeyName( childKey ) )
  {
    GEOSX_LOG_RANK_0( "Adding Solver of type " << childKey << ", named "
                                               << childName );
    rval = RegisterGroup(
      childName,
      SolverBase::CatalogInterface::Factory( childKey, childName, this ) );
  }
  return rval;
}

void
PhysicsSolverManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from SolverBase here
  for( auto & catalogIter : SolverBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

} /* namespace geosx */
