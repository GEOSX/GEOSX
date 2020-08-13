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

#include "FieldSpecificationManager.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include "constitutive/ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
FieldSpecificationManager::FieldSpecificationManager( string const & name,
                                                      Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

FieldSpecificationManager &
FieldSpecificationManager::get()
{
  static FieldSpecificationManager bcman( "FieldSpecifications", nullptr );
  return bcman;
}

FieldSpecificationManager::~FieldSpecificationManager()
{
  // TODO Auto-generated destructor stub
}

Group *
FieldSpecificationManager::CreateChild( string const & childKey,
                                        string const & childName )
{
  std::unique_ptr< FieldSpecificationBase > bc =
    FieldSpecificationBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup( childName, std::move( bc ) );
}

void
FieldSpecificationManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from BoundaryConditionBase here
  for( auto & catalogIter : FieldSpecificationBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

void
FieldSpecificationManager::ApplyInitialConditions( Group * domain ) const
{
  Apply( 0.0,
         domain,
         "",
         "",
         [&]( FieldSpecificationBase const * const bc,
              string const &,
              SortedArrayView< localIndex const > const & targetSet,
              Group * const targetGroup,
              string const fieldName ) {
           bc->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                           0.0,
                                                           targetGroup,
                                                           fieldName );
         } );
}

} /* namespace geosx */
