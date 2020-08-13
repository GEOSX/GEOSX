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
 * @file ConstitutiveManager.cpp
 */

#include "ConstitutiveManager.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
ConstitutiveManager::ConstitutiveManager( string const & name, Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );
}

ConstitutiveManager::~ConstitutiveManager()
{}

Group *
ConstitutiveManager::CreateChild( string const & childKey,
                                  string const & childName )
{
  std::unique_ptr< ConstitutiveBase > material =
    ConstitutiveBase::CatalogInterface::Factory( childKey, childName, this );
  return RegisterGroup< ConstitutiveBase >( childName, std::move( material ) );
}

void
ConstitutiveManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from ConstitutiveBase here
  for( auto & catalogIter : ConstitutiveBase::GetCatalog() )
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}

ConstitutiveBase *
ConstitutiveManager::HangConstitutiveRelation(
  string const & constitutiveRelationInstanceName,
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex ) const
{
  ConstitutiveBase const * const constitutiveRelation =
    GetConstitutiveRelation( constitutiveRelationInstanceName );

  std::unique_ptr< ConstitutiveBase > material;
  constitutiveRelation->DeliverClone( constitutiveRelationInstanceName,
                                      parent,
                                      material );

  material->AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  dataRepository::Group * constitutiveGroup =
    parent->GetGroup( groupKeyStruct::constitutiveModelsString );
  if( constitutiveGroup == nullptr )
  {
    constitutiveGroup =
      parent->RegisterGroup( groupKeyStruct::constitutiveModelsString )
        ->setSizedFromParent( 1 );
    constitutiveGroup->resize( parent->size() );
  }

  ConstitutiveBase * const rval =
    constitutiveGroup->RegisterGroup< ConstitutiveBase >(
      constitutiveRelationInstanceName,
      std::move( material ) );
  rval->setSizedFromParent( 1 );
  rval->resize( constitutiveGroup->size() );
  return rval;
}

}  // namespace constitutive

} /* namespace geosx */
