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
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
/**
 * @file ConstitutiveBase.cpp
 */

#include "ConstitutiveBase.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{
ConstitutiveBase::ConstitutiveBase( std::string const & name, Group * const parent ) :
  Group( name, parent ),
  m_numQuadraturePoints( 1 ),
  m_constitutiveDataGroup( nullptr )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

ConstitutiveBase::~ConstitutiveBase()
{}

ConstitutiveBase::CatalogInterface::CatalogType &
ConstitutiveBase::GetCatalog()
{
  static ConstitutiveBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void
ConstitutiveBase::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  m_numQuadraturePoints = numConstitutivePointsPerParentIndex;
  m_constitutiveDataGroup = parent;

  for( auto & group : this->GetSubGroups() )
  {
    for( auto & wrapper : group.second->wrappers() )
    {
      if( wrapper.second->sizedFromParent() )
      {
        string const wrapperName = wrapper.first;
        parent
          ->registerWrapper( makeFieldName( this->getName(), wrapperName ),
                             wrapper.second->clone( wrapperName, parent ) )
          ->setRestartFlags( RestartFlags::NO_WRITE );
      }
    }
  }

  for( auto & wrapper : this->wrappers() )
  {
    if( wrapper.second->sizedFromParent() )
    {
      string const wrapperName = wrapper.first;
      parent
        ->registerWrapper( makeFieldName( this->getName(), wrapperName ),
                           wrapper.second->clone( wrapperName, parent ) )
        ->setRestartFlags( RestartFlags::NO_WRITE );
    }
  }
}

void
ConstitutiveBase::resize( localIndex newsize )
{
  Group::resize( newsize );
}

void
ConstitutiveBase::DeliverClone( string const & GEOSX_UNUSED_PARAM( name ),
                                Group * const GEOSX_UNUSED_PARAM( parent ),
                                std::unique_ptr< ConstitutiveBase > & clone ) const
{
  GEOSX_ASSERT( clone );
  clone->forWrappers( [&]( WrapperBase & wrapper ) {
    wrapper.CopyWrapperAttributes( *( this->getWrapperBase( wrapper.getName() ) ) );
  } );
}

}  // namespace constitutive
} /* namespace geosx */
