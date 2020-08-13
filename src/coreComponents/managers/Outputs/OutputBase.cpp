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
 * @file OutputBase.cpp
 */

#include "OutputBase.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{
using namespace dataRepository;

OutputBase::OutputBase( std::string const & name, Group * const parent ) :
  ExecutableGroup( name, parent ),
  m_childDirectory(),
  m_parallelThreads( 1 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeysStruct::childDirectoryString, &m_childDirectory )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Child directory path" );

  registerWrapper( viewKeysStruct::parallelThreadsString, &m_parallelThreads )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Number of plot files." );
}

OutputBase::~OutputBase()
{}

OutputBase::CatalogInterface::CatalogType &
OutputBase::GetCatalog()
{
  static OutputBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void
OutputBase::InitializePreSubGroups( Group * const GEOSX_UNUSED_PARAM( group ) )
{
  // This command doesn't seem to work anymore
  // SetupDirectoryStructure();
}

void
OutputBase::SetupDirectoryStructure()
{
  string childDirectory = m_childDirectory;

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( rank == 0 )
  {
    if( !childDirectory.empty() )
    {
      makeDirsForPath( childDirectory );
    }
  }
}

} /* namespace geosx */
