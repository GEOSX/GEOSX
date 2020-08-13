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
 * @file SiloOutput.cpp
 */

#include "SiloOutput.hpp"

#include "common/TimingMacros.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/FunctionManager.hpp"

namespace geosx
{
using namespace dataRepository;

SiloOutput::SiloOutput( std::string const & name, Group * const parent ) :
  OutputBase( name, parent ),
  m_plotFileRoot( "plot" ),
  m_writeEdgeMesh( 0 ),
  m_writeFaceMesh( 0 ),
  m_writeCellElementMesh( 1 ),
  m_writeFaceElementMesh( 1 ),
  m_plotLevel()
{
  registerWrapper( viewKeysStruct::plotFileRoot, &m_plotFileRoot )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setApplyDefaultValue( "plot" )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::writeEdgeMesh, &m_writeEdgeMesh )
    ->setDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::writeFaceMesh, &m_writeFaceMesh )
    ->setDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::writeCellElementMesh, &m_writeCellElementMesh )
    ->setDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::writeFaceElementMesh, &m_writeFaceElementMesh )
    ->setDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );

  registerWrapper( viewKeysStruct::plotLevel, &m_plotLevel )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "" );
}

SiloOutput::~SiloOutput()
{}

void
SiloOutput::Execute( real64 const time_n,
                     real64 const dt,
                     integer const cycleNumber,
                     integer const eventCounter,
                     real64 const eventProgress,
                     Group * domain )
{
  GEOSX_MARK_FUNCTION;

  DomainPartition * domainPartition = Group::group_cast< DomainPartition * >( domain );
  SiloFile silo;

  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  MpiWrapper::Barrier( MPI_COMM_GEOSX );

  integer const numFiles = parallelThreads() == 0 ? size : parallelThreads();

  silo.setPlotLevel( m_plotLevel );
  silo.setWriteEdgeMesh( m_writeEdgeMesh );
  silo.setWriteFaceMesh( m_writeFaceMesh );
  silo.setWriteCellElementMesh( m_writeCellElementMesh );
  silo.setWriteFaceElementMesh( m_writeFaceElementMesh );
  silo.setPlotFileRoot( m_plotFileRoot );
  silo.Initialize( numFiles );
  silo.WaitForBatonWrite( rank, cycleNumber, eventCounter, false );
  silo.WriteDomainPartition( *domainPartition,
                             cycleNumber,
                             time_n + dt * eventProgress,
                             0 );
  silo.HandOffBaton();
  silo.ClearEmptiesFromMultiObjects( cycleNumber );
  silo.Finish();
}

REGISTER_CATALOG_ENTRY( OutputBase, SiloOutput, std::string const &, Group * const )
} /* namespace geosx */
