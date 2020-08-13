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
 * @file ChomboIO.cpp
 */
#include "ChomboIO.hpp"
#include "mesh/MeshLevel.hpp"
#include "managers/DomainPartition.hpp"
#include "fileIO/coupling/ChomboCoupler.hpp"
#include <string>
#include <fstream>
#include <chrono>

namespace geosx
{
using namespace dataRepository;

ChomboIO::ChomboIO( std::string const & name, Group * const parent ) :
  OutputBase( name, parent ),
  m_coupler( nullptr ),
  m_outputPath(),
  m_beginCycle( 0 ),
  m_inputPath( "/INVALID_INPUT_PATH" ),
  m_waitForInput(),
  m_useChomboPressures()
{
  registerWrapper( viewKeyStruct::outputPathString, &m_outputPath )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Path at which the geosx to chombo file will be written." );

  registerWrapper( viewKeyStruct::beginCycleString, &m_beginCycle )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Cycle at which the coupling will commence." );

  registerWrapper( viewKeyStruct::inputPathString, &m_inputPath )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDefaultValue( "/INVALID_INPUT_PATH" )
    ->setDescription( "Path at which the chombo to geosx file will be written." );

  registerWrapper( viewKeyStruct::waitForInputString, &m_waitForInput )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDefaultValue( 0 )
    ->setDescription(
      "True iff geosx should wait for chombo to write out a file. When true "
      "the inputPath must be set." );

  registerWrapper( viewKeyStruct::useChomboPressuresString, &m_useChomboPressures )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDefaultValue( 0 )
    ->setDescription(
      "True iff geosx should use the pressures chombo writes out." );
}

ChomboIO::~ChomboIO()
{
  delete m_coupler;
  m_coupler = nullptr;
}

void
ChomboIO::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                   real64 const dt,
                   integer const cycleNumber,
                   integer const GEOSX_UNUSED_PARAM( eventCounter ),
                   real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                   dataRepository::Group * const domain )
{
  if( m_coupler == nullptr )
  {
    GEOSX_ERROR_IF( m_waitForInput && m_inputPath == "/INVALID_INPUT_PATH",
                    "Waiting for input but no input path was specified." );

    DomainPartition * const domainPartition =
      Group::group_cast< DomainPartition * >( domain );
    MeshLevel * const meshLevel = domainPartition->getMeshBody( 0 )->getMeshLevel( 0 );
    m_coupler =
      new ChomboCoupler( MPI_COMM_GEOSX, m_outputPath, m_inputPath, *meshLevel );
  }

  if( cycleNumber < m_beginCycle )
  {
    return;
  }

  m_coupler->write( dt );

  if( m_waitForInput )
  {
    m_coupler->read( m_useChomboPressures );
  }
}

REGISTER_CATALOG_ENTRY( OutputBase, ChomboIO, std::string const &, Group * const )
} /* namespace geosx */
