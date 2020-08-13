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
 * @file RestartOutput.cpp
 */

#include "RestartOutput.hpp"
#include "fileIO/silo/SiloFile.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/Functions/FunctionManager.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"

namespace geosx
{
using namespace dataRepository;

RestartOutput::RestartOutput( std::string const & name, Group * const parent ) :
  OutputBase( name, parent )
{}

RestartOutput::~RestartOutput()
{}

void
RestartOutput::Execute( real64 const GEOSX_UNUSED_PARAM( time_n ),
                        real64 const GEOSX_UNUSED_PARAM( dt ),
                        integer const cycleNumber,
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        Group * domain )
{
  GEOSX_MARK_FUNCTION;

  DomainPartition * domainPartition = Group::group_cast< DomainPartition * >( domain );
  ProblemManager * problemManager =
    Group::group_cast< ProblemManager * >( domainPartition->getParent() );

  // Ignoring the eventProgress indicator for now to be compliant with the integrated test repo
  // integer const eventProgressPercent = static_cast<integer const>(eventProgress * 100.0);
  char fileName[200] = { 0 };
  sprintf( fileName,
           "%s_%s_%09d",
           problemManager->getProblemName().c_str(),
           "restart",
           cycleNumber );

  problemManager->prepareToWrite();
  FunctionManager::Instance().prepareToWrite();
  FieldSpecificationManager::get().prepareToWrite();
  writeTree( fileName );
  problemManager->finishWriting();
  FunctionManager::Instance().finishWriting();
  FieldSpecificationManager::get().finishWriting();
}

REGISTER_CATALOG_ENTRY( OutputBase, RestartOutput, std::string const &, Group * const )
} /* namespace geosx */
