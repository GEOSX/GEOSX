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

// Source includes
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// TPL includes
#ifdef GEOSX_USE_CALIPER
  #include <caliper/cali-manager.h>
#endif

// System includes
#include <cmath>
#include <iostream>
#include <sys/time.h>

using namespace geosx;

int
main( int argc, char * argv[] )
{
  basicSetup( argc, argv, true );

  GEOSX_MARK_FUNCTION_BEGIN;

  printTypeSummary();

  timeval tim;
  gettimeofday( &tim, nullptr );
  real64 t_start = tim.tv_sec + ( tim.tv_usec / 1000000.0 );

  std::string restartFileName;
  bool restart = ProblemManager::ParseRestart( restartFileName );
  if( restart )
  {
    GEOSX_LOG_RANK_0( "Loading restart file " << restartFileName );
    dataRepository::loadTree( restartFileName );
  }

  {
    ProblemManager problemManager( "Problem", nullptr );

    problemManager.InitializePythonInterpreter();
    problemManager.ParseCommandLineInput();

    if( !problemManager.getSchemaFileName().empty() )
    {
      problemManager.GenerateDocumentation();
    }
    else
    {
      problemManager.ParseInputFile();

      problemManager.ProblemSetup();

      if( restart )
      {
        problemManager.ReadRestartOverwrite();
      }

      MpiWrapper::Barrier( MPI_COMM_GEOSX );
      GEOSX_LOG_RANK_0( "Running simulation" );

      gettimeofday( &tim, nullptr );
      const real64 t_initialize = tim.tv_sec + ( tim.tv_usec / 1000000.0 );

      problemManager.RunSimulation();

      gettimeofday( &tim, nullptr );
      const real64 t_run = tim.tv_sec + ( tim.tv_usec / 1000000.0 );

      GEOSX_LOG_RANK_0( "\ninit time = "
                        << std::setprecision( 5 ) << t_initialize - t_start
                        << "s, run time = " << t_run - t_initialize << "s" );
    }

    problemManager.ClosePythonInterpreter();
  }

  GEOSX_MARK_FUNCTION_END;

  basicCleanup();

  return 0;
}
