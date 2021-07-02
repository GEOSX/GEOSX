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
#include "mainInterface/initialization.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include <unistd.h>

// System includes
#include <chrono>

using namespace geosx;


int main( int argc, char *argv[] )
{
  try
  {
    std::chrono::system_clock::time_point const startTime = std::chrono::system_clock::now();

    std::unique_ptr< CommandLineOptions > commandLineOptions = basicSetup( argc, argv, true );

    std::chrono::system_clock::duration initTime;
    std::chrono::system_clock::duration runTime;
    {
      GeosxState state( std::move( commandLineOptions ) );

      bool const problemToRun = state.initializeDataRepository();
      if( problemToRun )
      {
        state.applyInitialConditions();
        state.run();
        LVARRAY_WARNING_IF( state.getState() != State::COMPLETED, "Simulation exited early." );
      }

      initTime = state.getInitTime();
      runTime = state.getRunTime();
    }

    basicCleanup();

    std::chrono::system_clock::duration const totalTime = std::chrono::system_clock::now() - startTime;

    GEOSX_LOG_RANK_0( "total time          " << durationToString( totalTime ) );
    GEOSX_LOG_RANK_0( "initialization time " << durationToString( initTime ) );
    GEOSX_LOG_RANK_0( "run time            " << durationToString( runTime ) );

    return 0;
  }
  // A NotAnError is thrown if "-h" or "--help" option is used.
  catch( NotAnError const & )
  {
    return 0;
  }
  catch( std::exception const & e )
  {
    GEOSX_LOG( e.what() );
    LvArray::system::callErrorHandler();
    std::abort();
  }
}
