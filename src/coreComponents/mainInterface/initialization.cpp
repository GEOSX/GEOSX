/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-LiCense-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "initialization.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "common/Path.hpp"
#include "LvArray/src/system.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "common/MpiWrapper.hpp"

// TPL includes
#include <optionparser.h>
#include <umpire/ResourceManager.hpp>

#if defined( GEOSX_USE_CALIPER )
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#endif

// System includes
#include <iomanip>

#if defined( GEOSX_USE_OPENMP )
#include <omp.h>
#endif

#if defined( GEOSX_USE_CUDA )
#include <cuda.h>
#endif

#include <fenv.h>

namespace geosx
{

/**
 * @class Arg a class inheriting from option::Arg that can parse a command line argument.
 */
struct Arg : public option::Arg
{
  /**
   * @brief Parse an unknown option. Unknown options aren't supported so this throws an error.
   * @param option the option to parse.
   * @return option::ARG_ILLEGAL.
   */
  static option::ArgStatus unknown( option::Option const & option, bool )
  {
    GEOSX_LOG_RANK( "Unknown option: " << option.name );
    return option::ARG_ILLEGAL;
  }

  /**
   * @brief Parse a non-empty string option.
   * @param option the option to parse.
   * @return option::ARK_OK if the parse was successful, option::ARG_ILLEGAL otherwise.
   */
  static option::ArgStatus nonEmpty( const option::Option & option, bool )
  {
    if((option.arg != nullptr) && (option.arg[0] != 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK( "Error: " << option.name << " requires a non-empty argument!" );
    return option::ARG_ILLEGAL;
  }

  /**
   * @brief Parse a numeric string option.
   * @param option the option to parse.
   * @return option::ARK_OK if the parse was successful, option::ARG_ILLEGAL otherwise.
   */
  static option::ArgStatus numeric( const option::Option & option, bool )
  {
    char * endptr = nullptr;
    if((option.arg != nullptr) && strtol( option.arg, &endptr, 10 )) {}
    if((endptr != option.arg) && (*endptr == 0))
    {
      return option::ARG_OK;
    }

    GEOSX_LOG_RANK( "Error: " << option.name << " requires a long-int argument!" );
    return option::ARG_ILLEGAL;
  }
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr< CommandLineOptions > parseCommandLineOptions( int argc, char * * argv )
{
  std::unique_ptr< CommandLineOptions > commandLineOptions = std::make_unique< CommandLineOptions >();

  // Set the options structs and parse
  enum optionIndex
  {
    UNKNOWN,
    HELP,
    INPUT,
    RESTART,
    XPAR,
    YPAR,
    ZPAR,
    SCHEMA,
    NONBLOCKING_MPI,
    SUPPRESS_PINNED,
    PROBLEMNAME,
    OUTPUTDIR,
    TIMERS,
    SUPPRESS_MOVE_LOGGING,
  };

  const option::Descriptor usage[] =
  {
    { UNKNOWN, 0, "", "", Arg::unknown, "USAGE: geosx -i input.xml [options]\n\nOptions:" },
    { HELP, 0, "?", "help", Arg::None, "\t-?, --help" },
    { INPUT, 0, "i", "input", Arg::nonEmpty, "\t-i, --input, \t Input xml filename (required)" },
    { RESTART, 0, "r", "restart", Arg::nonEmpty, "\t-r, --restart, \t Target restart filename" },
    { XPAR, 0, "x", "xpartitions", Arg::numeric, "\t-x, --x-partitions, \t Number of partitions in the x-direction" },
    { YPAR, 0, "y", "ypartitions", Arg::numeric, "\t-y, --y-partitions, \t Number of partitions in the y-direction" },
    { ZPAR, 0, "z", "zpartitions", Arg::numeric, "\t-z, --z-partitions, \t Number of partitions in the z-direction" },
    { SCHEMA, 0, "s", "schema", Arg::nonEmpty, "\t-s, --schema, \t Name of the output schema" },
    { NONBLOCKING_MPI, 0, "b", "use-nonblocking", Arg::None, "\t-b, --use-nonblocking, \t Use non-blocking MPI communication" },
    { PROBLEMNAME, 0, "n", "name", Arg::nonEmpty, "\t-n, --name, \t Name of the problem, used for output" },
    { SUPPRESS_PINNED, 0, "s", "suppress-pinned", Arg::None, "\t-s, --suppress-pinned \t Suppress usage of pinned memory for MPI communication buffers" },
    { OUTPUTDIR, 0, "o", "output", Arg::nonEmpty, "\t-o, --output, \t Directory to put the output files" },
    { TIMERS, 0, "t", "timers", Arg::nonEmpty, "\t-t, --timers, \t String specifying the type of timer output." },
    { SUPPRESS_MOVE_LOGGING, 0, "", "suppress-move-logging", Arg::None, "\t--suppress-move-logging \t Suppress logging of host-device data migration" },
    { 0, 0, nullptr, nullptr, nullptr, nullptr }
  };

  argc -= ( argc > 0 );
  argv += ( argc > 0 );
  option::Stats stats( usage, argc, argv );
  option::Option options[ 100 ];//stats.options_max];
  option::Option buffer[ 100 ];//stats.buffer_max];
  option::Parser parse( usage, argc, argv, options, buffer );

  // Handle special cases
  bool const noXML = options[INPUT].count() == 0 && options[SCHEMA].count() == 0;
  if( parse.error() || options[HELP] || (argc == 0) || noXML )
  {
    int columns = getenv( "COLUMNS" ) ? atoi( getenv( "COLUMNS" )) : 120;
    option::printUsage( fwrite, stdout, usage, columns );

    if( options[HELP] )
    {
      throw NotAnError();
    }

    GEOSX_THROW( "Bad command line arguments.", InputError );
  }

  // Iterate over the remaining inputs
  for( int ii=0; ii<parse.optionsCount(); ++ii )
  {
    option::Option & opt = buffer[ii];
    switch( opt.index() )
    {
      case UNKNOWN:
      {}
      break;
      case HELP:
      {}
      break;
      case INPUT:
      {
        commandLineOptions->inputFileNames.emplace_back( opt.arg );
      }
      break;
      case RESTART:
      {
        commandLineOptions->restartFileName = opt.arg;
        commandLineOptions->beginFromRestart = 1;
      }
      break;
      case XPAR:
      {
        commandLineOptions->xPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case YPAR:
      {
        commandLineOptions->yPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case ZPAR:
      {
        commandLineOptions->zPartitionsOverride = std::stoi( opt.arg );
        commandLineOptions->overridePartitionNumbers = 1;
      }
      break;
      case NONBLOCKING_MPI:
      {
        commandLineOptions->useNonblockingMPI = true;
      }
      break;
      case SUPPRESS_PINNED:
      {
        commandLineOptions->suppressPinned = true;
      }
      break;
      case SCHEMA:
      {
        commandLineOptions->schemaName = opt.arg;
      }
      break;
      case PROBLEMNAME:
      {
        commandLineOptions->problemName = opt.arg;
      }
      break;
      case OUTPUTDIR:
      {
        commandLineOptions->outputDirectory = opt.arg;
      }
      break;
      case TIMERS:
      {
        commandLineOptions->timerOutput = opt.arg;
      }
      break;
      case SUPPRESS_MOVE_LOGGING:
      {
        commandLineOptions->suppressMoveLogging = true;
      }
      break;
    }
  }

  if( commandLineOptions->problemName == "" && options[INPUT].count() > 0 )
  {
    string & inputFileName = commandLineOptions->inputFileNames[0];
    if( inputFileName.length() > 4 && inputFileName.substr( inputFileName.length() - 4, 4 ) == ".xml" )
    {
      string::size_type start = inputFileName.find_last_of( '/' ) + 1;
      if( start >= inputFileName.length())
      {
        start = 0;
      }
      commandLineOptions->problemName.assign( inputFileName, start, inputFileName.length() - 4 - start );
    }
    else
    {
      commandLineOptions->problemName.assign( inputFileName );
    }
  }

  return commandLineOptions;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine )
{
  setupEnvironment( argc, argv );
  setupLAI();

  if( parseCommandLine )
  {
    return parseCommandLineOptions( argc, argv );
  }
  else
  {
    return std::make_unique< CommandLineOptions >();
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void basicCleanup()
{
  finalizeLAI();
  cleanupEnvironment();
}


} // namespace geosx
