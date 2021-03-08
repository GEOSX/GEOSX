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
#include "mpiCommunications/MpiWrapper.hpp"

// TPL includes
#include <optionparser.h>
#include <umpire/ResourceManager.hpp>

#if defined( GEOSX_USE_CALIPER )
#include <caliper/cali-manager.h>
#include <adiak.hpp>
#endif

// System includes
#include <iomanip>

#if defined( GEOSX_USE_MKL )
#include <mkl.h>
#endif

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
 * @brief For each Umpire::Allocator compute the total high water mark across all ranks
 *        and if using Adiak add statistics about the high water mark.
 */
static void addUmpireHighWaterMarks()
{
  umpire::ResourceManager & rm = umpire::ResourceManager::getInstance();

  // Get a list of all the allocators and sort it so that it's in the same order on each rank.
  std::vector< string > allocatorNames = rm.getAllocatorNames();
  std::sort( allocatorNames.begin(), allocatorNames.end() );

  // If each rank doesn't have the same number of allocators you can't aggregate them.
  std::size_t const numAllocators = allocatorNames.size();
  std::size_t const minNumAllocators = MpiWrapper::min( numAllocators );

  if( numAllocators != minNumAllocators )
  {
    GEOSX_WARNING( "Not all ranks have created the same number of umpire allocators, cannot compute high water marks." );
    return;
  }

  // Loop over the allocators.
  constexpr int MAX_NAME_LENGTH = 100;
  char allocatorNameBuffer[ MAX_NAME_LENGTH + 1 ];
  char allocatorNameMinCharsBuffer[ MAX_NAME_LENGTH + 1 ];
  for( string const & allocatorName : allocatorNames )
  {
    // Skip umpire internal allocators.
    if( allocatorName.rfind( "__umpire_internal", 0 ) == 0 )
      continue;

    GEOSX_ERROR_IF_GT( allocatorName.size(), MAX_NAME_LENGTH );

    memset( allocatorNameBuffer, '\0', sizeof( allocatorNameBuffer ) );
    memcpy( allocatorNameBuffer, allocatorName.data(), allocatorName.size() );

    memset( allocatorNameMinCharsBuffer, '\0', sizeof( allocatorNameMinCharsBuffer ) );

    // Make sure that each rank is looking at the same allocator.
    MpiWrapper::allReduce( allocatorNameBuffer, allocatorNameMinCharsBuffer, MAX_NAME_LENGTH, MPI_MIN, MPI_COMM_GEOSX );
    if( strcmp( allocatorNameBuffer, allocatorNameMinCharsBuffer ) != 0 )
    {
      GEOSX_WARNING( "Not all ranks have an allocator named " << allocatorNameBuffer << ", cannot compute high water mark." );
      continue;
    }

    // Get the total number of bytes allocated with this allocator across ranks.
    // This is a little redundant since
    std::size_t const mark = rm.getAllocator( allocatorName ).getHighWatermark();
    std::size_t const totalMark = MpiWrapper::sum( mark );
    GEOSX_LOG_RANK_0( "Umpire " << std::setw( 15 ) << allocatorName << " high water mark: " <<
                      std::setw( 9 ) << LvArray::system::calculateSize( totalMark ) );

    pushStatsIntoAdiak( allocatorName + " high water mark", mark );
  }
}

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
        commandLineOptions->inputFileName = opt.arg;
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

  if( commandLineOptions->problemName == "" )
  {
    string & inputFileName = commandLineOptions->inputFileName;
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

void printompinfo()
{
  int nthreads, tid;

/* Fork a team of threads giving them their own copies of variables */
#pragma omp parallel private(nthreads, tid) default(none)
  {

    /* Obtain thread number */
    tid = omp_get_thread_num();
    printf("Hello World from thread = %d\n", tid);

    /* Only master thread does this */
    if (tid == 0)
    {
      nthreads = omp_get_num_threads();
      printf("Number of threads = %d\n", nthreads);
    }

  } /* All threads join master thread and disband */

}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::unique_ptr< CommandLineOptions > basicSetup( int argc, char * argv[], bool const parseCommandLine )
{
  printompinfo();
  setupMPI( argc, argv );
  setupLogger();
  setupLvArray();
  setupOpenMP();  setupMKL();
  setupLAI( argc, argv );

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
  LvArray::system::resetSignalHandling();
  finalizeLAI();
  finalizeLogger();
  addUmpireHighWaterMarks();
  finalizeCaliper();
  finalizeMPI();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupLogger()
{
#ifdef GEOSX_USE_MPI
  logger::InitializeLogger( MPI_COMM_GEOSX );
#else
  logger::InitializeLogger();
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeLogger()
{
  logger::FinalizeLogger();
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupLvArray()
{
  LvArray::system::setErrorHandler( []()
  {
  #if defined( GEOSX_USE_MPI )
    int mpi = 0;
    MPI_Initialized( &mpi );
    if( mpi )
    {
      MPI_Abort( MPI_COMM_WORLD, EXIT_FAILURE );
    }
  #endif
    std::abort();
  } );

  LvArray::system::setSignalHandling( []( int const signal ) { LvArray::system::stackTraceHandler( signal, true ); } );

#if defined(GEOSX_USE_FPE)
  LvArray::system::setFPE();
#else
  LvArray::system::disableFloatingPointExceptions( FE_ALL_EXCEPT );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupMKL()
{
#ifdef GEOSX_USE_MKL
  GEOSX_LOG_RANK_0( "MKL max threads: " << mkl_get_max_threads() );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupOpenMP()
{
#ifdef GEOSX_USE_OPENMP
  GEOSX_LOG_RANK_0( "Max threads: " << omp_get_max_threads() );
#endif
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupMPI( int argc, char * argv[] )
{
  if( !MpiWrapper::initialized() )
  {
    MpiWrapper::init( &argc, &argv );
  }

  MPI_COMM_GEOSX = MpiWrapper::commDup( MPI_COMM_WORLD );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeMPI()
{
  MpiWrapper::commFree( MPI_COMM_GEOSX );
  MpiWrapper::finalize();
}

#if defined( GEOSX_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void setupCaliper( cali::ConfigManager & caliperManager,
                   CommandLineOptions const & commandLineOptions )
{
  caliperManager.add( commandLineOptions.timerOutput.c_str() );
  GEOSX_ERROR_IF( caliperManager.error(), "Caliper config error: " << caliperManager.error_msg() );
  caliperManager.start();

#if defined( GEOSX_USE_MPI )
  adiak::init( &MPI_COMM_GEOSX );
#else
  adiak::init( nullptr );
#endif

  GEOSX_WARNING_IF( !adiak::uid(), "Error getting the user info." );
  GEOSX_WARNING_IF( !adiak::launchdate(), "Error getting the launch date info." );
  GEOSX_WARNING_IF( !adiak::cmdline(), "Error getting the command line args." );
  GEOSX_WARNING_IF( !adiak::clustername(), "Error getting the clustername." );
  GEOSX_WARNING_IF( !adiak::walltime(), "Error getting the walltime." );
  GEOSX_WARNING_IF( !adiak::systime(), "Error getting the systime." );
  GEOSX_WARNING_IF( !adiak::cputime(), "Error getting the cputime." );

  adiak::value( "XML File", splitPath( commandLineOptions.inputFileName ).second );
  adiak::value( "Problem name", commandLineOptions.problemName );

  // MPI info
#if defined( GEOSX_USE_MPI )
  adiak::value( "MPI", "On" );
  adiak::value( "mpi ranks", MpiWrapper::commSize() );
#else
  adiak::value( "MPI", "Off" );
  adiak::value( "mpi ranks", 1 );
#endif

  // Build info
#if defined( __clang_version__ )
  adiak::value( "compiler", "clang" );
  adiak::value( "compiler version", adiak::version( "clang" __clang_version__ ) );
#elif defined( __INTEL_COMPILER )
  adiak::value( "compiler", "intel" );
  adiak::value( "compiler version", adiak::version( "intel" __INTEL_COMPILER ) );
#elif defined( __GNUC__ )
  adiak::value( "compiler", "gcc" );
  adiak::value( "compiler version", adiak::version( "gcc" __VERSION__ ) );
#else
  adiak::value( "compiler", "unknown" );
  adiak::value ( "compiler version", "unknown" );
#endif

  adiak::value( "build type", GEOSX_CMAKE_BUILD_TYPE );
  adiak::value( "compilation date", __DATE__ );

  // OpenMP info
#if defined( GEOSX_USE_OPENMP )
  std::int64_t const numThreads = omp_get_max_threads();
  adiak::value( "OpenMP", "On" );
#else
  std::int64_t const numThreads = 1;
  adiak::value( "OpenMP", "Off" );
#endif
  pushStatsIntoAdiak( "numThreads", numThreads );

  // CUDA info
  int cudaRuntimeVersion = 0;
  int cudaDriverVersion = 0;
#if defined( GEOSX_USE_CUDA )
  adiak::value( "CUDA", "On" );
  GEOSX_ERROR_IF_NE( cudaSuccess, cudaRuntimeGetVersion( &cudaRuntimeVersion ) );
  GEOSX_ERROR_IF_NE( cudaSuccess, cudaDriverGetVersion( &cudaDriverVersion ) );
#else
  adiak::value( "CUDA", "Off" );
#endif
  adiak::value( "CUDA runtime version", cudaRuntimeVersion );
  adiak::value( "CUDA driver version", cudaDriverVersion );
}
#endif // defined( GEOSX_USE_CALIPER )

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void finalizeCaliper()
{
#if defined( GEOSX_USE_CALIPER )
  adiak::fini();
#endif
}

} // namespace geosx
