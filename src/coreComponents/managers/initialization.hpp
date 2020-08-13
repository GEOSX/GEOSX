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

#ifndef GEOSX_MANAGERS_INITIALIZATION_HPP_
#define GEOSX_MANAGERS_INITIALIZATION_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "mpiCommunications/MpiWrapper.hpp"

// TPL includes
#ifdef GEOSX_USE_CALIPER
  #include <adiak.hpp>
#endif

namespace geosx
{
/**
 * CommandLineOptions class containing the parsed command line options.
 */
struct CommandLineOptions
{
  /// The path to the input xml.
  std::string inputFileName;

  /// True iff restarting from the middle of an existing run.
  bool beginFromRestart = false;

  /// The path to the restart file, if specified.
  std::string restartFileName;

  /// The number of partitions in the x direction.
  integer xPartitionsOverride;

  /// The number of partitions in the y direction.
  integer yPartitionsOverride;

  /// The number of partitions in the z direction.
  integer zPartitionsOverride;

  /// True if using the partition override.
  integer overridePartitionNumbers = false;

  /// True if processing mpi communications in any order.
  /// But leads to non-reproducible results.
  integer useNonblockingMPI = false;

  /// True iff supress the use of pinned memory buffers
  /// ( if available ) for MPI communication.
  /// Generally only used by the integration tests.
  integer suppressPinned = false;

  /// The name of the schema.
  std::string schemaName;

  /// The name of the problem being run.
  std::string problemName;

  /// The directory to put all output.
  std::string outputDirectory = ".";

  /// The string used to initialize caliper.
  std::string timerOutput = "";

  /// Suppress logging of host-device data migration.
  integer suppressMoveLogging = false;
};

/**
 * @brief Perform the basic GEOSX initialization and optionally parse the command line input.
 * @param [in] argc The number of command line arguments.
 * @param [in,out] argv The command line arguments.
 * @param [in] parseCommandLine True iff the command line options should be parsed.
 */
void
basicSetup( int argc, char * argv[], bool const parseCommandLine = false );

/**
 * @brief @return a struct containing all the parsed command line options.
 */
CommandLineOptions const &
getCommandLineOptions();

/**
 * @brief Override the input file name, useful only for tests.
 * @param inputFileName new input file name
 */
void
overrideInputFileName( std::string const & inputFileName );

/**
 * @brief Perform the basic GEOSX cleanup.
 */
void
basicCleanup();

/**
 * @brief Initialize the logger.
 */
void
setupLogger();

/**
 * @brief Finalize the logger.
 */
void
finalizeLogger();

/**
 * @brief Setup the LvArray library. This initializes signal handling
 *        and the floating point environment.
 */
void
setupCXXUtils();

/**
 * @brief Setup MKL if in use.
 */
void
setupMKL();

/**
 * @brief Setup OpenMP.
 */
void
setupOpenMP();

/**
 * @brief Setup MPI.
 * @param [in] argc the number of command line arguments.
 * @param [in,out] argv the command line arguments.
 */
void
setupMPI( int argc, char * argv[] );

/**
 * @brief Finalize MPI.
 */
void
finalizeMPI();

/**
 * @brief Compute the sum, mean, min, and max of @p value across ranks and push
 *        them into Adiak using @p name.
 * @tparam T The type of @p value.
 * @param name The name to use when adding the stats to Adiak.
 * @param value The value to compute the statistics of.
 */
template< typename T >
void
pushStatsIntoAdiak( std::string const & name, T const value )
{
#if defined( GEOSX_USE_CALIPER )
  T const total = MpiWrapper::Sum( value );
  adiak::value( name + " sum", total );
  adiak::value( name + " mean", double( total ) / MpiWrapper::Comm_size() );
  adiak::value( name + " min", MpiWrapper::Min( value ) );
  adiak::value( name + " max", MpiWrapper::Max( value ) );
#else
  GEOSX_UNUSED_VAR( name );
  GEOSX_UNUSED_VAR( value );
#endif
}

}  // namespace geosx

#endif  // GEOSX_MANAGERS_INITIALIZATION_HPP_
