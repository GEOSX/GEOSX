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
 * @file Logger.hpp
 */

#ifndef GEOSX_COMMON_LOGGER_HPP
#define GEOSX_COMMON_LOGGER_HPP

// Source incldes
#include "common/GeosxConfig.hpp"
#include "LvArray/src/Macros.hpp"

// System includes
#if defined( GEOSX_USE_MPI )
  #include <mpi.h>
#endif

/**
 * @brief Log a message on screen.
 * @details The expression to log must evaluate something that can be stream inserted.
 */
#define GEOSX_LOG( ... ) LVARRAY_LOG( __VA_ARGS__ )

/**
 * @brief Log an expression and its value on screen.
 * @details The expression to log must evaluate something that can be stream inserted.
 */
#define GEOSX_LOG_VAR( ... ) LVARRAY_LOG_VAR( __VA_ARGS__ )

/**
 * @brief Conditionally log a message on screen on rank 0.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_RANK_0_IF( EXP, msg )               \
  do                                                  \
  {                                                   \
    if( ::geosx::logger::internal::rank == 0 && EXP ) \
    {                                                 \
      std::ostringstream oss;                         \
      oss << msg;                                     \
      std::cout << oss.str() << std::endl;            \
    }                                                 \
  } while( false )

/**
 * @brief Log a message on screen on rank 0.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_RANK_0( msg ) GEOSX_LOG_RANK_0_IF( true, msg )

/**
 * @brief Conditionally log a message to the rank output stream.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_RANK_IF( EXP, msg )                                         \
  do                                                                          \
  {                                                                           \
    if( EXP )                                                                 \
    {                                                                         \
      std::ostringstream oss;                                                 \
      oss << "Rank " << ::geosx::logger::internal::rankString << ": " << msg; \
      *logger::internal::rankStream << oss.str() << std::endl;                \
    }                                                                         \
  } while( false )

/**
 * @brief Log a message to the rank output stream.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_RANK( msg ) GEOSX_LOG_RANK_IF( true, msg )

/**
 * @brief Log a variable/expression name and value on screen to the rank output stream.
 * @param var a variable or expression accessible from current scope that can be stream inserted
 */
#define GEOSX_LOG_RANK_VAR( var ) GEOSX_LOG_RANK( #var " = " << var )

/**
 * @brief Conditionally raise a hard error and terminate the program.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF( EXP, msg ) \
  LVARRAY_ERROR_IF(                \
    EXP,                           \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error and terminate the program.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR( msg ) GEOSX_ERROR_IF( true, msg )

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ASSERT_MSG( EXP, msg ) \
  LVARRAY_ASSERT_MSG(                \
    EXP,                             \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert a condition in debug builds.
 * @param EXP an expression that will be evaluated as a predicate
 */
#define GEOSX_ASSERT( EXP ) GEOSX_ASSERT_MSG( EXP, "" )

/**
 * @brief Conditionally report a warning.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_WARNING_IF( EXP, msg ) LVARRAY_WARNING_IF( EXP, msg )

/**
 * @brief Report a warning.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_WARNING( msg ) LVARRAY_WARNING( msg )

/**
 * @brief Conditionally log an info message.
 * @param EXP an expression that will be evaluated as a predicate
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_INFO_IF( EXP, msg ) LVARRAY_INFO_IF( EXP, msg )

/**
 * @brief Log an info message.
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_INFO( msg ) LVARRAY_INFO( msg )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_EQ_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_EQ_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if two values are equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_EQ( lhs, rhs ) GEOSX_ERROR_IF_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_NE_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_NE_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if two values are not equal.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_NE( lhs, rhs ) GEOSX_ERROR_IF_NE_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_GT_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_GT_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if one value compares greater than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_GT( lhs, rhs ) GEOSX_ERROR_IF_GT_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_GE_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_GE_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if one value compares greater than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_GE( lhs, rhs ) GEOSX_ERROR_IF_GE_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_LT_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_LT_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if one value compares less than the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_LT( lhs, rhs ) GEOSX_ERROR_IF_LT_MSG( lhs, rhs, "" )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ERROR_IF_LE_MSG( lhs, rhs, msg ) \
  LVARRAY_ERROR_IF_LE_MSG(                     \
    lhs,                                       \
    rhs,                                       \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Raise a hard error if one value compares less than or equal to the other.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ERROR_IF_LE( lhs, rhs ) GEOSX_ERROR_IF_LE_MSG( lhs, rhs, "" )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ASSERT_EQ_MSG( lhs, rhs, msg ) \
  LVARRAY_ASSERT_EQ_MSG(                     \
    lhs,                                     \
    rhs,                                     \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert that two values compare equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ASSERT_EQ( lhs, rhs ) GEOSX_ASSERT_EQ_MSG( lhs, rhs, "" )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ASSERT_NE_MSG( lhs, rhs, msg ) LVARRAY_ASSERT_NE_MSG( lhs, rhs, msg )

/**
 * @brief Assert that two values compare not equal in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ASSERT_NE( lhs, rhs ) LVARRAY_ASSERT_NE( lhs, rhs )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ASSERT_GT_MSG( lhs, rhs, msg ) \
  LVARRAY_ASSERT_GT_MSG(                     \
    lhs,                                     \
    rhs,                                     \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert that one value compares greater than the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ASSERT_GT( lhs, rhs ) GEOSX_ASSERT_GT_MSG( lhs, rhs, "" )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 * @param msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_ASSERT_GE_MSG( lhs, rhs, msg ) \
  LVARRAY_ASSERT_GE_MSG(                     \
    lhs,                                     \
    rhs,                                     \
    "***** Rank " << ::geosx::logger::internal::rankString << ": " << msg )

/**
 * @brief Assert that one value compares greater than or equal to the other in debug builds.
 * @param lhs expression to be evaluated and used as left-hand side in comparison
 * @param rhs expression to be evaluated and used as right-hand side in comparison
 */
#define GEOSX_ASSERT_GE( lhs, rhs ) GEOSX_ASSERT_GE_MSG( lhs, rhs, "" )

/**
 * @brief Macro used to turn on/off a function based on the log level.
 * @param[in] minLevel Minimum log level
 * @param[in] fn Function to filter
 */
#define GEOSX_LOG_LEVEL_FN( minLevel, fn ) \
  do                                       \
  {                                        \
    if( this->getLogLevel() >= minLevel )  \
    {                                      \
      fn;                                  \
    }                                      \
  } while( false )

/**
 * @brief Output messages based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_LEVEL( minLevel, msg ) \
  GEOSX_INFO_IF( this->getLogLevel() >= minLevel, msg );

/**
 * @brief Output messages (only on rank 0) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_LEVEL_RANK_0( minLevel, msg ) \
  GEOSX_LOG_RANK_0_IF( this->getLogLevel() >= minLevel, msg )

/**
 * @brief Output messages (with one line per rank) based on current Group's log level.
 * @param[in] minLevel minimum log level
 * @param[in] msg a message to log (any expression that can be stream inserted)
 */
#define GEOSX_LOG_LEVEL_BY_RANK( minLevel, msg ) \
  GEOSX_LOG_RANK_IF( this->getLogLevel() >= minLevel, msg )

namespace geosx
{
namespace logger
{
namespace internal
{
extern int rank;

extern std::string rankString;

extern int n_ranks;

extern std::ostream * rankStream;

#if defined( GEOSX_USE_MPI )
extern MPI_Comm comm;
#endif
}  // namespace internal

#if defined( GEOSX_USE_MPI )
/**
 * @brief Initialize the logger in a parallel build.
 * @param comm global MPI communicator
 * @param rank_output_dir output directory for rank log files
 */
void
InitializeLogger( MPI_Comm comm, const std::string & rank_output_dir = "" );
#endif

/**
 * @brief Initialize the logger in a serial build.
 * @param rank_output_dir output directory for rank log files
 */
void
InitializeLogger( const std::string & rank_output_dir = "" );

/**
 * @brief Finalize the logger and close the rank streams.
 */
void
FinalizeLogger();

}  // namespace logger

}  // namespace geosx

#endif /* GEOSX_COMMON_LOGGER_HPP */
