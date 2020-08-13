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

/* clang-format off */

/**
 * @file TimingMacros.hpp
 *
 * A collection of timing-related macros that wrap Caliper.
 */

#ifndef GEOSX_COMMON_TIMINGMACROS_HPP_
#define GEOSX_COMMON_TIMINGMACROS_HPP_

#include "common/GeosxConfig.hpp"
#include "common/GeosxMacros.hpp"

#ifdef GEOSX_USE_CALIPER
#include <caliper/cali.h>
#include <sys/time.h>
#include <string>
#include <iostream>

namespace timingHelpers
{
  inline std::string stripPF( char const * prettyFunction )
  {
    std::string input(prettyFunction);
    std::string::size_type const end = input.find_first_of( '(' );
    std::string::size_type const beg = input.find_last_of( ' ', end)+1;
    return input.substr( beg, end-beg );
  }
}

/// Mark a function or scope for timing with a given name
#define GEOSX_MARK_SCOPE(name) cali::Function __cali_ann##__LINE__(STRINGIZE_NX(name))

/// Mark a function for timing using a compiler-provided name
#define GEOSX_MARK_FUNCTION cali::Function __cali_ann##__func__(timingHelpers::stripPF(__PRETTY_FUNCTION__).c_str())

/// Mark the beginning of timed statement group
#define GEOSX_MARK_BEGIN(name) CALI_MARK_BEGIN(STRINGIZE(name))

/// Mark the end of timed statements group
#define GEOSX_MARK_END(name) CALI_MARK_END(STRINGIZE(name))

/// Mark the beginning of function, only useful when you don't want to or can't mark the whole function.
#define GEOSX_MARK_FUNCTION_BEGIN CALI_MARK_FUNCTION_BEGIN

/// Mark the end of function, only useful when you don't want to or can't mark the whole function.
#define GEOSX_MARK_FUNCTION_END CALI_MARK_FUNCTION_END

#else // GEOSX_USE_CALIPER

/// @cond DO_NOT_DOCUMENT
#define GEOSX_MARK_SCOPE(name)
#define GEOSX_MARK_FUNCTION_SCOPED
#define GEOSX_MARK_FUNCTION

#define GEOSX_MARK_BEGIN(name)
#define GEOSX_MARK_END(name)

#define GEOSX_MARK_FUNCTION_BEGIN
#define GEOSX_MARK_FUNCTION_END
/// @endcond

#endif // GEOSX_USE_CALIPER

/// Get current time of day as a floating point number of seconds in a variable @p time.
#ifdef GEOSX_USE_TIMERS
#define GEOSX_GET_TIME( time )                                                 \
  real64 time;                                                                 \
  do                                                                           \
  {                                                                            \
    timeval tim;                                                               \
    gettimeofday(&tim, nullptr);                                               \
    time = tim.tv_sec + (tim.tv_usec / 1000000.0);                             \
  } while (false)
#else
#define GEOSX_GET_TIME( time )
#endif

#endif // GEOSX_COMMON_TIMINGMACROS_HPP_
