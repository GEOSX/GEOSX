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
 * @file Stopwatch.hpp
 */

#ifndef GEOSX_STOPWATCH_HPP
#define GEOSX_STOPWATCH_HPP

#include <chrono>

namespace geosx
{
/**
 * @brief Class defining a simple stopwatch for interval timing.
 * @note Current implementation relies on std::chrono::steady_clock
 */
class Stopwatch
{
public:
  /**
   * @brief Constructor.
   */
  Stopwatch()
  {
    zero();
  }

  /**
   * @brief Zero out the timer.
   */
  void
  zero()
  {
    m_start = std::chrono::steady_clock::now();
  }

  /**
   * @brief Return elapsed time in seconds since zero() was last called.
   * @return elapsed time
   */
  real64
  elapsedTime()
  {
    std::chrono::steady_clock::time_point const end =
      std::chrono::steady_clock::now();
    std::chrono::duration< real64 > const diff = end - m_start;
    return diff.count();
  }

private:
  /// Time point of the last timer restart
  std::chrono::steady_clock::time_point m_start;
};

}  // namespace geosx

#endif
