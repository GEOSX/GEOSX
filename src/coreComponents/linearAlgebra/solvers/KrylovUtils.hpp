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
 * @file KrylovUtils.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_

#include "codingUtilities/Utilities.hpp"

/// Tolerance for division by zero in Krylov solvers
#define GEOSX_KRYLOV_MIN_DIV ::geosx::NumericTraits< real64 >::eps

#ifndef GEOSX_KRYLOV_BREAKDOWN_IF_ZERO
  /**
 * @brief Exit solver iteration and report a breakdown if value too close to zero.
 * @param VAR the variable or expression
 */
  #define GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( VAR )                              \
    do                                                                       \
    {                                                                        \
      if( isZero( VAR, GEOSX_KRYLOV_MIN_DIV ) )                              \
      {                                                                      \
        GEOSX_LOG_LEVEL_RANK_0(                                              \
          1,                                                                 \
          "Breakdown in " << methodName() << ": " << #VAR << " = " << VAR ); \
        m_result.status = LinearSolverResult::Status::Breakdown;             \
        break;                                                               \
      }                                                                      \
    } while( false )
#endif

#endif  //GEOSX_LINEARALGEBRA_SOLVERS_KRYLOVUTILS_HPP_
