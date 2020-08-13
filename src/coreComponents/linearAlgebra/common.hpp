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
 * @file common.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_

#include "common/DataTypes.hpp"

/**
 * Whether to check preconditions at runtime in LAI functions
 */
#define GEOSX_LAI_RUNTIME_ASSERT 1
#if GEOSX_LAI_RUNTIME_ASSERT

  /**
 * Assert expression is true
 * @param expr expression
 */
  #define GEOSX_LAI_ASSERT( expr ) GEOSX_ERROR_IF( !( expr ), "" )

  /**
 * Assert expression and output message if false
 * @param expr expression
 * @param msg message
 */
  #define GEOSX_LAI_ASSERT_MSG( expr, msg ) GEOSX_ERROR_IF( !( expr ), msg )

  /**
 * Assert lhs equals rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_EQ( lhs, rhs ) GEOSX_ERROR_IF_NE( lhs, rhs )

  /**
 * Assert lhs not equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_NE( lhs, rhs ) GEOSX_ERROR_IF_EQ( lhs, rhs )

  /**
 * Assert lhs greater than rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_GT( lhs, rhs ) GEOSX_ERROR_IF_GE( rhs, lhs )

  /**
 * Assert lhs greater than or equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_GE( lhs, rhs ) GEOSX_ERROR_IF_GT( rhs, lhs )
#else

  /**
 * Assert expression is true
 * @param expr expression
 */
  #define GEOSX_LAI_ASSERT( expr ) GEOSX_ASSERT( expr )

  /**
 * Assert expression and output message if false
 * @param expr expression
 * @param msg message
 */
  #define GEOSX_LAI_ASSERT_MSG( expr, msg ) GEOSX_ASSERT_MSG( expr, msg )

  /**
 * Assert lhs equals rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_EQ( lhs, rhs ) GEOSX_ASSERT_EQ( lhs, rhs )

  /**
 * Assert lhs not equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_NE( lhs, rhs ) GEOSX_ASSERT_NE( lhs, rhs )

  /**
 * Assert lhs greater rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_GT( lhs, rhs ) GEOSX_ASSERT_GT( lhs, rhs )

  /**
 * Assert lhs greater than or equal to rhs
 * @param lhs left hand side
 * @param rhs right hand side
 */
  #define GEOSX_LAI_ASSERT_GE( lhs, rhs ) GEOSX_ASSERT_GE( lhs, rhs )
#endif

/**
 * Macro for checking and reporting error codes from TPL packages
 * @param call call to check function
 */
#define GEOSX_LAI_CHECK_ERROR( call )                               \
  do                                                                \
  {                                                                 \
    auto const ierr = call;                                         \
    GEOSX_ERROR_IF_NE_MSG( ierr, 0, "Error in call to " << #call ); \
  } while( false )

/**
 * Macro for checking and reporting non-negative error codes from TPL packages
 * @param call call to check function
 */
#define GEOSX_LAI_CHECK_ERROR_NNEG( call )                          \
  do                                                                \
  {                                                                 \
    auto const ierr = call;                                         \
    GEOSX_ERROR_IF_GT_MSG( 0, ierr, "Error in call to " << #call ); \
  } while( false )

namespace geosx
{
/**
 *  Enumeration of available output formats for LAI objects
 */
enum class LAIOutputFormat
{
  NATIVE_ASCII,
  NATIVE_BINARY,
  MATLAB_ASCII,
  MATLAB_BINARY,
  MATRIX_MARKET
};

/**
 * Permutation info for row major and column major layouts
 */
struct MatrixLayout
{
  /// typedef LAI row major permutation consistent with RAJA order
  using ROW_MAJOR_PERM = RAJA::PERM_IJ;

  /// typedef LAI col major permutation consistent with RAJA order
  using COL_MAJOR_PERM = RAJA::PERM_JI;

  /// row major data unit stride dim
  constexpr static int const ROW_MAJOR =
    LvArray::typeManipulation::getStrideOneDimension( ROW_MAJOR_PERM {} );

  /// column major unit stride dim
  constexpr static int const COL_MAJOR =
    LvArray::typeManipulation::getStrideOneDimension( COL_MAJOR_PERM {} );
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_COMMON_HPP_
