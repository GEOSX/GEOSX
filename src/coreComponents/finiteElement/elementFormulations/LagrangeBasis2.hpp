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

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_
/**
 * @file LagrangeBasis2.hpp
 */

#include "common/DataTypes.hpp"

/**
 * This class contains the implementation for a second order (quadratic) Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o-------------o-------------o  ---> xi
 *  Index:         0             1             2
 *  Coordinate:   -1             0             1
 *
 */
class LagrangeBasis2
{
  /**
   * @brief The value of the basis function for support point 0.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value0( constexpr real64 xi )
  {
    constexpr real64 xi_div2 = 0.5 * xi;
    return -xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief The value of the basis function for support point 1.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value1( constexpr real64 xi )
  {
    return 1.0 - xi * xi;
  }

  /**
   * @brief The value of the basis function for support point 2.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value2( constexpr real64 xi )
  {
    constexpr real64 xi_div2 = 0.5 * xi;
    return xi_div2 + xi_div2 * xi;
  }

  /**
   * @brief The gradient of the basis function for support point 0 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient0( constexpr real64 xi )
  {
    return -0.5 + xi;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient1( constexpr real64 xi )
  {
    return -2 * xi;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient2( constexpr real64 xi )
  {
    return 0.5 + xi;
  }

  /**
   * @class TensorProduct3D
   *
   *                                                                  ____________________
   *                                                                 |Node   xi0  xi1  xi2|
   *                                                                 |=====  ===  ===  ===|
   *                                                                 |  0    -1   -1   -1 |
   *                                                                 |  1     0   -1   -1 |
   *                                                                 |  2     1   -1   -1 |
   *              24              25               26                |  3    -1    0   -1 |
   *                o--------------o--------------o                  |  4     0    0   -1 |
   *               /.                            /|                  |  5     1    0   -1 |
   *              / .                           / |                  |  6    -1    1   -1 |
   *          21 o  .           o 22        23 o  |                  |  7     0    1   -1 |
   *            /   .                         /   |                  |  8     1    1   -1 |
   *           /    .         19             /    |                  |  9    -1   -1    0 |
   *       18 o--------------o--------------o 20  |                  | 10     0   -1    0 |
   *          |     o              o        |     o                  | 11     1   -1    0 |
   *          |     .15             16      |     |17                | 12    -1    0    0 |
   *          |     .                       |     |                  | 13     0    0    0 |
   *          |  o  .           o           |  o  |                  | 14     1    0    0 |
   *          |   12.            13         |   14|                  | 15    -1    1    0 |
   *          |     .                       |     |                  | 16     0    1    0 |
   *        9 o     .        o 10           o 11  |                  | 17     1    1    0 |
   *          |     o..............o........|.....o                  | 18    -1   -1    1 |
   *          |    , 6              7       |    / 8                 | 19     0   -1    1 |
   *          |   ,                         |   /                    | 20     1   -1    1 |
   *          |  o              o           |  o         xi2         | 21    -1    0    1 |
   *          | , 3              4          | / 5        |           | 22     0    0    1 |
   *          |,                            |/           | / xi1     | 23     1    0    1 |
   *          o--------------o--------------o            |/          | 24    -1    1    1 |
   *         0                1              2           o----- xi0  | 25     0    1    1 |
   *                                                                 | 26     1    1    1 |
   *                                                                 |____________________|
   *
   */
  struct TensorProduct3D
  {
    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,1)
     * @param j The index in the xi1 direction (0,1)
     * @param k The index in the xi2 direction (0,1)
     * @return The linear index of the support/quadrature point (0-7)
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static int
    linearIndex( constexpr int i,
                 constexpr int j,
                 constexpr int k )
    {
      return i + 3 * j + 9 * k;
    }

    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param a The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     * @param i2 The Cartesian index of the support point in the xi2 direction.
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static void
    multiIndex( constexpr int linearIndex,
                int & i0,
                int & i1,
                int & i2 )
    {
      i2 = ( linearIndex * 29 ) >> 8;
      //i2 = a/9;

      i1 = ( ( linearIndex * 22 ) >> 6 ) - i2 * 3;
      //i1 = a/3 - i2 * 3;

      i0 = linearIndex - i1 * 3 - i2 * 9;
    }
  };
};

#endif /* GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS2_HPP_ */
