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

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_

/**
 * @file LagrangeBasis1.hpp
 */

#include "common/DataTypes.hpp"

namespace geosx
{
namespace finiteElement
{
/**
 * This class contains the implementation for a first order (linear) Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o-------------o  ---> xi
 *  Index:         0             1
 *  Coordinate:   -1             1
 *
 */
class LagrangeBasis1
{
public:
  /// The number of support points for the basis
  constexpr static localIndex numSupportPoints = 2;

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  parentSupportCoord( const localIndex supportPointIndex )
  {
    return -1.0 + 2.0 * ( supportPointIndex & 1 );
  }

  /**
   * @brief The value of the basis function for a support point evaluated at a
   *   point along the axes.
   * @param index The index of the support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of basis function.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value( const int index, const real64 xi )
  {
    return 0.5 + 0.5 * xi * parentSupportCoord( index );
  }

  /**
   * @brief The value of the basis function for the 0 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value0( const real64 xi )
  {
    return 0.5 - 0.5 * xi;
  }

  /**
   * @brief The value of the basis function for the 1 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  value1( const real64 xi )
  {
    return 0.5 + 0.5 * xi;
  }

  /**
   * @brief The gradient of the basis function for a support point evaluated at
   *   a point along the axes.
   * @param index The index of the support point associated with the basis
   *   function.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient( const int index, const real64 xi )
  {
    GEOSX_UNUSED_VAR( xi );
    return 0.5 * parentSupportCoord( index );
  }

  /**
   * @brief The gradient of the basis function for support point 0 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function (-0.5)
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient0( const real64 xi )
  {
    GEOSX_UNUSED_VAR( xi );
    return -0.5;
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function (0.5)
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64
  gradient1( const real64 xi )
  {
    GEOSX_UNUSED_VAR( xi );
    return 0.5;
  }

  /**
   * @struct TensorProduct3D
   *
   * A 3-dimensional basis formed from the tensor product of the 1d basis.
   *
   *                  6                   7                       ____________________
   *                   o-----------------o                       |Node   xi0  xi1  xi2|
   *                  /.                /|                       |=====  ===  ===  ===|
   *                 / .               / |                       | 0     -1   -1   -1 |
   *              4 o-----------------o 5|                       | 1      1   -1   -1 |
   *                |  .              |  |                       | 2     -1    1   -1 |
   *                |  .              |  |                       | 3      1    1   -1 |
   *                |  .              |  |                       | 4     -1   -1    1 |
   *                |  .              |  |                       | 5      1   -1    1 |
   *                |2 o..............|..o 3       xi2           | 6     -1    1    1 |
   *                | ,               | /          |             | 7      1    1    1 |
   *                |,                |/           | / xi1       |____________________|
   *                o-----------------o            |/
   *               0                   1           ------ xi0
   *
   */
  struct TensorProduct3D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 8;

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
    linearIndex( const int i, const int j, const int k )
    {
      return i + 2 * j + 4 * k;
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
    multiIndex( const int linearIndex,
                int & i0,
                int & i1,
                int & i2 )
    {
      i0 = ( linearIndex & 1 );
      i1 = ( linearIndex & 2 ) >> 1;
      i2 = ( linearIndex & 4 ) >> 2;
    }

    /**
     * @brief The value of the basis function for a support point evaluated at a
     *   point along the axes.
     *
     * @param coords The coordinates (in the parent frame) at which to evaluate the basis
     * @param N Array to hold the value of the basis functions at each support point.
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    static void
    value( const real64 ( &coords )[3], real64 ( &N )[numSupportPoints] )
    {
      for( int a = 0; a < 2; ++a )
      {
        for( int b = 0; b < 2; ++b )
        {
          for( int c = 0; c < 2; ++c )
          {
            const int lindex =
              LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );
            N[lindex] = LagrangeBasis1::value( a, coords[0] ) *
              LagrangeBasis1::value( b, coords[1] ) *
              LagrangeBasis1::value( c, coords[2] );
          }
        }
      }
    }

    /**
     * @brief The parent coordinates for a support point in the xi0 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64
    parentCoords0( localIndex const linearIndex )
    {
      return -1.0 + 2.0 * ( linearIndex & 1 );
    }

    /**
     * @brief The parent coordinates for a support point in the xi1 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64
    parentCoords1( localIndex const linearIndex )
    {
      return -1.0 + ( linearIndex & 2 );
    }

    /**
     * @brief The parent coordinates for a support point in the xi2 direction.
     * @param linearIndex The linear index of the support point
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static real64
    parentCoords2( localIndex const linearIndex )
    {
      return -1.0 + 0.5 * ( linearIndex & 4 );
    }
  };
};

}  // namespace finiteElement
}  // namespace geosx

#endif /* GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_LAGRANGEBASIS1_HPP_ */
