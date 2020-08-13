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
 * @file H1_Tetrahedron_Lagrange1_Gauss1.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_H1TETRAHEDRONLAGRANGE1GAUSS1
#define GEOSX_CORE_FINITEELEMENT_H1TETRAHEDRONLAGRANGE1GAUSS1

#include "FiniteElementBase.hpp"

namespace geosx
{
namespace finiteElement
{
/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal linear tetrahedron finite element with a 1-point
 * Gaussian quadrature rule.
 *
 *          3                                =====  ==  ==  ==
 *           +                               Node   r   s   t
 *           |\\_                            =====  ==  ==  ==
 *           ||  \_                          0      0   0   0
 *           | \   \_           t            1      1   0   1
 *           |  +__  \_         |   s        2      0   1   0
 *           | /2  \__ \_       |  /         3      0   0   1
 *           |/       \__\      | /          =====  ==  ==  ==
 *           +------------+     *------r
 *          0              1
 *
 */
class H1_Tetrahedron_Lagrange1_Gauss1 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 4;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 1;

  virtual ~H1_Tetrahedron_Lagrange1_Gauss1() override
  {}

  virtual localIndex
  getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  virtual localIndex
  getNumSupportPoints() const override
  {
    return numNodes;
  }

  /**
   * @brief Calculate shape functions values for each support point at a
   *   quadrature point.
   * @param q Index of the quadrature point.
   * @param N An array to pass back the shape function values for each support
   *   point.
   */
  GEOSX_HOST_DEVICE
  static void
  shapeFunctionValues( localIndex const q, real64 ( &N )[numNodes] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param dNdX Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64
  shapeFunctionDerivatives( localIndex const q,
                            real64 const ( &X )[numNodes][3],
                            real64 ( &dNdX )[numNodes][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64
  transformedQuadratureWeight( localIndex const q,
                               real64 const ( &X )[numNodes][3] );

private:
  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = 1.0 / 6.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /**
   * @brief Calculates the determinant of the Jacobian of the isoparametric
   *        mapping from the parent space to the physical space.
   * @param X Array containing the coordinates of the support points
   * @return determinant value
   */
  GEOSX_HOST_DEVICE
  static real64
  determinantJacobianTransformation( real64 const ( &X )[numNodes][3] );
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::determinantJacobianTransformation(
  real64 const ( &X )[numNodes][3] )
{
  return ( X[1][0] - X[0][0] ) *
    ( ( X[2][1] - X[0][1] ) * ( X[3][2] - X[0][2] ) -
      ( X[3][1] - X[0][1] ) * ( X[2][2] - X[0][2] ) ) +
    ( X[1][1] - X[0][1] ) *
    ( ( X[3][0] - X[0][0] ) * ( X[2][2] - X[0][2] ) -
      ( X[2][0] - X[0][0] ) * ( X[3][2] - X[0][2] ) ) +
    ( X[1][2] - X[0][2] ) *
    ( ( X[2][0] - X[0][0] ) * ( X[3][1] - X[0][1] ) -
      ( X[3][0] - X[0][0] ) * ( X[2][1] - X[0][1] ) );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Tetrahedron_Lagrange1_Gauss1::shapeFunctionValues( localIndex const q,
                                                      real64 ( &N )[numNodes] )
{
  GEOSX_UNUSED_VAR( q );

  // single quadrature point (centroid), i.e.  r = s = t = 1/4
  N[0] = 0.25;  // N0 = 1 - r - s - t
  N[1] = 0.25;  // N1 = r
  N[2] = 0.25;  // N2 = s
  N[3] = 0.25;  // N3 = t
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::shapeFunctionDerivatives(
  localIndex const q,
  real64 const ( &X )[numNodes][3],
  real64 ( &dNdX )[numNodes][3] )
{
  GEOSX_UNUSED_VAR( q );

  dNdX[0][0] = X[1][1] * ( X[3][2] - X[2][2] ) - X[2][1] * ( X[3][2] - X[1][2] ) +
    X[3][1] * ( X[2][2] - X[1][2] );
  dNdX[0][1] = -X[1][0] * ( X[3][2] - X[2][2] ) + X[2][0] * ( X[3][2] - X[1][2] ) -
    X[3][0] * ( X[2][2] - X[1][2] );
  dNdX[0][2] = X[1][0] * ( X[3][1] - X[2][1] ) - X[2][0] * ( X[3][1] - X[1][1] ) +
    X[3][0] * ( X[2][1] - X[1][1] );

  dNdX[1][0] = -X[0][1] * ( X[3][2] - X[2][2] ) + X[2][1] * ( X[3][2] - X[0][2] ) -
    X[3][1] * ( X[2][2] - X[0][2] );
  dNdX[1][1] = X[0][0] * ( X[3][2] - X[2][2] ) - X[2][0] * ( X[3][2] - X[0][2] ) +
    X[3][0] * ( X[2][2] - X[0][2] );
  dNdX[1][2] = -X[0][0] * ( X[3][1] - X[2][1] ) + X[2][0] * ( X[3][1] - X[0][1] ) -
    X[3][0] * ( X[2][1] - X[0][1] );

  dNdX[2][0] = X[0][1] * ( X[3][2] - X[1][2] ) - X[1][1] * ( X[3][2] - X[0][2] ) +
    X[3][1] * ( X[1][2] - X[0][2] );
  dNdX[2][1] = -X[0][0] * ( X[3][2] - X[1][2] ) + X[1][0] * ( X[3][2] - X[0][2] ) -
    X[3][0] * ( X[1][2] - X[0][2] );
  dNdX[2][2] = X[0][0] * ( X[3][1] - X[1][1] ) - X[1][0] * ( X[3][1] - X[0][1] ) +
    X[3][0] * ( X[1][1] - X[0][1] );

  dNdX[3][0] = -X[0][1] * ( X[2][2] - X[1][2] ) + X[1][1] * ( X[2][2] - X[0][2] ) -
    X[2][1] * ( X[1][2] - X[0][2] );
  dNdX[3][1] = X[0][0] * ( X[2][2] - X[1][2] ) - X[1][0] * ( X[2][2] - X[0][2] ) +
    X[2][0] * ( X[1][2] - X[0][2] );
  dNdX[3][2] = -X[0][0] * ( X[2][1] - X[1][1] ) + X[1][0] * ( X[2][1] - X[0][1] ) -
    X[2][0] * ( X[1][1] - X[0][1] );

  real64 detJ = determinantJacobianTransformation( X );
  real64 factor = 1.0 / ( detJ );

  for( int i = 0; i < numNodes; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      dNdX[i][j] *= factor;
    }
  }

  return detJ * weight;
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::transformedQuadratureWeight(
  localIndex const q,
  real64 const ( &X )[numNodes][3] )
{
  GEOSX_UNUSED_VAR( q );

  real64 detJ = determinantJacobianTransformation( X );

  return detJ * weight;
}

}  // namespace finiteElement
}  // namespace geosx

#endif  //GEOSX_CORE_FINITEELEMENT_H1TETRAHEDRONLAGRANGE1GAUSS1
