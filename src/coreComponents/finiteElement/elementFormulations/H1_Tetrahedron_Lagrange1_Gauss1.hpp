/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
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

/**
 * @file H1_Tetrahedron_Lagrange1_Gauss1.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1

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

  GEOSX_HOST_DEVICE
  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  virtual localIndex getNumSupportPoints() const override
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
  static void calcN( localIndex const q,
                     real64 ( &N )[numNodes] );

  /**
   * @brief Calculate the shape functions derivatives wrt the physical
   *   coordinates.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   * @return The determinant of the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           real64 const (&X)[numNodes][3],
                           real64 ( &gradN )[numNodes][3] );

  /**
   * @brief Calculate the integration weights for a quadrature point.
   * @param q Index of the quadrature point.
   * @param X Array containing the coordinates of the support points.
   * @return The product of the quadrature rule weight and the determinate of
   *   the parent/physical transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const q,
                                             real64 const (&X)[numNodes][3] );

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param q The quadrature point index in 3d space.
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   * @return The determinant of the Jacobian transformation matrix.
   */
  GEOSX_HOST_DEVICE
  static real64 invJacobianTransformation( int const q,
                                           real64 const (&X)[numNodes][3],
                                           real64 ( & J )[3][3] )
  {
    GEOSX_UNUSED_VAR( q, X );
//    jacobianTransformation( q, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }
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
  static real64 determinantJacobianTransformation( real64 const (&X)[numNodes][3] );

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  determinantJacobianTransformation( real64 const (&X)[numNodes][3] )
{
  return ( X[1][0] - X[0][0] )*( ( X[2][1] - X[0][1] )*( X[3][2] - X[0][2] ) - ( X[3][1] - X[0][1] )*( X[2][2] - X[0][2] ) )
         + ( X[1][1] - X[0][1] )*( ( X[3][0] - X[0][0] )*( X[2][2] - X[0][2] ) - ( X[2][0] - X[0][0] )*( X[3][2] - X[0][2] ) )
         + ( X[1][2] - X[0][2] )*( ( X[2][0] - X[0][0] )*( X[3][1] - X[0][1] ) - ( X[3][0] - X[0][0] )*( X[2][1] - X[0][1] ) );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Tetrahedron_Lagrange1_Gauss1::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{
  GEOSX_UNUSED_VAR( q );

  // single quadrature point (centroid), i.e.  r = s = t = 1/4
  N[0] = 0.25; // N0 = 1 - r - s - t
  N[1] = 0.25; // N1 = r
  N[2] = 0.25; // N2 = s
  N[3] = 0.25; // N3 = t
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             real64 (& gradN)[numNodes][3] )
{
  GEOSX_UNUSED_VAR( q );

  gradN[0][0] =  X[1][1]*( X[3][2] - X[2][2] ) - X[2][1]*( X[3][2] - X[1][2] ) + X[3][1]*( X[2][2] - X[1][2] );
  gradN[0][1] = -X[1][0]*( X[3][2] - X[2][2] ) + X[2][0]*( X[3][2] - X[1][2] ) - X[3][0]*( X[2][2] - X[1][2] );
  gradN[0][2] =  X[1][0]*( X[3][1] - X[2][1] ) - X[2][0]*( X[3][1] - X[1][1] ) + X[3][0]*( X[2][1] - X[1][1] );

  gradN[1][0] = -X[0][1]*( X[3][2] - X[2][2] ) + X[2][1]*( X[3][2] - X[0][2] ) - X[3][1]*( X[2][2] - X[0][2] );
  gradN[1][1] =  X[0][0]*( X[3][2] - X[2][2] ) - X[2][0]*( X[3][2] - X[0][2] ) + X[3][0]*( X[2][2] - X[0][2] );
  gradN[1][2] = -X[0][0]*( X[3][1] - X[2][1] ) + X[2][0]*( X[3][1] - X[0][1] ) - X[3][0]*( X[2][1] - X[0][1] );

  gradN[2][0] =  X[0][1]*( X[3][2] - X[1][2] ) - X[1][1]*( X[3][2] - X[0][2] ) + X[3][1]*( X[1][2] - X[0][2] );
  gradN[2][1] = -X[0][0]*( X[3][2] - X[1][2] ) + X[1][0]*( X[3][2] - X[0][2] ) - X[3][0]*( X[1][2] - X[0][2] );
  gradN[2][2] =  X[0][0]*( X[3][1] - X[1][1] ) - X[1][0]*( X[3][1] - X[0][1] ) + X[3][0]*( X[1][1] - X[0][1] );

  gradN[3][0] = -X[0][1]*( X[2][2] - X[1][2] ) + X[1][1]*( X[2][2] - X[0][2] ) - X[2][1]*( X[1][2] - X[0][2] );
  gradN[3][1] =  X[0][0]*( X[2][2] - X[1][2] ) - X[1][0]*( X[2][2] - X[0][2] ) + X[2][0]*( X[1][2] - X[0][2] );
  gradN[3][2] = -X[0][0]*( X[2][1] - X[1][1] ) + X[1][0]*( X[2][1] - X[0][1] ) - X[2][0]*( X[1][1] - X[0][1] );

  real64 detJ = determinantJacobianTransformation( X );
  real64 factor = 1.0 / ( detJ );

  for( int i = 0; i < numNodes; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      gradN[i][j] *= factor;
    }
  }

  return detJ * weight;
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Tetrahedron_Lagrange1_Gauss1::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  GEOSX_UNUSED_VAR( q );

  real64 detJ =  determinantJacobianTransformation( X );

  return detJ * weight;
}

}
}

#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1TETRAHEDRONLAGRANGE1GAUSS1
