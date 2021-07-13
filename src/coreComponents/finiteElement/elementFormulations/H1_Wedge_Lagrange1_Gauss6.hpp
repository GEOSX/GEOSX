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
 * @file H1_Wedge_Lagrange1_Gauss6.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1WEDGELAGRANGE1GAUSS6
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1WEDGELAGRANGE1GAUSS6

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal bilinear wedge finite element with a 6-point Gaussian
 * quadrature rule. It is assumed that the indexing for the quadrature points
 * mirrors that of the nodes. Also note that the assumed node ordering is not
 * the standard right-hand-rule used in the literature.
 *
 *                5
 *                /\                                =====  ==  ==  ===
 *               / .\                               Node   r   s   xi
 *              /  . \                              =====  ==  ==  ===
 *            1/______\ 3                           0      0   0   -1
 *             |   .   |                            1      0   0    1
 *             |   .   |                            2      1   0   -1
 *             |   .   |                            3      1   0    1
 *             |   .   |                            4      0   1   -1
 *             |   4   |        xi                  5      0   1    1
 *             |  . .  |          |   s             =====  ==  ==  ===
 *             | .   . |          |  /
 *             |.____ .|          | /
 *             0       2          |/____ r
 *
 */
class H1_Wedge_Lagrange1_Gauss6 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 6;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 6;

  virtual ~H1_Wedge_Lagrange1_Gauss6() override
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
   *
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
    jacobianTransformation( q, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }


private:
  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = 1.0;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /// The factor specifying the location of the quadrature points
  /// relative to the origin of the element in the parent space.
  constexpr static real64 quadratureCrossSectionCoord = 1.0 / 6.0;

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space.
  constexpr static real64 quadratureLongitudinalCoord = 1.0 / 1.732050807568877293528;

  /**
   * @brief Calculates the index for support/quadrature points from the
   *        pair (indexT, indexL).
   * @param indexT The index in the triangular r-s cross section (0,1,2)
   * @param indexL The index in the xi (longitudinal) direction (0,1)
   * @return The linear index of the support/quadrature point (0-5)
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const indexT, T const indexL )
  {
    return 2 * indexT + indexL;
  }

  /**
   * @brief Calculate the parent coordinates for the r direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the r direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return 0.5 * ( a & 2 );
  }

  /**
   * @brief Calculate the parent coordinates for the s direction, given the
   *   linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the s direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return 0.25 * ( a & 4 );
  }

  /**
   * @brief Calculate the parent coordinates for the xi direction, given the
   *        linear index of a support point.
   * @param q The linear index of support point
   * @return parent coordinate in the xi direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 2 * ( a & 1 );
  }

  /**
   * @brief Calculate the parent coordinates for the r direction, given the
   *        linear index of a quadrature point.
   * @param a The linear index of quadrature point
   * @return parent coordinate in the r direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return quadratureCrossSectionCoord + 0.5  * parentCoords0( q );
  }

  /**
   * @brief Calculate the parent coordinates for the s direction, given the
   *        linear index of a quadrature point.
   * @param q The linear index of quadrature point
   * @return parent coordinate in the r direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return quadratureCrossSectionCoord + 0.5  * parentCoords1( q );
  }

  /**
   * @brief Calculate the parent coordinates for the xi direction, given the
   *        linear index of a quadrature point.
   * @param q The linear index of quadrature point
   * @return parent coordinate in the xi direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return parentCoords2( q ) * quadratureLongitudinalCoord;
  }

  /**
   * @brief Calculates the Jacobian matrix of the isoparametric mapping.
   * @param q The linear index of quadrature point
   * @param X Array containing the coordinates of the support points
   * @param J Array to store the Jacobian matrix
   */
  GEOSX_HOST_DEVICE
  static void jacobianTransformation( int const q,
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );

  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param q The linear index of quadrature point
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *             support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void
    applyJacobianTransformationToShapeFunctionsDerivatives( int const q,
                                                            real64 const ( &invJ )[3][3],
                                                            real64 ( &gradN )[numNodes][3] );

};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Wedge_Lagrange1_Gauss6::
  jacobianTransformation( int const q,
                          real64 const (&X)[numNodes][3],
                          real64 ( & J )[3][3] )
{
  real64 const r  = quadratureParentCoords0( q );
  real64 const s  = quadratureParentCoords1( q );
  real64 const xi = quadratureParentCoords2( q );

  real64 const psiTRI[3] = { 1.0 - r - s, r, s };
  real64 const psiLIN[2] = { 0.5 - 0.5*xi, 0.5 + 0.5*xi };
  constexpr real64 dpsiTRI[2][3] = { { -1.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0 } };
  constexpr real64 dpsiLIN[2] = { -0.5, 0.5 };

  for( localIndex a=0; a<3; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsiTRI[0][a] * psiLIN[b],
                                dpsiTRI[1][a] * psiLIN[b],
                                psiTRI[a] * dpsiLIN[b] };
      localIndex const nodeIndex = linearMap( a, b );
      for( int i = 0; i < 3; ++i )
      {
        for( int j = 0; j < 3; ++j )
        {
          J[i][j] = J[i][j] + dNdXi[ j ] * X[nodeIndex][i];
        }
      }
    }
  }
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Wedge_Lagrange1_Gauss6::
  applyJacobianTransformationToShapeFunctionsDerivatives( int const q,
                                                          real64 const ( &invJ )[3][3],
                                                          real64 (& gradN)[numNodes][3] )
{
  real64 const r  = quadratureParentCoords0( q );
  real64 const s  = quadratureParentCoords1( q );
  real64 const xi = quadratureParentCoords2( q );

  real64 const psiTRI[3] = { 1.0 - r - s, r, s };
  real64 const psiLIN[2] = { 0.5 - 0.5*xi, 0.5 + 0.5*xi };
  constexpr real64 dpsiTRI[2][3] = { { -1.0, 1.0, 0.0 }, { -1.0, 0.0, 1.0 } };
  constexpr real64 dpsiLIN[2] = { -0.5, 0.5 };

  for( localIndex a=0; a<3; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsiTRI[0][a] * psiLIN[b],
                                dpsiTRI[1][a] * psiLIN[b],
                                psiTRI[a] * dpsiLIN[b] };
      localIndex const nodeIndex = linearMap( a, b );
      for( int i = 0; i < 3; ++i )
      {
        gradN[nodeIndex][i] = 0.0;
        for( int j = 0; j < 3; ++j )
        {
          gradN[nodeIndex][i] = gradN[nodeIndex][i] + dNdXi[ j ] * invJ[j][i];
        }
      }
    }
  }
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Wedge_Lagrange1_Gauss6::
  calcN( localIndex const q,
         real64 (& N)[numNodes] )
{
  real64 const r  = quadratureParentCoords0( q );
  real64 const s  = quadratureParentCoords1( q );
  real64 const xi = quadratureParentCoords2( q );

  N[0] = 0.5*( 1.0 - r - s ) * ( 1.0 - xi );
  N[1] = 0.5*( 1.0 - r - s ) * ( 1.0 + xi );
  N[2] = 0.5* r * ( 1.0 - xi );
  N[3] = 0.5* r * ( 1.0 + xi );
  N[4] = 0.5* s * ( 1.0 - xi );
  N[5] = 0.5* s * ( 1.0 + xi );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Wedge_Lagrange1_Gauss6::
  calcGradN( localIndex const q,
             real64 const (&X)[numNodes][3],
             real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  jacobianTransformation( q, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyJacobianTransformationToShapeFunctionsDerivatives( q, J, gradN );

  return detJ * weight;
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Wedge_Lagrange1_Gauss6::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  jacobianTransformation( q, X, J );

  return LvArray::tensorOps::determinant< 3 >( J ) * weight;
}

}
}

#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1WEDGELAGRANGE1GAUSS6
