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
 * @file H1_Pyramid_Lagrange1_Gauss5.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1PYRAMIDLAGRANGE1GAUSS5
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1PYRAMIDLAGRANGE1GAUSS5

#include "FiniteElementBase.hpp"


namespace geosx
{
namespace finiteElement
{

/**
 * This class contains the kernel accessible functions specific to the
 * H1-conforming nodal pyramid finite element with a 5-point Gaussian
 * quadrature rule [1]. It is assumed that the indexing for the quadrature
 * points mirrors that of the nodes. Also note that the assumed node ordering
 * is not the standard right-hand-rule used in the literature.
 *
 *               4                                  =====  ===  ===  ===
 *              /\.                                 Node   xi0  xi1  xi2
 *             /' \`.                               =====  ===  ===  ===
 *            /'   \ `.                             0      -1   -1   -1
 *           /'     \  `.                           1       1   -1   -1
 *          /'       \   `.                         2      -1    1   -1
 *         /' 2       \    `.   3                   3       1    1   -1
 *        /o.......... \.....`o      xi2            4       0    0    1
 *       /,             \    /       |              =====  ===  ===  ===
 *      /,               \  /        | / xi1
 *     o------------------o          |/
 *     0                   1          ----- xi0
 *
 * References
 * [1] Felippa CA (2004). A compendium of FEM integration formulas for symbolic work.
 *     Engineering Computations 21(8), 867--890.
 *
 *
 */
class H1_Pyramid_Lagrange1_Gauss5 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = 5;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 5;

  virtual ~H1_Pyramid_Lagrange1_Gauss5() override
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
  /// The weight for quadrature points paired with base nodes
  constexpr static real64 weight = 81.0 / 100.0;

  /// The weight increment for the quadrature point paired with the apex node
  constexpr static real64 weightDelta  = 125.0 / 27.0 - weight;

  /// The factor specifying the (xi0,xi1) location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space for quadrature points paired with base nodes
  constexpr static real64 quadratureCrossSectionCoord = 0.584237394672177;

  /// The xi2-coordinate of the quadrature points relative to the origin
  /// in the parent space for quadrature points paired with base nodes
  constexpr static real64 quadratureLongitudinalCoordNeg = -2.0 / 3.0;

  /// The xi2-coordinate increment for the quadrature point paired with
  /// the apex node
  constexpr static real64 quadratureLongitudinalCoordDelta = 16.0 / 15.0;

  /**
   * @brief Calculates the linear index for support/quadrature points paired
   *        with base nodes from ij coordinates.
   * @param i The index in the xi0 direction (0,1)
   * @param j The index in the xi1 direction (0,1)
   * @return The linear index of the support/quadrature point (0-3)
   */
  template< typename T >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static T linearMap( T const i, T const j )
  {
    return i + 2 * j;
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords0( localIndex const a )
  {
    return -1.0 + 2.0 * ( a & 1 ) + 0.25 * ( a & 4 );
  }

  /**
   * @brief Calculate the parent coordinates for the xi1 direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi1 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords1( localIndex const a )
  {
    return -1.0 + ( a & 2 ) + 0.25 * ( a & 4 );
  }

  /**
   * @brief Calculate the parent coordinates for the xi2 direction, given the
   *        linear index of a support point.
   * @param a The linear index of support point
   * @return parent coordinate in the xi2 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentCoords2( localIndex const a )
  {
    return -1.0 + 0.5 * ( a & 4 );
  }

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *        linear index of a quadrature point.
   * @param a The linear index of quadrature point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords0( localIndex const q )
  {
    return parentCoords0( q ) * quadratureCrossSectionCoord;
  }

  /**
   * @brief Calculate the parent coordinates for the xi1 direction, given the
   *        linear index of a quadrature point.
   * @param a The linear index of quadrature point
   * @return parent coordinate in the xi1 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords1( localIndex const q )
  {
    return parentCoords1( q ) * quadratureCrossSectionCoord;
  }

  /**
   * @brief Calculate the parent coordinates for the xi2 direction, given the
   *        linear index of a quadrature point.
   * @param a The linear index of quadrature point
   * @return parent coordinate in the xi2 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureParentCoords2( localIndex const q )
  {
    return quadratureLongitudinalCoordNeg + 0.5 * ( 1 + parentCoords2( q ) ) * quadratureLongitudinalCoordDelta;
  }

  /**
   * @brief Return the integration weights for a quadrature point in the parent
   *        element.
   * @return The quadrature rule weight.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 quadratureWeight( localIndex const q )
  {
    return weight + 0.5 * ( 1 + parentCoords2( q ) ) * weightDelta;
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
H1_Pyramid_Lagrange1_Gauss5::
  jacobianTransformation( int const q,
                          real64 const (&X)[numNodes][3],
                          real64 ( & J )[3][3] )
{
  real64 const quadratureCoords[3] = { quadratureParentCoords0( q ),
                                       quadratureParentCoords1( q ),
                                       quadratureParentCoords2( q ) };

  real64 const psi0[2] = { 0.5 - 0.5*quadratureCoords[0],
                           0.5 + 0.5*quadratureCoords[0] };
  real64 const psi1[2] = { 0.5 - 0.5*quadratureCoords[1],
                           0.5 + 0.5*quadratureCoords[1] };
  real64 const psi2 = 0.5 - 0.5*quadratureCoords[2];
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  // Contributions from basis functions paired with base nodes
  for( localIndex a=0; a<2; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2,
                                psi0[a] * dpsi[b] * psi2,
                                psi0[a] * psi1[b] * dpsi[0] };
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

  // Contribution from the basis function paired with the apex nodes
  for( int i = 0; i < 3; ++i )
  {
    J[i][2] = J[i][2] + dpsi[1] * X[4][i];
  }
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Pyramid_Lagrange1_Gauss5::
  applyJacobianTransformationToShapeFunctionsDerivatives( int const q,
                                                          real64 const ( &invJ )[3][3],
                                                          real64 (& gradN)[numNodes][3] )
{
  real64 const quadratureCoords[3] = { quadratureParentCoords0( q ),
                                       quadratureParentCoords1( q ),
                                       quadratureParentCoords2( q ) };

  real64 const psi0[2] = { 0.5*( 1.0 - quadratureCoords[0] ),
                           0.5*( 1.0 + quadratureCoords[0] ) };
  real64 const psi1[2] = { 0.5*( 1.0 - quadratureCoords[1] ),
                           0.5*( 1.0 + quadratureCoords[1] ) };
  real64 const psi2 = 0.5*( 1.0 - quadratureCoords[2]);
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  // Contributions from basis functions paired with base nodes
  for( localIndex a=0; a<2; ++a )
  {
    for( localIndex b=0; b<2; ++b )
    {
      real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2,
                                psi0[a] * dpsi[b] * psi2,
                                psi0[a] * psi1[b] * dpsi[0] };
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

  // Contribution from the basis function paired with the apex nodes
  for( int i = 0; i < 3; ++i )
  {
    gradN[4][i] = dpsi[1] * invJ[2][i];
  }
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Pyramid_Lagrange1_Gauss5::
  calcN( localIndex const q,
         real64 ( & N )[numNodes] )
{
  real64 const xi[3] = { quadratureParentCoords0( q ),
                         quadratureParentCoords1( q ),
                         quadratureParentCoords2( q ) };

  N[0] = 0.125*( 1.0 - xi[0] ) * ( 1.0 - xi[1] ) * ( 1.0 - xi[2] );
  N[1] = 0.125*( 1.0 + xi[0] ) * ( 1.0 - xi[1] ) * ( 1.0 - xi[2] );
  N[2] = 0.125*( 1.0 - xi[0] ) * ( 1.0 + xi[1] ) * ( 1.0 - xi[2] );
  N[3] = 0.125*( 1.0 + xi[0] ) * ( 1.0 + xi[1] ) * ( 1.0 - xi[2] );
  N[4] = 0.5*( 1.0 + xi[2] );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 H1_Pyramid_Lagrange1_Gauss5::calcGradN( localIndex const q,
                                               real64 const (&X)[numNodes][3],
                                               real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  jacobianTransformation( q, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyJacobianTransformationToShapeFunctionsDerivatives( q, J, gradN );

  return detJ * quadratureWeight( q );
}

//*************************************************************************************************

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Pyramid_Lagrange1_Gauss5::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  jacobianTransformation( q, X, J );

  return LvArray::tensorOps::determinant< 3 >( J ) * quadratureWeight( q );
}


}
}
#endif //GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_H1PYRAMIDLAGRANGE1GAUSS5
