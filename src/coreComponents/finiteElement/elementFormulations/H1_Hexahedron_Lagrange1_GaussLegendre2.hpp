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
 * @file H1_Hexahedron_Lagrange1_GaussLegendre2.hpp
 */

#ifndef GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON
#define GEOSX_CORE_FINITEELEMENT_TRILINEARHEXAHEDRON

#include "FiniteElementBase.hpp"
#include "LagrangeBasis1.hpp"

namespace geosx
{
namespace finiteElement
{
/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gaussian quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
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
 *
 */

class H1_Hexahedron_Lagrange1_GaussLegendre2 final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes =
    LagrangeBasis1::TensorProduct3D::numSupportPoints;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 8;

  virtual ~H1_Hexahedron_Lagrange1_GaussLegendre2() override
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
  GEOSX_FORCE_INLINE
  static void
  shapeFunctionValues( localIndex const q, real64 ( &N )[numNodes] )
  {
    int qa, qb, qc;
    LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );
    real64 const qCoords[3] = {
      quadratureFactor * LagrangeBasis1::parentSupportCoord( qa ),
      quadratureFactor * LagrangeBasis1::parentSupportCoord( qb ),
      quadratureFactor * LagrangeBasis1::parentSupportCoord( qc ) };

    LagrangeBasis1::TensorProduct3D::value( qCoords, N );
  }

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
  /// The length of one dimension of the parent element.
  constexpr static real64 parentLength = LagrangeBasis1::parentSupportCoord( 1 ) -
    LagrangeBasis1::parentSupportCoord( 0 );

  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume =
    parentLength * parentLength * parentLength;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space.
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

  /**
   * @brief Calculates the isoparametric "Jacobian" transformation
   *   matrix/mapping from the parent space to the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param X Array containing the coordinates of the support points.
   * @param J Array to store the Jacobian transformation.
   */
  GEOSX_HOST_DEVICE
  static void
  jacobianTransformation( int const qa,
                          int const qb,
                          int const qc,
                          real64 const ( &X )[numNodes][3],
                          real64 ( &J )[3][3] );

  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param dNdX Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void
  applyJacobianTransformationToShapeFunctionsDerivatives(
    int const qa,
    int const qb,
    int const qc,
    real64 const ( &invJ )[3][3],
    real64 ( &dNdX )[numNodes][3] );
};

#if 1
//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange1_GaussLegendre2::shapeFunctionDerivatives(
  localIndex const q,
  real64 const ( &X )[numNodes][3],
  real64 ( &dNdX )[numNodes][3] )
{
  real64 J[3][3] = { { 0 } };

  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  real64 const detJ = inverse( J );

  applyJacobianTransformationToShapeFunctionsDerivatives( qa, qb, qc, J, dNdX );

  return detJ * weight;
}

#else

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange1_GaussLegendre2::shapeFunctionDerivatives(
  localIndex const q,
  real64 const ( &X )[numNodes][3],
  real64 ( &dNdX )[numNodes][3] )
{
  real64 J[3][3] = { { 0 } };

  real64 const quadratureCoords[3] = {
    quadratureFactor * LagrangeBasis1::TensorProduct3D::parentCoords0( q ),
    quadratureFactor * LagrangeBasis1::TensorProduct3D::parentCoords1( q ),
    quadratureFactor * LagrangeBasis1::TensorProduct3D::parentCoords2( q ) };

  real64 const psi0[2] = { 0.5 - 0.5 * quadratureCoords[0],
                           0.5 + 0.5 * quadratureCoords[0] };
  real64 const psi1[2] = { 0.5 - 0.5 * quadratureCoords[1],
                           0.5 + 0.5 * quadratureCoords[1] };
  real64 const psi2[2] = { 0.5 - 0.5 * quadratureCoords[2],
                           0.5 + 0.5 * quadratureCoords[2] };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };

  for( localIndex a = 0; a < 2; ++a )
  {
    for( localIndex b = 0; b < 2; ++b )
    {
      for( localIndex c = 0; c < 2; ++c )
      {
        real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                  psi0[a] * dpsi[b] * psi2[c],
                                  psi0[a] * psi1[b] * dpsi[c] };
        localIndex const nodeIndex =
          LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + dNdXi[j] * X[nodeIndex][i];
          }
        }
      }
    }
  }

  real64 const detJ = inverse( J );

  for( localIndex a = 0; a < 2; ++a )
  {
    for( localIndex b = 0; b < 2; ++b )
    {
      for( localIndex c = 0; c < 2; ++c )
      {
        real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                  psi0[a] * dpsi[b] * psi2[c],
                                  psi0[a] * psi1[b] * dpsi[c] };
        localIndex const nodeIndex =
          LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );
        for( int i = 0; i < 3; ++i )
        {
          dNdX[nodeIndex][i] = 0.0;
          for( int j = 0; j < 3; ++j )
          {
            dNdX[nodeIndex][i] = dNdX[nodeIndex][i] + dNdXi[j] * J[j][i];
          }
        }
      }
    }
  }

  return detJ;
}

#endif

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange1_GaussLegendre2::jacobianTransformation(
  int const qa,
  int const qb,
  int const qc,
  real64 const ( &X )[numNodes][3],
  real64 ( &J )[3][3] )
{
  constexpr static real64 linearBasisAtQuadrature[2] = {
    0.5 + 0.5 * quadratureFactor,
    0.5 - 0.5 * quadratureFactor };
  constexpr static real64 psiProduct[3] = {
    0.5 * linearBasisAtQuadrature[0] * linearBasisAtQuadrature[0],
    0.5 * linearBasisAtQuadrature[0] * linearBasisAtQuadrature[1],
    0.5 * linearBasisAtQuadrature[1] * linearBasisAtQuadrature[1] };
  constexpr static int dpsi[2] = { -1, 1 };

  for( int a = 0; a < 2; ++a )
  {
    int const qaa = a ^ qa;
    for( int b = 0; b < 2; ++b )
    {
      int const qbb = b ^ qb;
      for( int c = 0; c < 2; ++c )
      {
        int const qcc = c ^ qc;
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[qbb + qcc],
                                  dpsi[b] * psiProduct[qaa + qcc],
                                  dpsi[c] * psiProduct[qaa + qbb] };

        localIndex const nodeIndex =
          LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + dNdXi[j] * X[nodeIndex][i];
          }
        }
      }
    }
  }
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange1_GaussLegendre2::applyJacobianTransformationToShapeFunctionsDerivatives(
  int const qa,
  int const qb,
  int const qc,
  real64 const ( &invJ )[3][3],
  real64 ( &dNdX )[numNodes][3] )
{
  constexpr static real64 linearBasisAtQuadrature[2] = {
    0.5 + 0.5 * quadratureFactor,
    0.5 - 0.5 * quadratureFactor };
  constexpr static real64 psiProduct[3] = {
    0.5 * linearBasisAtQuadrature[0] * linearBasisAtQuadrature[0],
    0.5 * linearBasisAtQuadrature[0] * linearBasisAtQuadrature[1],
    0.5 * linearBasisAtQuadrature[1] * linearBasisAtQuadrature[1] };
  constexpr static int dpsi[2] = { -1, 1 };

  for( int a = 0; a < 2; ++a )
  {
    int const qaa = ( a ^ qa );
    for( int b = 0; b < 2; ++b )
    {
      int const qbb = ( b ^ qb );
      for( int c = 0; c < 2; ++c )
      {
        int const qcc = ( c ^ qc );
        real64 const dNdXi[3] = { dpsi[a] * psiProduct[qbb + qcc],
                                  dpsi[b] * psiProduct[qaa + qcc],
                                  dpsi[c] * psiProduct[qaa + qbb] };
        localIndex const nodeIndex =
          LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );
        for( int i = 0; i < 3; ++i )
        {
          dNdX[nodeIndex][i] = 0.0;
          for( int j = 0; j < 3; ++j )
          {
            dNdX[nodeIndex][i] = dNdX[nodeIndex][i] + dNdXi[j] * invJ[j][i];
          }
        }
      }
    }
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange1_GaussLegendre2::transformedQuadratureWeight(
  localIndex const q,
  real64 const ( &X )[numNodes][3] )
{
  real64 J[3][3] = { { 0 } };

  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  return detJ( J );
}

}  // namespace finiteElement
}  // namespace geosx

#endif  //FINITE_ELEMENT_SHAPE_KERNEL
