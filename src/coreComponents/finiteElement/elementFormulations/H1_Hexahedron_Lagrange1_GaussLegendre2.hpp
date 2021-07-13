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
 * @file H1_Hexahedron_Lagrange1_GaussLegendre2.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_TRILINEARHEXAHEDRON
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_TRILINEARHEXAHEDRON

#include "FiniteElementBase.hpp"
#include "LagrangeBasis1.hpp"

#include <utility>



namespace geosx
{
namespace finiteElement
{

//constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 + 0.5 * 0.5773502691896257645092,
//                                                       0.5 - 0.5 * 0.5773502691896257645092 };
//__constant__ real64 psiProduct[3] = { 0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[0],
//                                          0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[1],
//                                          0.5 * linearBasisAtQuadrature[1]*linearBasisAtQuadrature[1] };

//__constant__ real64 psiProduct[3] = {  0.311004233964073108,
//                                       0.083333333333333333,
//                                       0.022329099369260226 };
//
//__constant__ short dpsi[2] = { -1, 1 };

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
  constexpr static localIndex numNodes = LagrangeBasis1::TensorProduct3D::numSupportPoints;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 8;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~H1_Hexahedron_Lagrange1_GaussLegendre2() override
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
  GEOSX_FORCE_INLINE
  static void calcN( localIndex const q,
                     real64 (& N)[numNodes] )
  {
    int qa, qb, qc;
    LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );
    real64 const qCoords[3] = { quadratureFactor * LagrangeBasis1::parentSupportCoord( qa ),
                                quadratureFactor * LagrangeBasis1::parentSupportCoord( qb ),
                                quadratureFactor * LagrangeBasis1::parentSupportCoord( qc ) };

    LagrangeBasis1::TensorProduct3D::value( qCoords, N );
  }

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
    int qa, qb, qc;
    LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );
    jacobianTransformation( qa, qb, qc, X, J );
    return LvArray::tensorOps::invert< 3 >( J );
  }


  /**
   * @brief Calculate the symmetric gradient of a vector valued support field
   *   at a quadrature point using the stored inverse of the Jacobian
   *   transformation matrix.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The symmetric gradient in Voigt notation.
   */
  GEOSX_HOST_DEVICE
  static void symmetricGradient( int const q,
                                 real64 const (&invJ)[3][3],
                                 real64 const (&var)[numNodes][3],
                                 real64 ( &grad )[6] );



  /**
   * @brief Calculate the gradient of a vector valued support field at a point
   *   using the stored basis function gradients for all support points.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The vector valued support field to apply the gradient
   *   operator on.
   * @param grad The gradient.
   *
   * More precisely, the operator is defined as:
   * \f[
   * grad_{ij}  = \sum_a^{nSupport} \left ( \frac{\partial N_a}{\partial X_j} var_{ai}\right ),
   * \f]
   *
   */
  GEOSX_HOST_DEVICE
  static void gradient( int const q,
                        real64 const (&invJ)[3][3],
                        real64 const (&var)[numNodes][3],
                        real64 ( &grad )[3][3] );


  /**
   * @brief Inner product of all basis function gradients and a rank-2
   *   symmetric tensor evaluated at a quadrature point.
   * @param q The linear index of the quadrature point.
   * @param invJ The inverse of the Jacobian transformation matrix.
   * @param var The rank-2 symmetric tensor at @p q.
   * @param R The vector resulting from the tensor contraction.
   *
   * More precisely, the operator is defined as:
   * \f[
   * R_i = \sum_a^{nSupport} \left( \frac{\partial N_a}{\partial X_j} var_{ij} \right),
   * \f]
   * where \f$\frac{\partial N_a}{\partial X_j}\f$ is the basis function gradient,
   *   \f$var_{ij}\f$ is the rank-2 symmetric tensor.
   */
  GEOSX_HOST_DEVICE
  static void plusGradNajAij( int const q,
                              real64 const (&invJ)[3][3],
                              real64 const (&var)[6],
                              real64 ( &R )[numNodes][3] );



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
  static void jacobianTransformation( int const qa,
                                      int const qb,
                                      int const qc,
                                      real64 const (&X)[numNodes][3],
                                      real64 ( &J )[3][3] );


  /**
   * @brief Apply a Jacobian transformation matrix from the parent space to the
   *   physical space on the parent shape function derivatives, producing the
   *   shape function derivatives in the physical space.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param invJ The Jacobian transformation from parent->physical space.
   * @param gradN Array to contain the shape function derivatives for all
   *   support points at the coordinates of the quadrature point @p q.
   */
  GEOSX_HOST_DEVICE
  static void
    applyTransformationToParentGradients( int const qa,
                                          int const qb,
                                          int const qc,
                                          real64 const ( &invJ )[3][3],
                                          real64 ( &gradN )[numNodes][3] );


private:
  /// The length of one dimension of the parent element.
  constexpr static real64 parentLength = LagrangeBasis1::parentSupportCoord( 1 ) - LagrangeBasis1::parentSupportCoord( 0 );

  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = parentLength*parentLength*parentLength;

  /// The weight of each quadrature point.
  constexpr static real64 weight = parentVolume / numQuadraturePoints;

  /// The scaling factor specifying the location of the quadrature points
  /// relative to the origin and the outer extent of the element in the
  /// parent space.
  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;

//  constexpr static real64 psiProduct0 = 0.311004233964073108;
//  constexpr static real64 psiProduct1 = 0.083333333333333333;
//  constexpr static real64 psiProduct2 = 0.022329099369260226;
//
//  constexpr static short dpsi0 = -1;
//  constexpr static short dpsi1 = 1;


  /**
   * @brief Applies a function inside a generic loop in over the tensor product
   *   indices.
   * @tparam FUNC The type of function to call within the support loop.
   * @tparam PARAMS The parameter pack types to pass through to @p FUNC.
   * @param qa The 1d quadrature point index in xi0 direction (0,1)
   * @param qb The 1d quadrature point index in xi1 direction (0,1)
   * @param qc The 1d quadrature point index in xi2 direction (0,1)
   * @param func The function to call within the support loop.
   * @param params The parameters to pass to @p func.
   */
  template< typename FUNC, typename ... PARAMS >
  GEOSX_HOST_DEVICE
  static void supportLoop( int const qa,
                           int const qb,
                           int const qc,
                           FUNC && func,
                           PARAMS &&... params );

};

//GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE real64
//psiProductFunc( int const n )
//{
//  // factor = 1./(2 + Sqrt[3]) = 0.267949192431122706
//
//  real64 rval = 0.311004233964073108;
//  for( int a=0; a<n; ++a )
//  {
//    rval *= 0.267949192431122706;
//  }
//  return rval;
//}

/// @cond Doxygen_Suppress

template< typename FUNC, typename ... PARAMS >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
H1_Hexahedron_Lagrange1_GaussLegendre2::supportLoop( int const qa,
                                                     int const qb,
                                                     int const qc,
                                                     FUNC && func,
                                                     PARAMS &&... params )
{

/// Options for how to calculate the parent gradients.
  #define PARENT_GRADIENT_METHOD 2
#if PARENT_GRADIENT_METHOD == 1
  // This option calculates the basis values at the quadrature point for each
  // linear basis index.

  real64 const quadratureCoords[3] = { -quadratureFactor + 1.154700538379252 * qa,
                                       -quadratureFactor + 1.154700538379252 * qb,
                                       -quadratureFactor + 1.154700538379252 * qc };

  real64 const psi0[2] = { 0.5 - 0.5 * quadratureCoords[0],
                           0.5 + 0.5 * quadratureCoords[0] };
  real64 const psi1[2] = { 0.5 - 0.5 * quadratureCoords[1],
                           0.5 + 0.5 * quadratureCoords[1] };
  real64 const psi2[2] = { 0.5 - 0.5 * quadratureCoords[2],
                           0.5 + 0.5 * quadratureCoords[2] };
  constexpr real64 dpsi[2] = { -0.5, 0.5 };
#elif PARENT_GRADIENT_METHOD == 2
  // This option calculates the product of linear basis prior to use.
  // The tensor product basis gradient may be expressed as a permutation of the
  // product between the two possible linear basis gradients. Thus the values
  // in the basis loop of qaa/qbb/qcc indicate which permutation to choose.
  // The quantities qaa, qbb, qcc are the difference in index between the
  // quadrature point and the basis. This is possible because there are 8
  // basis, and 8 quadrature points, which have correlated indices. So the
  // values of qaa/qbb/qcc are the "distance" from the quadrature point index
  // and the support point index.
  // THIS approach uses about 10 less registers than option 1, with no apparent
  // cost.

//  constexpr static real64 linearBasisAtQuadrature[2] = { 0.5 + 0.5 * quadratureFactor,
//                                                         0.5 - 0.5 * quadratureFactor };
//  constexpr static real64 psiProduct[3] = { 0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[0],
//                                            0.5 * linearBasisAtQuadrature[0]*linearBasisAtQuadrature[1],
//                                            0.5 * linearBasisAtQuadrature[1]*linearBasisAtQuadrature[1] };

  /// { 1/12 (2 + Sqrt[3]), 1/12, 1/12 (2 - Sqrt[3]) }
  constexpr static real64 psiProduct[3] = { 0.311004233964073108, 0.083333333333333333, 0.022329099369260226};
  constexpr static int dpsi[2] = { -1, 1 };

//  constexpr static real64 psiProduct[3] = { psiProduct0, psiProduct1, psiProduct2 };
//  constexpr short dpsi[2] = { dpsi0, dpsi1 };
#endif

  // Loop over the linear basis indices in each direction.
  for( int a=0; a<2; ++a )
  {
#if PARENT_GRADIENT_METHOD == 2
    int const qaa = ( a^qa ); // abs(a-qa)
#endif
    for( int b=0; b<2; ++b )
    {
#if PARENT_GRADIENT_METHOD == 2
      int const qbb = ( b^qb );
#endif
      for( int c=0; c<2; ++c )
      {
#if PARENT_GRADIENT_METHOD == 2
        int const qcc = ( c^qc );
#endif

#if PARENT_GRADIENT_METHOD == 1
        real64 const dNdXi[3] = { dpsi[a] * psi1[b] * psi2[c],
                                  psi0[a] * dpsi[b] * psi2[c],
                                  psi0[a] * psi1[b] * dpsi[c] };
#elif PARENT_GRADIENT_METHOD == 2


//        const real64 dNdXi[3] = { dpsi[a] * psiProductFunc( qbb + qcc ),
//                                  dpsi[b] * psiProductFunc( qaa + qcc ),
//                                  dpsi[c] * psiProductFunc( qaa + qbb ) };

        real64 const dNdXi[3] = { dpsi[a] * psiProduct[ qbb + qcc ],
                                  dpsi[b] * psiProduct[ qaa + qcc ],
                                  dpsi[c] * psiProduct[ qaa + qbb ] };
#endif
        localIndex const nodeIndex = LagrangeBasis1::TensorProduct3D::linearIndex( a, b, c );

        func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
      }
    }
  }
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange1_GaussLegendre2::calcGradN( localIndex const q,
                                                   real64 const (&X)[numNodes][3],
                                                   real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};


  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( qa, qb, qc, J, gradN );

  return detJ * weight;
}

//*************************************************************************************************
#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange1_GaussLegendre2::
  jacobianTransformation( int const qa,
                          int const qb,
                          int const qc,
                          real64 const (&X)[numNodes][3],
                          real64 ( & J )[3][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&X)[numNodes][3],
                                                  real64 (& J)[3][3] )
  {
    real64 const * const GEOSX_RESTRICT Xnode = X[nodeIndex];
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        J[i][j] = J[i][j] + dNdXi[ j ] * Xnode[i];
      }
    }

//    J[0][0] = J[0][0] + dNdXi[0] * Xnode[0];
//    J[0][1] = J[0][1] + dNdXi[1] * Xnode[0];
//    J[0][2] = J[0][2] + dNdXi[2] * Xnode[0];
//    J[1][0] = J[1][0] + dNdXi[0] * Xnode[1];
//    J[1][1] = J[1][1] + dNdXi[1] * Xnode[1];
//    J[1][2] = J[1][2] + dNdXi[2] * Xnode[1];
//    J[2][0] = J[2][0] + dNdXi[0] * Xnode[2];
//    J[2][1] = J[2][1] + dNdXi[1] * Xnode[2];
//    J[2][2] = J[2][2] + dNdXi[2] * Xnode[2];

  }, X, J );
}


//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
H1_Hexahedron_Lagrange1_GaussLegendre2::
  applyTransformationToParentGradients( int const qa,
                                        int const qb,
                                        int const qc,
                                        real64 const ( &invJ )[3][3],
                                        real64 (& gradN)[numNodes][3] )
{
  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 (& gradN)[numNodes][3] )
  {
//    for( int i = 0; i < 3; ++i )
//    {
//      gradN[nodeIndex][i] = 0.0;
//      for( int j = 0; j < 3; ++j )
//      {
//        gradN[nodeIndex][i] = gradN[nodeIndex][i] + dNdXi[ j ] * invJ[j][i];
//      }
//    }
    // smaller register footprint by manually unrolling the for loops.
    gradN[nodeIndex][0] = dNdXi[0] * invJ[0][0] + dNdXi[1] * invJ[1][0] + dNdXi[2] * invJ[2][0];
    gradN[nodeIndex][1] = dNdXi[0] * invJ[0][1] + dNdXi[1] * invJ[1][1] + dNdXi[2] * invJ[2][1];
    gradN[nodeIndex][2] = dNdXi[0] * invJ[0][2] + dNdXi[1] * invJ[1][2] + dNdXi[2] * invJ[2][2];


  }, invJ, gradN );
}


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
H1_Hexahedron_Lagrange1_GaussLegendre2::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  return LvArray::tensorOps::determinant< 3 >( J );
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange1_GaussLegendre2::symmetricGradient( int const q,
                                                                real64 const (&invJ)[3][3],
                                                                real64 const (&var)[numNodes][3],
                                                                real64 (& grad)[6] )
{
  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 const (&var)[numNodes][3],
                                                  real64 (& grad)[6] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }

    grad[0] = grad[0] + gradN[0] * var[ nodeIndex ][0];
    grad[1] = grad[1] + gradN[1] * var[ nodeIndex ][1];
    grad[2] = grad[2] + gradN[2] * var[ nodeIndex ][2];
    grad[3] = grad[3] + gradN[2] * var[ nodeIndex ][1] + gradN[1] * var[ nodeIndex ][2];
    grad[4] = grad[4] + gradN[2] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][2];
    grad[5] = grad[5] + gradN[1] * var[ nodeIndex ][0] + gradN[0] * var[ nodeIndex ][1];
  }, invJ, var, grad );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange1_GaussLegendre2::plusGradNajAij( int const q,
                                                             real64 const (&invJ)[3][3],
                                                             real64 const (&var)[6],
                                                             real64 (& R)[numNodes][3] )
{
  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc,
               [] GEOSX_HOST_DEVICE
                 ( real64 const (&dNdXi)[3],
                 int const nodeIndex,
                 real64 const (&invJ)[3][3],
                 real64 const (&var)[6],
                 real64 (& R)[numNodes][3] )
  {

    real64 gradN[3] = {0, 0, 0};
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        gradN[i] = gradN[i] + dNdXi[ j ] * invJ[j][i];
      }
    }
    R[ nodeIndex ][ 0 ] = R[ nodeIndex ][ 0 ] - var[ 0 ] * gradN[ 0 ] - var[ 5 ] * gradN[ 1 ] - var[ 4 ] * gradN[ 2 ];
    R[ nodeIndex ][ 1 ] = R[ nodeIndex ][ 1 ] - var[ 5 ] * gradN[ 0 ] - var[ 1 ] * gradN[ 1 ] - var[ 3 ] * gradN[ 2 ];
    R[ nodeIndex ][ 2 ] = R[ nodeIndex ][ 2 ] - var[ 4 ] * gradN[ 0 ] - var[ 3 ] * gradN[ 1 ] - var[ 2 ] * gradN[ 2 ];
  }, invJ, var, R );
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void H1_Hexahedron_Lagrange1_GaussLegendre2::gradient( int const q,
                                                       real64 const (&invJ)[3][3],
                                                       real64 const (&var)[numNodes][3],
                                                       real64 (& grad)[3][3] )
{
  int qa, qb, qc;
  LagrangeBasis1::TensorProduct3D::multiIndex( q, qa, qb, qc );

  supportLoop( qa, qb, qc, [] GEOSX_HOST_DEVICE ( real64 const (&dNdXi)[3],
                                                  int const nodeIndex,
                                                  real64 const (&invJ)[3][3],
                                                  real64 const (&var)[numNodes][3],
                                                  real64 (& grad)[3][3] )
  {
    for( int i = 0; i < 3; ++i )
    {
      real64 gradN=0.0;;
      for( int j = 0; j < 3; ++j )
      {
        gradN = gradN + dNdXi[ j ] * invJ[j][i];
      }
      for( int k = 0; k < 3; ++k )
      {
        grad[k][i] = grad[k][i] + gradN * var[ nodeIndex ][k];
      }
    }
  }, invJ, var, grad );
}

/// @endcond

#if __GNUC__
#pragma GCC diagnostic pop
#endif

#undef PARENT_GRADIENT_METHOD
}
}

#endif //FINITE_ELEMENT_SHAPE_KERNEL
