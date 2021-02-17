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
 * @file Q5_Hexahedron_Lagrange_GaussLobatto.hpp
 */

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_TRILINEARHEXAHEDRON
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_TRILINEARHEXAHEDRON

#include "FiniteElementBase.hpp"
#include "LagrangeBasis5GL.hpp"

#include <utility>



namespace geosx
{
namespace finiteElement
{


/**
 * This class contains the kernel accessible functions specific to the standard
 * Trilinear Hexahedron finite element with a Gauss Lobatto quadrature rule. It is
 * assumed that the indexing for the quadrature points mirrors that of the
 * nodes. Also note that the assumed node ordering is not the standard
 * right-hand-rule used in the literature. Here we use a Cartesian aligned
 * numbering in order to simplify the mapping to the parent coordinates and
 * tensor product indices.
 *
 *
 *                210      211      212      213      214      215                  _______________________________________________________________________
 *                   o--------o--------o--------o--------o--------o                |Node      xi0                        xi1                          xi2  |
 *                  /.                                           /|                |=====     ===                        ===                          ===  |
 *            204  / .  205      206      207      208      209 / |                |  0       -1                         -1                           -1   |
 *                o  .     o        o        o        o        o  |                |  1   -sqrt(1/21(7+/sqrt(7))         -1                           -1   |
 *               /   o                                        /   |                |  2    -sqrt(1/21(7-/sqrt(7))        -1                           -1   |
 *         198  /    .199     200      201      202      203 /    o                |  3    sqrt(1/21(7-/sqrt(7))         -1                           -1   |
 *             o     .  o        o        o        o        o     |                |  4    sqrt(1/21(7+/sqrt(7))         -1                           -1   |
 *            /      .                                     /      |                |  5       -1                         -1                           -1   |
 *      192  /   193 o     194      195      196      197 /    o  |                |  6       -1                 -sqrt(1/21(7+/sqrt(7))               -1   |
 *          o        o        o        o        o        o        o                |  7   -sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7))               -1   |
 *         /         .                                  /         |                |  8   -sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7))               -1   |
 *    186 /    187   .  188      189      190      191 /    o     |                |  9    sqrt(1/21(7-/sqrt(7)) -sqrt(1/21(7+/sqrt(7))               -1   |
 *       o        o  o     o        o        o        o        o  |                | 10    sqrt(1/21(7+/sqrt(7)) -sqrt(1/21(7+/sqrt(7))               -1   |
 *      /            .                               /            o                | 11       -1                 -sqrt(1/21(7+/sqrt(7))               -1   |
 * 180 /    181      . 182    183      184      185 /    o        |                | ..       ..                         ..                           ..   |
 *    o--------o--------o--------o--------o--------o        o     |                | ..       ..                         ..                           ..   |
 *    |           o  .                             |           o  |                | 204      -1                  sqrt(1/21(7+/sqrt(7))               1    |
 *    |  o           o        o        o        o  |  o  o        o                | 205  -sqrt(1/21(7+/sqrt(7))  sqrt(1/21(7+/sqrt(7))               1    |
 *    |     o        .                             |     o        |                | 206  -sqrt(1/21(7-/sqrt(7))  sqrt(1/21(7-/sqrt(7))               1    |
 *    |        o     .                             |        o     |                | 207  sqrt(1/21(7+/sqrt(7))   sqrt(1/21(7-/sqrt(7))               1    |
 *    o           o  .                             o           o  |                | 208  sqrt(1/21(7-/sqrt(7))   sqrt(1/21(7+/sqrt(7))               1    |
 *    |  o           .                             |  o           |                | 209       1                  sqrt(1/21(7+/sqrt(                  1    |
 *    |     o        o--------o--------o--------o--|-----o--------o                | 210      -1                          1                           1    |
 *    |        o    ,30       31      32        33 |     34 o    /35               | 211  -sqrt(1/21(7+/sqrt(7))          1                           1    |
 *    o            ,                               o            /                  | 212  -sqrt(1/21(7-/sqrt(7))          1                           1    |
 *    |  o        o        o         o       o     |  o        o                   | 213   sqrt(1/21(7-/sqrt(7))          1                           1    |
 *    |     o    ,24       25        26      27    |  28 o    /29                  | 214   sqrt(1/21(7+/sqrt(7))          1                           1    |
 *    |         ,                                  |         /                     | 215       1                          1                           1    |
 *    o        o        o         o       o     22 o        o                      |_______________________________________________________________________|
 *    |  o    ,18       19        20      21       |  o    /23
 *    |      ,                                     |      /
 *    |     o        o         o       o        o  |     o
 *    o    ,12       13        14      15       16 o    /17
 *    |   ,                                        |   /
 *    |  o        o        o        o        o     |  o               xi2
 *    | ,6        7        8        9        10    | /11               |
 *    |,                                           |/                  | / xi1
 *    o--------o--------o--------o--------o--------o                   |/
 *    0        1        2        3        4        5                   o----- xi0
 *
 *
 */
class Q5_Hexahedron_Lagrange_GaussLobatto final : public FiniteElementBase
{
public:
  /// The number of nodes/support points per element.
  constexpr static localIndex numNodes = LagrangeBasis5GL::TensorProduct3D::numSupportPoints;

  /// The number of quadrature points per element.
  constexpr static localIndex numQuadraturePoints = 216;

  /** @cond Doxygen_Suppress */
  USING_FINITEELEMENTBASE
  /** @endcond Doxygen_Suppress */

  virtual ~Q5_Hexahedron_Lagrange_GaussLobatto() override
  {}

  virtual localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

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
    LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

    real64 const qCoords[3] = { LagrangeBasis5GL::parentSupportCoord( qa ),
                                LagrangeBasis5GL::parentSupportCoord( qb ),
                                LagrangeBasis5GL::parentSupportCoord( qc ) };

    LagrangeBasis5GL::TensorProduct3D::value( qCoords, N );
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
    LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );
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
  static void plus_gradNajAij( int const q,
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
  constexpr static real64 parentLength = LagrangeBasis5GL::parentSupportCoord( 1 ) - LagrangeBasis5GL::parentSupportCoord( 0 );

  /// The volume of the element in the parent configuration.
  constexpr static real64 parentVolume = parentLength*parentLength*parentLength;

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

/// @cond Doxygen_Suppress
//MODFI2 : Changed the calcul of dNdXi implementing a new function gradient which calculate
// the derivative of the bases functions at the desired coord
template< typename FUNC, typename ... PARAMS >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
Q5_Hexahedron_Lagrange_GaussLobatto::supportLoop( int const qa,
                                                     int const qb,
                                                     int const qc,
                                                     FUNC && func,
                                                     PARAMS &&... params )
{
   
  real64 const qCoords[3] = { LagrangeBasis5GL::parentSupportCoord( qa ),
                              LagrangeBasis5GL::parentSupportCoord( qb ),
                              LagrangeBasis5GL::parentSupportCoord( qc ) };

  for( int c=0; c<6; ++c )
  {
    for( int b=0; b<6; ++b )
    {
      for( int a=0; a<6; ++a )
      {
        real64 const dNdXi[3] = { LagrangeBasis5GL::gradient( a, qCoords[0])*
                                  LagrangeBasis5GL::value( b, qCoords[1])*
                                  LagrangeBasis5GL::value( c, qCoords[2]),
                                  LagrangeBasis5GL::value( a, qCoords[0])*
                                  LagrangeBasis5GL::gradient( b, qCoords[1])*
                                  LagrangeBasis5GL::value( c, qCoords[2]),
                                  LagrangeBasis5GL::value( a, qCoords[0])*
                                  LagrangeBasis5GL::value( b, qCoords[1])*
                                  LagrangeBasis5GL::gradient( c, qCoords[2])};
        
        localIndex const nodeIndex = LagrangeBasis5GL::TensorProduct3D::linearIndex( a, b, c );

        func( dNdXi, nodeIndex, std::forward< PARAMS >( params )... );
      }
    }
  }

 
}

//*************************************************************************************************
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64
Q5_Hexahedron_Lagrange_GaussLobatto::calcGradN( localIndex const q,
                                                   real64 const (&X)[numNodes][3],
                                                   real64 (& gradN)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  //Define the Gauss Lobatto weights 
  real64 weight[6]={ 1.0/15.0,(1.0/30.0)*(14.0-sqrt(7.0)),(1.0/30.0)*(14.0+sqrt(7.0)),(1.0/30.0)*(14.0+sqrt(7.0)),(1.0/30.0)*(14.0-sqrt(7.0)),1.0/15.0};

  int qa, qb, qc;
  LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );
  
  real64 const detJ = LvArray::tensorOps::invert< 3 >( J );

  applyTransformationToParentGradients( qa, qb, qc, J, gradN );

//MODIF1 : Change the calcul of detJ multiplying by the right weight in each direction xi, yi, zi
  return detJ * weight[qa] * weight[qb] * weight[qc];
}

//*************************************************************************************************
#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
Q5_Hexahedron_Lagrange_GaussLobatto::
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
Q5_Hexahedron_Lagrange_GaussLobatto::
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
Q5_Hexahedron_Lagrange_GaussLobatto::
  transformedQuadratureWeight( localIndex const q,
                               real64 const (&X)[numNodes][3] )
{
  real64 J[3][3] = {{0}};

  int qa, qb, qc;
  LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

  jacobianTransformation( qa, qb, qc, X, J );

  return LvArray::tensorOps::determinant< 3 >( J );
}



GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void Q5_Hexahedron_Lagrange_GaussLobatto::symmetricGradient( int const q,
                                                                real64 const (&invJ)[3][3],
                                                                real64 const (&var)[numNodes][3],
                                                                real64 (& grad)[6] )
{
  int qa, qb, qc;
  LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

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
void Q5_Hexahedron_Lagrange_GaussLobatto::plus_gradNajAij( int const q,
                                                              real64 const (&invJ)[3][3],
                                                              real64 const (&var)[6],
                                                              real64 (& R)[numNodes][3] )
{
  int qa, qb, qc;
  LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

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
void Q5_Hexahedron_Lagrange_GaussLobatto::gradient( int const q,
                                                       real64 const (&invJ)[3][3],
                                                       real64 const (&var)[numNodes][3],
                                                       real64 (& grad)[3][3] )
{
  int qa, qb, qc;
  LagrangeBasis5GL::TensorProduct3D::multiIndex( q, qa, qb, qc );

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
