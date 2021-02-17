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

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_

/**
 * @file LagrangeBasis5GL.hpp
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
 *                 o-------------o------------------------o------------------------o------------------------o-------------------o  ---> xi
 *  Index:         0             1                        2                        3                        4                   5
 *  Coordinate:   -1    -sqrt(1/21(7+2*sqrt(7))) -sqrt(1/21(7-2*sqrt(7)))  sqrt(1/21(7-2*sqrt(7)))  sqrt(1/21(7-2*sqrt(7)))     1
 *
 */
class LagrangeBasis5GL
{
public:

  /**
   * @brief Calculate the parent coordinates for the xi0 direction, given the
   *   linear index of a support point. Here we hardcode GL points at order 5.
   * @param supportPointIndex The linear index of support point
   * @return parent coordinate in the xi0 direction.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 parentSupportCoord( const localIndex supportPointIndex )
  {

    real64 result=0.0;
   
    switch (supportPointIndex)
    {
      case 0:
        result = -1.0;
      break;

      case 1:
        result = -sqrt((1.0/21.0)*(7.0+2.0*sqrt(7.0)));
        break;

      case 2:
        result = -sqrt((1.0/21.0)*(7.0-2.0*sqrt(7.0)));
        break;

      case 3:
        result = sqrt((1.0/21.0)*(7.0-2.0*sqrt(7.0)));
        break;

      case 4:
        result = sqrt((1.0/21.0)*(7.0+2.0*sqrt(7.0)));
        break;

      case 5:
        result = 1.0;
        break;

      default:
        break;
    }

    return result;
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
  constexpr static real64 value( const int index,
                                 const real64 xi )
  {
     
    real64 result=0.0;

    switch (index)
    {
      case 0:
        return result = LagrangeBasis5GL::value0(xi);
        break;

      case 1:
        return result = LagrangeBasis5GL::value1(xi);
        break;

      case 2:
        return result = LagrangeBasis5GL::value2(xi);
        break;

      case 3:
        return result = LagrangeBasis5GL::value3(xi);
        break;

      case 4:
        return result = LagrangeBasis5GL::value4(xi);
        break;

      case 5:
        return result = LagrangeBasis5GL::value5(xi);
        break;

      default:
        break;
    }
    
    return result;

  }


  /**
   * @brief The value of the basis function for the 0 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value0( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 0. Here we need the points
       at index 3,4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (-21.0/16.0)*(xi*xi*xi*xi*xi-xi*xi*xi*xi-(lambda3*lambda3+lambda4*lambda4)*xi*xi*xi+(lambda3*lambda3+lambda4*lambda4)*xi*xi+
           lambda3*lambda3*lambda4*lambda4*xi-lambda3*lambda3*lambda4*lambda4);
  }


  /**
   * @brief The value of the basis function for the 1 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value1( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 1. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);

    return ((21.0/16.0)*sqrt(7.0-sqrt(7.0)/2.0))*(xi*xi*xi*xi*xi-lambda4*xi*xi*xi*xi-(lambda3*lambda3+1)*xi*xi*xi+lambda4*(lambda3*lambda3+1)*xi*xi+
           lambda3*lambda3*xi-lambda4*lambda3*lambda3);
  }

/**
   * @brief The value of the basis function for the 2 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value2( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 2. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return ((-21.0/16.0)*sqrt(7.0+sqrt(7.0)/2.0))*(xi*xi*xi*xi*xi-lambda3*xi*xi*xi*xi-(lambda4*lambda4+1)*xi*xi*xi+lambda3*(lambda4*lambda4+1)*xi*xi+
           lambda4*lambda4*xi-lambda3*lambda4*lambda4);
  }

/**
   * @brief The value of the basis function for the 3 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value3( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 3. Here we need the points
       at index 3 and 4 called lambda1, lambda2. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return ((21.0/16.0)*sqrt(7.0+sqrt(7.0)/2.0))*(xi*xi*xi*xi*xi+lambda3*xi*xi*xi*xi-(lambda4*lambda4+1)*xi*xi*xi-lambda3*(lambda4*lambda4+1)*xi*xi+
           lambda4*lambda4*xi+lambda3*lambda4*lambda4);
  }

/**
   * @brief The value of the basis function for the 4 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value4( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 4. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return ((-21.0/16.0)*sqrt(7.0-sqrt(7.0)/2.0))*(xi*xi*xi*xi*xi+lambda4*xi*xi*xi*xi-(lambda3*lambda3+1)*xi*xi*xi-lambda4*(lambda3*lambda3+1)*xi*xi+
           lambda3*lambda3*xi+lambda4*lambda3*lambda3);
  }

/**
   * @brief The value of the basis function for the 5 support point.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value5( const real64 xi )
  {
    /* Define the two GL points needed to compute the basis function at point index 5. Here we need the points
       at index 3 and 4 called lambda3, lambda4. */

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);

    return (21.0/16.0)*(xi*xi*xi*xi*xi+xi*xi*xi*xi-(lambda4*lambda4+lambda3*lambda3)*xi*xi*xi-(lambda3*lambda3+lambda4*lambda4)*xi*xi+
           lambda3*lambda3*lambda4*lambda4*xi+lambda4*lambda4*lambda3*lambda3);
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
  constexpr static real64 gradient( const int index,
                                    const real64 xi )
  {

    real64 result=0.0;

    switch (index)
    {
      case 0:
        result = LagrangeBasis5GL::gradient0(xi);
        break;

      case 1:
        result = LagrangeBasis5GL::gradient1(xi);
        break;

      case 2:
        result = LagrangeBasis5GL::gradient2(xi);
        break;

      case 3:
        result = LagrangeBasis5GL::gradient3(xi);
        break;

      case 4:
        result = LagrangeBasis5GL::gradient4(xi);
        break;

      case 5:
        result = LagrangeBasis5GL::gradient5(xi);
        break;

      default:
        break;
    }

    return result;

  }

  /**
   * @brief The gradient of the basis function for support point 0 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 0.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient0( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (-21.0/16.0)*(5.0*xi*xi*xi*xi-4.0*xi*xi*xi-3.0*(lambda3*lambda3+lambda4*lambda4)*xi*xi+2.0*(lambda3*lambda3+lambda4*lambda4)*xi+lambda3*lambda3*lambda4*lambda4);

  }

/**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 1.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient1( const real64 xi )
  {

    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);
    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);

    return (21.0/16.0)*sqrt(7.0-sqrt(7.0)/2.0)*(5.0*xi*xi*xi*xi-4.0*lambda4*xi*xi*xi-3.0*(lambda3*lambda3+1.0)*xi*xi+2.0*lambda4*(lambda3*lambda3+1.0)*xi+lambda3*lambda3);

  }

  /**
   * @brief The gradient of the basis function for support point 2 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 2.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient2( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (-21.0/16.0)*sqrt(7.0+sqrt(7.0)/2.0)*(5.0*xi*xi*xi*xi-4.0*lambda3*xi*xi*xi-3.0*(lambda4*lambda4+1.0)*xi*xi+2.0*lambda3*(lambda4*lambda4+1.0)*xi+lambda4*lambda4);

  }

  /**
   * @brief The gradient of the basis function for support point 3 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 3.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient3( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (21.0/16.0)*sqrt(7.0+sqrt(7.0)/2.0)*(5.0*xi*xi*xi*xi+4.0*lambda3*xi*xi*xi-3.0*(lambda4*lambda4+1.0)*xi*xi-2*lambda3*(lambda4*lambda4+1.0)*xi+lambda4*lambda4);

  }

  /**
   * @brief The gradient of the basis function for support point 4 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 0.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient4( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (-21.0/16.0)*sqrt(7.0-sqrt(7.0)/2.0)*(5.0*xi*xi*xi*xi+4.0*lambda4*xi*xi*xi-3.0*(lambda3*lambda3+1.0)*xi*xi-2.0*lambda4*(lambda3*lambda3+1.0)*xi+lambda3*lambda3);

  }

  /**
   * @brief The gradient of the basis function for support point 5 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function at point 5.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient5( const real64 xi )
  {

    real64 lambda4=LagrangeBasis5GL::parentSupportCoord(4);
    real64 lambda3=LagrangeBasis5GL::parentSupportCoord(3);

    return (21.0/16.0)*(5.0*xi*xi*xi*xi+4.0*xi*xi*xi-3.0*(lambda3*lambda3+lambda4*lambda4)*xi*xi-2.0*(lambda3*lambda3+lambda4*lambda4)*xi+lambda3*lambda3*lambda4*lambda4);

  }



  /**
   * @struct TensorProduct3D
   *
   * A 3-dimensional basis formed from the tensor product of the 1d basis.
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
   */

  struct TensorProduct3D
  {
    /// The number of support points in the basis.
    constexpr static localIndex numSupportPoints = 216;

    /**
     * @brief Calculates the linear index for support/quadrature points from ijk
     *   coordinates.
     * @param i The index in the xi0 direction (0,5)
     * @param j The index in the xi1 direction (0,5)
     * @param k The index in the xi2 direction (0,5)
     * @return The linear index of the support/quadrature point (0-125)
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 6 * j + 36 * k;
    }

    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param linearIndex The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     * @param i2 The Cartesian index of the support point in the xi2 direction.
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    constexpr static void multiIndex( const int linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {
       
       i2 = linearIndex/36;
  
       i1 = (linearIndex%36)/6;

       i0 = (linearIndex%36)%6;
    
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
    static void value( const real64 (& coords)[3],
                       real64 (& N)[numSupportPoints] )
    {
      for( int a=0; a<6; ++a )
      {
        for( int b=0; b<6; ++b )
        {
          for( int c=0; c<6; ++c )
          {
            const int lindex = LagrangeBasis5GL::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis5GL::value( a, coords[0] ) *
                          LagrangeBasis5GL::value( b, coords[1] ) *
                          LagrangeBasis5GL::value( c, coords[2] );
          }
        }
      }
    }

  };

};

}
}


#endif /* GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS5GL_HPP_ */
