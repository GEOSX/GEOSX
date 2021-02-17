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

#ifndef GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_
#define GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_

/**
 * @file LagrangeBasis3GL.hpp
 */

#include "common/DataTypes.hpp"

namespace geosx
{
namespace finiteElement
{

/**
 * This class contains the implementation for a second order (quadratic) Lagrange
 * polynomial basis. The parent space is defined by:
 *
 *                 o---------o--------o---------o  ---> xi
 *  Index:         0         1        2         3
 *  Coordinate:   -1    -1/sqrt(5) 1/sqrt(5)    1
 *
 */
class LagrangeBasis3GL
{
public:
    /// The number of support points for the basis
    constexpr static localIndex numSupportPoints = 4;

    /**
     * @brief Calculate the parent coordinates for the xi0 direction, given the
     *   linear index of a support point.
     * @param supportPointIndex The linear index of support point
     * @return parent coordinate in the xi0 direction.
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    // MODIF1 : Harcoding the Gauss-Lobatto coordinates and return the right one
    // depending on the supportPointIndex value
    //Switch case
    constexpr static real64 parentSupportCoord( const localIndex supportPointIndex )
    {
        real64 result=0.0;

        switch (supportPointIndex) {
        case 0:
            result = -1.0;
            break;
        case 1:
            result = -1.0/sqrt(5.0);
            break;
        case 2:
            result = 1.0/sqrt(5.0);
            break;
        case 3:
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
//MODIF3 : Change the  way to return the base function evaluated at the desired coord
    constexpr static real64 value( const int index,
                                   const real64 xi )
    {
      real64 result=0.0;

      switch (index) {
      case 0:
          result = LagrangeBasis3GL::value0( xi );
          break;
      case 1:
          result = LagrangeBasis3GL::value1( xi );
          break;
      case 2:
          result = LagrangeBasis3GL::value2( xi );
          break;
      case 3:
          result = LagrangeBasis3GL::value3( xi );
          break;
      default:
          break;
      }

      return result;
    }


  /**
   * @brief The value of the basis function for support point 0.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  //MODFI4 : Implemented new base functions and their derivative for Q3
  constexpr static real64 value0( const real64 xi )
  {
    return -(5.0/8.0)*(xi*xi*xi-xi*xi-(1.0/5.0)*xi+1.0/5.0);
  }

  /**
   * @brief The value of the basis function for support point 1.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value1( const real64 xi )
  {
    return (5.0*sqrt(5.0)/8.0)*(xi*xi*xi-(1.0/sqrt(5.0))*xi*xi-xi+1.0/sqrt(5.0));
  }

  /**
   * @brief The value of the basis function for support point 2.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value2( const real64 xi )
  {
    return -(5.0*sqrt(5.0)/8.0)*(xi*xi*xi+(1.0/sqrt(5.0))*xi*xi-xi-1.0/sqrt(5.0));
  }

  /**
   * @brief The value of the basis function for support point 3.
   * @param xi The coordinate at which to evaluate the basis.
   * @return The value of the basis.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 value3( const real64 xi )
  {
    return (5.0/8.0)*(xi*xi*xi+xi*xi-(1.0/5.0)*xi-1.0/5.0);
  }


  /**
    * @brief The gradient of the basis function for a support point evaluated at a
    *   point along the axes.
    * @param index The index of the support point.
    * @param xi The coordinate at which to evaluate the basis.
    * @return The value of basis function.
    */
   GEOSX_HOST_DEVICE
   GEOSX_FORCE_INLINE
   //MODIF5 : New function returning the derivated base function at desired coord
   constexpr static real64 gradient( const int index,
                                  const real64 xi )
   {
       real64 result=0.0;

       switch (index) {
       case 0:
           result = LagrangeBasis3GL::gradient0( xi );
           break;
       case 1:
           result = LagrangeBasis3GL::gradient1( xi );
           break;
       case 2:
           result = LagrangeBasis3GL::gradient2( xi );
           break;
       case 3:
           result = LagrangeBasis3GL::gradient3( xi );
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
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient0( const real64 xi )
  {
    return -(5.0/8.0)*(3.0*xi*xi-2.0*xi-(1.0/5.0));
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient1( const real64 xi )
  {
    return (5.0*sqrt(5.0)/8.0)*(3.0*xi*xi-(2.0/sqrt(5.0))*xi-1.0);
  }

  /**
   * @brief The gradient of the basis function for support point 1 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient2( const real64 xi )
  {
    return -(5.0*sqrt(5.0)/8.0)*(3.0*xi*xi+(2.0/sqrt(5.0))*xi-1.0);
  }

  /**
   * @brief The gradient of the basis function for support point 3 evaluated at
   *   a point along the axes.
   * @param xi The coordinate at which to evaluate the gradient.
   * @return The gradient of basis function
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  constexpr static real64 gradient3( const real64 xi )
  {
    return (5.0/8.0)*(3.0*xi*xi+2.0*xi-(1.0/5.0));;
  }

  /**
   * @class TensorProduct3D
 *                                                                  _____________________________________
 *                                                                 |Node      xi0         xi1         xi2|
 *                                                                 |=====     ===         ===         ===|
 *                                                                 |  0       -1          -1          -1 |
 *                                                                 |  1   -1/sqrt(5)      -1          -1 |
 *                                                                 |  2    1/sqrt(5)      -1          -1 |
 *              60       61         62        63                   |  3        1          -1          -1 |
 *                o---------o---------o---------o                  |  4       -1      -1/sqrt(5)      -1 |
 *            56 /.     57        58        59 /|                  |  5   -1/sqrt(5)  -1/sqrt(5)      -1 |
 *              o .       o         o         o |                  |  6    1/sqrt(5)  -1/sqrt(5)      -1 |
 *          52 /  .   53        54        55 /  |                  |  7        1      -1/sqrt(5)      -1 |
 *            o   .     o         o         o   |                  |  8       -1       1/sqrt(5)      -1 |
 *        48 /    o 49      o 50      o 51 /    o                  |  9   -1/sqrt(5)   1/sqrt(5)      -1 |
 *          o---------o---------o---------o     |                  | 10    1/sqrt(5)   1/sqrt(5)      -1 |
 *          |   o .       o         o     |   o |                  | 11        1       1/sqrt(5)      -1 |
 *          |     .                       |     |                  | 12       -1           1          -1 |
 *          | o   o     o   o     o   o   | o   o                  | 13   -1/sqrt(5)       1          -1 |
 *          |     .                       |     |                  | 14    1/sqrt(5)       1          -1 |
 *          o   o .   o   o     o   o     o   o |                  | 15        1           1          -1 |
 *          |     .                       |     |                  | ..       ..          ..          .. |
 *          | o   .     o         o       | o   |                  | ..       ..          ..          .. |
 *          |     o.........o.........o...|.....o                  | 55        1      -1/sqrt(5)       1 |
 *          o    ,12  o     13  o     14  o    /15                 | 56       -1       1/sqrt(5)       1 |
 *          |   o         o         o     |   o                    | 57   -1/sqrt(5)   1/sqrt(5)       1 |
 *          |  ,8         9         10    |  /11       xi2         | 58    1/sqrt(5)   1/sqrt(5)       1 |
 *          | o         o         o       | o          |           | 59        1       1/sqrt(5)       1 |
 *          |,4         5         6       |/7          | / xi1     | 60       -1           1           1 |
 *          o---------o---------o---------o            |/          | 61   -1/sqrt(5)       1           1 |
 *         0         1         2         3             o----- xi0  | 62    1/sqrt(5)       1           1 |
 *                                                                 | 63        1           1           1 |
 *                                                                 |_____________________________________|
 *
 */
  struct TensorProduct3D
  {
      constexpr static localIndex numSupportPoints = 64;

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
    //MODIF6 : Change linearIndex for 64 nodes
    constexpr static int linearIndex( const int i,
                                      const int j,
                                      const int k )
    {
      return i + 4 * j + 16 * k;
    }



    /**
     * @brief Calculate the Cartesian/TensorProduct index given the linear index
     *   of a support point.
     * @param linearIndex The linear index of support point
     * @param i0 The Cartesian index of the support point in the xi0 direction.
     * @param i1 The Cartesian index of the support point in the xi1 direction.
     * @param i2 The Cartesian index of the support point in the xi2 direction.
     * @return
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE
    //MODIF7 : Change calcul of multiIndex
    constexpr static void multiIndex( int const linearIndex,
                                      int & i0,
                                      int & i1,
                                      int & i2 )
    {

      i2 = linearIndex/16;

      i1 = (linearIndex%16)/4;

      i0 = (linearIndex%16)%4;

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
      for( int a=0; a<4; ++a )
      {
        for( int b=0; b<4; ++b )
        {
          for( int c=0; c<4; ++c )
          {
            const int lindex = LagrangeBasis3GL::TensorProduct3D::linearIndex( a, b, c );
            N[ lindex ] = LagrangeBasis3GL::value( a, coords[0] ) *
                          LagrangeBasis3GL::value( b, coords[1] ) *
                          LagrangeBasis3GL::value( c, coords[2] );
          }
        }
      }
    }
  };
};


}
}


#endif /* GEOSX_FINITEELEMENT_ELEMENTFORMULATIONS_ELEMENTFORMULATIONS_LAGRANGEBASIS3GL_HPP_ */
