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
 * @file H1_Pyramid_Lagrange1_Gauss5.hpp
 */

#include <finiteElement/elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp>
#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void
testKernelDriver()
{
  constexpr int numNodes = 5;
  constexpr int numQuadraturePoints = 5;
  constexpr static real64 weight = 81.0 / 100.0;
  constexpr static real64 weightDelta = 125.0 / 27.0 - weight;
  constexpr static real64 quadratureCrossSectionCoord = 0.584237394672177;
  constexpr static real64 quadratureLongitudinalCoordNeg = -2.0 / 3.0;
  constexpr static real64 quadratureLongitudinalCoordDelta = 16.0 / 15.0;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;

  constexpr real64 xCoords[numNodes][3] = { { 0.0, 0.0, 0.0 },
                                            { 1.0, 0.0, 0.0 },
                                            { 0.0, 1.0, 0.0 },
                                            { 1.0, 1.0, 0.0 },
                                            { 0.0, 0.0, 1.0 } };

  forAll< POLICY >( 1, [=] GEOSX_HOST_DEVICE( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = { 0 };
      H1_Pyramid_Lagrange1_Gauss5::shapeFunctionValues( q, N );
      for( localIndex a = 0; a < numNodes; ++a )
      {
        viewN( q, a ) = N[a];
      }
    }
  } );

  forAll< POLICY >( 1, [=] GEOSX_HOST_DEVICE( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 dNdX[numNodes][3] = { { 0 } };
      viewDetJ[q] =
        H1_Pyramid_Lagrange1_Gauss5::shapeFunctionDerivatives( q, xCoords, dNdX );

      for( localIndex a = 0; a < numNodes; ++a )
      {
        for( int i = 0; i < 3; ++i )
        {
          viewdNdX( q, a, i ) = dNdX[a][i];
        }
      }
    }
  } );

  constexpr real64 parentCoords[3][numNodes] = { { -1, 1, -1, 1, 0 },
                                                 { -1, -1, 1, 1, 0 },
                                                 { -1, -1, -1, -1, 1 } };

  forAll< serialPolicy >( 1, [=]( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 const xi[3] = {
        quadratureCrossSectionCoord * parentCoords[0][q],
        quadratureCrossSectionCoord * parentCoords[1][q],
        quadratureLongitudinalCoordNeg +
          quadratureLongitudinalCoordDelta * 0.5 * ( 1 + parentCoords[2][q] ) };

      for( localIndex a = 0; a < numNodes - 1; ++a )
      {
        real64 N = 0.125 * ( 1 + xi[0] * parentCoords[0][a] ) *
          ( 1 + xi[1] * parentCoords[1][a] ) * ( 1 + xi[2] * parentCoords[2][a] );
        EXPECT_FLOAT_EQ( N, viewN[q][a] );
      }

      {
        real64 N = 0.5 * ( 1 + xi[2] * parentCoords[2][numNodes - 1] );
        EXPECT_FLOAT_EQ( N, viewN[q][numNodes - 1] );
      }

      real64 J[3][3] = { { 0 } };
      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 dNdXi[3] = {
          0.5 * ( 1.0 - parentCoords[2][a] ) *
            ( 0.125 * parentCoords[0][a] * ( 1 + xi[1] * parentCoords[1][a] ) *
              ( 1 + xi[2] * parentCoords[2][a] ) ),
          0.5 * ( 1.0 - parentCoords[2][a] ) *
            ( 0.125 * ( 1 + xi[0] * parentCoords[0][a] ) * parentCoords[1][a] *
              ( 1 + xi[2] * parentCoords[2][a] ) ),
          0.5 * ( 1.0 - parentCoords[2][a] ) *
              ( 0.125 * ( 1 + xi[0] * parentCoords[0][a] ) *
                ( 1 + xi[1] * parentCoords[1][a] ) * parentCoords[2][a] ) +
            0.5 * ( 1.0 + parentCoords[2][a] ) * ( 0.5 * parentCoords[2][a] ) };

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            J[i][j] = J[i][j] + xCoords[a][i] * dNdXi[j];
          }
        }
      }
      real64 const detJ = FiniteElementBase::inverse( J );
      EXPECT_FLOAT_EQ( detJ * ( weight + 0.25 * ( q & 4 ) * weightDelta ),
                       viewDetJ[q] );

      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 dNdX[3] = { 0 };
        real64 dNdXi[3] = {
          0.5 * ( 1.0 - parentCoords[2][a] ) *
            ( 0.125 * parentCoords[0][a] * ( 1 + xi[1] * parentCoords[1][a] ) *
              ( 1 + xi[2] * parentCoords[2][a] ) ),
          0.5 * ( 1.0 - parentCoords[2][a] ) *
            ( 0.125 * ( 1 + xi[0] * parentCoords[0][a] ) * parentCoords[1][a] *
              ( 1 + xi[2] * parentCoords[2][a] ) ),
          0.5 * ( 1.0 - parentCoords[2][a] ) *
              ( 0.125 * ( 1 + xi[0] * parentCoords[0][a] ) *
                ( 1 + xi[1] * parentCoords[1][a] ) * parentCoords[2][a] ) +
            0.5 * ( 1.0 + parentCoords[2][a] ) * ( 0.5 * parentCoords[2][a] ) };

        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 3; ++j )
          {
            dNdX[i] += dNdXi[j] * J[j][i];
          }
        }

        EXPECT_FLOAT_EQ( dNdX[0], viewdNdX[q][a][0] );
        EXPECT_FLOAT_EQ( dNdX[1], viewdNdX[q][a][1] );
        EXPECT_FLOAT_EQ( dNdX[2], viewdNdX[q][a][2] );
      }
    }
  } );
}

#ifdef USE_CUDA
TEST( FiniteElementShapeFunctions, testKernelCuda )
{
  testKernelDriver< geosx::parallelDevicePolicy< 32 > >();
}
#endif
TEST( FiniteElementShapeFunctions, testKernelHost )
{
  testKernelDriver< serialPolicy >();
}

using namespace geosx;
int
main( int argc, char * argv[] )
{
  testing::InitGoogleTest();

  basicSetup( argc, argv, false );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
