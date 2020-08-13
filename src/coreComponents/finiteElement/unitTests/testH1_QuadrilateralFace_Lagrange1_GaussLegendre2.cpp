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
 * @file testH1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp
 */

#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void
testKernelDriver()
{
  constexpr int numNodes = 4;
  constexpr int numQuadraturePoints = 4;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;

  constexpr real64 xCoords[numNodes][3] = { { 0.1, 0.1, 0.9 },
                                            { 1.1, -0.1, 0.1 },
                                            { 0.0, 1.1, 0.0 },
                                            { 0.9, 1.1, -0.1 } };

  forAll< POLICY >( 1, [=] GEOSX_HOST_DEVICE( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = { 0 };
      H1_QuadrilateralFace_Lagrange1_GaussLegendre2::shapeFunctionValues( q, N );
      for( localIndex a = 0; a < numNodes; ++a )
      {
        viewN( q, a ) = N[a];
      }
    }
  } );

  forAll< POLICY >( 1, [=] GEOSX_HOST_DEVICE( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      viewDetJ[q] =
        H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight(
          q,
          xCoords );
    }
  } );

  constexpr real64 pCoords[2][numNodes] = { { -1, 1, -1, 1 }, { -1, -1, 1, 1 } };

  constexpr static real64 quadratureFactor = 1.0 / 1.732050807568877293528;
  constexpr static real64 weight = 1.0;

  forAll< serialPolicy >( 1, [=]( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 const xi[2] = { quadratureFactor * pCoords[0][q],
                             quadratureFactor * pCoords[1][q] };

      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 N =
          0.25 * ( 1 + xi[0] * pCoords[0][a] ) * ( 1 + xi[1] * pCoords[1][a] );
        EXPECT_FLOAT_EQ( N, viewN[q][a] );
      }

      real64 dXdXi[3][2] = { { 0 } };
      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 dNdXi[2] = { 0.25 * pCoords[0][a] * ( 1 + xi[1] * pCoords[1][a] ),
                            0.25 * ( 1 + xi[0] * pCoords[0][a] ) * pCoords[1][a] };
        for( int i = 0; i < 3; ++i )
        {
          for( int j = 0; j < 2; ++j )
          {
            dXdXi[i][j] = dXdXi[i][j] + xCoords[a][i] * dNdXi[j];
          }
        }
      }

      real64 n[3] = { dXdXi[1][0] * dXdXi[2][1] - dXdXi[2][0] * dXdXi[1][1],
                      dXdXi[2][0] * dXdXi[0][1] - dXdXi[0][0] * dXdXi[2][1],
                      dXdXi[0][0] * dXdXi[1][1] - dXdXi[1][0] * dXdXi[0][1] };
      real64 const detJ = sqrt( n[0] * n[0] + n[1] * n[1] + n[2] * n[2] );
      EXPECT_FLOAT_EQ( detJ * weight, viewDetJ[q] );
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
