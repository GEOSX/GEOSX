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
 * @file testLinearTriangleFaceShapeFunctionKernel.cpp
 */

#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>
#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void
testKernelDriver()
{
  constexpr int numNodes = 3;
  constexpr int numQuadraturePoints = 1;
  constexpr real64 weight = 0.5;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;

  constexpr real64 xCoords[numNodes][3] = { { 0.1, 0.1, 0.9 },
                                            { 1.1, -0.1, 0.1 },
                                            { 0.0, 1.1, 0.0 } };

  forAll< POLICY >( 1, [=] GEOSX_HOST_DEVICE( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = { 0 };
      H1_TriangleFace_Lagrange1_Gauss1::shapeFunctionValues( q, N );
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
        H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, xCoords );
    }
  } );

  constexpr real64 quadratureCoords[2][numQuadraturePoints] = { { 1.0 / 3.0 },
                                                                { 1.0 / 3.0 } };

  forAll< serialPolicy >( 1, [=]( localIndex const ) {
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      real64 const xi[2] = { quadratureCoords[0][q], quadratureCoords[1][q] };

      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 N = static_cast< real64 >( ( a | 0 ) < 1 ) +
          static_cast< real64 >( ( ( a ^ 1 ) < 1 ) - ( ( a ^ 1 ) == 1 ) ) * xi[0] +
          static_cast< real64 >( ( ( a ^ 2 ) < 1 ) - ( ( a ^ 2 ) == 2 ) ) * xi[1];
        EXPECT_FLOAT_EQ( N, viewN[q][a] );
      }

      real64 dXdXi[3][2] = { { 0 } };
      for( localIndex a = 0; a < numNodes; ++a )
      {
        real64 dNdXi[2] = { static_cast< real64 >( ( ( a ^ 1 ) < 1 ) - ( ( a ^ 1 ) == 1 ) ),
                            static_cast< real64 >( ( ( a ^ 2 ) < 1 ) - ( ( a ^ 2 ) == 2 ) ) };

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
