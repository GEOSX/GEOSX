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
 * @file testFiniteElementBase.cpp
 */

#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

#include "finiteElement/elementFormulations/FiniteElementBase.hpp"

using namespace geosx;
using namespace finiteElement;

static real64 inverse( real64 ( &J )[3][3] )
{
  real64 scratch[3][3];
  scratch[0][0] = J[1][1] * J[2][2] - J[1][2] * J[2][1];
  scratch[1][0] = J[0][2] * J[2][1] - J[0][1] * J[2][2];
  scratch[2][0] = J[0][1] * J[1][2] - J[0][2] * J[1][1];
  scratch[0][1] = J[1][2] * J[2][0] - J[1][0] * J[2][2];
  scratch[1][1] = J[0][0] * J[2][2] - J[0][2] * J[2][0];
  scratch[2][1] = J[0][2] * J[1][0] - J[0][0] * J[1][2];
  scratch[0][2] = J[1][0] * J[2][1] - J[1][1] * J[2][0];
  scratch[1][2] = J[0][1] * J[2][0] - J[0][0] * J[2][1];
  scratch[2][2] = J[0][0] * J[1][1] - J[0][1] * J[1][0];

  real64 const detJ =
    J[0][0] * scratch[0][0] + J[1][0] * scratch[1][0] + J[2][0] * scratch[2][0];
  real64 const invDet = 1 / detJ;

  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      J[i][j] = scratch[j][i] * invDet;
    }
  }

  return detJ;
}

struct Jacobian
{
  real64 data[3][3] = { {
                          1.19167,
                          -0.0372008,
                          -0.0766346,
                        },
                        {
                          0.0599679,
                          1.19167,
                          0.0205342,
                        },
                        {
                          -0.0438996,
                          0.00610042,
                          1.1378,
                        } };
};

template< typename POLICY >
void
testInverseDriver()
{
  Jacobian J;
  Jacobian invJ;
  real64 detJ;

  forAll< POLICY >( 1, [&]( localIndex const ) {
    detJ = FiniteElementBase::inverse( invJ.data );
  } );

  real64 const detJ_Ref = inverse( J.data );

  EXPECT_FLOAT_EQ( detJ, detJ_Ref );

  forAll< serialPolicy >( 1, [=]( localIndex const ) {
    for( int i = 0; i < 3; ++i )
    {
      for( int j = 0; j < 3; ++j )
      {
        EXPECT_FLOAT_EQ( J.data[i][j], invJ.data[i][j] );
      }
    }
  } );
}

template< typename POLICY >
void
testDetJDriver()
{
  Jacobian J;
  real64 detJ;

  forAll< POLICY >( 1, [&]( localIndex const ) {
    detJ = FiniteElementBase::detJ( J.data );
  } );

  real64 const detJ_Ref = inverse( J.data );

  EXPECT_FLOAT_EQ( detJ, detJ_Ref );
}

TEST( FiniteElementBase, testJacobianInverseHost )
{
  testInverseDriver< serialPolicy >();
}
TEST( FiniteElementBase, testDetJHost )
{
  testDetJDriver< serialPolicy >();
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
