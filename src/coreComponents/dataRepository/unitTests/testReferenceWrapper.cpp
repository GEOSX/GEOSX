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

#include <gtest/gtest.h>

#include "dataRepository/ReferenceWrapper.hpp"
#include "Array.hpp"

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace LvArray;

TEST( testReferenceWrapper, testIntWrapper )
{
  int var = 0;

  ReferenceWrapper< int > wrappedVar( var );

  wrappedVar = 5;

  int var2 = wrappedVar;

  EXPECT_TRUE( var == 5 );
  EXPECT_TRUE( var2 == var );

  EXPECT_TRUE( &var == &( wrappedVar.get() ) );
  EXPECT_TRUE( &var2 != &( wrappedVar.get() ) );
  EXPECT_TRUE( &var != &var2 );
}

TEST( testReferenceWrapper, testArrayWrapper )
{
  Array< int, 1, int > arr;
  arr.resize( 4 );

  ReferenceWrapper< Array< int, 1, int > > wrappedArr( arr );

  for( int i = 0; i < 4; ++i )
  {
    wrappedArr[i] = 2 * i;
  }

  EXPECT_TRUE( arr[0] == 0 );
  EXPECT_TRUE( arr[1] == 2 );
  EXPECT_TRUE( arr[2] == 4 );
  EXPECT_TRUE( arr[3] == 6 );
}

TEST( testReferenceWrapper, testArrayOfWrappedInts )
{
  int val0 = 0;
  int val1 = 1;

  std::vector< ReferenceWrapper< int > > arr2;
  arr2.resize( 2 );
  arr2[0].set( val0 );
  arr2[1].set( val1 );

  arr2[0] = 10;
  arr2[1] = 11;

  EXPECT_TRUE( val0 == 10 );
  EXPECT_TRUE( val1 == 11 );
}

TEST( testReferenceWrapper, testOperatorParen )
{
  using array2d = Array< int, 2, int >;
  array2d arr;
  arr.resize( 2, 3 );

  ReferenceWrapper< array2d > wrappedArr( arr );

  for( int i = 0; i < 2; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      wrappedArr( i, j ) = 3 * i + j;
    }
  }

  EXPECT_TRUE( arr[0][0] == 0 );
  EXPECT_TRUE( arr[0][1] == 1 );
  EXPECT_TRUE( arr[0][2] == 2 );
  EXPECT_TRUE( arr[1][0] == 3 );
  EXPECT_TRUE( arr[1][1] == 4 );
  EXPECT_TRUE( arr[1][2] == 5 );
}

TEST( testReferenceWrapper, testNestedOperatorSquare )
{
  using array2d = Array< int, 2, int >;
  array2d arr;
  arr.resize( 2, 3 );

  ReferenceWrapper< array2d > wrappedArr( arr );

  for( int i = 0; i < 2; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      wrappedArr[i][j] = 3 * i + j;
    }
  }

  EXPECT_TRUE( arr[0][0] == 0 );
  EXPECT_TRUE( arr[0][1] == 1 );
  EXPECT_TRUE( arr[0][2] == 2 );
  EXPECT_TRUE( arr[1][0] == 3 );
  EXPECT_TRUE( arr[1][1] == 4 );
  EXPECT_TRUE( arr[1][2] == 5 );
}
