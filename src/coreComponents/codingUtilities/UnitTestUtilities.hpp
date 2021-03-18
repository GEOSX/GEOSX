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

/*
 * @file UnitTestUtilities.hpp
 */

#ifndef GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_
#define GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_

#include "gtest/gtest.h"

#include "common/DataTypes.hpp"

#ifndef GTEST_SKIP
#define GTEST_SKIP() return
#endif

#define SKIP_TEST_IF( COND, REASON ) \
  do \
  { \
    if( COND ) \
    { \
      GEOSX_WARNING( "This test is currently known to fail when " #COND " because:\n" REASON "\n" \
                                                                                             "Therefore, we skip it entirely for this run (may show as PASSED or SKIPPED)" ); \
      GTEST_SKIP(); \
    } \
  } while(0)

#define SKIP_TEST_IN_SERIAL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize == 1, REASON ); \
  } while(0)

#define SKIP_TEST_IN_PARALLEL( REASON ) \
  do \
  { \
    int const mpiSize = MpiWrapper::commSize( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize > 1, REASON ); \
  } while(0)

namespace geosx
{

namespace testing
{

constexpr real64 DEFAULT_ABS_TOL = 1E-13;
constexpr real64 DEFAULT_REL_TOL = std::numeric_limits< real64 >::epsilon();

::testing::AssertionResult checkRelativeErrorFormat( const char *, const char *, const char *, const char *,
                                                     real64 const v1, real64 const v2, real64 const relTol, real64 const absTol )
{
  real64 const delta = std::abs( v1 - v2 );
  real64 const value = std::max( std::abs( v1 ), std::abs( v2 ) );
  if( delta > absTol && delta > relTol * value )
  {
    return ::testing::AssertionFailure() << std::scientific << std::setprecision( 5 )
                                         << " relative error: " << delta / value
                                         << " (" << v1 << " vs " << v2 << "),"
                                         << " exceeds " << relTol << std::endl;
  }
  return ::testing::AssertionSuccess();
}

::testing::AssertionResult checkRelativeErrorFormat( char const * a, char const * b, char const * c,
                                                     real64 const v1, real64 const v2, real64 const relTol )
{ return checkRelativeErrorFormat( a, b, c, "DEFAULT_ABS_TOL", v1, v2, relTol, DEFAULT_ABS_TOL ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol )
{ EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, real64 const absTol )
{ EXPECT_PRED_FORMAT4( checkRelativeErrorFormat, v1, v2, relTol, absTol ); }

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, string const & name )
{
  SCOPED_TRACE( name );
  EXPECT_PRED_FORMAT3( checkRelativeErrorFormat, v1, v2, relTol );
}

void checkRelativeError( real64 const v1, real64 const v2, real64 const relTol, real64 const absTol, string const & name )
{
  SCOPED_TRACE( name );
  EXPECT_PRED_FORMAT4( checkRelativeErrorFormat, v1, v2, relTol, absTol );
}

template< typename ROW_INDEX, typename COL_INDEX, typename VALUE >
void compareMatrixRow( ROW_INDEX const rowNumber, VALUE const relTol, VALUE const absTol,
                       localIndex const length1, COL_INDEX const * const indices1, VALUE const * const values1,
                       localIndex const length2, COL_INDEX const * const indices2, VALUE const * const values2 )
{
  SCOPED_TRACE( "Row " + std::to_string( rowNumber ));

  EXPECT_EQ( length1, length2 );

  for( localIndex j1 = 0, j2 = 0; j1 < length1 && j2 < length2; ++j1, ++j2 )
  {
    while( j1 < length1 && j2 < length2 && indices1[j1] != indices1[j2] )
    {
      while( j1 < length1 && indices1[j1] < indices2[j2] )
      {
        ADD_FAILURE() << "column " << indices1[j1] << ") in matrix 1 does not have a match";
      }
      while( j2 < length2 && indices2[j2] < indices1[j1] )
      {
        ADD_FAILURE() << "column " << indices2[j2] << ") in matrix 2 does not have a match";
      }
    }
    if( j1 < length1 && j2 < length2 )
    {
      SCOPED_TRACE( "Column " + std::to_string( indices1[j1] ));

      checkRelativeError( values1[j1], values2[j1], relTol, absTol );
    }
  }
}

template< typename ROW_INDEX, typename COL_INDEX, typename VALUE >
void compareMatrixRow( ROW_INDEX const rowNumber, VALUE const relTol, VALUE const absTol,
                       arraySlice1d< COL_INDEX const > indices1, arraySlice1d< VALUE const > values1,
                       arraySlice1d< COL_INDEX const > indices2, arraySlice1d< VALUE const > values2 )
{
  ASSERT_EQ( indices1.size(), values1.size() );
  ASSERT_EQ( indices2.size(), values2.size() );

  compareMatrixRow( rowNumber, relTol, absTol,
                    indices1.size(), indices1.dataIfContiguous(), values1.dataIfContiguous(),
                    indices2.size(), indices2.dataIfContiguous(), values2.dataIfContiguous() );
}

template< typename MATRIX >
void compareMatrices( MATRIX const & matrix1,
                      MATRIX const & matrix2,
                      real64 const relTol = DEFAULT_REL_TOL,
                      real64 const absTol = DEFAULT_ABS_TOL )
{
  ASSERT_EQ( matrix1.numGlobalRows(), matrix2.numGlobalRows() );
  ASSERT_EQ( matrix1.numGlobalCols(), matrix2.numGlobalCols() );

  ASSERT_EQ( matrix1.numLocalRows(), matrix2.numLocalRows() );
  ASSERT_EQ( matrix1.numLocalCols(), matrix2.numLocalCols() );

  array1d< globalIndex > indices1, indices2;
  array1d< real64 > values1, values2;

  // check the accuracy across local rows
  for( globalIndex i = matrix1.ilower(); i < matrix1.iupper(); ++i )
  {
    indices1.resize( matrix1.rowLength( i ) );
    values1.resize( matrix1.rowLength( i ) );
    matrix1.getRowCopy( i, indices1, values1 );

    indices2.resize( matrix2.rowLength( i ) );
    values2.resize( matrix2.rowLength( i ) );
    matrix2.getRowCopy( i, indices2, values2 );

    compareMatrixRow( i, relTol, absTol,
                      indices1.size(), indices1.data(), values1.data(),
                      indices2.size(), indices2.data(), values2.data() );
  }
}

template< typename T, typename COL_INDEX >
void compareLocalMatrices( CRSMatrixView< T const, COL_INDEX const > const & matrix1,
                           CRSMatrixView< T const, COL_INDEX const > const & matrix2,
                           real64 const relTol = DEFAULT_REL_TOL,
                           real64 const absTol = DEFAULT_ABS_TOL )
{
  ASSERT_EQ( matrix1.numRows(), matrix2.numRows() );
  ASSERT_EQ( matrix1.numColumns(), matrix2.numColumns() );

  matrix1.move( LvArray::MemorySpace::host, false );
  matrix2.move( LvArray::MemorySpace::host, false );

  // check the accuracy across local rows
  for( localIndex i = 0; i < matrix1.numRows(); ++i )
  {
    compareMatrixRow( i, relTol, absTol,
                      matrix1.getColumns( i ), matrix1.getEntries( i ),
                      matrix2.getColumns( i ), matrix2.getEntries( i ) );
  }
}

} // namespace testing

} // namespace geosx

#endif //GEOSX_CODINGUTILITIES_UNITTESTUTILITIES_HPP_
