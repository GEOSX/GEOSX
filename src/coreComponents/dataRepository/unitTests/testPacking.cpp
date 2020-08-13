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

#include "LvArray/src/Array.hpp"
#include "dataRepository/BufferOps.hpp"
#include "dataRepository/BufferOpsDevice.hpp"
#include "managers/initialization.hpp"
#include "mpiCommunications/CommunicationTools.hpp"

#include <ctime>
#include <cstdlib>

using namespace geosx;

#ifndef GTEST_SKIP
  #define GTEST_SKIP() return
#endif

#define SKIP_TEST_IF( COND, REASON )                                          \
  do                                                                          \
  {                                                                           \
    if( COND )                                                                \
    {                                                                         \
      GEOSX_WARNING( "This test is currently known to fail when " #COND       \
                     " because:\n" REASON                                     \
                     "\n"                                                     \
                     "Therefore, we skip it entirely for this run (may show " \
                     "as PASSED or SKIPPED)" );                               \
      GTEST_SKIP();                                                           \
    }                                                                         \
  } while( 0 )

#define SKIP_TEST_IN_SERIAL( REASON )                            \
  do                                                             \
  {                                                              \
    int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize == 1, REASON );                        \
  } while( 0 )

#define SKIP_TEST_IN_PARALLEL( REASON )                          \
  do                                                             \
  {                                                              \
    int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX ); \
    SKIP_TEST_IF( mpiSize != 1, REASON );                        \
  } while( 0 )

real64
drand( real64 min = 0.0, real64 max = 1.0 )
{
  real64 f = (real64)rand() / RAND_MAX;
  return min + f * ( max - min );
}

TEST( testPacking, testPacking )
{
  std::srand( std::time( nullptr ) );
  constexpr localIndex size = 10000;
  array1d< R1Tensor > veloc( size );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj ) veloc[ii][jj] = drand();

  buffer_unit_type * null_buf = NULL;
  localIndex calc_size = bufferOps::Pack< false >( null_buf, veloc );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::Pack< true >( buffer, veloc );
  buffer_unit_type const * cbuffer = &buf[0];
  bufferOps::Unpack( cbuffer, unpacked );
  for( localIndex ii = 0; ii < size; ++ii )
    EXPECT_TRUE( veloc[ii] == unpacked[ii] );
}

TEST( testPacking, testPackByIndex )
{
  std::srand( std::time( nullptr ) );
  constexpr localIndex size = 10000;
  localIndex pack_count = std::rand() % size;
  array1d< R1Tensor > veloc( size );
  array1d< localIndex > indices( pack_count );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj ) veloc[ii][jj] = drand();

  for( localIndex ii = 0; ii < pack_count; ++ii )
    indices[ii] = std::rand() % size;

  // std::sort(indices.begin(),indices.end())

  buffer_unit_type * null_buf = NULL;
  // [ ndim, din_i_stride, .. , num_tens_in_arr_i, tens_len_i, data_i, ... ]
  localIndex calc_size = bufferOps::PackByIndex< false >( null_buf, veloc, indices );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::PackByIndex< true >( buffer, veloc, indices );
  buffer_unit_type const * cbuffer = &buf[0];
  bufferOps::UnpackByIndex( cbuffer, unpacked, indices );
  for( localIndex ii = 0; ii < size; ++ii )
  {
    if( std::find( indices.begin(), indices.end(), ii ) != indices.end() )
    {
      EXPECT_EQ( veloc[ii], unpacked[ii] );
    }
  }
}

TEST( testPacking, testTensorPacking )
{
  std::srand( std::time( nullptr ) );
  array1d< R1Tensor > tns( 1 );
  for( localIndex ii = 0; ii < 3; ++ii ) tns[0][ii] = drand();

  buffer_unit_type * null_buf = nullptr;
  localIndex calc_size =
    bufferOps::PackDevice< false >( null_buf, tns.toViewConst() );
  buffer_type buf( calc_size );
  buffer_unit_type * b = &buf[0];
  bufferOps::PackDevice< true >( b, tns.toViewConst() );

  array1d< R1Tensor > unp( 1 );
  buffer_unit_type const * bc = &buf[0];
  bufferOps::UnpackDevice( bc, unp.toView() );
  unp.move( LvArray::MemorySpace::CPU );
  for( localIndex ii = 0; ii < 3; ++ii ) EXPECT_TRUE( tns[0][ii] = unp[0][ii] );
}

TEST( testPacking, testPackingDevice )
{
  std::srand( std::time( nullptr ) );
  constexpr localIndex size = 10000;
  array1d< R1Tensor > veloc( size );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj ) veloc[ii][jj] = drand();

  buffer_unit_type * null_buf = NULL;
  localIndex calc_size =
    bufferOps::PackDevice< false >( null_buf, veloc.toViewConst() );

  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  bufferOps::PackDevice< true >( buffer, veloc.toViewConst() );

  buffer_unit_type const * cbuffer = &buf[0];
  bufferOps::UnpackDevice( cbuffer, unpacked.toView() );
  unpacked.move( LvArray::MemorySpace::CPU );
  for( localIndex ii = 0; ii < size; ++ii ) EXPECT_EQ( veloc[ii], unpacked[ii] );
}

TEST( testPacking, testPackByIndexDevice )
{
  std::srand( std::time( nullptr ) );
  constexpr localIndex size = 10000;
  localIndex pack_count = std::rand() % size;
  array1d< R1Tensor > veloc( size );
  array1d< localIndex > indices( pack_count );
  array1d< R1Tensor > unpacked( size );

  for( localIndex ii = 0; ii < size; ++ii )
    for( localIndex jj = 0; jj < 3; ++jj ) veloc[ii][jj] = drand();

  for( localIndex ii = 0; ii < pack_count; ++ii )
    indices[ii] = std::rand() % size;

  // std::sort(indices.begin(),indices.end())

  buffer_unit_type * null_buf = NULL;
  // [ num_dim, stride_i.. , tensor_0, tensor_1, ..., tensor_n ]
  arrayView1d< R1Tensor const > const & veloc_view = veloc.toViewConst();
  localIndex calc_size =
    bufferOps::PackByIndexDevice< false >( null_buf,
                                           veloc_view,
                                           indices.toViewConst() );
  buffer_type buf( calc_size );
  buffer_unit_type * buffer = &buf[0];
  localIndex packed_size =
    bufferOps::PackByIndexDevice< true >( buffer, veloc_view, indices.toViewConst() );
  EXPECT_EQ( calc_size, packed_size );
  buffer_unit_type const * cbuffer = &buf[0];
  localIndex unpacked_size =
    bufferOps::UnpackByIndexDevice( cbuffer,
                                    unpacked.toView(),
                                    indices.toViewConst() );
  EXPECT_EQ( unpacked_size, packed_size );
  unpacked.move( LvArray::MemorySpace::CPU );
  for( localIndex ii = 0; ii < size; ++ii )
  {
    if( std::find( indices.begin(), indices.end(), ii ) != indices.end() )
    {
      EXPECT_EQ( veloc[ii], unpacked[ii] );
    }
  }
}

int
main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
