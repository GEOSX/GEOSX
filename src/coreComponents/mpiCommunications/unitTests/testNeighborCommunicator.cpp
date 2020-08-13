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

#include "managers/initialization.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "LvArray/src/Array.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#ifdef UMPIRE_ENABLE_CUDA
  #include "umpire/alloc/CudaPinnedAllocator.hpp"
#endif

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

char
crand()
{
  return 'a' + rand() % 26;
}

TEST( TestNeighborComms, testBuffers )
{
  {
    size_t sz = 100;
    std::vector< char > sd( sz );
    for( size_t ii = 0; ii < sz; ++ii ) sd[ii] = crand();

    auto nc = NeighborCommunicator();
    nc.resizeSendBuffer( 0, sz );

    auto & sb = nc.SendBuffer( 0 );
    for( size_t ii = 0; ii < sz; ++ii ) sb[ii] = sd[ii];

    for( size_t ii = 0; ii < sz; ++ii ) EXPECT_EQ( sb[ii], sd[ii] );
  }
}

#if defined( UMPIRE_ENABLE_CUDA ) && defined( USE_CHAI )
void
pack( buffer_unit_type * buf, arrayView1d< const int > & veloc_view, localIndex size )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex ii ) {
    reinterpret_cast< int * >( buf )[ii] = veloc_view.data()[ii];
  } );
}

void
unpack( buffer_unit_type * buf, arrayView1d< int > & veloc_view, localIndex size )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex ii ) {
    veloc_view.data()[ii] = reinterpret_cast< int * >( buf )[ii];
  } );
}

TEST( TestNeighborComms, testMPICommunication_fromPinnedSetOnDevice )
{
  SKIP_TEST_IN_SERIAL( "Parallel test" );
  {
    int rnk = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

    constexpr localIndex size = 1000;
    constexpr localIndex byte_size = 1000 * sizeof( int );
    array1d< int > veloc( size );
    for( int ii = 0; ii < size; ++ii ) veloc[ii] = rnk == 0 ? ii : 0;

    auto nc = NeighborCommunicator();
    nc.resizeSendBuffer( 0, byte_size );
    auto & sb = nc.SendBuffer( 0 );
    buffer_unit_type * buf = &sb[0];

    MPI_Request request;
    if( rnk == 0 )
    {
      veloc.move( LvArray::MemorySpace::GPU );
      auto veloc_view = veloc.toViewConst();
      pack( buf, veloc_view, size );
      MpiWrapper::iSend( buf, byte_size, 1, 0, MPI_COMM_GEOSX, &request );
    }
    else
    {
      int err = MpiWrapper::iRecv( buf, byte_size, 0, 0, MPI_COMM_GEOSX, &request );
      EXPECT_EQ( err, MPI_SUCCESS );
      MPI_Status status;
      err = MpiWrapper::Wait( &request, &status );
      EXPECT_EQ( err, MPI_SUCCESS );
      auto veloc_view = veloc.toView();
      unpack( buf, veloc_view, size );
      veloc.move( LvArray::MemorySpace::CPU );
      for( int ii = 0; ii < size; ++ii ) EXPECT_EQ( veloc[ii], ii );
    }
  }
}
#endif

int
main( int ac, char * av[] )
{
  ::testing::InitGoogleTest( &ac, av );
  geosx::basicSetup( ac, av );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
