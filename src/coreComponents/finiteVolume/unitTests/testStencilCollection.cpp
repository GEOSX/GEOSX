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

// Source includes
#include "managers/initialization.hpp"
#include "finiteVolume/FluxStencil.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

// TPL includes
#include <gtest/gtest.h>

// System includes
#include <chrono>
#include <iostream>

#define TEST_SIZE 1000000

using namespace geosx;

struct Cell
{
  GEOSX_HOST_DEVICE
  Cell() :
    er( 0 ),
    esr( 0 ),
    ei( 0 )
  {}
  GEOSX_HOST_DEVICE
  Cell( localIndex i ) :
    er( 0 ),
    esr( 0 ),
    ei( i )
  {}

  localIndex er;
  localIndex esr;
  localIndex ei;
};

template< typename INDEX, typename T >
void
makeStencilTPFA( localIndex size, FluxStencil< INDEX, T > & stencil )
{
  stencil.reserve( size, 2 );

  stackArray1d< INDEX, 2 > cells( 2 );
  stackArray1d< real64, 2 > weights( 2 );

  for( localIndex kf = 0; kf < size; ++kf )
  {
    cells[0] = INDEX( kf );
    cells[1] = INDEX( kf + 1 );
    weights[0] = 1.0 / 2.0;
    weights[1] = 1.0 / 2.0;

    stencil.add( 2, cells.data(), weights.data(), kf );
  }
  stencil.compress();
}

template< typename LAMBDA >
void
testStencilLoop( LAMBDA && compute )
{
  FluxStencil< Cell, double > stencil;
  makeStencilTPFA( TEST_SIZE, stencil );

  array1d< double > in( TEST_SIZE + 1 );
  array1d< double > out( TEST_SIZE );
  for( localIndex i = 0; i < TEST_SIZE; ++i )
  {
    in[i] = double( i );
  }

  auto start = std::chrono::high_resolution_clock::now();
  compute( in, out, stencil );
  auto end = std::chrono::high_resolution_clock::now();

  std::chrono::duration< double > diff = end - start;
  std::cout << "Time elapsed: " << diff.count() << " s" << std::endl;

  for( localIndex kf = 0; kf < TEST_SIZE; ++kf )
  {
    EXPECT_FLOAT_EQ( out[kf], ( in[kf] / 2 + in[kf + 1] / 2 ) );
  }
}

TEST( testStencilCollection, noAcessorIterationTPFA )
{
  testStencilLoop( [&]( arrayView1d< double const > const & dataIn,
                        arrayView1d< double > const & dataOut,
                        FluxStencil< Cell, double > const & stencil ) {
    ArrayOfArraysView< FluxStencil< Cell, double >::Entry const, true > const & connections =
      stencil.getConnections();

    forAll< serialPolicy >( connections.size(), [=]( localIndex iconn ) {
      for( localIndex i = 0; i < connections.sizeOfArray( iconn ); ++i )
      {
        FluxStencil< Cell, double >::Entry const & entry = connections( iconn, i );
        dataOut[iconn] += entry.weight * dataIn[entry.index.ei];
      }
    } );
  } );
}

int
main( int argc, char * argv[] )
{
  geosx::basicSetup( argc, argv );

  int result = 0;
  testing::InitGoogleTest( &argc, argv );
  result = RUN_ALL_TESTS();

  geosx::basicCleanup();
  return result;
}
