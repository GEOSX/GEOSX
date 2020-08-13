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

#include "dataRepository/BufferOps.hpp"

using namespace geosx;
using namespace bufferOps;

TEST( testGeosxTraits, test_is_noncontainer_type_packable )
{
  static_assert( is_noncontainer_type_packable< int >, "Should be true." );
  static_assert( is_noncontainer_type_packable< double >, "Should be true." );
  static_assert( is_noncontainer_type_packable< R1Tensor >, "Should be true." );
  static_assert( is_noncontainer_type_packable< string >, "Should be true." );

  static_assert( !is_noncontainer_type_packable< void >, "Should be false." );
  static_assert( !is_noncontainer_type_packable< array1d< double > >,
                 "Should be false." );
  static_assert( !is_noncontainer_type_packable< SortedArray< double > >,
                 "Should be false." );
  static_assert( !is_noncontainer_type_packable< map< string, int > >,
                 "Should be false." );
  static_assert( !is_noncontainer_type_packable< std::pair< string, int > >,
                 "Should be false." );
}

TEST( testGeosxTraits, test_is_array_packable )
{
  static_assert( is_packable_array< array2d< real64, RAJA::PERM_IJ > >,
                 "Should be true." );
  static_assert( is_packable_array< array2d< real64, RAJA::PERM_JI > >,
                 "Should be true." );

  static_assert( !is_packable_array< int >, "Should be false." );
  static_assert( !is_packable_array< double >, "Should be false." );
  static_assert( !is_packable_array< void >, "Should be false." );
}

TEST( testGeosxTraits, test_is_packable_map )
{
  static_assert( is_packable_map< map< string, int > >, "Should be true." );
  static_assert( is_packable_map< map< string, array1d< int > > >, "Should be true." );
  static_assert( !is_packable_map< map< string, std::pair< int, int > > >,
                 "Should be false" );
}
