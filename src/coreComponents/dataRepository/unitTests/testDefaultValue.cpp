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
#include "dataRepository/DefaultValue.hpp"

// TPL includes
#include <gtest/gtest.h>

// System includes
#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace dataRepository;

TEST( testDefaultValue, testScalar )
{
  // The comparison with true is to avoid a linker error because has_default_value is constexpr.
  EXPECT_TRUE( DefaultValue< int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< long int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< long long int >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< double >::has_default_value == true );
}

TEST( testDefaultValue, testArray )
{
  // The comparison with true is to avoid a linker error because has_default_value is constexpr.
  EXPECT_TRUE( DefaultValue< array1d< double > >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< array2d< double > >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< array3d< double > >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< array1d< int > >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< array1d< long int > >::has_default_value == true );
  EXPECT_TRUE( DefaultValue< array1d< long long int > >::has_default_value == true );
}
