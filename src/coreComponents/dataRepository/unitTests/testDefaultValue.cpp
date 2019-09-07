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

#include "Array.hpp"
#include "dataRepository/DefaultValue.hpp"

#include <functional>
#include <string>
#include <typeindex>
#include <vector>

using namespace geosx;
using namespace LvArray;
using namespace dataRepository;


TEST( testDefaultValue, testScalar )
{
EXPECT_TRUE( DefaultValue< int >::has_default_value == true );
EXPECT_TRUE( DefaultValue< long int >::has_default_value == true );
EXPECT_TRUE( DefaultValue< long long int >::has_default_value == true );
EXPECT_TRUE( DefaultValue< double >::has_default_value == true );
}

TEST( testDefaultValue, testArray )
{
using array1 = Array< double, 1, int >;
using array2 = Array< double, 2, int >;
using array3 = Array< double, 3, int >;
using array4 = Array< int, 1, int >;
using array5 = Array< long int, 1, long int >;
using array6 = Array< long long int, 1, long long int >;
EXPECT_TRUE( DefaultValue< array1 >::has_default_value==true );
EXPECT_TRUE( DefaultValue< array2 >::has_default_value==true );
EXPECT_TRUE( DefaultValue< array3 >::has_default_value==true );
EXPECT_TRUE( DefaultValue< array4 >::has_default_value==true );
EXPECT_TRUE( DefaultValue< array5 >::has_default_value==true );
EXPECT_TRUE( DefaultValue< array6 >::has_default_value==true );
}
