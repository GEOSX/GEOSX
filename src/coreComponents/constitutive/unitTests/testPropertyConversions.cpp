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

#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PropertyConversions.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( PropertyConversionTests, testElasticConversions )
{
  real64 bulk, poisson, young, shear;
  real64 bulk_out, poisson_out, young_out, shear_out;

  // bulk and poisson
  bulk        = 5.4321;
  poisson     = 0.2345;
  young      = conversions::BulkModAndPoissonRatio::toYoungMod( bulk, poisson );
  shear       = conversions::BulkModAndPoissonRatio::toShearMod( bulk, poisson );

  EXPECT_DOUBLE_EQ( young, 3*bulk*(1-2*poisson) );
  EXPECT_DOUBLE_EQ( shear, 3*bulk*(1-2*poisson)/(2+2*poisson) );

  // bulk and shear
  poisson_out = conversions::BulkModAndShearMod::toPoissonRatio( bulk, shear );
  young_out = conversions::BulkModAndShearMod::toYoungMod( bulk, shear );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // bulk and young
  shear_out = conversions::BulkModAndYoungMod::toShearMod( bulk, young );
  poisson_out = conversions::BulkModAndYoungMod::toPoissonRatio( bulk, young );

  EXPECT_DOUBLE_EQ( poisson, poisson_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // young and shear
  bulk_out    = conversions::ShearModAndYoungMod::toBulkMod( shear, young );
  poisson_out = conversions::ShearModAndYoungMod::toPoissonRatio( shear, young );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( poisson, poisson_out );

  // young and poisson
  bulk_out = conversions::YoungModAndPoissonRatio::toBulkMod( young, poisson );
  shear_out = conversions::YoungModAndPoissonRatio::toShearMod( young, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( shear, shear_out );

  // shear and poisson
  bulk_out = conversions::ShearModAndPoissonRatio::toBulkMod( shear, poisson );
  young_out =conversions::ShearModAndPoissonRatio::toYoungMod( shear, poisson );

  EXPECT_DOUBLE_EQ( bulk, bulk_out );
  EXPECT_DOUBLE_EQ( young, young_out );

  // first lame (bulk and shear input only)
  real64 lame = conversions::BulkModAndShearMod::toFirstLame( bulk, shear );

  EXPECT_DOUBLE_EQ( lame, bulk-2*shear/3 );
}
