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
#include "constitutive/solid/TriaxialDriver.hpp"

#include "dataRepository/xmlWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

template< typename POLICY >
void testTriaxialDriver()
{
  TriaxialDriver driver;
}


#ifdef USE_CUDA
TEST( TriaxialTests, testTriaxialDevice )
{
  testCamClayDriver< geosx::parallelDevicePolicy< > >();
}
#endif
TEST( TriaxialTests, testTriaxialHost )
{
  testCamClayDriver< serialPolicy >();
}
