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

// Source includes
#include "mainInterface/initialization.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::SinglePhaseBaseKernels;

// Sphinx start after test mobility

TEST( SinglePhaseBaseKernels, mobility )
{
  int constexpr NTEST = 3;

  real64 const dens[NTEST]        = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST] = { 1e-5, 1e-10, 0.0    };
  real64 const visc[NTEST]        = { 5.0, 2.0, 1.0    };
  real64 const dVisc_dPres[NTEST] = { 1e-7, 0.0, 0.0    };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    real64 mob;
    real64 dMob_dPres;

    MobilityKernel::compute( dens[i], dDens_dPres[i], visc[i], dVisc_dPres[i], mob, dMob_dPres );

    // compute etalon
    real64 const mob_et = dens[i] / visc[i];
    real64 const dMob_dPres_et = mob_et * (dDens_dPres[i] / dens[i] - dVisc_dPres[i] / visc[i]);

    EXPECT_DOUBLE_EQ( mob, mob_et );
    EXPECT_DOUBLE_EQ( dMob_dPres, dMob_dPres_et );
  }
}

// Sphinx end before test mobility

/**
 * @brief Test the assembly of accumulation contribution
 *
 * @note This only tests uncoupled version.
 */
TEST( SinglePhaseBaseKernels, accumulation )
{
  int constexpr NTEST = 3;

  real64 const densOld[NTEST]       = { 700.0, 1200.0, 1500.0 };
  real64 const densNew[NTEST]       = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST]   = { 1e-5, 1e-10, 0.0   };
  real64 const poroOld[NTEST]       = { 6e-2, 1e-1, 2e-1   };
  real64 const poroNew[NTEST]       = { 7e-2, 2e-1, 3e-1   };
  real64 const dPoro_dPres[NTEST]   = { 0.0, 1e-1, 2e-1   };
  real64 const volume               = 1.0;


  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    real64 accum;
    real64 accumJacobian;

    real64 const poreVolumeNew = volume * poroNew[i];
    real64 const poreVolumeOld = volume * poroOld[i];
    real64 const dPoreVol_dPres = dPoro_dPres[i] * volume;


    AccumulationKernel::compute( densNew[i], densOld[i], dDens_dPres[i],
                                 poreVolumeNew, poreVolumeOld, dPoreVol_dPres,
                                 accum, accumJacobian );

    // compute etalon
    real64 const accum_et = poroNew[i] * densNew[i] * volume - poroOld[i] * densOld[i] * volume;
    real64 const accumJacobian_et = dPoreVol_dPres * densNew[i] + dDens_dPres[i] * poreVolumeNew;

    EXPECT_DOUBLE_EQ( accum, accum_et );
    EXPECT_DOUBLE_EQ( accumJacobian, accumJacobian_et );
  }
}


int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
