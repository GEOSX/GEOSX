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
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::SinglePhaseBaseKernels;

TEST( SinglePhaseBaseKernels, mobility )
{
  int constexpr NTEST = 3;

  real64 const dens[NTEST] = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST] = { 1e-5, 1e-10, 0.0 };
  real64 const visc[NTEST] = { 5.0, 2.0, 1.0 };
  real64 const dVisc_dPres[NTEST] = { 1e-7, 0.0, 0.0 };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    real64 mob;
    real64 dMob_dPres;

    MobilityKernel::Compute( dens[i],
                             dDens_dPres[i],
                             visc[i],
                             dVisc_dPres[i],
                             mob,
                             dMob_dPres );

    // compute etalon
    real64 const mob_et = dens[i] / visc[i];
    real64 const dMob_dPres_et =
      mob_et * ( dDens_dPres[i] / dens[i] - dVisc_dPres[i] / visc[i] );

    EXPECT_DOUBLE_EQ( mob, mob_et );
    EXPECT_DOUBLE_EQ( dMob_dPres, dMob_dPres_et );
  }
}

/**
 * @brief Test the assembly of accumulation contribution
 *
 * @note This only tests uncoupled version.
 * In future, porosity update will be elsewhere and this will be simplified.
 */
TEST( SinglePhaseBaseKernels, accumulation )
{
  int constexpr NTEST = 3;

  real64 const densOld[NTEST] = { 700.0, 1200.0, 1500.0 };
  real64 const densNew[NTEST] = { 800.0, 1000.0, 1500.0 };
  real64 const dDens_dPres[NTEST] = { 1e-5, 1e-10, 0.0 };
  real64 const dVol[NTEST] = { 1e-2, 2e-2, 1e-1 };
  real64 const poroOld[NTEST] = { 6e-2, 1e-1, 2e-1 };
  real64 const poroRef[NTEST] = { 5e-2, 2e-1, 3e-1 };
  real64 const pvMult[NTEST] = { 1.20, 1.10, 1.05 };
  real64 const dPvMult_dPres[NTEST] = { 1e-5, 1e-7, 0.0 };
  real64 const volume = 1.0;

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );

    real64 accum;
    real64 accumJacobian;
    real64 poroNew;

    AccumulationKernel< CellElementSubRegion >::Compute< false >( 0.0,
                                                                  densNew[i],
                                                                  densOld[i],
                                                                  dDens_dPres[i],
                                                                  volume,
                                                                  dVol[i],
                                                                  poroRef[i],
                                                                  poroOld[i],
                                                                  pvMult[i],
                                                                  dPvMult_dPres[i],
                                                                  0.0,
                                                                  0.0,
                                                                  0.0,
                                                                  0.0,
                                                                  poroNew,
                                                                  accum,
                                                                  accumJacobian );

    // compute etalon
    real64 const poroNew_et = poroRef[i] * pvMult[i];
    real64 const accum_et = poroNew_et * densNew[i] * ( volume + dVol[i] ) -
      poroOld[i] * densOld[i] * volume;
    real64 const accumJacobian_et =
      dPvMult_dPres[i] * poroRef[i] * densNew[i] * ( volume + dVol[i] ) +
      poroNew_et * dDens_dPres[i] * ( volume + dVol[i] );

    EXPECT_DOUBLE_EQ( poroNew, poroNew_et );
    EXPECT_DOUBLE_EQ( accum, accum_et );
    EXPECT_DOUBLE_EQ( accumJacobian, accumJacobian_et );
  }
}

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
