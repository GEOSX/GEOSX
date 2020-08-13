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

// Source inclues
#include "managers/initialization.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

void
testNumericalDerivatives( CapillaryPressureBase & capPressure,
                          arraySlice1d< real64 const > const & saturation,
                          real64 perturbParameter,
                          real64 relTol )
{
  localIndex const NP = capPressure.numFluidPhases();
  auto const & phases = capPressure.phaseNames();

  // create a clone of the capillary pressure to run updates on
  std::unique_ptr< ConstitutiveBase > capPressureCopyPtr;
  capPressure.DeliverClone( "fluidCopy", nullptr, capPressureCopyPtr );
  CapillaryPressureBase & capPressureCopy =
    *capPressureCopyPtr->group_cast< CapillaryPressureBase * >();

  capPressure.AllocateConstitutiveData( capPressure.getParent(), 1 );
  capPressureCopy.AllocateConstitutiveData( capPressure.getParent(), 1 );

  arraySlice1d< real64 const > phaseCapPressure =
    capPressure.phaseCapPressure()[0][0];
  arraySlice2d< real64 const > dPhaseCapPressure_dSat =
    capPressure.dPhaseCapPressure_dPhaseVolFraction()[0][0];
  arraySlice1d< real64 const > phaseCapPressureCopy =
    capPressureCopy.phaseCapPressure()[0][0];

  // set the fluid state to current
  constitutive::constitutiveUpdatePassThru( capPressure, [&]( auto & castedCapPres ) {
    typename TYPEOFREF( castedCapPres )::KernelWrapper capPresWrapper =
      castedCapPres.createKernelWrapper();
    capPresWrapper.Update( 0, 0, saturation );
  } );

  // update saturation and check derivatives
  auto dPhaseCapPressure_dS = invertLayout( dPhaseCapPressure_dSat, NP, NP );

  array1d< real64 > satNew( NP );
  for( localIndex jp = 0; jp < NP; ++jp )
  {
    real64 const dS = perturbParameter * ( saturation[jp] + perturbParameter );
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      satNew[ip] = saturation[ip];
    }
    satNew[jp] += dS;

    constitutive::constitutiveUpdatePassThru(
      capPressureCopy,
      [&]( auto & castedCapPres ) {
        typename TYPEOFREF( castedCapPres )::KernelWrapper capPresWrapper =
          castedCapPres.createKernelWrapper();
        capPresWrapper.Update( 0, 0, satNew );
      } );

    string var = "phaseVolFrac[" + phases[jp] + "]";
    checkDerivative( phaseCapPressureCopy.toSliceConst(),
                     phaseCapPressure.toSliceConst(),
                     dPhaseCapPressure_dS[jp].toSliceConst(),
                     dS,
                     relTol,
                     "phaseCapPressure",
                     var,
                     phases );
  }
}

CapillaryPressureBase *
makeBrooksCoreyCapPressureTwoPhase( string const & name,
                                    Group * parent )
{
  auto capPressure = parent->RegisterGroup< BrooksCoreyCapillaryPressure >( name );

  auto & phaseNames = capPressure->getReference< string_array >(
    CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "water";
  phaseNames[1] = "oil";

  auto & phaseMinSat = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.1;
  phaseMinSat[1] = 0.05;

  auto & phaseCapPressureExpInv = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 2;
  phaseCapPressureExpInv[1] = 4;

  auto & phaseEntryPressure = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString );
  phaseEntryPressure.resize( 2 );
  phaseEntryPressure[0] = 1;
  phaseEntryPressure[1] = 1;

  auto & capPressureEpsilon = capPressure->getReference< real64 >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4;

  capPressure->PostProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase *
makeBrooksCoreyCapPressureThreePhase( string const & name,
                                      Group * parent )
{
  auto capPressure = parent->RegisterGroup< BrooksCoreyCapillaryPressure >( name );

  auto & phaseNames = capPressure->getReference< string_array >(
    CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "water";
  phaseNames[1] = "oil";
  phaseNames[2] = "gas";

  auto & phaseMinSat = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04;
  phaseMinSat[1] = 0.02;
  phaseMinSat[2] = 0.1;

  auto & phaseCapPressureExpInv = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 2;
  phaseCapPressureExpInv[1] = -3;
  phaseCapPressureExpInv[2] = 2.5;

  auto & phaseEntryPressure = capPressure->getReference< array1d< real64 > >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::phaseEntryPressureString );
  phaseEntryPressure.resize( 3 );
  phaseEntryPressure[0] = 1;
  phaseEntryPressure[1] = -1;
  phaseEntryPressure[2] = 2;

  auto & capPressureEpsilon = capPressure->getReference< real64 >(
    BrooksCoreyCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-7;

  capPressure->PostProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase *
makeVanGenuchtenCapPressureTwoPhase( string const & name,
                                     Group * parent )
{
  auto capPressure = parent->RegisterGroup< VanGenuchtenCapillaryPressure >( name );

  auto & phaseNames = capPressure->getReference< string_array >(
    CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil";
  phaseNames[1] = "gas";

  auto & phaseMinSat = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 2 );
  phaseMinSat[0] = 0.04;
  phaseMinSat[1] = 0.1;

  auto & phaseCapPressureExpInv = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 2 );
  phaseCapPressureExpInv[0] = 0.4;
  phaseCapPressureExpInv[1] = 0.5;

  auto & phaseCapPressureMultiplier = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString );
  phaseCapPressureMultiplier.resize( 2 );
  phaseCapPressureMultiplier[0] = 0.5;
  phaseCapPressureMultiplier[1] = 1;

  auto & capPressureEpsilon = capPressure->getReference< real64 >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4;

  capPressure->PostProcessInputRecursive();
  return capPressure;
}

CapillaryPressureBase *
makeVanGenuchtenCapPressureThreePhase( string const & name,
                                       Group * parent )
{
  auto capPressure = parent->RegisterGroup< VanGenuchtenCapillaryPressure >( name );

  auto & phaseNames = capPressure->getReference< string_array >(
    CapillaryPressureBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil";
  phaseNames[1] = "gas";
  phaseNames[2] = "water";

  auto & phaseMinSat = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseMinVolumeFractionString );
  phaseMinSat.resize( 3 );
  phaseMinSat[0] = 0.04;
  phaseMinSat[1] = 0.1;
  phaseMinSat[2] = 0.;

  auto & phaseCapPressureExpInv = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureExponentInvString );
  phaseCapPressureExpInv.resize( 3 );
  phaseCapPressureExpInv[0] = 0.33;
  phaseCapPressureExpInv[1] = 0.4;
  phaseCapPressureExpInv[2] = 0.5;

  auto & phaseCapPressureMultiplier = capPressure->getReference< array1d< real64 > >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::phaseCapPressureMultiplierString );
  phaseCapPressureMultiplier.resize( 3 );
  phaseCapPressureMultiplier[0] = 0.5;
  phaseCapPressureMultiplier[1] = 1;
  phaseCapPressureMultiplier[2] = 0.2;

  auto & capPressureEpsilon = capPressure->getReference< real64 >(
    VanGenuchtenCapillaryPressure::viewKeyStruct::capPressureEpsilonString );
  capPressureEpsilon = 1e-4;

  capPressure->PostProcessInputRecursive();
  return capPressure;
}

TEST( testCapPressure, numericalDerivatives_brooksCoreyCapPressureTwoPhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid =
    makeBrooksCoreyCapPressureTwoPhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = 1.0 - sat[0];
  while( sat[0] <= end_sat )
  {
    testNumericalDerivatives( *fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}

TEST( testCapPressure, numericalDerivatives_brooksCoreyCapPressureThreePhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid =
    makeBrooksCoreyCapPressureThreePhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5 * ( 1 - sat[0] );
  sat[2] = 1.0 - sat[0] - sat[1];
  while( sat[0] <= end_sat )
  {
    testNumericalDerivatives( *fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1 - sat[0] );
    sat[2] = 1.0 - sat[0] - sat[1];
  }
}

TEST( testCapPressure, numericalDerivatives_vanGenuchtenCapPressureTwoPhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid =
    makeVanGenuchtenCapPressureTwoPhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 2 );
  sat[0] = start_sat;
  sat[1] = 1 - sat[1];
  while( sat[0] <= end_sat )
  {
    testNumericalDerivatives( *fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 1 - sat[0];
  }
}

TEST( testCapPressure, numericalDerivatives_vanGenuchtenCapPressureThreePhase )
{
  auto parent = std::make_unique< Group >( "parent", nullptr );
  parent->resize( 1 );

  CapillaryPressureBase * fluid =
    makeVanGenuchtenCapPressureThreePhase( "capPressure", parent.get() );

  parent->Initialize( parent.get() );
  parent->InitializePostInitialConditions( parent.get() );

  real64 const eps = std::sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const tol = 1e-4;

  real64 const start_sat = 0.4;
  real64 const end_sat = 0.6;
  real64 const dS = 1e-1;
  array1d< real64 > sat( 3 );
  sat[0] = start_sat;
  sat[1] = 0.5 * ( 1 - sat[0] );
  sat[2] = 1.0 - sat[0] - sat[1];
  while( sat[0] <= end_sat )
  {
    testNumericalDerivatives( *fluid, sat, eps, tol );
    sat[0] += dS;
    sat[1] = 0.5 * ( 1 - sat[0] );
    sat[2] = 1 - sat[0] - sat[1];
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
