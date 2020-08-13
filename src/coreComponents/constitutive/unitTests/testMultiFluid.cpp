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
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/fluid/MultiFluidUtils.hpp"
#include "physicsSolvers/fluidFlow/unitTests/testCompFlowUtils.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::testing;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

/// Black-oil tables written into temporary files during testing

static const char * pvtg_str =
  "#\tPg(Pa)\t\tRv(sm3/sm3)\tBg(m3/sm3)\tVisc(Pa.s)\n"
  "\n"
  "\t3000000\t\t0.000132\t0.04234\t    0.00001344\n"
  "\t\t\t\t0\t\t\t0.04231\t    0.00001389\n"
  "\t6000000\t\t0.000124\t0.02046\t    0.0000142\n"
  "\t\t\t\t0\t\t\t0.02043\t    0.0000145\n"
  "\t9000000\t\t0.000126\t0.01328\t    0.00001526\n"
  "\t\t\t\t0\t\t\t0.01325\t    0.00001532\n"
  "   12000000\t\t0.000135\t0.00977\t    0.0000166\n"
  "\t\t\t\t0\t\t\t0.00973\t    0.00001634\n"
  "   15000000\t\t0.000149\t0.00773\t    0.00001818\n"
  "\t\t\t\t0\t\t\t0.00769\t    0.00001752\n"
  "   18000000\t\t0.000163\t0.006426\t0.00001994\n"
  "\t\t\t\t0\t\t\t0.006405\t0.00001883\n"
  "   21000000\t\t0.000191\t0.005541\t0.00002181\n"
  "\t\t\t\t0\t\t\t0.005553\t0.00002021\n"
  "   24000000\t\t0.000225\t0.004919\t0.0000237\n"
  "\t\t\t\t0\t\t\t0.004952\t0.00002163\n"
  "   27000000\t\t0.000272\t0.004471\t0.00002559\n"
  "\t\t\t\t0\t\t\t0.004511\t0.00002305\n"
  "   29500000\t\t0.000354\t0.004194\t0.00002714\n"
  "\t\t\t\t0\t\t\t0.004225\t0.00002423\n"
  "   31000000\t\t0.000403\t0.004031\t0.00002806\n"
  "\t\t\t\t0.000354\t0.004059\t0.00002768\n"
  "   33000000\t\t0.000354\t0.00391\t    0.00002832\n"
  "\t\t\t\t0\t\t\t0.003913\t0.00002583\n"
  "   53000000\t\t0.000479\t0.003868\t0.00002935\n"
  "\t\t\t\t0.000354\t0.0039\t\t0.00002842\n"
  "\t\t\t\t0\t\t\t0.003903\t0.00002593";

static const char * pvto_str =
  "# Rs[sm3/sm3]\tPbub[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
  "\n"
  "  2\t            2000000\t    1.02\t    0.000975\n"
  "  5\t            5000000\t    1.03\t    0.00091\n"
  " 10\t            10000000\t1.04\t    0.00083\n"
  " 15\t            20000000\t1.05\t    0.000695\n"
  "                90000000\t1.03\t    0.000985  -- some line comment\n"
  " 30\t            30000000\t1.07\t    0.000594\n"
  " 40\t            40000000\t1.08\t    0.00051\n"
  "                50000000\t1.07\t    0.000549  -- another one\n"
  "                90000000\t1.06\t    0.00074\n"
  " 50\t            50000000.7\t1.09\t    0.000449\n"
  "                90000000.7\t1.08\t    0.000605";

static const char * pvtw_str =
  "#\tPref[bar]\tBw[m3/sm3]\tCp[1/bar]\t    Visc[cP]\n"
  "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

/// Dead-oil tables written into temporary files during testing

static const char * pvdg_str =
  "# Pg(Pa)\tBg(m3/sm3)\tVisc(Pa.s)\n"
  "\n"
  "3000000\t\t0.04234\t\t0.00001344\n"
  "6000000\t\t0.02046\t\t0.0000142\n"
  "9000000\t\t0.01328\t\t0.00001526\n"
  "12000000\t0.00977\t\t0.0000166\n"
  "15000000\t0.00773\t\t0.00001818\n"
  "18000000\t0.006426\t0.00001994\n"
  "21000000\t0.005541\t0.00002181\n"
  "24000000\t0.004919\t0.0000237\n"
  "27000000\t0.004471\t0.00002559\n"
  "29500000\t0.004194\t0.00002714\n"
  "31000000\t0.004031\t0.00002806\n"
  "33000000\t0.00391\t\t0.00002832\n"
  "53000000\t0.003868\t0.00002935";

static const char * pvdo_str =
  "#P[Pa]\tBo[m3/sm3]\tVisc(Pa.s)\n"
  "\n"
  "2000000\t\t1.02\t0.000975\n"
  "5000000\t\t1.03\t0.00091\n"
  "10000000\t1.04\t0.00083\n"
  "20000000\t1.05\t0.000695\n"
  "30000000\t1.07\t0.000594\n"
  "40000000\t1.08\t0.00051\n"
  "50000000.7\t1.09\t0.000449";

static const char * pvdw_str =
  "#\tPref[bar]\tBw[m3/sm3]\tCp[1/bar]\t    Visc[cP]\n"
  "\t30600000.1\t1.03\t\t0.00000000041\t0.0003";

void
testNumericalDerivatives(
  MultiFluidBase & fluid,
  real64 const P,
  real64 const T,
  arraySlice1d< real64 > const & composition,
  real64 const perturbParameter,
  real64 const relTol,
  real64 const absTol = std::numeric_limits< real64 >::max() )
{
  localIndex const NC = fluid.numFluidComponents();
  localIndex const NP = fluid.numFluidPhases();

  auto const & components = fluid.getReference< string_array >(
    MultiFluidBase::viewKeyStruct::componentNamesString );
  auto const & phases = fluid.getReference< string_array >(
    MultiFluidBase::viewKeyStruct::phaseNamesString );

  // create a clone of the fluid to run updates on
  std::unique_ptr< ConstitutiveBase > fluidCopyPtr;
  fluid.DeliverClone( "fluidCopy", nullptr, fluidCopyPtr );
  MultiFluidBase & fluidCopy = *fluidCopyPtr->group_cast< MultiFluidBase * >();

  fluid.AllocateConstitutiveData( fluid.getParent(), 1 );
  fluidCopy.AllocateConstitutiveData( fluid.getParent(), 1 );

// extract data views from both fluids
#define GET_FLUID_DATA( FLUID, DIM, KEY )     \
  FLUID.getReference< Array< real64, DIM > >( \
    MultiFluidBase::viewKeyStruct::KEY )[0][0]

  CompositionalVarContainer< 1 > phaseFrac {
    GET_FLUID_DATA( fluid, 3, phaseFractionString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseFraction_dGlobalCompFractionString ) };

  CompositionalVarContainer< 1 > phaseDens {
    GET_FLUID_DATA( fluid, 3, phaseDensityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseDensity_dGlobalCompFractionString ) };

  CompositionalVarContainer< 1 > phaseVisc {
    GET_FLUID_DATA( fluid, 3, phaseViscosityString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dPressureString ),
    GET_FLUID_DATA( fluid, 3, dPhaseViscosity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseViscosity_dGlobalCompFractionString ) };

  CompositionalVarContainer< 2 > phaseCompFrac {
    GET_FLUID_DATA( fluid, 4, phaseCompFractionString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dPressureString ),
    GET_FLUID_DATA( fluid, 4, dPhaseCompFraction_dTemperatureString ),
    GET_FLUID_DATA( fluid, 5, dPhaseCompFraction_dGlobalCompFractionString ) };

  CompositionalVarContainer< 0 > totalDens {
    GET_FLUID_DATA( fluid, 2, totalDensityString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dPressureString ),
    GET_FLUID_DATA( fluid, 2, dTotalDensity_dTemperatureString ),
    GET_FLUID_DATA( fluid, 3, dTotalDensity_dGlobalCompFractionString ) };

  auto const & phaseFracCopy = GET_FLUID_DATA( fluidCopy, 3, phaseFractionString );
  auto const & phaseDensCopy = GET_FLUID_DATA( fluidCopy, 3, phaseDensityString );
  auto const & phaseViscCopy = GET_FLUID_DATA( fluidCopy, 3, phaseViscosityString );
  auto const & phaseCompFracCopy =
    GET_FLUID_DATA( fluidCopy, 4, phaseCompFractionString );
  auto const & totalDensCopy = GET_FLUID_DATA( fluidCopy, 2, totalDensityString );

#undef GET_FLUID_DATA

  // set the original fluid state to current
  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid ) {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
      castedFluid.createKernelWrapper();
    fluidWrapper.Update( 0, 0, P, T, composition );
  } );

  // now perturb variables and update the copied fluid's state
  constitutive::constitutiveUpdatePassThru( fluidCopy, [&]( auto & castedFluid ) {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
      castedFluid.createKernelWrapper();

    // update pressure and check derivatives
    {
      real64 const dP = perturbParameter * ( P + perturbParameter );
      fluidWrapper.Update( 0, 0, P + dP, T, composition );

      checkDerivative( phaseFracCopy,
                       phaseFrac.value,
                       phaseFrac.dPres,
                       dP,
                       relTol,
                       absTol,
                       "phaseFrac",
                       "Pres",
                       phases );
      checkDerivative( phaseDensCopy,
                       phaseDens.value,
                       phaseDens.dPres,
                       dP,
                       relTol,
                       absTol,
                       "phaseDens",
                       "Pres",
                       phases );
      checkDerivative( phaseViscCopy,
                       phaseVisc.value,
                       phaseVisc.dPres,
                       dP,
                       relTol,
                       absTol,
                       "phaseVisc",
                       "Pres",
                       phases );
      checkDerivative( totalDensCopy,
                       totalDens.value,
                       totalDens.dPres,
                       dP,
                       relTol,
                       absTol,
                       "totalDens",
                       "Pres" );
      checkDerivative( phaseCompFracCopy.toSliceConst(),
                       phaseCompFrac.value.toSliceConst(),
                       phaseCompFrac.dPres.toSliceConst(),
                       dP,
                       relTol,
                       absTol,
                       "phaseCompFrac",
                       "Pres",
                       phases,
                       components );
    }

    // update temperature and check derivatives
    {
      real64 const dT = perturbParameter * ( T + perturbParameter );
      fluidWrapper.Update( 0, 0, P, T + dT, composition );

      checkDerivative( phaseFracCopy,
                       phaseFrac.value,
                       phaseFrac.dTemp,
                       dT,
                       relTol,
                       absTol,
                       "phaseFrac",
                       "Temp",
                       phases );
      checkDerivative( phaseDensCopy,
                       phaseDens.value,
                       phaseDens.dTemp,
                       dT,
                       relTol,
                       absTol,
                       "phaseDens",
                       "Temp",
                       phases );
      checkDerivative( phaseViscCopy,
                       phaseVisc.value,
                       phaseVisc.dTemp,
                       dT,
                       relTol,
                       absTol,
                       "phaseVisc",
                       "Temp",
                       phases );
      checkDerivative( totalDensCopy,
                       totalDens.value,
                       totalDens.dTemp,
                       dT,
                       relTol,
                       absTol,
                       "totalDens",
                       "Temp" );
      checkDerivative( phaseCompFracCopy.toSliceConst(),
                       phaseCompFrac.value.toSliceConst(),
                       phaseCompFrac.dTemp.toSliceConst(),
                       dT,
                       relTol,
                       absTol,
                       "phaseCompFrac",
                       "Temp",
                       phases,
                       components );
    }

    // update composition and check derivatives
    auto dPhaseFrac_dC = invertLayout( phaseFrac.dComp, NP, NC );
    auto dPhaseDens_dC = invertLayout( phaseDens.dComp, NP, NC );
    auto dPhaseVisc_dC = invertLayout( phaseVisc.dComp, NP, NC );
    auto dTotalDens_dC = invertLayout( totalDens.dComp, NC );
    auto dPhaseCompFrac_dC = invertLayout( phaseCompFrac.dComp, NP, NC, NC );

    array1d< real64 > compNew( NC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      real64 const dC = perturbParameter * ( composition[jc] + perturbParameter );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compNew[ic] = composition[ic];
      }
      compNew[jc] += dC;

      // renormalize
      real64 sum = 0.0;
      for( localIndex ic = 0; ic < NC; ++ic ) sum += compNew[ic];
      for( localIndex ic = 0; ic < NC; ++ic ) compNew[ic] /= sum;

      fluidWrapper.Update( 0, 0, P, T, compNew );

      string const var = "compFrac[" + components[jc] + "]";
      checkDerivative( phaseFracCopy,
                       phaseFrac.value,
                       dPhaseFrac_dC[jc],
                       dC,
                       relTol,
                       absTol,
                       "phaseFrac",
                       var,
                       phases );
      checkDerivative( phaseDensCopy,
                       phaseDens.value,
                       dPhaseDens_dC[jc],
                       dC,
                       relTol,
                       absTol,
                       "phaseDens",
                       var,
                       phases );
      checkDerivative( phaseViscCopy,
                       phaseVisc.value,
                       dPhaseVisc_dC[jc],
                       dC,
                       relTol,
                       absTol,
                       "phaseVisc",
                       var,
                       phases );
      checkDerivative( totalDensCopy,
                       totalDens.value,
                       dTotalDens_dC[jc],
                       dC,
                       relTol,
                       absTol,
                       "totalDens",
                       var );
      checkDerivative( phaseCompFracCopy.toSliceConst(),
                       phaseCompFrac.value.toSliceConst(),
                       dPhaseCompFrac_dC[jc].toSliceConst(),
                       dC,
                       relTol,
                       absTol,
                       "phaseCompFrac",
                       var,
                       phases,
                       components );
    }
  } );
}

MultiFluidBase *
makeCompositionalFluid( string const & name, Group & parent )
{
  auto fluid = parent.RegisterGroup< CompositionalMultiphaseFluid >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::componentNamesString );
  compNames.resize( 4 );
  compNames[0] = "N2";
  compNames[1] = "C10";
  compNames[2] = "C20";
  compNames[3] = "H20";

  auto & molarWgt = fluid->getReference< array1d< real64 > >(
    MultiFluidBase::viewKeyStruct::componentMolarWeightString );
  molarWgt.resize( 4 );
  molarWgt[0] = 28e-3;
  molarWgt[1] = 134e-3;
  molarWgt[2] = 275e-3;
  molarWgt[3] = 18e-3;

  auto & phaseNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 2 );
  phaseNames[0] = "oil";
  phaseNames[1] = "gas";

  auto & eqnOfState = fluid->getReference< string_array >(
    CompositionalMultiphaseFluid::viewKeyStruct::equationsOfStateString );
  eqnOfState.resize( 2 );
  eqnOfState[0] = "PR";
  eqnOfState[1] = "PR";

  auto & critPres = fluid->getReference< array1d< real64 > >(
    CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalPressureString );
  critPres.resize( 4 );
  critPres[0] = 34e5;
  critPres[1] = 25.3e5;
  critPres[2] = 14.6e5;
  critPres[3] = 220.5e5;

  auto & critTemp = fluid->getReference< array1d< real64 > >(
    CompositionalMultiphaseFluid::viewKeyStruct::componentCriticalTemperatureString );
  critTemp.resize( 4 );
  critTemp[0] = 126.2;
  critTemp[1] = 622.0;
  critTemp[2] = 782.0;
  critTemp[3] = 647.0;

  auto & acFactor = fluid->getReference< array1d< real64 > >(
    CompositionalMultiphaseFluid::viewKeyStruct::componentAcentricFactorString );
  acFactor.resize( 4 );
  acFactor[0] = 0.04;
  acFactor[1] = 0.443;
  acFactor[2] = 0.816;
  acFactor[3] = 0.344;

  fluid->PostProcessInputRecursive();
  return fluid;
}

class CompositionalFluidTest : public ::testing::Test
{
protected:
  virtual void
  SetUp() override
  {
    parent = std::make_unique< Group >( "parent", nullptr );
    parent->resize( 1 );

    fluid = makeCompositionalFluid( "fluid", *parent );

    parent->Initialize( parent.get() );
    parent->InitializePostInitialConditions( parent.get() );
  }

  std::unique_ptr< Group > parent;
  MultiFluidBase * fluid;
};

TEST_F( CompositionalFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099;
  comp[1] = 0.3;
  comp[2] = 0.6;
  comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol );
}

TEST_F( CompositionalFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 4 );
  comp[0] = 0.099;
  comp[1] = 0.3;
  comp[2] = 0.6;
  comp[3] = 0.001;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-2;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol );
}

MultiFluidBase *
makeLiveOilFluid( string const & name, Group * parent )
{
  auto fluid = parent->RegisterGroup< BlackOilFluid >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::componentNamesString );
  compNames.resize( 3 );
  compNames[0] = "oil";
  compNames[1] = "gas";
  compNames[2] = "water";

  auto & molarWgt = fluid->getReference< array1d< real64 > >(
    MultiFluidBase::viewKeyStruct::componentMolarWeightString );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3;
  molarWgt[1] = 16e-3;
  molarWgt[2] = 18e-3;

  auto & phaseNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil";
  phaseNames[1] = "gas";
  phaseNames[2] = "water";

  auto & surfaceDens = fluid->getReference< array1d< real64 > >(
    BlackOilFluid::viewKeyStruct::surfaceDensitiesString );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0;
  surfaceDens[1] = 0.9907;
  surfaceDens[2] = 1022.0;

  auto & tableNames = fluid->getReference< path_array >(
    BlackOilFluid::viewKeyStruct::tableFilesString );
  tableNames.resize( 3 );
  tableNames[0] = "pvto.txt";
  tableNames[1] = "pvtg.txt";
  tableNames[2] = "pvtw.txt";

  auto & fluidType =
    fluid->getReference< string >( BlackOilFluid::viewKeyStruct::fluidTypeString );
  fluidType = "LiveOil";

  fluid->PostProcessInputRecursive();
  return fluid;
}

MultiFluidBase *
makeDeadOilFluid( string const & name, Group * parent )
{
  auto fluid = parent->RegisterGroup< BlackOilFluid >( name );

  // TODO we should actually create a fake XML node with data, but this seemed easier...

  auto & compNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::componentNamesString );
  compNames.resize( 3 );
  compNames[0] = "oil";
  compNames[1] = "gas";
  compNames[2] = "water";

  auto & molarWgt = fluid->getReference< array1d< real64 > >(
    MultiFluidBase::viewKeyStruct::componentMolarWeightString );
  molarWgt.resize( 3 );
  molarWgt[0] = 114e-3;
  molarWgt[1] = 16e-3;
  molarWgt[2] = 18e-3;

  auto & phaseNames = fluid->getReference< string_array >(
    MultiFluidBase::viewKeyStruct::phaseNamesString );
  phaseNames.resize( 3 );
  phaseNames[0] = "oil";
  phaseNames[1] = "gas";
  phaseNames[2] = "water";

  auto & surfaceDens = fluid->getReference< array1d< real64 > >(
    BlackOilFluid::viewKeyStruct::surfaceDensitiesString );
  surfaceDens.resize( 3 );
  surfaceDens[0] = 800.0;
  surfaceDens[1] = 0.9907;
  surfaceDens[2] = 1022.0;

  auto & tableNames = fluid->getReference< path_array >(
    BlackOilFluid::viewKeyStruct::tableFilesString );
  tableNames.resize( 3 );
  tableNames[0] = "pvdo.txt";
  tableNames[1] = "pvdg.txt";
  tableNames[2] = "pvdw.txt";

  auto & fluidType =
    fluid->getReference< string >( BlackOilFluid::viewKeyStruct::fluidTypeString );
  fluidType = "DeadOil";

  fluid->PostProcessInputRecursive();
  return fluid;
}

void
writeTableToFile( std::string const & filename, char const * str )
{
  std::ofstream os( filename );
  ASSERT_TRUE( os.is_open() );
  os << str;
  os.close();
}

void
removeFile( std::string const & filename )
{
  int const ret = std::remove( filename.c_str() );
  ASSERT_TRUE( ret == 0 );
}

class LiveOilFluidTest : public ::testing::Test
{
protected:
  virtual void
  SetUp() override
  {
    writeTableToFile( "pvto.txt", pvto_str );
    writeTableToFile( "pvtg.txt", pvtg_str );
    writeTableToFile( "pvtw.txt", pvtw_str );

    parent = std::make_unique< Group >( "parent", nullptr );
    parent->resize( 1 );
    fluid = makeLiveOilFluid( "fluid", parent.get() );

    parent->Initialize( parent.get() );
    parent->InitializePostInitialConditions( parent.get() );
  }

  virtual void
  TearDown() override
  {
    removeFile( "pvto.txt" );
    removeFile( "pvtg.txt" );
    removeFile( "pvtw.txt" );
  }

  std::unique_ptr< Group > parent;
  MultiFluidBase * fluid;
};

TEST_F( LiveOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1;
  comp[1] = 0.3;
  comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol );
}

TEST_F( LiveOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1;
  comp[1] = 0.3;
  comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-2;
  real64 const absTol = 1e-14;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol, absTol );
}

class DeadOilFluidTest : public ::testing::Test
{
protected:
  virtual void
  SetUp() override
  {
    writeTableToFile( "pvdo.txt", pvdo_str );
    writeTableToFile( "pvdg.txt", pvdg_str );
    writeTableToFile( "pvdw.txt", pvdw_str );

    parent = std::make_unique< Group >( "parent", nullptr );
    parent->resize( 1 );
    fluid = makeDeadOilFluid( "fluid", parent.get() );

    parent->Initialize( parent.get() );
    parent->InitializePostInitialConditions( parent.get() );
  }

  virtual void
  TearDown() override
  {
    removeFile( "pvdo.txt" );
    removeFile( "pvdg.txt" );
    removeFile( "pvdw.txt" );
  }

  std::unique_ptr< Group > parent;
  MultiFluidBase * fluid;
};

TEST_F( DeadOilFluidTest, numericalDerivativesMolar )
{
  fluid->setMassFlag( false );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1;
  comp[1] = 0.3;
  comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-4;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol );
}

TEST_F( DeadOilFluidTest, numericalDerivativesMass )
{
  fluid->setMassFlag( true );

  // TODO test over a range of values
  real64 const P = 5e6;
  real64 const T = 297.15;
  array1d< real64 > comp( 3 );
  comp[0] = 0.1;
  comp[1] = 0.3;
  comp[2] = 0.6;

  real64 const eps = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const relTol = 1e-2;
  real64 const absTol = 1e-14;

  testNumericalDerivatives( *fluid, P, T, comp, eps, relTol, absTol );
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
