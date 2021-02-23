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

/**
 * @file MultiFluidBase.cpp
 */

#include "MultiFluidBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

MultiFluidBase::MultiFluidBase( string const & name, Group * const parent )
  : ConstitutiveBase( name, parent ),
  m_useMass( false )
{
  // We make base inputs optional here, since derived classes may want to predefine/hardcode
  // components/phases. Models that do need these inputs should change input flags accordingly.

  registerWrapper( viewKeyStruct::componentNamesString(), &m_componentNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of component names" );

  registerWrapper( viewKeyStruct::componentMolarWeightString(), &m_componentMolarWeight ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Component molar weights" );

  registerWrapper( viewKeyStruct::phaseNamesString(), &m_phaseNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseFractionString(), &m_phaseFraction ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseFraction_dPressureString(), &m_dPhaseFraction_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseFraction_dTemperatureString(), &m_dPhaseFraction_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseFraction_dGlobalCompFractionString(), &m_dPhaseFraction_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseDensityString(), &m_phaseDensity ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseDensity_dPressureString(), &m_dPhaseDensity_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseDensity_dTemperatureString(), &m_dPhaseDensity_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseDensity_dGlobalCompFractionString(), &m_dPhaseDensity_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseMassDensityString(), &m_phaseMassDensity ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dPressureString(), &m_dPhaseMassDensity_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dTemperatureString(), &m_dPhaseMassDensity_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseMassDensity_dGlobalCompFractionString(), &m_dPhaseMassDensity_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseViscosityString(), &m_phaseViscosity ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dPressureString(), &m_dPhaseViscosity_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dTemperatureString(), &m_dPhaseViscosity_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseViscosity_dGlobalCompFractionString(), &m_dPhaseViscosity_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseCompFractionString(), &m_phaseCompFraction );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dPressureString(), &m_dPhaseCompFraction_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dTemperatureString(), &m_dPhaseCompFraction_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseCompFraction_dGlobalCompFractionString(), &m_dPhaseCompFraction_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::totalDensityString(), &m_totalDensity ).setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dTotalDensity_dPressureString(), &m_dTotalDensity_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dTotalDensity_dTemperatureString(), &m_dTotalDensity_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dTotalDensity_dGlobalCompFractionString(), &m_dTotalDensity_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::useMassString(), &m_useMass ).
    setRestartFlags( RestartFlags::NO_WRITE );

  registerWrapper( viewKeyStruct::phaseEnthalpyString(), &m_phaseEnthalpy ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseEnthalpy_dPressureString(), &m_dPhaseEnthalpy_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseEnthalpy_dTemperatureString(), &m_dPhaseEnthalpy_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseEnthalpy_dGlobalCompFractionString(), &m_dPhaseEnthalpy_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );


  registerWrapper( viewKeyStruct::phaseInternalEnergyString(), &m_phaseInternalEnergy ).
    setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseInternalEnergy_dPressureString(), &m_dPhaseInternalEnergy_dPressure ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseInternalEnergy_dTemperatureString(), &m_dPhaseInternalEnergy_dTemperature ).
    setRestartFlags( RestartFlags::NO_WRITE );
  registerWrapper( viewKeyStruct::dPhaseInternalEnergy_dGlobalCompFractionString(), &m_dPhaseInternalEnergy_dGlobalCompFraction ).
    setRestartFlags( RestartFlags::NO_WRITE );
  
}

void MultiFluidBase::resizeFields( localIndex const size, localIndex const numPts )
{
  localIndex const NP = numFluidPhases();
  localIndex const NC = numFluidComponents();

  m_phaseFraction.resize( size, numPts, NP );
  m_dPhaseFraction_dPressure.resize( size, numPts, NP );
  m_dPhaseFraction_dTemperature.resize( size, numPts, NP );
  m_dPhaseFraction_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseDensity.resize( size, numPts, NP );
  m_dPhaseDensity_dPressure.resize( size, numPts, NP );
  m_dPhaseDensity_dTemperature.resize( size, numPts, NP );
  m_dPhaseDensity_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseMassDensity.resize( size, numPts, NP );
  m_dPhaseMassDensity_dPressure.resize( size, numPts, NP );
  m_dPhaseMassDensity_dTemperature.resize( size, numPts, NP );
  m_dPhaseMassDensity_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseViscosity.resize( size, numPts, NP );
  m_dPhaseViscosity_dPressure.resize( size, numPts, NP );
  m_dPhaseViscosity_dTemperature.resize( size, numPts, NP );
  m_dPhaseViscosity_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseCompFraction.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dPressure.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dTemperature.resize( size, numPts, NP, NC );
  m_dPhaseCompFraction_dGlobalCompFraction.resize( size, numPts, NP, NC, NC );

  m_totalDensity.resize( size, numPts );
  m_dTotalDensity_dPressure.resize( size, numPts );
  m_dTotalDensity_dTemperature.resize( size, numPts );
  m_dTotalDensity_dGlobalCompFraction.resize( size, numPts, NC );

  m_phaseEnthalpy.resize( size, numPts, NP );
  m_dPhaseEnthalpy_dPressure.resize( size, numPts, NP );
  m_dPhaseEnthalpy_dTemperature.resize( size, numPts, NP );
  m_dPhaseEnthalpy_dGlobalCompFraction.resize( size, numPts, NP, NC );

  m_phaseInternalEnergy.resize( size, numPts, NP );
  m_dPhaseInternalEnergy_dPressure.resize( size, numPts, NP );
  m_dPhaseInternalEnergy_dTemperature.resize( size, numPts, NP );
  m_dPhaseInternalEnergy_dGlobalCompFraction.resize( size, numPts, NP, NC );
  

}

void MultiFluidBase::allocateConstitutiveData( dataRepository::Group & parent,
                                               localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  resizeFields( parent.size(), numConstitutivePointsPerParentIndex );
}

MultiFluidBase::~MultiFluidBase()
{}


void MultiFluidBase::postProcessInput()
{
  ConstitutiveBase::postProcessInput();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF( NC == 0, "MultiFluidBase: No fluid components specified" );

  GEOSX_ERROR_IF( NC > MAX_NUM_COMPONENTS,
                  "MultiFluidBase: Number of fluid components exceeds the maximum of " << MAX_NUM_COMPONENTS );

  GEOSX_ERROR_IF( NP == 0, "MultiFluidBase: No fluid phases specified" );

  GEOSX_ERROR_IF( NP > MAX_NUM_PHASES,
                  "MultiFluidBase: Number of fluid phases exceeds the maximum of " << MAX_NUM_PHASES );

  #define MULTIFLUID_CHECK_INPUT_LENGTH( data, expected, attr ) \
    if( LvArray::integerConversion< localIndex >((data).size()) != LvArray::integerConversion< localIndex >( expected )) \
    { \
      GEOSX_ERROR( "MultiFluidBase: invalid number of entries in " \
                   << (attr) << " attribute (" \
                   << (data).size() << "given, " \
                   << (expected) << " expected)" ); \
    }

  MULTIFLUID_CHECK_INPUT_LENGTH( m_componentMolarWeight, NC, viewKeyStruct::componentMolarWeightString() )

  // call to correctly set member array tertiary sizes on the 'main' material object
  resizeFields( 0, 0 );
}

bool MultiFluidBase::getMassFlag() const
{
  return m_useMass;
}

void MultiFluidBase::setMassFlag( bool const flag )
{
  m_useMass = flag;
}

} //namespace constitutive

} //namespace geosx
