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

/**
 * @file CapillaryPressureBase.cpp
 */

#include "CapillaryPressureBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
constexpr integer CapillaryPressureBase::PhaseType::GAS;
constexpr integer CapillaryPressureBase::PhaseType::OIL;
constexpr integer CapillaryPressureBase::PhaseType::WATER;

namespace
{
std::unordered_map< string, integer > const phaseDict = {
  { "gas", CapillaryPressureBase::PhaseType::GAS },
  { "oil", CapillaryPressureBase::PhaseType::OIL },
  { "water", CapillaryPressureBase::PhaseType::WATER } };

}

CapillaryPressureBase::CapillaryPressureBase( std::string const & name,
                                              Group * const parent ) :
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseTypesString, &m_phaseTypes )
    ->setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseOrderString, &m_phaseOrder )
    ->setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseCapPressureString, &m_phaseCapPressure )
    ->setPlotLevel( PlotLevel::LEVEL_0 );

  registerWrapper( viewKeyStruct::dPhaseCapPressure_dPhaseVolFractionString,
                   &m_dPhaseCapPressure_dPhaseVolFrac );
}

CapillaryPressureBase::~CapillaryPressureBase()
{}

void
CapillaryPressureBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF(
    NP < 2,
    "CapillaryPressureBase: number of fluid phases should be at least 2" );

  GEOSX_ERROR_IF(
    NP > PhaseType::MAX_NUM_PHASES,
    "CapillaryPressureBase: number of fluid phases exceeds the maximum of "
      << PhaseType::MAX_NUM_PHASES );

  m_phaseTypes.resize( NP );
  m_phaseOrder.resize( PhaseType::MAX_NUM_PHASES );
  m_phaseOrder.setValues< serialPolicy >( -1 );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    auto it = phaseDict.find( m_phaseNames[ip] );
    GEOSX_ERROR_IF(
      it == phaseDict.end(),
      "CapillaryPressureBase: phase not supported: " << m_phaseNames[ip] );
    integer const phaseIndex = it->second;
    GEOSX_ERROR_IF( phaseIndex >= PhaseType::MAX_NUM_PHASES,
                    "CapillaryPressureBase: invalid phase index " << phaseIndex );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = LvArray::integerConversion< integer >( ip );
  }

  GEOSX_ERROR_IF( m_phaseOrder[CapillaryPressureBase::REFERENCE_PHASE] < 0,
                  "CapillaryPressureBase: reference oil phase has not been "
                  "defined and should be included in model" );

  // call to correctly set member array tertiary sizes on the 'main' material object
  ResizeFields( 0, 0 );
}

void
CapillaryPressureBase::ResizeFields( localIndex const size,
                                     localIndex const numPts )
{
  localIndex const NP = numFluidPhases();

  m_phaseCapPressure.resize( size, numPts, NP );
  m_dPhaseCapPressure_dPhaseVolFrac.resize( size, numPts, NP, NP );
}

void
CapillaryPressureBase::DeliverClone( string const & name,
                                     dataRepository::Group * const parent,
                                     std::unique_ptr< ConstitutiveBase > & clone ) const
{
  ConstitutiveBase::DeliverClone( name, parent, clone );
  CapillaryPressureBase & relPerm = dynamicCast< CapillaryPressureBase & >( *clone );

  relPerm.m_phaseNames = m_phaseNames;
  relPerm.m_phaseTypes = m_phaseTypes;
  relPerm.m_phaseOrder = m_phaseOrder;
}

void
CapillaryPressureBase::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent,
                                              numConstitutivePointsPerParentIndex );
  ResizeFields( parent->size(), numConstitutivePointsPerParentIndex );
}

}  // namespace constitutive

}  // namespace geosx
