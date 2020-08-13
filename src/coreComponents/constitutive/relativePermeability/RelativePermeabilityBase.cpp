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
 * @file RelativePermeabilityBase.cpp
 */

#include "RelativePermeabilityBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
constexpr integer RelativePermeabilityBase::PhaseType::GAS;
constexpr integer RelativePermeabilityBase::PhaseType::OIL;
constexpr integer RelativePermeabilityBase::PhaseType::WATER;

namespace
{
std::unordered_map< string, integer > const phaseDict = {
  { "gas", RelativePermeabilityBase::PhaseType::GAS },
  { "oil", RelativePermeabilityBase::PhaseType::OIL },
  { "water", RelativePermeabilityBase::PhaseType::WATER } };

}

RelativePermeabilityBase::RelativePermeabilityBase( std::string const & name,
                                                    Group * const parent ) :
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseNamesString, &m_phaseNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "List of fluid phases" );

  registerWrapper( viewKeyStruct::phaseTypesString, &m_phaseTypes )
    ->setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseOrderString, &m_phaseOrder )
    ->setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::phaseRelPermString, &m_phaseRelPerm )
    ->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dPhaseRelPerm_dPhaseVolFractionString,
                   &m_dPhaseRelPerm_dPhaseVolFrac );
}

RelativePermeabilityBase::~RelativePermeabilityBase()
{}

void
RelativePermeabilityBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF(
    NP < 2,
    "RelativePermeabilityBase: number of fluid phases should be at least 2" );

  GEOSX_ERROR_IF(
    NP > PhaseType::MAX_NUM_PHASES,
    "RelativePermeabilityBase: number of fluid phases exceeds the maximum of "
      << PhaseType::MAX_NUM_PHASES );

  m_phaseTypes.resize( NP );
  m_phaseOrder.resize( PhaseType::MAX_NUM_PHASES );
  m_phaseOrder.setValues< serialPolicy >( -1 );

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    auto it = phaseDict.find( m_phaseNames[ip] );
    GEOSX_ERROR_IF(
      it == phaseDict.end(),
      "RelativePermeabilityBase: phase not supported: " << m_phaseNames[ip] );
    integer const phaseIndex = it->second;
    GEOSX_ERROR_IF(
      phaseIndex >= PhaseType::MAX_NUM_PHASES,
      "RelativePermeabilityBase: invalid phase index " << phaseIndex );

    m_phaseTypes[ip] = phaseIndex;
    m_phaseOrder[phaseIndex] = LvArray::integerConversion< integer >( ip );
  }

  // call to correctly set member array tertiary sizes on the 'main' material object
  ResizeFields( 0, 0 );
}

void
RelativePermeabilityBase::ResizeFields( localIndex const size,
                                        localIndex const numPts )
{
  localIndex const NP = numFluidPhases();

  m_phaseRelPerm.resize( size, numPts, NP );
  m_dPhaseRelPerm_dPhaseVolFrac.resize( size, numPts, NP, NP );
}

void
RelativePermeabilityBase::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent,
                                              numConstitutivePointsPerParentIndex );
  ResizeFields( parent->size(), numConstitutivePointsPerParentIndex );
}

void
RelativePermeabilityBase::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  ConstitutiveBase::DeliverClone( name, parent, clone );
  RelativePermeabilityBase & relPerm =
    dynamicCast< RelativePermeabilityBase & >( *clone );

  relPerm.m_phaseNames = m_phaseNames;
  relPerm.m_phaseTypes = m_phaseTypes;
  relPerm.m_phaseOrder = m_phaseOrder;
}

}  // namespace constitutive

}  // namespace geosx
