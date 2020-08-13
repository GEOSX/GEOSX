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
 * @file VanGenuchtenBakerRelativePermeability.cpp
 */

#include "VanGenuchtenBakerRelativePermeability.hpp"

#include <cmath>

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
VanGenuchtenBakerRelativePermeability::VanGenuchtenBakerRelativePermeability(
  std::string const & name,
  Group * const parent ) :
  RelativePermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString,
                   &m_phaseMinVolumeFraction )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::waterOilRelPermExponentInvString,
                   &m_waterOilRelPermExponentInv )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Rel perm power law exponent inverse for the pair (water phase, oil "
      "phase) at residual gas saturation" );

  registerWrapper( viewKeyStruct::waterOilRelPermMaxValueString,
                   &m_waterOilRelPermMaxValue )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Maximum rel perm value for the pair (water phase, oil phase) at "
      "residual gas saturation" );

  registerWrapper( viewKeyStruct::gasOilRelPermExponentInvString,
                   &m_gasOilRelPermExponentInv )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Rel perm power law exponent inverse for the pair (gas phase, oil phase) "
      "at residual water saturation" );

  registerWrapper( viewKeyStruct::gasOilRelPermMaxValueString,
                   &m_gasOilRelPermMaxValue )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Maximum rel perm value for the pair (gas phase, oil phase) at residual "
      "water saturation" );
}

VanGenuchtenBakerRelativePermeability::~VanGenuchtenBakerRelativePermeability()
{}

void
VanGenuchtenBakerRelativePermeability::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< VanGenuchtenBakerRelativePermeability >( name, parent );
  }

  RelativePermeabilityBase::DeliverClone( name, parent, clone );
  VanGenuchtenBakerRelativePermeability & relPerm =
    dynamicCast< VanGenuchtenBakerRelativePermeability & >( *clone );

  relPerm.m_phaseMinVolumeFraction = m_phaseMinVolumeFraction;
  relPerm.m_waterOilRelPermExponentInv = m_waterOilRelPermExponentInv;
  relPerm.m_waterOilRelPermMaxValue = m_waterOilRelPermMaxValue;
  relPerm.m_gasOilRelPermExponentInv = m_gasOilRelPermExponentInv;
  relPerm.m_gasOilRelPermMaxValue = m_gasOilRelPermMaxValue;
  relPerm.m_volFracScale = m_volFracScale;
}

void
VanGenuchtenBakerRelativePermeability::PostProcessInput()
{
  RelativePermeabilityBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

  GEOSX_ERROR_IF( m_phaseOrder[PhaseType::OIL] < 0,
                  "VanGenuchtenBakerRelativePermeability: reference oil phase "
                  "has not been defined and must be included in model" );

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr )                             \
  if( LvArray::integerConversion< localIndex >( ( data ).size() ) !=                 \
      LvArray::integerConversion< localIndex >( expected ) )                         \
  {                                                                                  \
    GEOSX_ERROR(                                                                     \
      "VanGenuchtenBakerRelativePermeability: invalid number of entries in "         \
      << ( attr ) << " attribute (" << ( data ).size() << " given, " << ( expected ) \
      << " expected)" );                                                             \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction,
                            NP,
                            viewKeyStruct::phaseMinVolumeFractionString )

  if( m_phaseOrder[PhaseType::WATER] >= 0 )
  {
    COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermExponentInv,
                              2,
                              viewKeyStruct::waterOilRelPermExponentInvString )
    COREY_CHECK_INPUT_LENGTH( m_waterOilRelPermMaxValue,
                              2,
                              viewKeyStruct::waterOilRelPermMaxValueString )
  }

  if( m_phaseOrder[PhaseType::GAS] >= 0 )
  {
    COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermExponentInv,
                              2,
                              viewKeyStruct::gasOilRelPermExponentInvString )
    COREY_CHECK_INPUT_LENGTH( m_gasOilRelPermMaxValue,
                              2,
                              viewKeyStruct::gasOilRelPermMaxValueString )
  }

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    GEOSX_ERROR_IF(
      m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0,
      "VanGenuchtenBakerRelativePermeability: invalid phase min volume "
      "fraction value: "
        << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];
  }
  GEOSX_ERROR_IF( m_volFracScale < 0.0,
                  "VanGenuchtenBakerRelativePermeability: sum of min volume "
                  "fractions exceeds 1.0" );

  for( localIndex ip = 0; ip < 2; ++ip )
  {
    if( m_phaseOrder[PhaseType::WATER] >= 0 )
    {
      GEOSX_ERROR_IF( m_waterOilRelPermExponentInv[ip] < 0.0,
                      "VanGenuchtenBakerRelativePermeability: invalid water-oil "
                      "exponent inverse value: "
                        << m_waterOilRelPermExponentInv[ip] );
      GEOSX_ERROR_IF(
        m_waterOilRelPermMaxValue[ip] < 0.0 || m_waterOilRelPermMaxValue[ip] > 1.0,
        "VanGenuchtenBakerRelativePermeability: invalid maximum value: "
          << m_waterOilRelPermMaxValue[ip] );
    }

    if( m_phaseOrder[PhaseType::GAS] >= 0 )
    {
      GEOSX_ERROR_IF( m_gasOilRelPermExponentInv[ip] < 0.0,
                      "VanGenuchtenBakerRelativePermeability: invalid gas-oil "
                      "exponent inverse value: "
                        << m_gasOilRelPermExponentInv[ip] );
      GEOSX_ERROR_IF(
        m_gasOilRelPermMaxValue[ip] < 0.0 || m_gasOilRelPermMaxValue[ip] > 1.0,
        "VanGenuchtenBakerRelativePermeability: invalid maximum value: "
          << m_gasOilRelPermMaxValue[ip] );
    }
  }

  if( m_phaseOrder[PhaseType::WATER] >= 0 && m_phaseOrder[PhaseType::GAS] >= 0 )
  {
    real64 const mean = 0.5 *
      ( m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL] +
        m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] );
    m_gasOilRelPermMaxValue[GasOilPairPhaseType::OIL] = mean;
    m_waterOilRelPermMaxValue[WaterOilPairPhaseType::OIL] = mean;
  }
}

VanGenuchtenBakerRelativePermeability::KernelWrapper
VanGenuchtenBakerRelativePermeability::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_waterOilRelPermExponentInv,
                        m_waterOilRelPermMaxValue,
                        m_gasOilRelPermExponentInv,
                        m_gasOilRelPermMaxValue,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseRelPerm,
                        m_dPhaseRelPerm_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        VanGenuchtenBakerRelativePermeability,
                        std::string const &,
                        Group * const )
}  // namespace constitutive

}  // namespace geosx
