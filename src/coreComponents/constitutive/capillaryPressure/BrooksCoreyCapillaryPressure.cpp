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
 * @file BrooksCoreyCapillaryPressure.cpp
 */

#include "BrooksCoreyCapillaryPressure.hpp"

#include <cmath>

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
BrooksCoreyCapillaryPressure::BrooksCoreyCapillaryPressure( std::string const & name,
                                                            Group * const parent ) :
  CapillaryPressureBase( name, parent )
{
  registerWrapper( viewKeyStruct::phaseMinVolumeFractionString,
                   &m_phaseMinVolumeFraction )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Minimum volume fraction value for each phase" );

  registerWrapper( viewKeyStruct::phaseCapPressureExponentInvString,
                   &m_phaseCapPressureExponentInv )
    ->setApplyDefaultValue( 2.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Inverse of capillary power law exponent for each phase" );

  registerWrapper( viewKeyStruct::phaseEntryPressureString, &m_phaseEntryPressure )
    ->setApplyDefaultValue( 1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Entry pressure value for each phase" );

  registerWrapper( viewKeyStruct::capPressureEpsilonString, &m_capPressureEpsilon )
    ->setApplyDefaultValue( 1e-6 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Wetting-phase saturation at which the max cap. pressure is attained; "
      "used to avoid infinite cap. pressure values for saturations close to "
      "zero" );
}

BrooksCoreyCapillaryPressure::~BrooksCoreyCapillaryPressure()
{}

void
BrooksCoreyCapillaryPressure::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< BrooksCoreyCapillaryPressure >( name, parent );
  }

  CapillaryPressureBase::DeliverClone( name, parent, clone );
  BrooksCoreyCapillaryPressure & relPerm =
    dynamicCast< BrooksCoreyCapillaryPressure & >( *clone );

  relPerm.m_phaseMinVolumeFraction = m_phaseMinVolumeFraction;
  relPerm.m_phaseCapPressureExponentInv = m_phaseCapPressureExponentInv;
  relPerm.m_phaseEntryPressure = m_phaseEntryPressure;
  relPerm.m_capPressureEpsilon = m_capPressureEpsilon;
  relPerm.m_volFracScale = m_volFracScale;
}

void
BrooksCoreyCapillaryPressure::PostProcessInput()
{
  CapillaryPressureBase::PostProcessInput();

  localIndex const NP = numFluidPhases();

#define COREY_CHECK_INPUT_LENGTH( data, expected, attr )                       \
  if( LvArray::integerConversion< localIndex >( ( data ).size() ) !=           \
      LvArray::integerConversion< localIndex >( expected ) )                   \
  {                                                                            \
    GEOSX_ERROR( "BrooksCoreyCapillaryPressure: invalid number of entries in " \
                 << ( attr ) << " attribute (" << ( data ).size() << "given, " \
                 << ( expected ) << " expected)" );                            \
  }

  COREY_CHECK_INPUT_LENGTH( m_phaseMinVolumeFraction,
                            NP,
                            viewKeyStruct::phaseMinVolumeFractionString )
  COREY_CHECK_INPUT_LENGTH( m_phaseCapPressureExponentInv,
                            NP,
                            viewKeyStruct::phaseCapPressureExponentInvString )
  COREY_CHECK_INPUT_LENGTH( m_phaseEntryPressure,
                            NP,
                            viewKeyStruct::phaseEntryPressureString )

#undef COREY_CHECK_INPUT_LENGTH

  m_volFracScale = 1.0;
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    GEOSX_ERROR_IF(
      ( m_phaseMinVolumeFraction[ip] < 0.0 || m_phaseMinVolumeFraction[ip] > 1.0 ),
      "BrooksCoreyCapillaryPressure: invalid min volume fraction value: "
        << m_phaseMinVolumeFraction[ip] );
    m_volFracScale -= m_phaseMinVolumeFraction[ip];

    GEOSX_ERROR_IF(
      ( m_phaseCapPressureExponentInv[ip] < 1.0 ) &&
        ( m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE ),
      "BrooksCoreyCapillaryPressure: invalid exponent inverse value: "
        << m_phaseCapPressureExponentInv[ip] );

    GEOSX_ERROR_IF( ( m_phaseEntryPressure[ip] < 0.0 ) &&
                      ( m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE ),
                    "BrooksCoreyCapillaryPressure: invalid entry pressure: "
                      << m_phaseEntryPressure[ip] );

    GEOSX_ERROR_IF(
      ( m_capPressureEpsilon < 0.0 || m_capPressureEpsilon > 0.2 ) &&
        ( m_phaseTypes[ip] != CapillaryPressureBase::REFERENCE_PHASE ),
      "BrooksCoreyCapillaryPressure: invalid epsilon: " << m_capPressureEpsilon );
  }

  GEOSX_ERROR_IF(
    m_volFracScale < 0.0,
    "BrooksCoreyCapillaryPressure: sum of min volume fractions exceeds 1.0" );
}

BrooksCoreyCapillaryPressure::KernelWrapper
BrooksCoreyCapillaryPressure::createKernelWrapper()
{
  return KernelWrapper( m_phaseMinVolumeFraction,
                        m_phaseCapPressureExponentInv,
                        m_phaseEntryPressure,
                        m_capPressureEpsilon,
                        m_volFracScale,
                        m_phaseTypes,
                        m_phaseOrder,
                        m_phaseCapPressure,
                        m_dPhaseCapPressure_dPhaseVolFrac );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        BrooksCoreyCapillaryPressure,
                        std::string const &,
                        Group * const )
}  // namespace constitutive

}  // namespace geosx
