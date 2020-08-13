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
 * @file CompositionalMultiphaseFluid.cpp
 */

#include "CompositionalMultiphaseFluid.hpp"

#include "codingUtilities/Utilities.hpp"

// PVTPackage includes
#include "MultiphaseSystem/CompositionalMultiphaseSystem.hpp"

#include <map>

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
namespace
{
PVTPackage::EOS_TYPE
getCompositionalEosType( string const & name )
{
  static std::map< string, PVTPackage::EOS_TYPE > const eosTypes = {
    { "PR", PVTPackage::EOS_TYPE::PENG_ROBINSON },
    { "SRK", PVTPackage::EOS_TYPE::REDLICH_KWONG_SOAVE } };
  auto const it = eosTypes.find( name );
  GEOSX_ERROR_IF( it == eosTypes.end(),
                  "Compositional EOS type not supported by PVTPackage: " << name );
  return it->second;
}

}  // namespace

CompositionalMultiphaseFluid::CompositionalMultiphaseFluid( std::string const & name,
                                                            Group * const parent ) :
  MultiFluidPVTPackageWrapper( name, parent )
{
  getWrapperBase( viewKeyStruct::componentNamesString )
    ->setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::componentMolarWeightString )
    ->setInputFlag( InputFlags::REQUIRED );
  getWrapperBase( viewKeyStruct::phaseNamesString )->setInputFlag( InputFlags::REQUIRED );

  registerWrapper( viewKeyStruct::equationsOfStateString, &m_equationsOfState )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "List of equation of state types for each phase" );

  registerWrapper( viewKeyStruct::componentCriticalPressureString,
                   &m_componentCriticalPressure )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Component critical pressures" );

  registerWrapper( viewKeyStruct::componentCriticalTemperatureString,
                   &m_componentCriticalTemperature )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Component critical temperatures" );

  registerWrapper( viewKeyStruct::componentAcentricFactorString,
                   &m_componentAcentricFactor )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Component acentric factors" );

  registerWrapper( viewKeyStruct::componentVolumeShiftString,
                   &m_componentVolumeShift )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Component volume shifts" );

  registerWrapper( viewKeyStruct::componentBinaryCoeffString,
                   &m_componentBinaryCoeff )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Table of binary interaction coefficients" );
}

CompositionalMultiphaseFluid::~CompositionalMultiphaseFluid()
{}

void
CompositionalMultiphaseFluid::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< CompositionalMultiphaseFluid >( name, parent );
  }

  MultiFluidPVTPackageWrapper::DeliverClone( name, parent, clone );
  CompositionalMultiphaseFluid & fluid =
    dynamicCast< CompositionalMultiphaseFluid & >( *clone );

  fluid.m_equationsOfState = m_equationsOfState;
  fluid.m_componentCriticalPressure = m_componentCriticalPressure;
  fluid.m_componentCriticalTemperature = m_componentCriticalTemperature;
  fluid.m_componentAcentricFactor = m_componentAcentricFactor;
  fluid.m_componentVolumeShift = m_componentVolumeShift;
  fluid.m_componentBinaryCoeff = m_componentBinaryCoeff;

  fluid.createFluid();
}

void
CompositionalMultiphaseFluid::PostProcessInput()
{
  MultiFluidPVTPackageWrapper::PostProcessInput();

  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

#define COMPFLUID_CHECK_INPUT_LENGTH( data, expected, attr )                   \
  if( LvArray::integerConversion< localIndex >( ( data ).size() ) !=           \
      LvArray::integerConversion< localIndex >( expected ) )                   \
  {                                                                            \
    GEOSX_ERROR( "CompositionalMultiphaseFluid: invalid number of entries in " \
                 << ( attr ) << " attribute (" << ( data ).size() << "given, " \
                 << ( expected ) << " expected)" );                            \
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_equationsOfState,
                                NP,
                                viewKeyStruct::equationsOfStateString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentCriticalPressure,
                                NC,
                                viewKeyStruct::componentCriticalPressureString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentCriticalTemperature,
                                NC,
                                viewKeyStruct::componentCriticalTemperatureString )
  COMPFLUID_CHECK_INPUT_LENGTH( m_componentAcentricFactor,
                                NC,
                                viewKeyStruct::componentAcentricFactorString )

  if( m_componentVolumeShift.empty() )
  {
    m_componentVolumeShift.resize( NC );
    m_componentVolumeShift.setValues< serialPolicy >( 0.0 );
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_componentVolumeShift,
                                NC,
                                viewKeyStruct::componentVolumeShiftString )
  //if (m_componentBinaryCoeff.empty()) TODO needs reading of 2D arrays
  {
    m_componentBinaryCoeff.resize( NC, NC );
    m_componentBinaryCoeff.setValues< serialPolicy >( 0.0 );
  }

  COMPFLUID_CHECK_INPUT_LENGTH( m_componentBinaryCoeff,
                                NC * NC,
                                viewKeyStruct::componentBinaryCoeffString )

#undef COMPFLUID_CHECK_INPUT_LENGTH
}

void
CompositionalMultiphaseFluid::createFluid()
{
  localIndex const NC = numFluidComponents();
  localIndex const NP = numFluidPhases();

  std::vector< PVTPackage::EOS_TYPE > eos( NP );
  std::transform(
    m_equationsOfState.begin(),
    m_equationsOfState.end(),
    eos.begin(),
    []( string const & name ) { return getCompositionalEosType( name ); } );

  std::vector< PVTPackage::PHASE_TYPE > phases( m_phaseTypes.begin(),
                                                m_phaseTypes.end() );
  std::vector< std::string > const components( m_componentNames.begin(),
                                               m_componentNames.end() );
  std::vector< double > const Pc( m_componentCriticalPressure.begin(),
                                  m_componentCriticalPressure.end() );
  std::vector< double > const Tc( m_componentCriticalTemperature.begin(),
                                  m_componentCriticalTemperature.end() );
  std::vector< double > const Mw( m_componentMolarWeight.begin(),
                                  m_componentMolarWeight.end() );
  std::vector< double > const Omega( m_componentAcentricFactor.begin(),
                                     m_componentAcentricFactor.end() );

  ComponentProperties const compProps( NC, components, Mw, Tc, Pc, Omega );

  m_fluid = std::make_unique< PVTPackage::CompositionalMultiphaseSystem >(
    phases,
    eos,
    PVTPackage::COMPOSITIONAL_FLASH_TYPE::NEGATIVE_OIL_GAS,
    compProps );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompositionalMultiphaseFluid,
                        std::string const &,
                        Group * const )
}  // namespace constitutive

}  // namespace geosx
