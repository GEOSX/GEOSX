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
 * @file CompressibleSinglePhaseFluid.cpp
 */

#include "CompressibleSinglePhaseFluid.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
CompressibleSinglePhaseFluid::CompressibleSinglePhaseFluid( std::string const & name,
                                                            Group * const parent ) :
  SingleFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Fluid compressibility" );

  registerWrapper( viewKeyStruct::viscosibilityString, &m_viscosibility )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Fluid viscosity exponential coefficient" );

  registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Reference pressure" );

  registerWrapper( viewKeyStruct::referenceDensityString, &m_referenceDensity )
    ->setApplyDefaultValue( 1000.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Reference fluid density" );

  registerWrapper( viewKeyStruct::referenceViscosityString, &m_referenceViscosity )
    ->setApplyDefaultValue( 0.001 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Reference fluid viscosity" );

  registerWrapper( viewKeyStruct::densityModelString, &m_densityModelString )
    ->setApplyDefaultValue( "linear" )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Type of density model (linear, quadratic, exponential)" );

  registerWrapper( viewKeyStruct::viscosityModelString, &m_viscosityModelString )
    ->setApplyDefaultValue( "linear" )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Type of viscosity model (linear, quadratic, exponential)" );
}

CompressibleSinglePhaseFluid::~CompressibleSinglePhaseFluid() = default;

void
CompressibleSinglePhaseFluid::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  SingleFluidBase::AllocateConstitutiveData( parent,
                                             numConstitutivePointsPerParentIndex );

  m_density.setValues< serialPolicy >( m_referenceDensity );
  m_viscosity.setValues< serialPolicy >( m_referenceViscosity );
}

void
CompressibleSinglePhaseFluid::DeliverClone(
  string const & name,
  Group * const parent,
  std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< CompressibleSinglePhaseFluid >( name, parent );
  }
  SingleFluidBase::DeliverClone( name, parent, clone );
  CompressibleSinglePhaseFluid & fluid =
    dynamicCast< CompressibleSinglePhaseFluid & >( *clone );

  fluid.m_compressibility = m_compressibility;
  fluid.m_viscosibility = m_viscosibility;
  fluid.m_referencePressure = m_referencePressure;
  fluid.m_referenceDensity = m_referenceDensity;
  fluid.m_referenceViscosity = m_referenceViscosity;
  fluid.m_densityModelString = m_densityModelString;
  fluid.m_viscosityModelString = m_viscosityModelString;
  fluid.m_densityModelType = m_densityModelType;
  fluid.m_viscosityModelType = m_viscosityModelType;
}

void
CompressibleSinglePhaseFluid::PostProcessInput()
{
  SingleFluidBase::PostProcessInput();

  GEOSX_ERROR_IF_LT_MSG(
    m_compressibility,
    0.0,
    getName() << ": invalid value of " << viewKeyStruct::compressibilityString );

  GEOSX_ERROR_IF_LT_MSG(
    m_viscosibility,
    0.0,
    getName() << ": invalid value of " << viewKeyStruct::viscosibilityString );

  GEOSX_ERROR_IF_LE_MSG( m_referenceDensity,
                         0.0,
                         getName() << ": invalid value of "
                                   << viewKeyStruct::referenceDensityString );

  GEOSX_ERROR_IF_LE_MSG( m_referenceViscosity,
                         0.0,
                         getName() << ": invalid value of "
                                   << viewKeyStruct::referenceViscosityString );

  m_densityModelType = stringToExponentType( m_densityModelString );
  m_viscosityModelType = stringToExponentType( m_viscosityModelString );

  // Due to the way update wrapper is currently implemented, we can only support one model type

  GEOSX_ERROR_IF( m_densityModelType != ExponentApproximationType::Linear,
                  getName() << ": model type currently not supported: "
                            << m_densityModelString );

  GEOSX_ERROR_IF( m_viscosityModelType != ExponentApproximationType::Linear,
                  getName() << ": model type currently not supported: "
                            << m_viscosityModelString );

  // Set default values for derivatives (cannot be done in base class)
  // TODO: reconsider the necessity of this

  real64 dRho_dP;
  real64 dVisc_dP;
  createKernelWrapper().Compute( m_referencePressure,
                                 m_referenceDensity,
                                 dRho_dP,
                                 m_referenceViscosity,
                                 dVisc_dP );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::dDens_dPresString )
    ->setDefaultValue( dRho_dP );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::dVisc_dPresString )
    ->setDefaultValue( dVisc_dP );
}

CompressibleSinglePhaseFluid::KernelWrapper
CompressibleSinglePhaseFluid::createKernelWrapper()
{
  return KernelWrapper( KernelWrapper::DensRelationType( m_referencePressure,
                                                         m_referenceDensity,
                                                         m_compressibility ),
                        KernelWrapper::ViscRelationType( m_referencePressure,
                                                         m_referenceViscosity,
                                                         m_viscosibility ),
                        m_density,
                        m_dDensity_dPressure,
                        m_viscosity,
                        m_dViscosity_dPressure );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        CompressibleSinglePhaseFluid,
                        std::string const &,
                        Group * const )

} /* namespace constitutive */

} /* namespace geosx */
