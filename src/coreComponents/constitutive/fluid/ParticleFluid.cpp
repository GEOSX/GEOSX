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
 * @file ParticleFluid.cpp
 */

#include "ParticleFluid.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
namespace
{
ParticleSettlingModel
stringToParticleSettlingModel( string const & str )
{
  if( str == "Stokes" )
  {
    return ParticleSettlingModel::Stokes;
  }
  else if( str == "Intermediate" )
  {
    return ParticleSettlingModel::Intermediate;
  }
  else if( str == "Turbulence" )
  {
    return ParticleSettlingModel::Turbulence;
  }
  else
  {
    GEOSX_ERROR( "Unrecognized particle settling velocity model: " << str );
  }
  return ParticleSettlingModel::Stokes;
}

}  // namespace

ParticleFluid::ParticleFluid( std::string const & name, Group * const parent ) :
  ParticleFluidBase( name, parent )
{
  registerWrapper( viewKeyStruct::particleSettlingModelString,
                   &m_particleSettlingModelString )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Particle settling velocity model" );

  registerWrapper( viewKeyStruct::proppantDensityString, &m_proppantDensity )
    ->setApplyDefaultValue( 1400.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::fluidViscosityString, &m_fluidViscosity )
    ->setApplyDefaultValue( 0.001 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Fluid viscosity" );

  registerWrapper( viewKeyStruct::proppantDiameterString, &m_proppantDiameter )
    ->setApplyDefaultValue( 200e-6 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::hinderedSettlingCoefficientString,
                   &m_hinderedSettlingCoefficient )
    ->setApplyDefaultValue( 5.9 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Hindered settling coefficient" );

  registerWrapper( viewKeyStruct::collisionAlphaString, &m_collisionAlpha )
    ->setApplyDefaultValue( 1.27 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Collision alpha coefficient" );

  registerWrapper( viewKeyStruct::slipConcentrationString, &m_slipConcentration )
    ->setApplyDefaultValue( 0.1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Slip concentration" );

  registerWrapper( viewKeyStruct::collisionBetaString, &m_collisionBeta )
    ->setApplyDefaultValue( 1.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Collision beta coefficient" );

  registerWrapper( viewKeyStruct::sphericityString, &m_sphericity )
    ->setApplyDefaultValue( 1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Sphericity" );
}

ParticleFluid::~ParticleFluid() = default;

void
ParticleFluid::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ParticleFluidBase::AllocateConstitutiveData(
    parent,
    numConstitutivePointsPerParentIndex );
}

void
ParticleFluid::DeliverClone( string const & name,
                             Group * const parent,
                             std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< ParticleFluid >( name, parent );
  }
  ParticleFluidBase::DeliverClone( name, parent, clone );
  ParticleFluid & fluid = dynamicCast< ParticleFluid & >( *clone );

  fluid.m_proppantDensity = m_proppantDensity;
  fluid.m_fluidViscosity = m_fluidViscosity;
  fluid.m_proppantDiameter = m_proppantDiameter;

  fluid.m_hinderedSettlingCoefficient = m_hinderedSettlingCoefficient;
  fluid.m_collisionAlpha = m_collisionAlpha;
  fluid.m_slipConcentration = m_slipConcentration;
  fluid.m_collisionBeta = m_collisionBeta;

  fluid.m_sphericity = m_sphericity;

  fluid.m_particleSettlingModelString = m_particleSettlingModelString;
  fluid.m_particleSettlingModel = m_particleSettlingModel;
}

void
ParticleFluid::PostProcessInput()
{
  ParticleFluidBase::PostProcessInput();

  GEOSX_ERROR_IF( m_particleSettlingModelString.empty(),
                  "Invalid particle settling model in ParticleFluid " );

  GEOSX_ERROR_IF(
    m_proppantDensity < 500.0,
    "Invalid proppantDensity in ParticleFluid, which must >= 500.0 " );

  GEOSX_ERROR_IF(
    m_proppantDiameter < 10e-6,
    "Invalid proppantDiameter in ParticleFluid, which must >= 10e-6 " );

  GEOSX_ERROR_IF(
    m_hinderedSettlingCoefficient < 0.0 || m_hinderedSettlingCoefficient > 10.0,
    "Invalid hinderedSettlingCoefficient in ParticleFluid, which must between "
    "0 and 10 " );

  GEOSX_ERROR_IF( m_collisionAlpha < 1.0,
                  "Invalid collisionAlpha in ParticleFluid, which must >= 1 " );

  GEOSX_ERROR_IF( m_collisionBeta < 0.0,
                  "Invalid collisionBeta in ParticleFluid, which must >= 0" );

  GEOSX_ERROR_IF(
    m_slipConcentration > 0.3,
    "Invalid slipConcentration in ParticleFluid, which must <= 0.3" );

  m_packPermeabilityCoef = pow( m_sphericity * m_proppantDiameter, 2.0 ) / 180.0;

  m_particleSettlingModel =
    stringToParticleSettlingModel( m_particleSettlingModelString );
}

ParticleFluid::KernelWrapper
ParticleFluid::createKernelWrapper() const
{
  return KernelWrapper( m_particleSettlingModel,
                        m_proppantDensity,
                        m_proppantDiameter,
                        m_hinderedSettlingCoefficient,
                        m_collisionAlpha,
                        m_slipConcentration,
                        m_collisionBeta,
                        m_isCollisionalSlip,
                        m_maxProppantConcentration,
                        m_settlingFactor,
                        m_dSettlingFactor_dPressure,
                        m_dSettlingFactor_dProppantConcentration,
                        m_dSettlingFactor_dComponentConcentration,
                        m_collisionFactor,
                        m_dCollisionFactor_dProppantConcentration,
                        m_proppantPackPermeability );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        ParticleFluid,
                        std::string const &,
                        Group * const )

} /* namespace constitutive */

} /* namespace geosx */
