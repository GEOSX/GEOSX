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
 * @file ParticleFluidBase.cpp
 */

#include "ParticleFluidBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
ParticleFluidBase::ParticleFluidBase( std::string const & name, Group * const parent ) :
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::settlingFactorString, &m_settlingFactor );
  registerWrapper( viewKeyStruct::dSettlingFactor_dPressureString,
                   &m_dSettlingFactor_dPressure );
  registerWrapper( viewKeyStruct::dSettlingFactor_dProppantConcentrationString,
                   &m_dSettlingFactor_dProppantConcentration );
  registerWrapper( viewKeyStruct::dSettlingFactor_dComponentConcentrationString,
                   &m_dSettlingFactor_dComponentConcentration );

  registerWrapper( viewKeyStruct::collisionFactorString, &m_collisionFactor );
  registerWrapper( viewKeyStruct::dCollisionFactor_dProppantConcentrationString,
                   &m_dCollisionFactor_dProppantConcentration );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString,
                   &m_maxProppantConcentration )
    ->setApplyDefaultValue( 0.6 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Max proppant concentration" );

  registerWrapper( viewKeyStruct::isCollisionalSlipString, &m_isCollisionalSlip )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Whether the collisional component of the slip velocity is considered" );

  registerWrapper( viewKeyStruct::proppantPackPermeabilityString,
                   &m_proppantPackPermeability );
}

ParticleFluidBase::~ParticleFluidBase() = default;

void
ParticleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
}

void
ParticleFluidBase::AllocateConstitutiveData(
  Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent,
                                              numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_dSettlingFactor_dComponentConcentration.resize( parent->size(),
                                                    MAX_NUM_COMPONENTS );
}

void
ParticleFluidBase::DeliverClone( string const & name,
                                 Group * const parent,
                                 std::unique_ptr< ConstitutiveBase > & clone ) const
{
  GEOSX_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  ParticleFluidBase & fluid = dynamicCast< ParticleFluidBase & >( *clone );

  fluid.m_settlingFactor = m_settlingFactor;
  fluid.m_dSettlingFactor_dPressure = m_dSettlingFactor_dPressure;
  fluid.m_dSettlingFactor_dProppantConcentration =
    m_dSettlingFactor_dProppantConcentration;
  fluid.m_dSettlingFactor_dComponentConcentration =
    m_dSettlingFactor_dComponentConcentration;

  fluid.m_collisionFactor = m_collisionFactor;
  fluid.m_dCollisionFactor_dProppantConcentration =
    m_dCollisionFactor_dProppantConcentration;

  fluid.m_maxProppantConcentration = m_maxProppantConcentration;
  fluid.m_isCollisionalSlip = m_isCollisionalSlip;
  fluid.m_proppantPackPermeability = m_proppantPackPermeability;
}

}  //namespace constitutive

}  //namespace geosx
