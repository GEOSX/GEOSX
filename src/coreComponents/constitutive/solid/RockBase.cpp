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
 * @file RockBase.cpp
 */

#include "RockBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


RockBase::RockBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_newPorosity(),
  m_oldPorosity(),
  m_dPorosity_dPressure(),
  m_referencePorosity(),
  m_defaultReferencePorosity(),
  m_grainBulkModulus(),
  m_grainDensity()
{
  registerWrapper( viewKeyStruct::newPorosityString(), &m_newPorosity ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setApplyDefaultValue( 0.2 ); // will be overwritten

  registerWrapper( viewKeyStruct::oldPorosityString(), &m_oldPorosity ).
    setApplyDefaultValue( -1.0 );// will be overwritten

  registerWrapper( viewKeyStruct::dPorosity_dPressureString(), &m_dPorosity_dPressure ).
    setApplyDefaultValue( 0.0 );// will be overwritten

  registerWrapper( viewKeyStruct::defaultRefererencePorosityString(), &m_defaultReferencePorosity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default value of the reference porosity" );

  registerWrapper( viewKeyStruct::referencePorosityString(), &m_referencePorosity ).
    setApplyDefaultValue( 1.0 );

  registerWrapper( viewKeyStruct::grainBulkModulusString(), &m_grainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Grain bulk modulus" );

  registerWrapper( viewKeyStruct::grainDensityString(), &m_grainDensity ).
    setDescription( "Grain density" );

  registerWrapper( viewKeyStruct::defaultGrainDensityString(), &m_defaultGrainDensity ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Default grain density. It's used to set the default value of the grain density." );
}

RockBase::~RockBase() = default;

std::unique_ptr< ConstitutiveBase >
RockBase::deliverClone( string const & name,
                        Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void RockBase::allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldPorosity.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPorosity_dPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_grainDensity.resize( 0, numConstitutivePointsPerParentIndex );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void RockBase::postProcessInput()
{
  this->getWrapper< array2d< real64 > >( viewKeyStruct::grainDensityString() ).
    setApplyDefaultValue( m_defaultGrainDensity );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
    setApplyDefaultValue( m_defaultReferencePorosity );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, RockBase, string const &, Group * const )
}
} /* namespace geosx */
