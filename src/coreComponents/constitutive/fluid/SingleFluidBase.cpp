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
 * @file SingleFluidBase.cpp
 */

#include "SingleFluidBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
SingleFluidBase::SingleFluidBase( std::string const & name, Group * const parent ) :
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Default value for density." );

  registerWrapper( viewKeyStruct::defaultViscosityString, &m_defaultViscosity )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Default value for viscosity." );

  registerWrapper( viewKeyStruct::densityString, &m_density )
    ->setPlotLevel( PlotLevel::LEVEL_0 );

  registerWrapper( viewKeyStruct::dDens_dPresString, &m_dDensity_dPressure );

  registerWrapper( viewKeyStruct::viscosityString, &m_viscosity )
    ->setPlotLevel( PlotLevel::LEVEL_0 );
  registerWrapper( viewKeyStruct::dVisc_dPresString, &m_dViscosity_dPressure );
}

SingleFluidBase::~SingleFluidBase() = default;

void
SingleFluidBase::PostProcessInput()
{
  ConstitutiveBase::PostProcessInput();
  this->getWrapper< array2d< real64 > >( viewKeyStruct::densityString )
    ->setApplyDefaultValue( m_defaultDensity );
  this->getWrapper< array2d< real64 > >( viewKeyStruct::viscosityString )
    ->setApplyDefaultValue( m_defaultViscosity );
}

void
SingleFluidBase::AllocateConstitutiveData(
  Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent,
                                              numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );

  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dDensity_dPressure.resize( parent->size(),
                               numConstitutivePointsPerParentIndex );

  m_viscosity.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_dViscosity_dPressure.resize( parent->size(),
                                 numConstitutivePointsPerParentIndex );
}

void
SingleFluidBase::DeliverClone( string const & name,
                               Group * const parent,
                               std::unique_ptr< ConstitutiveBase > & clone ) const
{
  GEOSX_ERROR_IF( !clone, "clone not allocated" );

  ConstitutiveBase::DeliverClone( name, parent, clone );
  SingleFluidBase & newConstitutiveRelation =
    dynamicCast< SingleFluidBase & >( *clone );

  newConstitutiveRelation.m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation.m_defaultViscosity = m_defaultViscosity;
  newConstitutiveRelation.m_density = m_density;
  newConstitutiveRelation.m_viscosity = m_viscosity;
  newConstitutiveRelation.m_dDensity_dPressure = m_dDensity_dPressure;
  newConstitutiveRelation.m_dViscosity_dPressure = m_dViscosity_dPressure;
}
}  //namespace constitutive

}  //namespace geosx
