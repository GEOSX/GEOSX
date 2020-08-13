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
 * @file SolidBase.cpp
 */

#include "SolidBase.hpp"

namespace geosx
{
using namespace dataRepository;

namespace constitutive
{
SolidBase::SolidBase( string const & name, Group * const parent ) :
  ConstitutiveBase( name, parent ),
  m_defaultDensity( 0 ),
  m_density(),
  m_stress( 0, 0, 6 ),
  m_postProcessed( false )
{
  registerWrapper( viewKeyStruct::defaultDensityString, &m_defaultDensity )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Default Material Density" );

  registerWrapper( viewKeyStruct::densityString, &m_density )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "Material Density" );

  registerWrapper( viewKeyStruct::stressString, &m_stress )
    ->setPlotLevel( PlotLevel::LEVEL_0 )
    ->setDescription( "Material Stress" );
}

SolidBase::~SolidBase()
{}

void
SolidBase::DeliverClone( string const & GEOSX_UNUSED_PARAM( name ),
                         Group * const GEOSX_UNUSED_PARAM( parent ),
                         std::unique_ptr< ConstitutiveBase > & clone ) const
{
  SolidBase * const newConstitutiveRelation =
    dynamic_cast< SolidBase * >( clone.get() );

  newConstitutiveRelation->m_defaultDensity = m_defaultDensity;
  newConstitutiveRelation->m_density = m_density;

  newConstitutiveRelation->m_stress = m_stress;
}

void
SolidBase::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::AllocateConstitutiveData( parent,
                                              numConstitutivePointsPerParentIndex );

  this->resize( parent->size() );
  m_density.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_density.setValues< serialPolicy >( m_defaultDensity );

  m_stress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
}

}  // namespace constitutive
} /* namespace geosx */
