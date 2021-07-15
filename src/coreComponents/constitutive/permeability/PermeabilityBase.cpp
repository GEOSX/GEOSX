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
 * @file PermeabilityBase.cpp
 */

#include "../permeability/PermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PermeabilityBase::PermeabilityBase( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent ),
  m_permeability(),
  m_dPerm_dPressure()
{
  registerWrapper( viewKeyStruct::permeabilityString(), &m_permeability ).
    setPlotLevel( PlotLevel::LEVEL_0 ).
    setDescription( " permeability of the rock." ).
    setApplyDefaultValue( -1.0 );   // will be overwritten

  registerWrapper( viewKeyStruct::dPerm_dPressureString(), &m_dPerm_dPressure ).
    setPlotLevel( PlotLevel::LEVEL_3 ).
    setDescription( " dPerm_dPressure of the rock." ).
    setApplyDefaultValue( 0.0 );     // will be overwritten
}

PermeabilityBase::~PermeabilityBase() = default;

std::unique_ptr< ConstitutiveBase >
PermeabilityBase::deliverClone( string const & name,
                                Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void PermeabilityBase::allocateConstitutiveData( dataRepository::Group & parent,
                                                 localIndex const numConstitutivePointsPerParentIndex )
{
  m_permeability.resize( 0, 3 );
  m_dPerm_dPressure.resize( 0, 3 );

  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void PermeabilityBase::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PermeabilityBase, string const &, Group * const )
}
} /* namespace geosx */
