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
 * @file ConstantPermeability.cpp
 */

#include "ConstantPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


ConstantPermeability::ConstantPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{
  registerWrapper( viewKeyStruct::permeabilityComponentsString(), &m_permeabilityComponents ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "xx, yy and zz components of a diagonal permeability tensor." );
}

ConstantPermeability::~ConstantPermeability() = default;

std::unique_ptr< ConstitutiveBase >
ConstantPermeability::deliverClone( string const & name,
                                    Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void ConstantPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                     localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  for( localIndex ei=0; ei < parent.size(); ei++ )
  {
    for( localIndex q=0; q < numConstitutivePointsPerParentIndex; q++ )
    {
      m_permeability[ei][q][0] =  m_permeabilityComponents[0];
      m_permeability[ei][q][1] =  m_permeabilityComponents[1];
      m_permeability[ei][q][2] =  m_permeabilityComponents[2];
    }
  }
}

void ConstantPermeability::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ConstantPermeability, string const &, Group * const )

}
} /* namespace geosx */
