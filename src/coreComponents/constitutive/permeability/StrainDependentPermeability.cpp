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
 * @file StrainDependentPermeability.cpp
 */

#include "StrainDependentPermeability.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


StrainDependentPermeability::StrainDependentPermeability( string const & name, Group * const parent ):
  PermeabilityBase( name, parent )
{}

StrainDependentPermeability::~StrainDependentPermeability() = default;

std::unique_ptr< ConstitutiveBase >
StrainDependentPermeability::deliverClone( string const & name,
                                           Group * const parent ) const
{
  return ConstitutiveBase::deliverClone( name, parent );
}

void StrainDependentPermeability::allocateConstitutiveData( dataRepository::Group & parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  PermeabilityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void StrainDependentPermeability::postProcessInput()
{}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, StrainDependentPermeability, string const &, Group * const )

}
} /* namespace geosx */
