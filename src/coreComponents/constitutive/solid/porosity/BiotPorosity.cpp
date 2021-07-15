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
 * @file BiotPorosity.cpp
 */

#include "BiotPorosity.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


BiotPorosity::BiotPorosity( string const & name, Group * const parent ):
  PorosityBase( name, parent )
{
  registerWrapper( viewKeyStruct::grainBulkModulusString(), &m_grainBulkModulus ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Grain bulk modulus" );

  registerWrapper( viewKeyStruct::biotCoefficientString(), &m_biotCoefficient ).
    setDescription( "Biot coefficient." );
}

BiotPorosity::~BiotPorosity() = default;

std::unique_ptr< ConstitutiveBase >
BiotPorosity::deliverClone( string const & name,
                            Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  return clone;
}

void BiotPorosity::allocateConstitutiveData( dataRepository::Group & parent,
                                             localIndex const numConstitutivePointsPerParentIndex )
{
  PorosityBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

void BiotPorosity::postProcessInput()
{
  PorosityBase::postProcessInput();
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, BiotPorosity, string const &, Group * const )
}
} /* namespace geosx */
