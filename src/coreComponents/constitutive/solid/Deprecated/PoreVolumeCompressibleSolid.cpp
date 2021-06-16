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
 * @file PoreVolumeCompressibleSolid.cpp
 */

#include "PoreVolumeCompressibleSolid.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{


PoreVolumeCompressibleSolid::PoreVolumeCompressibleSolid( string const & name, Group * const parent ):
  ConstitutiveBase( name, parent )
{
  registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Solid compressibility" );

  registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Reference pressure for fluid compressibility" );

  registerWrapper( viewKeyStruct::poreVolumeMultiplierString(), &m_poreVolumeMultiplier ).
    setDefaultValue( 1.0 );

  registerWrapper( viewKeyStruct::dPVMult_dPresString(), &m_dPVMult_dPressure );
}

PoreVolumeCompressibleSolid::~PoreVolumeCompressibleSolid() = default;

std::unique_ptr< ConstitutiveBase >
PoreVolumeCompressibleSolid::deliverClone( string const & name,
                                           Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = ConstitutiveBase::deliverClone( name, parent );

  PoreVolumeCompressibleSolid * const
  newConstitutiveRelation = dynamic_cast< PoreVolumeCompressibleSolid * >(clone.get());

  newConstitutiveRelation->m_poreVolumeRelation = this->m_poreVolumeRelation;

  return clone;
}

void PoreVolumeCompressibleSolid::allocateConstitutiveData( dataRepository::Group & parent,
                                                            localIndex const numConstitutivePointsPerParentIndex )
{
  ConstitutiveBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  this->resize( parent.size() );

  m_poreVolumeMultiplier.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent.size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier.setValues< serialPolicy >( 1.0 );
}

void PoreVolumeCompressibleSolid::postProcessInput()
{
  if( m_compressibility < 0.0 )
  {
    string const message = "An invalid value of fluid bulk modulus (" + std::to_string( m_compressibility ) + ") is specified";
    GEOSX_ERROR( message );
  }
  m_poreVolumeRelation.setCoefficients( m_referencePressure, 1.0, m_compressibility );
}

void PoreVolumeCompressibleSolid::stateUpdatePointPressure( real64 const & pres,
                                                            localIndex const k,
                                                            localIndex const q )
{
  m_poreVolumeRelation.compute( pres, m_poreVolumeMultiplier[k][q], m_dPVMult_dPressure[k][q] );
}

void PoreVolumeCompressibleSolid::stateUpdateBatchPressure( arrayView1d< real64 const > const & pres,
                                                            arrayView1d< real64 const > const & dPres )
{
  localIndex const numElems = m_poreVolumeMultiplier.size( 0 );
  localIndex const numQuad  = m_poreVolumeMultiplier.size( 1 );

  GEOSX_ASSERT_EQ( pres.size(), numElems );
  GEOSX_ASSERT_EQ( dPres.size(), numElems );

  ExponentialRelation< real64, ExponentApproximationType::Linear > const relation = m_poreVolumeRelation;

  arrayView2d< real64 > const & pvmult = m_poreVolumeMultiplier;
  arrayView2d< real64 > const & dPVMult_dPres = m_dPVMult_dPressure;

  forAll< parallelDevicePolicy<> >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQuad; ++q )
    {
      relation.compute( pres[k] + dPres[k], pvmult[k][q], dPVMult_dPres[k][q] );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoreVolumeCompressibleSolid, string const &, Group * const )
}
} /* namespace geosx */
