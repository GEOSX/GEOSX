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
 * @file PoroElastic.cpp
 */

#include "PoroElastic.hpp"

#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
PoroElastic< BASE >::PoroElastic( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_compressibility(),
  m_referencePressure(),
  m_biotCoefficient(),
  m_poreVolumeMultiplier(),
  m_dPVMult_dPressure(),
  m_poreVolumeRelation()
{
  this->registerWrapper( viewKeyStruct::biotCoefficientString(), &m_biotCoefficient ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Biot's coefficient" );

  this->registerWrapper( viewKeyStruct::compressibilityString(), &m_compressibility ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Pore volume compressibilty" );

  this->registerWrapper( viewKeyStruct::referencePressureString(), &m_referencePressure ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "ReferencePressure" );


  this->registerWrapper( viewKeyStruct::poreVolumeMultiplierString(), &m_poreVolumeMultiplier ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "" );

  this->registerWrapper( viewKeyStruct::dPVMult_dPresString(), &m_dPVMult_dPressure ).
    setApplyDefaultValue( -1 ).
    setDescription( "" );
}

template< typename BASE >
PoroElastic< BASE >::~PoroElastic()
{}

template< typename BASE >
void PoroElastic< BASE >::postProcessInput()
{
  BASE::postProcessInput();
  m_poreVolumeRelation.setCoefficients( m_referencePressure, 1.0, m_compressibility );
}

template< typename BASE >
std::unique_ptr< ConstitutiveBase >
PoroElastic< BASE >::deliverClone( string const & name,
                                   dataRepository::Group * const parent ) const
{
  std::unique_ptr< ConstitutiveBase > clone = BASE::deliverClone( name, parent );
  PoroElastic< BASE > & castedClone = dynamic_cast< PoroElastic< BASE > & >( *clone );
  castedClone.m_poreVolumeRelation = m_poreVolumeRelation;

  return clone;
}

template< typename BASE >
void PoroElastic< BASE >::allocateConstitutiveData( dataRepository::Group & parent,
                                                    localIndex const numConstitutivePointsPerParentIndex )
{
  m_poreVolumeMultiplier.resize( 0, numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( 0, numConstitutivePointsPerParentIndex );
  BASE::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}

template< typename BASE >
void PoroElastic< BASE >::stateUpdateBatchPressure( arrayView1d< real64 const > const & pres,
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

typedef PoroElastic< ElasticIsotropic > PoroElasticIsotropic;
typedef PoroElastic< ElasticTransverseIsotropic > PoroElasticTransverseIsotropic;
typedef PoroElastic< DruckerPrager > PoroDruckerPrager;
typedef PoroElastic< DruckerPragerExtended > PoroDruckerPragerExtended;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PoroDruckerPragerExtended, string const &, Group * const )

}
} /* namespace geosx */
