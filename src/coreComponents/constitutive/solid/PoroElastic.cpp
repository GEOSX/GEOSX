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
 * @file PoroElastic.cpp
 */

#include "PoroElastic.hpp"

#include "LinearElasticAnisotropic.hpp"
#include "LinearElasticIsotropic.hpp"
#include "LinearElasticTransverseIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{
template< typename BASE >
PoroElastic< BASE >::PoroElastic( string const & name, Group * const parent ) :
  BASE( name, parent ),
  m_compressibility(),
  m_referencePressure(),
  m_biotCoefficient(),
  m_poreVolumeMultiplier(),
  m_dPVMult_dPressure(),
  m_poreVolumeRelation()
{
  this
    ->registerWrapper( viewKeyStruct::biotCoefficientString, &m_biotCoefficient )
    ->setApplyDefaultValue( 1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Biot's coefficient" );

  this
    ->registerWrapper( viewKeyStruct::compressibilityString, &m_compressibility )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Pore volume compressibilty" );

  this
    ->registerWrapper( viewKeyStruct::referencePressureString, &m_referencePressure )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "ReferencePressure" );

  this
    ->registerWrapper( viewKeyStruct::poreVolumeMultiplierString,
                       &m_poreVolumeMultiplier )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "" );

  this
    ->registerWrapper( viewKeyStruct::dPVMult_dPresString, &m_dPVMult_dPressure )
    ->setApplyDefaultValue( -1 )
    ->setDescription( "" );
}

template< typename BASE >
PoroElastic< BASE >::~PoroElastic()
{}

template< typename BASE >
void
PoroElastic< BASE >::PostProcessInput()
{
  //    m_compressibility = 1 / K;

  if( m_compressibility <= 0 )
  {
    //    string const message = std::to_string( numConstantsSpecified ) + " Elastic Constants Specified. Must specify 2
    // constants!";
    //    GEOSX_ERROR( message );
  }
  m_poreVolumeRelation.SetCoefficients( m_referencePressure, 1.0, m_compressibility );
}

template< typename BASE >
void
PoroElastic< BASE >::DeliverClone( string const & name,
                                   Group * const parent,
                                   std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< PoroElastic< BASE > >( name, parent );
  }
  BASE::DeliverClone( name, parent, clone );
  PoroElastic< BASE > * const newConstitutiveRelation =
    dynamic_cast< PoroElastic< BASE > * >( clone.get() );

  newConstitutiveRelation->m_compressibility = m_compressibility;
  newConstitutiveRelation->m_referencePressure = m_referencePressure;
  newConstitutiveRelation->m_biotCoefficient = m_biotCoefficient;
  newConstitutiveRelation->m_poreVolumeMultiplier = m_poreVolumeMultiplier;
  newConstitutiveRelation->m_dPVMult_dPressure = m_dPVMult_dPressure;
  newConstitutiveRelation->m_poreVolumeRelation = m_poreVolumeRelation;
}

template< typename BASE >
void
PoroElastic< BASE >::AllocateConstitutiveData(
  dataRepository::Group * const parent,
  localIndex const numConstitutivePointsPerParentIndex )
{
  BASE::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_poreVolumeMultiplier.resize( parent->size(),
                                 numConstitutivePointsPerParentIndex );
  m_dPVMult_dPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_poreVolumeMultiplier.setValues< serialPolicy >( 1.0 );
}

template< typename BASE >
void
PoroElastic< BASE >::StateUpdateBatchPressure(
  arrayView1d< real64 const > const & pres,
  arrayView1d< real64 const > const & dPres )
{
  localIndex const numElems = m_poreVolumeMultiplier.size( 0 );
  localIndex const numQuad = m_poreVolumeMultiplier.size( 1 );

  GEOSX_ASSERT_EQ( pres.size(), numElems );
  GEOSX_ASSERT_EQ( dPres.size(), numElems );

  ExponentialRelation< real64, ExponentApproximationType::Linear > const relation =
    m_poreVolumeRelation;

  arrayView2d< real64 > const & pvmult = m_poreVolumeMultiplier;
  arrayView2d< real64 > const & dPVMult_dPres = m_dPVMult_dPressure;

  forAll< parallelDevicePolicy<> >(
    numElems,
    [=] GEOSX_HOST_DEVICE( localIndex const k ) {
      for( localIndex q = 0; q < numQuad; ++q )
      {
        relation.Compute( pres[k] + dPres[k], pvmult[k][q], dPVMult_dPres[k][q] );
      }
    } );
}

typedef PoroElastic< LinearElasticIsotropic > PoroLinearElasticIsotropic;
typedef PoroElastic< LinearElasticAnisotropic > PoroLinearElasticAnisotropic;
typedef PoroElastic< LinearElasticTransverseIsotropic > PoroLinearElasticTransverseIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        PoroLinearElasticIsotropic,
                        string const &,
                        Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        PoroLinearElasticAnisotropic,
                        string const &,
                        Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase,
                        PoroLinearElasticTransverseIsotropic,
                        string const &,
                        Group * const )

}  // namespace constitutive
} /* namespace geosx */
