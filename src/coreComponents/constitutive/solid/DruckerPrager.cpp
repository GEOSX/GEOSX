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
 *  @file DruckerPrager.cpp
 */

#include "DruckerPrager.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

DruckerPrager::DruckerPrager( std::string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultBulkModulus(),
  m_defaultShearModulus(),
  m_defaultTanFrictionAngle(),
  m_defaultTanDilationAngle(),
  m_defaultCohesion(),
  m_defaultHardeningRate(),
  m_bulkModulus(),
  m_shearModulus(),
  m_tanFrictionAngle(),
  m_tanDilationAngle(),
  m_hardeningRate(),
  m_newCohesion(),
  m_oldCohesion(),
  m_newStress(),
  m_oldStress()
{
  // register default values
  
  registerWrapper( viewKeyStruct::defaultBulkModulusString, &m_defaultBulkModulus )->
    setApplyDefaultValue( 1e9 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic bulk modulus parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusString, &m_defaultShearModulus )->
    setApplyDefaultValue( 0.6e9 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic shear modulus parameter" );
  
  registerWrapper( viewKeyStruct::defaultTanFrictionAngleString, &m_defaultTanFrictionAngle )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Yield surface slope parameter tan(phi)" );
  
  registerWrapper( viewKeyStruct::defaultTanDilationAngleString, &m_defaultTanDilationAngle )->
    setApplyDefaultValue( 0.5 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Plastic potential slope parameter tan(psi)" );
  
  registerWrapper( viewKeyStruct::defaultHardeningRateString, &m_defaultHardeningRate )->
    setApplyDefaultValue( 1e8 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Cohesion hardening/softening rate parameter" );
  
  registerWrapper( viewKeyStruct::defaultCohesionString, &m_defaultCohesion )->
    setApplyDefaultValue( 5e6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial cohesion parameter" );

  // register fields
  
  registerWrapper( viewKeyStruct::bulkModulusString, &m_bulkModulus )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic bulk modulus field" );

  registerWrapper( viewKeyStruct::shearModulusString, &m_shearModulus )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic shear modulus field" );
  
  registerWrapper( viewKeyStruct::tanFrictionAngleString, &m_tanFrictionAngle )->
    setApplyDefaultValue( -1 )->
    setDescription( "Yield surface slope tan(phi) field" );
  
  registerWrapper( viewKeyStruct::tanDilationAngleString, &m_tanDilationAngle )->
    setApplyDefaultValue( -1 )->
    setDescription( "Plastic potential slope tan(psi) field" );
  
  registerWrapper( viewKeyStruct::hardeningRateString, &m_hardeningRate )->
    setApplyDefaultValue( -1 )->
    setDescription( "Hardening rate field" );
  
  registerWrapper( viewKeyStruct::newCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "New cohesion field" );
  
  registerWrapper( viewKeyStruct::oldCohesionString, &m_newCohesion )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old cohesion field" );
  
  registerWrapper( viewKeyStruct::newStressString, &m_newStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "New stress field" );
  
  registerWrapper( viewKeyStruct::oldStressString, &m_oldStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old stress field" );
}


DruckerPrager::~DruckerPrager()
{}


void
DruckerPrager::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< DruckerPrager >( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  DruckerPrager * const newConstitutiveRelation = dynamic_cast< DruckerPrager * >(clone.get());

  newConstitutiveRelation->m_defaultBulkModulus      = m_defaultBulkModulus;
  newConstitutiveRelation->m_defaultShearModulus     = m_defaultShearModulus;
  newConstitutiveRelation->m_defaultTanFrictionAngle = m_defaultTanFrictionAngle;
  newConstitutiveRelation->m_defaultTanDilationAngle = m_defaultTanDilationAngle;
  newConstitutiveRelation->m_defaultCohesion         = m_defaultCohesion;
  newConstitutiveRelation->m_defaultHardeningRate    = m_defaultHardeningRate;
  
  newConstitutiveRelation->m_bulkModulus = m_bulkModulus;
  newConstitutiveRelation->m_shearModulus = m_shearModulus;
  newConstitutiveRelation->m_tanFrictionAngle = m_tanFrictionAngle;
  newConstitutiveRelation->m_tanDilationAngle = m_tanDilationAngle;
  newConstitutiveRelation->m_hardeningRate = m_hardeningRate;
  newConstitutiveRelation->m_newCohesion = m_newCohesion;
  newConstitutiveRelation->m_oldCohesion = m_oldCohesion;
  newConstitutiveRelation->m_newStress = m_newStress;
  newConstitutiveRelation->m_oldStress = m_oldStress;
}

void DruckerPrager::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  this->resize( parent->size() );
  
  m_bulkModulus.resize( parent->size() );
  m_shearModulus.resize( parent->size() );
  m_tanFrictionAngle.resize( parent->size() );
  m_tanDilationAngle.resize( parent->size() );
  m_hardeningRate.resize( parent->size() );
  
  m_newCohesion.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_oldCohesion.resize( parent->size(), numConstitutivePointsPerParentIndex );
  
  m_newStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 ); // TODO: figure out how to set initial stress
  
  // set arrays to default values
  m_bulkModulus.setValues< serialPolicy >( m_defaultBulkModulus );
  m_shearModulus.setValues< serialPolicy >( m_defaultShearModulus );
  m_tanFrictionAngle.setValues< serialPolicy >( m_defaultTanFrictionAngle );
  m_tanDilationAngle.setValues< serialPolicy >( m_defaultTanDilationAngle );
  m_hardeningRate.setValues< serialPolicy >( m_defaultHardeningRate );
  m_newCohesion.setValues< serialPolicy >( m_defaultCohesion );
  m_oldCohesion.setValues< serialPolicy >( m_defaultCohesion );
}

void DruckerPrager::PostProcessInput()
{
  GEOSX_ASSERT_MSG(m_defaultCohesion >= 0, "Negative cohesion value detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= 0, "Negative friction angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanDilationAngle >= 0, "Negative dilation angle detected");
  GEOSX_ASSERT_MSG(m_defaultTanFrictionAngle >= m_defaultTanDilationAngle, "Friction angle should exceed dilation angle");
  
  m_postProcessed = true; // TODO: add parameter conversion helper class for more flexible input
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DruckerPrager, std::string const &, Group * const )
}
} /* namespace geosx */
