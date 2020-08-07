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
 *  @file DelftEgg.cpp
 */

#include "DelftEgg.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

DelftEgg::DelftEgg( std::string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultRefPressure(),
  m_defaultRefStrainVol(),
  m_defaultRecompressionIndex(),
  m_defaultVirginCompressionIndex(),
  m_defaultCslSlope(),
  m_defaultShapeParameter(),
  m_defaultPreConsolidationPressure(),
  m_shearModulus(),
  m_recompressionIndex(),
  m_virginCompressionIndex(),
  m_cslSlope(),
  m_shapeParameter(),
  m_newPreConsolidationPressure(),
  m_oldPreConsolidationPressure(),
  m_newStress(),
  m_oldStress()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultShearModulusString, &m_defaultShearModulus )->
    setApplyDefaultValue( 0.6e9 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Elastic shear modulus parameter" );
  
  registerWrapper( viewKeyStruct::defaultRefPressureString, &m_defaultRefPressure )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference pressure" );
  
  registerWrapper( viewKeyStruct::defaultRefStrainVolString, &m_defaultRefStrainVol )->
    setApplyDefaultValue( 1e-6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Reference volumetric strain" );
  
  registerWrapper( viewKeyStruct::defaultRecompressionIndexString, &m_defaultRecompressionIndex )->
    setApplyDefaultValue( 2e-3 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Recompresion index" );
  
  registerWrapper( viewKeyStruct::defaultVirginCompressionIndexString, &m_defaultVirginCompressionIndex )->
    setApplyDefaultValue( 5e-3 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Virgin Compression index" );

  registerWrapper( viewKeyStruct::defaultCslSlopeString, &m_defaultCslSlope )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::defaultShapeParameterString, &m_defaultShapeParameter )->
    setApplyDefaultValue( 1.0 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::defaultPreConsolidationPressureString, &m_defaultPreConsolidationPressure )->
    setApplyDefaultValue( 5e6 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Initial preconsolidation pressure" );

  // register fields

  registerWrapper( viewKeyStruct::shearModulusString, &m_shearModulus )->
    setApplyDefaultValue( -1 )->
    setDescription( "Elastic shear modulus field" );
  
  registerWrapper( viewKeyStruct::refPressureString, &m_refPressure )->
    setApplyDefaultValue( -1 )->
    setDescription( "Reference pressure" );
  
  registerWrapper( viewKeyStruct::refStrainVolString, &m_refStrainVol )->
    setApplyDefaultValue( -1 )->
    setDescription( "Reference volumetric strain" );
  
  registerWrapper( viewKeyStruct::recompressionIndexString, &m_recompressionIndex )->
    setApplyDefaultValue( -1 )->
    setDescription( "Recompression index" );

  registerWrapper( viewKeyStruct::virginCompressionIndexString, &m_virginCompressionIndex )->
    setApplyDefaultValue( -1 )->
    setDescription( "Virgin compression index" );

  registerWrapper( viewKeyStruct::cslSlopeString, &m_cslSlope )->
    setApplyDefaultValue( -1 )->
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::shapeParameterString, &m_shapeParameter )->
    setApplyDefaultValue( -1 )->
    setDescription( "Shape parameter for the yield surface" );
  
  registerWrapper( viewKeyStruct::newPreConsolidationPressureString, &m_newPreConsolidationPressure )->
    setApplyDefaultValue( -1 )->
    setDescription( "New preconsolidation pressure" );
  
  registerWrapper( viewKeyStruct::oldPreConsolidationPressureString, &m_oldPreConsolidationPressure )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old preconsolidation pressure" );
  
  registerWrapper( viewKeyStruct::newStressString, &m_newStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "New stress field" );
  
  registerWrapper( viewKeyStruct::oldStressString, &m_oldStress )->
    setApplyDefaultValue( -1 )->
    setDescription( "Old stress field" );
}


DelftEgg::~DelftEgg()
{}


void
DelftEgg::DeliverClone( string const & name,
                                      Group * const parent,
                                      std::unique_ptr< ConstitutiveBase > & clone ) const
{
  if( !clone )
  {
    clone = std::make_unique< DelftEgg >( name, parent );
  }
  SolidBase::DeliverClone( name, parent, clone );
  DelftEgg * const newConstitutiveRelation = dynamic_cast< DelftEgg * >(clone.get());

  newConstitutiveRelation->m_defaultShearModulus     = m_defaultShearModulus;
  newConstitutiveRelation->m_defaultRefPressure = m_defaultRefPressure;
  newConstitutiveRelation->m_defaultRefStrainVol = m_defaultRefStrainVol;
  newConstitutiveRelation->m_defaultRecompressionIndex         = m_defaultRecompressionIndex;
  newConstitutiveRelation->m_defaultVirginCompressionIndex    = m_defaultVirginCompressionIndex;
  newConstitutiveRelation->m_defaultCslSlope    = m_defaultCslSlope;
  newConstitutiveRelation->m_defaultShapeParameter    = m_defaultShapeParameter;
  newConstitutiveRelation->m_defaultPreConsolidationPressure    = m_defaultPreConsolidationPressure;

  newConstitutiveRelation->m_shearModulus = m_shearModulus;
  newConstitutiveRelation->m_refPressure = m_refPressure;
  newConstitutiveRelation->m_refStrainVol = m_refStrainVol;
  newConstitutiveRelation->m_recompressionIndex = m_recompressionIndex;
  newConstitutiveRelation->m_virginCompressionIndex = m_virginCompressionIndex;
  newConstitutiveRelation->m_cslSlope = m_cslSlope;
  newConstitutiveRelation->m_shapeParameter = m_shapeParameter;
  newConstitutiveRelation->m_newPreConsolidationPressure = m_newPreConsolidationPressure;
  newConstitutiveRelation->m_oldPreConsolidationPressure = m_oldPreConsolidationPressure;
  newConstitutiveRelation->m_newStress = m_newStress;
  newConstitutiveRelation->m_oldStress = m_oldStress;
}

void DelftEgg::AllocateConstitutiveData( dataRepository::Group * const parent,
                                                       localIndex const numConstitutivePointsPerParentIndex )
{
  SolidBase::AllocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
  this->resize( parent->size() );
  
  m_shearModulus.resize( parent->size() );
  m_recompressionIndex.resize( parent->size() );
  m_virginCompressionIndex.resize( parent->size() );
  m_cslSlope.resize( parent->size() );
  m_shapeParameter.resize( parent->size() );
  
  m_newPreConsolidationPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_oldPreConsolidationPressure.resize( parent->size(), numConstitutivePointsPerParentIndex );
  m_refPressure.resize( parent->size() , numConstitutivePointsPerParentIndex);
  m_refStrainVol.resize( parent->size() , numConstitutivePointsPerParentIndex);
  
  m_newStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 );
  m_oldStress.resize( parent->size(), numConstitutivePointsPerParentIndex, 6 ); // TODO: figure out how to set initial stress
  
  // set arrays to default values
  m_shearModulus.setValues< serialPolicy >( m_defaultShearModulus );
  m_refPressure.setValues< serialPolicy >( m_defaultRefPressure );
  m_refStrainVol.setValues< serialPolicy >( m_defaultRefStrainVol );
  m_recompressionIndex.setValues< serialPolicy >( m_defaultRecompressionIndex );
  m_virginCompressionIndex.setValues< serialPolicy >( m_defaultVirginCompressionIndex);
  m_cslSlope.setValues< serialPolicy >( m_defaultCslSlope );
  m_shapeParameter.setValues< serialPolicy >( m_defaultShapeParameter );
  m_newPreConsolidationPressure.setValues< serialPolicy >( m_defaultPreConsolidationPressure );
  m_oldPreConsolidationPressure.setValues< serialPolicy >( m_defaultPreConsolidationPressure );
}

void DelftEgg::PostProcessInput()
{
  GEOSX_ASSERT_MSG(m_defaultCslSlope >= 0, "Negative slope of critical state line detected");
  GEOSX_ASSERT_MSG(m_defaultRecompressionIndex >= 0, "Negative recompresion index detected");
  GEOSX_ASSERT_MSG(m_defaultVirginCompressionIndex >= 0, "Negative virgin compression index");
  GEOSX_ASSERT_MSG(m_defaultVirginCompressionIndex>= m_defaultRecompressionIndex, "Recompression index should exceed virgin recompression index");
   
  m_postProcessed = true; // TODO: add parameter conversion helper class for more flexible input
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, DelftEgg, std::string const &, Group * const )
}
} /* namespace geosx */
