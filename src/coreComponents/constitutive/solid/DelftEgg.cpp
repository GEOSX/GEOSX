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
 *  @file ModifiedCamClay.cpp
 */

#include "ModifiedCamClay.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ModifiedCamClay::ModifiedCamClay( string const & name, Group * const parent ):
  ElasticIsotropicPressureDependent( name, parent ),
  m_defaultVirginCompressionIndex(),
  m_defaultCslSlope(),
  m_defaultShapeParameter(),
  m_defaultPreConsolidationPressure(),
  m_virginCompressionIndex(),
  m_cslSlope(),
  m_shapeParameter(),
  m_newPreConsolidationPressure(),
  m_oldPreConsolidationPressure()
{
  // register default values

  registerWrapper( viewKeyStruct::defaultVirginCompressionIndexString(), &m_defaultVirginCompressionIndex ).
    setApplyDefaultValue( 5e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Virgin compression index" );

  registerWrapper( viewKeyStruct::defaultCslSlopeString(), &m_defaultCslSlope ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::defaultShapeParameterString(), &m_defaultShapeParameter ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::defaultPreConsolidationPressureString(), &m_defaultPreConsolidationPressure ).
    setApplyDefaultValue( -1.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial preconsolidation pressure" );

  // register fields

  registerWrapper( viewKeyStruct::virginCompressionIndexString(), &m_virginCompressionIndex ).
    setApplyDefaultValue( -1 ).
    setDescription( "Virgin compression index" );

  registerWrapper( viewKeyStruct::cslSlopeString(), &m_cslSlope ).
    setApplyDefaultValue( -1 ).
    setDescription( "Slope of the critical state line" );

  registerWrapper( viewKeyStruct::shapeParameterString(), &m_shapeParameter ).
    setApplyDefaultValue( -1 ).
    setDescription( "Shape parameter for the yield surface" );

  registerWrapper( viewKeyStruct::newPreConsolidationPressureString(), &m_newPreConsolidationPressure ).
    setApplyDefaultValue( -1 ).
    setPlotLevel( dataRepository::PlotLevel::LEVEL_3 ).
    setDescription( "New preconsolidation pressure" );

  registerWrapper( viewKeyStruct::oldPreConsolidationPressureString(), &m_oldPreConsolidationPressure ).
    setApplyDefaultValue( -1 ).
    setDescription( "Old preconsolidation pressure" );
}


ModifiedCamClay::~ModifiedCamClay()
{}


void ModifiedCamClay::allocateConstitutiveData( dataRepository::Group & parent,
                                                localIndex const numConstitutivePointsPerParentIndex )
{
  m_newPreConsolidationPressure.resize( 0, numConstitutivePointsPerParentIndex );
  m_oldPreConsolidationPressure.resize( 0, numConstitutivePointsPerParentIndex );

  ElasticIsotropicPressureDependent::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );
}


void ModifiedCamClay::postProcessInput()
{
  ElasticIsotropicPressureDependent::postProcessInput();

  GEOSX_THROW_IF( m_defaultCslSlope <= 0, "Non-positive slope of critical state line detected", InputError );
  GEOSX_THROW_IF( m_defaultShapeParameter < 1., "Shape parameter for yield surface must be greater than or equal to one", InputError );
  GEOSX_THROW_IF( m_defaultVirginCompressionIndex <= 0, "Non-positive virgin compression index detected", InputError );
  GEOSX_THROW_IF( m_defaultVirginCompressionIndex <= m_defaultRecompressionIndex, "Recompression index should exceed virgin recompression index", InputError );

//  real64 poisson = conversions::BulkModAndShearMod::toPoissonRatio( -1*m_defaultRefPressure/m_defaultRecompressionIndex,
// m_defaultShearModulus );
//  GEOSX_THROW_IF( poisson < 0, "Elastic parameters lead to negative Poisson ratio at reference pressure", InputError );

  // set results as array default values

  getWrapper< array2d< real64 > >( viewKeyStruct::oldPreConsolidationPressureString() ).
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getWrapper< array2d< real64 > >( viewKeyStruct::newPreConsolidationPressureString() ).
    setApplyDefaultValue( m_defaultPreConsolidationPressure );

  getWrapper< array1d< real64 > >( viewKeyStruct::virginCompressionIndexString() ).
    setApplyDefaultValue( m_defaultVirginCompressionIndex );

  getWrapper< array1d< real64 > >( viewKeyStruct::cslSlopeString() ).
    setApplyDefaultValue( m_defaultCslSlope );

  getWrapper< array1d< real64 > >( viewKeyStruct::shapeParameterString() ).
    setApplyDefaultValue( m_defaultShapeParameter );
}


void ModifiedCamClay::saveConvergedState() const
{
  SolidBase::saveConvergedState(); // TODO: not ideal, as we have separate loops for base and derived data

  localIndex const numE = numElem();
  localIndex const numQ = numQuad();

  arrayView2d< real64 const > newPreConsolidationPressure = m_newPreConsolidationPressure;
  arrayView2d< real64 > oldPreConsolidationPressure = m_oldPreConsolidationPressure;

  forAll< parallelDevicePolicy<> >( numE, [=] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    for( localIndex q = 0; q < numQ; ++q )
    {
      oldPreConsolidationPressure( k, q ) = newPreConsolidationPressure( k, q );
    }
  } );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ModifiedCamClay, std::string const &, Group * const )
}
} /* namespace geosx */
