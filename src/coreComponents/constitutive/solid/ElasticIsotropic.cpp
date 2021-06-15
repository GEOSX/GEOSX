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
 *  @file ElasticIsotropic.cpp
 */

#include "ElasticIsotropic.hpp"

namespace geosx
{
using namespace dataRepository;
namespace constitutive
{

ElasticIsotropic::ElasticIsotropic( string const & name, Group * const parent ):
  SolidBase( name, parent ),
  m_defaultBulkModulus(),
  m_defaultShearModulus(),
  m_bulkModulus(),
  m_shearModulus()
{
  registerWrapper( viewKeyStruct::defaultBulkModulusString(), &m_defaultBulkModulus ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Bulk Modulus Parameter" );

  registerWrapper( viewKeyStruct::defaultShearModulusString(), &m_defaultShearModulus ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Shear Modulus Parameter" );

  registerWrapper< real64 >( viewKeyStruct::defaultYoungsModulusString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Elastic Young's Modulus." );

  registerWrapper< real64 >( viewKeyStruct::defaultPoissonRatioString() ).
    setApplyDefaultValue( -1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Poisson's ratio" );

  registerWrapper( viewKeyStruct::bulkModulusString(), &m_bulkModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Bulk Modulus Field" );

  registerWrapper( viewKeyStruct::shearModulusString(), &m_shearModulus ).
    setApplyDefaultValue( -1 ).
    setDescription( "Elastic Shear Modulus" );
}


ElasticIsotropic::~ElasticIsotropic()
{}


void ElasticIsotropic::postProcessInput()
{
  // check what constants the user actually input, and do conversions as needed

  SolidBase::postProcessInput();

  real64 & nu = getReference< real64 >( viewKeyStruct::defaultPoissonRatioString() );
  real64 & E  = getReference< real64 >( viewKeyStruct::defaultYoungsModulusString() );
  real64 & K  = m_defaultBulkModulus;
  real64 & G  = m_defaultShearModulus;

  string errorCheck( "( " );
  int numConstantsSpecified = 0;
  if( nu >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "nu, ";
  }
  if( E >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "E, ";
  }
  if( K >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "K, ";
  }
  if( G >= 0.0 )
  {
    ++numConstantsSpecified;
    errorCheck += "G, ";
  }
  errorCheck += ")";

  GEOSX_ERROR_IF( numConstantsSpecified != 2,
                  "A specific pair of elastic constants is required. Either (K,G) or (E,nu). "<<
                  "You have specified "<<errorCheck );

  if( nu >= 0.0 && E >= 0.0 )
  {
    K = BulkModulus().
          setYoungModulus( E ).
          setPoissonRatio( nu ).
          getValue();

    G = ShearModulus().
          setYoungModulus( E ).
          setPoissonRatio( nu ).
          getValue();
  }
  else if( nu >= 0.0 && G >= 0.0 )
  {
    E = YoungModulus().
          setShearModulus( G ).
          setPoissonRatio( nu ).
          getValue();

    K = BulkModulus().
          setShearModulus( G ).
          setPoissonRatio( nu ).
          getValue();
  }
  else if( nu >= 0 && K >= 0.0 )
  {
    E = YoungModulus().
          setBulkModulus( K ).
          setPoissonRatio( nu ).
          getValue();

    G = ShearModulus().
          setBulkModulus( K ).
          setPoissonRatio( nu ).
          getValue();
  }
  else if( E >= 0.0 && K >=0 )
  {
    nu = PoissonRatio( BulkModulus( K ), YoungModulus( E ) ).value;

    G = ShearModulus().
          setYoungModulus( E ).
          setBulkModulus( K ).
          getValue();
  }
  else if( E >= 0.0 && G >= 0 )
  {
    nu = PoissonRatio( ShearModulus( G ), YoungModulus( E ) ).value;

    K = BulkModulus().
          setYoungModulus( E ).
          setShearModulus( G ).
          getValue();
  }
  else if( K >= 0.0 && G >= 0.0 )
  {

    E = YoungModulus().
          setBulkModulus( K ).
          setShearModulus( G ).
          getValue();
    nu = PoissonRatio( BulkModulus( K ), ShearModulus( G ) ).value;
  }
  else
  {
    GEOSX_ERROR( "invalid specification for default elastic constants. "<<errorCheck<<" has been specified." );
  }

  // set results as array default values
  this->getWrapper< array1d< real64 > >( viewKeyStruct::bulkModulusString() ).
    setApplyDefaultValue( m_defaultBulkModulus );

  this->getWrapper< array1d< real64 > >( viewKeyStruct::shearModulusString() ).
    setApplyDefaultValue( m_defaultShearModulus );
}

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ElasticIsotropic, string const &, Group * const )
}
} /* namespace geosx */
