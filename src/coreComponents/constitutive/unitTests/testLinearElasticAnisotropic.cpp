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

#include "gtest/gtest.h"

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/LinearElasticAnisotropic.hpp"

#include "dataRepository/xmlWrapper.hpp"
using namespace geosx;
using namespace ::geosx::constitutive;

TEST( LinearElasticAnisotropicTests, testAllocation )
{
  LinearElasticAnisotropic cm( "model", nullptr );

  localIndex constexpr numElems = 2;
  localIndex constexpr numQuadraturePoints = 3;

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElems );
  cm.AllocateConstitutiveData( &disc, numQuadraturePoints );

  EXPECT_EQ( cm.size(), numElems );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuadraturePoints );

  //  arrayView1d<LinearElasticAnisotropic::StiffnessTensor const> const &
  //  stiffness = cm.stiffness() ;

  arrayView3d< real64 const, solid::STRESS_USD > const & stress = cm.getStress();

  //  EXPECT_EQ( stiffness.size(), numElems );
  EXPECT_EQ( stress.size( 0 ), numElems );
  EXPECT_EQ( stress.size( 1 ), numQuadraturePoints );
  EXPECT_EQ( stress.size( 2 ), 6 );
}

void
stressCalc( real64 const c[6][6],
            real64 const ( &Ddt )[6],
            real64 stressVoigt[6] )
{
  real64 const DdtVoigt[6] =
    { Ddt[0], Ddt[1], Ddt[2], 2 * Ddt[3], 2 * Ddt[4], 2 * Ddt[5] };

  for( int i = 0; i < 6; ++i )
  {
    stressVoigt[i] = 0.0;
    for( int j = 0; j < 6; ++j )
    {
      stressVoigt[i] += c[i][j] * DdtVoigt[j];
    }
  }
}

void
voigtStrain( real64 strainVoigt[6], real64 const Ddt[6] )
{
  strainVoigt[0] = Ddt[0];
  strainVoigt[1] = Ddt[1];
  strainVoigt[2] = Ddt[2];
  strainVoigt[3] = Ddt[3] * 2;
  strainVoigt[4] = Ddt[4] * 2;
  strainVoigt[5] = Ddt[5] * 2;
}

void
stressSliceCheck( arrayView3d< real64 const, solid::STRESS_USD > const & stress,
                  real64 const stressV[6] )
{
  ASSERT_DOUBLE_EQ( stress( 0, 0, 0 ), stressV[0] );
  ASSERT_DOUBLE_EQ( stress( 0, 0, 1 ), stressV[1] );
  ASSERT_DOUBLE_EQ( stress( 0, 0, 2 ), stressV[2] );
  ASSERT_DOUBLE_EQ( stress( 0, 0, 3 ), stressV[3] );
  ASSERT_DOUBLE_EQ( stress( 0, 0, 4 ), stressV[4] );
  ASSERT_DOUBLE_EQ( stress( 0, 0, 5 ), stressV[5] );
}

void
stressCheck( real64 const stressV[6], real64 const stressV2[6] )
{
  ASSERT_DOUBLE_EQ( stressV2[0], stressV[0] );
  ASSERT_DOUBLE_EQ( stressV2[1], stressV[1] );
  ASSERT_DOUBLE_EQ( stressV2[2], stressV[2] );
  ASSERT_DOUBLE_EQ( stressV2[3], stressV[3] );
  ASSERT_DOUBLE_EQ( stressV2[4], stressV[4] );
  ASSERT_DOUBLE_EQ( stressV2[5], stressV[5] );
}

TEST( LinearElasticAnisotropicTests, testStateUpdatePoint )
{
  LinearElasticAnisotropic cm( "model", nullptr );

  real64 c[6][6] = { { 1.0e11, 0.1e10, 0.2e10, 0.3e10, 0.4e10, 0.5e10 },
                     { 1.0e10, 1.1e11, 1.2e10, 1.3e10, 1.4e10, 1.5e10 },
                     { 2.0e10, 2.1e10, 2.2e11, 2.3e10, 2.4e10, 2.5e10 },
                     { 3.0e10, 3.1e10, 3.2e10, 3.3e10, 3.4e10, 3.5e10 },
                     { 4.0e10, 4.1e10, 4.2e10, 4.3e10, 4.4e10, 4.5e10 },
                     { 5.0e10, 5.1e10, 5.2e10, 5.3e10, 5.4e10, 5.5e10 } };

  cm.setDefaultStiffness( c );

  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( 2 );
  cm.AllocateConstitutiveData( &disc, 2 );

  auto cw = cm.createKernelUpdates();

  auto const & stateStress = cw.m_stress;

  real64 const strain = 0.1;
  real64 Ddt[6] = { 0 };
  real64 strainV[6] = { 0 };
  real64 stressV[6] = { 0 };
  real64 stressV2[6] = { 0 };
  real64 Rot[3][3] = { { 0 } };

  {
    Ddt[0] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );

    stateStress.setValues< serialPolicy >( 0 );
    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressSliceCheck( stateStress, stressV );

    voigtStrain( strainV, Ddt );
    stateStress.setValues< serialPolicy >( 0 );
    cw.SmallStrain( 0, 0, strainV );
    stressSliceCheck( stateStress, stressV );

    stateStress.setValues< serialPolicy >( 0 );
    cw.SmallStrainNoState( 0, strainV, stressV2 );
    stressCheck( stressV, stressV2 );
  }

  {
    stateStress.setValues< serialPolicy >( 0 );
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[1] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );
  }

  {
    stateStress.setValues< serialPolicy >( 0 );
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[2] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );
  }

  {
    stateStress.setValues< serialPolicy >( 0 );
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[5] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );
  }

  {
    stateStress.setValues< serialPolicy >( 0 );
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[4] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );
  }

  {
    stateStress.setValues< serialPolicy >( 0 );
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[3] = strain;
    Rot[0][0] = 1;
    Rot[1][1] = 1;
    Rot[2][2] = 1;

    cw.HypoElastic( 0, 0, Ddt, Rot );
    stressCalc( c, Ddt, stressV );
    stressSliceCheck( stateStress, stressV );
  }
}

TEST( LinearElasticAnisotropicTests, testXML )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  string const inputStream =
    "<Constitutive>"
    "  <LinearElasticAnisotropic name=\"granite\""
    "                            defaultDensity=\"2700\""
    "                            defaultStiffness=\"{ { 1.0e10,  1.1e9,  "
    "1.2e9, 1.3e9, 1.4e9, 1.5e9},"
    "                                                 {  2.0e9, 2.1e10,  "
    "2.2e9, 2.3e9, 2.4e9, 2.5e9},"
    "                                                 {  3.0e9,  3.1e9, "
    "3.2e10, 3.3e9, 3.4e9, 3.5e9},"
    "                                                 {  4.0e9,  4.1e9,  "
    "4.2e9, 4.3e9, 4.4e9, 4.5e9},"
    "                                                 {  5.0e9,  5.1e9,  "
    "5.2e9, 5.3e9, 5.4e9, 5.5e9},"
    "                                                 {  6.0e9,  6.1e9,  "
    "6.2e9, 6.3e9, 6.4e9, 6.5e9} }\" />"
    "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult =
    xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description() );
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();

  LinearElasticAnisotropic * const model =
    constitutiveManager.GetConstitutiveRelation< LinearElasticAnisotropic >(
      "granite" );
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( 1 );
  model->AllocateConstitutiveData( &disc, 1 );

  LinearElasticAnisotropicUpdates kernelWrapper = model->createKernelUpdates();

  real64 stiffness[6][6];
  kernelWrapper.GetStiffness( 0, stiffness );

  real64 c[6][6] = { { 1.0e10, 1.1e9, 1.2e9, 1.3e9, 1.4e9, 1.5e9 },
                     { 2.0e9, 2.1e10, 2.2e9, 2.3e9, 2.4e9, 2.5e9 },
                     { 3.0e9, 3.1e9, 3.2e10, 3.3e9, 3.4e9, 3.5e9 },
                     { 4.0e9, 4.1e9, 4.2e9, 4.3e9, 4.4e9, 4.5e9 },
                     { 5.0e9, 5.1e9, 5.2e9, 5.3e9, 5.4e9, 5.5e9 },
                     { 6.0e9, 6.1e9, 6.2e9, 6.3e9, 6.4e9, 6.5e9 } };

  for( int i = 0; i < 6; ++i )
  {
    for( int j = 0; j < 6; ++j )
    {
      ASSERT_DOUBLE_EQ( stiffness[i][j], c[i][j] );
    }
  }
}
