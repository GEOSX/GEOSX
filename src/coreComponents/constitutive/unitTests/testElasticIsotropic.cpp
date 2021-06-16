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

// Source includes
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "dataRepository/xmlWrapper.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( ElasticIsotropicTests, testAllocation )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ElasticIsotropic cm( "model", &rootGroup );

  localIndex constexpr numElems = 2;
  localIndex constexpr numQuadraturePoints = 3;

  dataRepository::Group disc( "discretization", &rootGroup );
  disc.resize( numElems );
  cm.allocateConstitutiveData( disc, numQuadraturePoints );

  EXPECT_EQ( cm.size(), numElems );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuadraturePoints );

  arrayView1d< real64 const > const bulkModulus = cm.bulkModulus();
  arrayView1d< real64 const > const shearModulus = cm.shearModulus();
  arrayView3d< real64 const, solid::STRESS_USD > const stress = cm.getStress();

  EXPECT_EQ( bulkModulus.size(), numElems );
  EXPECT_EQ( shearModulus.size(), numElems );
  EXPECT_EQ( stress.size( 0 ), numElems );
  EXPECT_EQ( stress.size( 1 ), numQuadraturePoints );
  EXPECT_EQ( stress.size( 2 ), 6 );
}

TEST( ElasticIsotropicTests, testStateUpdatePoint )
{
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  real64 constexpr K = 2e10;
  real64 constexpr G = 1e10;

  string const inputStream =
    "<Constitutive>"
    "   <ElasticIsotropic"
    "      name=\"granite\" "
    "      defaultReferencePorosity=\"0.2\" "
    "      grainBulkModulus=\"1e15\" "
    "      defaultGrainDensity=\"2700\" "
    "      defaultBulkModulus=\"" + std::to_string( K ) + "\" "
                                                          "      defaultShearModulus=\"" + std::to_string( G ) + "\"/>"
                                                                                                                 "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.postProcessInputRecursive();

  dataRepository::Group disc( "discretization", &rootGroup );
  disc.resize( 2 );

  ElasticIsotropic & cm = constitutiveManager.getConstitutiveRelation< ElasticIsotropic >( "granite" );

  cm.allocateConstitutiveData( disc, 2 );
  ElasticIsotropic::KernelWrapper cmw = cm.createKernelUpdates();

  arrayView3d< real64, solid::STRESS_USD > const & stress = cm.getStress();

  real64 const strain = 0.1;
  real64 Ddt[ 6 ] = { 0 };
  real64 Rot[ 3 ][ 3 ] = { { 0 } };

  real64 pointStress[6];
  {
    Ddt[ 0 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), (2.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), 0 );
  }

  {
    stress.zero();
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[ 1 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), (2.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), 0 );
  }

  {
    stress.zero();
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[ 2 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), (-1.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), (2.0/3.0*strain)*2*G + strain*K );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), 0 );
  }

  {
    stress.zero();
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[ 5 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), strain*G );
  }

  {
    stress.zero();
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[ 4 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), strain*G );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), 0 );
  }

  {
    stress.zero();
    LvArray::tensorOps::fill< 6 >( Ddt, 0 );

    Ddt[ 3 ] = strain;
    Rot[ 0 ][ 0 ] = 1;
    Rot[ 1 ][ 1 ] = 1;
    Rot[ 2 ][ 2 ] = 1;

    cmw.hypoUpdate_StressOnly( 0, 0, Ddt, Rot, pointStress );

    EXPECT_DOUBLE_EQ( stress( 0, 0, 0 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 1 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 2 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 3 ), strain*G );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 4 ), 0 );
    EXPECT_DOUBLE_EQ( stress( 0, 0, 5 ), 0 );
  }
}
