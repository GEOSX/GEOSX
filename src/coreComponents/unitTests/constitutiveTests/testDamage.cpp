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

#include <gtest/gtest.h>

#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/ElasticIsotropic.hpp"
#include "constitutive/solid/DamageSpectral.hpp"

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;

TEST( DamageTests, testDamageSpectral )
{
  // create a model and test xml input

  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  string const inputStream =
    "<Constitutive>"
    "   <DamageSpectralElasticIsotropic"
    "      name=\"shale\" "
    "      defaultDensity=\"2700\" "
    "      defaultBulkModulus=\"1.7e5\" "
    "      defaultShearModulus=\"8.0e4\" "
    "      lengthScale=\"0.25\" "
    "      criticalFractureEnergy=\"54.0\" "
    "      criticalStrainEnergy=\"15.0\" "
    "   />"
    "</Constitutive>";

  xmlWrapper::xmlDocument xmlDocument;
  xmlWrapper::xmlResult xmlResult = xmlDocument.load_buffer( inputStream.c_str(),
                                                             inputStream.size() );
  if( !xmlResult )
  {
    GEOSX_LOG_RANK_0( "XML parsed with errors!" );
    GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description());
    GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
  }

  xmlWrapper::xmlNode xmlConstitutiveNode = xmlDocument.child( "Constitutive" );
  constitutiveManager.processInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.postProcessInputRecursive();

  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;

  dataRepository::Group disc( "discretization", &rootGroup );
  disc.resize( numElem );

  DamageSpectral< ElasticIsotropic > & cm = constitutiveManager.getConstitutiveRelation< DamageSpectral< ElasticIsotropic > >( "shale" );
  cm.allocateConstitutiveData( disc, numQuad );

  // confirm allocation sizes

  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );

  // use updates class to run a uniaxial load test

  DamageSpectral< ElasticIsotropic >::KernelWrapper cmw = cm.createKernelUpdates();

  real64 strainIncrement[6] = {1e-4, 0, 0, 0, 0, 0};
  real64 stress[6];
  real64 stiffness[6][6];

  for( localIndex loadstep=0; loadstep < 50; ++loadstep )
  {
    cmw.smallStrainUpdate( 0, 0, strainIncrement, stress, stiffness );
    cm.saveConvergedState();
  }

  // we now use a finite-difference check of tangent stiffness to confirm
  // the analytical form is working properly.
  //cmw.checkSmallStrainStiffness( 0, 0, strainIncrement, true );
  EXPECT_TRUE( cmw.checkSmallStrainStiffness( 0, 0, strainIncrement ) );
}
