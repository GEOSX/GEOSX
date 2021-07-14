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
#include "constitutive/solid/ModifiedCamClay.hpp"
#include "constitutive/solid/InvariantDecompositions.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;


struct StrainData
{
  real64 strainIncrement[6] = { };
};


template< typename CMW >
void getStress( CMW const cmw,
                real64 (& stress)[6] )
{
  forAll< serialPolicy >( 1, [&stress, cmw] ( localIndex const k )
  {
    stress[0] = cmw.m_newStress( k, 0, 0 );
    stress[1] = cmw.m_newStress( k, 0, 1 );
    stress[2] = cmw.m_newStress( k, 0, 2 );
    stress[3] = cmw.m_newStress( k, 0, 3 );
    stress[4] = cmw.m_newStress( k, 0, 4 );
    stress[5] = cmw.m_newStress( k, 0, 5 );
  } );
}


template< typename CMW >
void setStress( CMW cmw,
                real64 const (&stress)[6] )
{
  forAll< serialPolicy >( 1, [&stress, cmw] ( localIndex const k )
  {
    cmw.m_oldStress( k, 0, 0 ) = stress[0];
    cmw.m_oldStress( k, 0, 1 ) = stress[1];
    cmw.m_oldStress( k, 0, 2 ) = stress[2];
    cmw.m_oldStress( k, 0, 3 ) = stress[3];
    cmw.m_oldStress( k, 0, 4 ) = stress[4];
    cmw.m_oldStress( k, 0, 5 ) = stress[5];
  } );
}


template< typename POLICY >
void testModifiedCamClayDriver()
{
  // create a Cam-Clay model, and test xml input
  conduit::Node node;
  dataRepository::Group rootGroup( "root", node );
  ConstitutiveManager constitutiveManager( "constitutive", &rootGroup );

  string const inputStream =
    "<Constitutive>"
    "   <ModifiedCamClay"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultRefPressure=\"-1.0\" "
    "      defaultRefStrainVol=\"0.0\" "
    "      defaultShearModulus=\"200.0\" "
    "      defaultPreConsolidationPressure =\"-1.5\" "
    "      defaultCslSlope=\"1.0\" "
    "      defaultVirginCompressionIndex=\"0.003\" "
    "      defaultRecompressionIndex=\"0.002\"/>"
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

  ModifiedCamClay & cm = constitutiveManager.getConstitutiveRelation< ModifiedCamClay >( "granite" );
  cm.allocateConstitutiveData( disc, numQuad );

  // confirm allocation sizes

  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );

  // use updates class to run a uniaxial compression test.  the material
  // begins with some cohesion, but it quickly degrades to zero under
  // plastic loading.  at the end of the loading, we confirm the
  // stress point lies on the correctly-positioned yield surface.

  ModifiedCamClay::KernelWrapper cmw = cm.createKernelUpdates();

  StrainData data;
  data.strainIncrement[0] = -1e-4;

  // set initial stress

  real64 stress[6] = {0};
  stress[0] = -1.0;
  stress[1] = -1.0;
  stress[2] = -1.0;

  setStress( cmw, stress );

  // run loading

  for( localIndex loadstep=0; loadstep < 500; ++loadstep )
  {
    forAll< parallelDevicePolicy<> >( 1, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      real64 stressLocal[6] = {0};
      real64 stiffnessLocal[6][6] = {{0}};
      cmw.smallStrainUpdate( k, 0, data.strainIncrement, stressLocal, stiffnessLocal );
      //    std::cout<< stressLocal[0] <<std::endl;
    } );
    cm.saveConvergedState();
  }

  // get final stress state in p-q space

  getStress( cmw, stress );

  real64 invariantP, invariantQ;
  real64 deviator[6];

  twoInvariant::stressDecomposition( stress,
                                     invariantP,
                                     invariantQ,
                                     deviator );

  //EXPECT_TRUE( ... to be determined ... );

  // we now use a finite-difference check of tangent stiffness to confirm
  // the analytical form is working properly.

  cmw.checkSmallStrainStiffness( 0, 0, data.strainIncrement, true );

  EXPECT_TRUE( cmw.checkSmallStrainStiffness( 0, 0, data.strainIncrement ) );
}


#ifdef USE_CUDA
TEST( ModifiedCamClayTests, testModifiedCamClayDevice )
{
  testModifiedCamClayDriver< geosx::parallelDevicePolicy< > >();
}
#endif
TEST( ModifiedCamClayTests, testModifiedCamClayHost )
{
  testModifiedCamClayDriver< serialPolicy >();
}
