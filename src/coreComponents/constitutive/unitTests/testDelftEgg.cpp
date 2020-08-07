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
#include "constitutive/solid/DelftEgg.hpp"

#include "dataRepository/xmlWrapper.hpp"

using namespace geosx;
using namespace ::geosx::constitutive;


// NOTE: using this for debugging, will set up proper unit tests later

TEST( DelftEggTests, testModel )
{
  ConstitutiveManager constitutiveManager( "constitutive", nullptr );

  string const inputStream =
    "<Constitutive>"
    "   <DelftEgg"
    "      name=\"granite\" "
    "      defaultDensity=\"2700\" "
    "      defaultRefPressure=\"-3.0\" "
    "      defaultRefStrainVol=\"-1e-4\" "
    "      defaultShearModulus=\"100.0\" "
    "      defaultPreConsolidationPressure =\"-4.0\" "
    "      defaultShapeParameter=\"1.0\" "
    "      defaultCslSlope=\"1.2\" "
    "      defaultVirginCompressionIndex=\"0.01\" "
    "      defaultRecompressionIndex=\"0.1\"/>"
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
  constitutiveManager.ProcessInputFileRecursive( xmlConstitutiveNode );
  constitutiveManager.PostProcessInputRecursive();
  
  localIndex constexpr numElem = 2;
  localIndex constexpr numQuad = 4;
  
  dataRepository::Group disc( "discretization", nullptr );
  disc.resize( numElem );
  
  DelftEgg & cm = *(constitutiveManager.GetConstitutiveRelation<DelftEgg>("granite"));
  cm.AllocateConstitutiveData( &disc, numQuad );
  
  EXPECT_EQ( cm.size(), numElem );
  EXPECT_EQ( cm.numQuadraturePoints(), numQuad );
  
  DelftEgg::KernelWrapper cmw = cm.createKernelWrapper();
  
  real64 inc = -1e-4; // compression
  real64 total = 0;
  
  //array2d< real64 > strainIncrement(1,6);
  //                  strainIncrement.setValues< serialPolicy >( 0 );
  //                  strainIncrement[0][0] = inc;

  real64 strainIncrement[6] = {inc, 0, 0, 0, 0, 0};
  real64 stress[6] = {0, 0, 0, 0, 0, 0};
  real64 stiffness[6][6];

  //array2d< real64 > stress(1,6);
  //array3d< real64 > stiffness(1,6,6);
  
  for(localIndex loadstep=0; loadstep < 100; ++loadstep)
  {
    cmw.SmallStrainUpdate(0,0,strainIncrement,stress,stiffness);
    cmw.SaveConvergedState(0,0);
    
    total += inc;
    
    real64 mean = (stress[0]+stress[1]+stress[2])/3;
    real64 deviator = stress[0]-stress[1];
    
    std::cout << mean << " " << deviator << " " << total << std::endl;
  }
  
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      std::cout << stiffness[i][j] << " ";
    }
    std::cout << "\n";
  }
  
  // finite-difference check of tangent stiffness
  
  array2d< real64 > fd_stiffness(6,6);
  //array2d< real64 > pstress(1,6);
  real64 pstress[6] = {0, 0, 0, 0, 0, 0};
  real64 pstiffness[6][6];

  //array3d< real64 > pstiffness(1,6,6);
  
  real64 eps = 1e-12;
  
  cmw.SmallStrainUpdate(0,0,strainIncrement,stress,stiffness);
  
  for(localIndex i=0; i<6; ++i)
  {
    strainIncrement[i] += eps;
    
    if(i>0)
    {
      strainIncrement[i-1] -= eps;
    }
    
    cmw.SmallStrainUpdate(0,0,strainIncrement,pstress,pstiffness);
    
    for(localIndex j=0; j<6; ++j)
    {
      fd_stiffness[j][i] = (pstress[j]-stress[j])/eps;
    }
  }
  
  for(localIndex i=0; i<6; ++i)
  {
    for(localIndex j=0; j<6; ++j)
    {
      std::cout << fd_stiffness[i][j] << " ";
    }
    std::cout << "\n";
  }
  
}



