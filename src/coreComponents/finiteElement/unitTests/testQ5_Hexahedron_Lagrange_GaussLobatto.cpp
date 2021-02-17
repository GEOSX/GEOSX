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
 * @file testQ5_Hexahedron_Lagrange_GaussLobatto
 */

#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include "gtest/gtest.h"

#include <chrono>

#include "finiteElement/elementFormulations/Q5_Hexahedron_Lagrange_GaussLobatto.hpp"

using namespace geosx;
using namespace finiteElement;

template< typename POLICY >
void testKernelDriver()
{
  constexpr int numNodes = 216;
  constexpr int numQuadraturePoints = 216;

  array1d< real64 > arrDetJ( numQuadraturePoints );
  array2d< real64 > arrN( numQuadraturePoints, numNodes );
  array3d< real64 > arrdNdX( numQuadraturePoints, numNodes, 3 );

  arrayView1d< real64 > const & viewDetJ = arrDetJ;
  arrayView2d< real64 > const & viewN = arrN;
  arrayView3d< real64 > const & viewdNdX = arrdNdX;

  
  forAll< POLICY >( 1,
                    [=] GEOSX_HOST_DEVICE ( localIndex const )
  {
    real64 Ntest[numQuadraturePoints][numNodes]={{0.0}};
    for( localIndex q=0; q<numQuadraturePoints; ++q )
    {
      real64 N[numNodes] = {0};
      Ntest[q][q]=1.0;
      Q5_Hexahedron_Lagrange_GaussLobatto::calcN( q, N );
      for( localIndex a=0; a<numNodes; ++a )
      {
        if (abs(N[a])<1e-9){
           viewN(q,a)=0;
        }
        else {
           viewN( q, a ) = N[a];
        }
        EXPECT_FLOAT_EQ( Ntest[q][a], viewN[q][a] );
      }
    }
  } );



  real64 xCoords[numNodes][3];
  xCoords[0][0]=-1.0;
  xCoords[1][0]=-sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[2][0]=-sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[3][0]=sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[4][0]=sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[5][0]=1.0;

  xCoords[0][1]=-1.0;
  xCoords[6][1]=-sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[12][1]=-sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[18][1]=sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[24][1]=sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[30][1]=1.0;

  xCoords[0][2]=-1.0;
  xCoords[36][2]=-sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[72][2]=-sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[108][2]=sqrt(1.0/21.0*(7.0-2.0*sqrt(7.0)));
  xCoords[144][2]=sqrt(1.0/21.0*(7.0+2.0*sqrt(7.0)));
  xCoords[180][2]=1.0;

  for( localIndex k=0; k<6; ++k )
  {
      for( localIndex j=0; j<6; ++j )
      {
          for( localIndex i=0; i<6; ++i )
          {
              xCoords[i+6*j+36*k][0]=xCoords[i][0];
              xCoords[i+6*j+36*k][1]=xCoords[6*j][1];
              xCoords[i+6*j+36*k][2]=xCoords[36*k][2];
          }
      }
  }
  
  real64 gradNxtest[numNodes][numQuadraturePoints] = {{0}};
  real64 gradNytest[numNodes][numQuadraturePoints] = {{0}};
  real64 gradNztest[numNodes][numQuadraturePoints] = {{0}};

  gradNxtest[0][0]=-15.0/2.0;
  gradNxtest[0][1]=1.0/6.0*(-2.0-sqrt(7.0)-sqrt(21.0+6.0*sqrt(7.0)));
  gradNxtest[0][2]=1.0/6.0*(-2.0+sqrt(7.0)+sqrt(21.0-6.0*sqrt(7.0)));
  gradNxtest[0][3]=1.0/6.0*(-2.0+sqrt(7.0)-sqrt(21.0-6.0*sqrt(7.0)));
  gradNxtest[0][4]=-1.0/6.0*(2.0+sqrt(7.0)-sqrt(21.0+6.0*sqrt(7.0)));
  gradNxtest[0][5]= -1.0/2.0;

  gradNxtest[1][0]=1.0/4.0*(7.0+4.0*sqrt(7.0)+sqrt(343.0+70.0*sqrt(7.0)));
  gradNxtest[1][2]=-1.0/24.0*sqrt(28.0-2.0*sqrt(7.0))*(5.0+sqrt(3.0)-sqrt(7.0)+sqrt(21.0));
  gradNxtest[1][3]=1.0/24.0*sqrt(28.0-2.0*sqrt(7.0))*(-5.0+sqrt(3.0)+sqrt(7.0)+sqrt(21.0));
  gradNxtest[1][4]=-1.0/2.0*sqrt(7-2.0*sqrt(7.0));
  gradNxtest[1][5]=1.0/4.0*(-7.0-4.0*sqrt(7.0)+sqrt(343.0+70.0*sqrt(7.0)));

  gradNxtest[2][0]=1.0/4.0*(7.0-4.0*sqrt(7.0)-sqrt(343.0-70.0*sqrt(7.0)));
  gradNxtest[2][1]=1.0/12.0*sqrt(7.0+sqrt(7.0)/2.0)*(5.0-sqrt(3.0)+sqrt(7.0)+sqrt(21.0));
  gradNxtest[2][3]=-1.0/2.0*sqrt(7.0+2.0*sqrt(7.0));
  gradNxtest[2][4]=1.0/12.0*sqrt(7.0+sqrt(7.0)/2.0)*(5.0+sqrt(3.0)+sqrt(7.0)-sqrt(21.0));
  gradNxtest[2][5]=1.0/4.0*(-7.0+4.0*sqrt(7.0)-sqrt(343.0-70.0*sqrt(7.0)));

  gradNxtest[3][0]=-gradNxtest[2][5];
  gradNxtest[3][1]=-gradNxtest[2][4];
  gradNxtest[3][2]=-gradNxtest[2][3];
  gradNxtest[3][4]=-gradNxtest[2][1];
  gradNxtest[3][5]=-gradNxtest[2][0];
 
  gradNxtest[4][0]=-gradNxtest[1][5];
  gradNxtest[4][1]=-gradNxtest[1][4];
  gradNxtest[4][2]=-gradNxtest[1][3];
  gradNxtest[4][3]=-gradNxtest[1][2];
  gradNxtest[4][5]=-gradNxtest[1][0];

  gradNxtest[5][0]=-gradNxtest[0][5];
  gradNxtest[5][1]=-gradNxtest[0][4];
  gradNxtest[5][2]=-gradNxtest[0][3];
  gradNxtest[5][3]=-gradNxtest[0][2];
  gradNxtest[5][4]=-gradNxtest[0][1];
  gradNxtest[5][5]=-gradNxtest[0][0];

  for( localIndex k=0; k<6; ++k )
  {
      for( localIndex j=0; j<6; ++j )
      {
          for( localIndex i=0; i<6; ++i )
          {   
              for (localIndex l=0; l<6; ++l )
              {
                gradNxtest[i+6*j+36*k][l+6*j+36*k]=gradNxtest[i][l];
                gradNxtest[i+6*j+36*k][l+6*j+36*k]=gradNxtest[i][l];

                gradNytest[i+6*j+36*k][i+6*l+36*k]=gradNxtest[j][l];
                gradNytest[i+6*j+36*k][i+6*l+36*k]=gradNxtest[j][l];

                gradNztest[i+6*j+36*k][i+6*j+36*l]=gradNxtest[k][l];
                gradNztest[i+6*j+36*k][i+6*j+36*l]=gradNxtest[k][l];
              }
          }
      }
  }

    forAll< POLICY >( 1,
                    [=] GEOSX_HOST_DEVICE ( localIndex const )                   
  { 
  
  for( localIndex q=0; q<numQuadraturePoints; ++q )
  { 
    
    real64 dNdX[numNodes][3] = {{0}};
    
    viewDetJ[q] = Q5_Hexahedron_Lagrange_GaussLobatto::calcGradN( q,
                                                                  xCoords,
                                                                  dNdX );

    for( localIndex a=0; a<numNodes; ++a )
    {         
       for( int i = 0; i < 3; ++i )
          {
            if (abs(dNdX[a][i])<1e-9)
            {
              viewdNdX(q,a,i)=0;
            }
            else 
            {
              viewdNdX( q, a,i ) = dNdX[a][i];
            }
          }
      EXPECT_FLOAT_EQ( gradNxtest[a][q], viewdNdX(q, a, 0) );
      EXPECT_FLOAT_EQ( gradNytest[a][q], viewdNdX(q, a, 1) );
      EXPECT_FLOAT_EQ( gradNztest[a][q], viewdNdX(q, a, 2) );
    } 
  }
 } );

}  

#ifdef USE_CUDA
TEST( FiniteElementShapeFunctions, testKernelCuda )
{
  testKernelDriver< geosx::parallelDevicePolicy< 32 > >();
}
#endif
TEST( FiniteElementShapeFunctions, testKernelHost )
{
  testKernelDriver< serialPolicy >();
}



using namespace geosx;
int main( int argc, char * argv[] )
{
  testing::InitGoogleTest();

  basicSetup( argc, argv, false );

  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
