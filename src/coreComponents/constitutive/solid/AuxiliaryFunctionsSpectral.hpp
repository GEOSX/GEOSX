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
 * @file AuxiliaryFunctionsSpectral.hpp
 * @brief Helper functions to perform spectral decomposition of stresses.
 *
 * A detailed description of the calculations performed here can be found in
 * Jiang, Wen et al. "Three-dimensional phase-field modeling of porosity dependent intergranular fracture in UO2"
 * Computational Materials Science 171 (2020): 109269
 *
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_
#include <algorithm>
#include <cmath>
#include <iostream>
#include "LvArray/src/output.hpp"
#include "LvArray/src/tensorOps.hpp"
namespace geosx
{

/**
 *@brief Get the positive part only of a symmetric, 3x3, tensor T using spectral split. The eigenvectors and eigenvalues of T must be passed in.
 *@param[in] eigs Array of eigenvalues of T.
 *@param[in] eigevecs 3x3 array with the eigenvectors of T as rows.
 *@param[out] positivePart Array that stores the positive part of T in Voigt Notation.
 */
GEOSX_HOST_DEVICE inline
void PositivePartOfTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& positivePart)[6] )
{
  real64 positiveEigs[6]={};
  for( int i=0; i < 3; i++ )
  {
    positiveEigs[i] = fmax( 0.0, eigs[i] );
  }
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( positivePart, eigvecs, positiveEigs );
}

/**
 *@brief Get the negative part only of a symmetric, 3x3, tensor T using spectral split. The eigenvectors and eigenvalues of T must be passed in.
 *@param[in] eigs Array of eigenvalues of T.
 *@param[in] eigevecs 3x3 array with the eigenvectors of T as rows.
 *@param[out] negativePart Array that stores the negative part of T in Voigt Notation.
 */
GEOSX_HOST_DEVICE inline
void NegativePartOfTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& negativePart)[6] )
{
  real64 negativeEigs[6]={};
  for( int i=0; i < 3; i++ )
  {
    negativeEigs[i] = fmin( 0.0, eigs[i] );
  }
  LvArray::tensorOps::Rij_eq_AikSymBklAjl< 3 >( negativePart, eigvecs, negativeEigs );
}

/**
 *@brief Implements the : (double contraction) operator between two symmetric, 3x3, second-order tensors in voigt form.
 *@param[in] A First operand in Voigt notation.
 *@param[in] B Second operand in Voigt notation.
 *@return Result of the operation.
 */
GEOSX_HOST_DEVICE inline
real64 doubleContraction( real64 (& A)[6], real64 (& B)[6] )
{
  real64 ans = 0;
  for( int i=0; i < 6; i++ )
  {
    if( i < 3 )
    {
      ans = ans + A[i]*B[i];
    }
    else
    {
      ans = ans + 2*A[i]*B[i];
    }
  }
  return ans;
}

/**
 *@brief Compute strain from stresses using SSLE.
 *@tparam USD
 *@param[in] stress The stress tensor.
 *@param[in] K The bulk modulus.
 *@param[in] mu The shear modulus.
 *@param[out] strain The strain tensor in Voigt notaion.
 */
template< int USD >
GEOSX_HOST_DEVICE inline
void recoverStrainFromStress( arraySlice1d< real64 const, USD > const & stress,
                              real64 const K,
                              real64 const mu,
			      real64 (& strain)[6])
{
  real64 E = 9*K*mu / (3*K + mu);
  real64 nu = (3*K - 2*mu) / (6*K + 2*mu);
  strain[0] = (stress[0] - nu*(stress[1] + stress[2]))/E;
  strain[1] = (stress[1] - nu*(stress[0] + stress[2]))/E;
  strain[2] = (stress[2] - nu*(stress[0] + stress[1]))/E;
  strain[3] = (1 + nu)*stress[3]/E;
  strain[4] = (1 + nu)*stress[4]/E;
  strain[5] = (1 + nu)*stress[5]/E;
}

/**
 *@brief Heaviside function, return 1 for positive, 0 for negatives and 0.5 if argument is zero or very close.
 *@param[in] x The argument (real number).
 *@return The result.
 */
GEOSX_HOST_DEVICE inline
real64 heaviside( real64 x )
{
  if( fabs( x ) < 1e-12 )
  {
    return 0.5;
  }
  if( x > 0 )
  {
    return 1;
  }
  if( x <= 0 )
  {
    return 0;
  }
  GEOSX_ERROR( "This function was not supposed to reach this line" );
  return 1000000;
}

/**
 *@brief Computes a tensor that enters the calculation of the Jacobian of the Spectral Split - check reference paper for more details.
 *@param[in] eigvector One of the eigenvectors of the tensor being studied.
 *@param[out] Q The 4th-order tensor that results of the operation, in Voigt form.
 */
GEOSX_HOST_DEVICE inline
void QTensor( real64 const (&eigvector)[3], real64 (& Q)[6][6] )
{
  real64 M[6]={0};
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M, eigvector );
  for( int i = 0; i<6; i++ )
  {
    for( int j = 0; j<6; j++ )
    {
      Q[i][j] = M[ i ] * M[ j ];
    }
  }
}

/**
 *@brief Computes another tensor that enters the calculation of the Jacobian of the Spectral Split - check reference paper for more details.
 *@param[in] eigvec1 One of the eigenvectors of the tensor being studied.
 *@param[in] eigvec2 Another eigenvector.
 *@param[out] G The 4th-order tensor that results of the operation, in Voigt form.
 */
GEOSX_HOST_DEVICE inline
void GTensor( real64 (& eigvec1)[3], real64 (& eigvec2)[3], real64 (& G)[6][6] )
{
  GEOSX_UNUSED_VAR( eigvec1, eigvec2, G );
  real64 M1[6]={0};
  real64 M2[6]={0};
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M1, eigvec1 );
  LvArray::tensorOps::symRij_eq_AiAj< 3 >( M2, eigvec2 );

  G[0][0] = M1[0]*M2[0] + M1[0]*M2[0];
  G[0][1] = M1[5]*M2[5] + M1[5]*M2[5];
  G[0][2] = M1[4]*M2[4] + M1[4]*M2[4];
  G[0][3] = M1[5]*M2[4] + M1[4]*M2[5];
  G[0][4] = M1[0]*M2[4] + M1[4]*M2[0];
  G[0][5] = M1[0]*M2[5] + M1[5]*M2[0];
  G[1][0] = M1[5]*M2[5] + M1[5]*M2[5];
  G[1][1] = M1[1]*M2[1] + M1[1]*M2[1];
  G[1][2] = M1[3]*M2[3] + M1[3]*M2[3];
  G[1][3] = M1[1]*M2[3] + M1[3]*M2[1];
  G[1][4] = M1[5]*M2[3] + M1[3]*M2[5];
  G[1][5] = M1[5]*M2[1] + M1[1]*M2[5];
  G[2][0] = M1[4]*M2[4] + M1[4]*M2[4];
  G[2][1] = M1[3]*M2[3] + M1[3]*M2[3];
  G[2][2] = M1[2]*M2[2] + M1[2]*M2[2];
  G[2][3] = M1[3]*M2[2] + M1[2]*M2[3];
  G[2][4] = M1[4]*M2[2] + M1[2]*M2[4];
  G[2][5] = M1[4]*M2[3] + M1[3]*M2[4];
  G[3][0] = M1[5]*M2[4] + M1[5]*M2[4];
  G[3][1] = M1[1]*M2[3] + M1[1]*M2[3];
  G[3][2] = M1[3]*M2[2] + M1[3]*M2[2];
  G[3][3] = M1[1]*M2[2] + M1[3]*M2[3];
  G[3][4] = M1[5]*M2[2] + M1[3]*M2[4];
  G[3][5] = M1[5]*M2[3] + M1[1]*M2[4];
  G[4][0] = M1[0]*M2[4] + M1[0]*M2[4];
  G[4][1] = M1[5]*M2[3] + M1[5]*M2[3];
  G[4][2] = M1[4]*M2[2] + M1[4]*M2[2];
  G[4][3] = M1[5]*M2[2] + M1[4]*M2[3];
  G[4][4] = M1[0]*M2[2] + M1[4]*M2[4];
  G[4][5] = M1[0]*M2[3] + M1[5]*M2[4];
  G[5][0] = M1[0]*M2[5] + M1[0]*M2[5];
  G[5][1] = M1[5]*M2[1] + M1[5]*M2[1];
  G[5][2] = M1[4]*M2[3] + M1[4]*M2[3];
  G[5][3] = M1[5]*M2[3] + M1[4]*M2[1];
  G[5][4] = M1[0]*M2[3] + M1[4]*M2[5];
  G[5][5] = M1[0]*M2[1] + M1[5]*M2[5];
}

/**
 *@brief This function takes the eigen-decomposition of a tensor and builds the 4th order Positive Projector associated with it
 *@param[in] eigs array with the 3 eigenvalues of a tensor
 *@param[in] eigvecs 3x3 array with the 3 eigenvectors of a tensor (in rows)
 *@param[out] PositiveProjector empty array that will be populated with the voigt form of the Positive Projector
 *
 * Given a symmetric tensor T, the positive projector is defined as P+ = variation(T+)/variation(T). That is, if we
 * define the function f+ to be the positive spectral part of T, then, P+ is just the variational derivative of f+.
 * Note that we don't take the tensor T as a parameter, only its eigenvectors and eigenvalues.
 *
 */
GEOSX_HOST_DEVICE inline
void PositiveProjectorTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& PositiveProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-12;
  if( fabs( eigs[0] - eigs[1] ) < tol || fabs( eigs[0]-eigs[2] ) < tol || fabs( eigs[1]-eigs[2] ) < tol )
  {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6] = {};
  //init GVoigt
  real64 Gsym[6][6] = {};
  real64 Gji[6][6] = {};

  //compute projector
  for( int i = 0; i < 3; i++ )
  {
    real64 ithEigenVector[3] = {};
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor( ithEigenVector, Qi );
    //do update
    LvArray::tensorOps::scale< 6, 6 >( Qi, heaviside( eigs[i] ));
    LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Qi );
    if( !repeatedEigenvalues )
    {

      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3]={};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        //Do update
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (fmax( eigs[i], 0.0 ) - fmax( eigs[j], 0.0 ))/(2*(eigs[i]-eigs[j])));
        LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Gsym );
      }
    }
    else
    {
      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3]={};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (heaviside( eigs[i] ) + heaviside( eigs[j] ))/4 );
        //do update
        LvArray::tensorOps::add< 6, 6 >( PositiveProjector, Gsym );
      }
    }
  }

}

//This is the negative projector, check the documentation of the positive projector for more details.
/**
 *@brief This function takes the eigen-decomposition of a tensor and builds the 4th order Negative Projector associated with it
 *@param[in] eigs array with the 3 eigenvalues of a tensor
 *@param[in] eigvecs 3x3 array with the 3 eigenvectors of a tensor (in rows)
 *@param[out] NegativeProjector empty array that will be populated with the voigt form of the Negative Projector
 *
 * Given a symmetric tensor T, the negative projector is defined as P- = variation(T-)/variation(T). That is, if we
 * define the function f- to be the negative spectral part of T, then, P- is just the variational derivative of f-.
 * Note that we don't take the tensor T as a parameter, only its eigenvectors and eigenvalues.
 *
 */
GEOSX_HOST_DEVICE inline
void NegativeProjectorTensor( real64 (& eigs)[3], real64 (& eigvecs)[3][3], real64 (& NegativeProjector)[6][6] )
{
  //test for repeated eigenvalues
  bool repeatedEigenvalues = false;
  real64 tol = 1e-12;
  if( fabs( eigs[0] - eigs[1] ) < tol || fabs( eigs[0]-eigs[2] ) < tol || fabs( eigs[1]-eigs[2] ) < tol )
  {
    repeatedEigenvalues = true;
  }

  //init QVoigt
  real64 Qi[6][6] = {};
  //init GVoigt
  real64 Gsym[6][6] = {};
  real64 Gji[6][6] = {};

  //compute projector
  for( int i = 0; i < 3; i++ )
  {
    real64 ithEigenVector[3] = {};
    ithEigenVector[0]=eigvecs[0][i];
    ithEigenVector[1]=eigvecs[1][i];
    ithEigenVector[2]=eigvecs[2][i];
    //First Part
    //compute Qi
    QTensor( ithEigenVector, Qi );
    //do update
    LvArray::tensorOps::scale< 6, 6 >( Qi, heaviside( -eigs[i] ));
    LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Qi );
    if( !repeatedEigenvalues )
    {

      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3] = {};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        //Do update
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (fmin( eigs[i], 0.0 ) - fmin( eigs[j], 0.0 ))/(2*(eigs[i]-eigs[j])));
        LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Gsym );
      }
    }
    else
    {
      for( int j = 0; j < 3; j++ )
      {
        real64 jthEigenVector[3] = {};
        jthEigenVector[0]=eigvecs[0][j];
        jthEigenVector[1]=eigvecs[1][j];
        jthEigenVector[2]=eigvecs[2][j];
        if( i == j )
        {
          continue;
        }
        //compute Gij and Gji
        GTensor( ithEigenVector, jthEigenVector, Gsym );
        GTensor( jthEigenVector, ithEigenVector, Gji );
        LvArray::tensorOps::add< 6, 6 >( Gsym, Gji );
        LvArray::tensorOps::scale< 6, 6 >( Gsym, 0.5 * (heaviside( -eigs[i] ) + heaviside( -eigs[j] ))/4 );
        //do update
        LvArray::tensorOps::add< 6, 6 >( NegativeProjector, Gsym );
      }
    }
  }

}

/**
 *@brief This function tests the GetStiffness function from DamageSpectral.hpp.
 *@param[in] c 4th order constitutive tensor in Voigt form.
 *@param[in] strain Strain tensor in Voigt form.
 *@param[in] damage Scalar damage value.
 */
GEOSX_HOST_DEVICE inline
void GetStiffnessTest( real64 (& c)[6][6], real64 (& strain)[6], real64 damage )
{

  //Spectral Split
  real64 const damageFactor = (1-damage)*(1-damage);
  real64 const mu = 1;
  real64 const lambda = 1;
  //get strain tensor in voigt form
  real64 traceOfStrain = strain[0] + strain[1] + strain[2];
  //get eigenvalues and eigenvectors
  real64 eigenValues[3]={};
  real64 eigenVectors[3][3]={};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
  //construct 4th order IxI tensor
  real64 IxITensor[6][6]={};
  for( int i=0; i < 3; i++ )
  {
    for( int j=0; j < 3; j++ )
    {
      IxITensor[i][j] = 1.0;
    }
  }

  //construct positive part
  real64 cPositive[6][6]={};
  real64 positiveProjector[6][6]={};
  PositiveProjectorTensor( eigenValues, eigenVectors, positiveProjector );
  LvArray::tensorOps::scaledCopy< 6, 6 >( cPositive, IxITensor, lambda*heaviside( traceOfStrain ));
  LvArray::tensorOps::scale< 6, 6 >( positiveProjector, 2*mu );
  LvArray::tensorOps::add< 6, 6 >( cPositive, positiveProjector );

  //construct negative part
  real64 negativeProjector[6][6]={};
  NegativeProjectorTensor( eigenValues, eigenVectors, negativeProjector );
  LvArray::tensorOps::scaledCopy< 6, 6 >( c, IxITensor, lambda*heaviside( -traceOfStrain ));
  LvArray::tensorOps::scale< 6, 6 >( negativeProjector, 2*mu );
  LvArray::tensorOps::add< 6, 6 >( c, negativeProjector );
  //finish up
  LvArray::tensorOps::scale< 6, 6 >( cPositive, damageFactor );
  LvArray::tensorOps::add< 6, 6 >( c, cPositive );

}

/**
 *@brief This function tests the GetStress function from DamageSpectral.hpp.
 *@param[in] strain Strain tensor in Voigt form.
 *@param[in] stress Stress tensor in Voigt form.
 */
GEOSX_HOST_DEVICE inline
void getTestStress( real64 (& strain)[6], real64 (& stress)[6] )
{

  //Spectral split
  real64 const damageFactor = 0.25;
  real64 const mu = 1;
  real64 const lambda = 1;
  //get strain tensor in voigt form
  real64 traceOfStrain = strain[0] + strain[1] + strain[2];
  //get eigenvalues and eigenvectors
  real64 eigenValues[3] = {};
  real64 eigenVectors[3][3] = {};
  LvArray::tensorOps::symEigenvectors< 3 >( eigenValues, eigenVectors, strain );
  //transpose eigenVectors matrix to match convention
  real64 temp[3][3] = {};
  LvArray::tensorOps::transpose< 3, 3 >( temp, eigenVectors );
  LvArray::tensorOps::copy< 3, 3 >( eigenVectors, temp );
  real64 tracePlus = fmax( traceOfStrain, 0.0 );
  real64 traceMinus = fmin( traceOfStrain, 0.0 );
  //build symmetric matrices of positive and negative eigenvalues
  real64 Itensor[6] = {};
  for( int i = 0; i < 3; i++ )
  {
    Itensor[i] = 1;
  }
  real64 positivePartOfStrain[6] = {};
  real64 negativePartOfStrain[6] = {};
  PositivePartOfTensor( eigenValues, eigenVectors, positivePartOfStrain );
  NegativePartOfTensor( eigenValues, eigenVectors, negativePartOfStrain );
  real64 positiveStress[6] = {};
  real64 negativeStress[6] = {};
  LvArray::tensorOps::scaledCopy< 6 >( positiveStress, Itensor, lambda*tracePlus );
  LvArray::tensorOps::scaledCopy< 6 >( negativeStress, Itensor, lambda*traceMinus );
  LvArray::tensorOps::scaledAdd< 6 >( positiveStress, positivePartOfStrain, 2*mu );
  LvArray::tensorOps::scaledAdd< 6 >( negativeStress, negativePartOfStrain, 2*mu );
  LvArray::tensorOps::copy< 6 >( stress, negativeStress );
  LvArray::tensorOps::scaledAdd< 6 >( stress, positiveStress, damageFactor );
}

} //namespacegeosx

#endif /* GEOSX_CONSTITUTIVE_SOLID_AUXFUNSPECTRAL_HPP_ */
