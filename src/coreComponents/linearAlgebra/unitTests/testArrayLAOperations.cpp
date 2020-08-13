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
 * @file testDenseLAOperations.cpp
 */

#include "gtest/gtest.h"

#include <numeric>

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

using namespace geosx;

using INDEX_TYPE = std::ptrdiff_t;

static real64 const machinePrecision =
  20.0 * std::numeric_limits< real64 >::epsilon();
static real64 const pi = std::atan( 1.0 ) * 4.0;

template< typename LAI >
void
vector_norm1_test()
{
  INDEX_TYPE N = 24;
  array1d< real64 > vec( N );

  // Populate vector with random coefficients
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute norm1
  real64 norm1 = 0.0;
  for( INDEX_TYPE i = 0; i < vec.size(); ++i )
  {
    norm1 += std::abs( vec( i ) );
  }

  // Check
  EXPECT_NEAR( norm1, LAI::vectorNorm1( vec ), norm1 * machinePrecision );
}

template< typename LAI >
void
vector_norm2_test()
{
  INDEX_TYPE N = 24;
  array1d< real64 > vec( N );

  // Populate vector with random coefficients
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute norm2
  real64 norm2 = 0.0;
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    norm2 += std::pow( vec( i ), 2 );
  }
  norm2 = std::sqrt( norm2 );

  // Check
  EXPECT_NEAR( norm2, LAI::vectorNorm2( vec ), norm2 * machinePrecision );
}

template< typename LAI >
void
vector_normInf_test()
{
  INDEX_TYPE N = 24;
  array1d< real64 > vec( N );

  // Populate vector with random coefficients
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute normInf
  real64 normInf = std::abs( vec( 0 ) );

  if( N > 1 )
  {
    for( INDEX_TYPE i = 1; i < N; ++i )
    {
      normInf = std::max( normInf, std::abs( vec( i ) ) );
    }
  }

  // Check
  EXPECT_NEAR( normInf, LAI::vectorNormInf( vec ), normInf * machinePrecision );
}

template< typename LAI >
void
determinant_test()
{
  array2d< real64 > Laplacian1d;
  real64 determinant;

  real64 theta;
  real64 lambda_max;
  real64 lambda_min;
  real64 N_real;

  for( INDEX_TYPE N = 1; N <= 100; ++N )
  {
    // Construct 1d discrete Laplacian with Dirichlet boundary conditions
    // at both ends
    Laplacian1d.resize( N, N );
    Laplacian1d.setValues< serialPolicy >( 0.0 );
    for( INDEX_TYPE i = 0; i < N; ++i )
    {
      for( INDEX_TYPE j = 0; j < N; ++j )
      {
        if( i == j )
        {
          Laplacian1d( i, i ) = 2;
        }
        else if( std::abs( i - j ) == 1 )
        {
          Laplacian1d( i, j ) = -1;
        }
      }
    }

    // Check
    N_real = static_cast< real64 >( N );
    theta = pi * 0.5 / ( 2.0 * N_real + 1 );
    lambda_min = std::pow( 2.0 * std::sin( theta ), 2 );
    theta = pi * ( N_real - 0.5 ) / ( 2 * N_real + 1 );
    lambda_max = std::pow( 2.0 * std::sin( theta ), 2 );

    determinant = N_real + 1.0;  // Exact determinant
    EXPECT_NEAR( determinant,
                 LAI::determinant( Laplacian1d ),
                 ( lambda_max / lambda_min ) * machinePrecision );
  }
}

template< typename LAI >
void
matrix_normInf_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > mat( M, N );

  // Populate matrix with random coefficients
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute normInf
  real64 normInf = 0.0;
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    real64 tmp = 0.0;
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      tmp += std::abs( mat( i, j ) );
    }
    normInf = std::max( normInf, tmp );
  }

  // Check
  EXPECT_NEAR( normInf, LAI::matrixNormInf( mat ), normInf * machinePrecision );
}

template< typename LAI >
void
matrix_norm1_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > mat( M, N );

  // Populate matrix with random coefficients
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute norm1
  array1d< real64 > tmp( N );
  tmp.setValues< serialPolicy >( 0 );
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      tmp( j ) += std::abs( mat( i, j ) );
    }
  }
  real64 * norm1 = std::max_element( tmp.begin(), tmp.end() );

  // Check
  EXPECT_NEAR( *norm1, LAI::matrixNorm1( mat ), *norm1 * machinePrecision );
}

template< typename LAI >
void
matrix_normFrobenius_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > mat( M, N );

  // Populate matrix with random coefficients
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute normFrobeniusm
  real64 normFrobenius = 0.0;
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      normFrobenius += std::pow( mat( i, j ), 2 );
    }
  }
  normFrobenius = std::sqrt( normFrobenius );

  // Check
  EXPECT_NEAR( normFrobenius,
               LAI::matrixNormFrobenius( mat ),
               normFrobenius * machinePrecision );
}

template< typename LAI >
void
vector_vector_add_test()
{
  INDEX_TYPE N = 24;
  array1d< real64 > vec1( N );
  array1d< real64 > vec2( N );
  array1d< real64 > vecSum( N );

  // Populate vectors with random coefficients
  LAI::vectorRand( vec1, LAI::RandomNumberDistribution::UNIFORM_01 );
  LAI::vectorRand( vec2, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute vector sum vecSum = alpha*vec1 + vec2
  real64 alpha = 3.0;
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    vecSum( i ) = alpha * vec1( i ) + vec2( i );
  }

  // Compute v2 = alpha*vec1 + vec2
  LAI::vectorVectorAdd( vec1, vec2, alpha );

  // Check
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    EXPECT_NEAR( vec2( i ), vecSum( i ), machinePrecision );
  }
}

template< typename LAI >
void
matrix_matrix_add_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > mat1( M, N );
  array2d< real64 > mat2( M, N );
  array2d< real64 > matSum( M, N );

  // Populate vectors with random coefficients
  // Populate matrix with random coefficients
  LAI::matrixRand( mat1, LAI::RandomNumberDistribution::UNIFORM_01 );
  // Populate matrix with random coefficients
  LAI::matrixRand( mat2, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute vector sum matSum = alpha*mat1 + mat2
  real64 alpha = 3.0;
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      matSum( i, j ) = alpha * mat1( i, j ) + mat2( i, j );
    }
  }

  // Compute mat2 = alpha*mat1 + mat2
  LAI::matrixMatrixAdd( mat1, mat2, alpha );

  // Check
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      EXPECT_NEAR( mat2( i, j ), matSum( i, j ), machinePrecision );
    }
  }
}

template< typename LAI >
void
vector_scale_test()
{
  INDEX_TYPE N = 24;
  array1d< real64 > vec( N );
  array1d< real64 > vecScaled( N );

  // Populate vectors with random coefficients
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute vector vecScaled = alpha*vec
  real64 alpha = 3.0;
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    vecScaled( i ) = alpha * vec( i );
  }

  // Compute vec = alpha*vec
  LAI::vectorScale( alpha, vec );

  // Check
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    EXPECT_NEAR( vec( i ), vecScaled( i ), machinePrecision );
  }
}

template< typename LAI >
void
matrix_scale_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > mat( M, N );
  array2d< real64 > matScaled( M, N );

  // Populate vectors with random coefficients
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute vector matScaled = alpha*mat
  real64 alpha = 3.0;
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      matScaled( i, j ) = alpha * mat( i, j );
    }
  }

  // Compute mat = alpha*mat
  LAI::matrixScale( alpha, mat );

  // Check
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      EXPECT_NEAR( mat( i, j ), matScaled( i, j ), machinePrecision );
    }
  }
}

template< typename LAI >
void
vector_dot_test()
{
  INDEX_TYPE N = 6;
  array1d< real64 > vec1( N );
  array1d< real64 > vec2( N );

  // Populate vectors with random coefficients
  LAI::vectorRand( vec1, LAI::RandomNumberDistribution::UNIFORM_01 );
  LAI::vectorRand( vec2, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

  // Compute vector vec1_dot_v2 = vec1^t*vec2
  real64 vec1DotVec2 = 0.0;
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    vec1DotVec2 += vec1( i ) * vec2( i );
  }

  // Check
  EXPECT_NEAR( vec1DotVec2,
               LAI::vectorDot( vec1, vec2 ),
               static_cast< real64 >( N ) * machinePrecision );
}

template< typename LAI >
void
matrix_vector_multiply_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > A( M, N );
  array1d< real64 > X( N );
  array1d< real64 > Y( M );
  array1d< real64 > vecResult( M );

  // Populate matrix and vectors with random coefficients
  LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
  LAI::vectorRand( X, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
  LAI::vectorRand( Y, LAI::RandomNumberDistribution::NORMAL_01 );

  // Compute vector vecResults = alpha*A*X + beta*Y
  real64 alpha = 3.0;
  real64 beta = 7.0;
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    vecResult( i ) = beta * Y( i );
    for( INDEX_TYPE j = 0; j < N; ++j )
    {
      vecResult( i ) += alpha * A( i, j ) * X( j );
    }
  }

  // Compute Y = alpha*A*X + beta*Y
  LAI::matrixVectorMultiply( A, X, Y, alpha, beta );

  // Check
  for( INDEX_TYPE i = 0; i < M; ++i )
  {
    EXPECT_NEAR( vecResult( i ), Y( i ), static_cast< real64 >( N + 1 ) * machinePrecision );
  }
}

template< typename LAI >
void
matrixT_vector_multiply_test()
{
  INDEX_TYPE M = 6;
  INDEX_TYPE N = 24;
  array2d< real64 > A( M, N );
  array1d< real64 > X( M );
  array1d< real64 > Y( N );
  array1d< real64 > vecResult( N );

  // Populate matrix and vectors with random coefficients
  LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
  LAI::vectorRand( X, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
  LAI::vectorRand( Y, LAI::RandomNumberDistribution::NORMAL_01 );

  // Compute vector vecResults = alpha*transpose(A)*X + beta*Y
  real64 alpha = 3.0;
  real64 beta = 7.0;
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    vecResult( i ) = beta * Y( i );
    for( INDEX_TYPE j = 0; j < M; ++j )
    {
      vecResult( i ) += alpha * A( j, i ) * X( j );
    }
  }

  // Compute Y = alpha*A*X + beta*Y
  LAI::matrixTVectorMultiply( A, X, Y, alpha, beta );

  // Check
  for( INDEX_TYPE i = 0; i < N; ++i )
  {
    EXPECT_NEAR( vecResult( i ), Y( i ), static_cast< real64 >( N + 1 ) * machinePrecision );
  }
}

template< typename LAI >
void
matrix_matrix_multiply_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 24 );
  M_indices.emplace_back( 100 );
  array1d< INDEX_TYPE > N_indices( M_indices );
  array1d< INDEX_TYPE > K_indices( M_indices );

  array2d< real64 > A;
  array2d< real64 > B;
  array2d< real64 > C;
  array2d< real64 > matResult;
  real64 alpha = 3.0;
  real64 beta = 7.0;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      for( INDEX_TYPE K : K_indices )
      {
        // Resize matrices
        A.resize( M, K );
        B.resize( K, N );
        C.resize( M, N );
        matResult.resize( M, N );

        // Populate matrices with random coefficients
        LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
        LAI::matrixRand( B, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
        LAI::matrixRand( C, LAI::RandomNumberDistribution::NORMAL_01 );

        // Compute matrix matResult = alpha*A*B + beta*C
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            matResult( i, j ) = beta * C( i, j );
            for( INDEX_TYPE l = 0; l < K; ++l )
            {
              matResult( i, j ) += alpha * A( i, l ) * B( l, j );
            }
          }
        }

        // Compute C = alpha*A*B + beta*C
        LAI::matrixMatrixMultiply( A, B, C, alpha, beta );

        // Check
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            EXPECT_NEAR( C( i, j ),
                         matResult( i, j ),
                         static_cast< real64 >( K + 1 ) * machinePrecision );
          }
        }
      }
    }
  }
}

template< typename LAI >
void
matrixT_matrix_multiply_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 24 );
  M_indices.emplace_back( 100 );
  array1d< INDEX_TYPE > N_indices( M_indices );
  array1d< INDEX_TYPE > K_indices( M_indices );

  array2d< real64 > A;
  array2d< real64 > B;
  array2d< real64 > C;
  array2d< real64 > matResult;
  real64 alpha = 3.0;
  real64 beta = 7.0;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      for( INDEX_TYPE K : K_indices )
      {
        // Resize matrices
        A.resize( K, M );
        B.resize( K, N );
        C.resize( M, N );
        matResult.resize( M, N );

        // Populate matrices with random coefficients
        LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
        LAI::matrixRand( B, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
        LAI::matrixRand( C, LAI::RandomNumberDistribution::NORMAL_01 );

        // Compute matrix matResult = alpha*A*B + beta*C
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            matResult( i, j ) = beta * C( i, j );
            for( INDEX_TYPE l = 0; l < K; ++l )
            {
              matResult( i, j ) += alpha * A( l, i ) * B( l, j );
            }
          }
        }

        // Compute C = alpha*A*B + beta*C
        LAI::matrixTMatrixMultiply( A, B, C, alpha, beta );

        // Check
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            EXPECT_NEAR( C( i, j ),
                         matResult( i, j ),
                         static_cast< real64 >( K + 1 ) * machinePrecision );
          }
        }
      }
    }
  }
}

template< typename LAI >
void
matrix_matrixT_multiply_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 24 );
  M_indices.emplace_back( 100 );
  array1d< INDEX_TYPE > N_indices( M_indices );
  array1d< INDEX_TYPE > K_indices( M_indices );

  array2d< real64 > A;
  array2d< real64 > B;
  array2d< real64 > C;
  array2d< real64 > matResult;
  real64 alpha = 3.0;
  real64 beta = 7.0;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      for( INDEX_TYPE K : K_indices )
      {
        // Resize matrices
        A.resize( M, K );
        B.resize( N, K );
        C.resize( M, N );
        matResult.resize( M, N );

        // Populate matrices with random coefficients
        LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
        LAI::matrixRand( B, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
        LAI::matrixRand( C, LAI::RandomNumberDistribution::NORMAL_01 );

        // Compute matrix matResult = alpha*A*B + beta*C
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            matResult( i, j ) = beta * C( i, j );
            for( INDEX_TYPE l = 0; l < K; ++l )
            {
              matResult( i, j ) += alpha * A( i, l ) * B( j, l );
            }
          }
        }

        // Compute C = alpha*A*B + beta*C
        LAI::matrixMatrixTMultiply( A, B, C, alpha, beta );

        // Check
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            EXPECT_NEAR( C( i, j ),
                         matResult( i, j ),
                         static_cast< real64 >( K + 1 ) * machinePrecision );
          }
        }
      }
    }
  }
}

template< typename LAI >
void
matrixT_matrixT_multiply_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 24 );
  M_indices.emplace_back( 100 );
  array1d< INDEX_TYPE > N_indices( M_indices );
  array1d< INDEX_TYPE > K_indices( M_indices );

  array2d< real64 > A;
  array2d< real64 > B;
  array2d< real64 > C;
  array2d< real64 > matResult;
  real64 alpha = 3.0;
  real64 beta = 7.0;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      for( INDEX_TYPE K : K_indices )
      {
        // Resize matrices
        A.resize( K, M );
        B.resize( N, K );
        C.resize( M, N );
        matResult.resize( M, N );

        // Populate matrices with random coefficients
        LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );
        LAI::matrixRand( B, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
        LAI::matrixRand( C, LAI::RandomNumberDistribution::NORMAL_01 );

        // Compute matrix matResult = alpha*A*B + beta*C
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            matResult( i, j ) = beta * C( i, j );
            for( INDEX_TYPE l = 0; l < K; ++l )
            {
              matResult( i, j ) += alpha * A( l, i ) * B( j, l );
            }
          }
        }

        // Compute C = alpha*A*B + beta*C
        LAI::matrixTMatrixTMultiply( A, B, C, alpha, beta );

        // Check
        for( INDEX_TYPE i = 0; i < M; ++i )
        {
          for( INDEX_TYPE j = 0; j < N; ++j )
          {
            EXPECT_NEAR( C( i, j ),
                         matResult( i, j ),
                         static_cast< real64 >( K + 1 ) * machinePrecision );
          }
        }
      }
    }
  }
}

template< typename LAI >
void
matrix_inverse_test()
{
  array2d< real64 > Laplacian1d;
  array2d< real64 > Laplacian1dInv;
  real64 exact_entry;
  real64 theta;
  real64 lambda_max;
  real64 lambda_min;
  real64 N_real;

  for( INDEX_TYPE N = 1; N <= 100; ++N )
  {
    // Construct 1d discrete Laplacian with Dirichlet boundary conditions
    // at both ends
    Laplacian1d.resize( N, N );
    Laplacian1d.setValues< serialPolicy >( 0 );
    for( INDEX_TYPE i = 0; i < N; ++i )
    {
      for( INDEX_TYPE j = 0; j < N; ++j )
      {
        if( i == j )
        {
          Laplacian1d( i, i ) = 2;
        }
        else if( std::abs( i - j ) == 1 )
        {
          Laplacian1d( i, j ) = -1;
        }
      }
    }

    // Compute inverse of Laplacian1d
    Laplacian1dInv.resize( N, N );
    LAI::matrixInverse( Laplacian1d, Laplacian1dInv );

    // Check
    N_real = static_cast< real64 >( N );
    theta = pi * 0.5 / ( 2.0 * N_real + 1 );
    lambda_min = std::pow( 2.0 * std::sin( theta ), 2 );
    theta = pi * ( N_real - 0.5 ) / ( 2 * N_real + 1 );
    lambda_max = std::pow( 2.0 * std::sin( theta ), 2 );

    for( INDEX_TYPE i = 0; i < N; ++i )
    {
      for( INDEX_TYPE j = 0; j < N; ++j )
      {
        exact_entry = static_cast< real64 >( ( 1 + std::min( i, j ) ) *
                                             ( N + 1 - ( 1 + std::max( i, j ) ) ) ) /
          static_cast< real64 >( ( N + 1 ) );

        EXPECT_NEAR( Laplacian1dInv( i, j ),
                     exact_entry,
                     ( lambda_max / lambda_min ) * machinePrecision );
      }
    }
  }
}

template< typename LAI >
void
vector_copy_test()
{
  array1d< INDEX_TYPE > N_indices;
  N_indices.emplace_back( 1 );
  N_indices.emplace_back( 6 );
  N_indices.emplace_back( 24 );
  N_indices.emplace_back( 100 );

  array1d< real64 > src;
  array1d< real64 > dst;

  for( INDEX_TYPE N : N_indices )
  {
    src.resize( N );
    dst.resize( N );

    // Populate src vector with random coefficients
    LAI::vectorRand( src, LAI::RandomNumberDistribution::NORMAL_01 );

    // Copy vector, dst = src
    LAI::vectorCopy( src, dst );

    // Check
    for( INDEX_TYPE i = 0; i < N; ++i )
    {
      EXPECT_NEAR( src( i ), dst( i ), machinePrecision );
    }
  }
}

template< typename LAI >
void
matrix_copy_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 24 );
  M_indices.emplace_back( 100 );
  array1d< INDEX_TYPE > N_indices( M_indices );

  array2d< real64 > src;
  array2d< real64 > dst;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      src.resize( M, N );
      dst.resize( M, N );

      // Populate src matrix with random coefficients
      LAI::matrixRand( src, LAI::RandomNumberDistribution::NORMAL_01 );

      // Copy matrix, dst = src
      LAI::matrixCopy( src, dst );

      // Check
      for( INDEX_TYPE i = 0; i < M; ++i )
      {
        for( INDEX_TYPE j = 0; j < N; ++j )
        {
          EXPECT_NEAR( src( i, j ), dst( i, j ), machinePrecision );
        }
      }
    }
  }
}

template< typename LAI >
void
vector_rand_test()
{
  INDEX_TYPE N = 10000;
  array1d< real64 > vec( N );

  // Populate vector with random coefficients

  // --- uniform distribution (0,1);
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_01 );
  real64 * v_max = std::max_element( vec.begin(), vec.end() );
  real64 * v_min = std::min_element( vec.begin(), vec.end() );
  EXPECT_TRUE( 0.0 <= *v_min && *v_max <= 1.0 );

  // --- uniform distribution (-1,1);
  LAI::vectorRand( vec, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
  v_max = std::max_element( vec.begin(), vec.end() );
  v_min = std::min_element( vec.begin(), vec.end() );
  EXPECT_TRUE( -1.0 <= *v_min && *v_max <= 1.0 );

  // --- normal distribution (0,1);
  // TODO: Add normality test
}

template< typename LAI >
void
matrix_rand_test()
{
  INDEX_TYPE M = 99;
  INDEX_TYPE N = 101;
  array2d< real64 > mat( M, N );

  // Populate vector with random coefficients

  // --- uniform distribution (0,1);
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_01 );
  real64 * A_max = std::max_element( mat.begin(), mat.end() );
  real64 * A_min = std::min_element( mat.begin(), mat.end() );
  EXPECT_TRUE( 0.0 <= *A_min && *A_max <= 1.0 );

  // --- uniform distribution (-1,1);
  LAI::matrixRand( mat, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
  A_max = std::max_element( mat.begin(), mat.end() );
  A_min = std::min_element( mat.begin(), mat.end() );
  EXPECT_TRUE( -1.0 <= *A_min && *A_max <= 1.0 );

  // --- normal distribution (0,1);
  // TODO: Add normality test
}

template< typename LAI >
void
set_get_random_number_generator_seed_test()
{
  array1d< int > seedSet( 4 );
  seedSet( 0 ) = 1;
  seedSet( 1 ) = 2;
  seedSet( 2 ) = 3;
  seedSet( 3 ) = 7;  // must be odd
  LAI::setRandomNumberGeneratorSeed( seedSet );

  array1d< int > seedGet( 4 );
  LAI::getRandomNumberGeneratorSeed( seedGet );

  // Check
  for( INDEX_TYPE i = 0; i < seedSet.size(); ++i )
  {
    EXPECT_EQ( seedSet( i ), seedGet( i ) );
  }
}

template< typename LAI >
void
performance_test()
{
  constexpr localIndex MAX_SIZE = 1024 * 8 * 4;

  array2d< double > a( MAX_SIZE, MAX_SIZE );
  array2d< double > b( MAX_SIZE, MAX_SIZE );
  array2d< double > c( MAX_SIZE, MAX_SIZE );

  for( localIndex size = 3; size <= MAX_SIZE; size = localIndex( size * 1.5 + 1 ) )
  {
    a.resize( size, size );
    b.resize( size, size );
    c.resize( size, size );

    LAI::matrixRand( a, LAI::RandomNumberDistribution::UNIFORM_m1p1 );
    LAI::matrixRand( b, LAI::RandomNumberDistribution::UNIFORM_m1p1 );

    std::chrono::duration< double > multiplicationTime {};
    localIndex const nIter = MAX_SIZE / ( size + 1 ) + 1;
    for( localIndex iter = 0; iter < nIter; ++iter )
    {
      auto start = std::chrono::high_resolution_clock::now();
      LAI::matrixMatrixMultiply( a, b, c, 1.0, 1.0 );
      auto end = std::chrono::high_resolution_clock::now();
      multiplicationTime += end - start;
    }

    GEOSX_LOG( std::setiosflags( std::ios::scientific )
               << std::setprecision( 5 ) << "size = " << size << ",\tniter : " << nIter
               << ",\ttime : " << multiplicationTime.count() << ",\tscaled time : "
               << multiplicationTime.count() / ( size * size * size ) );
  }
}

template< typename LAI >
void
matrix_svd_test()
{
  array1d< INDEX_TYPE > M_indices;
  M_indices.emplace_back( 1 );
  M_indices.emplace_back( 2 );
  M_indices.emplace_back( 3 );
  M_indices.emplace_back( 4 );
  M_indices.emplace_back( 5 );
  M_indices.emplace_back( 6 );
  M_indices.emplace_back( 7 );
  M_indices.emplace_back( 8 );
  array1d< INDEX_TYPE > N_indices( M_indices );

  array2d< real64 > A;
  array2d< real64 > U;
  array1d< real64 > S_vec;
  array2d< real64 > VT;

  array2d< real64 > A_svd;
  array2d< real64 > UTU;
  array2d< real64 > U_extended;
  array2d< real64 > S_extended;
  array2d< real64 > VT_extended;
  array2d< real64 > work0;

  for( INDEX_TYPE M : M_indices )
  {
    for( INDEX_TYPE N : N_indices )
    {
      int const minDim = ( M > N ) ? N : M;

      // Resize matrices
      A.resize( M, N );
      A_svd.resize( M, N );
      U.resize( M, minDim );
      S_vec.resize( minDim );
      VT.resize( minDim, N );
      UTU.resize( minDim, minDim );
      U_extended.resize( M, M );
      S_extended.resize( M, N );
      VT_extended.resize( N, N );
      work0.resize( M, N );

      // Populate matrix A with random coefficients
      LAI::matrixRand( A, LAI::RandomNumberDistribution::UNIFORM_01 );

      // Compute the SVD of A
      LAI::matrixSVD( A, U, S_vec, VT );

      // 1) Check that we recover the matrix with the decomposition

      // fill U_extended
      for( int i = 0; i < U_extended.size( 0 ); ++i )
      {
        for( int j = 0; j < U_extended.size( 1 ); ++j )
        {
          U_extended( i, j ) = 0.0;
        }
      }
      for( int i = 0; i < U.size( 0 ); ++i )
      {
        for( int j = 0; j < U.size( 1 ); ++j )
        {
          U_extended( i, j ) = U( i, j );
        }
      }

      // fill VT_extended
      for( int i = 0; i < VT_extended.size( 0 ); ++i )
      {
        for( int j = 0; j < VT_extended.size( 1 ); ++j )
        {
          VT_extended( i, j ) = 0.0;
        }
      }
      for( int i = 0; i < VT.size( 0 ); ++i )
      {
        for( int j = 0; j < VT.size( 1 ); ++j )
        {
          VT_extended( i, j ) = VT( i, j );
        }
      }

      // fill S_extended
      for( int i = 0; i < S_extended.size( 0 ); ++i )
      {
        for( int j = 0; j < S_extended.size( 1 ); ++j )
        {
          S_extended( i, j ) = 0.0;
        }
      }
      for( int i = 0; i < S_vec.size(); ++i )
      {
        S_extended( i, i ) = S_vec( i );
      }

      LAI::matrixMatrixMultiply( S_extended, VT_extended, work0 );
      LAI::matrixMatrixMultiply( U_extended, work0, A_svd );

      for( INDEX_TYPE i = 0; i < M; ++i )
      {
        for( INDEX_TYPE j = 0; j < N; ++j )
        {
          EXPECT_NEAR( A_svd( i, j ), A( i, j ), machinePrecision );
        }
      }

      // 2) Check that U is orthonormal
      LAI::matrixTMatrixMultiply( U, U, UTU );
      for( int i = 0; i < UTU.size( 0 ); ++i )
      {
        UTU( i, i ) = 1 - UTU( i, i );
      }

      EXPECT_NEAR( 0.0, LAI::matrixNormFrobenius( UTU ), machinePrecision );
    }
  }
}

TEST( Array1D, vectorNorm1 )
{
  vector_norm1_test< BlasLapackLA >();
}

TEST( Array1D, vectorNorm2 )
{
  vector_norm2_test< BlasLapackLA >();
}

TEST( Array1D, vectorNormInf )
{
  vector_normInf_test< BlasLapackLA >();
}

TEST( Array2D, determinant )
{
  determinant_test< BlasLapackLA >();
}

TEST( Array2D, matrixNormInf )
{
  matrix_normInf_test< BlasLapackLA >();
}

TEST( Array2D, matrixNorm1 )
{
  matrix_norm1_test< BlasLapackLA >();
}

TEST( Array2D, matrixNormFrobenius )
{
  matrix_normFrobenius_test< BlasLapackLA >();
}

TEST( Array1D, vectorVectorAdd )
{
  vector_vector_add_test< BlasLapackLA >();
}

TEST( Array2D, matrixMatrixAdd )
{
  matrix_matrix_add_test< BlasLapackLA >();
}

TEST( Array1D, vectorScale )
{
  vector_scale_test< BlasLapackLA >();
}

TEST( Array2D, matrixScale )
{
  matrix_scale_test< BlasLapackLA >();
}

TEST( Array1D, vectorDot )
{
  vector_dot_test< BlasLapackLA >();
}

TEST( Array2D, matrixVectorMultiply )
{
  matrix_vector_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixTVectorMultiply )
{
  matrixT_vector_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixMatrixMultiply )
{
  matrix_matrix_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixTMatrixMultiply )
{
  matrixT_matrix_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixMatrixTMultiply )
{
  matrix_matrixT_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixTMatrixTMultiply )
{
  matrixT_matrixT_multiply_test< BlasLapackLA >();
}

TEST( Array2D, matrixInverse )
{
  matrix_inverse_test< BlasLapackLA >();
}

TEST( Array1D, vectorCopy )
{
  vector_copy_test< BlasLapackLA >();
}

TEST( Array2D, matrixCopy )
{
  matrix_copy_test< BlasLapackLA >();
}

TEST( Array1D, vectorRand )
{
  vector_rand_test< BlasLapackLA >();
}

TEST( Array2D, matrixRand )
{
  matrix_rand_test< BlasLapackLA >();
}

TEST( DenseLAInterface, setGetRandomNumberGeneratorSeed )
{
  set_get_random_number_generator_seed_test< BlasLapackLA >();
}

TEST( DenseLAInterface, performanceTest )
{
  // performance_test<BlasLapackLA>();
  SUCCEED();
}

TEST( DenseLAInterface, matrixSVD )
{
  matrix_svd_test< BlasLapackLA >();
}

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
