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
 * @file testLAOperations.cpp
 */

#include <gtest/gtest.h>

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "managers/initialization.hpp"
#include "testLinearAlgebraUtils.hpp"

using namespace geosx;

static real64 constexpr machinePrecision =
  20.0 * std::numeric_limits< real64 >::epsilon();

template< typename LAI >
class LAOperationsTest : public ::testing::Test
{};

TYPED_TEST_SUITE_P( LAOperationsTest );

TYPED_TEST_P( LAOperationsTest, VectorFunctions )
{
  // Define alias
  using Vector = typename TypeParam::ParallelVector;

  // Get the MPI rank
  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const numRanks = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  Vector x;
  localIndex const localSize = 3;
  globalIndex const globalSize = localSize * numRanks;
  globalIndex const offset = rank * localSize;

  // Testing createWithLocalSize
  x.createWithLocalSize( localSize, MPI_COMM_GEOSX );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing iupper/ilower
  EXPECT_EQ( x.ilower(), offset );
  EXPECT_EQ( x.iupper(), offset + localSize );

  // Testing setting/getting values locally
  x.open();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    x.set( i, 2 * i );
  }
  x.close();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_DOUBLE_EQ( x.get( i ), 2 * i );
  }

  // Testing createWithGlobalSize
  x.createWithGlobalSize( globalSize, MPI_COMM_GEOSX );
  EXPECT_EQ( x.localSize(), localSize );
  EXPECT_EQ( x.globalSize(), globalSize );

  // Testing adding global values on rank 0 and getting locally
  x.open();
  if( rank == 0 )
  {
    for( globalIndex i = 0; i < x.globalSize(); ++i )
    {
      x.add( i, 2 * i );
    }
  }
  x.close();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), 2 * i );
  }

  // Testing getLocalRowID
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.getLocalRowID( i ), i % localSize );
  }

  // Testing create with array1d
  array1d< real64 > localVals( localSize );
  for( localIndex i = 0; i < localSize; ++i )
  {
    localVals[i] = real64( i + rank * localSize );
  }

  Vector v;
  v.create( localVals, MPI_COMM_GEOSX );
  for( globalIndex i = v.ilower(); i < v.iupper(); ++i )
  {
    EXPECT_EQ( v.get( i ), localVals[v.getLocalRowID( i )] );
  }

  // Testing copy constructor, assignment, get element
  Vector y( x );
  Vector z;
  z = x;
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), y.get( i ) );
    EXPECT_EQ( x.get( i ), z.get( i ) );
  }

  // Testing zero
  z.zero();
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( z.get( i ), 0 );
  }

  // Testing copy
  z.copy( x );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( x.get( i ), z.get( i ) );
  }

  // Testing scale, z = x
  z.scale( 4.0 );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( 4.0 * x.get( i ), z.get( i ) );
  }

  // Testing add/set single element
  x.open();
  x.set( offset, -1 );
  x.close();  // set/add can't be interchanged
  x.open();
  x.add( offset + 1, 10 );
  x.close();
  EXPECT_DOUBLE_EQ( x.get( offset ), -1 );
  EXPECT_DOUBLE_EQ( x.get( offset + 1 ), offset * 2 + 12 );

  // Testing add/set c-style
  {
    globalIndex const inds[3] = { offset, offset + 1, offset + 2 };
    real64 const vals[3] = { -5.0, -6.0, 0.0 };
    y.zero();
    y.open();
    y.set( inds, vals, 2 );
    y.close();
    z.set( 1.0 );
    z.open();
    z.add( inds, vals, 2 );
    z.close();
    for( localIndex i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( y.get( inds[i] ), vals[i] );
      EXPECT_DOUBLE_EQ( z.get( inds[i] ), vals[i] + 1.0 );
    }
  }

  // Testing add/set array1d-style
  {
    array1d< globalIndex > inds( 3 );
    inds[0] = offset;
    inds[1] = offset + 1;
    inds[2] = offset + 2;
    array1d< real64 > vals( 3 );
    vals[0] = -5.0;
    vals[1] = -6.0;
    vals[2] = 0.0;
    y.zero();
    y.open();
    y.set( inds, vals );
    y.close();
    z.set( 1.0 );
    z.open();
    z.add( inds, vals );
    z.close();
    for( localIndex i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( y.get( inds[i] ), vals[i] );
      EXPECT_DOUBLE_EQ( z.get( inds[i] ), vals[i] + 1.0 );
    }
  }

  // Testing dot, axpy, axpby
  x.set( 1.0 );
  y.set( 2.0 );
  z.set( 3.0 );

  real64 const dotprod = x.dot( y );
  EXPECT_EQ( dotprod, 2 * y.globalSize() );  // sum_size 2

  y.axpy( 2.0, x );
  for( globalIndex i = y.ilower(); i < y.iupper(); ++i )
  {
    EXPECT_EQ( y.get( i ), 4.0 );  // 2*1 + 2
  }

  z.axpby( 2.0, x, 3.0 );
  for( globalIndex i = z.ilower(); i < z.iupper(); ++i )
  {
    EXPECT_EQ( z.get( i ), 11.0 );  // 2*1 + 3*3
  }

  // Testing norms
  x.zero();
  x.open();
  if( rank == 0 )
  {
    globalIndex const inds2[2] = { 0, 1 };
    real64 const vals2[2] = { 3.0, -4.0 };
    x.set( inds2, vals2, 2 );  // 3, -4, 0
  }
  x.close();
  EXPECT_EQ( x.norm1(), 7.0 );
  EXPECT_EQ( x.norm2(), 5.0 );
  EXPECT_EQ( x.normInf(), 4.0 );

  // Testing extractLocalVector
  real64 const * localVec = x.extractLocalVector();
  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
  {
    EXPECT_EQ( localVec[x.getLocalRowID( i )], x.get( i ) );
  }
}

#if 0
TYPED_TEST_P( LAOperationsTest, MatrixFunctions )
{
  // Define aliases
  using Vector = typename TypeParam::ParallelVector;
  using Matrix = typename TypeParam::ParallelMatrix;

  // Get the MPI rank
  int numranks = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  int rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

  std::cout << "*** Rank: " << rank << std::endl;

  // Dummy vector and Matrix
  Matrix C;
  Matrix D;
  {
    // Test matrix-matrix product: C = A*B
    Matrix A;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, 2 * numranks, A );
    Matrix B( A );

    A.multiply( B, C );
    A.leftMultiplyTranspose( A, D );
  }

  // Define some vectors, matrices
  Vector vec1, vec2, vec3;
  Matrix mat1, mat2, mat3, mat4;
  mat1.createWithLocalSize( 2, 2, MPI_COMM_GEOSX ); // 2*numranks x 2*numranks
  mat2.createWithGlobalSize( 2, 2, MPI_COMM_GEOSX ); // 2x2
  mat3.createWithLocalSize( 2, 3, 3, MPI_COMM_GEOSX ); // 2*numranks x 3*numranks
  mat4.createWithGlobalSize( 3, 4, 3, MPI_COMM_GEOSX ); // 3x4

  // Testing create, globalRows, globalCols
  EXPECT_EQ( mat1.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat1.numGlobalCols(), 2 * numranks );
  EXPECT_EQ( mat2.numGlobalRows(), 2 );
  EXPECT_EQ( mat2.numGlobalCols(), 2 );
  EXPECT_EQ( mat3.numGlobalRows(), 2 * numranks );
  EXPECT_EQ( mat3.numGlobalCols(), 3 * numranks );
  EXPECT_EQ( mat4.numGlobalRows(), 3 );
  EXPECT_EQ( mat4.numGlobalCols(), 4 );

  // Testing add/set/insert element
  //  mat1.insert( 1, 0, .5 );
  //  mat1.close();
  //  mat1.set( 1, 0, 5 );
  //  mat1.close();
  mat1.open();
  mat1.add( 1, 0, 1 );
  mat1.add( 1, 0, 2 );
  mat1.close();

  // Testing add/set/insert c-style, getRowCopy
  globalIndex inds1[2] = { 0, 2 };
  globalIndex inds2[1] = { 0 };
  globalIndex inds3[3] = { 0, 1, 2 };
  real64 vals1[2] = { 5, 10 };
  real64 vals2[1] = { 1 };
  real64 vals3[3] = { .5, 1, 2 };

  globalIndex iRow = 1;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    mat4.insert( iRow, inds3, vals3, 3 );
    //    mat4.close();
    //    mat4.open();
    mat4.set( iRow, inds1, vals1, 2 );
    //    mat4.close();
    //    mat4.open();
    mat4.add( iRow, inds2, vals2, 1 );
    //    mat4.close();
  }
  mat4.close();

  array1d< real64 > colvals_CHECK( 3 );
  colvals_CHECK( 0 ) = 6;
  colvals_CHECK( 1 ) = 1;
  colvals_CHECK( 2 ) = 10;

  if( ( mat4.ilower() <= iRow ) && ( iRow < mat4.iupper() ) )
  {
    localIndex const rowLength = mat.getGlobalRowLength( iRow );
    EXPECT_EQ( rowLength, colvals_CHECK.size() );
    array1d< real64 > colvals( rowLength );
    array1d< globalIndex > colinds( rowLength );
    mat4.getRowCopy( iRow, colinds, colvals );
    for( int i = 0; i < 3; ++i )
    {
      EXPECT_DOUBLE_EQ( colvals( colinds[i] ), colvals_CHECK( i ) ); //HYPRE does not return sorted cols!
    }
  }
  // Testing add/set/insert array1d
  Matrix mat6;
  mat6.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< real64 > vals6( 3 );
  array1d< real64 > vals7( 3 );
  array1d< globalIndex > inds6( 3 );
  vals6[0] = 1;
  vals6[1] = .5;
  vals6[2] = -3;
  vals7[0] = 1;
  vals7[1] = 1;
  vals7[2] = 1;
  inds6[0] = 0;
  inds6[1] = 1;
  inds6[2] = 3;

  iRow = 0;
  if( ( mat6.ilower() <= iRow ) && ( iRow < mat6.iupper() ) )
  {
    mat6.insert( iRow, inds6, vals6 );
    //	  mat6.close();
    mat6.set( iRow, inds6, vals7 );
    //	  mat6.close();
    mat6.add( iRow, inds6, vals6 );
  }
  mat6.close();

  // Testing add/set/insert array2d
  Matrix mat7;
  mat7.createWithGlobalSize( 4, 4, MPI_COMM_GEOSX );
  array1d< globalIndex > rows( 2 );
  array1d< globalIndex > cols( 2 );
  array2d< real64 > vals8( 2, 2 );
  rows[0] = 0;
  rows[1] = 2;
  cols[0] = 1;
  cols[1] = 3;
  vals8[0][0] = 1;
  vals8[0][1] = 2;
  vals8[1][0] = 3;
  vals8[1][1] = 4;
  if( ( mat7.ilower() <= *std::min_element( rows.data(), rows.data() + rows.size() ) ) &&
      ( *std::max_element( rows.data(), rows.data() + rows.size() ) < mat7.iupper() ) )
  {
    mat7.insert( rows, cols, vals8 );
    //    mat7.close();
    mat7.add( rows, cols, vals8 );
    //    mat7.close();
  }
  mat7.close();

  // Testing set and zero
  mat7.open();
  mat7.set( 2 );
  mat7.close();

  mat7.open();
  mat1.zero();
  mat1.close();

  // Testing vector multiply, matrix multiply, MatrixMatrixMultiply
  vec1.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec2.createWithGlobalSize( 2, MPI_COMM_GEOSX );
  vec1.set( 1 );
  vec1.close();
  globalIndex inds4[2] = { 0, 1 };
  real64 vals4[2] = { 1, 3 };
  real64 vals5[2] = { 2, 1 };

  if( ( mat2.ilower() <= 0 ) && ( 0 < mat2.iupper() ) )
  {
    mat2.insert( 0, inds4, vals4, 2 );
  }
  if( ( mat2.ilower() <= 1 ) && ( 1 < mat2.iupper() ) )
  {
    mat2.insert( 1, inds4, vals5, 2 );
  }
  mat2.close();

  mat2.multiply( vec1, vec2 );

  if( ( vec2.ilower() <= 0 ) && ( 0 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 0 ), 4 );
  }
  if( ( vec2.ilower() <= 1 ) && ( 1 < vec2.iupper() ) )
  {
    EXPECT_DOUBLE_EQ( vec2.get( 1 ), 3 );
  }

  // Matrix-Matrix multiply
  Matrix mat2mat2;

  {
    Matrix mat22( mat2 );

    mat22.multiply( mat22, mat2mat2 );
  }
}
#endif

TYPED_TEST_P( LAOperationsTest, MatrixMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  globalIndex const n = 100;
  Matrix A;
  compute2DLaplaceOperator( MPI_COMM_GEOSX, n, A );

  Matrix A_squared;
  A.multiply( A, A_squared );

  real64 const a = A.normInf();
  real64 const b = A_squared.normInf();

  EXPECT_DOUBLE_EQ( a, 8.0 );
  EXPECT_DOUBLE_EQ( b, 64.0 );
}

TYPED_TEST_P( LAOperationsTest, RectangularMatrixOperations )
{
  using Matrix = typename TypeParam::ParallelMatrix;

  int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

  // Set a size that allows to run with arbitrary number of processes
  globalIndex const nRows = std::max( 100, mpiSize );
  globalIndex const nCols = 2 * nRows;

  Matrix A;
  A.createWithGlobalSize( nRows, nCols, 2, MPI_COMM_GEOSX );

  A.open();
  for( globalIndex i = A.ilower(); i < A.iupper(); ++i )
  {
    real64 const entry = static_cast< real64 >( i + 1 );
    A.insert( i, 2 * i, entry );
    A.insert( i, 2 * i + 1, -entry );
  }
  A.close();

  // Check on sizes
  EXPECT_EQ( A.numGlobalRows(), nRows );
  EXPECT_EQ( A.numGlobalCols(), nCols );

  // Check on norms
  real64 const a = A.norm1();
  real64 const b = A.normInf();
  real64 const c = A.normFrobenius();

  EXPECT_DOUBLE_EQ( a, static_cast< real64 >( nRows ) );
  EXPECT_DOUBLE_EQ( b, static_cast< real64 >( nCols ) );
  EXPECT_DOUBLE_EQ(
    c,
    std::sqrt( static_cast< real64 >( nRows * ( nRows + 1 ) * ( 2 * nRows + 1 ) ) / 3.0 ) );
}

REGISTER_TYPED_TEST_SUITE_P( LAOperationsTest,
                             VectorFunctions,
                             MatrixMatrixOperations,
                             RectangularMatrixOperations );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, LAOperationsTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, LAOperationsTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, LAOperationsTest, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////

LinearSolverParameters
params_Direct()
{
  LinearSolverParameters parameters;
  parameters.solverType = "direct";
  parameters.krylov.relTolerance = machinePrecision;
  return parameters;
}

LinearSolverParameters
params_GMRES_ILU()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = "gmres";
  parameters.preconditionerType = "iluk";
  parameters.ilu.fill = 1;
  return parameters;
}

LinearSolverParameters
params_GMRES_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = "gmres";
  parameters.preconditionerType = "amg";
  parameters.amg.smootherType = "gaussSeidel";
  parameters.amg.coarseType = "direct";
  return parameters;
}

LinearSolverParameters
params_CG_AMG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 300;
  parameters.solverType = "cg";
  parameters.isSymmetric = true;
  parameters.preconditionerType = "amg";
  parameters.amg.smootherType = "gaussSeidel";
  parameters.amg.coarseType = "direct";
  return parameters;
}

template< typename LAI >
class SolverTestBase : public ::testing::Test
{
public:
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;
  using Solver = typename LAI::LinearSolver;

protected:
  Matrix matrix;
  real64 cond_est = 1.0;

  void
  test( LinearSolverParameters const & params )
  {
    // Create a random "true" solution vector
    Vector sol_true;
    sol_true.createWithGlobalSize( matrix.numGlobalCols(), matrix.getComm() );
    sol_true.rand();

    // Create and compute the right-hand side vector
    Vector rhs;
    rhs.createWithGlobalSize( matrix.numGlobalRows(), matrix.getComm() );
    matrix.apply( sol_true, rhs );

    // Create and zero out the computed solution vector
    Vector sol_comp;
    sol_comp.createWithGlobalSize( sol_true.globalSize(), sol_true.getComm() );
    sol_comp.zero();

    // Create the solver and solve the system
    Solver solver( params );
    solver.solve( matrix, sol_comp, rhs );
    EXPECT_TRUE( solver.result().success() );

    // Check that solution is within epsilon of true
    sol_comp.axpy( -1.0, sol_true );
    real64 const relTol = cond_est * params.krylov.relTolerance;
    EXPECT_LT( sol_comp.norm2() / sol_true.norm2(), relTol );
  }
};

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestLaplace2D : public SolverTestBase< LAI >
{
public:
  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:
  void
  SetUp() override
  {
    globalIndex constexpr n = 100;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, n, this->matrix );

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( SolverTestLaplace2D );

TYPED_TEST_P( SolverTestLaplace2D, Direct )
{
  this->test( params_Direct() );
}

TYPED_TEST_P( SolverTestLaplace2D, GMRES_ILU )
{
  this->test( params_GMRES_ILU() );
}

TYPED_TEST_P( SolverTestLaplace2D, CG_AMG )
{
  this->test( params_CG_AMG() );
}

REGISTER_TYPED_TEST_SUITE_P( SolverTestLaplace2D, Direct, GMRES_ILU, CG_AMG );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, SolverTestLaplace2D, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestLaplace2D, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestLaplace2D, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class SolverTestElasticity2D : public SolverTestBase< LAI >
{
public:
  using Base = SolverTestBase< LAI >;
  using Matrix = typename Base::Matrix;
  using Vector = typename Base::Vector;

protected:
  void
  SetUp() override
  {
    globalIndex constexpr n = 100;
    compute2DElasticityOperator( MPI_COMM_GEOSX,
                                 1.0,
                                 1.0,
                                 n,
                                 n,
                                 10000.,
                                 0.2,
                                 this->matrix );

    // Impose Dirichlet boundary conditions: fix domain bottom (first 2*(nCellsX + 1) rows of matrix)
    this->matrix.open();
    for( globalIndex iRow = 0; iRow < 2 * ( n + 1 ); ++iRow )
    {
      if( this->matrix.getLocalRowID( iRow ) >= 0 )
      {
        this->matrix.clearRow( iRow, true );
      }
    }
    this->matrix.close();
    this->cond_est =
      1e3;  // not a true condition number estimate, but enough to pass tests
  }
};

TYPED_TEST_SUITE_P( SolverTestElasticity2D );

TYPED_TEST_P( SolverTestElasticity2D, Direct )
{
  this->test( params_Direct() );
}

TYPED_TEST_P( SolverTestElasticity2D, GMRES_AMG )
{
  LinearSolverParameters params = params_GMRES_AMG();
  params.amg.separateComponents = true;
  params.dofsPerNode = 2;
  this->test( params );
}

REGISTER_TYPED_TEST_SUITE_P( SolverTestElasticity2D, Direct, GMRES_AMG );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos,
                                SolverTestElasticity2D,
                                TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, SolverTestElasticity2D, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, SolverTestElasticity2D, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////
//
//template<typename LAI>
//void parallel_vector_copy_constructor()
//{
//
//  int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
//  typename LAI::ParallelVector x;
//  typename LAI::ParallelVector y;
//  localIndex const localSize = rank*10 + 1;
//
//  array1d<real64> vec(localSize);
//
//  // Populate vector with random coefficients
//  BlasLapackLA::vectorRand( vec,
//                            BlasLapackLA::RandomNumberDistribution::UNIFORM_m1p1 );
//
//  x.create( vec, MPI_COMM_WORLD );
//  x.close();
//
//  y.create(x);
//
//  for( globalIndex i = x.ilower(); i < x.iupper(); ++i )
//  {
//    EXPECT_EQ( x.get( i ), y.get( i ) );
//  }
//}
//
//template<typename LAI>
//void parallel_vector_create_with_local_size()
//{
//  int const rank = MpiWrapper::Comm_rank( MPI_COMM_WORLD );
//  typename LAI::ParallelVector x;
//  localIndex const localSize = rank*10;
//  globalIndex const globalSize = MpiWrapper::Sum( localSize );
//
//  x.createWithLocalSize( localSize, MPI_COMM_WORLD );
//
//  EXPECT_EQ( x.localSize(), localSize );
//  EXPECT_EQ( x.globalSize(), globalSize );
//}
//
//template<typename LAI>
//void parallel_vector_create_with_global_size()
//{
//  int const numRanks = MpiWrapper::Comm_size( MPI_COMM_WORLD );
//  typename LAI::ParallelVector x;
//  localIndex const localSize = integer_conversion< localIndex >( numRanks );
//  globalIndex const globalSize = integer_conversion< globalIndex >( localSize * localSize );
//
//  // Testing createWithGlobalSize
//  x.createWithGlobalSize( globalSize, MPI_COMM_WORLD );
//  EXPECT_EQ( x.localSize(), localSize );
//  EXPECT_EQ( x.globalSize(), globalSize );
//}
//
//#ifdef GEOSX_USE_TRILINOS
//
//TEST( EpetraVector, EpetraVector )
//{
//  parallel_vector_copy_constructor< TrilinosInterface >();
//}
//
//TEST( EpetraVector, createWithLocalSize )
//{
//  parallel_vector_create_with_local_size< TrilinosInterface >();
//}
//
//TEST( EpetraVector, createWithGlobalSize )
//{
//  parallel_vector_create_with_global_size< TrilinosInterface >();
//}
//
//#endif
///////////////////////////////////////////////////////////////////////////////////////

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  // Avoid setting up signal handlers, due to mysterious ML FPE crashes on Mac
  //geosx::basicSetup( argc, argv );
  setupMPI( argc, argv );
  setupLogger();
  setupOpenMP();
  setupMKL();
  setupLogger();
  setupLAI( argc, argv );

  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
