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
 * @file testKrylovSolvers.cpp
 */

#include <gtest/gtest.h>

#include "testLinearAlgebraUtils.hpp"

#include "common/DataTypes.hpp"
#include "managers/initialization.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/BlockOperatorWrapper.hpp"
#include "linearAlgebra/solvers/PreconditionerIdentity.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"

using namespace geosx;

LinearSolverParameters
params_CG()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 500;
  parameters.solverType = "cg";
  parameters.isSymmetric = true;
  return parameters;
}

LinearSolverParameters
params_BiCGSTAB()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 500;
  parameters.solverType = "bicgstab";
  return parameters;
}

LinearSolverParameters
params_GMRES()
{
  LinearSolverParameters parameters;
  parameters.krylov.relTolerance = 1e-8;
  parameters.krylov.maxIterations = 500;
  parameters.solverType = "gmres";
  return parameters;
}

template< typename OPERATOR, typename PRECOND, typename VECTOR >
class KrylovSolverTestBase : public ::testing::Test
{
public:
  using Base = ::testing::Test;

  // Constructor that allows subclass to pass initialization args to matrix/precond
  // Unfortunately, only one set of args can be passed to both and none to vectors
  template< typename... ARGS >
  KrylovSolverTestBase( ARGS &&... args ) :
    Base(),
    matrix( std::forward< ARGS >( args )... ),
    precond( std::forward< ARGS >( args )... )
  {}

protected:
  OPERATOR matrix;
  PRECOND precond;
  VECTOR sol_true;
  VECTOR sol_comp;
  VECTOR rhs_true;
  real64 cond_est = 1.0;

  void
  test( LinearSolverParameters const & params )
  {
    sol_true.rand();
    sol_comp.zero();
    matrix.apply( sol_true, rhs_true );

    // Create the solver and solve the system
    using Vector = typename OPERATOR::Vector;
    std::unique_ptr< KrylovSolver< Vector > > const solver =
      KrylovSolver< Vector >::Create( params, matrix, precond );
    solver->solve( rhs_true, sol_comp );
    EXPECT_TRUE( solver->result().success() );

    // Check that solution is within epsilon of true
    sol_comp.axpy( -1.0, sol_true );
    real64 const relTol = cond_est * params.krylov.relTolerance;
    EXPECT_LT( sol_comp.norm2() / sol_true.norm2(), relTol );
  }
};

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class KrylovSolverTest
  : public KrylovSolverTestBase< typename LAI::ParallelMatrix,
                                 PreconditionerIdentity< LAI >,
                                 typename LAI::ParallelVector >
{
public:
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

  using Base = KrylovSolverTestBase< typename LAI::ParallelMatrix,
                                     PreconditionerIdentity< LAI >,
                                     typename LAI::ParallelVector >;

  KrylovSolverTest() :
    Base()
  {}

protected:
  void
  SetUp()
  {
    // Compute matrix and preconditioner
    globalIndex constexpr n = 100;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, n, this->matrix );
    this->precond.compute( this->matrix );

    // Set up vectors
    this->sol_true.createWithGlobalSize( this->matrix.numGlobalCols(),
                                         MPI_COMM_GEOSX );
    this->sol_comp.createWithGlobalSize( this->matrix.numGlobalCols(),
                                         MPI_COMM_GEOSX );
    this->rhs_true.createWithGlobalSize( this->matrix.numGlobalRows(),
                                         MPI_COMM_GEOSX );

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 1.5 * 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( KrylovSolverTest );

TYPED_TEST_P( KrylovSolverTest, CG )
{
  this->test( params_CG() );
}

TYPED_TEST_P( KrylovSolverTest, BiCGSTAB )
{
  this->test( params_BiCGSTAB() );
}

TYPED_TEST_P( KrylovSolverTest, GMRES )
{
  this->test( params_GMRES() );
}

REGISTER_TYPED_TEST_SUITE_P( KrylovSolverTest, CG, BiCGSTAB, GMRES );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos, KrylovSolverTest, TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, KrylovSolverTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, KrylovSolverTest, PetscInterface, );
#endif

///////////////////////////////////////////////////////////////////////////////////////

template< typename LAI >
class KrylovSolverBlockTest
  : public KrylovSolverTestBase<
      BlockOperatorWrapper< typename LAI::ParallelVector, typename LAI::ParallelMatrix >,
      BlockOperatorWrapper< typename LAI::ParallelVector >,
      BlockVector< typename LAI::ParallelVector > >
{
public:
  using Matrix = typename LAI::ParallelMatrix;
  using Vector = typename LAI::ParallelVector;

  using Base = KrylovSolverTestBase<
    BlockOperatorWrapper< typename LAI::ParallelVector, typename LAI::ParallelMatrix >,
    BlockOperatorWrapper< typename LAI::ParallelVector >,
    BlockVector< typename LAI::ParallelVector > >;

  KrylovSolverBlockTest() :
    Base( 2, 2 )
  {}

protected:
  Matrix laplace2D;
  PreconditionerIdentity< LAI > identity;

  void
  SetUp()
  {
    globalIndex constexpr n = 100;
    compute2DLaplaceOperator( MPI_COMM_GEOSX, n, laplace2D );

    // We are going to assembly the following dummy system
    // [L 0] [x_true] = [b_0]
    // [0 L] [x_true] = [b_1]
    this->matrix.set( 0, 0, laplace2D );
    this->matrix.set( 1, 1, laplace2D );

    // Set up block identity preconditioning operator
    identity.compute( laplace2D );
    this->precond.set( 0, 0, identity );
    this->precond.set( 1, 1, identity );

    // Setup vectors
    this->sol_true.resize( 2 );
    this->sol_comp.resize( 2 );
    this->rhs_true.resize( 2 );

    for( localIndex i = 0; i < 2; ++i )
    {
      this->sol_true.block( i ).createWithGlobalSize( laplace2D.numGlobalCols(),
                                                      MPI_COMM_GEOSX );
      this->sol_comp.block( i ).createWithGlobalSize( laplace2D.numGlobalCols(),
                                                      MPI_COMM_GEOSX );
      this->rhs_true.block( i ).createWithGlobalSize( laplace2D.numGlobalRows(),
                                                      MPI_COMM_GEOSX );
    }

    // Condition number for the Laplacian matrix estimate: 4 * n^2 / pi^2
    this->cond_est = 3 * 4.0 * n * n / std::pow( M_PI, 2 );
  }
};

TYPED_TEST_SUITE_P( KrylovSolverBlockTest );

TYPED_TEST_P( KrylovSolverBlockTest, CG )
{
  this->test( params_CG() );
}

TYPED_TEST_P( KrylovSolverBlockTest, BiCGSTAB )
{
  this->test( params_BiCGSTAB() );
}

TYPED_TEST_P( KrylovSolverBlockTest, GMRES )
{
  this->test( params_GMRES() );
}

REGISTER_TYPED_TEST_SUITE_P( KrylovSolverBlockTest, CG, BiCGSTAB, GMRES );

#ifdef GEOSX_USE_TRILINOS
INSTANTIATE_TYPED_TEST_SUITE_P( Trilinos,
                                KrylovSolverBlockTest,
                                TrilinosInterface, );
#endif

#ifdef GEOSX_USE_HYPRE
INSTANTIATE_TYPED_TEST_SUITE_P( Hypre, KrylovSolverBlockTest, HypreInterface, );
#endif

#ifdef GEOSX_USE_PETSC
INSTANTIATE_TYPED_TEST_SUITE_P( Petsc, KrylovSolverBlockTest, PetscInterface, );
#endif

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
