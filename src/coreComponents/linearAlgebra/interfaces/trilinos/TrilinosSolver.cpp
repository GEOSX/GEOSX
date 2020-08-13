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
 * @file TrilinosSolver.cpp
 */

#include "TrilinosSolver.hpp"
#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraMatrix.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosPreconditioner.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <AztecOO.h>
#include <Amesos.h>
#include <ml_MultiLevelPreconditioner.h>

#include <memory>

namespace geosx
{
TrilinosSolver::TrilinosSolver( LinearSolverParameters parameters ) :
  m_parameters( std::move( parameters ) )
{}

TrilinosSolver::~TrilinosSolver() = default;

void
TrilinosSolver::solve( EpetraMatrix & mat,
                       EpetraVector & sol,
                       EpetraVector & rhs,
                       DofManager const * const dofManager )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  GEOSX_UNUSED_VAR( dofManager );

  if( m_parameters.scaling.useRowScaling )
  {
    Epetra_FECrsMatrix & mat_raw = mat.unwrapped();
    Epetra_MultiVector & rhs_raw = rhs.unwrapped();

    Epetra_Vector scaling( mat_raw.RowMap() );
    mat_raw.InvRowSums( scaling );
    mat_raw.LeftScale( scaling );

    Epetra_MultiVector tmp( rhs_raw );
    rhs_raw.Multiply( 1.0, scaling, tmp, 0.0 );
  }

  if( m_parameters.solverType == "direct" )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs );
  }
}

void
TrilinosSolver::solve_direct( EpetraMatrix & mat,
                              EpetraVector & sol,
                              EpetraVector & rhs )
{
  // Time setup and solve
  Stopwatch watch;

  // Create Epetra linear problem and instantiate solver.
  Epetra_LinearProblem problem( &mat.unwrapped(),
                                &sol.unwrapped(),
                                &rhs.unwrapped() );

  // Instantiate the Amesos solver.
  Amesos Factory;

  // Select KLU solver (only one available as of 9/20/2018)
  std::unique_ptr< Amesos_BaseSolver > solver( Factory.Create( "Klu", problem ) );

  // Factorize the matrix
  GEOSX_LAI_CHECK_ERROR( solver->SymbolicFactorization() );
  GEOSX_LAI_CHECK_ERROR( solver->NumericFactorization() );
  m_result.setupTime = watch.elapsedTime();

  // Solve the system
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( solver->Solve() );
  m_result.solveTime = watch.elapsedTime();

  // Basic output
  if( m_parameters.logLevel > 0 )
  {
    solver->PrintStatus();
    solver->PrintTiming();
  }

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;
  m_result.residualReduction = NumericTraits< real64 >::eps;
}

namespace
{
void
CreateTrilinosKrylovSolver( LinearSolverParameters const & params,
                            AztecOO & solver )
{
  // Choose the solver type
  if( params.solverType == "gmres" )
  {
    GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_gmres ) );
    GEOSX_LAI_CHECK_ERROR(
      solver.SetAztecOption( AZ_kspace, params.krylov.maxRestart ) );
  }
  else if( params.solverType == "bicgstab" )
  {
    GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_bicgstab ) );
  }
  else if( params.solverType == "cg" )
  {
    GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_solver, AZ_cg ) );
  }
  else
  {
    GEOSX_ERROR( "Unsupported linear solver type: " << params.solverType );
  }

  // Ask for a convergence normalized by the right hand side
  GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_conv, AZ_rhs ) );

  // Control output
  switch( params.logLevel )
  {
    case 1:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_summary ) );
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_diagnostics, AZ_all ) );
      break;
    }
    case 2:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_all ) );
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_diagnostics, AZ_all ) );
      break;
    }
    default:
    {
      GEOSX_LAI_CHECK_ERROR( solver.SetAztecOption( AZ_output, AZ_none ) );
    }
  }
}

}  // namespace

void
TrilinosSolver::solve_krylov( EpetraMatrix & mat,
                              EpetraVector & sol,
                              EpetraVector & rhs )
{
  // Time setup and solve
  Stopwatch watch;

  // Create Epetra linear problem.
  Epetra_LinearProblem problem( &mat.unwrapped(),
                                &sol.unwrapped(),
                                &rhs.unwrapped() );

  // Instantiate the AztecOO solver.
  AztecOO solver( problem );

  // Choose the solver type
  CreateTrilinosKrylovSolver( m_parameters, solver );

  // Deal with separate component approximation
  EpetraMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat,
                                                 separateComponentMatrix,
                                                 m_parameters.dofsPerNode );
  }
  EpetraMatrix & precondMat =
    m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // Create and compute preconditioner
  TrilinosPreconditioner precond( m_parameters );
  precond.compute( precondMat );
  GEOSX_LAI_CHECK_ERROR( solver.SetPrecOperator( &precond.unwrapped() ) );

  m_result.setupTime = watch.elapsedTime();

  // Actually solve
  watch.zero();
  int const result = solver.Iterate( m_parameters.krylov.maxIterations,
                                     m_parameters.krylov.relTolerance );
  m_result.solveTime = watch.elapsedTime();

  m_result.status = result ? LinearSolverResult::Status::NotConverged
                           : LinearSolverResult::Status::Success;
  m_result.numIterations = solver.NumIters();
  m_result.residualReduction = solver.ScaledResidual();
}

}  // namespace geosx
