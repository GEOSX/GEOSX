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
 * @file PetscSolver.cpp
 */

#include "PetscSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"
#include "linearAlgebra/interfaces/petsc/PetscPreconditioner.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>

// Put everything under the geosx namespace.
namespace geosx
{
PetscSolver::PetscSolver( LinearSolverParameters parameters ) :
  m_parameters( std::move( parameters ) )
{}

void
PetscSolver::solve( PetscMatrix & mat,
                    PetscVector & sol,
                    PetscVector & rhs,
                    DofManager const * const dofManager )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  GEOSX_UNUSED_VAR( dofManager );

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
PetscSolver::solve_direct( PetscMatrix & mat, PetscVector & sol, PetscVector & rhs )
{
  MPI_Comm const comm = mat.getComm();

  // create linear solver
  KSP ksp;
  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetOperators( ksp, mat.unwrapped(), mat.unwrapped() ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPPREONLY ) );

  // use direct solve preconditioner SUPERLU DIST
  Stopwatch watch;
  PC prec;
  GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp, &prec ) );
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCLU ) );
  GEOSX_LAI_CHECK_ERROR( PCFactorSetMatSolverType( prec, MATSOLVERSUPERLU_DIST ) );
  GEOSX_LAI_CHECK_ERROR( PCSetUp( prec ) );
  m_result.setupTime = watch.elapsedTime();

  // solve system
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );
  m_result.solveTime = watch.elapsedTime();

  KSPConvergedReason reason;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &reason ) );

  m_result.status = reason >= 0 ? LinearSolverResult::Status::Success
                                : LinearSolverResult::Status::Breakdown;
  m_result.numIterations = 1.0;
  m_result.residualReduction = NumericTraits< real64 >::eps;

  // destroy solver
  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
}

namespace
{
void
CreatePetscKrylovSolver( LinearSolverParameters const & params,
                         MPI_Comm const comm,
                         KSP & ksp )
{
  GEOSX_LAI_CHECK_ERROR( KSPCreate( comm, &ksp ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetNormType( ksp, KSP_NORM_UNPRECONDITIONED ) );
  GEOSX_LAI_CHECK_ERROR( KSPSetTolerances( ksp,
                                           params.krylov.relTolerance,
                                           PETSC_DEFAULT,
                                           PETSC_DEFAULT,
                                           params.krylov.maxIterations ) );

  // pick the solver type
  if( params.solverType == "gmres" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPGMRES ) );
    GEOSX_LAI_CHECK_ERROR( KSPGMRESSetRestart( ksp, params.krylov.maxRestart ) );
  }
  else if( params.solverType == "bicgstab" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPBCGS ) );
  }
  else if( params.solverType == "cg" )
  {
    GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp, KSPCG ) );
  }
  else
  {
    GEOSX_ERROR( "Unsupported Petsc solver type: " << params.solverType );
  }
}

}  // namespace

void
PetscSolver::solve_krylov( PetscMatrix & mat, PetscVector & sol, PetscVector & rhs )
{
  Stopwatch watch;

  // create linear solver
  KSP ksp;
  CreatePetscKrylovSolver( m_parameters, mat.getComm(), ksp );

  // This can be used to extract residual norm history:
  //array1d< real64 > residualNorms( m_parameters.krylov.maxIterations + 1 );
  //GEOSX_LAI_CHECK_ERROR( KSPSetResidualHistory( ksp, residualNorms.data(), residualNorms.size(), PETSC_TRUE ) );

  // Deal with separate component approximation
  PetscMatrix separateComponentMatrix;
  if( m_parameters.amg.separateComponents )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat,
                                                 separateComponentMatrix,
                                                 m_parameters.dofsPerNode );
  }
  PetscMatrix & precondMat =
    m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // create and compute a preconditioner and set into KSP
  PetscPreconditioner precond( m_parameters );
  precond.compute( precondMat );
  KSPSetPC( ksp, precond.unwrapped() );

  GEOSX_LAI_CHECK_ERROR(
    KSPSetOperators( ksp, mat.unwrapped(), precondMat.unwrapped() ) );

  // display output
  if( m_parameters.logLevel > 0 )
  {
    GEOSX_LAI_CHECK_ERROR( PetscOptionsSetValue( nullptr, "-ksp_monitor", nullptr ) );
  }
  GEOSX_LAI_CHECK_ERROR( KSPSetFromOptions( ksp ) );

  m_result.setupTime = watch.elapsedTime();

  // Actually solve
  watch.zero();
  GEOSX_LAI_CHECK_ERROR( KSPSolve( ksp, rhs.unwrapped(), sol.unwrapped() ) );
  m_result.solveTime = watch.elapsedTime();

  // Get status indicator
  KSPConvergedReason result;
  GEOSX_LAI_CHECK_ERROR( KSPGetConvergedReason( ksp, &result ) );
  GEOSX_WARNING_IF( result < 0, "PetscSolver: Krylov convergence not achieved" );

  switch( result )
  {
    case KSP_CONVERGED_RTOL:
    case KSP_CONVERGED_ATOL:
    case KSP_CONVERGED_ITS:
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }
    case KSP_DIVERGED_ITS:
    case KSP_DIVERGED_DTOL:
    {
      m_result.status = LinearSolverResult::Status::NotConverged;
      break;
    }
    default:
    {
      m_result.status = LinearSolverResult::Status::Breakdown;
    }
  }

  // Get number of iterations performed
  PetscInt numIter;
  GEOSX_LAI_CHECK_ERROR( KSPGetIterationNumber( ksp, &numIter ) );
  m_result.numIterations = numIter;

  // reset verbosity option
  GEOSX_LAI_CHECK_ERROR( PetscOptionsClearValue( nullptr, "-ksp_monitor" ) );

  GEOSX_LAI_CHECK_ERROR( KSPDestroy( &ksp ) );
}

}  // namespace geosx
