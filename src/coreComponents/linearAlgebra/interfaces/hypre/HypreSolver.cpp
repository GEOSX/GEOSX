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
 * @file HypreSolver.cpp
 */

#include "HypreSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/interfaces/hypre/HyprePreconditioner.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <krylov.h>

namespace geosx
{
typedef HYPRE_Int ( *HYPRE_PtrToSolverDestroyFcn )( HYPRE_Solver );

HypreSolver::HypreSolver( LinearSolverParameters parameters ) :
  m_parameters( std::move( parameters ) )
{}

//// ----------------------------
//// Top-Level Solver
//// ----------------------------
//// We switch between different solverTypes here
//
void
HypreSolver::solve( HypreMatrix & mat,
                    HypreVector & sol,
                    HypreVector & rhs,
                    DofManager const * const dofManager )
{
  GEOSX_LAI_ASSERT( mat.ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );

  if( m_parameters.solverType == "direct" )
  {
    solve_direct( mat, sol, rhs );
  }
  else
  {
    solve_krylov( mat, sol, rhs, dofManager );
  }
}

void
HypreSolver::solve_direct( HypreMatrix & mat, HypreVector & sol, HypreVector & rhs )
{
  // Instantiate solver
  HYPRE_Solver solver;

  Stopwatch watch;
  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistSetup( &solver, mat.unwrapped(), 0 ) );
  m_result.setupTime = watch.elapsedTime();
  watch.zero();
  GEOSX_LAI_CHECK_ERROR(
    hypre_SLUDistSolve( solver, rhs.unwrapped(), sol.unwrapped() ) );
  m_result.solveTime = watch.elapsedTime();

  m_result.status = LinearSolverResult::Status::Success;
  m_result.numIterations = 1;
  m_result.residualReduction = NumericTraits< real64 >::eps;

  GEOSX_LAI_CHECK_ERROR( hypre_SLUDistDestroy( solver ) );
}

namespace
{
void
CreateHypreGMRES( LinearSolverParameters const & params,
                  MPI_Comm const comm,
                  HYPRE_Solver & solver,
                  HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRGMRESSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRGMRESSetKDim( solver, params.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRGMRESSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRGMRESSetPrintLevel( solver, params.logLevel ) );  // print iteration info
  //GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver, 0 ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRGMRESSetLogging( solver, 1 ) ); /* needed to get run info later */

  solverFuncs.setPrecond = HYPRE_ParCSRGMRESSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRGMRESSetup;
  solverFuncs.solve = HYPRE_ParCSRGMRESSolve;
  solverFuncs.getNumIter = HYPRE_GMRESGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_GMRESGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRGMRESDestroy;
}

void
CreateHypreFlexGMRES( LinearSolverParameters const & params,
                      MPI_Comm const comm,
                      HYPRE_Solver & solver,
                      HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRFlexGMRESSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRFlexGMRESSetKDim( solver, params.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRFlexGMRESSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetPrintLevel(
    solver,
    params.logLevel ) );  // print iteration info
  //GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetPrintLevel( solver, 0 ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRFlexGMRESSetLogging( solver, 1 ) ); /* needed to get run info later */

  solverFuncs.setPrecond = HYPRE_ParCSRFlexGMRESSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRFlexGMRESSetup;
  solverFuncs.solve = HYPRE_ParCSRFlexGMRESSolve;
  solverFuncs.getNumIter = HYPRE_FlexGMRESGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_FlexGMRESGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRFlexGMRESDestroy;
}

void
CreateHypreBiCGSTAB( LinearSolverParameters const & params,
                     MPI_Comm const comm,
                     HYPRE_Solver & solver,
                     HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRBiCGSTABSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRBiCGSTABSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrintLevel(
    solver,
    params.logLevel ) );  // print iteration info
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParCSRBiCGSTABSetLogging( solver, 1 ) );  // needed to get run info later

  solverFuncs.setPrecond = HYPRE_ParCSRBiCGSTABSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRBiCGSTABSetup;
  solverFuncs.solve = HYPRE_ParCSRBiCGSTABSolve;
  solverFuncs.getNumIter = HYPRE_BiCGSTABGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_BiCGSTABGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRBiCGSTABDestroy;
}

void
CreateHypreCG( LinearSolverParameters const & params,
               MPI_Comm const comm,
               HYPRE_Solver & solver,
               HypreSolverFuncs & solverFuncs )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( comm, &solver ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( solver, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTol( solver, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_PCGSetPrintLevel( solver, params.logLevel ) ); /* print the iteration info */
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_PCGSetLogging( solver, 1 ) ); /* needed to get run info later */
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_PCGSetTwoNorm( solver, 1 ) ); /* use the two norm as the stopping criteria */

  solverFuncs.setPrecond = HYPRE_ParCSRPCGSetPrecond;
  solverFuncs.setup = HYPRE_ParCSRPCGSetup;
  solverFuncs.solve = HYPRE_ParCSRPCGSolve;
  solverFuncs.getNumIter = HYPRE_PCGGetNumIterations;
  solverFuncs.getFinalNorm = HYPRE_PCGGetFinalRelativeResidualNorm;
  solverFuncs.destroy = HYPRE_ParCSRPCGDestroy;
}

void
CreateHypreKrylovSolver( LinearSolverParameters const & params,
                         MPI_Comm const comm,
                         HYPRE_Solver & solver,
                         HypreSolverFuncs & solverFuncs )
{
  if( params.solverType == "gmres" )
  {
    CreateHypreGMRES( params, comm, solver, solverFuncs );
  }
  else if( params.solverType == "fgmres" )
  {
    CreateHypreFlexGMRES( params, comm, solver, solverFuncs );
  }
  else if( params.solverType == "bicgstab" )
  {
    CreateHypreBiCGSTAB( params, comm, solver, solverFuncs );
  }
  else if( params.solverType == "cg" )
  {
    CreateHypreCG( params, comm, solver, solverFuncs );
  }
  else
  {
    GEOSX_ERROR( "Unsupported Hypre solver type: " << params.solverType );
  }
}

}  // namespace

void
HypreSolver::solve_krylov( HypreMatrix & mat,
                           HypreVector & sol,
                           HypreVector & rhs,
                           DofManager const * const dofManager )
{
  Stopwatch watch;

  // Create the preconditioner, but don't compute (this is done by solver setup)
  HyprePreconditioner precond( m_parameters, dofManager );

  // Deal with separate component approximation
  // TODO: preliminary version for separate displacement components
  HypreMatrix separateComponentMatrix;
  HYPRE_Solver uu_amg_solver = {};  //TODO: this is a quick and dirty first implementation

  if( m_parameters.amg.separateComponents &&
      m_parameters.preconditionerType != "mgr" )
  {
    LAIHelperFunctions::SeparateComponentFilter( mat,
                                                 separateComponentMatrix,
                                                 m_parameters.dofsPerNode );
  }
  else if( m_parameters.preconditionerType == "mgr" &&
           m_parameters.mgr.separateComponents )
  {
    // Extract displacement block
    HypreMatrix Pu;
    HypreMatrix scr_mat;
    dofManager->makeRestrictor( { { m_parameters.mgr.displacementFieldName, 0, 3 } },
                                mat.getComm(),
                                true,
                                Pu );
    mat.multiplyPtAP( Pu, scr_mat );
    LAIHelperFunctions::SeparateComponentFilter( scr_mat,
                                                 separateComponentMatrix,
                                                 m_parameters.dofsPerNode );

    HYPRE_BoomerAMGCreate( &uu_amg_solver );
    HYPRE_BoomerAMGSetTol( uu_amg_solver, 0.0 );
    HYPRE_BoomerAMGSetMaxIter( uu_amg_solver, 1 );
    HYPRE_BoomerAMGSetPrintLevel( uu_amg_solver, 0 );
    HYPRE_BoomerAMGSetRelaxOrder( uu_amg_solver, 1 );
    HYPRE_BoomerAMGSetAggNumLevels( uu_amg_solver, 1 );
    HYPRE_BoomerAMGSetNumFunctions( uu_amg_solver, 3 );

    HYPRE_BoomerAMGSetup( uu_amg_solver,
                          separateComponentMatrix.unwrapped(),
                          nullptr,
                          nullptr );

    HYPRE_MGRSetFSolver( precond.unwrapped(),
                         HYPRE_BoomerAMGSolve,
                         HYPRE_BoomerAMGSetup,
                         uu_amg_solver );
  }
  HypreMatrix & precondMat =
    m_parameters.amg.separateComponents ? separateComponentMatrix : mat;

  // Instantiate the solver
  HYPRE_Solver solver {};
  HypreSolverFuncs solverFuncs;
  CreateHypreKrylovSolver( m_parameters, mat.getComm(), solver, solverFuncs );

  // Set the preconditioner
  GEOSX_LAI_CHECK_ERROR( solverFuncs.setPrecond( solver,
                                                 precond.unwrappedFuncs().apply,
                                                 precond.unwrappedFuncs().setup,
                                                 precond.unwrapped() ) );

  // Setup
  GEOSX_LAI_CHECK_ERROR( solverFuncs.setup( solver,
                                            precondMat.unwrapped(),
                                            rhs.unwrapped(),
                                            sol.unwrapped() ) );
  m_result.setupTime = watch.elapsedTime();

  // Solve
  watch.zero();
  HYPRE_Int const result =
    solverFuncs.solve( solver, mat.unwrapped(), rhs.unwrapped(), sol.unwrapped() );
  m_result.solveTime = watch.elapsedTime();

  // Set result status based on return value
  m_result.status = result ? LinearSolverResult::Status::NotConverged
                           : LinearSolverResult::Status::Success;

  // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
  GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

  // Get final residual norm
  HYPRE_Real finalNorm;
  GEOSX_LAI_CHECK_ERROR( solverFuncs.getFinalNorm( solver, &finalNorm ) );
  m_result.residualReduction = finalNorm;

  // Get number of iterations
  HYPRE_Int numIter;
  GEOSX_LAI_CHECK_ERROR( solverFuncs.getNumIter( solver, &numIter ) );
  m_result.numIterations = numIter;

  if( m_parameters.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( "\t\tLinear Solver | Iter = "
                      << numIter << " | Final Relative Tol " << finalNorm
                      << " | SetupTime " << m_result.setupTime << " | SolveTime "
                      << m_result.solveTime );
  }

  // Destroy solver
  GEOSX_LAI_CHECK_ERROR( solverFuncs.destroy( solver ) );
  if( m_parameters.preconditionerType == "mgr" &&
      m_parameters.mgr.separateComponents )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGDestroy( uu_amg_solver ) );
  }
}

}  // namespace geosx
