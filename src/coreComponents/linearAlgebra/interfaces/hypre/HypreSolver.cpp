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
 * @file HypreSolver.cpp
 */

#include "HypreSolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <krylov.h>

namespace geosx
{

/**
 * @brief Container for hypre Krylov solver function pointers.
 */
struct HypreSolverWrapper : public HyprePrecWrapper
{
  /// Alias for set preconditioner function type
  using SetPrecondFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_Solver );

  /// Alias for get number of iterations function type
  using GetNumIter = HYPRE_Int ( * )( HYPRE_Solver solver,
                                      HYPRE_Int * num_iterations );

  /// Alias for get final residual norm function type
  using GetFinalNorm = HYPRE_Int ( * )( HYPRE_Solver solver,
                                        HYPRE_Real * norm );

  SetPrecondFunc setPrecond{}; ///< pointer to set preconditioner function
  GetNumIter getNumIter{};     ///< pointer to get number of iterations function
  GetFinalNorm getFinalNorm{}; ///< pointer to get final residual norm function
};

HypreSolver::HypreSolver( LinearSolverParameters parameters )
  : Base( std::move( parameters ) ),
  m_precond( m_params )
{ }

HypreSolver::~HypreSolver()
{
  HypreSolver::clear();
}

namespace
{

void createHypreGMRES( LinearSolverParameters const & params,
                       MPI_Comm const comm,
                       HypreSolverWrapper & solver )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESCreate( comm, &solver.ptr ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetKDim( solver.ptr, params.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver.ptr, params.logLevel ) ); // print iteration info
  //GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetPrintLevel( solver, 0 ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRGMRESSetLogging( solver.ptr, 1 ) ); /* needed to get run info later */

  solver.setPrecond = HYPRE_ParCSRGMRESSetPrecond;
  solver.setup = HYPRE_ParCSRGMRESSetup;
  solver.solve = HYPRE_ParCSRGMRESSolve;
  solver.getNumIter = HYPRE_GMRESGetNumIterations;
  solver.getFinalNorm = HYPRE_GMRESGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRGMRESDestroy;
}

void createHypreFlexGMRES( LinearSolverParameters const & params,
                           MPI_Comm const comm,
                           HypreSolverWrapper & solver )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESCreate( comm, &solver.ptr ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetKDim( solver.ptr, params.krylov.maxRestart ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetPrintLevel( solver.ptr, params.logLevel ) ); // print iteration info
  //GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetPrintLevel( solver, 0 ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRFlexGMRESSetLogging( solver.ptr, 1 ) ); /* needed to get run info later */

  solver.setPrecond = HYPRE_ParCSRFlexGMRESSetPrecond;
  solver.setup = HYPRE_ParCSRFlexGMRESSetup;
  solver.solve = HYPRE_ParCSRFlexGMRESSolve;
  solver.getNumIter = HYPRE_FlexGMRESGetNumIterations;
  solver.getFinalNorm = HYPRE_FlexGMRESGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRFlexGMRESDestroy;
}

void createHypreBiCGSTAB( LinearSolverParameters const & params,
                          MPI_Comm const comm,
                          HypreSolverWrapper & solver )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABCreate( comm, &solver.ptr ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetPrintLevel( solver.ptr, params.logLevel ) ); // print iteration info
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRBiCGSTABSetLogging( solver.ptr, 1 ) ); // needed to get run info later

  solver.setPrecond = HYPRE_ParCSRBiCGSTABSetPrecond;
  solver.setup = HYPRE_ParCSRBiCGSTABSetup;
  solver.solve = HYPRE_ParCSRBiCGSTABSolve;
  solver.getNumIter = HYPRE_BiCGSTABGetNumIterations;
  solver.getFinalNorm = HYPRE_BiCGSTABGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRBiCGSTABDestroy;
}

void createHypreCG( LinearSolverParameters const & params,
                    MPI_Comm const comm,
                    HypreSolverWrapper & solver )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParCSRPCGCreate( comm, &solver.ptr ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetMaxIter( solver.ptr, params.krylov.maxIterations ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTol( solver.ptr, params.krylov.relTolerance ) );

  // Default for now
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetPrintLevel( solver.ptr, params.logLevel ) ); /* print the iteration info */
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetLogging( solver.ptr, 1 ) );    /* needed to get run info later */
  GEOSX_LAI_CHECK_ERROR( HYPRE_PCGSetTwoNorm( solver.ptr, 1 ) );    /* use the two norm as the stopping criteria */

  solver.setPrecond = HYPRE_ParCSRPCGSetPrecond;
  solver.setup = HYPRE_ParCSRPCGSetup;
  solver.solve = HYPRE_ParCSRPCGSolve;
  solver.getNumIter = HYPRE_PCGGetNumIterations;
  solver.getFinalNorm = HYPRE_PCGGetFinalRelativeResidualNorm;
  solver.destroy = HYPRE_ParCSRPCGDestroy;
}

void createHypreKrylovSolver( LinearSolverParameters const & params,
                              MPI_Comm const comm,
                              HypreSolverWrapper & solver )
{
  switch( params.solverType )
  {
    case LinearSolverParameters::SolverType::gmres:
    {
      createHypreGMRES( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::fgmres:
    {
      createHypreFlexGMRES( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::bicgstab:
    {
      createHypreBiCGSTAB( params, comm, solver );
      break;
    }
    case LinearSolverParameters::SolverType::cg:
    {
      createHypreCG( params, comm, solver );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Solver type not supported in hypre interface: " << params.solverType );
    }
  }
}

} // namespace

void HypreSolver::setup( HypreMatrix const & mat )
{
  clear();
  Base::setup( mat );
  Stopwatch timer( m_result.setupTime );

  m_precond.setup( mat );
  m_componentFilterTime = m_precond.componentFilterTime();

  m_solver = std::make_unique< HypreSolverWrapper >();
  createHypreKrylovSolver( m_params, mat.getComm(), *m_solver );

  // Set the preconditioner
  GEOSX_LAI_CHECK_ERROR( m_solver->setPrecond( m_solver->ptr,
                                               m_precond.unwrapped().solve,
                                               hypre::DummySetup,
                                               m_precond.unwrapped().ptr ) );

  // Setup the solver (need a dummy vector for rhs/sol to avoid hypre segfaulting in setup)
  HypreVector dummy;
  dummy.createWithLocalSize( mat.numLocalRows(), mat.getComm() );
  GEOSX_LAI_CHECK_ERROR( m_solver->setup( m_solver->ptr,
                                          mat.unwrapped(),
                                          dummy.unwrapped(),
                                          dummy.unwrapped() ) );
}

int HypreSolver::doSolve( HypreVector const & rhs,
                          HypreVector & sol ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( sol.ready() );
  GEOSX_LAI_ASSERT( rhs.ready() );
  return m_solver->solve( m_solver->ptr, matrix().unwrapped(), rhs.unwrapped(), sol.unwrapped() );
}

void HypreSolver::apply( HypreVector const & rhs,
                         HypreVector & sol ) const
{
  doSolve( rhs, sol );
}

void HypreSolver::solve( HypreVector const & rhs,
                         HypreVector & sol ) const
{
  if( isZero( rhs.norm2(), 0.0 ) )
  {
    sol.zero();
    m_result.numIterations = 0;
    m_result.residualReduction = 0.0;
    m_result.solveTime = 0.0;
    m_result.status = LinearSolverResult::Status::Success;
    return;
  }

  {
    Stopwatch timer( m_result.solveTime );
    int const result = doSolve( rhs, sol );
    m_result.status = result ? LinearSolverResult::Status::NotConverged : LinearSolverResult::Status::Success;
  }

  // Clear error code to avoid GEOSX from crashing if Krylov method did not converge
  GEOSX_LAI_CHECK_ERROR( HYPRE_ClearAllErrors() );

  // Get final residual norm
  GEOSX_LAI_CHECK_ERROR( m_solver->getFinalNorm( m_solver->ptr, &m_result.residualReduction ) );

  // Get number of iterations
  HYPRE_Int numIter;
  GEOSX_LAI_CHECK_ERROR( m_solver->getNumIter( m_solver->ptr, &numIter ) );
  m_result.numIterations = numIter;

  if( m_params.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( "\t\tLinear Solver | " << m_result.status <<
                      " | Iterations: " << m_result.numIterations <<
                      " | Final Rel Res: " << m_result.residualReduction <<
                      " | SC Filter Time: " << m_componentFilterTime <<
                      " | Setup Time: " << m_result.setupTime << " s" <<
                      " | Solve Time: " << m_result.solveTime << " s" );
  }
}

void HypreSolver::clear()
{
  Base::clear();
  if( m_solver )
  {
    GEOSX_LAI_CHECK_ERROR( m_solver->destroy( m_solver->ptr ) );
    m_solver = nullptr;
  }
  m_solver.reset();
}

} // end geosx namespace
