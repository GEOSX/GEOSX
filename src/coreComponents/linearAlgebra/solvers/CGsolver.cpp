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
 * @file CGsolver.cpp
 *
 */

#include "CGsolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/utilities/BlockVectorView.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geosx
{
// BEGIN_RST_NARRATIVE CGsolver.rst
// ==============================
// CG Solver
// ==============================
// Implementation of the Preconditioned Conjugate Gradient (CG) algorithm.
// The notation is consistent with "Iterative Methods for
// Linear and Non-Linear Equations" from C.T. Kelley (1995)
// and "Iterative Methods for Sparse Linear Systems" from Y. Saad (2003).

// ----------------------------
// Constructor
// ----------------------------
// Empty constructor.
template< typename VECTOR >
CGsolver< VECTOR >::CGsolver( LinearOperator< Vector > const & A,
                              LinearOperator< Vector > const & M,
                              real64 const tolerance,
                              localIndex const maxIterations,
                              integer const verbosity ) :
  KrylovSolver< VECTOR >( A, M, tolerance, maxIterations, verbosity )
{}

// ----------------------------
// Destructor
// ----------------------------
template< typename VECTOR >
CGsolver< VECTOR >::~CGsolver() = default;

// ----------------------------
// Monolithic CG solver
// ----------------------------
template< typename VECTOR >
void
CGsolver< VECTOR >::solve( Vector const & b, Vector & x ) const

{
  Stopwatch watch;

  // Compute the target absolute tolerance
  real64 const absTol = b.norm2() * m_tolerance;

  // Define residual vector
  VectorTemp r = createTempVector( b );

  // Compute initial rk =  b - Ax
  m_operator.residual( x, b, r );

  // Preconditioning
  VectorTemp z = createTempVector( x );

  // Search direction
  VectorTemp p = createTempVector( z );
  VectorTemp Ap = createTempVector( z );

  // Keep old value of preconditioned residual norm
  real64 tau_old = 0.0;

  p.zero();
  m_result.status = LinearSolverResult::Status::NotConverged;
  m_result.numIterations = 0;
  m_residualNorms.resize( m_maxIterations + 1 );

  localIndex k;
  real64 rnorm = 0.0;

  for( k = 0; k <= m_maxIterations; ++k )
  {
    rnorm = r.norm2();
    logProgress( k, rnorm );

    // Convergence check on ||rk||/||b||
    if( rnorm < absTol )
    {
      m_result.status = LinearSolverResult::Status::Success;
      break;
    }

    // Update z = Mr
    m_precond.apply( r, z );

    // Compute beta
    real64 const tau = z.dot( r );
    real64 const beta = k > 0 ? tau / tau_old : 0.0;

    // Update p = z + beta*p
    p.axpby( 1.0, z, beta );

    // Compute Ap
    m_operator.apply( p, Ap );

    // compute alpha
    real64 const pAp = p.dot( Ap );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( pAp );
    real64 const alpha = tau / pAp;

    // Update x = x + alpha*p
    x.axpby( alpha, p, 1.0 );

    // Update rk = rk - alpha*Ap
    r.axpby( -alpha, Ap, 1.0 );

    // Keep the old tau value
    tau_old = tau;
  }

  m_result.numIterations = k;
  m_result.residualReduction = rnorm / absTol * m_tolerance;
  m_result.solveTime = watch.elapsedTime();

  logResult();
  m_residualNorms.resize( m_result.numIterations + 1 );
}

// END_RST_NARRATIVE

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class CGsolver< TrilinosInterface::ParallelVector >;
template class CGsolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class CGsolver< HypreInterface::ParallelVector >;
template class CGsolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class CGsolver< PetscInterface::ParallelVector >;
template class CGsolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}  //namespace geosx
