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
 * @file BiCGSTABsolver.cpp
 *
 */

#include "BiCGSTABsolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geosx
{
template< typename VECTOR >
BiCGSTABsolver< VECTOR >::BiCGSTABsolver( LinearOperator< Vector > const & A,
                                          LinearOperator< Vector > const & M,
                                          real64 const tolerance,
                                          localIndex const maxIterations,
                                          integer const verbosity ) :
  KrylovSolver< VECTOR >( A, M, tolerance, maxIterations, verbosity )
{}

template< typename VECTOR >
BiCGSTABsolver< VECTOR >::~BiCGSTABsolver() = default;

template< typename VECTOR >
void
BiCGSTABsolver< VECTOR >::solve( Vector const & b, Vector & x ) const
{
  Stopwatch watch;

  // Compute the target absolute tolerance
  real64 const absTol = b.norm2() * m_tolerance;

  // Define vectors
  VectorTemp r( x );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Define vectors
  VectorTemp r0( r );

  // Define scalars and reinitialize some
  real64 rho_old = r.dot( r0 );
  real64 alpha = 1.0;
  real64 omega = 1.0;

  // Define temporary vectors
  VectorTemp v = createTempVector( r );
  VectorTemp p = createTempVector( r );
  VectorTemp y = createTempVector( r );
  VectorTemp z = createTempVector( r );
  VectorTemp t = createTempVector( r );
  VectorTemp s = createTempVector( r );
  VectorTemp q = createTempVector( x );

  v.zero();
  p.zero();
  m_result.status = LinearSolverResult::Status::NotConverged;
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

    // Compute r0.rk
    real64 const rho = r.dot( r0 );

    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( rho_old );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( omega );

    // Compute beta
    real64 const beta = rho / rho_old * alpha / omega;

    // Update p = r + beta*(p - omega*v)
    p.axpy( -omega, v );
    p.axpby( 1., r, beta );

    // Update vk = MApk
    m_precond.apply( p, y );
    m_operator.apply( y, v );

    // Compute alpha
    real64 const vr0 = v.dot( r0 );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( vr0 );
    alpha = rho / vr0;

    // compute x = x + alpha*y
    x.axpy( alpha, y );

    // Compute s = rk - alpha*vk
    s.copy( r );
    s.axpy( -alpha, v );

    // Compute z = Ms
    m_precond.apply( s, z );

    // Compute t = Az
    m_operator.apply( z, t );

    // Compute t = Mt
    m_precond.apply( t, q );

    // Update omega
    real64 const q2 = q.dot( q );
    GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( q2 );
    omega = q.dot( z ) / q2;

    // Update x = x + omega*z
    x.axpy( omega, z );

    // Update rk = s - omega*t
    s.axpy( -omega, t );
    r.copy( s );

    // Keep the old value of rho
    rho_old = rho;
  }

  m_result.numIterations = k;
  m_result.residualReduction = rnorm / absTol * m_tolerance;
  m_result.solveTime = watch.elapsedTime();

  logResult();
  m_residualNorms.resize( m_result.numIterations + 1 );
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BiCGSTABsolver< TrilinosInterface::ParallelVector >;
template class BiCGSTABsolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class BiCGSTABsolver< HypreInterface::ParallelVector >;
template class BiCGSTABsolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class BiCGSTABsolver< PetscInterface::ParallelVector >;
template class BiCGSTABsolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}  //namespace geosx
