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
 * @file GMRESsolver.cpp
 */

#include "GMRESsolver.hpp"

#include "common/Stopwatch.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/solvers/KrylovUtils.hpp"

namespace geosx
{
template< typename VECTOR >
GMRESsolver< VECTOR >::GMRESsolver( LinearOperator< Vector > const & A,
                                    LinearOperator< Vector > const & M,
                                    real64 tolerance,
                                    localIndex maxIterations,
                                    integer verbosity,
                                    localIndex maxRestart ) :
  KrylovSolver< VECTOR >( A, M, tolerance, maxIterations, verbosity ),
  m_maxRestart( maxRestart ),
  m_kspace( m_maxRestart + 1 ),
  m_kspaceInitialized( false )
{
  GEOSX_ERROR_IF_LE_MSG(
    m_maxRestart,
    0,
    "GMRES: max number of restart iterations must be positive." );
}

template< typename VECTOR >
GMRESsolver< VECTOR >::~GMRESsolver() = default;

namespace
{
void
ComputeGivensRotation( real64 const x, real64 const y, real64 & c, real64 & s )
{
  if( isZero( y ) )
  {
    c = 1.0;
    s = 0.0;
  }
  else if( std::fabs( y ) > std::fabs( x ) )
  {
    real64 const nu = x / y;
    s = 1.0 / std::sqrt( 1.0 + nu * nu );
    c = nu * s;
  }
  else
  {
    real64 const nu = y / x;
    c = 1.0 / std::sqrt( 1.0 + nu * nu );
    s = nu * c;
  }
}

void
ApplyGivensRotation( real64 const c, real64 const s, real64 & dx, real64 & dy )
{
  real64 const temp = c * dx + s * dy;
  dy = -s * dx + c * dy;
  dx = temp;
}

void
Backsolve( localIndex const k,
           arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & H,
           arraySlice1d< real64 > const & g )
{
  for( localIndex j = k - 1; j >= 0; --j )
  {
    g[j] /= H( j, j );
    for( localIndex i = j - 1; i >= 0; --i )
    {
      g[i] -= H( i, j ) * g[j];
    }
  }
}

}  // namespace

template< typename VECTOR >
void
GMRESsolver< VECTOR >::solve( Vector const & b, Vector & x ) const
{
  // We create Krylov subspace vectors once using the size and partitioning of b.
  // It is assumed that on every repeated call to solve() input vectors will keep
  // the same (or at least compatible) size and partitioning.
  if( !m_kspaceInitialized )
  {
    for( localIndex i = 0; i < m_maxRestart + 1; ++i )
    {
      m_kspace[i] = createTempVector( b );
    }
  }

  Stopwatch watch;

  // Compute the target absolute tolerance
  real64 const absTol = b.norm2() * m_tolerance;

  // Define vectors
  VectorTemp r = createTempVector( b );
  VectorTemp w = createTempVector( b );
  VectorTemp z = createTempVector( b );

  // Compute initial rk
  m_operator.residual( x, b, r );

  // Create upper Hessenberg matrix
  array2d< real64, MatrixLayout::COL_MAJOR_PERM > H( m_maxRestart + 1, m_maxRestart );

  // Create plane rotation storage
  array1d< real64 > c( m_maxRestart + 1 );
  array1d< real64 > s( m_maxRestart + 1 );
  array1d< real64 > g( m_maxRestart + 1 );

  m_result.status = LinearSolverResult::Status::NotConverged;
  m_residualNorms.resize( m_maxIterations + 1 );

  localIndex k;
  real64 rnorm = 0.0;

  for( k = 0; k <= m_maxIterations &&
       m_result.status == LinearSolverResult::Status::NotConverged; )
  {
    // Re-initialize Krylov subspace
    g.setValues< serialPolicy >( 0.0 );
    g[0] = r.norm2();
    m_kspace[0].axpby( 1.0 / g[0], r, 0.0 );

    localIndex j;
    for( j = 0; j < m_maxRestart && k <= m_maxIterations; ++j, ++k )
    {
      // Record iteration progress
      rnorm = std::fabs( g[j] );
      logProgress( k, rnorm );

      // Convergence check
      if( rnorm < absTol )
      {
        m_result.status = LinearSolverResult::Status::Success;
        break;
      }

      // Compute the new vector
      m_precond.apply( m_kspace[j], z );
      m_operator.apply( z, w );

      // Orthogonalization
      for( localIndex i = 0; i <= j; ++i )
      {
        H( i, j ) = w.dot( m_kspace[i] );
        w.axpby( -H( i, j ), m_kspace[i], 1.0 );
      }

      H( j + 1, j ) = w.norm2();
      GEOSX_KRYLOV_BREAKDOWN_IF_ZERO( H( j + 1, j ) );
      m_kspace[j + 1].axpby( 1.0 / H( j + 1, j ), w, 0.0 );

      // Apply all previous rotations to the new column
      for( localIndex i = 0; i < j; ++i )
      {
        ApplyGivensRotation( c[i], s[i], H( i, j ), H( i + 1, j ) );
      }

      // Compute and apply the new rotation to eliminate subdiagonal element
      ComputeGivensRotation( H( j, j ), H( j + 1, j ), c[j], s[j] );
      ApplyGivensRotation( c[j], s[j], H( j, j ), H( j + 1, j ) );
      ApplyGivensRotation( c[j], s[j], g[j], g[j + 1] );
    }

    // Regardless of how we quit out of inner loop, j is the actual size of H
    Backsolve( j, H, g );
    w.zero();
    for( localIndex i = 0; i < j; ++i )
    {
      w.axpy( g[i], m_kspace[i] );
    }
    m_precond.apply( w, z );

    // Update the solution vector and recompute residual
    x.axpy( 1.0, z );
    m_operator.residual( x, b, r );
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
template class GMRESsolver< TrilinosInterface::ParallelVector >;
template class GMRESsolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class GMRESsolver< HypreInterface::ParallelVector >;
template class GMRESsolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class GMRESsolver< PetscInterface::ParallelVector >;
template class GMRESsolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}  // namespace geosx
