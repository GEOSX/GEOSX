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
 * @file KrylovSolver.cpp
 */

#include "KrylovSolver.hpp"
#include "linearAlgebra/solvers/BiCGSTABsolver.hpp"
#include "linearAlgebra/solvers/CGsolver.hpp"
#include "linearAlgebra/solvers/GMRESsolver.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{
template< typename VECTOR >
KrylovSolver< VECTOR >::KrylovSolver( LinearOperator< Vector > const & matrix,
                                      LinearOperator< Vector > const & precond,
                                      real64 const tolerance,
                                      localIndex const maxIterations,
                                      integer const verbosity ) :
  Base(),
  m_operator( matrix ),
  m_precond( precond ),
  m_tolerance( tolerance ),
  m_maxIterations( maxIterations ),
  m_logLevel( verbosity )
{
  GEOSX_ERROR_IF_LE_MSG(
    m_maxIterations,
    0,
    "Krylov solver: max number of iteration must be positive." );
  GEOSX_LAI_ASSERT_EQ( m_operator.numGlobalRows(), m_precond.numGlobalRows() );
  GEOSX_LAI_ASSERT_EQ( m_operator.numGlobalCols(), m_precond.numGlobalCols() );
}

template< typename VECTOR >
std::unique_ptr< KrylovSolver< VECTOR > >
KrylovSolver< VECTOR >::Create(
  LinearSolverParameters const & parameters,
  LinearOperator< VECTOR > const & matrix,
  LinearOperator< VECTOR > const & precond )
{
  if( parameters.solverType == "cg" )
  {
    GEOSX_ERROR_IF( !parameters.isSymmetric,
                    "Cannot use CG solver with a non-symmetric system" );
    return std::make_unique< CGsolver< Vector > >( matrix,
                                                   precond,
                                                   parameters.krylov.relTolerance,
                                                   parameters.krylov.maxIterations,
                                                   parameters.logLevel );
  }
  else if( parameters.solverType == "bicgstab" )
  {
    return std::make_unique< BiCGSTABsolver< Vector > >(
      matrix,
      precond,
      parameters.krylov.relTolerance,
      parameters.krylov.maxIterations,
      parameters.logLevel );
  }
  else if( parameters.solverType == "gmres" )
  {
    return std::make_unique< GMRESsolver< Vector > >( matrix,
                                                      precond,
                                                      parameters.krylov.relTolerance,
                                                      parameters.krylov.maxIterations,
                                                      parameters.logLevel,
                                                      parameters.krylov.maxRestart );
  }
  else
  {
    GEOSX_ERROR( "Unsupported linear solver type: " << parameters.solverType );
    return {};
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class KrylovSolver< TrilinosInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< TrilinosInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_HYPRE
template class KrylovSolver< HypreInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< HypreInterface::ParallelVector > >;
#endif

#ifdef GEOSX_USE_PETSC
template class KrylovSolver< PetscInterface::ParallelVector >;
template class KrylovSolver< BlockVectorView< PetscInterface::ParallelVector > >;
#endif

}  // namespace geosx
