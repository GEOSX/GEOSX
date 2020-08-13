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
 * @file HypreSolver.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_HYPRESOLVER_HPP_
#define GEOSX_LINEARALGEBRA_HYPRESOLVER_HPP_

#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LinearSolverResult.hpp"

namespace geosx
{
class DofManager;
class HypreVector;
class HypreMatrix;

/**
 * @brief This class creates and provides basic support for Hypre solvers.
 */
class HypreSolver
{
public:
  /**
   * @brief Solver constructor, with parameter list reference
   * @param[in] parameters structure containing linear solver parameters
   */
  HypreSolver( LinearSolverParameters parameters );

  /**
   * @brief Virtual destructor.
   */
  virtual ~HypreSolver() = default;

  /**
   * @brief Solve system with an iterative solver (HARD CODED PARAMETERS, GMRES).
   * @param[in,out] mat the matrix
   * @param[in,out] sol the solution
   * @param[in,out] rhs the right-hand side
   * @param dofManager the Degree-of-Freedom manager associated with matrix
   *
   * Solve Ax=b with A an HypreMatrix, x and b HypreVector.
   */
  void
  solve( HypreMatrix & mat,
         HypreVector & sol,
         HypreVector & rhs,
         DofManager const * const dofManager = nullptr );

  /**
   * @brief Get the result of previous solve.
   * @return struct with last solve stats
   */
  LinearSolverResult const &
  result()
  {
    return m_result;
  }

private:
  LinearSolverParameters m_parameters;
  LinearSolverResult m_result;

  void
  solve_direct( HypreMatrix & mat, HypreVector & sol, HypreVector & rhs );

  void
  solve_krylov( HypreMatrix & mat,
                HypreVector & sol,
                HypreVector & rhs,
                DofManager const * const dofManager );
};

}  // namespace geosx

#endif /* HYPRESOLVER_HPP_ */
