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
 * @file HypreInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_

#include "linearAlgebra/interfaces/hypre/HypreSolver.hpp"
#include "linearAlgebra/interfaces/hypre/HypreMatrix.hpp"
#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"

#include <memory>

namespace geosx
{
/**
 * @class HypreInterface
 * @brief This class holds aliases based on the Hypre library.
 */
struct HypreInterface
{
  /**
   * @brief Initializes the MPI environment for the Hypre library
   *
   * @param[in] argc standard argc as in any C main
   * @param[in] argv standard argv as in any C main
   */
  static void
  initialize( int & argc, char **& argv );

  /**
   * @brief Finalizes the MPI environment for the Hypre library
   */
  static void
  finalize();

  /**
   * @brief Create a hypre-based preconditioner object.
   * @param params the preconditioner parameters
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< HypreInterface > >
  createPreconditioner(
    LinearSolverParameters params );

  /// Alias for HypreMatrix
  using ParallelMatrix = HypreMatrix;
  /// Alias for HypreVector
  using ParallelVector = HypreVector;
  /// Alias for HypreSolver
  using LinearSolver = HypreSolver;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREINTERFACE_HPP_*/
