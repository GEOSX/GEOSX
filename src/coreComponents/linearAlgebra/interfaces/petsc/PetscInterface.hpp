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
 * @file PetscInterface.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_

#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/interfaces/petsc/PetscMatrix.hpp"
#include "linearAlgebra/interfaces/petsc/PetscSolver.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"

#include <memory>

namespace geosx
{
/**
 * @class PetscInterface
 * @brief This class holds aliases based on the Petsc library.
 */
struct PetscInterface
{
  /**
   * @brief Initializes the MPI environment for the Petsc library
   *
   * @param[in] argc standard argc as in any C main
   * @param[in] argv standard argv as in any C main
   *
   * Essentially, it is a wrapper for PetscInitialize
   */
  static void
  initialize( int & argc, char **& argv );

  /**
   * @brief Finalizes the MPI environment for the Petsc library
   *
   * Essentially, it is a wrapper for PetscFinalize
   */
  static void
  finalize();

  /**
   * @brief Create a PETSc-based preconditioner object.
   * @param params the parameters for preconditioner
   * @return owning pointer to the newly created preconditioner
   */
  static std::unique_ptr< PreconditionerBase< PetscInterface > >
  createPreconditioner(
    LinearSolverParameters params );

  /// Alias for PetscMatrix
  using ParallelMatrix = PetscMatrix;
  /// Alias for PetscVector
  using ParallelVector = PetscVector;
  /// Alias for PetscSolver
  using LinearSolver = PetscSolver;
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCINTERFACE_HPP_*/
