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
 * @file PetscPreconditioner.hpp
 */

#ifndef GEOSX_PETSCPRECONDITIONER_HPP
#define GEOSX_PETSCPRECONDITIONER_HPP

#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/petsc/PetscInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

/**
 * @name Hypre forward declarations.
 *
 * Forward declare PETSc's solver structs and pointer aliases in order
 * to avoid including PETSc headers and leaking into the rest of GEOSX.
 */
///@{

/// PETSc preconditioner struct forward declaration
extern "C" struct _p_PC;

/// Preconditioner pointer alias
using PC = _p_PC *;

///@}

namespace geosx
{
/**
 * @brief Wrapper around PETSc-based preconditioners.
 */
class PetscPreconditioner final : public PreconditionerBase< PetscInterface >
{
public:
  /// Alias for base type
  using Base = PreconditionerBase< PetscInterface >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /// Allow for partial overload of Base::compute()
  using Base::compute;

  /**
   * @brief Constructor.
   * @param params preconditioner parameters
   */
  explicit PetscPreconditioner( LinearSolverParameters params );

  /**
   * @brief Destructor.
   */
  virtual ~PetscPreconditioner() override;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void
  compute( Matrix const & mat ) override;

  /**
   * @brief Apply operator to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   *
   * @warning @p src and @p dst cannot alias the same vector.
   */
  virtual void
  apply( Vector const & src, Vector & dst ) const override;

  virtual void
  clear() override;

  /**
   * @brief Access the underlying implementation.
   * @return the wrapped PETSc preconditioner
   */
  PC const &
  unwrapped() const;

private:
  /// Parameters for all preconditioners
  LinearSolverParameters m_parameters;

  /// Pointer to the PETSc implementation
  PC m_precond;
};

}  // namespace geosx

#endif  //GEOSX_PETSCPRECONDITIONER_HPP
