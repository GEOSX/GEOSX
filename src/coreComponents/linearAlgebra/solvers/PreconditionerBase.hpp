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

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_

#include "linearAlgebra/common.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"

namespace geosx
{
class DofManager;

/**
 * @brief Common interface for preconditioning operators
 * @tparam LAI linear algebra interface providing vectors, matrices and solvers
 */
template< typename LAI >
class PreconditionerBase : public LinearOperator< typename LAI::ParallelVector >
{
public:
  PreconditionerBase() :
    m_mat {}
  {}

  virtual ~PreconditionerBase() = default;

  /// Alias for base type
  using Base = LinearOperator< typename LAI::ParallelVector >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Compute the preconditioner from a matrix.
   * @param mat the matrix to precondition.
   */
  virtual void
  compute( Matrix const & mat )
  {
    GEOSX_LAI_ASSERT( mat.ready() );
    m_mat = &mat;
  }

  /**
   * @brief Compute the preconditioner from a matrix
   * @param mat the matrix to precondition
   * @param dofManager the Degree-of-Freedom manager associated with matrix
   */
  virtual void
  compute( Matrix const & mat, DofManager const & dofManager )
  {
    GEOSX_UNUSED_VAR( dofManager );
    compute( mat );
  }

  /**
   * @brief Clean up the preconditioner setup.
   *
   * Releases memory used and allows the matrix to be deleted cleanly.
   * This method should be called before the matrix used to compute the preconditioner
   * goes out of scope or is re-created. Some implementations require the matrix
   * to outlive the preconditioner (for example, Trilinos/ML may crash the program if
   * deleted after the matrix).
   *
   * @note Should be properly overridden in derived classes, which may call this method.
   */
  virtual void
  clear()
  {
    m_mat = nullptr;
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex
  numGlobalRows() const override
  {
    return m_mat->numGlobalRows();
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex
  numGlobalCols() const override
  {
    return m_mat->numGlobalCols();
  }

  /**
   * @brief Chech if preconditioner is ready to use
   * @return @p true if compute() has been called but not clear().
   */
  bool
  ready() const
  {
    return m_mat != nullptr;
  }

  /**
   * @brief Access the matrix the preconditioner was computed from
   * @return reference to the matrix (user's repsonsibility to ensure it's still valid)
   */
  Matrix const &
  matrix() const
  {
    GEOSX_LAI_ASSERT( ready() );
    return *m_mat;
  }

private:
  /// Pointer to the matrix
  Matrix const * m_mat;
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_SOLVERS_PRECONDITIONERBASE_HPP_
