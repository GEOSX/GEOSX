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
 * @file TrilinosPreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_

#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

class Epetra_Operator;

namespace Teuchos
{
class ParameterList;
}

namespace geosx
{
/**
 * @brief Wrapper around Trilinos-based preconditioners.
 */
class TrilinosPreconditioner final : public PreconditionerBase< TrilinosInterface >
{
public:
  /// Alias for base type
  using Base = PreconditionerBase< TrilinosInterface >;

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
  explicit TrilinosPreconditioner( LinearSolverParameters params );

  /**
   * @brief Destructor.
   */
  virtual ~TrilinosPreconditioner() override;

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
   * @brief Access the underlying Epetra preconditioning operator.
   * @return reference to the Epetra operator
   */
  Epetra_Operator const &
  unwrapped() const;

  /**
   * @copydoc unwrapped() const
   */
  Epetra_Operator &
  unwrapped();

private:
  /// Parameters for all preconditioners
  LinearSolverParameters m_parameters;

  /// Pointer to the Trilinos implementation
  std::unique_ptr< Epetra_Operator > m_precond;
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_TRILINOSAMG_HPP_
