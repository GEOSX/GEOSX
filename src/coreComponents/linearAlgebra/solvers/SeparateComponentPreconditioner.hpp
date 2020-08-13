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
 * @file SeparateComponentPreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_

#include "PreconditionerBase.hpp"

#include <memory>

namespace geosx
{
/**
 * @brief Separate component filter implemented as a compound preconditioner.
 * @tparam LAI linear algebra interface to use
 */
template< typename LAI >
class SeparateComponentPreconditioner : public PreconditionerBase< LAI >
{
public:
  /// Alias for base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param numComp number of components in the field
   * @param precond the actual preconditioner to apply to filtered matrix (ownership transferred)
   */
  SeparateComponentPreconditioner( localIndex const numComp,
                                   std::unique_ptr< PreconditionerBase< LAI > > precond );

  /**
   * @brief Destructor.
   */
  virtual ~SeparateComponentPreconditioner() override;

  using PreconditionerBase< LAI >::compute;

  virtual void
  compute( Matrix const & mat, DofManager const & dofManager ) override;

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
   * @brief Access the preconditioning matrix.
   * @return reference to the filtered matrix
   */
  Matrix const &
  getPrecondMatrix() const
  {
    return m_matSC;
  }

  /**
   * @brief Access to the nested preconditioner.
   * @return reference to the preconditioner passed at construction
   */
  PreconditionerBase< LAI > const &
  getNestedPrecond() const
  {
    return *m_precond;
  }

private:
  /// Number of components in the matrix
  localIndex m_numComp;

  /// Separate component approximation of the original matrix
  Matrix m_matSC;

  /// Actual preconditioner
  std::unique_ptr< PreconditionerBase< LAI > > m_precond;
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_SOLVERS_SEPARATECOMPONENTPRECONDITIONER_HPP_
