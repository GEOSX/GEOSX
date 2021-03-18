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
 * @file TransposeOperator.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_TRANSPOSEMATRIXOPERATOR_HPP_
#define GEOSX_LINEARALGEBRA_TRANSPOSEMATRIXOPERATOR_HPP_

#include "linearAlgebra/common/LinearOperator.hpp"

namespace geosx
{

/**
 * @brief Simple class that wraps a matrix and represents its transpose as a linear operator.
 * @tparam LAI the linear algebra interface
 */
template< typename LAI >
class TransposeOperator : public LinearOperator< typename LAI::ParallelVector >
{
public:

  /// Alias for base type
  using Base = LinearOperator< typename LAI::ParallelVector >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /**
   * @brief Constructor.
   * @param mat the underlying matrix
   */
  explicit TransposeOperator( Matrix const & mat )
    : Base(),
    m_matrix( mat )
  { }

  /**
   * @brief Destructor.
   */
  virtual ~TransposeOperator() override = default;

  /**
   * @brief Apply operator to a vector.
   * @param src input vector (x)
   * @param dst output vector (b)
   *
   * @warning @p src and @p dst cannot alias the same vector (some implementations may allow this).
   */
  virtual void apply( Vector const & src, Vector & dst ) const override
  {
    m_matrix.applyTranspose( src, dst );
  }

  /**
   * @brief Get the number of global rows.
   * @return Number of global rows in the operator.
   */
  virtual globalIndex numGlobalRows() const override
  {
    return m_matrix.numGlobalCols();
  }

  /**
   * @brief Get the number of global columns.
   * @return Number of global columns in the operator.
   */
  virtual globalIndex numGlobalCols() const override
  {
    return m_matrix.numGlobalRows();
  }

  /**
   * @brief Get the number of local rows.
   * @return Number of local rows in the operator.
   */
  virtual localIndex numLocalRows() const override
  {
    return m_matrix.numLocalCols();
  }

  /**
   * @brief Get the number of local columns.
   * @return Number of local columns in the operator.
   */
  virtual localIndex numLocalCols() const override
  {
    return m_matrix.numLocalRows();
  }

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator of the underlying matrix
   */
  virtual MPI_Comm getComm() const override
  {
    return m_matrix.getComm();
  }

private:

  Matrix const & m_matrix;
};

}

#endif //GEOSX_LINEARALGEBRA_TRANSPOSEMATRIXOPERATOR_HPP_
