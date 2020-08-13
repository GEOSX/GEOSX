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
 * @file BlockOperatorWrapper.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_

#include "linearAlgebra/utilities/BlockOperatorView.hpp"
#include "linearAlgebra/common.hpp"

namespace geosx
{
/**
 * @brief "Shallow" representation of a block operator.
 * @tparam VECTOR type of vector that sub-blocks of this view can operate on
 * @tparam OPERATOR type of operator that can operate on @p VECTOR
 *                  (can be base class or a more specialized derived class)
 *
 * This extends BlockOperatorView class by providing a way to assign sub-block pointers.
 * The sub-blocks themselves must be stored elsewhere.
 * Therefore, it's an easy way to assemble a block operator representation from pre-existing blocks.
 */
template< typename VECTOR, typename OPERATOR = LinearOperator< VECTOR > >
class BlockOperatorWrapper : public BlockOperatorView< VECTOR, OPERATOR >
{
public:
  /// Alias for base type
  using Base = BlockOperatorView< VECTOR, OPERATOR >;

  /// Alias for vector type
  using Vector = typename Base::Vector;

  /**
   * @brief Create a vector wrapper of @p nBlocks blocks.
   * @param nRows number of block rows
   * @param nCols number of block columns
   */
  explicit BlockOperatorWrapper( localIndex const nRows, localIndex const nCols ) :
    Base( nRows, nCols )
  {}

  /**
   * @brief Deleted copy constructor.
   * @param rhs the block operator to copy
   */
  BlockOperatorWrapper( BlockOperatorWrapper const & rhs ) = delete;

  /**
   * @brief Deleted move constructor.
   * @param rhs the block operator to move from
   */
  BlockOperatorWrapper( BlockOperatorWrapper && rhs ) = delete;

  /**
   * @brief Destructor.
   */
  virtual ~BlockOperatorWrapper() override = default;

  /**
   * @brief Set a single block of the operator.
   * @param blockRowIndex block row index
   * @param blockColIndex block column index
   * @param op            reference to the operator (which must not go out of scope before the block wrapper)
   */
  void
  set( localIndex const blockRowIndex,
       localIndex const blockColIndex,
       OPERATOR & op )
  {
    this->setPointer( blockRowIndex, blockColIndex, &op );
  }
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKOPERATORWRAPPER_HPP_
