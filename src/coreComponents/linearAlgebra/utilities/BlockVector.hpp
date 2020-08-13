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
 * @file BlockVector.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_

#include "linearAlgebra/utilities/BlockVectorView.hpp"

namespace geosx
{
/**
 * @brief Concrete representation of a block vector.
 * @tparam VECTOR type of sub-vectors
 *
 * This extends BlockVectorView class by providing storage for sub-block vectors.
 * The @p VECTOR type needs to be default-constructible.
 */
template< typename VECTOR >
class BlockVector : public BlockVectorView< VECTOR >
{
public:
  /// Alias for base type
  using Base = BlockVectorView< VECTOR >;

  /**
   * @brief Create a vector of @p nBlocks blocks.
   * @param nBlocks Number of blocks
   */
  explicit BlockVector( localIndex const nBlocks ) :
    Base( nBlocks ),
    m_vectorStorage( nBlocks )
  {
    setPointers();
  }

  /**
   * @brief Create a vector of @p nBlocks blocks.
   */
  explicit BlockVector() :
    BlockVector( 0 )
  {}

  /**
   * @brief Copy constructor that performs a deep copy of each sub-vector.
   * @param rhs the block vector to copy
   */
  BlockVector( BlockVector const & rhs ) :
    Base( rhs ),
    m_vectorStorage( rhs.m_vectorStorage )
  {
    setPointers();
  }

  /**
   * @brief Move constructor.
   * @param rhs the block vector to move from
   */
  BlockVector( BlockVector && rhs ) :
    Base( std::move( rhs ) ),
    m_vectorStorage( std::move( rhs.m_vectorStorage ) )
  {
    setPointers();
  }

  /**
   * @brief Conversion constructor from a compatible view with a deep copy of each sub-vector.
   * @param rhs the block vector view to copy from
   * @note declared explicit to avoid unintended deep copying
   */
  explicit BlockVector( BlockVectorView< VECTOR > const & rhs ) :
    Base( rhs.blockSize() )
  {
    for( localIndex i = 0; i < rhs.blockSize(); ++i )
    {
      m_vectorStorage.emplace_back( rhs.block( i ) );
    }
    setPointers();
  }

  /**
   * @brief Copy assignment
   * @param x the vector to copy
   * @return reference to @p this
   */
  BlockVector &
  operator=( BlockVector const & x )
  {
    m_vectorStorage = x.m_vectorStorage;
    setPointers();
    return *this;
  }

  /**
   * @brief Move assignment
   * @param x the vector to move from
   * @return reference to @p this
   */
  BlockVector &
  operator=( BlockVector && x ) noexcept
  {
    m_vectorStorage = std::move( x.m_vectorStorage );
    setPointers();
    return *this;
  }

  /**
   * @brief Destructor.
   */
  virtual ~BlockVector() override = default;

  /**
   * @brief Resize to a different number of blocks.
   * @param nBlocks the new number of blocks
   *
   * @note If the new number of blocks is larger than the previous, new vectors
   *       will not be initialized. It is the user's responsibility to do that.
   */
  void
  resize( localIndex const nBlocks )
  {
    m_vectorStorage.resize( nBlocks );
    setPointers();
  }

private:
  void
  setPointers()
  {
    Base::resize( m_vectorStorage.size() );
    for( localIndex i = 0; i < m_vectorStorage.size(); ++i )
    {
      this->setPointer( i, &m_vectorStorage[i] );
    }
  }

  /// Storage for actual vectors
  array1d< VECTOR > m_vectorStorage;
};

}  //namespace geosx

#endif  //GEOSX_LINEARALGEBRA_UTILITIES_BLOCKVECTOR_HPP_
