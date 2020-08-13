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
 * @file VectorBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_

#include "linearAlgebra/common.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
/**
 * @brief Common base template for all vector wrapper types.
 * @tparam VECTOR derived vector type
 *
 * This class template provides a common interface for all derived vector
 * wrapper types. Most methods are pure abstract in order to get the compiler
 * to enforce a common interface in all derived classes; however there is no
 * runtime polymorphism or virtual dispatch - derived classes are not related
 * between each other, and pointers/references to base should never be used -
 * the destructor is protected to enforce that.
 *
 * As an added benefit, the documentation for vector interface also lives here
 * and does not need to be duplicated across vector wrappers. Derived classes
 * should still document specific functions in case a particular LA package's
 * behavior deviates from expectations or has unexpected performance impacts.
 * In that case, @c \@copydoc tag can be used to copy over the documentation.
 */
template< typename VECTOR >
class VectorBase
{
protected:
  /// Alias for VECTOR
  using Vector = VECTOR;

  /**
   * @name Constructors/destructor/assignment
   */
  ///@{

  /**
   * @brief Constructs a vector in default state
   */
  VectorBase() :
    m_closed( true )
  {}

  /**
   * @brief Copy constructor.
   */
  VectorBase( VectorBase const & ) = default;

  /**
   * @brief Move constructor.
   */
  VectorBase( VectorBase && ) = default;

  /**
   * @brief Copy assignment.
   * @return reference to this object
   */
  VectorBase &
  operator=( VectorBase const & ) = default;

  /**
   * @brief Move assignment.
   * @return reference to this object
   */
  VectorBase &
  operator=( VectorBase && ) = default;

  /**
   * @brief Destructor.
   */
  ~VectorBase() = default;

  ///@}

  /**
   * @name Status query methods
   */
  ///@{

  /**
   * @brief Query vector closed status
   * @return @p true if vector has been opened and has not been closed since; @p false otherwise
   */
  inline bool
  closed() const
  {
    return m_closed;
  }

  /**
   * @brief Query vector creation status
   * @return @p true if vector has been created
   */
  virtual bool
  created() const = 0;

  /**
   * @brief Query vector ready status
   * @return @p true if vector has been created and is currently closed
   */
  inline bool
  ready() const
  {
    return created() && closed();
  }

  ///@}

  /**
   * @name Create Methods
   */
  ///@{

  /**
   * @brief Create a vector based on local number of elements.
   * @param localSize local number of elements
   * @param comm MPI communicator to use
   *
   * Create a vector based on local number of elements.  Global size is
   * the sum across processors.  For specifying a global size and having
   * automatic partitioning, see createWithGlobalSize().
   */
  virtual void
  createWithLocalSize( localIndex const localSize,
                       MPI_Comm const & comm ) = 0;

  /**
   * @brief Create a vector based on global number of elements.
   * @param globalSize Global number of elements
   * @param comm MPI communicator to use
   *
   * Create a vector based on global number of elements. Every processors
   * gets the same number of local elements except proc 0, which gets any
   * remainder elements as well if the split can't be done evenly.
   */
  virtual void
  createWithGlobalSize( globalIndex const globalSize,
                        MPI_Comm const & comm ) = 0;

  /**
   * @brief Construct parallel vector from a local array.
   * @param localValues local data to put into vector
   * @param comm MPI communicator to use
   */
  virtual void
  create( arrayView1d< real64 const > const & localValues,
          MPI_Comm const & comm ) = 0;

  ///@}

  /**
   * @name Open/close methods
   */
  ///@{

  /**
   * @brief Open the vector for modifying entries
   */
  virtual void
  open() = 0;

  /**
   * @brief Assemble vector
   *
   * Performs parallel communication to scatter assembled entries to appropriate locations
   */
  virtual void
  close() = 0;

  /**
   * @brief Reset the matrix to default state
   */
  virtual void
  reset()
  {
    m_closed = true;
  };

  ///@}

  /**
   * @name Add/Set Methods
   */
  ///@{

  /**
   * @brief Set vector value.
   * @param globalRow global row index
   * @param value Value to add at given row
   *
   * Set vector value at given element.
   */
  virtual void
  set( globalIndex const globalRow, real64 const value ) = 0;

  /**
   * @brief Add into vector value.
   * @param globalRow global row
   * @param value Values to add in given row
   *
   * Add into vector value at given row.
   */
  virtual void
  add( globalIndex const globalRow, real64 const value ) = 0;

  /**
   * @brief Set vector values.
   * @param globalIndices global row indices
   * @param values Values to add in given rows
   * @param size Number of elements
   *
   * Set vector values at given elements.
   */
  virtual void
  set( globalIndex const * globalIndices,
       real64 const * values,
       localIndex const size ) = 0;

  /**
   * @brief Add vector values.
   * @param globalIndices global row indices
   * @param values values to add in given rows
   * @param size number of elements
   *
   * Add vector values at given elements.
   */
  virtual void
  add( globalIndex const * globalIndices,
       real64 const * values,
       localIndex const size ) = 0;

  /**
   * @brief Set vector values using array1d
   * @param globalIndices global row indices
   * @param values values to add in given rows
   *
   * Set vector values at given elements.
   */
  virtual void
  set( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 const > const & values ) = 0;

  /**
   * @brief Add into vector values using array1d
   * @param globalIndices global rows indices
   * @param values values to add in given rows
   *
   * Add into vector values at given rows.
   */
  virtual void
  add( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 const > const & values ) = 0;

  /**
   * @brief Set all elements to a constant value.
   * @param value value to set vector elements to
   */
  virtual void
  set( real64 const value ) = 0;

  /**
   * @brief Set vector elements to zero.
   */
  virtual void
  zero() = 0;

  /**
   * @brief Set vector elements to random entries.
   * @param seed the random number seed to use
   */
  virtual void
  rand( unsigned const seed = 1984 ) = 0;

  ///@}

  /**
   * @name Algebraic Operations
   */
  ///@{

  /**
   * @brief Multiply all elements by scalingFactor.
   * @param scalingFactor scaling Factor
   */
  virtual void
  scale( real64 const scalingFactor ) = 0;

  /**
   * @brief Replace vector elements by their reciprocals
   * @note No guarding is done against division by zero.
   */
  virtual void
  reciprocal() = 0;

  /**
   * @brief Dot product with the vector vec.
   * @param vec vector to dot-product with
   * @return dot product
   */
  virtual real64
  dot( Vector const & vec ) const = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>x</tt>.
   * @param x vector to copy
   */
  virtual void
  copy( Vector const & x ) = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + y</tt>.
   * @param alpha scaling factor for added vector
   * @param x vector to add
   */
  virtual void
  axpy( real64 const alpha, Vector const & x ) = 0;

  /**
   * @brief Update vector <tt>y</tt> as <tt>y</tt> = <tt>alpha*x + beta*y</tt>.
   * @param alpha scaling factor for added vector
   * @param x vector to add
   * @param beta scaling factor for self vector
   */
  virtual void
  axpby( real64 const alpha, Vector const & x, real64 const beta ) = 0;

  /**
   * @brief 1-norm of the vector.
   * @return the 1-norm value
   */
  virtual real64
  norm1() const = 0;

  /**
   * @brief 2-norm of the vector.
   * @return the 2-norm value
   */
  virtual real64
  norm2() const = 0;

  /**
   * @brief Infinity-norm of the vector.
   * @return the inf-norm value
   */
  virtual real64
  normInf() const = 0;

  ///@}

  /**
   * @name Accessor Methods
   */
  ///@{

  /**
   * @brief Returns the global of the vector.
   * @return the global size
   */
  virtual globalIndex
  globalSize() const = 0;

  /**
   * @brief Returns the local size of the vector.
   * @return the local size (on this processor)
   */
  virtual localIndex
  localSize() const = 0;

  /**
   * @brief Get lower bound of local partition
   * @return index of the first global row owned by this processor
   */
  virtual globalIndex
  ilower() const = 0;

  /**
   * @brief Get upper bound of local partition
   * @brief @return next index after last global row owned by that processor
   * @note The intention is for [ilower; iupper) to be used as a half-open index range
   */
  virtual globalIndex
  iupper() const = 0;

  /**
   * @brief Get a value by index.
   * @param globalRow global row index
   * @return value at global index @p globalRow
   */
  virtual real64
  get( globalIndex globalRow ) const = 0;

  /**
   * @brief Get a sequence of values by index.
   * @param[in] globalIndices array of global row indices
   * @param[out] values array of vector values
   */
  virtual void
  get( arraySlice1d< globalIndex const > const & globalIndices,
       arraySlice1d< real64 > const & values ) const = 0;

  /**
   * @brief Map a global row index to local row index.
   * @param[in] globalRow the global row index
   * @return global row index corresponding to @p globalRow
   */
  virtual localIndex
  getLocalRowID( globalIndex const globalRow ) const = 0;

  /**
   * @brief Map a local row index to global row index.
   * @param[in] localRow the local row index
   * @return global row index corresponding to @p localRow
   */
  virtual globalIndex
  getGlobalRowID( localIndex const localRow ) const = 0;

  /**
   * @brief Extract a view of the local portion of the array.
   * @return pointer to local vector data
   */
  virtual real64 const *
  extractLocalVector() const = 0;

  /**
   * @brief Extract a view of the local portion of the array.
   * @return pointer to local vector data
   */
  virtual real64 *
  extractLocalVector() = 0;

  /**
   * @brief Extract local solution by copying into a user-provided array
   * @param localVector the array view to write to (must be properly sized)
   */
  virtual void
  extract( arrayView1d< real64 > const & localVector ) const
  {
    real64 const * const data = extractLocalVector();
    forAll< parallelHostPolicy >( localSize(), [=]( localIndex const k ) {
      localVector[k] = data[k];
    } );
  }

  /**
   * @brief Get the communicator used by this vector
   * @return the MPI communicator
   */
  virtual MPI_Comm
  getComm() const = 0;

  ///@}

  /**
   * @name I/O Methods
   */
  ///@{

  /**
   * @brief Print the vector in Trilinos format to the terminal.
   * @param os the output stream to print to
   */
  virtual void
  print( std::ostream & os = std::cout ) const = 0;

  /**
   * @brief Write the vector to a file.
   * @param filename name of the output file
   * @param[in] format output format
   */
  virtual void
  write(
    string const & filename,
    LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const = 0;

  ///@}

  /**
   * @brief Stream insertion operator for all vector types
   * @param os the output stream to print to
   * @param vec the vector to print
   * @return reference to @p os
   */
  friend std::ostream &
  operator<<( std::ostream & os, Vector const & vec )
  {
    vec.print( os );
    return os;
  }

  /// Flag indicating whether the vector is closed
  bool m_closed;
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_VECTORBASE_HPP_
