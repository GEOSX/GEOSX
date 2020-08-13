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
 * @file MatrixBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP_

#include "linearAlgebra/common.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
//#include "LvArray/src/streamIO.hpp"

namespace geosx
{
/**
 * @brief Common base template for all matrix wrapper types.
 * @tparam MATRIX derived matrix type
 * @tparam VECTOR compatible vector type
 *
 * This class template provides a common interface for all derived matrix
 * wrapper types. Most methods are pure abstract in order to get the compiler
 * to enforce a common interface in all derived classes; the class should be
 * inherited from privately by implementations to avoid the possibility of
 * accidentally invoking virtual function calls. Some basic facilities are
 * provided related to matrix lifecycle. These include status flags and
 * corresponding status query functions.
 *
 * As an added benefit, the documentation for matrix interface also lives here
 * and does not need to be duplicated across matrix wrappers. Derived classes
 * should still document specific functions in case a particular LA package's
 * behavior deviates from expectations or has unexpected performance impacts.
 * In that case, @c \@copydoc tag can be used to copy over the documentation.
 */
template< typename MATRIX, typename VECTOR >
class MatrixBase : public virtual LinearOperator< VECTOR >
{
public:
  /// Type alias for actual derived matrix class
  using Matrix = MATRIX;

  /// Type alias for a compatible vector class
  using Vector = VECTOR;

protected:
  /**
   * @name Constructors/destructor/assignment
   */
  ///@{

  /**
   * @brief Constructs a matrix in default state
   */
  MatrixBase() :
    m_closed( true ),
    m_assembled( false )
  {}

  /**
   * @brief Copy constructor.
   */
  MatrixBase( MatrixBase const & ) = default;

  /**
   * @brief Move constructor.
   */
  MatrixBase( MatrixBase && ) = default;

  /**
   * @brief Copy assignment.
   * @return reference to this object
   */
  MatrixBase &
  operator=( MatrixBase const & ) = default;

  /**
   * @brief Move assignment.
   * @return reference to this object
   */
  MatrixBase &
  operator=( MatrixBase && ) = default;

  /**
   * @brief Destructor.
   */
  ~MatrixBase() = default;

  ///@}

  /**
   * @name Status query methods
   */
  ///@{

  /**
   * @brief Query matrix closed status
   * @return @p true if matrix has been opened and has not been closed since; @p false otherwise
   */
  inline bool
  closed() const
  {
    return m_closed;
  }

  /**
   * @brief Query matrix assembled status
   * @return @p true if matrix has been opened and closed since creation; @p false otherwise
   */
  inline bool
  assembled() const
  {
    return m_assembled;
  }

  /**
   * @brief Query matrix ready status
   * @return @p true if matrix has been assembled and is currently closed;
   *         this implies it's ready to be used or re-opened for adding/setting values
   */
  inline bool
  ready() const
  {
    return closed() && assembled();
  }

  /**
   * @brief Query matrix status
   * @return @p true if matrix has been assembled and is currently open;
   *         this implies individual entries within existing sparsity pattern
   *         can be altered via set()/add() methods.
   */
  inline bool
  modifiable() const
  {
    return !closed() && assembled();
  }

  /**
   * @brief Query matrix status
   * @return @p true if matrix has NOT been assembled yet (not closed since
   *         last create() call) and is currently open for insertion of new entries
   */
  inline bool
  insertable() const
  {
    return !closed() && !assembled();
  }

  /**
   * @brief Query matrix creation status
   * @return @p true if matrix has been created
   */
  virtual bool
  created() const = 0;

  ///@}

  /**
   * @name Create Methods
   */
  ///@{

  /**
   * @brief Create a square matrix from local number of rows.
   *
   * @param localSize local number of rows for square matrix.
   * @param maxEntriesPerRow Maximum number of non-zero entries per row.
   * @param comm MPI communicator.
   *
   */
  virtual void
  createWithLocalSize( localIndex const localSize,
                       localIndex const maxEntriesPerRow,
                       MPI_Comm const & comm )
  {
    createWithLocalSize( localSize, localSize, maxEntriesPerRow, comm );
  }

  /**
   * @brief Create a square matrix from global number of rows.
   *
   * Create a square matrix with an (approximately) even partitioning of rows.
   *
   * @param globalSize Global dimensions for a square matrix.
   * @param maxEntriesPerRow Maximum number of non-zero entries per row.
   * @param comm MPI communicator.
   *
   */
  virtual void
  createWithGlobalSize( globalIndex const globalSize,
                        localIndex const maxEntriesPerRow,
                        MPI_Comm const & comm )
  {
    createWithGlobalSize( globalSize, globalSize, maxEntriesPerRow, comm );
  }

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * @param comm MPI communicator.
   * @param localRows Local number of rows.
   * @param localCols Local number of columns.
   * @param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  virtual void
  createWithLocalSize( localIndex const localRows,
                       localIndex const localCols,
                       localIndex const maxEntriesPerRow,
                       MPI_Comm const & comm ) = 0;

  /**
   * @brief Create a rectangular matrix from number of rows/columns.
   *
   * @param comm MPI communicator.
   * @param globalRows Global number of rows.
   * @param globalCols Global number of columns.
   * @param maxEntriesPerRow Maximum number of entries per row (hint).
   */
  virtual void
  createWithGlobalSize( globalIndex const globalRows,
                        globalIndex const globalCols,
                        localIndex const maxEntriesPerRow,
                        MPI_Comm const & comm ) = 0;

  /**
   * @brief Create parallel matrix from a local CRS matrix.
   * @param localMatrix The input local matrix.
   * @param comm The MPI communicator to use.
   *
   * @note Copies values, so that @p localMatrix does not need to retain its values after the call.
   * @todo Replace generic implementation with more efficient ones in each package.
   */
  virtual void
  create( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
          MPI_Comm const & comm )
  {
    localMatrix.move( LvArray::MemorySpace::CPU, false );

    localIndex maxEntriesPerRow = 0;
    for( localIndex i = 0; i < localMatrix.numRows(); ++i )
    {
      maxEntriesPerRow = std::max( maxEntriesPerRow, localMatrix.numNonZeros( i ) );
    }

    createWithLocalSize( localMatrix.numRows(),
                         localMatrix.numRows(),
                         maxEntriesPerRow,
                         comm );

    globalIndex const rankOffset = ilower();

    open();
    for( localIndex localRow = 0; localRow < localMatrix.numRows(); ++localRow )
    {
      insert( localRow + rankOffset,
              localMatrix.getColumns( localRow ),
              localMatrix.getEntries( localRow ) );
    }
    close();
  }

  ///@}

  /**
   * @name Open/close methods
   */
  ///@{

  /**
   * @brief Open matrix for adding new entries.
   *
   * @note Adding entries that result in modifications of sparsity pattern may not be allowed
   *       by most implementations. An error will be raised in that case.
   */
  virtual void
  open() = 0;

  /**
   * @brief Assemble and compress the matrix.
   *
   * Compresses the matrix to CSR format with contiguous memory on each processor. Prevents from
   * adding new entries in the sparsity pattern but allows for modification of existing entries.
   */
  virtual void
  close() = 0;

  /**
   * @brief Reset the matrix to default state
   */
  virtual void
  reset()
  {
    m_assembled = false;
    m_closed = true;
  }

  ///@}

  /**
   * @name Global modification methods
   */
  ///@{

  /**
   * @brief Set all non-zero elements to a value.
   * @param value the value to set all elements to
   */
  virtual void
  set( real64 const value ) = 0;

  /**
   * @brief Set all elements to zero.
   */
  virtual void
  zero() = 0;

  ///@}

  /**
   * @name Add/Set/Insert Methods
   *
   * @p insert method can only be used while the matrix is being filled for the first time,
   * e.g. during sparsity pattern construction. @p add and @p set methods can only be used
   * after the matrix has been assembled, and both only accept entries that already exist
   * in the sparsity pattern.
   *
   * @note Caution: these methods are not thread-safe.
   */
  ///@{

  /**
   * @brief Add to one element.
   * @param rowIndex Global row index
   * @param colIndex Global column index
   * @param value Value to add to prescribed location
   */
  virtual void
  add( globalIndex const rowIndex,
       globalIndex const colIndex,
       real64 const value ) = 0;

  /**
   * @brief Set one element.
   * @param rowIndex Global row index
   * @param colIndex Global column index
   * @param value Value to set at prescribed location
   */
  virtual void
  set( globalIndex const rowIndex,
       globalIndex const colIndex,
       real64 const value ) = 0;

  /**
   * @brief Insert one element.
   * @param rowIndex Global row index
   * @param colIndex Global column index
   * @param value Value to insert at prescribed location
   */
  virtual void
  insert( globalIndex const rowIndex,
          globalIndex const colIndex,
          real64 const value ) = 0;

  /**
   * @brief Add elements to one row using c-style arrays
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   * @param size Number of elements
   */
  virtual void
  add( globalIndex const rowIndex,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const size ) = 0;

  /**
   * @brief Set elements to one row using c-style arrays
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   * @param size Number of elements
   */
  virtual void
  set( globalIndex const rowIndex,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const size ) = 0;

  /**
   * @brief Insert elements to one row using c-style arrays
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   * @param size Number of elements
   */
  virtual void
  insert( globalIndex const rowIndex,
          globalIndex const * colIndices,
          real64 const * values,
          localIndex const size ) = 0;

  /**
   * @brief Add elements to one row using array1d
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   */
  virtual void
  add( globalIndex const rowIndex,
       arraySlice1d< globalIndex const > const & colIndices,
       arraySlice1d< real64 const > const & values ) = 0;

  /**
   * @brief Set elements of one row using array1d
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   */
  virtual void
  set( globalIndex const rowIndex,
       arraySlice1d< globalIndex const > const & colIndices,
       arraySlice1d< real64 const > const & values ) = 0;

  /**
   * @brief Insert elements of one row using array1d
   * @param rowIndex Global row index
   * @param colIndices Global column indices
   * @param values Values to add to prescribed locations
   */
  virtual void
  insert( globalIndex const rowIndex,
          arraySlice1d< globalIndex const > const & colIndices,
          arraySlice1d< real64 const > const & values ) = 0;

  /**
   * @brief Add a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  add(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) = 0;

  /**
   * @brief Set a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  set(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) = 0;

  /**
   * @brief Insert a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  insert(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) = 0;

  /**
   * @brief Add a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  add(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) = 0;

  /**
   * @brief Set a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  set(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) = 0;

  /**
   * @brief Insert a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   */
  virtual void
  insert(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) = 0;

  /**
   * @brief Add a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   * @param numRows Number of row indices
   * @param numCols Number of column indices
   *
   * @note Row major layout assumed in values
   */
  virtual void
  add( globalIndex const * rowIndices,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const numRows,
       localIndex const numCols ) = 0;

  /**
   * @brief Set a dense block of values.
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   * @param numRows Number of row indices
   * @param numCols Number of column indices
   *
   * @note Row major layout assumed in values
   */
  virtual void
  set( globalIndex const * rowIndices,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const numRows,
       localIndex const numCols ) = 0;

  /**
   * @brief Insert dense matrix.
   *
   * @param rowIndices Global row indices
   * @param colIndices Global col indices
   * @param values Dense local matrix of values
   * @param numRows Number of row indices
   * @param numCols Number of column indices
   * @note Row major layout assumed in values
   */
  virtual void
  insert( globalIndex const * rowIndices,
          globalIndex const * colIndices,
          real64 const * values,
          localIndex const numRows,
          localIndex const numCols ) = 0;

  ///@}

  /**
   * @name Linear Algebra Methods
   */
  ///@{

  virtual void
  apply( Vector const & src, Vector & dst ) const override = 0;

  /**
   * @brief Apply transpose of the matrix to a vector
   * @param src Input vector (x).
   * @param dst Output vector (b).
   */
  virtual void
  applyTranspose( Vector const & src, Vector & dst ) const = 0;

  /**
   * @brief Compute residual <tt>r = Ax - b</tt>.
   *
   * Overrides LinearOperator::residual().
   *
   * @param x Input solution.
   * @param b Input right hand side.
   * @param r Output residual.
   *
   * @note @p x and @p r cannot alias the same vector.
   */
  virtual void
  residual( Vector const & x, Vector const & b, Vector & r ) const override
  {
    if( &b != &r )
    {
      r.copy( b );
    }
    gemv( -1.0, x, 1.0, r );
  }

  /**
   * @brief Matrix/Matrix multiplication.
   *
   * Compute <tt>this * B = C</tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void
  multiply( Matrix const & src, Matrix & dst ) const = 0;

  /**
   * @brief Matrix/Matrix transpose multiplication.
   *
   * Compute <tt>this^T * B = C</tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void
  leftMultiplyTranspose( Matrix const & src, Matrix & dst ) const = 0;

  /**
   * @brief Matrix/Matrix transpose multiplication.
   *
   * Compute <tt>B * this^T = C</tt>.
   *
   * @param src Input matrix (B).
   * @param dst Output matrix (C).
   *
   * Note that the output matrix C should have the same
   * row-map as this.  If close() has already been called
   * on C, then C's sparsity pattern must already contain
   * the nonzero entries produced by the product this*B.
   */
  virtual void
  rightMultiplyTranspose( Matrix const & src, Matrix & dst ) const = 0;

  /**
   * @brief Compute the triple product <tt>dst = R * this * P</tt>
   * @param R the "restriction" matrix
   * @param P the "prolongation" matrix
   * @param dst the resulting product matrix (will be re-created as needed)
   */
  virtual void
  multiplyRAP( Matrix const & R, Matrix const & P, Matrix & dst ) const
  {
    GEOSX_LAI_ASSERT( ready() );
    GEOSX_LAI_ASSERT( R.ready() );
    GEOSX_LAI_ASSERT( P.ready() );
    GEOSX_LAI_ASSERT_EQ( numGlobalRows(), R.numGlobalCols() );
    GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

    Matrix AP;
    multiply( P, AP );
    R.multiply( AP, dst );
  }

  /**
   * @brief Compute the triple product <tt>dst = P^T * this * P</tt>
   * @param P the "prolongation" matrix
   * @param dst the resulting product matrix (will be re-created as needed)
   */
  virtual void
  multiplyPtAP( Matrix const & P, Matrix & dst ) const
  {
    GEOSX_LAI_ASSERT( ready() );
    GEOSX_LAI_ASSERT( P.ready() );
    GEOSX_LAI_ASSERT_EQ( numGlobalRows(), P.numGlobalRows() );
    GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

    Matrix AP;
    multiply( P, AP );
    P.leftMultiplyTranspose( AP, dst );
  }

  /**
   * @brief Compute gemv <tt>y = alpha*A*x + beta*y</tt>.
   *
   * @note The naming convention follows the BLAS library.
   *
   * @param alpha Scalar factor for added matvec product.
   * @param x Input vector.
   * @param beta Scalar factor for right hand side.
   * @param y Output vector.
   * @param useTranspose Boolean, set to true to use <tt>A^T</tt>.
   *
   * @warning @p x and @p y cannot alias the same vector.
   */
  virtual void
  gemv( real64 const alpha,
        Vector const & x,
        real64 const beta,
        Vector & y,
        bool useTranspose = false ) const = 0;

  /**
   * @brief Multiply all elements by scalingFactor.
   * @param scalingFactor Scaling factor.
   */
  virtual void
  scale( real64 const scalingFactor ) = 0;

  /**
   * @brief Pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   * @param vec Vector to pre-multiply with.
   */
  virtual void
  leftScale( Vector const & vec ) = 0;

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vec.
   * @param vec Vector to post-multiply with.
   */
  virtual void
  rightScale( Vector const & vec ) = 0;

  /**
   * @brief Post-multiplies (right) with diagonal matrix consisting of the values in vecRight
   * and pre-multiplies (left) with diagonal matrix consisting of the values in vec.
   * @param vecLeft vec to pre-multiply with.
   * @param vecRight vec to post-multiply with.
   */
  virtual void
  leftRightScale( Vector const & vecLeft, Vector const & vecRight ) = 0;

  /**
   * @brief Matrix transposition.
   *
   * Compute <tt>B = this^T</tt>.
   *
   * @param dst Output matrix (B).
   *
   */
  virtual void
  transpose( Matrix & dst ) const = 0;

  /**
   * @brief Clear a row, and optionally set diagonal element to <tt>diagValue</tt>.
   * @param row globalIndex of the row to be cleared.
   * @param diagValue (Optional) set diagonal element to desired value.
   * @param keepDiag if @p true, @p diagValue is ignored and original diagonal is preserved
   * @return original diagonal value if matrix is square; zero otherwise
   *
   * @note @p diagValue and @p keepDiag are ignored if the matrix is not square
   */
  virtual real64
  clearRow( globalIndex const row,
            bool const keepDiag = false,
            real64 const diagValue = 0.0 ) = 0;

  /**
   * @brief Add entries of another matrix to this.
   * @param src   the source matrix
   * @param scale factor to scale entries of @p src by
   *
   * @note Sparsity pattern of @p this must be a superset of sparsity of @p src.
   *       @p this and @p src must have the same parallel row distribution.
   */
  virtual void
  addEntries( Matrix const & src, real64 const scale = 1.0 ) = 0;

  /**
   * @brief Add entries of a vector to the diagonal of this matrix.
   * @param src the source vector
   *
   * @note @p this must be square and have a (possibly zero) diagonal entry in every row.
   *       @p this and @p src must have the same parallel row distribution.
   */
  virtual void
  addDiagonal( Vector const & src ) = 0;

  ///@}

  /**
   * @name Accessors Methods
   */
  ///@{

  /**
   * @brief Returns the number of nonzero entries in the longest row of the matrix.
   * @return the max length of a row
   *
   * Collective.
   */
  virtual localIndex
  maxRowLength() const = 0;

  /**
   * @brief Get row length via local row index.
   * @param[in] localRowIndex the local row index
   * @return the number of nonzero entries in the row
   *
   * TODO: Breaks the goal of hiding local row indexing from user.
   *       Revise use cases to use ilower() and iupper().
   */
  virtual localIndex
  localRowLength( localIndex localRowIndex ) const = 0;

  /**
   * @brief Get row length via global row index.
   * @param[in] globalRowIndex the global row index
   * @return the number of nonzero entries in the row
   */
  virtual localIndex
  globalRowLength( globalIndex const globalRowIndex ) const = 0;

  /**
   * @brief Returns a copy of the data in row @p globalRow.
   * @param[in]  globalRow  the index of global row to extract
   * @param[out] colIndices the output array of global column indices (must have a large enough size)
   * @param[out] values     the output array of values (must have a large enough size)
   */
  virtual void
  getRowCopy( globalIndex const globalRow,
              arraySlice1d< globalIndex > const & colIndices,
              arraySlice1d< real64 > const & values ) const = 0;

  /**
   * @brief get diagonal element value on a given row
   * @param globalRow global row index
   * @return value of diagonal element on the row
   */
  virtual real64
  getDiagValue( globalIndex globalRow ) const = 0;

  /**
   * @brief Extract diagonal values into a vector.
   * @param dst the target vector, must have the same row partitioning as @p this
   */
  virtual void
  extractDiagonal( Vector & dst ) const = 0;

  /**
   * @brief Returns the number of global rows.
   * @return number of global rows
   */
  virtual globalIndex
  numGlobalRows() const override = 0;

  /**
   * @brief Returns the number of global columns.
   * @return number of global columns
   */
  virtual globalIndex
  numGlobalCols() const override = 0;

  /**
   * @brief Return the local number of columns on each processor.
   * @return number of local columns
   */
  virtual localIndex
  numLocalRows() const = 0;

  /**
   * @brief Return the local number of columns on each processor.
   * @return number of local rows
   *
   * @note Matrix implementations don't physically "own" columns ranges as they do rows.
   * Instead, the local column range refers to the "diagonal" block of columns which would
   * correspond to the local range of entries of a source vector created with the same
   * local/global size as the number of matrix columns.
   */
  virtual localIndex
  numLocalCols() const = 0;

  /**
   * @brief Returns the index of the first global row owned by that processor.
   * @return the index of the first global row owned by that processor
   */
  virtual globalIndex
  ilower() const = 0;

  /**
   * @brief Returns index one past the last global row owned by that processor.
   * @return the next index after last global row owned by that processor
   *
   * @note The intention is for [ilower; iupper) to be used as a half-open index range
   */
  virtual globalIndex
  iupper() const = 0;

  /**
   * @brief Returns the index of the first global col owned by that processor.
   * @return index of the first owned global col
   *
   * @note Matrix implementations don't physically "own" column ranges the same way
   * they do row ranges. Instead, the column range refers to the "diagonal" block of
   * columns which would correspond to the local range of entries of a vector created
   * with the same local/global size as the number of matrix columns.
   */
  virtual globalIndex
  jlower() const = 0;

  /**
   * @brief Returns index one past the last global col owned by that processor.
   * @return index one past the last owned global col
   *
   * @note The intention is for [jlower; jupper) to be used as a half-open index range.
   * @note Also see note for @p jlower() about the meaning of "owned" columns.
   */
  virtual globalIndex
  jupper() const = 0;

  /**
   * @brief Returns the number of nonzeros in the local portion of the matrix
   * @return the number of nonzeros in the local portion of the matrix
   */
  virtual localIndex
  numLocalNonzeros() const = 0;

  /**
   * @brief Returns the total number of nonzeros in the matrix
   * @return the total number of nonzeros in the matrix
   */
  virtual globalIndex
  numGlobalNonzeros() const = 0;

  /**
   * @brief Returns the infinity norm of the matrix.
   * @return the value of infinity norm
   */
  virtual real64
  normInf() const = 0;

  /**
   * @brief Returns the one norm of the matrix.
   * @return the value of 1-norm
   */
  virtual real64
  norm1() const = 0;

  /**
   * @brief Returns the Frobenius norm of the matrix.
   * @return the value of Frobenius norm
   */
  virtual real64
  normFrobenius() const = 0;

  /**
   * @brief Map a global row index to local row index
   * @param index the global row index
   * @return the local row index corresponding to @p index, or -1 if not a local row
   */
  virtual localIndex
  getLocalRowID( globalIndex const index ) const = 0;

  /**
   * @brief Map a local row index to global row index
   * @param index the local row index (between 0 and number of local rows)
   * @return the global row index corresponding to @p index
   */
  virtual globalIndex
  getGlobalRowID( localIndex const index ) const = 0;

  /**
   * @brief Get the MPI communicator the matrix was created with
   * @return MPI communicator passed in @p create...()
   *
   * @note when build without MPI, may return anything
   *       (MPI_Comm will be a mock type defined in MpiWrapper)
   */
  virtual MPI_Comm
  getComm() const = 0;

  ///@}

  /**
   * @name I/O Methods
   */
  ///@{

  /**
   * @brief Print the matrix in Trilinos format to a stream.
   * @param os the output stream
   */
  virtual void
  print( std::ostream & os = std::cout ) const = 0;

  /**
   * @brief Write the matrix to filename in a matlab-compatible format.
   * @param filename name of the output file
   * @param format   output format
   *
   * Within octave / matlab:
   * >> load filename
   * >> M = spconvert(filename_root)
   */
  virtual void
  write(
    string const & filename,
    LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const = 0;

  ///@}

  /**
   * @brief Stream insertion operator for all matrix types
   * @param os the output stream
   * @param matrix the matrix to be printed
   * @return reference to the output stream
   */
  friend std::ostream &
  operator<<( std::ostream & os, Matrix const & matrix )
  {
    matrix.print( os );
    return os;
  }

  /// Flag indicating whether the matrix is currently open for adding new entries
  bool m_closed;

  /// Flag indicating whether the matrix (sparsity pattern) has been assembled
  bool m_assembled;
};

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_MATRIXBASE_HPP_
