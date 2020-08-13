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
 * @file EpetraMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraVector.hpp"
#include "linearAlgebra/interfaces/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

class Epetra_Map;

class Epetra_CrsMatrix;

class Epetra_FECrsMatrix;

namespace geosx
{
/**
 * @brief Wrapper class for Epetra's CrsMatrix.
 */
class EpetraMatrix final : public virtual LinearOperator< EpetraVector >,
                           private MatrixBase< EpetraMatrix, EpetraVector >
{
public:
  /// Compatible vector type
  using Vector = EpetraVector;

  /**
   * @name Constructor/Destructor methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  EpetraMatrix();

  /**
   * @brief Copy constructor.
   * @param[in] src the matrix to be copied
   */
  EpetraMatrix( EpetraMatrix const & src );

  /**
   * @brief Virtual destructor.
   */
  virtual ~EpetraMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::assembled;
  using MatrixBase::closed;
  using MatrixBase::create;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::createWithLocalSize;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;
  using MatrixBase::residual;

  virtual void
  createWithLocalSize( localIndex const localRows,
                       localIndex const localCols,
                       localIndex const maxEntriesPerRow,
                       MPI_Comm const & comm ) override;

  virtual void
  createWithGlobalSize( globalIndex const globalRows,
                        globalIndex const globalCols,
                        localIndex const maxEntriesPerRow,
                        MPI_Comm const & comm ) override;

  virtual void
  open() override;

  virtual void
  close() override;

  virtual bool
  created() const override;

  virtual void
  reset() override;

  virtual void
  set( real64 const value ) override;

  virtual void
  zero() override;

  virtual void
  add( globalIndex const rowIndex,
       globalIndex const colIndex,
       real64 const value ) override;

  virtual void
  set( globalIndex const rowIndex,
       globalIndex const colIndex,
       real64 const value ) override;

  virtual void
  insert( globalIndex const rowIndex,
          globalIndex const colIndex,
          real64 const value ) override;

  virtual void
  add( globalIndex const rowIndex,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const size ) override;

  virtual void
  set( globalIndex const rowIndex,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const size ) override;

  virtual void
  insert( globalIndex const rowIndex,
          globalIndex const * colIndices,
          real64 const * values,
          localIndex const size ) override;

  virtual void
  add( globalIndex const rowIndex,
       arraySlice1d< globalIndex const > const & colIndices,
       arraySlice1d< real64 const > const & values ) override;

  virtual void
  set( globalIndex const rowIndex,
       arraySlice1d< globalIndex const > const & colIndices,
       arraySlice1d< real64 const > const & values ) override;

  virtual void
  insert( globalIndex const rowIndex,
          arraySlice1d< globalIndex const > const & colIndices,
          arraySlice1d< real64 const > const & values ) override;

  virtual void
  add(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void
  set(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void
  insert(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values ) override;

  virtual void
  add(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void
  set(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void
  insert(
    arraySlice1d< globalIndex const > const & rowIndices,
    arraySlice1d< globalIndex const > const & colIndices,
    arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values ) override;

  virtual void
  add( globalIndex const * rowIndices,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const numRows,
       localIndex const numCols ) override;

  virtual void
  set( globalIndex const * rowIndices,
       globalIndex const * colIndices,
       real64 const * values,
       localIndex const numRows,
       localIndex const numCols ) override;

  virtual void
  insert( globalIndex const * rowIndices,
          globalIndex const * colIndices,
          real64 const * values,
          localIndex const numRows,
          localIndex const numCols ) override;

  virtual void
  apply( EpetraVector const & src, EpetraVector & dst ) const override;

  virtual void
  applyTranspose( EpetraVector const & src,
                  EpetraVector & dst ) const override;

  virtual void
  multiply( EpetraMatrix const & src, EpetraMatrix & dst ) const override;

  virtual void
  leftMultiplyTranspose( EpetraMatrix const & src,
                         EpetraMatrix & dst ) const override;

  virtual void
  rightMultiplyTranspose( EpetraMatrix const & src,
                          EpetraMatrix & dst ) const override;

  virtual void
  multiplyRAP( EpetraMatrix const & R,
               EpetraMatrix const & P,
               EpetraMatrix & dst ) const override;

  virtual void
  multiplyPtAP( EpetraMatrix const & P,
                EpetraMatrix & dst ) const override;

  virtual void
  gemv( real64 const alpha,
        EpetraVector const & x,
        real64 const beta,
        EpetraVector & y,
        bool useTranspose = false ) const override;

  virtual void
  scale( real64 const scalingFactor ) override;

  virtual void
  leftScale( EpetraVector const & vec ) override;

  virtual void
  rightScale( EpetraVector const & vec ) override;

  virtual void
  leftRightScale( EpetraVector const & vecLeft,
                  EpetraVector const & vecRight ) override;

  virtual void
  transpose( EpetraMatrix & dst ) const override;

  virtual real64
  clearRow( globalIndex const row,
            bool const keepDiag = false,
            real64 const diagValue = 0.0 ) override;

  virtual void
  addEntries( EpetraMatrix const & src,
              real64 const scale = 1.0 ) override;

  virtual void
  addDiagonal( EpetraVector const & src ) override;

  virtual localIndex
  maxRowLength() const override;

  virtual localIndex
  localRowLength( localIndex localRowIndex ) const override;

  virtual localIndex
  globalRowLength( globalIndex globalRowIndex ) const override;

  virtual void
  getRowCopy( globalIndex globalRow,
              arraySlice1d< globalIndex > const & colIndices,
              arraySlice1d< real64 > const & values ) const override;

  virtual real64
  getDiagValue( globalIndex globalRow ) const override;

  virtual void
  extractDiagonal( EpetraVector & dst ) const override;

  virtual globalIndex
  numGlobalRows() const override;

  virtual globalIndex
  numGlobalCols() const override;

  virtual localIndex
  numLocalRows() const override;

  virtual localIndex
  numLocalCols() const override;

  virtual globalIndex
  ilower() const override;

  virtual globalIndex
  iupper() const override;

  virtual globalIndex
  jlower() const override;

  virtual globalIndex
  jupper() const override;

  virtual localIndex
  numLocalNonzeros() const override;

  virtual globalIndex
  numGlobalNonzeros() const override;

  virtual real64
  normInf() const override;

  virtual real64
  norm1() const override;

  virtual real64
  normFrobenius() const override;

  virtual localIndex
  getLocalRowID( globalIndex const index ) const override;

  virtual globalIndex
  getGlobalRowID( localIndex const index ) const override;

  virtual MPI_Comm
  getComm() const override;

  virtual void
  print( std::ostream & os = std::cout ) const override;

  virtual void
  write( string const & filename,
         LAIOutputFormat const format =
           LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying matrix.
   * @return const pointer to the underlying matrix
   */
  Epetra_FECrsMatrix const &
  unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   * @return non-const pointer to the underlying matrix
   */
  Epetra_FECrsMatrix &
  unwrapped();

private:
  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void
  multiply( bool const transA,
            EpetraMatrix const & B,
            bool const transB,
            EpetraMatrix & C ) const;

  /**
   * @brief Create the matrix by copying data from an Epetra_CrsMatrix
   * @param src the source matrix
   */
  void
  create( Epetra_CrsMatrix const & src );

  /// Pointer to the underlying Epetra_CrsMatrix.
  std::unique_ptr< Epetra_FECrsMatrix > m_matrix;

  /// Map representing the parallel partitioning of a source vector (x in y=Ax)
  std::unique_ptr< Epetra_Map > m_src_map;

  /// Map representing the parallel partitioning of a destination vector (y in y=Ax)
  std::unique_ptr< Epetra_Map > m_dst_map;
};

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_EPETRAMATRIX_HPP_*/
