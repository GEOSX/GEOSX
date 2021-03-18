/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PetscMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/petsc/PetscVector.hpp"
#include "linearAlgebra/interfaces/petsc/PetscExport.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

/**
 * @name PETSc forward declarations.
 *
 * Forward declare PETSc's matrix struct and pointer aliases in order
 * to avoid including PETSc headers and leaking into the rest of GEOSX.
 */
///@{

/// Mat struct forward declaration
extern "C" struct _p_Mat;

///@}

namespace geosx
{

/**
 * @brief This class creates and provides basic support for the Mat
 *        matrix object type used in PETSc.
 */
class PetscMatrix final : public virtual LinearOperator< PetscVector >,
  private MatrixBase< PetscMatrix, PetscVector >
{
public:

  /// Compatible vector type
  using Vector = PetscVector;

  /// Associated exporter type
  using Export = PetscExport;

  /// Alias for PETSc matrix struct pointer
  using Mat = struct _p_Mat *;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   */
  PetscMatrix();

  /**
   * @brief Copy constructor.
   * @param[in] src the matrix to be copied
   */
  PetscMatrix( PetscMatrix const & src );

  /**
   * @brief Move constructor.
   * @param[in] src the matrix to be copied
   */
  PetscMatrix( PetscMatrix && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src matrix to be copied.
   * @return the new vector.
   */
  PetscMatrix & operator=( PetscMatrix const & src );

  /**
   * @brief Move assignment.
   * @param src matrix to be moved from.
   * @return the new matrix.
   */
  PetscMatrix & operator=( PetscMatrix && src ) noexcept;

  /**
   * @brief Destructor.
   */
  ~PetscMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::create;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;
  using MatrixBase::residual;
  using MatrixBase::setDofManager;
  using MatrixBase::dofManager;
  using MatrixBase::getInvRowSums;

  virtual void createWithLocalSize( localIndex const localRows,
                                    localIndex const localCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm ) override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numGlobalRows
   */
  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

  virtual void open() override;

  virtual void close() override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const colIndex,
                    real64 const value ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value ) override;

  virtual void add( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void set( globalIndex const rowIndex,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const size ) override;

  virtual void insert( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const size ) override;

  virtual void add( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void set( globalIndex const rowIndex,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice1d< real64 const > const & values ) override;

  virtual void insert( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values ) override;

  virtual void add( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const > const & values ) override;

  virtual void set( arraySlice1d< globalIndex const > const & rowIndices,
                    arraySlice1d< globalIndex const > const & colIndices,
                    arraySlice2d< real64 const > const & values ) override;

  virtual void insert( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values ) override;

  virtual void add( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void set( globalIndex const * rowIndices,
                    globalIndex const * colIndices,
                    real64 const * values,
                    localIndex const numRows,
                    localIndex const numCols ) override;

  virtual void insert( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols ) override;

  virtual void apply( PetscVector const & src,
                      PetscVector & dst ) const override;

  virtual void multiply( PetscMatrix const & src,
                         PetscMatrix & dst ) const override;

  virtual void applyTranspose( Vector const & src,
                               Vector & dst ) const override;

  virtual void leftMultiplyTranspose( PetscMatrix const & src,
                                      PetscMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( PetscMatrix const & src,
                                       PetscMatrix & dst ) const override;

  virtual void multiplyRAP( PetscMatrix const & R,
                            PetscMatrix const & P,
                            PetscMatrix & dst ) const override;

  virtual void multiplyPtAP( PetscMatrix const & P,
                             PetscMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     PetscVector const & x,
                     real64 const beta,
                     PetscVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( PetscVector const & vec ) override;

  virtual void rightScale( PetscVector const & vec ) override;

  virtual void leftRightScale( PetscVector const & vecLeft,
                               PetscVector const & vecRight ) override;

  virtual void transpose( PetscMatrix & dst ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual void addEntries( PetscMatrix const & src,
                           real64 const scale = 1.0,
                           bool const samePattern = true ) override;

  virtual void addDiagonal( PetscVector const & src ) override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::maxRowLength
   */
  virtual localIndex maxRowLength() const override;

  virtual localIndex rowLength( globalIndex const globalRowIndex ) const override;

  virtual void getRowLengths( arrayView1d< localIndex > const & lengths ) const override;

  virtual void extractDiagonal( PetscVector & dst ) const override;

  virtual void getRowSums( PetscVector & dst ) const override;

  virtual void getRowCopy( globalIndex globalRow,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numGlobalRows
   */
  virtual globalIndex numGlobalRows() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numGlobalCols
   */
  virtual globalIndex numGlobalCols() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numLocalRows
   */
  virtual localIndex numLocalRows() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numLocalCols
   */
  virtual localIndex numLocalCols() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::iupper
   */
  virtual globalIndex iupper() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::jlower
   */
  virtual globalIndex jlower() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::jupper
   */
  virtual globalIndex jupper() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numLocalNonzeros
   */
  virtual localIndex numLocalNonzeros() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::numGlobalNonzeros
   */
  virtual globalIndex numGlobalNonzeros() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::normFrobenius
   */
  virtual real64 normFrobenius() const override;

  virtual localIndex getLocalRowID( globalIndex const index ) const override;

  virtual globalIndex getGlobalRowID( localIndex const index ) const override;

  /**
   * @copydoc MatrixBase<PetscMatrix,PetscVector>::getComm
   */
  virtual MPI_Comm getComm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a const pointer to the underlying matrix.
   * @return the const pointer to the underlying matrix.
   */
  const Mat & unwrapped() const;

  /**
   * @brief Returns a non-const pointer to the underlying matrix.
   * @return the non-const pointer to the underlying matrix.
   */
  Mat & unwrapped();

private:

  /// Underlying Petsc object.
  Mat m_mat{};

  /// Indices of rows to be cleared on next close()
  array1d< globalIndex > m_rowsToClear;

  /// Diagonal values of rows to be set on next close()
  array1d< real64 > m_diagValues;

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_PETSCSPARSEMATRIX_HPP_*/
