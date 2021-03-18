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
 * @file HypreMatrix.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/hypre/HypreVector.hpp"
#include "linearAlgebra/interfaces/hypre/HypreExport.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/interfaces/MatrixBase.hpp"

/**
 * @name Hypre forward declarations.
 *
 * Forward declare hypre's matrix structs and pointer aliases in order
 * to avoid including hypre headers and leaking into the rest of GEOSX.
 */
///@{

/// IJMatrix struct forward declaration
extern "C" struct hypre_IJMatrix_struct;

/// ParCSRMatrix struct forward declaration
extern "C" struct hypre_ParCSRMatrix_struct;

///@}

namespace geosx
{

/**
 * @brief Wrapper class for hypre's ParCSRMatrix.
 *
 * This class creates and provides basic support for the HYPRE_ParCSRMatrix object
 * type used in Hypre using the linear-algebraic system interface (IJ interface).
 */
class HypreMatrix final : public virtual LinearOperator< HypreVector >,
  private MatrixBase< HypreMatrix, HypreVector >
{
public:

  /// Compatible vector type
  using Vector = HypreVector;

  /// Associated exporter type
  using Export = HypreExport;

  /// IJMatrix pointer alias
  using HYPRE_IJMatrix = hypre_IJMatrix_struct *;

  /// ParCSRMatrix pointer alias
  using HYPRE_ParCSRMatrix = hypre_ParCSRMatrix_struct *;

  /**
   * @name Constructor/Destructor Methods
   */
  ///@{

  /**
   * @brief Empty matrix constructor.
   *
   * Create an empty (distributed) matrix.
   */
  HypreMatrix();

  /**
   * @brief Copy constructor.
   * @param[in] src matrix to be copied
   */
  HypreMatrix( HypreMatrix const & src );

  /**
   * @brief Move constructor.
   * @param src matrix to be moved from
   */
  HypreMatrix( HypreMatrix && src ) noexcept;

  /**
   * @brief Copy assignment.
   * @param src matrix to be copied
   * @return the new vector
   */
  HypreMatrix & operator=( HypreMatrix const & src );

  /**
   * @brief Move assignment.
   * @param src matrix to be moved from
   * @return the new vector
   */
  HypreMatrix & operator=( HypreMatrix && src ) noexcept;

  /**
   * @brief Virtual destructor.
   */
  ~HypreMatrix() override;

  ///@}

  /**
   * @name MatrixBase interface
   */
  ///@{

  using MatrixBase::createWithLocalSize;
  using MatrixBase::createWithGlobalSize;
  using MatrixBase::closed;
  using MatrixBase::assembled;
  using MatrixBase::insertable;
  using MatrixBase::modifiable;
  using MatrixBase::ready;
  using MatrixBase::residual;
  using MatrixBase::setDofManager;
  using MatrixBase::dofManager;
  using MatrixBase::getInvRowSums;
  using MatrixBase::create;

  virtual void create( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                       localIndex const numLocalColumns,
                       MPI_Comm const & comm ) override;

  virtual void createWithLocalSize( localIndex const localRows,
                                    localIndex const localCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm ) override;

  virtual void createWithGlobalSize( globalIndex const globalRows,
                                     globalIndex const globalCols,
                                     localIndex const maxEntriesPerRow,
                                     MPI_Comm const & comm ) override;

  virtual void open() override;

  virtual void close() override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreMatrix>::created
   */
  virtual bool created() const override;

  virtual void reset() override;

  virtual void set( real64 const value ) override;

  virtual void zero() override;

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

  virtual void apply( HypreVector const & src,
                      HypreVector & dst ) const override;

  virtual void applyTranspose( Vector const & src,
                               Vector & dst ) const override;

  virtual void multiply( HypreMatrix const & src,
                         HypreMatrix & dst ) const override;

  virtual void leftMultiplyTranspose( HypreMatrix const & src,
                                      HypreMatrix & dst ) const override;

  virtual void rightMultiplyTranspose( HypreMatrix const & src,
                                       HypreMatrix & dst ) const override;

  virtual void multiplyRAP( HypreMatrix const & R,
                            HypreMatrix const & P,
                            HypreMatrix & dst ) const override;

  virtual void multiplyPtAP( HypreMatrix const & P,
                             HypreMatrix & dst ) const override;

  virtual void gemv( real64 const alpha,
                     HypreVector const & x,
                     real64 const beta,
                     HypreVector & y,
                     bool useTranspose = false ) const override;

  virtual void scale( real64 const scalingFactor ) override;

  virtual void leftScale( HypreVector const & vec ) override;

  virtual void rightScale( HypreVector const & vec ) override;

  virtual void leftRightScale( HypreVector const & vecLeft,
                               HypreVector const & vecRight ) override;

  virtual void transpose( HypreMatrix & dst ) const override;

  virtual real64 clearRow( globalIndex const row,
                           bool const keepDiag = false,
                           real64 const diagValue = 0.0 ) override;

  virtual void addEntries( HypreMatrix const & src,
                           real64 const scale = 1.0,
                           bool const samePattern = true ) override;

  virtual void addDiagonal( HypreVector const & src ) override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::maxRowLength
   */
  virtual localIndex maxRowLength() const override;

  virtual localIndex rowLength( globalIndex const globalRowIndex ) const override;

  virtual void getRowLengths( arrayView1d< localIndex > const & lengths ) const override;

  virtual void getRowCopy( globalIndex globalRowIndex,
                           arraySlice1d< globalIndex > const & colIndices,
                           arraySlice1d< real64 > const & values ) const override;

  virtual void extractDiagonal( HypreVector & dst ) const override;

  virtual void getRowSums( HypreVector & dst ) const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numGlobalRows
   */
  virtual globalIndex numGlobalRows() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numGlobalCols
   */
  virtual globalIndex numGlobalCols() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numLocalRows
   */
  virtual localIndex numLocalRows() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numLocalCols
   */
  virtual localIndex numLocalCols() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::ilower
   */
  virtual globalIndex ilower() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::iupper
   */
  virtual globalIndex iupper() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::jlower
   */
  virtual globalIndex jlower() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::jupper
   */
  virtual globalIndex jupper() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numLocalNonzeros
   */
  virtual localIndex numLocalNonzeros() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::numGlobalNonzeros
   */
  virtual globalIndex numGlobalNonzeros() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::normInf
   */
  virtual real64 normInf() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::norm1
   */
  virtual real64 norm1() const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::normFrobenius
   */
  virtual real64 normFrobenius() const override;

  virtual localIndex getLocalRowID( globalIndex const index ) const override;

  virtual globalIndex getGlobalRowID( localIndex const index ) const override;

  /**
   * @copydoc MatrixBase<HypreMatrix,HypreVector>::getComm
   */
  virtual MPI_Comm getComm() const override;

  virtual void print( std::ostream & os = std::cout ) const override;

  virtual void write( string const & filename,
                      LAIOutputFormat const format = LAIOutputFormat::MATRIX_MARKET ) const override;

  ///@}

  /**
   * @brief Returns a pointer to implementation.
   * @return the underlying HYPRE_ParCSRMatrix object.
   */
  HYPRE_ParCSRMatrix const & unwrapped() const;

  /**
   * @brief Returns a pointer to implementation.
   * @return the underlying HYPRE_IJMatrix object.
   */
  HYPRE_IJMatrix const & unwrappedIJ() const;

private:

  /**
   * @brief Perform a matrix matrix product with Parallel Matrix
   */
  void parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix );

  /**
   * Pointer to underlying HYPRE_IJMatrix type.
   */
  HYPRE_IJMatrix m_ij_mat{};

  /**
   * Pointer to underlying HYPRE_ParCSRMatrix type.
   */
  HYPRE_ParCSRMatrix m_parcsr_mat{};

};

} // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREMATRIX_HPP_*/
