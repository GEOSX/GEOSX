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
 * @file HypreMatrix.cpp
 */

#include "HypreMatrix.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "LvArray/src/output.hpp"

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

#include <iomanip>
#include <numeric>

namespace geosx
{

// Helper function that performs the following sequence of IJMatrix
// call: Create, SetObjectType, Initialize.
static void initialize( MPI_Comm const & comm,
                        HYPRE_BigInt const & ilower,
                        HYPRE_BigInt const & iupper,
                        HYPRE_BigInt const & jlower,
                        HYPRE_BigInt const & jupper,
                        arrayView1d< HYPRE_Int const > const & ncols,
                        HYPRE_IJMatrix & ij_matrix )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixCreate( comm,
                                               ilower,
                                               iupper,
                                               jlower,
                                               jupper,
                                               &ij_matrix ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetObjectType( ij_matrix, HYPRE_PARCSR ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetRowSizes( ij_matrix, ncols.data() ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( ij_matrix ) );
}

HypreMatrix::HypreMatrix()
  : LinearOperator(),
  MatrixBase()
{}

HypreMatrix::HypreMatrix( HypreMatrix const & src )
  : HypreMatrix()
{
  *this = src;
}

HypreMatrix::HypreMatrix( HypreMatrix && src ) noexcept
  : HypreMatrix()
{
  *this = std::move( src );
}

HypreMatrix & HypreMatrix::operator=( HypreMatrix const & src )
{
  if( &src != this )
  {
    reset();
    if( src.ready() )
    {
      // Copy parcsr matrix
      HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixClone( src.m_parcsr_mat, 1 );
      // Create IJ layer (with matrix closed)
      parCSRtoIJ( dst_parcsr );
    }
    m_dofManager = src.dofManager();
  }
  return *this;
}

HypreMatrix & HypreMatrix::operator=( HypreMatrix && src ) noexcept
{
  if( &src != this )
  {
    std::swap( m_ij_mat, src.m_ij_mat );
    std::swap( m_parcsr_mat, src.m_parcsr_mat );
    MatrixBase::operator=( std::move( src ) );
  }
  return *this;
}

HypreMatrix::~HypreMatrix()
{
  reset();
}

void HypreMatrix::createWithGlobalSize( globalIndex const globalRows,
                                        globalIndex const globalCols,
                                        localIndex const maxEntriesPerRow,
                                        MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_Int const rank  = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commRank( comm ) );
  HYPRE_Int const nproc = LvArray::integerConversion< HYPRE_Int >( MpiWrapper::commSize( comm ) );

  HYPRE_Int const localRowSize = LvArray::integerConversion< HYPRE_Int >( globalRows / nproc );
  HYPRE_Int const rowResidual = LvArray::integerConversion< HYPRE_Int >( globalRows % nproc );

  HYPRE_Int const localColSize = LvArray::integerConversion< HYPRE_Int >( globalCols / nproc );
  HYPRE_Int const colResidual = LvArray::integerConversion< HYPRE_Int >( globalCols % nproc );

  HYPRE_BigInt const ilower = rank * localRowSize + ( rank == 0 ? 0 : rowResidual );
  HYPRE_BigInt const iupper = ilower + localRowSize + ( rank == 0 ? rowResidual : 0 ) - 1;
  HYPRE_BigInt const jlower = rank * localColSize + ( rank == 0 ? 0 : colResidual );
  HYPRE_BigInt const jupper = jlower + localColSize + ( rank == 0 ? colResidual : 0 ) - 1;

  array1d< HYPRE_Int > row_sizes;
  row_sizes.resizeDefault( LvArray::integerConversion< localIndex >( iupper - ilower + 1 ),
                           LvArray::integerConversion< HYPRE_Int >( maxEntriesPerRow ) );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::create( CRSMatrixView< real64 const, globalIndex const > const & localMatrix,
                          localIndex const numLocalColumns,
                          MPI_Comm const & comm )
{
  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, localIndex > maxRowEntries( 0 );
  forAll< hypre::execPolicy >( localMatrix.numRows(),
                               [localMatrix, maxRowEntries] GEOSX_HYPRE_HOST_DEVICE ( localIndex const row )
  {
    maxRowEntries.max( localMatrix.numNonZeros( row ) );
  } );

  createWithLocalSize( localMatrix.numRows(), numLocalColumns, maxRowEntries.get(), comm );
  globalIndex const rankOffset = ilower();

  array1d< HYPRE_BigInt > rows;
  rows.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  array1d< HYPRE_Int > sizes;
  sizes.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  array1d< HYPRE_Int > offsets;
  offsets.resizeWithoutInitializationOrDestruction( hypre::memorySpace, localMatrix.numRows() );

  forAll< hypre::execPolicy >( localMatrix.numRows(),
                               [localMatrix, rankOffset,
                                rowsV = rows.toView(),
                                sizesV = sizes.toView(),
                                offsetsV = offsets.toView()] GEOSX_HYPRE_HOST_DEVICE ( localIndex const row )
  {
    rowsV[row] = LvArray::integerConversion< HYPRE_BigInt >( row + rankOffset );
    sizesV[row] = LvArray::integerConversion< HYPRE_Int >( localMatrix.numNonZeros( row ) );
    offsetsV[row] = LvArray::integerConversion< HYPRE_Int >( localMatrix.getOffsets()[row] );
  } );

  // This is necessary so that localMatrix.getColumns() and localMatrix.getEntries() return the device pointers.
  localMatrix.move( hypre::memorySpace, false );

  open();
  GEOSX_HYPRE_CHECK_DEVICE_ERRORS( "Cuda error detected before HYPRE_IJMatrixAddToValues2" );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues2( m_ij_mat,
                                                     localMatrix.numRows(),
                                                     sizes.data(),
                                                     rows.data(),
                                                     offsets.data(),
                                                     localMatrix.getColumns(),
                                                     localMatrix.getEntries() ) );
  close();
}

void HypreMatrix::createWithLocalSize( localIndex const localRows,
                                       localIndex const localCols,
                                       localIndex const maxEntriesPerRow,
                                       MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  HYPRE_BigInt const ilower = MpiWrapper::prefixSum< HYPRE_BigInt >( localRows );
  HYPRE_BigInt const iupper = ilower + localRows - 1;

  HYPRE_BigInt const jlower = MpiWrapper::prefixSum< HYPRE_BigInt >( localCols );
  HYPRE_BigInt const jupper = jlower + localCols - 1;

  array1d< HYPRE_Int > row_sizes;
  row_sizes.resizeDefault( localRows, LvArray::integerConversion< HYPRE_Int >( maxEntriesPerRow ) );

  initialize( comm,
              ilower,
              iupper,
              jlower,
              jupper,
              row_sizes,
              m_ij_mat );
}

void HypreMatrix::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  open();
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetConstantValues( m_ij_mat, static_cast< HYPRE_Real >( value ) ) );
  close();
}

void HypreMatrix::reset()
{
  MatrixBase::reset();
  if( m_ij_mat )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixDestroy( m_ij_mat ) );
    m_ij_mat = nullptr;
    m_parcsr_mat = nullptr;
  }
}

void HypreMatrix::zero()
{
  GEOSX_LAI_ASSERT( ready() );
  open();
  GEOSX_LAI_CHECK_ERROR( hypre_IJMatrixSetConstantValuesParCSR( m_ij_mat, 0.0 ) );
  close();
}

void HypreMatrix::open()
{
  GEOSX_LAI_ASSERT( created() && closed() );
  if( m_assembled )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixInitialize( m_ij_mat ) );
  }
  m_closed = false;
}

void HypreMatrix::close()
{
  GEOSX_LAI_ASSERT( !closed() );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAssemble( m_ij_mat ) );

  // Get a reference to the constructed matrix object. Done only on the first
  // assembly call when the sparsity pattern of the matrix is defined.
  if( !m_assembled )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetObject( m_ij_mat, (void * *) &m_parcsr_mat ) );
    // Compute column partitioning if needed
    if( hypre_IJMatrixRowPartitioning( m_ij_mat ) !=
        hypre_IJMatrixColPartitioning( m_ij_mat ) )
    {
      if( !hypre_ParCSRMatrixCommPkg( m_parcsr_mat ) )
      {
        GEOSX_LAI_CHECK_ERROR( hypre_MatvecCommPkgCreate( m_parcsr_mat ) );
      }
    }
  }

  m_closed = true;
  m_assembled = true;
}

bool HypreMatrix::created() const
{
  return m_ij_mat != nullptr;
}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    hypre::toHypreBigInt( &rowIndex ),
                                                    hypre::toHypreBigInt( &colIndex ),
                                                    &value ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const colIndex,
                       real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  hypre::toHypreBigInt( &rowIndex ),
                                                  hypre::toHypreBigInt( &colIndex ),
                                                  &value ) );

}

void HypreMatrix::insert( globalIndex const rowIndex,
                          globalIndex const colIndex,
                          real64 const value )
{
  GEOSX_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = 1;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    hypre::toHypreBigInt( &rowIndex ),
                                                    hypre::toHypreBigInt( &colIndex ),
                                                    &value ) );
}

void HypreMatrix::add( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( size );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    hypre::toHypreBigInt( &rowIndex ),
                                                    hypre::toHypreBigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( size );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  hypre::toHypreBigInt( &rowIndex ),
                                                  hypre::toHypreBigInt( colIndices ),
                                                  values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex0,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex size )
{
  GEOSX_LAI_ASSERT( insertable() );

#if defined(GEOSX_USE_HYPRE_CUDA)
  array1d< globalIndex > rowIndexDevice( 1 );
  array1d< HYPRE_Int > ncolsDevice( 1 );

  rowIndexDevice[0] = rowIndex0;
  ncolsDevice[0] = LvArray::integerConversion< HYPRE_Int >( size );

  rowIndexDevice.move( LvArray::MemorySpace::cuda, false );
  ncolsDevice.move( LvArray::MemorySpace::cuda, false );

  globalIndex const * const rowIndex = rowIndexDevice.data();
  HYPRE_Int * const ncols = ncolsDevice.data();
#else
  globalIndex const * const rowIndex = &rowIndex0;
  HYPRE_Int hypreSize = size;
  HYPRE_Int * const ncols = &hypreSize;
#endif

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    ncols,
                                                    rowIndex,
                                                    hypre::toHypreBigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::add( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    hypre::toHypreBigInt( &rowIndex ),
                                                    hypre::toHypreBigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::set( globalIndex const rowIndex,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( rowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), rowIndex );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixSetValues( m_ij_mat,
                                                  1,
                                                  &ncols,
                                                  hypre::toHypreBigInt( &rowIndex ),
                                                  hypre::toHypreBigInt( colIndices ),
                                                  values ) );
}

void HypreMatrix::insert( globalIndex const rowIndex,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );

  HYPRE_Int ncols = LvArray::integerConversion< HYPRE_Int >( colIndices.size() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixAddToValues( m_ij_mat,
                                                    1,
                                                    &ncols,
                                                    hypre::toHypreBigInt( &rowIndex ),
                                                    hypre::toHypreBigInt( colIndices ),
                                                    values ) );
}

void HypreMatrix::add( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    add( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::set( arraySlice1d< globalIndex const > const & rowIndices,
                       arraySlice1d< globalIndex const > const & colIndices,
                       arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( rowIndices.size() ); ++i )
  {
    set( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::insert( arraySlice1d< globalIndex const > const & rowIndices,
                          arraySlice1d< globalIndex const > const & colIndices,
                          arraySlice2d< real64 const > const & values )
{
  for( localIndex i = 0; i < rowIndices.size(); ++i )
  {
    insert( rowIndices[i], colIndices, values[i] );
  }
}

void HypreMatrix::add( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    add( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::set( globalIndex const * rowIndices,
                       globalIndex const * colIndices,
                       real64 const * values,
                       localIndex const numRows,
                       localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    set( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::insert( globalIndex const * rowIndices,
                          globalIndex const * colIndices,
                          real64 const * values,
                          localIndex const numRows,
                          localIndex const numCols )
{
  for( localIndex i = 0; i < numRows; ++i )
  {
    insert( rowIndices[i], colIndices, values + numCols * i, numCols );
  }
}

void HypreMatrix::apply( HypreVector const & src,
                         HypreVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( 1.0,
                                                   m_parcsr_mat,
                                                   src.unwrapped(),
                                                   0.0,
                                                   dst.unwrapped() ) );
}

void HypreMatrix::applyTranspose( Vector const & src,
                                  Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvecT( 1.0,
                                                    m_parcsr_mat,
                                                    src.unwrapped(),
                                                    0.0,
                                                    dst.unwrapped() ) );
}

void HypreMatrix::multiply( HypreMatrix const & src,
                            HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::leftMultiplyTranspose( HypreMatrix const & src,
                                         HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.numGlobalRows() );

  // Compute product
  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParTMatmul( m_parcsr_mat, src.m_parcsr_mat );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::rightMultiplyTranspose( HypreMatrix const & src,
                                          HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.numGlobalCols() );

  // Transpose this
  HypreMatrix tmp;
  transpose( tmp );

  // Compute product
  src.multiply( tmp, dst );
}

void HypreMatrix::multiplyRAP( HypreMatrix const & R,
                               HypreMatrix const & P,
                               HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( R.ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), R.numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  HypreMatrix Rt;
  R.transpose( Rt );

  HYPRE_Int const Rt_owns_its_col_starts = hypre_ParCSRMatrixOwnsColStarts( Rt.unwrapped() );
  HYPRE_Int const P_owns_its_col_starts = hypre_ParCSRMatrixOwnsColStarts( P.unwrapped() );

  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixRAP( Rt.unwrapped(),
                                                               m_parcsr_mat,
                                                               P.unwrapped() );

  hypre_ParCSRMatrixSetRowStartsOwner( dst_parcsr, 0 );
  hypre_ParCSRMatrixSetColStartsOwner( dst_parcsr, 0 );

  hypre_ParCSRMatrixSetColStartsOwner( Rt.unwrapped(), Rt_owns_its_col_starts );
  hypre_ParCSRMatrixSetColStartsOwner( P.unwrapped(), P_owns_its_col_starts );

  dst.parCSRtoIJ( dst_parcsr );
  Rt.reset();
}

void HypreMatrix::multiplyPtAP( HypreMatrix const & P,
                                HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), P.numGlobalRows() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  HYPRE_Int const P_owns_its_col_starts = hypre_ParCSRMatrixOwnsColStarts( P.unwrapped() );

  HYPRE_ParCSRMatrix const dst_parcsr = hypre_ParCSRMatrixRAPKT( P.unwrapped(),
                                                                 m_parcsr_mat,
                                                                 P.unwrapped(),
                                                                 0 );

  hypre_ParCSRMatrixSetRowStartsOwner( dst_parcsr, 0 );
  hypre_ParCSRMatrixSetColStartsOwner( dst_parcsr, 0 );

  hypre_ParCSRMatrixSetColStartsOwner( P.unwrapped(), P_owns_its_col_starts );

  dst.parCSRtoIJ( dst_parcsr );
}

void HypreMatrix::parCSRtoIJ( HYPRE_ParCSRMatrix const & parCSRMatrix )
{
  reset();
  m_closed = false;

  hypre_IJMatrix * const ijmatrix = hypre_CTAlloc( hypre_IJMatrix, 1, HYPRE_MEMORY_HOST );

  hypre_IJMatrixComm( ijmatrix ) = hypre_ParCSRMatrixComm( parCSRMatrix );
  hypre_IJMatrixObject( ijmatrix ) = parCSRMatrix;
  hypre_IJMatrixTranslator( ijmatrix ) = nullptr;
  hypre_IJMatrixAssumedPart( ijmatrix ) = hypre_ParCSRMatrixAssumedPartition( parCSRMatrix );
  hypre_ParCSRMatrixOwnsAssumedPartition( parCSRMatrix ) = 0;

  hypre_IJMatrixAssembleFlag( ijmatrix ) = 1;

  hypre_IJMatrixObjectType( ijmatrix ) = HYPRE_PARCSR;
#ifdef HYPRE_USING_OPENMP
  hypre_IJMatrixOMPFlag( ijmatrix ) = 1;
#else
  hypre_IJMatrixOMPFlag( ijmatrix ) = 0;
#endif
  hypre_IJMatrixPrintLevel( ijmatrix ) = 0;

  array1d< HYPRE_BigInt > info( 2 );
  if( MpiWrapper::commRank( hypre_IJMatrixComm( ijmatrix ) ) == 0 )
  {
    info( 0 ) = hypre_ParCSRMatrixFirstRowIndex( parCSRMatrix );
    info( 1 ) = hypre_ParCSRMatrixFirstColDiag( parCSRMatrix );
  }
  MpiWrapper::bcast( info.data(), 2, 0, hypre_IJMatrixComm( ijmatrix ) );
  hypre_IJMatrixGlobalFirstRow( ijmatrix ) = info( 0 );
  hypre_IJMatrixGlobalFirstCol( ijmatrix ) = info( 1 );

  hypre_IJMatrixGlobalNumRows( ijmatrix ) = hypre_ParCSRMatrixGlobalNumRows( parCSRMatrix );
  hypre_IJMatrixGlobalNumCols( ijmatrix ) = hypre_ParCSRMatrixGlobalNumCols( parCSRMatrix );

  // Set row partitioning
  if( hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix ) )
  {
    hypre_IJMatrixRowPartitioning( ijmatrix ) = hypre_ParCSRMatrixRowStarts( parCSRMatrix );
  }
  else
  {
    HYPRE_BigInt * const row_partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
    row_partitioning[0] = hypre_ParCSRMatrixFirstRowIndex( parCSRMatrix );
    row_partitioning[1] = hypre_ParCSRMatrixLastRowIndex( parCSRMatrix ) + 1;
    hypre_IJMatrixRowPartitioning( ijmatrix ) = row_partitioning;
    hypre_ParCSRMatrixRowStarts( parCSRMatrix ) = row_partitioning;
  }
  hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix ) = 0;

  if( hypre_IJMatrixGlobalNumRows( ijmatrix ) != hypre_IJMatrixGlobalNumCols( ijmatrix ) )
  {
    // Rectangular matrix
    // Set column partitioning
    if( hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) )
    {
      hypre_IJMatrixColPartitioning( ijmatrix ) = hypre_ParCSRMatrixColStarts( parCSRMatrix );
    }
    else
    {
      HYPRE_BigInt * const col_partitioning = hypre_CTAlloc( HYPRE_BigInt, 2, HYPRE_MEMORY_HOST );
      col_partitioning[0] = hypre_ParCSRMatrixFirstColDiag( parCSRMatrix );
      col_partitioning[1] = hypre_ParCSRMatrixLastColDiag( parCSRMatrix ) + 1;
      hypre_IJMatrixColPartitioning( ijmatrix ) = col_partitioning;
      hypre_ParCSRMatrixColStarts( parCSRMatrix ) = col_partitioning;
    }
    hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) = 0;
  }
  else
  {
    // Square matrix
    hypre_IJMatrixColPartitioning( ijmatrix ) = hypre_IJMatrixRowPartitioning( ijmatrix );
    hypre_ParCSRMatrixOwnsColStarts( parCSRMatrix ) = hypre_ParCSRMatrixOwnsRowStarts( parCSRMatrix );
  }

  m_ij_mat = (HYPRE_IJMatrix) ijmatrix;
  close();
}

void HypreMatrix::gemv( real64 const alpha,
                        HypreVector const & x,
                        real64 const beta,
                        HypreVector & y,
                        bool useTranspose ) const
{
  GEOSX_LAI_ASSERT( ready() );

  if( !useTranspose )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvec( alpha,
                                                     m_parcsr_mat,
                                                     x.unwrapped(),
                                                     beta,
                                                     y.unwrapped() ) );
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixMatvecT( alpha,
                                                      m_parcsr_mat,
                                                      x.unwrapped(),
                                                      beta,
                                                      y.unwrapped() ) );
  }
}

void HypreMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  hypre_CSRMatrix const * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );
  HYPRE_Real * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

  forAll< hypre::execPolicy >( diag_nnz, [=] GEOSX_HYPRE_HOST_DEVICE ( HYPRE_Int const i )
  {
    ptr_diag_data[i] *= scalingFactor;
  } );

  hypre_CSRMatrix const * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );
  HYPRE_Real * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

  forAll< hypre::execPolicy >( offdiag_nnz, [=] GEOSX_HYPRE_HOST_DEVICE ( HYPRE_Int const i )
  {
    ptr_offdiag_data[i] *= scalingFactor;
  } );
}

void HypreMatrix::leftScale( HypreVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_EQ( vec.localSize(), numLocalRows() );

  hypre_Vector const * const vec_local = hypre_ParVectorLocalVector( vec.unwrapped() );
  HYPRE_Real const * const vec_data = hypre_VectorData( vec_local );

  hypre_CSRMatrix const * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  hypre_CSRMatrix const * const csr_offdiag = hypre_ParCSRMatrixOffd( m_parcsr_mat );

  HYPRE_Int const * const ia_diag = hypre_CSRMatrixI( csr_diag );
  HYPRE_Int const * const ia_offdiag = hypre_CSRMatrixI( csr_offdiag );

  HYPRE_Real * const va_diag = hypre_CSRMatrixData( csr_diag );
  HYPRE_Real * const va_offdiag = hypre_CSRMatrixData( csr_offdiag );

  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    for( HYPRE_Int j = ia_diag[i]; j < ia_diag[i + 1]; ++j )
    {
      va_diag[j] *= vec_data[i];
    }
    for( HYPRE_Int j = ia_offdiag[i]; j < ia_offdiag[i + 1]; ++j )
    {
      va_offdiag[j] *= vec_data[i];
    }
  } );
}

void HypreMatrix::addEntries( HypreMatrix const & src, real64 const scale, bool samePattern )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == src.numGlobalRows() );
  GEOSX_LAI_ASSERT( numGlobalCols() == src.numGlobalCols() );
  GEOSX_UNUSED_VAR( samePattern );

  HYPRE_ParCSRMatrix parCSRMatrix;
  GEOSX_LAI_CHECK_ERROR( hypre_ParcsrAdd( 1.0,
                                          unwrapped(),
                                          scale,
                                          src.unwrapped(),
                                          &parCSRMatrix ) );

  parCSRtoIJ( parCSRMatrix );
}

void HypreMatrix::addDiagonal( HypreVector const & src )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == numGlobalCols() );
  GEOSX_LAI_ASSERT( numLocalRows() == src.localSize() );

  open();
  real64 const * values = src.extractLocalVector();
  for( globalIndex i = ilower(), ii = 0; i < iupper(); ++i, ++ii )
  {
    add( i, i, values[ii] );
  }
  close();
}

localIndex HypreMatrix::maxRowLength() const
{
  GEOSX_LAI_ASSERT( assembled() );

  hypre_CSRMatrix const * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const * const ia_diag        = hypre_CSRMatrixI( csr_diag );

  hypre_CSRMatrix const * const csr_offd = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const * const ia_offd        = hypre_CSRMatrixI( csr_offd );

  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, localIndex > localMaxRowLength( 0 );
  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const localRow )
  {
    localMaxRowLength.max( (ia_diag[localRow + 1] - ia_diag[localRow]) + (ia_offd[localRow + 1] - ia_offd[localRow] ) );
  } );

  return MpiWrapper::max( localMaxRowLength.get(), getComm() );
}

localIndex HypreMatrix::rowLength( globalIndex const globalRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );

  HYPRE_BigInt row = LvArray::integerConversion< HYPRE_BigInt >( globalRowIndex );
  HYPRE_Int ncols;

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                                     1,
                                                     &row,
                                                     &ncols ) );

  return LvArray::integerConversion< localIndex >( ncols );
}

void HypreMatrix::getRowLengths( arrayView1d< localIndex > const & lengths ) const
{
  GEOSX_LAI_ASSERT( assembled() );

  hypre_CSRMatrix const * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const * const ia_diag        = hypre_CSRMatrixI( csr_diag );

  hypre_CSRMatrix const * const csr_offd = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const * const ia_offd        = hypre_CSRMatrixI( csr_offd );

  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const localRow )
  {
    lengths[localRow] = (ia_diag[localRow + 1] - ia_diag[localRow]) + (ia_offd[localRow + 1] - ia_offd[localRow]);
  } );
}

void HypreMatrix::getRowCopy( globalIndex const globalRowIndex,
                              arraySlice1d< globalIndex > const & colIndices,
                              arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRowIndex, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRowIndex );

  HYPRE_BigInt row = LvArray::integerConversion< HYPRE_BigInt >( globalRowIndex );
  HYPRE_Int numEntries;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetRowCounts( m_ij_mat,
                                                     1,
                                                     &row,
                                                     &numEntries ) );

  GEOSX_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOSX_LAI_ASSERT_GE( values.size(), numEntries );

  GEOSX_LAI_CHECK_ERROR( hypre_IJMatrixGetValuesParCSR( m_ij_mat,
                                                        -1,
                                                        &numEntries,
                                                        &row,
                                                        hypre::toHypreBigInt( colIndices ),
                                                        values ) );
}

void HypreMatrix::extractDiagonal( HypreVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  dst.zero();
  real64 * const values = dst.extractLocalVector();

  hypre_CSRMatrix const * const csr = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const * const ia        = hypre_CSRMatrixI( csr );
  HYPRE_Int const * const ja        = hypre_CSRMatrixJ( csr );
  HYPRE_Real const * const va       = hypre_CSRMatrixData( csr );

  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const localRow )
  {
    for( HYPRE_Int j = ia[localRow]; j < ia[localRow + 1]; ++j )
    {
      if( ja[j] == localRow )
      {
        values[localRow] = va[j];
      }
    }
  } );
}

void HypreMatrix::getRowSums( HypreMatrix::Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  dst.zero();
  real64 * const values = dst.extractLocalVector();

  hypre_CSRMatrix const * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const * const ia_diag        = hypre_CSRMatrixI( csr_diag );
  HYPRE_Real const * const va_diag       = hypre_CSRMatrixData( csr_diag );

  hypre_CSRMatrix const * const csr_offd = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const * const ia_offd        = hypre_CSRMatrixI( csr_offd );
  HYPRE_Real const * const va_offd       = hypre_CSRMatrixData( csr_offd );

  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const localRow )
  {
    for( HYPRE_Int j = ia_diag[localRow]; j < ia_diag[localRow + 1]; ++j )
    {
      values[localRow] += va_diag[j];
    }
    for( HYPRE_Int j = ia_offd[localRow]; j < ia_offd[localRow + 1]; ++j )
    {
      values[localRow] += va_offd[j];
    }
  } );
}

real64 HypreMatrix::clearRow( globalIndex const globalRow,
                              bool const keepDiag,
                              real64 const diagValue )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  // Get local row index
  HYPRE_Int const localRow = LvArray::integerConversion< HYPRE_Int >( getLocalRowID( globalRow ) );

  // Clear row in diagonal block
  hypre_CSRMatrix * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const * const ia_diag  = hypre_CSRMatrixI( csr_diag );
  HYPRE_Int const * const ja_diag  = hypre_CSRMatrixJ( csr_diag );
  HYPRE_Real * const va_diag       = hypre_CSRMatrixData( csr_diag );

  bool const square = numGlobalRows() == numGlobalCols();

  real64 oldDiag = 0.0;
  for( HYPRE_Int j = ia_diag[localRow]; j < ia_diag[localRow + 1]; ++j )
  {
    if( square && ja_diag[j] == localRow )
    {
      oldDiag = va_diag[j];
    }
    va_diag[j] = 0.0;
  }

  // Clear row in off-diagonal block
  hypre_CSRMatrix * const csr_offdiag = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const * const ia_offdiag  = hypre_CSRMatrixI( csr_offdiag );
  HYPRE_Real * const va_offdiag       = hypre_CSRMatrixData( csr_offdiag );

  for( HYPRE_Int j = ia_offdiag[localRow]; j < ia_offdiag[localRow + 1]; ++j )
  {
    va_offdiag[j] = 0.0;
  }

  // Set diagonal value
  real64 const newDiag = keepDiag ? oldDiag : diagValue;
  if( square && std::fabs( newDiag ) > 0.0 )
  {
    set( globalRow, globalRow, newDiag );
  }
  return oldDiag;
}

HYPRE_ParCSRMatrix const & HypreMatrix::unwrapped() const
{
  return m_parcsr_mat;
}

HYPRE_IJMatrix const & HypreMatrix::unwrappedIJ() const
{
  return m_ij_mat;
}

localIndex HypreMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return (index >= ilower && index <= iupper) ? LvArray::integerConversion< localIndex >( index - ilower ) : -1;
}

globalIndex HypreMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( index, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), index );
  return ilower() + index;
}

globalIndex HypreMatrix::numGlobalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumRows( m_ij_mat );
}

globalIndex HypreMatrix::numGlobalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixGlobalNumCols( m_ij_mat );
}

localIndex HypreMatrix::numLocalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< localIndex >( iupper - ilower + 1 );
}

localIndex HypreMatrix::numLocalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< localIndex >( jupper - jlower + 1 );
}

globalIndex HypreMatrix::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< globalIndex >( ilower );
}

globalIndex HypreMatrix::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< globalIndex >( iupper + 1 );
}

globalIndex HypreMatrix::jlower() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< globalIndex >( jlower );
}

globalIndex HypreMatrix::jupper() const
{
  GEOSX_LAI_ASSERT( created() );
  HYPRE_BigInt ilower, iupper, jlower, jupper;
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJMatrixGetLocalRange( m_ij_mat, &ilower, &iupper, &jlower, &jupper ) );
  return LvArray::integerConversion< globalIndex >( jupper + 1 );
}

localIndex HypreMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );

  hypre_CSRMatrix const * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  HYPRE_Int const diag_nnz = hypre_CSRMatrixNumNonzeros( prt_diag_CSR );

  hypre_CSRMatrix const * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
  HYPRE_Int const offdiag_nnz = hypre_CSRMatrixNumNonzeros( prt_offdiag_CSR );

  return LvArray::integerConversion< localIndex >( diag_nnz + offdiag_nnz );
}

globalIndex HypreMatrix::numGlobalNonzeros() const
{
  return MpiWrapper::sum( LvArray::integerConversion< globalIndex >( numLocalNonzeros() ), getComm() );
}

void HypreMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );

  int const this_mpi_process = MpiWrapper::commRank( getComm() );
  int const n_mpi_process = MpiWrapper::commSize( getComm() );
  char str[77];

  if( this_mpi_process == 0 )
  {
    os << "MPI_Process         GlobalRowID         GlobalColID                   Value" << std::endl;
  }

  for( int iRank = 0; iRank < n_mpi_process; iRank++ )
  {
    MpiWrapper::barrier( getComm() );
    if( iRank == this_mpi_process )
    {
      globalIndex const firstRowID = ilower();
      globalIndex const firstDiagColID = jlower();

      hypre_CSRMatrix const * const prt_diag_CSR = hypre_ParCSRMatrixDiag( m_parcsr_mat );
      HYPRE_Int const * const diag_IA = hypre_CSRMatrixI( prt_diag_CSR );
      HYPRE_Int const * const diag_JA = hypre_CSRMatrixJ( prt_diag_CSR );
      HYPRE_Real const * const ptr_diag_data = hypre_CSRMatrixData( prt_diag_CSR );

      hypre_CSRMatrix const * const prt_offdiag_CSR = hypre_ParCSRMatrixOffd( m_parcsr_mat );
      HYPRE_Int const * const offdiag_IA = hypre_CSRMatrixI( prt_offdiag_CSR );
      HYPRE_Int const * const offdiag_JA = hypre_CSRMatrixJ( prt_offdiag_CSR );
      HYPRE_BigInt const * const col_map_offdiag = hypre_ParCSRMatrixColMapOffd( m_parcsr_mat );
      HYPRE_Real const * const ptr_offdiag_data = hypre_CSRMatrixData( prt_offdiag_CSR );

      for( HYPRE_Int i = 0; i < hypre_CSRMatrixNumRows( prt_diag_CSR ); ++i )
      {
        for( HYPRE_Int j = diag_IA[i]; j < diag_IA[i + 1]; ++j )
        {

          sprintf( str,
#ifdef GEOSX_USE_HYPRE_CUDA
                   "%i%20i%20i%24.10e\n",
#else
                   "%i%20lli%20lli%24.10e\n",
#endif
                   iRank,
                   firstRowID + LvArray::integerConversion< globalIndex >( i ),
                   firstDiagColID + LvArray::integerConversion< globalIndex >( diag_JA[j] ),
                   ptr_diag_data[j] );
          os << str;
        }
        for( HYPRE_Int j = offdiag_IA[i]; j < offdiag_IA[i + 1]; ++j )
        {
          sprintf( str,
#ifdef GEOSX_USE_HYPRE_CUDA
                   "%i%20i%20i%24.10e\n",
#else
                   "%i%20lli%20lli%24.10e\n",
#endif
                   iRank,
                   firstRowID + LvArray::integerConversion< globalIndex >( i ),
                   col_map_offdiag[ offdiag_JA[j] ],
                   ptr_offdiag_data[j] );
          os << str;
        }
      }
    }
  }
}

void HypreMatrix::write( string const & filename,
                         LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );

  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixPrintIJ( m_parcsr_mat, 1, 1, filename.c_str() ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      MPI_Comm const comm = getComm();
      int const rank = MpiWrapper::commRank( comm );

      if( numGlobalRows() * numGlobalCols() == 0 )
      {
        if( rank == 0 )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );
          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix coordinate real general\n" );
          hypre_fprintf( fp, "%lld %lld %d\n", numGlobalRows(), numGlobalCols(), 0 );
          std::fclose( fp );
        }
      }
      else
      {
        // Copy distributed parcsr matrix in a local CSR matrix on every process with at least one row
        // Warning: works for a parcsr matrix that is smaller than 2^31-1
        hypre_CSRMatrix * const fullMatrix = hypre_ParCSRMatrixToCSRMatrixAll( m_parcsr_mat );

        // Identify the smallest process where CSRmatrix exists
        int const printID = MpiWrapper::min( fullMatrix ? rank : MpiWrapper::commSize( comm ), comm );

        // Write to file CSRmatrix on one rank
        if( rank == printID )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );

          HYPRE_Int const num_rows = hypre_CSRMatrixNumRows( fullMatrix );
          HYPRE_Int const num_cols = hypre_CSRMatrixNumCols( fullMatrix );
          HYPRE_Int const num_nnz  = hypre_CSRMatrixNumNonzeros( fullMatrix );

          HYPRE_Real const * const matrix_data = hypre_CSRMatrixData( fullMatrix );
          HYPRE_Int const * const matrix_i     = hypre_CSRMatrixI( fullMatrix );
          HYPRE_Int const * const matrix_j     = hypre_CSRMatrixJ( fullMatrix );

          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix coordinate real general\n" );
          hypre_fprintf( fp, "%d %d %d\n", num_rows, num_cols, num_nnz );
          for( HYPRE_Int i = 0; i < num_rows; i++ )
          {
            for( HYPRE_Int j = matrix_i[i]; j < matrix_i[i+1]; j++ )
            {
              hypre_fprintf( fp, "%d %d %.16e\n", i + 1, matrix_j[j] + 1, matrix_data[j] );
            }
          }
          std::fclose( fp );
        }

        // Destroy CSRmatrix
        if( fullMatrix )
        {
          GEOSX_LAI_CHECK_ERROR( hypre_CSRMatrixDestroy( fullMatrix ) );
        }
      }
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported matrix output format" );
    }
  }
}

real64 HypreMatrix::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  HypreMatrix matT;
  transpose( matT );
  return matT.normInf();
}

real64 HypreMatrix::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );

  hypre_CSRMatrix const * const csr_diag = hypre_ParCSRMatrixDiag( m_parcsr_mat );
  hypre_CSRMatrix const * const csr_offdiag = hypre_ParCSRMatrixOffd( m_parcsr_mat );

  HYPRE_Int const * const ia_diag = hypre_CSRMatrixI( csr_diag );
  HYPRE_Int const * const ia_offdiag = hypre_CSRMatrixI( csr_offdiag );
  HYPRE_Real const * const va_diag = hypre_CSRMatrixData( csr_diag );
  HYPRE_Real const * const va_offdiag = hypre_CSRMatrixData( csr_offdiag );

  RAJA::ReduceMax< ReducePolicy< hypre::execPolicy >, HYPRE_Real > maxRowAbsSum( 0.0 );
  forAll< hypre::execPolicy >( numLocalRows(), [=] GEOSX_HYPRE_HOST_DEVICE ( localIndex const i )
  {
    HYPRE_Real rowAbsSum = 0.0;
    for( HYPRE_Int j = ia_diag[i]; j < ia_diag[i + 1]; ++j )
    {
      rowAbsSum += fabs( va_diag[j] );
    }
    for( HYPRE_Int j = ia_offdiag[i]; j < ia_offdiag[i + 1]; ++j )
    {
      rowAbsSum += fabs( va_offdiag[j] );
    }
    maxRowAbsSum.max( rowAbsSum );
  } );

  return MpiWrapper::max( maxRowAbsSum.get(), getComm() );

}

real64 HypreMatrix::normFrobenius() const
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_ParCSRMatrixFnorm( m_parcsr_mat );
}

void HypreMatrix::rightScale( HypreVector const & vec )
{
  GEOSX_UNUSED_VAR( vec )
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_ERROR( "Not implemented" );
}

void HypreMatrix::leftRightScale( HypreVector const & vecLeft,
                                  HypreVector const & vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

void HypreMatrix::transpose( HypreMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );

  // Transpose this->m_parcsr_mat
  HYPRE_ParCSRMatrix dst_parcsr;
  GEOSX_LAI_CHECK_ERROR( hypre_ParCSRMatrixTranspose( m_parcsr_mat, &dst_parcsr, 1 ) );

  // Create IJ layer (with matrix closed)
  dst.parCSRtoIJ( dst_parcsr );
}

MPI_Comm HypreMatrix::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJMatrixComm( m_ij_mat );
}

}// end namespace geosx
