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
 * @file EpetraMatrix.cpp
 */

#include "EpetraMatrix.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_Map.h>
#include <Epetra_FECrsGraph.h>
#include <Epetra_FECrsMatrix.h>
#include <Epetra_FEVector.h>
#include <Epetra_Vector.h>
#include <EpetraExt_MatrixMatrix.h>
#include <EpetraExt_RowMatrixOut.h>
#include <EpetraExt_Transpose_RowMatrix.h>

#ifdef GEOSX_USE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
using Epetra_MpiComm = Epetra_SerialComm;
#endif

namespace geosx
{
EpetraMatrix::EpetraMatrix() :
  LinearOperator(),
  MatrixBase()
{}

EpetraMatrix::EpetraMatrix( EpetraMatrix const & src ) :
  EpetraMatrix()
{
  GEOSX_LAI_ASSERT( src.ready() );
  m_matrix = std::make_unique< Epetra_FECrsMatrix >( *src.m_matrix );
  m_src_map = std::make_unique< Epetra_Map >( m_matrix->DomainMap() );
  m_dst_map = std::make_unique< Epetra_Map >( m_matrix->RangeMap() );
  m_assembled = true;
}

EpetraMatrix::~EpetraMatrix() = default;

void
EpetraMatrix::createWithGlobalSize( globalIndex const globalRows,
                                    globalIndex const globalCols,
                                    localIndex const maxEntriesPerRow,
                                    MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalRows, 0 );
  GEOSX_LAI_ASSERT_GE( globalCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  m_dst_map =
    std::make_unique< Epetra_Map >( globalRows, 0, Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_src_map =
    std::make_unique< Epetra_Map >( globalCols, 0, Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_matrix = std::make_unique< Epetra_FECrsMatrix >(
    Copy,
    *m_dst_map,
    LvArray::integerConversion< int >( maxEntriesPerRow ),
    false );
}

void
EpetraMatrix::createWithLocalSize( localIndex const localRows,
                                   localIndex const localCols,
                                   localIndex const maxEntriesPerRow,
                                   MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localRows, 0 );
  GEOSX_LAI_ASSERT_GE( localCols, 0 );
  GEOSX_LAI_ASSERT_GE( maxEntriesPerRow, 0 );

  reset();

  m_dst_map =
    std::make_unique< Epetra_Map >( LvArray::integerConversion< globalIndex >( -1 ),
                                    LvArray::integerConversion< int >( localRows ),
                                    0,
                                    Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_src_map =
    std::make_unique< Epetra_Map >( LvArray::integerConversion< globalIndex >( -1 ),
                                    LvArray::integerConversion< int >( localCols ),
                                    0,
                                    Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_matrix = std::make_unique< Epetra_FECrsMatrix >(
    Copy,
    *m_dst_map,
    LvArray::integerConversion< int >( maxEntriesPerRow ),
    false );
}

bool
EpetraMatrix::created() const
{
  return bool( m_matrix );
}

void
EpetraMatrix::reset()
{
  MatrixBase::reset();
  m_matrix.reset();
  m_dst_map.reset();
  m_src_map.reset();
}

void
EpetraMatrix::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->PutScalar( value ) );
}

void
EpetraMatrix::zero()
{
  set( 0 );
}

void
EpetraMatrix::open()
{
  GEOSX_LAI_ASSERT( created() && closed() );
  m_closed = false;
}

void
EpetraMatrix::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->GlobalAssemble( *m_src_map, *m_dst_map ) );
  m_assembled = true;
  m_closed = true;
}

void
EpetraMatrix::add( globalIndex const rowIndex,
                   globalIndex const colIndex,
                   real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->SumIntoGlobalValues( rowIndex, 1, &value, toEpetraLongLong( &colIndex ) ) );
}

void
EpetraMatrix::set( globalIndex const rowIndex,
                   globalIndex const colIndex,
                   real64 const value )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->ReplaceGlobalValues( rowIndex, 1, &value, toEpetraLongLong( &colIndex ) ) );
}

void
EpetraMatrix::insert( globalIndex const rowIndex,
                      globalIndex const colIndex,
                      real64 const value )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG(
    m_matrix->InsertGlobalValues( rowIndex, 1, &value, toEpetraLongLong( &colIndex ) ) );
}

void
EpetraMatrix::add( globalIndex const rowIndex,
                   globalIndex const * colIndices,
                   real64 const * values,
                   localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->SumIntoGlobalValues( rowIndex,
                                   LvArray::integerConversion< int >( size ),
                                   values,
                                   toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::set( globalIndex const rowIndex,
                   globalIndex const * colIndices,
                   real64 const * values,
                   localIndex size )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->ReplaceGlobalValues( rowIndex,
                                   LvArray::integerConversion< int >( size ),
                                   values,
                                   toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::insert( globalIndex const rowIndex,
                      globalIndex const * colIndices,
                      real64 const * values,
                      localIndex size )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG(
    m_matrix->InsertGlobalValues( rowIndex,
                                  LvArray::integerConversion< int >( size ),
                                  values,
                                  toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::add( globalIndex const rowIndex,
                   arraySlice1d< globalIndex const > const & colIndices,
                   arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->SumIntoGlobalValues(
    rowIndex,
    LvArray::integerConversion< int >( colIndices.size() ),
    values,
    toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::set( globalIndex const rowIndex,
                   arraySlice1d< globalIndex const > const & colIndices,
                   arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->ReplaceGlobalValues(
    rowIndex,
    LvArray::integerConversion< int >( colIndices.size() ),
    values,
    toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::insert( globalIndex const rowIndex,
                      arraySlice1d< globalIndex const > const & colIndices,
                      arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG( m_matrix->InsertGlobalValues(
    rowIndex,
    LvArray::integerConversion< int >( colIndices.size() ),
    values,
    toEpetraLongLong( colIndices ) ) );
}

void
EpetraMatrix::add(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->SumIntoGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::set(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->ReplaceGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::insert(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::ROW_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG( m_matrix->InsertGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::add(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->SumIntoGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::COLUMN_MAJOR ) );
}

void
EpetraMatrix::set(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->ReplaceGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::COLUMN_MAJOR ) );
}

void
EpetraMatrix::insert(
  arraySlice1d< globalIndex const > const & rowIndices,
  arraySlice1d< globalIndex const > const & colIndices,
  arraySlice2d< real64 const, MatrixLayout::COL_MAJOR > const & values )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG( m_matrix->InsertGlobalValues(
    LvArray::integerConversion< int >( rowIndices.size() ),
    toEpetraLongLong( rowIndices ),
    LvArray::integerConversion< int >( colIndices.size() ),
    toEpetraLongLong( colIndices ),
    values.dataIfContiguous(),
    Epetra_FECrsMatrix::COLUMN_MAJOR ) );
}

void
EpetraMatrix::add( globalIndex const * rowIndices,
                   globalIndex const * colIndices,
                   real64 const * values,
                   localIndex const numRows,
                   localIndex const numCols )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->SumIntoGlobalValues( LvArray::integerConversion< int >( numRows ),
                                   toEpetraLongLong( rowIndices ),
                                   LvArray::integerConversion< int >( numCols ),
                                   toEpetraLongLong( colIndices ),
                                   values,
                                   Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::set( globalIndex const * rowIndices,
                   globalIndex const * colIndices,
                   real64 const * values,
                   localIndex const numRows,
                   localIndex const numCols )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->ReplaceGlobalValues( LvArray::integerConversion< int >( numRows ),
                                   toEpetraLongLong( rowIndices ),
                                   LvArray::integerConversion< int >( numCols ),
                                   toEpetraLongLong( colIndices ),
                                   values,
                                   Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::insert( globalIndex const * rowIndices,
                      globalIndex const * colIndices,
                      real64 const * values,
                      localIndex const numRows,
                      localIndex const numCols )
{
  GEOSX_LAI_ASSERT( insertable() );
  GEOSX_LAI_CHECK_ERROR_NNEG(
    m_matrix->InsertGlobalValues( LvArray::integerConversion< int >( numRows ),
                                  toEpetraLongLong( rowIndices ),
                                  LvArray::integerConversion< int >( numCols ),
                                  toEpetraLongLong( colIndices ),
                                  values,
                                  Epetra_FECrsMatrix::ROW_MAJOR ) );
}

void
EpetraMatrix::apply( EpetraVector const & src, EpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR(
    m_matrix->Multiply( false, src.unwrapped(), dst.unwrapped() ) );
}

void
EpetraMatrix::applyTranspose( EpetraVector const & src, EpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), dst.globalSize() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), src.globalSize() );

  GEOSX_LAI_CHECK_ERROR(
    m_matrix->Multiply( true, src.unwrapped(), dst.unwrapped() ) );
}

void
EpetraMatrix::multiply( EpetraMatrix const & src, EpetraMatrix & dst ) const
{
  this->multiply( false, src, false, dst );
}

void
EpetraMatrix::leftMultiplyTranspose( EpetraMatrix const & src,
                                     EpetraMatrix & dst ) const
{
  this->multiply( true, src, false, dst );
}

void
EpetraMatrix::rightMultiplyTranspose( EpetraMatrix const & src,
                                      EpetraMatrix & dst ) const
{
  src.multiply( false, *this, true, dst );
}

void
EpetraMatrix::create( Epetra_CrsMatrix const & src )
{
  GEOSX_LAI_ASSERT( closed() );
  reset();

  m_matrix = std::make_unique< Epetra_FECrsMatrix >( Copy, src.Graph() );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->FillComplete( src.DomainMap(), src.RangeMap(), true ) );
  m_src_map = std::make_unique< Epetra_Map >( m_matrix->DomainMap() );
  m_dst_map = std::make_unique< Epetra_Map >( m_matrix->RangeMap() );

  // We have to copy the data using the "expert" interface,
  // since there is no direct way to convert CrsMatrix to FECrsMatrix.
  // This works, because we initialized dst.m_matrix with the graph (sparsity) of result

  GEOSX_LAI_ASSERT_EQ( src.NumMyRows(), m_matrix->NumMyRows() );
  GEOSX_LAI_ASSERT_EQ( src.NumMyNonzeros(), m_matrix->NumMyNonzeros() );

  int * offsets_src;
  int * indices_src;
  double * values_src;
  GEOSX_LAI_CHECK_ERROR(
    src.ExtractCrsDataPointers( offsets_src, indices_src, values_src ) );

  int * offsets_dst;
  int * indices_dst;
  double * values_dst;
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->ExtractCrsDataPointers( offsets_dst, indices_dst, values_dst ) );

  std::copy( offsets_src, offsets_src + src.NumMyRows() + 1, offsets_dst );
  std::copy( indices_src, indices_src + src.NumMyNonzeros(), indices_dst );
  std::copy( values_src, values_src + src.NumMyNonzeros(), values_dst );

  m_assembled = true;
}

void
EpetraMatrix::multiplyRAP( EpetraMatrix const & R,
                           EpetraMatrix const & P,
                           EpetraMatrix & dst ) const
{
  // TODO: ML_Epetra_RAP does not work with long long indices, find a workaround?
#if 0
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( R.ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), R.numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  Epetra_CrsMatrix * result = nullptr;
  GEOSX_LAI_CHECK_ERROR( ML_Epetra::ML_Epetra_RAP( unwrapped(), P.unwrapped(), R.unwrapped(), result, false ) );

  dst.create( *result );
  delete result;
#else
  MatrixBase::multiplyRAP( R, P, dst );
#endif
}

void
EpetraMatrix::multiplyPtAP( EpetraMatrix const & P, EpetraMatrix & dst ) const
{
  // TODO: ML_Epetra_PtAP does not work with long long indices, find a workaround?
#if 0
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( P.ready() );
  GEOSX_LAI_ASSERT_EQ( numGlobalRows(), P.numGlobalRows() );
  GEOSX_LAI_ASSERT_EQ( numGlobalCols(), P.numGlobalRows() );

  Epetra_CrsMatrix * result = nullptr;
  GEOSX_LAI_CHECK_ERROR( ML_Epetra::ML_Epetra_PtAP( unwrapped(), P.unwrapped(), result, false ) );

  dst.create( *result );
  delete result;
#else
  MatrixBase::multiplyPtAP( P, dst );
#endif
}

void
EpetraMatrix::gemv( real64 const alpha,
                    EpetraVector const & x,
                    real64 const beta,
                    EpetraVector & y,
                    bool useTranspose ) const
{
  GEOSX_LAI_ASSERT( ready() );
  EpetraVector Ax( y );
  GEOSX_LAI_CHECK_ERROR(
    m_matrix->Multiply( useTranspose, x.unwrapped(), Ax.unwrapped() ) );
  y.axpby( alpha, Ax, beta );
}

void
EpetraMatrix::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  GEOSX_LAI_CHECK_ERROR( m_matrix->Scale( scalingFactor ) );
}

void
EpetraMatrix::leftScale( EpetraVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->LeftScale( *vec.unwrapped()( 0 ) ) );
}

void
EpetraMatrix::rightScale( EpetraVector const & vec )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_matrix->RightScale( *vec.unwrapped()( 0 ) ) );
}

void
EpetraMatrix::leftRightScale( EpetraVector const & vecLeft,
                              EpetraVector const & vecRight )
{
  leftScale( vecLeft );
  rightScale( vecRight );
}

void
EpetraMatrix::transpose( EpetraMatrix & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );

  // XXX: Transposer always creates a CrsMatrix, but doesn't expose this fact in the API:
  // https://docs.trilinos.org/dev/packages/epetraext/doc/html/EpetraExt__Transpose__RowMatrix_8cpp_source.html#l00248
  EpetraExt::RowMatrix_Transpose transposer( m_src_map.get() );
  Epetra_CrsMatrix & trans =
    dynamic_cast< Epetra_CrsMatrix & >( transposer( *m_matrix ) );
  dst.create( trans );
}

real64
EpetraMatrix::clearRow( globalIndex const globalRow,
                        bool const keepDiag,
                        real64 const diagValue )
{
  GEOSX_LAI_ASSERT( modifiable() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  int length;
  int * colIndices;
  double * values;
  GEOSX_LAI_CHECK_ERROR( m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ),
                                                     length,
                                                     values,
                                                     colIndices ) );

  bool const square = numGlobalRows() == numGlobalCols();

  real64 oldDiag = 0.0;
  for( int j = 0; j < length; ++j )
  {
    if( square && m_matrix->GCID64( colIndices[j] ) == globalRow )
    {
      oldDiag = values[j];
    }
    values[j] = 0.0;
  }

  // Set diagonal value
  real64 const newDiag = keepDiag ? oldDiag : diagValue;
  if( square && std::fabs( newDiag ) > 0.0 )
  {
    set( globalRow, globalRow, newDiag );
  }
  return oldDiag;
}

void
EpetraMatrix::addEntries( EpetraMatrix const & src, real64 const scale )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == src.numGlobalRows() );
  GEOSX_LAI_ASSERT( numGlobalCols() == src.numGlobalCols() );

  GEOSX_LAI_CHECK_ERROR(
    EpetraExt::MatrixMatrix::Add( src.unwrapped(), false, scale, *m_matrix, 1.0 ) );
}

void
EpetraMatrix::addDiagonal( EpetraVector const & src )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( numGlobalRows() == numGlobalCols() );
  GEOSX_LAI_ASSERT( numLocalRows() == src.localSize() );

  // XXX: It's not clear if this is better or worse than element-wise add()
  Epetra_Vector diag( *m_dst_map, false );
  GEOSX_LAI_CHECK_ERROR( m_matrix->ExtractDiagonalCopy( diag ) );
  GEOSX_LAI_CHECK_ERROR( diag.Update( 1.0, src.unwrapped(), 1.0 ) );
  GEOSX_LAI_CHECK_ERROR( m_matrix->ReplaceDiagonalValues( diag ) );
}

localIndex
EpetraMatrix::maxRowLength() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->MaxNumEntries();
}

localIndex
EpetraMatrix::localRowLength( localIndex localRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->NumMyEntries( LvArray::integerConversion< int >( localRowIndex ) );
}

localIndex
EpetraMatrix::globalRowLength( globalIndex globalRowIndex ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->NumGlobalEntries( globalRowIndex );
}

void
EpetraMatrix::getRowCopy( globalIndex globalRow,
                          arraySlice1d< globalIndex > const & colIndices,
                          arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  int numEntries;
  int * indices_ptr;
  double * values_ptr;

  GEOSX_LAI_CHECK_ERROR( m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ),
                                                     numEntries,
                                                     values_ptr,
                                                     indices_ptr ) );

  GEOSX_LAI_ASSERT_GE( colIndices.size(), numEntries );
  GEOSX_LAI_ASSERT_GE( values.size(), numEntries );

  localIndex const length = LvArray::integerConversion< localIndex >( numEntries );
  for( localIndex i = 0; i < length; ++i )
  {
    colIndices[i] =
      LvArray::integerConversion< globalIndex >( m_matrix->GCID64( indices_ptr[i] ) );
    values[i] = values_ptr[i];
  }
}

real64
EpetraMatrix::getDiagValue( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( assembled() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  int length;
  int * indices_ptr;
  double * values_ptr;

  GEOSX_LAI_CHECK_ERROR( m_matrix->ExtractMyRowView( m_matrix->LRID( globalRow ),
                                                     length,
                                                     values_ptr,
                                                     indices_ptr ) );

  for( int j = 0; j < length; ++j )
  {
    if( m_matrix->GCID64( indices_ptr[j] ) == globalRow )
    {
      return values_ptr[j];
    }
  }

  return 0.0;
}

void
EpetraMatrix::extractDiagonal( EpetraVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  // This doesn't work because ExtractDiagonalCopy takes an Epetra_Vector,
  // not Epetra_FEVector which is unrelated (ugh):
  // GEOSX_LAI_CHECK_ERROR( m_matrix->ExtractDiagonalCopy( dst.unwrapped() ) );

  dst.zero();
  real64 * const values = dst.extractLocalVector();

  int length;
  int * indices_ptr;
  double * values_ptr;

  for( localIndex localRow = 0; localRow < numLocalRows(); ++localRow )
  {
    GEOSX_LAI_CHECK_ERROR(
      m_matrix->ExtractMyRowView( localRow, length, values_ptr, indices_ptr ) );
    for( int j = 0; j < length; ++j )
    {
      if( indices_ptr[j] == localRow )
      {
        values[localRow] = values_ptr[j];
      }
    }
  }
}

Epetra_FECrsMatrix const &
EpetraMatrix::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return *m_matrix;
}

Epetra_FECrsMatrix &
EpetraMatrix::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return *m_matrix;
}

globalIndex
EpetraMatrix::numGlobalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->NumGlobalElements64();
}

globalIndex
EpetraMatrix::numGlobalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->NumGlobalElements64();
}

globalIndex
EpetraMatrix::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->MinMyGID64();
}

globalIndex
EpetraMatrix::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_dst_map->MaxMyGID64() + 1;
}

globalIndex
EpetraMatrix::jlower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->MinMyGID64();
}

globalIndex
EpetraMatrix::jupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->MaxMyGID64() + 1;
}

localIndex
EpetraMatrix::numLocalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->NumMyNonzeros();
}

globalIndex
EpetraMatrix::numGlobalNonzeros() const
{
  GEOSX_LAI_ASSERT( assembled() );
  return m_matrix->NumGlobalNonzeros64();
}

real64
EpetraMatrix::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_matrix->NormInf();
}

real64
EpetraMatrix::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_matrix->NormOne();
}

real64
EpetraMatrix::normFrobenius() const
{
  GEOSX_LAI_ASSERT( ready() );
  return m_matrix->NormFrobenius();
}

localIndex
EpetraMatrix::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  return m_matrix->LRID( index );
}

globalIndex
EpetraMatrix::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( index, 0 );
  GEOSX_LAI_ASSERT_GT( numLocalRows(), index );
  return m_matrix->GRID64( LvArray::integerConversion< int >( index ) );
}

localIndex
EpetraMatrix::numLocalCols() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_src_map->NumMyElements();
}

localIndex
EpetraMatrix::numLocalRows() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_matrix->RowMap().NumMyElements();
}

MPI_Comm
EpetraMatrix::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
#ifdef GEOSX_USE_MPI
  return dynamic_cast< Epetra_MpiComm const & >( m_matrix->RowMap().Comm() ).Comm();
#else
  return MPI_COMM_GEOSX;
#endif
}

void
EpetraMatrix::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  m_matrix->Print( os );
}

void
EpetraMatrix::write( string const & filename, LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      std::ofstream ofs( filename );
      print( ofs );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      GEOSX_LAI_CHECK_ERROR(
        EpetraExt::RowMatrixToMatrixMarketFile( filename.c_str(), *m_matrix ) );
      break;
    }
    case LAIOutputFormat::MATLAB_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR(
        EpetraExt::RowMatrixToMatlabFile( filename.c_str(), *m_matrix ) );
      break;
    }
    default:
      GEOSX_ERROR( "Unsupported matrix output format" );
  }
}

void
EpetraMatrix::multiply( bool const transA,
                        EpetraMatrix const & B,
                        bool const transB,
                        EpetraMatrix & C ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( B.ready() );

  C.createWithLocalSize( transA ? numLocalCols() : numLocalRows(),
                         transB ? B.numLocalRows() : B.numLocalCols(),
                         1,  // TODO: estimate entries per row?
                         getComm() );

  GEOSX_LAI_CHECK_ERROR( EpetraExt::MatrixMatrix::Multiply( unwrapped(),
                                                            transA,
                                                            B.unwrapped(),
                                                            transB,
                                                            C.unwrapped(),
                                                            true ) );
  C.m_assembled = true;
}

}  // namespace geosx
