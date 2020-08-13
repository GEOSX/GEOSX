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
 * @file EpetraVector.cpp
 */

#include "EpetraVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/trilinos/EpetraUtils.hpp"

#include <Epetra_FEVector.h>
#include <Epetra_Map.h>
#include <EpetraExt_MultiVectorOut.h>

#ifdef GEOSX_USE_MPI
  #include <Epetra_MpiComm.h>
#else
  #include <Epetra_SerialComm.h>
typedef Epetra_SerialComm Epetra_MpiComm;
#endif

namespace geosx
{
// Check matching requirements on index/value types between GEOSX and PETSc

static_assert( sizeof( long long ) == sizeof( globalIndex ),
               "long long and geosx::globalIndex must have the same size" );

static_assert(
  std::is_signed< long long >::value == std::is_signed< globalIndex >::value,
  "long long and geosx::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< double, real64 >::value,
               "double and geosx::real64 must be the same type" );

EpetraVector::EpetraVector() :
  VectorBase(),
  m_vector {}
{}

EpetraVector::EpetraVector( EpetraVector const & src ) :
  EpetraVector()
{
  *this = src;
}

EpetraVector::EpetraVector( EpetraVector && src ) noexcept :
  EpetraVector()
{
  *this = std::move( src );
}

EpetraVector &
EpetraVector::operator=( EpetraVector const & src )
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  if( m_vector )
  {
    *m_vector = *src.m_vector;
  }
  else
  {
    m_vector = std::make_unique< Epetra_FEVector >( *src.m_vector );
  }
  return *this;
}

EpetraVector &
EpetraVector::operator=( EpetraVector && src ) noexcept
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  m_vector = std::move( src.m_vector );
  return *this;
}

EpetraVector::~EpetraVector() = default;

bool
EpetraVector::created() const
{
  return bool( m_vector );
}

void
EpetraVector::createWithLocalSize( localIndex const localSize,
                                   MPI_Comm const & MPI_PARAM( comm ) )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );
  Epetra_Map const map( LvArray::integerConversion< long long >( -1 ),
                        LvArray::integerConversion< int >( localSize ),
                        0,
                        Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_vector = std::make_unique< Epetra_FEVector >( map );
}

void
EpetraVector::createWithGlobalSize( globalIndex const globalSize,
                                    MPI_Comm const & MPI_PARAM( comm ) )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );
  Epetra_Map const map( LvArray::integerConversion< long long >( globalSize ),
                        0,
                        Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_vector = std::make_unique< Epetra_FEVector >( map );
}

void
EpetraVector::create( arrayView1d< real64 const > const & localValues,
                      MPI_Comm const & MPI_PARAM( comm ) )
{
  GEOSX_LAI_ASSERT( closed() );

  localValues.move( LvArray::MemorySpace::CPU, false );

  int const localSize = LvArray::integerConversion< int >( localValues.size() );
  Epetra_Map const map( LvArray::integerConversion< long long >( -1 ),
                        localSize,
                        0,
                        Epetra_MpiComm( MPI_PARAM( comm ) ) );
  m_vector =
    std::make_unique< Epetra_FEVector >( Copy,
                                         map,
                                         const_cast< real64 * >( localValues.data() ),
                                         localSize,
                                         1 );
}

void
EpetraVector::set( globalIndex const globalRowIndex, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    m_vector->ReplaceGlobalValues( 1, toEpetraLongLong( &globalRowIndex ), &value ) );
}

void
EpetraVector::add( globalIndex const globalRowIndex, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    m_vector->SumIntoGlobalValues( 1, toEpetraLongLong( &globalRowIndex ), &value ) );
}

void
EpetraVector::set( globalIndex const * globalRowIndices,
                   real64 const * values,
                   localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    m_vector->ReplaceGlobalValues( LvArray::integerConversion< int >( size ),
                                   toEpetraLongLong( globalRowIndices ),
                                   values ) );
}

void
EpetraVector::add( globalIndex const * globalRowIndices,
                   real64 const * values,
                   localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    m_vector->SumIntoGlobalValues( LvArray::integerConversion< int >( size ),
                                   toEpetraLongLong( globalRowIndices ),
                                   values ) );
}

void
EpetraVector::set( arraySlice1d< globalIndex const > const & globalRowIndices,
                   arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( m_vector->ReplaceGlobalValues(
    LvArray::integerConversion< int >( values.size() ),
    toEpetraLongLong( globalRowIndices.dataIfContiguous() ),
    values.dataIfContiguous() ) );
}

void
EpetraVector::add( arraySlice1d< globalIndex const > const & globalRowIndices,
                   arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( m_vector->SumIntoGlobalValues(
    LvArray::integerConversion< int >( values.size() ),
    toEpetraLongLong( globalRowIndices.dataIfContiguous() ),
    values.dataIfContiguous() ) );
}

void
EpetraVector::set( real64 value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vector->PutScalar( value ) );
}

void
EpetraVector::zero()
{
  set( 0.0 );
}

void
EpetraVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vector->SetSeed( seed ) );
  GEOSX_LAI_CHECK_ERROR( m_vector->Random() );
}

void
EpetraVector::open()
{
  GEOSX_LAI_ASSERT( ready() );
  m_closed = false;
}

void
EpetraVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( m_vector->GlobalAssemble() );
  m_closed = true;
}

void
EpetraVector::reset()
{
  VectorBase::reset();
  m_vector.reset();
}

void
EpetraVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  GEOSX_LAI_CHECK_ERROR( m_vector->Scale( scalingFactor ) );
}

void
EpetraVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( m_vector->Reciprocal( *m_vector ) );
}

real64
EpetraVector::dot( EpetraVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  real64 tmp;
  GEOSX_LAI_CHECK_ERROR( m_vector->Dot( vec.unwrapped(), &tmp ) );
  return tmp;
}

void
EpetraVector::copy( EpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vector->Update( 1., x.unwrapped(), 0. ) );
}

void
EpetraVector::axpy( real64 const alpha, EpetraVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vector->Update( alpha, x.unwrapped(), 1. ) );
}

void
EpetraVector::axpby( real64 const alpha,
                     EpetraVector const & x,
                     real64 const beta )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( m_vector->Update( alpha, x.unwrapped(), beta ) );
}

real64
EpetraVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 tmp;
  m_vector->Norm1( &tmp );
  return tmp;
}

real64
EpetraVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 tmp;
  m_vector->Norm2( &tmp );
  return tmp;
}

real64
EpetraVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 tmp;
  m_vector->NormInf( &tmp );
  return tmp;
}

void
EpetraVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  m_vector->Print( os );
}

void
EpetraVector::write( string const & filename, LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::MATLAB_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR(
        EpetraExt::MultiVectorToMatlabFile( filename.c_str(), *m_vector ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      GEOSX_LAI_CHECK_ERROR(
        EpetraExt::MultiVectorToMatrixMarketFile( filename.c_str(), *m_vector ) );
      break;
    }
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }
}

real64
EpetraVector::get( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  return extractLocalVector()[getLocalRowID( globalRow )];
}

void
EpetraVector::get( arraySlice1d< globalIndex const > const & globalIndices,
                   arraySlice1d< real64 > const & values ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( values.size(), globalIndices.size() );
  GEOSX_LAI_ASSERT_GE(
    *std::min_element( globalIndices.dataIfContiguous(),
                       globalIndices.dataIfContiguous() + globalIndices.size() ),
    ilower() );
  GEOSX_LAI_ASSERT_GT(
    iupper(),
    *std::max_element( globalIndices.dataIfContiguous(),
                       globalIndices.dataIfContiguous() + globalIndices.size() ) );

  real64 const * const localVector = extractLocalVector();
  for( localIndex i = 0; i < globalIndices.size(); ++i )
  {
    values[i] = localVector[getLocalRowID( globalIndices[i] )];
  }
}

Epetra_FEVector const &
EpetraVector::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vector;
}

Epetra_FEVector &
EpetraVector::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return *m_vector;
}

globalIndex
EpetraVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->GlobalLength64();
}

localIndex
EpetraVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->MyLength();
}

localIndex
EpetraVector::getLocalRowID( globalIndex const globalRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->Map().LID( LvArray::integerConversion< long long >( globalRow ) );
}

globalIndex
EpetraVector::getGlobalRowID( localIndex const localRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( localRow, 0 );
  GEOSX_LAI_ASSERT_GT( localSize(), localRow );
  return m_vector->Map().GID64( LvArray::integerConversion< int >( localRow ) );
}

real64 const *
EpetraVector::extractLocalVector() const
{
  GEOSX_LAI_ASSERT( ready() );
  int dummy;
  double * localVector;
  GEOSX_LAI_CHECK_ERROR( m_vector->ExtractView( &localVector, &dummy ) );
  return localVector;
}

real64 *
EpetraVector::extractLocalVector()
{
  GEOSX_LAI_ASSERT( ready() );
  int dummy;
  double * localVector;
  GEOSX_LAI_CHECK_ERROR( m_vector->ExtractView( &localVector, &dummy ) );
  return localVector;
}

globalIndex
EpetraVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->Map().MinMyGID64();
}

globalIndex
EpetraVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vector->Map().MaxMyGID64() + 1;
}

MPI_Comm
EpetraVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
#ifdef GEOSX_USE_MPI
  return dynamic_cast< Epetra_MpiComm const & >( m_vector->Map().Comm() ).Comm();
#else
  return MPI_COMM_GEOSX;
#endif
}

}  // namespace geosx
