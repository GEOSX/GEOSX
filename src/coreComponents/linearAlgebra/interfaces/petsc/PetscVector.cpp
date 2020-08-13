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
 * @file PetscVector.cpp
 */

#include "PetscVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/petsc/PetscUtils.hpp"

#include <petscvec.h>

namespace geosx
{
// Check matching requirements on index/value types between GEOSX and PETSc

static_assert( sizeof( PetscInt ) == sizeof( globalIndex ),
               "PetscInt and geosx::globalIndex must have the same size" );

static_assert(
  std::is_signed< PetscInt >::value == std::is_signed< globalIndex >::value,
  "PetscInt and geoex::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< PetscScalar, real64 >::value,
               "PetscScalar and geosx::real64 must be the same type" );

PetscVector::PetscVector() :
  VectorBase(),
  m_vec {}
{}

PetscVector::PetscVector( PetscVector const & src ) :
  PetscVector()
{
  *this = src;
}

PetscVector::PetscVector( PetscVector && src ) noexcept :
  PetscVector()
{
  *this = std::move( src );
}

PetscVector &
PetscVector::operator=( PetscVector const & src )
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  reset();
  GEOSX_LAI_CHECK_ERROR( VecDuplicate( src.m_vec, &m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecCopy( src.m_vec, m_vec ) );
  return *this;
}

PetscVector &
PetscVector::operator=( PetscVector && src ) noexcept
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  std::swap( m_vec, src.m_vec );
  return *this;
}

PetscVector::~PetscVector()
{
  reset();
}

void
PetscVector::reset()
{
  VectorBase::reset();
  if( m_vec != nullptr )
  {
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &m_vec ) );
    m_vec = nullptr;
  }
}

void
PetscVector::createWithLocalSize( localIndex const localSize,
                                  MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );
  reset();
  GEOSX_LAI_CHECK_ERROR( VecCreate( comm, &m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecSetType( m_vec, VECMPI ) );
  GEOSX_LAI_CHECK_ERROR( VecSetSizes( m_vec, localSize, PETSC_DETERMINE ) );
}

void
PetscVector::createWithGlobalSize( globalIndex const globalSize,
                                   MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );
  reset();
  GEOSX_LAI_CHECK_ERROR( VecCreate( comm, &m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecSetType( m_vec, VECMPI ) );
  GEOSX_LAI_CHECK_ERROR( VecSetSizes( m_vec, PETSC_DECIDE, globalSize ) );
}

void
PetscVector::create( arrayView1d< real64 const > const & localValues,
                     MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  reset();
  PetscInt const size = localValues.size();
  PetscScalar * values;

  localValues.move( LvArray::MemorySpace::CPU, false );

  GEOSX_LAI_CHECK_ERROR( VecCreate( comm, &m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecSetType( m_vec, VECMPI ) );
  GEOSX_LAI_CHECK_ERROR( VecSetSizes( m_vec, size, PETSC_DETERMINE ) );
  GEOSX_LAI_CHECK_ERROR( VecGetArray( m_vec, &values ) );

  // set vector values
  for( int i = 0; i < size; i++ )
  {
    values[i] = localValues[i];
  }

  GEOSX_LAI_CHECK_ERROR( VecRestoreArray( m_vec, &values ) );
}

bool
PetscVector::created() const
{
  return m_vec != nullptr;
}

void
PetscVector::set( globalIndex const globalRow, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( VecSetValue( m_vec, globalRow, value, INSERT_VALUES ) );
}

void
PetscVector::add( globalIndex const globalRow, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( VecSetValue( m_vec, globalRow, value, ADD_VALUES ) );
}

void
PetscVector::set( globalIndex const * globalIndices,
                  real64 const * values,
                  localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    VecSetValues( m_vec, size, toPetscInt( globalIndices ), values, INSERT_VALUES ) );
}

void
PetscVector::add( globalIndex const * globalIndices,
                  real64 const * values,
                  localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    VecSetValues( m_vec, size, toPetscInt( globalIndices ), values, ADD_VALUES ) );
}

void
PetscVector::set( arraySlice1d< globalIndex const > const & globalIndices,
                  arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( VecSetValues( m_vec,
                                       values.size(),
                                       toPetscInt( globalIndices.dataIfContiguous() ),
                                       values.dataIfContiguous(),
                                       INSERT_VALUES ) );
}

void
PetscVector::add( arraySlice1d< globalIndex const > const & globalIndices,
                  arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( VecSetValues( m_vec,
                                       values.size(),
                                       toPetscInt( globalIndices.dataIfContiguous() ),
                                       values.dataIfContiguous(),
                                       ADD_VALUES ) );
}

void
PetscVector::set( real64 const value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( VecSet( m_vec, value ) );
}

void
PetscVector::zero()
{
  set( 0.0 );
}

void
PetscVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );

  // create random context
  PetscRandom ran;
  GEOSX_LAI_CHECK_ERROR( PetscRandomCreate( PETSC_COMM_WORLD, &ran ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomSetInterval( ran, -1.0, 1.0 ) );

  // set random seed
  GEOSX_LAI_CHECK_ERROR( PetscRandomSetSeed( ran, seed ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomSeed( ran ) );

  // create random vector
  GEOSX_LAI_CHECK_ERROR( VecSetRandom( m_vec, ran ) );
  GEOSX_LAI_CHECK_ERROR( PetscRandomDestroy( &ran ) );
}

void
PetscVector::open()
{
  GEOSX_LAI_ASSERT( ready() );
  m_closed = false;
}

void
PetscVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  // assemble the vector after setting values
  GEOSX_LAI_CHECK_ERROR( VecAssemblyBegin( m_vec ) );
  GEOSX_LAI_CHECK_ERROR( VecAssemblyEnd( m_vec ) );
  m_closed = true;
}

void
PetscVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  GEOSX_LAI_CHECK_ERROR( VecScale( m_vec, scalingFactor ) );
}

void
PetscVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( VecReciprocal( m_vec ) );
}

real64
PetscVector::dot( PetscVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  real64 dot;
  GEOSX_LAI_CHECK_ERROR( VecDot( m_vec, vec.m_vec, &dot ) );
  return dot;
}

void
PetscVector::copy( PetscVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( VecSet( m_vec, 0 ) );
  GEOSX_LAI_CHECK_ERROR( VecAXPY( m_vec, 1.0, x.m_vec ) );
}

void
PetscVector::axpy( real64 const alpha, PetscVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( VecAXPY( m_vec, alpha, x.m_vec ) );
}

void
PetscVector::axpby( real64 const alpha, PetscVector const & x, real64 const beta )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( VecScale( m_vec, beta ) );
  GEOSX_LAI_CHECK_ERROR( VecAXPY( m_vec, alpha, x.m_vec ) );
}

real64
PetscVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_1, &result ) );
  return result;
}

real64
PetscVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_2, &result ) );
  return result;
}

real64
PetscVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );
  real64 result;
  GEOSX_LAI_CHECK_ERROR( VecNorm( m_vec, NORM_INFINITY, &result ) );
  return result;
}

void
PetscVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );
  if( &os == &std::cout )
  {
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, PETSC_VIEWER_STDOUT_WORLD ) );
  }
  else if( &os == &std::cerr )
  {
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, PETSC_VIEWER_STDERR_WORLD ) );
  }
  else
  {
    GEOSX_ERROR( "Output to a generic stream not implemented" );
  }
}

void
PetscVector::write( string const & filename, LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  PetscViewerFormat petscFormat = PETSC_VIEWER_DEFAULT;

  bool useMatrixMarket = false;
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_COMMON;
    }
    break;
    case LAIOutputFormat::MATLAB_ASCII:
    {
      petscFormat = PETSC_VIEWER_ASCII_MATLAB;
    }
    break;
    case LAIOutputFormat::MATRIX_MARKET:
    {
      useMatrixMarket = true;
    }
    break;
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }

  if( !useMatrixMarket )
  {
    PetscViewer viewer;
    GEOSX_LAI_CHECK_ERROR(
      PetscViewerASCIIOpen( getComm(), filename.c_str(), &viewer ) );
    GEOSX_LAI_CHECK_ERROR( PetscViewerPushFormat( viewer, petscFormat ) );
    GEOSX_LAI_CHECK_ERROR( VecView( m_vec, viewer ) );
    GEOSX_LAI_CHECK_ERROR( PetscViewerDestroy( &viewer ) );
  }
  else
  {
    VecScatter scatter;
    Vec globalVec;
    GEOSX_LAI_CHECK_ERROR( VecScatterCreateToAll( m_vec, &scatter, &globalVec ) );
    GEOSX_LAI_CHECK_ERROR(
      VecScatterBegin( scatter, m_vec, globalVec, INSERT_VALUES, SCATTER_FORWARD ) );
    GEOSX_LAI_CHECK_ERROR(
      VecScatterEnd( scatter, m_vec, globalVec, INSERT_VALUES, SCATTER_FORWARD ) );
    if( MpiWrapper::Comm_rank( getComm() ) == 0 )
    {
      PetscScalar * v;
      GEOSX_LAI_CHECK_ERROR( VecGetArray( globalVec, &v ) );

      FILE * fp = std::fopen( filename.c_str(), "w" );
      fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
      fprintf( fp, "%lld %d\n", globalSize(), 1 );
      for( globalIndex i = 0; i < globalSize(); i++ )
      {
        fprintf( fp, "%.16e\n", v[i] );
      }
      std::fclose( fp );

      GEOSX_LAI_CHECK_ERROR( VecRestoreArray( globalVec, &v ) );
    }
    GEOSX_LAI_CHECK_ERROR( VecScatterDestroy( &scatter ) );
    GEOSX_LAI_CHECK_ERROR( VecDestroy( &globalVec ) );
  }
}

real64
PetscVector::get( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  real64 value;
  GEOSX_LAI_CHECK_ERROR( VecGetValues( m_vec, 1, toPetscInt( &globalRow ), &value ) );
  return value;
}

void
PetscVector::get( arraySlice1d< globalIndex const > const & globalIndices,
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

  GEOSX_LAI_CHECK_ERROR( VecGetValues( m_vec,
                                       globalIndices.size(),
                                       toPetscInt( globalIndices.dataIfContiguous() ),
                                       values.dataIfContiguous() ) );
}

Vec const &
PetscVector::unwrapped() const
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec;
}

Vec &
PetscVector::unwrapped()
{
  GEOSX_LAI_ASSERT( created() );
  return m_vec;
}

MPI_Comm
PetscVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  MPI_Comm comm;
  GEOSX_LAI_CHECK_ERROR(
    PetscObjectGetComm( reinterpret_cast< PetscObject >( m_vec ), &comm ) );
  return comm;
}

globalIndex
PetscVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt size;
  GEOSX_LAI_CHECK_ERROR( VecGetSize( m_vec, &size ) );
  return size;
}

localIndex
PetscVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt size;
  GEOSX_LAI_CHECK_ERROR( VecGetLocalSize( m_vec, &size ) );
  return size;
}

localIndex
PetscVector::getLocalRowID( globalIndex const globalRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( VecGetOwnershipRange( m_vec, &low, &high ) );
  return ( globalRow >= low && globalRow < high )
    ? LvArray::integerConversion< localIndex >( globalRow - low )
    : -1;
}

globalIndex
PetscVector::getGlobalRowID( localIndex const localRow ) const
{
  GEOSX_LAI_ASSERT( created() );
  GEOSX_LAI_ASSERT_GE( localRow, 0 );
  GEOSX_LAI_ASSERT_GT( localSize(), localRow );
  return ilower() + localRow;
}

real64 const *
PetscVector::extractLocalVector() const
{
  GEOSX_LAI_ASSERT( ready() );
  PetscScalar * avec;
  GEOSX_LAI_CHECK_ERROR( VecGetArray( m_vec, &avec ) );
  real64 const * const localVector = avec;
  GEOSX_LAI_CHECK_ERROR( VecRestoreArray( m_vec, &avec ) );
  return localVector;
}

real64 *
PetscVector::extractLocalVector()
{
  GEOSX_LAI_ASSERT( ready() );
  PetscScalar * avec;
  GEOSX_LAI_CHECK_ERROR( VecGetArray( m_vec, &avec ) );
  real64 * const localVector = avec;
  GEOSX_LAI_CHECK_ERROR( VecRestoreArray( m_vec, &avec ) );
  return localVector;
}

globalIndex
PetscVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( VecGetOwnershipRange( m_vec, &low, &high ) );
  return low;
}

globalIndex
PetscVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  PetscInt low, high;
  GEOSX_LAI_CHECK_ERROR( VecGetOwnershipRange( m_vec, &low, &high ) );
  return high;
}

}  // namespace geosx
