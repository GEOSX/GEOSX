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
 * @file HypreVector.cpp
 */

#include "HypreVector.hpp"

#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"

#include "HYPRE.h"
#include "_hypre_IJ_mv.h"
#include "_hypre_parcsr_mv.h"

#include <iomanip>

namespace geosx
{
// Check matching requirements on index/value types between GEOSX and Hypre

static_assert( sizeof( HYPRE_BigInt ) == sizeof( globalIndex ),
               "HYPRE_BigInt and geosx::globalIndex must have the same size" );

static_assert(
  std::is_signed< HYPRE_BigInt >::value == std::is_signed< globalIndex >::value,
  "HYPRE_BigInt and geoex::globalIndex must both be signed or unsigned" );

static_assert( std::is_same< HYPRE_Real, real64 >::value,
               "HYPRE_Real and geosx::real64 must be the same type" );

// Helper function that performs the following sequence of IJVEctor
// call: Create, SetObjectType, Initialize.
static void
initialize( MPI_Comm const & comm,
            HYPRE_BigInt const & jlower,
            HYPRE_BigInt const & jupper,
            HYPRE_IJVector & ij_vector )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorCreate( comm, jlower, jupper, &ij_vector ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetObjectType( ij_vector, HYPRE_PARCSR ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorInitialize( ij_vector ) );
}

// Helper function that performs the following sequence of IJVEctor
// call: Assemble, GetObject.
static void
finalize( HYPRE_IJVector & ij_vector, HYPRE_ParVector & par_vector )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( ij_vector ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetObject( ij_vector, (void **)&par_vector ) );
}

HypreVector::HypreVector() :
  VectorBase(),
  m_ij_vector {},
  m_par_vector {}
{}

// Copy constructor
HypreVector::HypreVector( HypreVector const & src ) :
  HypreVector()
{
  *this = src;
}

// Move constructor
HypreVector::HypreVector( HypreVector && src ) noexcept :
  HypreVector()
{
  *this = std::move( src );
}

HypreVector &
HypreVector::operator=( HypreVector const & src )
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );

  reset();

  HYPRE_BigInt jlower, jupper;
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorGetLocalRange( src.m_ij_vector, &jlower, &jupper ) );

  initialize( src.getComm(), jlower, jupper, m_ij_vector );
  finalize( m_ij_vector, m_par_vector );
  copy( src );

  return *this;
}

HypreVector &
HypreVector::operator=( HypreVector && src ) noexcept
{
  GEOSX_LAI_ASSERT( &src != this );
  GEOSX_LAI_ASSERT( src.ready() );
  std::swap( m_ij_vector, src.m_ij_vector );
  std::swap( m_par_vector, src.m_par_vector );
  return *this;
}

void
HypreVector::reset()
{
  VectorBase::reset();
  if( m_ij_vector )
  {
    HYPRE_IJVectorDestroy( m_ij_vector );
    m_ij_vector = nullptr;
  }
}

HypreVector::~HypreVector()
{
  reset();
}

void
HypreVector::createWithLocalSize( localIndex const localSize,
                                  MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( localSize, 0 );

  HYPRE_BigInt const jlower = MpiWrapper::PrefixSum< HYPRE_BigInt >(
    LvArray::integerConversion< HYPRE_BigInt >( localSize ) );
  HYPRE_BigInt const jupper =
    jlower + LvArray::integerConversion< HYPRE_BigInt >( localSize ) - 1;

  initialize( comm, jlower, jupper, m_ij_vector );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
  finalize( m_ij_vector, m_par_vector );
}

void
HypreVector::createWithGlobalSize( globalIndex const globalSize,
                                   MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );
  GEOSX_LAI_ASSERT_GE( globalSize, 0 );

  HYPRE_Int const rank =
    LvArray::integerConversion< HYPRE_Int >( MpiWrapper::Comm_rank( comm ) );
  HYPRE_Int const nproc =
    LvArray::integerConversion< HYPRE_Int >( MpiWrapper::Comm_size( comm ) );

  HYPRE_Int const localSize =
    LvArray::integerConversion< HYPRE_Int >( globalSize / nproc );
  HYPRE_Int const residual =
    LvArray::integerConversion< HYPRE_Int >( globalSize % nproc );

  HYPRE_BigInt const ilower = LvArray::integerConversion< HYPRE_BigInt >(
    rank * localSize + ( rank == 0 ? 0 : residual ) );
  HYPRE_BigInt const iupper = LvArray::integerConversion< HYPRE_BigInt >(
    ilower + localSize + ( rank == 0 ? residual : 0 ) - 1 );

  initialize( comm, ilower, iupper, m_ij_vector );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
  finalize( m_ij_vector, m_par_vector );
}

void
HypreVector::create( arrayView1d< real64 const > const & localValues,
                     MPI_Comm const & comm )
{
  GEOSX_LAI_ASSERT( closed() );

  localValues.move( LvArray::MemorySpace::CPU, false );

  HYPRE_BigInt const localSize =
    LvArray::integerConversion< HYPRE_BigInt >( localValues.size() );

  HYPRE_BigInt const jlower = MpiWrapper::PrefixSum< HYPRE_BigInt >( localSize );
  HYPRE_BigInt const jupper = jlower + localSize - 1;

  initialize( comm, jlower, jupper, m_ij_vector );
  finalize( m_ij_vector, m_par_vector );

  HYPRE_Real * const local_data = extractLocalVector();
  for( localIndex i = 0; i < localValues.size(); ++i )
  {
    local_data[i] = localValues[i];
  }
}

bool
HypreVector::created() const
{
  return m_ij_vector != nullptr && m_par_vector != nullptr;
}

void
HypreVector::set( globalIndex const globalRow, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorSetValues( m_ij_vector, 1, &globalRow, &value ) );
}

void
HypreVector::add( globalIndex const globalRow, real64 const value )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorAddToValues( m_ij_vector, 1, &globalRow, &value ) );
}

void
HypreVector::set( globalIndex const * globalIndices,
                  real64 const * values,
                  localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE( *std::min_element( globalIndices, globalIndices + size ),
                       ilower() );
  GEOSX_LAI_ASSERT_GE(
    iupper(),
    getLocalRowID( *std::max_element( globalIndices, globalIndices + size ) ) );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorSetValues( m_ij_vector,
                             LvArray::integerConversion< HYPRE_Int >( size ),
                             toHYPRE_BigInt( globalIndices ),
                             values ) );
}

void
HypreVector::add( globalIndex const * globalIndices,
                  real64 const * values,
                  localIndex size )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorAddToValues( m_ij_vector,
                               LvArray::integerConversion< HYPRE_Int >( size ),
                               toHYPRE_BigInt( globalIndices ),
                               values ) );
}

void
HypreVector::set( arraySlice1d< globalIndex const > const & globalIndices,
                  arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_ASSERT_GE(
    *std::min_element( globalIndices.dataIfContiguous(),
                       globalIndices.dataIfContiguous() + globalIndices.size() ),
    ilower() );
  GEOSX_LAI_ASSERT_GT(
    iupper(),
    *std::max_element( globalIndices.dataIfContiguous(),
                       globalIndices.dataIfContiguous() + globalIndices.size() ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorSetValues(
    m_ij_vector,
    LvArray::integerConversion< HYPRE_Int >( values.size() ),
    toHYPRE_BigInt( globalIndices.dataIfContiguous() ),
    values.dataIfContiguous() ) );
}

void
HypreVector::add( arraySlice1d< globalIndex const > const & globalIndices,
                  arraySlice1d< real64 const > const & values )
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAddToValues(
    m_ij_vector,
    LvArray::integerConversion< HYPRE_Int >( values.size() ),
    toHYPRE_BigInt( globalIndices.dataIfContiguous() ),
    values.dataIfContiguous() ) );
}

void
HypreVector::set( real64 value )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorSetConstantValues( m_par_vector, value ) );
}

void
HypreVector::zero()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( hypre_IJVectorZeroValues( m_ij_vector ) );
}

void
HypreVector::rand( unsigned const seed )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorSetRandomValues( m_par_vector, seed ) );
}

void
HypreVector::open()
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorInitialize( m_ij_vector ) );
  m_closed = false;
}

void
HypreVector::close()
{
  GEOSX_LAI_ASSERT( !closed() );
  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorAssemble( m_ij_vector ) );
  m_closed = true;
}

void
HypreVector::scale( real64 const scalingFactor )
{
  GEOSX_LAI_ASSERT( ready() );

  if( isEqual( scalingFactor, 1.0 ) )
  {
    return;
  }

  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorScale( scalingFactor, m_par_vector ) );
}

void
HypreVector::reciprocal()
{
  GEOSX_LAI_ASSERT( ready() );
  real64 * const values = extractLocalVector();
  for( localIndex i = 0; i < localSize(); ++i )
  {
    values[i] = 1.0 / values[i];
  }
}

real64
HypreVector::dot( HypreVector const & vec ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( vec.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), vec.globalSize() );

  HYPRE_Real result;
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_ParVectorInnerProd( m_par_vector, vec.m_par_vector, &result ) );
  return result;
}

void
HypreVector::copy( HypreVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorCopy( x.m_par_vector, m_par_vector ) );
}

void
HypreVector::axpy( real64 const alpha, HypreVector const & x )
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( x.ready() );
  GEOSX_LAI_ASSERT_EQ( globalSize(), x.globalSize() );

  GEOSX_LAI_CHECK_ERROR( HYPRE_ParVectorAxpy( alpha, x.m_par_vector, m_par_vector ) );
}

void
HypreVector::axpby( real64 const alpha, HypreVector const & x, real64 const beta )
{
  scale( beta );
  axpy( alpha, x );
}

real64
HypreVector::norm1() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 const * const local_data = extractLocalVector();
  real64 loc_norm1 = 0.0;
  for( HYPRE_Int i = 0; i < localSize(); ++i )
  {
    loc_norm1 += std::fabs( local_data[i] );
  }
  return MpiWrapper::Sum( loc_norm1, getComm() );
}

real64
HypreVector::norm2() const
{
  GEOSX_LAI_ASSERT( ready() );
  return std::sqrt( dot( *this ) );
}

real64
HypreVector::normInf() const
{
  GEOSX_LAI_ASSERT( ready() );

  real64 const * const local_data = extractLocalVector();
  real64 loc_normInf = 0.0;
  for( HYPRE_Int i = 0; i < localSize(); ++i )
  {
    loc_normInf = std::max( loc_normInf, std::fabs( local_data[i] ) );
  }
  return MpiWrapper::Max( loc_normInf, getComm() );
}

void
HypreVector::print( std::ostream & os ) const
{
  GEOSX_LAI_ASSERT( ready() );

  int const this_mpi_process = MpiWrapper::Comm_rank( getComm() );
  int const n_mpi_process = MpiWrapper::Comm_size( getComm() );
  char str[77];

  if( this_mpi_process == 0 )
  {
    os << "MPI_Process         GlobalRowID         GlobalColID                 "
          "  Value"
       << std::endl;
  }

  for( int iRank = 0; iRank < n_mpi_process; iRank++ )
  {
    MpiWrapper::Barrier( getComm() );
    if( iRank == this_mpi_process )
    {
      real64 const * const local_data = extractLocalVector();
      globalIndex const firstRowID = ilower();
      for( localIndex i = 0; i < localSize(); ++i )
      {
        sprintf( str,
                 "%11i%20lli%24.10e\n",
                 iRank,
                 firstRowID + LvArray::integerConversion< globalIndex >( i ),
                 local_data[i] );
        os << str;
      }
    }
  }
}

void
HypreVector::write( string const & filename, LAIOutputFormat const format ) const
{
  GEOSX_LAI_ASSERT( ready() );
  switch( format )
  {
    case LAIOutputFormat::NATIVE_ASCII:
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorPrint( m_ij_vector, filename.c_str() ) );
      break;
    }
    case LAIOutputFormat::MATRIX_MARKET:
    {
      if( globalSize() == 0 )
      {
        if( MpiWrapper::Comm_rank( getComm() ) == 0 )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );
          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
          hypre_fprintf( fp, "%d %d\n", 0, 1 );
          std::fclose( fp );
        }
      }
      else
      {
        // Copy distributed parVector in a local vector on every process
        // with at least one component
        // Warning: works for a parVector that is smaller than 2^31-1
        hypre_Vector * vector;
        vector = hypre_ParVectorToVectorAll( m_par_vector );

        // Identify the smallest process where vector exists
        int myID = MpiWrapper::Comm_rank( getComm() );
        if( vector == 0 )
        {
          myID = MpiWrapper::Comm_size( getComm() );
        }
        int printID = MpiWrapper::Min( myID, getComm() );

        // Write to file vector
        if( MpiWrapper::Comm_rank( getComm() ) == printID )
        {
          FILE * fp = std::fopen( filename.c_str(), "w" );
          HYPRE_Real * data = hypre_VectorData( vector );
          HYPRE_Int size = hypre_VectorSize( vector );

          hypre_fprintf( fp, "%s", "%%MatrixMarket matrix array real general\n" );
          hypre_fprintf( fp, "%d %d\n", size, 1 );

          for( HYPRE_Int i = 0; i < size; i++ )
          {
            hypre_fprintf( fp, "%.16e\n", data[i] );
          }

          std::fclose( fp );
        }

        // Destroy vector
        if( vector )
        {
          GEOSX_LAI_CHECK_ERROR( hypre_SeqVectorDestroy( vector ) );
        }
      }
      break;
    }
    default:
      GEOSX_ERROR( "Unsupported vector output format" );
  }
}

real64
HypreVector::get( globalIndex globalRow ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT_GE( globalRow, ilower() );
  GEOSX_LAI_ASSERT_GT( iupper(), globalRow );

  HYPRE_Real value;
  GEOSX_LAI_CHECK_ERROR(
    HYPRE_IJVectorGetValues( m_ij_vector, 1, toHYPRE_BigInt( &globalRow ), &value ) );
  return value;
}

void
HypreVector::get( arraySlice1d< globalIndex const > const & globalIndices,
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

  GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetValues(
    m_ij_vector,
    LvArray::integerConversion< HYPRE_Int >( globalIndices.size() ),
    toHYPRE_BigInt( globalIndices.dataIfContiguous() ),
    values.dataIfContiguous() ) );
}

HYPRE_ParVector const &
HypreVector::unwrapped() const
{
  return m_par_vector;
}

HYPRE_IJVector const &
HypreVector::unwrappedIJ() const
{
  return m_ij_vector;
}

globalIndex
HypreVector::globalSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJVectorGlobalNumRows( m_ij_vector );
}

localIndex
HypreVector::localSize() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_ParVectorActualLocalSize( m_par_vector );
}

localIndex
HypreVector::getLocalRowID( globalIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  return ( index >= ilower() && index < iupper() )
    ? LvArray::integerConversion< localIndex >( index - ilower() )
    : -1;
}

globalIndex
HypreVector::getGlobalRowID( localIndex const index ) const
{
  GEOSX_LAI_ASSERT( created() );
  return ilower() + index;
}

real64 const *
HypreVector::extractLocalVector() const
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_VectorData( hypre_ParVectorLocalVector( m_par_vector ) );
}

real64 *
HypreVector::extractLocalVector()
{
  GEOSX_LAI_ASSERT( ready() );
  return hypre_VectorData( hypre_ParVectorLocalVector( m_par_vector ) );
}

globalIndex
HypreVector::ilower() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >(
    hypre_ParVectorFirstIndex( m_par_vector ) );
}

globalIndex
HypreVector::iupper() const
{
  GEOSX_LAI_ASSERT( created() );
  return LvArray::integerConversion< globalIndex >(
           hypre_ParVectorLastIndex( m_par_vector ) ) +
    1;
}

MPI_Comm
HypreVector::getComm() const
{
  GEOSX_LAI_ASSERT( created() );
  return hypre_IJVectorComm( m_ij_vector );
}

}  // end namespace geosx
