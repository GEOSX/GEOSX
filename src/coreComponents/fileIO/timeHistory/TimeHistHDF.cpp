#include "TimeHistHDF.hpp"

#include "fileIO/Outputs/TimeHistoryOutput.hpp"

#include "common/MpiWrapper.hpp"

#include <hdf5.h>

namespace geosx
{

/**
 * @brief Get the HDF data type for the specified type
 * @tparam T the type to get info for
 * @return The HDF-specific data type info for T.
 */
template< typename T >
inline hid_t GetHDFDataType();
template<>
inline hid_t GetHDFDataType< char >() { return H5T_NATIVE_CHAR; }
template<>
inline hid_t GetHDFDataType< signed char >() { return H5T_NATIVE_CHAR; }
template<>
inline hid_t GetHDFDataType< float >() { return H5T_NATIVE_FLOAT; }
template<>
inline hid_t GetHDFDataType< double >() { return H5T_NATIVE_DOUBLE; }
template<>
inline hid_t GetHDFDataType< int >() { return H5T_NATIVE_INT; }
template<>
inline hid_t GetHDFDataType< long >() { return H5T_NATIVE_LONG; }
template<>
inline hid_t GetHDFDataType< long long >() { return H5T_NATIVE_LLONG; }

/**
 * @brief Get the HDF data type from the type_index of the type.
 * @param type the std::type_index(typeid(T)) of the type T
 * @return The HDF data type.
 */
inline hid_t GetHDFDataType( std::type_index const & type )
{
  if( type == std::type_index( typeid(char)) )
  {
    return GetHDFDataType< char >();
  }
  else if( type == std::type_index( typeid(signed char)) )
  {
    return GetHDFDataType< signed char >();
  }
  else if( type == std::type_index( typeid(real32)) )
  {
    return GetHDFDataType< real32 >();
  }
  else if( type == std::type_index( typeid(real64)) )
  {
    return GetHDFDataType< real64 >();
  }
  else if( type == std::type_index( typeid(integer)) )
  {
    return GetHDFDataType< integer >();
  }
  else if( type == std::type_index( typeid(localIndex)) )
  {
    return GetHDFDataType< localIndex >();
  }
  else if( type == std::type_index( typeid(globalIndex)) )
  {
    return GetHDFDataType< globalIndex >();
  }
  else
  {
    return GetHDFDataType< char >();
  }
}

/**
 * @brief Get an HDF array data type.
 * @param type The std::type_index(typeid(T)) of the type T inside the array
 * @param rank The rank of the array.
 * @param dims The extent of each dimension of the array.
 */
inline hid_t GetHDFArrayDataType( std::type_index const & type, hsize_t const rank, hsize_t const * dims )
{
  return H5Tarray_create( GetHDFDataType( type ), rank, dims );
}

HDFFile::HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, MPI_Comm comm ):
  m_filename( ),
  m_fileId( 0 ),
  m_faplId( 0 ),
  m_mpioFapl( parallelAccess ),
  m_comm( comm )
{
  int rnk = MpiWrapper::commRank( comm );
#ifdef GEOSX_USE_MPI
  if( m_mpioFapl )
  {
    m_faplId = H5Pcreate( H5P_FILE_ACCESS );
    H5Pset_fapl_mpio( m_faplId, m_comm, MPI_INFO_NULL );
    m_filename = fnm + ".hdf5";
  }
  else
  {
    m_faplId = H5P_DEFAULT;
    m_filename = fnm + "." + std::to_string( rnk ) + ".hdf5";
  }
#else
  m_faplId = H5P_DEFAULT;
  m_filename = fnm + ".hdf5";
#endif
  // check if file already exists
  htri_t exists = 0;
  H5E_BEGIN_TRY {
    exists = H5Fis_hdf5( m_filename.c_str() );
  } H5E_END_TRY
  // if there is an non-hdf file with the same name,
  // and we're either not using parallel access or we're rank 0
  if( exists == 0 && ( !m_mpioFapl || rnk == 0 ) )
  {
    remove( m_filename.c_str() );
  }
  if( exists > 0 && !deleteExisting )
  {
    m_fileId = H5Fopen( m_filename.c_str(), H5F_ACC_RDWR, m_faplId );
  }
  else if( exists >= 0 && deleteExisting )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, m_faplId );
  }
  else if( exists < 0 )
  {
    m_fileId = H5Fcreate( m_filename.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, m_faplId );
  }
}

HDFFile::~HDFFile()
{
  if( m_mpioFapl )
  {
    H5Pclose( m_faplId );
  }
  H5Fclose( m_fileId );
}

HDFHistIO::HDFHistIO( string const & filename,
                      localIndex rank,
                      std::vector< localIndex > const & dims,
                      string const & name,
                      std::type_index typeId,
                      localIndex writeHead,
                      localIndex initAlloc,
                      localIndex overallocMultiple,
                      MPI_Comm comm ):
  BufferedHistoryIO(),
  m_filename( filename ),
  m_overallocMultiple( overallocMultiple ),
  m_globalIdxOffset( 0 ),
  m_globalIdxCount( 0 ),
  m_globalIdxHighwater( 0 ),
  m_chunkSize( 0 ),
  m_writeLimit( initAlloc ),
  m_writeHead( writeHead ),
  m_hdfType( GetHDFDataType( typeId )),
  m_typeSize( H5Tget_size( m_hdfType )),
  m_typeCount( 1 ),
  m_rank( LvArray::integerConversion< hsize_t >( rank )),
  m_dims( rank ),
  m_localIdxCounts_buffered( ),
  m_name( name ),
  m_comm( comm ),
  m_subcomm( MPI_COMM_NULL ),
  m_sizeChanged( true )
{
  for( hsize_t dd = 0; dd < m_rank; ++dd )
  {
    m_dims[dd] = LvArray::integerConversion< hsize_t >( dims[dd] );
    m_typeCount *= m_dims[dd];
  }
  m_dataBuffer.resize( initAlloc * m_typeSize * m_typeCount );
  m_bufferHead = &m_dataBuffer[0];
}

void HDFHistIO::setupPartition( globalIndex localIdxCount )
{
  int size = MpiWrapper::commSize( m_comm );
  int rank = MpiWrapper::commRank( m_comm );
  int color = MPI_UNDEFINED;

  if( localIdxCount > 0 )
  {
    color = 0;
  }

  std::vector< globalIndex > counts( size );
  MpiWrapper::allgather( &localIdxCount, 1, &counts[0], 1, m_comm );

  m_chunkSize = std::numeric_limits< hsize_t >::max( );
  globalIndex globalIdxCount = 0;
  int key = 0;
  for( int ii = 0; ii < size; ii++ )
  {
    if( counts[ii] != 0 && ii < rank )
    {
      key++;
    }
    if( ii == rank )
    {
      m_globalIdxOffset = globalIdxCount;
    }
    globalIdxCount += counts[ii];
    if( counts[ii] > 0 && LvArray::integerConversion< hsize_t >( counts[ii] ) < m_chunkSize )
    {
      m_chunkSize = LvArray::integerConversion< hsize_t >( counts[ii] );
    }
  }

  // update global index count, account for no indices prior to init
  m_globalIdxCount = globalIdxCount;
  if( m_globalIdxCount == 0 )
  {
    m_chunkSize = 1;
  }

  // update global index highwater
  if( m_globalIdxCount > m_globalIdxHighwater )
  {
    m_globalIdxHighwater = m_globalIdxCount;
  }

  // free the previous partition comm
  if( m_subcomm != MPI_COMM_NULL )
  {
    MpiWrapper::commFree( m_subcomm );
  }
  m_subcomm = MpiWrapper::commSplit( m_comm, color, key );
}

void HDFHistIO::init( bool existsOkay )
{
  setupPartition( m_dims[0] );
  int rank = MpiWrapper::commRank( m_comm );
  MPI_Comm subcomm = m_subcomm;
  if( m_globalIdxCount == 0 && rank == 0 )
  {
    subcomm = MPI_COMM_SELF;
  }
  // create a dataset in the file if needed, don't erase file
  if( subcomm != MPI_COMM_NULL )
  {

    std::vector< hsize_t > historyFileDims( m_rank+1 );
    historyFileDims[0] = LvArray::integerConversion< hsize_t >( m_writeLimit );

    std::vector< hsize_t > dimChunks( m_rank+1 );
    dimChunks[0] = 1;

    for( hsize_t dd = 1; dd < m_rank+1; ++dd )
    {
      // hdf5 doesn't like chunk size 0, hence the subcomm
      dimChunks[dd] = m_dims[dd-1];
      historyFileDims[dd] = m_dims[dd-1];
    }
    dimChunks[1] = m_chunkSize;
    historyFileDims[1] = LvArray::integerConversion< hsize_t >( m_globalIdxCount );

    HDFFile target( m_filename, false, true, subcomm );
    bool inTarget = target.checkInTarget( m_name );
    if( !inTarget )
    {
      hid_t dcplId = 0;
      std::vector< hsize_t > maxFileDims( historyFileDims );
      // chunking is required to create an extensible dataset
      dcplId = H5Pcreate( H5P_DATASET_CREATE );
      H5Pset_chunk( dcplId, m_rank + 1, &dimChunks[0] );
      maxFileDims[0] = H5S_UNLIMITED;
      maxFileDims[1] = H5S_UNLIMITED;
      hid_t space = H5Screate_simple( m_rank+1, &historyFileDims[0], &maxFileDims[0] );
      hid_t dataset = H5Dcreate( target, m_name.c_str(), m_hdfType, space, H5P_DEFAULT, dcplId, H5P_DEFAULT );
      H5Dclose( dataset );
      H5Sclose( space );
    }
    else if( existsOkay )
    {
      updateDatasetExtent( m_writeLimit );
    }
    GEOSX_ERROR_IF( inTarget && !existsOkay, "Dataset (" + m_name + ") already exists in output file: " + m_filename );
  }
}

void HDFHistIO::write( )
{
  // check if the size has changed on any process in the primary comm
  bool anyChanged = false;
  MpiWrapper::allReduce( &m_sizeChanged, &anyChanged, 1, MPI_LOR, m_comm );
  m_sizeChanged = anyChanged;

  // this will set the first dim large enough to hold all the rows we're about to write
  resizeFileIfNeeded( m_bufferedCount );
  if( m_bufferedCount > 0 )
  {
    buffer_unit_type * dataBuffer = nullptr;
    if( m_dataBuffer.size() > 0 )
    {
      dataBuffer = &m_dataBuffer[0];
    }
    for( localIndex row = 0; row < m_bufferedCount; ++row )
    {
      // if the size changed at all, update the partitioning and dataset extent before each row is to be written
      //  to ensure the correct mpi ranks participate and that there is enough room to write the largest row during execution
      if( m_sizeChanged )
      {
        // since the highwater might change (the max # of indices / 2nd dimension) when updating the partitioning
        setupPartition( m_localIdxCounts_buffered[ row ] );
        // keep the write limit the same (will only change in resizeFileIfNeeded call above)
        updateDatasetExtent( m_writeLimit );
      }

      if( m_subcomm != MPI_COMM_NULL )
      {
        HDFFile target( m_filename, false, true, m_subcomm );

        hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
        hid_t filespace = H5Dget_space( dataset );

        std::vector< hsize_t > fileOffset( m_rank+1 );
        fileOffset[0] = LvArray::integerConversion< hsize_t >( m_writeHead );
        // the m_globalIdxOffset will be updated for each row during the partition setup if the size has changed during buffered collection
        fileOffset[1] = LvArray::integerConversion< hsize_t >( m_globalIdxOffset );

        std::vector< hsize_t > bufferedCounts( m_rank+1 );
        bufferedCounts[0] = LvArray::integerConversion< hsize_t >( 1 );
        bufferedCounts[1] = LvArray::integerConversion< hsize_t >( m_localIdxCounts_buffered[ row ] );
        for( hsize_t dd = 2; dd < m_rank+1; ++dd )
        {
          bufferedCounts[dd] = m_dims[dd-1];
        }
        hid_t memspace = H5Screate_simple( m_rank+1, &bufferedCounts[0], nullptr );

        hid_t fileHyperslab = filespace;
        H5Sselect_hyperslab( fileHyperslab, H5S_SELECT_SET, &fileOffset[0], nullptr, &bufferedCounts[0], nullptr );

        H5Dwrite( dataset, m_hdfType, memspace, fileHyperslab, H5P_DEFAULT, dataBuffer );

        // forward the data buffer pointer to the start of the next row
        if( dataBuffer )
        {
          dataBuffer += m_localIdxCounts_buffered[ row ] * m_typeSize;
        }

        // unfortunately have to close/open the file for each row since the accessing mpi ranks and extents can change over time
        H5Sclose( memspace );
        H5Sclose( filespace );
        H5Dclose( dataset );
      }

      m_writeHead++;
    }
  }
  m_sizeChanged = false;
  m_localIdxCounts_buffered.clear( );
  emptyBuffer( );
}

void HDFHistIO::compressInFile( )
{
  // set the write limit in the file to the current write head
  updateDatasetExtent( m_writeHead );
  m_writeLimit = m_writeHead;
}

inline void HDFHistIO::resizeFileIfNeeded( localIndex buffered_count )
{
  if( m_writeHead + buffered_count > m_writeLimit )
  {
    while( m_writeHead + buffered_count > m_writeLimit )
    {
      m_writeLimit *= m_overallocMultiple;
    }
    updateDatasetExtent( m_writeLimit );
  }
}

void HDFHistIO::updateDatasetExtent( hsize_t rowLimit )
{
  if( m_subcomm != MPI_COMM_NULL )
  {
    HDFFile target( m_filename, false, true, m_subcomm );
    std::vector< hsize_t > maxFileDims( m_rank+1 );
    maxFileDims[0] = rowLimit;
    maxFileDims[1] = LvArray::integerConversion< hsize_t >( m_globalIdxHighwater );
    for( hsize_t dd = 2; dd < m_rank+1; ++dd )
    {
      maxFileDims[dd] = m_dims[dd-1];
    }
    hid_t dataset = H5Dopen( target, m_name.c_str(), H5P_DEFAULT );
    H5Dset_extent( dataset, &maxFileDims[0] );
    H5Dclose( dataset );
  }
}

size_t HDFHistIO::getRowBytes( )
{
  return m_typeCount * m_typeSize;
}

void HDFHistIO::resizeBuffer( )
{
  // need to store the count every time we collect (which calls this)
  //  regardless of whether the size changes
  m_localIdxCounts_buffered.emplace_back( m_dims[0] );

  size_t capacity = m_dataBuffer.size();
  size_t inUse = m_bufferHead - &m_dataBuffer[0];
  size_t nextRow = getRowBytes( );
  // if needed, resize the buffer
  if( inUse + nextRow > capacity )
  {
    // resize based on the ammount currently in use rather than on the capacity ( less aggressive w/ changing sizes )
    m_dataBuffer.resize( inUse + ( nextRow * m_overallocMultiple ) );
  }
  // reset the buffer head and advance based on count in case the underlying data buffer moves during a resize
  m_bufferHead = &m_dataBuffer[0] + inUse;
}

void HDFHistIO::updateCollectingCount( localIndex count )
{
  if( LvArray::integerConversion< hsize_t >( count ) != m_dims[0] )
  {
    m_sizeChanged = true;
    m_dims[0] = count;
    m_typeCount = count;
    for( hsize_t dd = 1; dd < m_rank; ++dd )
    {
      m_typeCount *= m_dims[dd];
    }
  }
}

}
