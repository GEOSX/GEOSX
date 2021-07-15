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
 * @file TimeHistHDF.hpp
 */

#ifndef GEOSX_FILEIO_TIMEHISTORY_HDFFILE_HPP_
#define GEOSX_FILEIO_TIMEHISTORY_HDFFILE_HPP_

#include "LvArray/src/Array.hpp"
#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"
#include "common/MpiWrapper.hpp"

#include "dataRepository/HistoryDataSpec.hpp"
#include "HistoryIO.hpp"
#include <hdf5.h>



namespace geosx
{

/**
 * @class HDFTarget
 * @brief An abstract class representing an HDF output target.
 */
class HDFTarget
{
public:
  /**
   * @brief Get the HDF hid_t of the target.
   * @return The hid_t of the target.
   */
  virtual operator hid_t() { return 0; }

  /**
   * @brief Whether a dataset/group with the specified name exists in the target.
   * @param name The dataset/group name to check for.
   * @return Whether the dataset/group exists in the target.
   */
  virtual bool checkInTarget( const string & name )
  {
    htri_t exists = 0;
    H5E_BEGIN_TRY {
      exists = H5Gget_objinfo( this->operator hid_t(), name.c_str(), 0, NULL );
    } H5E_END_TRY
    return (exists == 0);
  }
};

/**
 * @class HDFFile
 * A class used to control access to an HDF file target.
 */
class HDFFile : public HDFTarget
{
public:
  /**
   * @brief Constructor -- this creates/opens the target file for read/write.
   * @param fnm The filename.
   * @param deleteExisting Whether to remove/recreate if a file with the same name exists.
   * @param parallelAccess Whether to access one file in parallel or one file per rank in the comm.
   * @param comm An MPI communicator where each rank in the communicator will be accesing the target file.
   */
  HDFFile( string const & fnm, bool deleteExisting, bool parallelAccess, MPI_Comm comm );

  /**
   * Destructor -- Close the file and acccessors.
   */
  ~HDFFile();

  /**
   * @brief Get the HDF hid_t file identifier.
   * @return the HDF hid_t file id.
   */
  virtual operator hid_t() final { return m_fileId; }
private:
  /// The filename
  string m_filename;
  /// The hdf file id
  hid_t m_fileId;
  /// The hdf file access properties list id
  hid_t m_faplId;
  /// Whether to open the same file in parallel, or open a file per process
  bool m_mpioFapl;
  /// The comminator to operate on the file collectively over
  MPI_Comm m_comm;
};

/**
 * @class HDFHistIO
 * @brief Perform buffered history I/O for a single type(really just output) on using HDF5.
 */
class HDFHistIO : public BufferedHistoryIO
{
public:
  /**
   * @brief Constructor
   * @param filename The filename to perform history output to.
   * @param rank The rank of the history data being collected.
   * @param dims The dimensional extent for each dimension of the history data being collected.
   * @param name The name to use to create/modify the dataset for the history data.
   * @param typeId The std::type_index(typeid(T)) of the underlying data type.
   * @param writeHead How many time history states have been written to the file (used on restart and to compress data on exit).
   * @param initAlloc How many states to preallocate the internal buffer to hold.
   * @param overallocMultiple Integer to scale the internal buffer when we fill the existing space.
   * @param comm A communicator where every rank will participate in writting to the output file.
   */
  HDFHistIO( string const & filename,
             localIndex rank,
             std::vector< localIndex > const & dims,
             string const & name,
             std::type_index typeId,
             localIndex writeHead = 0,
             localIndex initAlloc = 1,
             localIndex overallocMultiple = 2,
             MPI_Comm comm = MPI_COMM_GEOSX );

  /**
   * @brief Constructor
   * @param filename The filename to perform history output to.
   * @param spec HistoryMetadata to use to call the other constructor.
   * @param writeHead How many time states have been written to the file (used on restart and to compress data on exit).
   * @param initAlloc How many states to preallocate the internal buffer to hold.
   * @param overallocMultiple Integer to scale the internal buffer when we fill the existing space.
   * @param comm A communicator where every rank will participate in writing to the output file.
   */
  HDFHistIO( string const & filename,
             const HistoryMetadata & spec,
             localIndex writeHead = 0,
             localIndex initAlloc = 1,
             localIndex overallocMultiple = 2,
             MPI_Comm comm = MPI_COMM_GEOSX ):
    HDFHistIO( filename,
               spec.getRank(),
               spec.getDims(),
               spec.getName(),
               spec.getType(),
               writeHead,
               initAlloc,
               overallocMultiple,
               comm )
  { }

  /// Destructor
  virtual ~HDFHistIO() { }

  /// @copydoc geosx::BufferedHistoryIO::init
  virtual void init( bool existsOkay ) override;

  /// @copydoc geosx::BufferedHistoryIO::write
  virtual void write( ) override;

  /// @copydoc geosx::BufferedHistoryIO::compressInFile
  virtual void compressInFile( ) override;

  /// @copydoc geosx::BufferedHistoryIO::updateCollectingCount
  virtual void updateCollectingCount( localIndex count ) override;

  /// @copydoc geosx::BufferedHistoryIO::getRowBytes
  virtual size_t getRowBytes( ) override;

protected:

  /**
   * @brief Setup the parallel 'partitioning' of the data to allow dynamically sized output over time
   * @param localIdxCount The number of pieces of data associated with the local rank
   * @note This is collective over the communicator provided in the constructor
   */
  void setupPartition( globalIndex localIdxCount );

  /**
   * @brief Update the extent of the dataset in the target file.
   * @param rowLimit The new discrete 'row' count (the first output dim) to extend the extent to.
   * @note The second dimension is set to the global index highwater ( the largest number of output
   *        pieces of data encountered during execution ).
   */
  void updateDatasetExtent( hsize_t rowLimit );

  /**
   * @brief Resize the dataspace in the target file if needed to perform the current write of buffered states.
   * @param bufferedCount The number of buffered states to use to determine if the file needs to be resized.
   */
  void resizeFileIfNeeded( localIndex bufferedCount );

  virtual void resizeBuffer( ) override;

private:
  // file io params
  /// The filename to write to
  string m_filename;
  /// How much to scale the internal and file allocations by when room runs out
  const localIndex m_overallocMultiple;
  /// The global index offset for this mpi rank for this data set
  globalIndex m_globalIdxOffset;
  /// The global index count for this mpi rank for this data set
  globalIndex m_globalIdxCount;
  /// The largest index count encountered (globally) during execution (this determines the 2nd dimension of
  ///   the output data array)
  globalIndex m_globalIdxHighwater;
  /// The current chunk size for writing to file. This is the smallest (necessarily nonzero) local index
  ///   count from the last partition setup.
  hsize_t m_chunkSize;
  /// The current limit in discrete history counts for this data set in the file
  localIndex m_writeLimit;
  /// The current history count for this data set in the file
  localIndex m_writeHead;
  // history metadata
  /// The underlying data type for this history data set
  hsize_t m_hdfType;
  /// The size in byte of the data type
  size_t m_typeSize;
  /// The number of variables of data type in this data set
  hsize_t m_typeCount;   // prod(dims[0:n])
  /// The rank of the data set
  hsize_t m_rank;
  /// The dimensions of the data set
  std::vector< hsize_t > m_dims;
  ///
  std::vector< globalIndex > m_localIdxCounts_buffered;
  /// The name of the data set
  string m_name;
  /// The communicator across which the data set is distributed
  MPI_Comm m_comm;
  /// The communicator with only members of the m_comm comm which have nonzero ammounts of local data (required for chunking output ->
  /// growing the data size in the file)
  MPI_Comm m_subcomm;
  /// Whether the size of the collected data has changed between writes to file
  bool m_sizeChanged;
};

}

#endif
