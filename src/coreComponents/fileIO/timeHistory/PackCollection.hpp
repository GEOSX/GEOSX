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
 * @file PackCollection.hpp
 */

#ifndef GEOSX_FILEIO_TIMEHISTORY_PACKCOLLECTION_HPP_
#define GEOSX_FILEIO_TIMEHISTORY_PACKCOLLECTION_HPP_

#include "TimeHistoryCollection.hpp"

namespace geosx
{

class DomainPartition;

/**
 * @class PackCollection
 *
 * A task class for serializing history information using the MPI communication packing routines.
 */
class PackCollection : public HistoryCollection
{
public:
  /**
   * @brief Constructor
   * @copydetails dataRepository::Group::Group( string const & name, Group * parent );
   */
  PackCollection ( string const & name, Group * parent );

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string catalogName() { return "PackCollection"; }

  virtual void initializePostSubGroups() override;

  /// @copydoc geosx::HistoryCollection::getMetadata
  virtual HistoryMetadata getMetadata( DomainPartition const & domain, localIndex collectionIdx ) override;

  /// @copydoc geosx::HistoryCollection::getTargetName
  virtual const string & getTargetName( ) const override
  {
    return m_fieldName;
  }

  /**
   * @brief Update the indices related to the sets being collected.
   * @param domain The domain partition.
   * @note This is only required because we don't want to copy/move the
   *       indices each collection execution, because that causes data movement
   *       when collecting data from the device.
   * @note Refactoring the packing functions to allow direct usage of set indices
   *       from SortedArrayView instead of only ArrayViews will remove this
   *       duplication.
   */
  virtual void updateSetsIndices( DomainPartition & domain ) override final;

  virtual localIndex getNumMetaCollectors( ) const override final;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct
  {
    static constexpr char const * objectPathString() { return "objectPath"; }
    static constexpr char const * fieldNameString() { return "fieldName"; }
    static constexpr char const * setNamesString() { return "setNames"; }
    static constexpr char const * onlyOnSetChangeString() { return "onlyOnSetChange"; }

    dataRepository::ViewKey objectPath = { "objectPath" };
    dataRepository::ViewKey fieldName = { "fieldName" };
    dataRepository::ViewKey setNames = { "setNames" };
    dataRepository::ViewKey onlyOnSetChange = { "onlyOnSetChange" };
  } viewKeys;
  /// @endcond

protected:
  /// Construct the metadata collectors for this collector.
  void buildMetaCollectors( );

  /// Do not construct metadata collectors to collect coordinate information.
  ///   ( Prevents reccuring initialization of coordinate collection for coordinate collectors ).
  void disableCoordCollection( )
  {
    m_disableCoordCollection = true;
  }

  /// @copydoc geosx::HistoryCollection::collect
  virtual void collect( DomainPartition & domain,
                        real64 const time_n,
                        real64 const dt,
                        localIndex const collectionIdx,
                        buffer_unit_type * & buffer ) override;

private:
  // todo : replace this with a vector of references to the actual set sortedarrays (after packing rework to allow sorted arrays to be used
  // for indexing)
  /// The indices for the specified sets to pack
  std::vector< array1d< localIndex > > m_setsIndices;
  /// The dataRepository name/path to get history data from, relative paths are assumed to be relative to mesh body 0, mesh level 0
  string m_objectPath;
  /// The (packable) field associated with the specified object to get data from
  string m_fieldName;
  /// The names of the sets to collect history info from
  string_array m_setNames;
  /// Whether any of the collected set(s) have changed since the last collection
  bool m_setChanged;
  /// Whether to only pack when the collected set(s) change size (mostly used for collecting metadata)
  localIndex m_onlyOnSetChange;
  /// Whether to create coordinate meta-collectors if collected objects are mesh objects (set to true for coordinate meta-collectors to
  /// avoid init recursion)
  bool m_disableCoordCollection;
  /// Whether initializePostSubGroups has been called, since we only wan't to execute it once
  ///  It is called explicitly by the output to ensure this is in a valid state to collect info from to perform setup
  ///  It is also called by the normal initialization process
  bool m_initialized;
};

}
#endif
