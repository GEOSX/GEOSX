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
 * @file ObjectManagerBase.hpp
 */

#ifndef GEOSX_MESH_OBJECTMANAGERBASE_HPP_
#define GEOSX_MESH_OBJECTMANAGERBASE_HPP_

#include "dataRepository/Group.hpp"
#include "common/TimingMacros.hpp"
#include "mpiCommunications/NeighborData.hpp"

namespace geosx
{
class SiloFile;
class NodeManager;

/**
 * @brief The ObjectManagerBase is the base object of all object managers in the mesh data hierachy.
 */
class ObjectManagerBase : public dataRepository::Group
{
public:
  ObjectManagerBase() = delete;

  /**
   * @brief Constructor.
   * @param[in] name Name of this object manager
   * @param[in] parent Parent Group
   */
  explicit ObjectManagerBase( string const & name,
                              dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ObjectManagerBase() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  /**
   * @brief Nested type for the `factory` pattern, defining the base class (ObjectManagerBase)
   *        and the builder arguments (string const &, dataRepository::Group * const) of the derived products.
   */
  using CatalogInterface = dataRepository::CatalogInterface< ObjectManagerBase, string const &, dataRepository::Group * const >;

  /**
   * @brief Acessing the unique instance of this catalog.
   * @return A reference to the singleton.
   */
  static CatalogInterface::CatalogType & getCatalog();

  /**
   * @brief Get the name of the catalog.
   * @return The name.
   */
  virtual const string getCatalogName() const = 0;
  ///@}

  using dataRepository::Group::packSize;
  using dataRepository::Group::pack;

  virtual localIndex packSize( string_array const & wrapperNames,
                               arrayView1d< localIndex const > const & packList,
                               integer const recursive,
                               bool onDevice,
                               parallelDeviceEvents & events ) const override;

  virtual localIndex pack( buffer_unit_type * & buffer,
                           string_array const & wrapperNames,
                           arrayView1d< localIndex const > const & packList,
                           integer const recursive,
                           bool onDevice,
                           parallelDeviceEvents & events ) const override;

  virtual localIndex unpack( buffer_unit_type const * & buffer,
                             arrayView1d< localIndex > & packList,
                             integer const recursive,
                             bool onDevice,
                             parallelDeviceEvents & events ) override;

  /**
   * @brief Packs the elements of each set that actually are in @p packList.
   * @tparam DOPACK Template parameter that decides at compile time whether one should actually pack or not.
   * @param buffer The buffer that will store the packed data.
   * @param packList The elements of each set that should be packed.
   * @return The size of the (potentially) packed data.
   *
   * Note that the returned value does not depend on parameter @p DOPACK.
   */
  template< bool DOPACK >
  localIndex packSets( buffer_unit_type * & buffer,
                       arrayView1d< localIndex const > const & packList ) const;

  /**
   * @brief Unpack the content of @p buffer into the sets of the instance.
   * @param buffer The buffer containing the packed data.
   * @return The unpacked size.
   */
  localIndex unpackSets( buffer_unit_type const * & buffer );

  /**
   * @brief Inserts in @p exclusionList the data that shall not be packed.
   * @param[in,out] exclusionList Will receive the wrapper indices of the data that should not be packed.
   *
   * Note that data will be inserted into @p exclusionList
   * and that data previously present in @p exclusionList may remain.
   */
  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const;

  /**
   * @brief Computes the pack size of the global maps elements in the @ packList.
   * @param packList The element we want packed.
   * @param recursive Boolean like integer for sub-groups packing.
   * @return The packed size.
   */
  virtual localIndex packGlobalMapsSize( arrayView1d< localIndex const > const & packList,
                                         integer const recursive ) const;

  /**
   * @brief Packs the global maps elements in the @ packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @param recursive Boolean like integer for sub-groups packing.
   * @return The packed size.
   */
  virtual localIndex packGlobalMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList,
                                     integer const recursive ) const;

  /**
   * @brief Clear and redefines the ghosts to receive.
   */
  void setReceiveLists();

  /**
   * @brief Computes the pack size of the specific elements in the @ packList.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
  {
    GEOSX_UNUSED_VAR( packList );
    return 0;
  }

  /**
   * @brief Packs the specific elements in the @ packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const
  {
    GEOSX_UNUSED_VAR( buffer );
    GEOSX_UNUSED_VAR( packList );
    return 0;
  }

  /**
   * @brief Unpacks the specific elements in the @ packList.
   * @param buffer The buffer containing the packed data.
   * @param packList The (un)packed element.
   * @param overwriteUpMaps Clear the up maps provided.
   * @param overwriteDownMaps Clear the down maps provided.
   * @return The packed size.
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps )
  {
    GEOSX_UNUSED_VAR( buffer );
    GEOSX_UNUSED_VAR( packList );
    GEOSX_UNUSED_VAR( overwriteUpMaps );
    GEOSX_UNUSED_VAR( overwriteDownMaps );
    return 0;
  }

  /**
   * @brief Unpacks the global maps from @p buffer.
   * @param buffer The buffer containing the packed data.
   * @param packList The (un)packed element.
   * @param recursive Boolean like integer for sub-groups unpacking.
   * @return The unpacked size.
   */
  virtual localIndex unpackGlobalMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       integer const recursive );

  /**
   * @brief Computes the pack size of the parent/child relations in @p packList.
   * @param packList The indices we want packed.
   * @return The packed size.
   */
  localIndex packParentChildMapsSize( arrayView1d< localIndex const > const & packList ) const
  {
    buffer_unit_type * buffer = nullptr;
    return packParentChildMapsPrivate< false >( buffer, packList );
  }

  /**
   * @brief Packs the parent/child relations in @p packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The indices we want packed.
   * @return The packed size.
   */
  localIndex packParentChildMaps( buffer_unit_type * & buffer,
                                  arrayView1d< localIndex const > const & packList ) const
  {
    return packParentChildMapsPrivate< true >( buffer, packList );
  }

  /**
   * @brief Unacks the parent/child relations in @p packList.
   * @param buffer The buffer containing the packed data.
   * @param packList The unpacked indices.
   * @return
   */
  localIndex unpackParentChildMaps( buffer_unit_type const * & buffer,
                                    localIndex_array & packList );

private:
  /**
   * @brief Concrete implementation of the packing method.
   * @tparam DOPACK A template parameter to discriminate between actually packing or only computing the packing size.
   * @param buffer The buffer that will receive the packed data.
   * @param wrapperNames
   * @param packList The element we want packed.
   * @param recursive recursive pack or not.
   * @param onDevice Whether to use device-based packing functions
   *                  (buffer must be either pinned or a device pointer)
   * @return The packed size.
   */
  template< bool DOPACK >
  localIndex packPrivate( buffer_unit_type * & buffer,
                          string_array const & wrapperNames,
                          arrayView1d< localIndex const > const & packList,
                          integer const recursive,
                          bool onDevice,
                          parallelDeviceEvents & events ) const;

  /**
   * @brief Packing global maps.
   * @tparam DOPACK A template parameter to discriminate between actually packing or only computing the packing size.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @param recursive recursive pack or not.
   * @return The packed size.
   */
  template< bool DOPACK >
  localIndex packGlobalMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList,
                                    integer const recursive ) const;

  /**
   * @brief Pack parent and child maps.
   * @tparam DOPACK A template parameter to discriminate between actually packing or only computing the packing size.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @return The packed size.
   */
  template< bool DOPACK >
  localIndex packParentChildMapsPrivate( buffer_unit_type * & buffer,
                                         arrayView1d< localIndex const > const & packList ) const;

  //**********************************************************************************************************************
  // functions for compatibility with old data structure
  // TODO Deprecate or modernize all these suckers

public:

  /**
   * @brief Manually move all sets to a memory space.
   * @param targetSpace The memory space to move sets to.
   */
  void moveSets( LvArray::MemorySpace const targetSpace );

  /**
   * @copydoc geosx::dataRepository::Group::resize(indexType const)
   * @return Always 0, whatever the new size is.
   */
  localIndex resize( localIndex const newSize,
                     const bool /*assignGlobals*/ )
  {
    dataRepository::Group::resize( newSize );
    return 0;
  }

  using dataRepository::Group::resize;

  /**
   * @brief Creates a new set.
   * @param newSetName The set name.
   * @return reference to the set
   */
  SortedArray< localIndex > & createSet( const string & newSetName );

  /**
   * @brief Builds a new set on this instance given another objects set and the map between them.
   * @param inputSet The input set.
   * @param map The map between the newly created set and the @p inputSet.
   * @param setName The newly created set name.
   */
  void constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                  const array2d< localIndex > & map,
                                  const string & setName );
  /**
   * @brief Builds a new set on this instance given another objects set and the map between them.
   * @param inputSet The input set.
   * @param map The map between the newly created set and the @p inputSet.
   * @param setName The newly created set name.
   */
  void constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                  const array1d< localIndex_array > & map,
                                  const string & setName );
  /**
   * @brief Builds a new set on this instance given another objects set and the map between them.
   * @param inputSet The input set.
   * @param map The map between the newly created set and the @p inputSet.
   * @param setName The newly created set name.
   */
  void constructSetFromSetAndMap( SortedArrayView< localIndex const > const & inputSet,
                                  ArrayOfArraysView< localIndex const > const & map,
                                  const string & setName );

  /**
   * @brief Constructs the global to local map.
   */
  void constructGlobalToLocalMap();

  /**
   * @brief Computes the (local) index list that are domain boundaries.
   * @param[in,out] objectList Container that is filled with the local indices.
   *
   * Note that @p objectList is not cleared and domain boundary indices are only appended.
   */
  void constructLocalListOfBoundaryObjects( localIndex_array & objectList ) const;

  /**
   * @brief Computes the (global) index list that are domain boundaries.
   * @param[in,out] objectList Sorted container that is filled with the global indices.
   *
   * Note that @p objectList is not cleared and domain boundary indices are only appended.
   */
  void constructGlobalListOfBoundaryObjects( globalIndex_array & objectList ) const;

  /**
   * @brief Extract map from object and assign global indices.
   * @param nodeManager The node manager.
   * @param map The map.
   *
   * Dummy version, needs to be specialised by derived classes.
   */
  virtual void extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                                std::vector< std::vector< globalIndex > > & map )
  {
    GEOSX_UNUSED_VAR( nodeManager );
    GEOSX_UNUSED_VAR( map );
  }

  /**
   * @brief Defines @p neighborRank ownership for ghost objects.
   * @param neighborRank The rank that owns the ghost objects.
   */
  void setGhostRankForSenders( int const neighborRank )
  {
    arrayView1d< localIndex const > const ghostsToSend = getNeighborData( neighborRank ).ghostsToSend();
    array1d< std::pair< globalIndex, int > > & nonLocalGhosts = getNeighborData( neighborRank ).nonLocalGhosts();
    nonLocalGhosts.clear();

    for( localIndex const index : ghostsToSend )
    {
      integer & owningRank = m_ghostRank[ index ];
      if( owningRank >= 0 )
      {
        nonLocalGhosts.emplace_back( m_localToGlobalMap[ index ], owningRank );
      }
      else
      {
        owningRank = -1;
      }
    }
  }

  /**
   * @brief Get the number of ghost objects.
   * @return The number of ghost objects.
   */
  localIndex getNumberOfGhosts() const;

  /**
   * @brief Get the number of locally owned objects.
   * @return The number of locally owned objects.
   */
  localIndex getNumberOfLocalIndices() const;

  /**
   * @brief Split object to deal with topology changes.
   * @param indexToSplit Where to split.
   * @param rank MPI rank.
   * @param newIndex At which index the we'll have the new split instance.
   * @return Always 1.
   */
  integer splitObject( localIndex const indexToSplit,
                       int const rank,
                       localIndex & newIndex );

  /**
   * @brief sets the value of m_ghostRank to the value of the objects parent.
   * @param indices the list of indices for which to set the ghost rank
   *
   * This function takes a list of indices, and then sets the value of m_ghostRank for those indices to be equal to the
   * value of the "parent" index. This assumes that "parentIndex" is allocated and filled correctly.
   */
  void inheritGhostRankFromParent( std::set< localIndex > const & indices );

  /**
   * @brief Copy object from @p source to @ destination
   * @param source The source index.
   * @param destination The destination index.
   */
  void copyObject( localIndex const source, localIndex const destination );

  /**
   * @brief Computes the maximum global index allong all the MPI ranks.
   */
  virtual void setMaxGlobalIndex();

  /**
   * @brief Fixing the up/down maps by mapping the unmapped indices.
   * @tparam TYPE_RELATION Some InterObjectRelation template class instance.
   * @param relation Global to local indices relations.
   * @param unmappedIndices Unmapped indices we will map during this function call.
   * @param clearIfUnmapped Shall we clear the unmapped indices. Here unused.
   */
  template< typename TYPE_RELATION >
  static void fixUpDownMaps( TYPE_RELATION & relation,
                             map< localIndex, array1d< globalIndex > > & unmappedIndices,
                             bool const clearIfUnmapped );

  /**
   * @brief Fixing the up/down maps by mapping the unmapped indices.
   * @tparam TYPE_RELATION Some InterObjectRelation template class instance.
   * @param relation Global to local indices relations.
   * @param unmappedIndices Unmapped indices we will map during this function call.
   * @param clearIfUnmapped Shall we clear the unmapped indices.
   */
  template< typename TYPE_RELATION >
  static void fixUpDownMaps( TYPE_RELATION & relation,
                             map< localIndex, SortedArray< globalIndex > > & unmappedIndices,
                             bool const clearIfUnmapped );

  /**
   * @brief Fixing the up/down maps by mapping the unmapped indices.
   * @param relation Global to local indices relations.
   * @param globalToLocal This map is used instead of the map provided by @p relation.
   * @param unmappedIndices Unmapped indices we will map during this function call.
   * @param clearIfUnmapped Shall we clear the unmapped indices. Here unused.
   */
  static void fixUpDownMaps( ArrayOfSets< localIndex > & relation,
                             unordered_map< globalIndex, localIndex > const & globalToLocal,
                             map< localIndex, SortedArray< globalIndex > > & unmappedIndices,
                             bool const clearIfUnmapped );

  /**
   * @brief Removes from the list of arrays of @p upmap all the elements
   * for which the "mirror target array" of @p downmap does not contain the proper target index.
   * @param targetIndices The indices we want to keep.
   * @param[in,out] upmap The map to be filtered
   * @param downmap The map used to check for target availability.
   */
  static void cleanUpMap( std::set< localIndex > const & targetIndices,
                          array1d< SortedArray< localIndex > > & upmap,
                          arrayView2d< localIndex const > const & downmap );

  /**
   * @brief Removes from the list of sets of @p upmap all the elements
   * for which the "mirror target array" of @p downmap does not contain the proper target index.
   * @param targetIndices The indices we want to keep.
   * @param[in,out] upmap The map to be filtered
   * @param downmap The map used to check for target availability.
   */
  static void cleanUpMap( std::set< localIndex > const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          arrayView2d< localIndex const > const & downmap );

  /**
   * @brief Removes from the list of arrays of @p upmap all the elements
   * for which the "mirror target array" of @p downmap does not contain the proper target index.
   * @param targetIndices The indices we want to keep.
   * @param[in,out] upmap The map to be filtered
   * @param downmap The map used to check for target availability.
   */
  static void cleanUpMap( std::set< localIndex > const & targetIndices,
                          array1d< SortedArray< localIndex > > & upmap,
                          arrayView1d< arrayView1d< localIndex const > const > const & downmap );
  /**
   * @brief Removes from the list of sets of @p upmap all the elements
   * for which the "mirror target array" of @p downmap does not contain the proper target index.
   * @param targetIndices The indices we want to keep.
   * @param[in,out] upmap The map to be filtered
   * @param downmap The map used to check for target availability.
   */
  static void cleanUpMap( std::set< localIndex > const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          arrayView1d< arrayView1d< localIndex const > const > const & downmap );
  /**
   * @brief Removes from the list of sets of @p upmap all the elements
   * for which the "mirror target array" of @p downmap does not contain the proper target index.
   * @param targetIndices The indices we want to keep.
   * @param[in,out] upmap The map to be filtered
   * @param downmap The map used to check for target availability.
   */
  static void cleanUpMap( std::set< localIndex > const & targetIndices,
                          ArrayOfSetsView< localIndex > const & upmap,
                          ArrayOfArraysView< localIndex const > const & downmap );

  /**
   * @brief Updates the child and target indices after a topology change.
   * @param targetIndices The indices top update.
   */
  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices );

  /**
   * @brief Get the upmost parent.
   * @param parentIndices The list of parent indices.
   * @param lookup The index for which we are looking for the parent
   * @return The upmost parent index.
   *
   * Get the upmost parent (i.e. parent of parent of parent...) of @p lookup
   * that has no more valid parent in @p parentIndices.
   */
  static localIndex getParentRecusive( arraySlice1d< localIndex const > const & parentIndices,
                                       localIndex const lookup )
  {
    localIndex rval = lookup;

    while( parentIndices[rval] != -1 )
    {
      rval = parentIndices[rval];
    }

    return rval;
  }

  /**
   * @brief Register data with this ObjectManagerBase using a
   *   dataRepository::Wrapper.
   * @tparam MESH_DATA_TRAIT The trait struct that holds the information for
   *   the data being registered with the repository.
   * @param nameOfRegisteringObject The name of the object that is requesting
   *   that this data be registered.
   * @return A wrapper to the type specified by @p MESH_DATA_TRAIT.
   */
  template< typename MESH_DATA_TRAIT >
  dataRepository::Wrapper< typename MESH_DATA_TRAIT::type > &
  registerExtrinsicData( string const & nameOfRegisteringObject )
  {
    // These are required to work-around the need for instantiation of
    // the static constexpr trait components. This will not be required once
    // we move to c++17.

    //constexpr typename MESH_DATA_TRAIT::DataType defaultValue = MESH_DATA_TRAIT::defaultValue;
    // This is required for the Tensor classes.
    typename MESH_DATA_TRAIT::dataType defaultValue( MESH_DATA_TRAIT::defaultValue );

    return this->registerWrapper< typename MESH_DATA_TRAIT::type >( MESH_DATA_TRAIT::key() ).
             setApplyDefaultValue( defaultValue ).
             setPlotLevel( MESH_DATA_TRAIT::plotLevel ).
             setRestartFlags( MESH_DATA_TRAIT::restartFlag ).
             setDescription( MESH_DATA_TRAIT::description ).
             setRegisteringObjects( nameOfRegisteringObject );
  }

  /**
   * @brief Register a collection of data with this ObjectManagerBase using a
   *   dataRepository::Wrapper.
   * @tparam MESH_DATA_TRAIT0 The first of the trait structs that holds the
   *   information for the data being registered with the repository.
   * @tparam MESH_DATA_TRAIT1 The second of the trait structs that holds the
   *   information for the data being registered with the repository.
   * @tparam MESH_DATA_TRAITS The parameter pack of trait structs that holds
   *   the information for the data being registered with the repository.
   * @param nameOfRegisteringObject The name of the object that is requesting
   *   that this data be registered.
   */
  template< typename MESH_DATA_TRAIT0, typename MESH_DATA_TRAIT1, typename ... MESH_DATA_TRAITS >
  void registerExtrinsicData( string const & nameOfRegisteringObject )
  {
    registerExtrinsicData< MESH_DATA_TRAIT0 >( nameOfRegisteringObject );
    registerExtrinsicData< MESH_DATA_TRAIT1, MESH_DATA_TRAITS... >( nameOfRegisteringObject );
  }

  /**
   * @brief Get a view to the data associated with a trait from this
   *   ObjectManagerBase.
   * @tparam MESH_DATA_TRAIT The trait that holds the type and key of the data
   *   to be retrieved from this ObjectManagerBase.
   * @return A const reference to a view to const data.
   */
  template< typename MESH_DATA_TRAIT >
  GEOSX_DECLTYPE_AUTO_RETURN getExtrinsicData() const
  {
    return this->getWrapper< typename MESH_DATA_TRAIT::type >( MESH_DATA_TRAIT::key() ).reference();
  }

  /**
   * @brief Get the data associated with a trait from this ObjectManagerBase.
   * @tparam MESH_DATA_TRAIT The trait that holds the type and key of the data
   *   to be retrieved from this ObjectManagerBase.
   * @return A reference to the data.
   */
  template< typename MESH_DATA_TRAIT >
  GEOSX_DECLTYPE_AUTO_RETURN getExtrinsicData()
  {
    return this->getWrapper< typename MESH_DATA_TRAIT::type >( MESH_DATA_TRAIT::key() ).reference();
  }

  /**
   * @brief Checks if an extrinsic data has been registered.
   * @tparam MESH_DATA_TRAIT The trait that holds the type and key of the data
   *   to be retrieved from this ObjectManagerBase.
   * @return @p true if the data has been registered, @p false otherwise.
   */
  template< typename MESH_DATA_TRAIT >
  bool hasExtrinsicData() const
  {
    return this->hasWrapper( MESH_DATA_TRAIT::key() );
  }

#if 0
  template< typename MESH_DATA_TRAIT >
  dataRepository::Wrapper< typename MESH_DATA_TRAIT::Type > &
  registerExtrinsicData( string const & nameOfRegisteringObject,
                         MESH_DATA_TRAIT const & extrinisicDataTrait )
  {
    // These are required to work-around the need for instantiation of
    // the static constexpr trait components. This will not be required once
    // we move to c++17.

    //constexpr typename MESH_DATA_TRAIT::DataType defaultValue = MESH_DATA_TRAIT::defaultValue;
    constexpr dataRepository::PlotLevel plotLevel = MESH_DATA_TRAIT::plotLevel;
    string const description = MESH_DATA_TRAIT::description;

    // This is required for the Tensor classes.
    typename MESH_DATA_TRAIT::DataType defaultValue( MESH_DATA_TRAIT::defaultValue );

    return *(this->registerWrapper< typename MESH_DATA_TRAIT::Type >( extrinisicDataTrait.viewKey ).
               setApplyDefaultValue( defaultValue ).
               setPlotLevel( plotLevel ).
               setDescription( description ).
               setRegisteringObjects( nameOfRegisteringObject ) );
  }

  template< typename MESH_DATA_TRAIT0, typename MESH_DATA_TRAIT1, typename ... MESH_DATA_TRAITS >
  void registerExtrinsicData( string const & nameOfRegisteringObject,
                              MESH_DATA_TRAIT0 const & extrinisicDataTrait0,
                              MESH_DATA_TRAIT1 const & extrinisicDataTrait1,
                              MESH_DATA_TRAITS && ... extrinisicDataTraits )
  {
    registerExtrinsicData< MESH_DATA_TRAIT0 >( nameOfRegisteringObject, extrinisicDataTrait0 );
    registerExtrinsicData< MESH_DATA_TRAIT1,
                           MESH_DATA_TRAITS... >( nameOfRegisteringObject,
                                                  extrinisicDataTrait1,
                                                  std::forward< MESH_DATA_TRAITS >( extrinisicDataTraits )... );
  }

  template< typename MESH_DATA_TRAIT >
  auto const & getExtrinsicData( MESH_DATA_TRAIT const & extrinisicDataTrait ) const
  {
    return this->getWrapper< typename MESH_DATA_TRAIT::Type >( extrinisicDataTrait.viewKey ).referenceAsView();
  }

  template< typename MESH_DATA_TRAIT >
  auto & getExtrinsicData( MESH_DATA_TRAIT const & extrinisicDataTrait )
  {
    return this->getWrapper< typename MESH_DATA_TRAIT::Type >( extrinisicDataTrait.viewKey ).referenceAsView();
  }
#endif

  //**********************************************************************************************************************

  /**
   * @struct viewKeyStruct
   * @brief struct to serve as a container for variable strings and keys
   */
  struct viewKeyStruct
  {
    /// @return String key to adjacency list
    static constexpr char const * adjacencyListString() { return "adjacencyList"; }

    /// @return String key to domain boundary indicator
    static constexpr char const * domainBoundaryIndicatorString() { return "domainBoundaryIndicator"; }

    /// @return String key to external set
    static constexpr char const * externalSetString() { return "externalSet"; }

    /// @return String key to ghost ranks
    static constexpr char const * ghostRankString() { return "ghostRank"; }

    /// @return String key to ghosts to receive
    static constexpr char const * ghostsToReceiveString() { return "ghostsToReceive"; }

    /// @return String key to global->local mao
    static constexpr char const * globalToLocalMapString() { return "globalToLocalMap"; }

    /// @return String key to the 'is external' vector
    static constexpr char const * isExternalString() { return "isExternal"; }

    /// @return String key to the local->global map
    static constexpr char const * localToGlobalMapString() { return "localToGlobalMap"; }

    /// View key to external set
    dataRepository::ViewKey externalSet = { externalSetString() };

    /// View key to ghost ranks
    dataRepository::ViewKey ghostRank = { ghostRankString() };

    /// View key to global->local map
    dataRepository::ViewKey globalToLocalMap = { globalToLocalMapString() };

    /// View key to the local->global map
    dataRepository::ViewKey localToGlobalMap = { localToGlobalMapString() };
  }
  /// viewKey struct for the ObjectManagerBase class
  m_ObjectManagerBaseViewKeys;

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct
  {
    /// @return String key to the Group holding the object sets
    static constexpr char const * setsString() { return "sets"; }

    /// @return String key to the Groupholding all the NeighborData objects
    static constexpr char const * neighborDataString() { return "neighborData"; }

    /// View key to the Group holding the object sets
    dataRepository::GroupKey sets = { setsString() };

    /// View key to the Group holding the neighbor data
    dataRepository::GroupKey neighborData{ neighborDataString() };
  }
  /// groupKey struct for the ObjectManagerBase class
  m_ObjectManagerBaseGroupKeys;

  /**
   * @brief Get the view keys for Group access.
   * @return The keys.
   */
  virtual viewKeyStruct & viewKeys() { return m_ObjectManagerBaseViewKeys; }

  /**
   * @brief Get the view keys for Group access, const version.
   * @return The keys.
   */
  virtual viewKeyStruct const & viewKeys() const { return m_ObjectManagerBaseViewKeys; }
  /**
   * @brief Get the group keys for Group access.
   * @return The keys.
   */
  virtual groupKeyStruct & groupKeys() { return m_ObjectManagerBaseGroupKeys; }
  /**
   * @brief Get the group keys for Group access, const version.
   * @return The keys.
   */
  virtual groupKeyStruct const & groupKeys() const { return m_ObjectManagerBaseGroupKeys; }

  /**
   * @brief Get the group holding the object sets.
   * @return The Group intance holding the sets.
   */
  Group & sets()
  { return m_sets; }

  /**
   * @brief Get the group holding the object sets, const version.
   * @return The Group intance holding the sets.
   */
  Group const & sets() const
  { return m_sets; }

  /**
   * @brief Get a set by name.
   * @param setName Name of the set.
   * @return Sorted array indices.
   */
  SortedArray< localIndex > & getSet( string const & setName )
  { return m_sets.getReference< SortedArray< localIndex > >( setName ); }

  /**
   * @brief Get a set by name, const version.
   * @param setName Name of the set.
   * @return Sorted array indices.
   */
  SortedArrayView< localIndex const > getSet( string const & setName ) const
  { return m_sets.getReference< SortedArray< localIndex > >( setName ).toViewConst(); }

  /**
   * @brief Get the external set.
   * @return Sorted array indices.
   */
  SortedArray< localIndex > & externalSet()
  { return m_sets.getReference< SortedArray< localIndex > >( m_ObjectManagerBaseViewKeys.externalSet ); }

  /**
   * @brief Get the external set, const version.
   * @return Sorted array indices.
   */
  SortedArrayView< localIndex const > externalSet() const
  { return m_sets.getReference< SortedArray< localIndex > >( m_ObjectManagerBaseViewKeys.externalSet ); }

  /**
   * @brief Updates (if needed) the global index for local index @p lid.
   * @param lid The local index
   */
  void updateGlobalToLocalMap( localIndex const lid )
  {
    globalIndex const gid = m_localToGlobalMap[ lid ];
    m_localMaxGlobalIndex = std::max( m_localMaxGlobalIndex, gid );
    m_globalToLocalMap[ gid ] = lid;
  }

  /**
   * @brief Get local to global map.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex > localToGlobalMap()
  { return m_localToGlobalMap; }

  /**
   * @brief Get local to global map, const version.
   * @return The mapping relationship as a array.
   */
  arrayView1d< globalIndex const > localToGlobalMap() const
  { return m_localToGlobalMap; }

  /**
   * @brief Get global to local map.
   * @return The mapping relationship as a array.
   */
  unordered_map< globalIndex, localIndex > const & globalToLocalMap() const
  { return m_globalToLocalMap; }

  /**
   * @brief Retrieves the local index for given global index.
   * @param gid The global index.
   * @return The local index.
   */
  localIndex globalToLocalMap( globalIndex const gid ) const
  { return m_globalToLocalMap.at( gid ); }

  /**
   * @brief Get the locality information of the objects.
   * @return The information is stored as an array of integers, see #m_isExternal
   */
  array1d< integer > const & isExternal()
  { return this->m_isExternal; }

  /**
   * @brief Get the locality information of the objects.
   * @return The information is stored as an array of integers, see #m_isExternal
   */
  arrayView1d< integer const > isExternal() const
  { return this->m_isExternal; }

  /**
   * @brief Get the ghost information of each object.
   * @return See @see #m_ghostRank
   */
  array1d< integer > const & ghostRank()
  { return this->m_ghostRank; }

  /**
   * @brief Get the ghost information of each object, const version.
   * @return See @see #m_ghostRank
   */
  arrayView1d< integer const > ghostRank() const
  { return this->m_ghostRank; }

  /**
   * @brief Get neighbor data for given @p rank.
   * @param rank The rank
   * @return The neighbor data.
   */
  NeighborData & getNeighborData( int const rank )
  { return m_neighborData.at( rank ); }

  /**
   * @brief Get neighbor data for given @p rank, const version.
   * @param rank The rank
   * @return The neighbor data.
   */
  NeighborData const & getNeighborData( int const rank ) const
  { return m_neighborData.at( rank ); }

  /**
   * @brief Add a neighbor for @p rank.
   * @param rank The rank of the new neighbor.
   */
  void addNeighbor( int const rank )
  {
    string const & rankString = std::to_string( rank );
    m_neighborData.emplace( std::piecewise_construct, std::make_tuple( rank ), std::make_tuple( rankString, &m_neighborGroup ) );
    m_neighborGroup.registerGroup( rankString, &getNeighborData( rank ) );
  }

  /**
   * @brief Remove neighbor for @p rank.
   * @param rank The rank of the removed neighbor.
   */
  void removeNeighbor( int const rank )
  {
    m_neighborGroup.deregisterGroup( getNeighborData( rank ).getName() );
    m_neighborData.erase( rank );
  }

  /**
   * @brief Get the local maximum global index on this rank.
   * @return The index.
   */
  globalIndex localMaxGlobalIndex() const
  { return m_localMaxGlobalIndex; }

  /**
   * @brief Get the maximum global index of all objects across all rank. See @see #m_maxGlobalIndex
   * @return The index.
   */
  globalIndex maxGlobalIndex() const
  { return m_maxGlobalIndex; }


  /**
   * @brief Get the domain boundary indicator
   * @return The information in an array of integers, mainly treated as booleans
   *         (1 meaning the "index" is on the boundary).
   */
  arrayView1d< integer > getDomainBoundaryIndicator()
  {
    return m_domainBoundaryIndicator.toView();
  }

  /// @copydoc getDomainBoundaryIndicator()
  arrayView1d< integer const > getDomainBoundaryIndicator() const
  {
    return m_domainBoundaryIndicator.toViewConst();
  }

protected:
  /// Group that holds object sets.
  Group m_sets;

  /// Group that holds all the NeighborData objects.
  Group m_neighborGroup;

  /// Contains the global index of each object.
  array1d< globalIndex > m_localToGlobalMap;

  /// Map from object global index to the local index.
  unordered_map< globalIndex, localIndex > m_globalToLocalMap;

  /// Array that holds if an object is external.
  array1d< integer > m_isExternal;

  /// Domain boundary indicator: 1 means the "index" is on the boundary.
  array1d< integer > m_domainBoundaryIndicator;

  /**
   * @brief Array that holds the ghost information about each object.
   *
   * A value of -2 means that the object is owned locally and not communicated.
   * A value of -1 means that the object is owned locally and is communicated.
   * A positive value means that the object is a ghost and is owned by that rank.
   */
  array1d< integer > m_ghostRank;

  /// A map from rank to the associated NeighborData object.
  unordered_map< int, NeighborData > m_neighborData;

  /// Factor by which to overallocate when adding objects.
  real64 m_overAllocationFactor = 1.1;

  /// The maximum global index of all objects across all rank.
  globalIndex m_maxGlobalIndex = -1;

  /// The maximum global index of any object of all objects on this rank.
  globalIndex m_localMaxGlobalIndex = -1;
};


template< typename TYPE_RELATION >
void ObjectManagerBase::fixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, array1d< globalIndex > > & unmappedIndices,
                                       bool const )
{
  GEOSX_MARK_FUNCTION;

  bool allValuesMapped = true;
  unordered_map< globalIndex, localIndex > const & globalToLocal = relation.relatedObjectGlobalToLocal();
  for( map< localIndex, array1d< globalIndex > >::iterator iter = unmappedIndices.begin();
       iter != unmappedIndices.end();
       ++iter )
  {
    localIndex const li = iter->first;
    array1d< globalIndex > const & globalIndices = iter->second;
    for( localIndex a=0; a<globalIndices.size(); ++a )
    {
      if( globalIndices[a] != unmappedLocalIndexValue )
      {
        if( relation[li][a] == unmappedLocalIndexValue )
        {
          relation[li][a] = globalToLocal.at( globalIndices[a] );
        }
        else
        {
          allValuesMapped = false;
        }
      }
      GEOSX_ERROR_IF( relation[li][a]==unmappedLocalIndexValue, "Index not set" );
    }
  }
  GEOSX_ERROR_IF( !allValuesMapped, "some values of unmappedIndices were not used" );
  unmappedIndices.clear();
}


template< typename TYPE_RELATION >
void ObjectManagerBase::fixUpDownMaps( TYPE_RELATION & relation,
                                       map< localIndex, SortedArray< globalIndex > > & unmappedIndices,
                                       bool const clearIfUnmapped )
{
  GEOSX_MARK_FUNCTION;

  unordered_map< globalIndex, localIndex > const & globalToLocal = relation.RelatedObjectGlobalToLocal();
  for( map< localIndex, SortedArray< globalIndex > >::iterator iter = unmappedIndices.begin();
       iter != unmappedIndices.end();
       ++iter )
  {
    localIndex const li = iter->first;
    if( clearIfUnmapped )
    {
      relation[li].clear();
    }
    else
    {
      SortedArray< globalIndex > const & globalIndices = iter->second;
      for( auto const newGlobalIndex : globalIndices )
      {
        // NOTE: This simply ignores if newGlobalIndex is not found. This is OK if this function is
        // used for an upmap and the object shouldn't exist on this rank. There should be a better
        // way to check this.
        auto iterG2L = globalToLocal.find( newGlobalIndex );
        if( iterG2L != globalToLocal.end() )
        {
          relation[li].insert( iterG2L->second );
        }
      }
    }
  }
  unmappedIndices.clear();
}

inline
void ObjectManagerBase::fixUpDownMaps( ArrayOfSets< localIndex > & relation,
                                       unordered_map< globalIndex, localIndex > const & globalToLocal,
                                       map< localIndex, SortedArray< globalIndex > > & unmappedIndices,
                                       bool const clearIfUnmapped )
{
  GEOSX_MARK_FUNCTION;

  for( map< localIndex, SortedArray< globalIndex > >::iterator iter = unmappedIndices.begin();
       iter != unmappedIndices.end();
       ++iter )
  {
    localIndex const li = iter->first;
    if( clearIfUnmapped )
    {
      relation.clearSet( li );
    }
    else
    {
      SortedArray< globalIndex > const & globalIndices = iter->second;
      for( globalIndex const newGlobalIndex : globalIndices )
      {
        // NOTE: This simply ignores if newGlobalIndex is not found. This is OK if this function is
        // used for an upmap and the object shouldn't exist on this rank. There should be a better
        // way to check this.
        auto iterG2L = globalToLocal.find( newGlobalIndex );
        if( iterG2L != globalToLocal.end() )
        {
          relation.insertIntoSet( li, iterG2L->second );
        }
      }
    }
  }
  unmappedIndices.clear();
}

} /* namespace geosx */


/**
 * @brief Alias to ObjectManagerBase
 */
typedef geosx::ObjectManagerBase ObjectDataStructureBaseT;

#endif /* GEOSX_MESH_OBJECTMANAGERBASE_HPP_ */
