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

#ifndef GEOSX_MESH_WELLELEMENTSUBREGION_HPP_
#define GEOSX_MESH_WELLELEMENTSUBREGION_HPP_

#include "mesh/ElementSubRegionBase.hpp"
#include "mesh/InterObjectRelation.hpp"
#include "mesh/PerforationData.hpp"

namespace geosx
{

/**
 * @class WellElementSubRegion
 * @brief This class describes a collection of local well elements and perforations.
 */
class WellElementSubRegion : public ElementSubRegionBase
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = FixedOneToManyRelation; // unused but needed in MeshLevel::generateAdjacencyLists
  /// Alias for the type of the element-to-face map
  using FaceMapType = FixedOneToManyRelation; // unused but needed in MeshLevel::generateAdjacencyLists

  /**
   * @brief enumeration for values in segmentStatusList parameter of Generate()
   */
  enum WellElemStatus : unsigned
  {
    UNOWNED = 0,             // there are no perforations on this element
    REMOTE = 1,              // all perforations are remote
    LOCAL  = 2,              // all perforations are local
    SHARED = REMOTE | LOCAL  // both remote and local perforations
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  WellElementSubRegion( string const & name,
                        Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~WellElementSubRegion() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static const string catalogName() { return "wellElementSubRegion"; }

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override { return WellElementSubRegion::catalogName(); }

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void calculateElementGeometricQuantities( NodeManager const &,
                                                    FaceManager const & ) override
  {}

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the element-to-edge map.
   * @return a reference to element-to-edge map
   */
  FixedOneToManyRelation & edgeList()
  {
    return m_toEdgesRelation;
  }

  /**
   * @copydoc edgeList()
   */
  FixedOneToManyRelation const & edgeList() const
  {
    return m_toEdgesRelation;
  }

  /**
   * @brief Get the element-to-face map.
   * @return a reference to the element to face map
   */
  FixedOneToManyRelation & faceList()
  {
    return m_toFacesRelation;
  }

  /**
   * @copydoc faceList()
   */
  FixedOneToManyRelation const & faceList() const
  {
    return m_toFacesRelation;
  }

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  NodeMapType & nodeList()
  {
    return m_toNodesRelation;
  }

  /**
   * @copydoc nodeList()
   */
  NodeMapType const & nodeList() const
  {
    return m_toNodesRelation;
  }

  /**
   * @brief Get for the top element index.
   * @return local index of well's top element or -1 if it is not on current rank
   */
  localIndex getTopWellElementIndex() const
  {
    return m_topWellElementIndex;
  }

  /**
   * @brief Set the name of the WellControls object of this well.
   * @param[in] name the name of the WellControls object
   */
  void setWellControlsName( string const & name )
  {
    m_wellControlsName = name;
  }

  /**
   * @brief Get the name of the WellControls object of this well.
   * @return a string containing the name of the WellControls object
   */
  string const & getWellControlsName() const
  {
    return m_wellControlsName;
  }

  /**
   * @brief Get all the local perforations.
   * @return a pointer to the PerforationData object
   */
  PerforationData * getPerforationData()
  {
    return &m_perforationData;
  }

  /**
   * @copydoc getPerforationData()
   */
  PerforationData const * getPerforationData() const
  {
    return &m_perforationData;
  }

  /**
   * @brief Set for the MPI rank that owns this well (i.e. the top segment).
   * @param[in] rank MPI rank of the owner process
   */
  void setTopRank( int rank )
  {
    m_topRank = rank;
  }

  /**
   * @brief Check if well is owned by current rank
   * @return true if the well is owned by current rank, false otherwise
   */
  bool isLocallyOwned() const;

  ///@}

  /**
   * @name Construction of the well connectivity
   */
  ///@{

  /**
   * @brief Build the local well elements from global well element data.
   * @param[in] mesh the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] elemStatus list of well element status, as determined by perforations connected
   *                       to local or remote mesh partitions. Status values are defined in
   *                       enum SegmentStatus. They are used to partition well elements.
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void generate( MeshLevel & mesh,
                 InternalWellGenerator const & wellGeometry,
                 arrayView1d< integer > & elemStatus,
                 globalIndex nodeOffsetGlobal,
                 globalIndex elemOffsetGlobal );

  /**
   * @brief For each perforation, find the reservoir element that contains the perforation.
   * @param[in] mesh the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   */
  void connectPerforationsToMeshElements( MeshLevel & mesh,
                                          InternalWellGenerator const & wellGeometry );

  /**
   * @brief Reconstruct the (local) map nextWellElemId using nextWellElemIdGlobal after the ghost exchange.
   */
  void reconstructLocalConnectivity();

  ///@}

  /**
   * @name Overriding packing/unpacking functions
   */
  ///@{

  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpacks the specific elements in the @ packList.
   * @param[in] buffer The buffer containing the packed data.
   * @param[in] packList The (un)packed element.
   * @param[in] overwriteUpMaps Clear the up maps provided.
   * @param[in] overwriteDownMaps Clear the down maps provided.
   * @return The packed size.
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  virtual void fixUpDownMaps( bool const clearIfUnmapped ) final override;

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementSubRegionBase::viewKeyStruct
  {
    /// @return String key for the well control name
    static constexpr char const * wellControlsString() { return "wellControlsName"; }
    /// @return String key for the well element-to-node list
    static constexpr char const * wellNodeListString() { return "nodeList"; }
    /// @return String key for the local indices of the next well element (used in solvers)
    static constexpr char const * nextWellElementIndexString() { return "nextWellElementIndex"; }
    /// @return String key for the global indices of the next well element (to reconstruct maps)
    static constexpr char const * nextWellElementIndexGlobalString() { return "nextWellElementIndexGlobal"; }
    /// @return String key for the top well element index
    static constexpr char const * topWellElementIndexString() { return "topWellElementIndex"; }
    /// @return String key for the rank owning the top element
    static constexpr char const * topRankString() { return "topRank"; }
    /// @return String key for the well radius
    static constexpr char const * radiusString() { return "radius"; }

    /// ViewKey for the well control name
    dataRepository::ViewKey wellControlsName     = { wellControlsString() };
    /// ViewKey for the well element-to-node list
    dataRepository::ViewKey wellNodeList         = { wellNodeListString() };
    /// ViewKey for the local indices of the next well element (used in solvers)
    dataRepository::ViewKey nextWellElementIndex = { nextWellElementIndexString() };
    /// ViewKey for the global indices of the next well element (to reconstruct maps)
    dataRepository::ViewKey nextWellElementIndexGlobal = { nextWellElementIndexGlobalString() };
    /// ViewKey for the top well element index
    dataRepository::ViewKey topWellElementIndex = { topWellElementIndexString() };
    /// ViewKey for the rank owning the top element
    dataRepository::ViewKey topRank            = { topRankString() };
    /// ViewKey for the well radius
    dataRepository::ViewKey radius             = { radiusString() };
  }
  /// ViewKey struct for the WellElementSubRegion class
  viewKeysWellElementSubRegion;

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ElementSubRegionBase::groupKeyStruct
  {
    /// @return String key for the PerforationData object
    static constexpr char const * perforationDataString() { return "wellElementSubRegion"; }

    /// GroupKey for the PerforationData object
    dataRepository::GroupKey perforationData = { perforationDataString() };

  }
  /// groupKey struct for the WellElementSubRegion class
  groupKeysWellElementSubRegion;


private:

  /**
   * @brief Assign the unowned well elements (= well elem without perforation ) that are
            in the reservoir (and that can therefore be matched with a reservoir element) to an MPI rank.
   * @param[in] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] unownedElems set of unowned well elems.
   * @param[out] localElems set of local well elems. It contains the perforated well elements
                            connected to local mesh elements before the call, and is filled
                            with the newly assigned well elements in this function.
   * @param[out] wellElemStatus list of current well element status. Status values are defined in
   *                            enum SegmentStatus. They are used to partition well elements.
   */
  void assignUnownedElementsInReservoir( MeshLevel & mesh,
                                         InternalWellGenerator const & wellGeometry,
                                         SortedArray< globalIndex >           const & unownedElems,
                                         SortedArray< globalIndex > & localElems,
                                         arrayView1d< integer > & elemStatusGlobal ) const;

  /**
   * @brief Check that all the well elements have been assigned to a single rank.
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[out] localElems set of local well elems.
   * @param[out] wellElemStatus list of current well element status. Status values are defined in
   *                            enum SegmentStatus. They are used to partition well elements.
   *
   * This function also checks that if two ranks are neighbors in the well, they are also neighbors in the mesh.
   */
  void checkPartitioningValidity( InternalWellGenerator const & wellGeometry,
                                  SortedArray< globalIndex > & localElems,
                                  arrayView1d< integer > & elemStatusGlobal ) const;

  /**
   * @brief Add the well nodes to the nodeManager (properly resized).
   * @param[inout] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] localNodes set of local well nodes (includes boundary nodes). At this point all the nodes have been
   * collected
   * @param[in] boundaryNodes set of local well nodes that are at the boundary between this rank and another rank
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   *
   * The function WellElementSubRegion::CollectLocalAndBoundaryNodes must have been called before this function.
   */
  void updateNodeManagerSize( MeshLevel & mesh,
                              InternalWellGenerator const & wellGeometry,
                              SortedArray< globalIndex > const & localNodes,
                              SortedArray< globalIndex > const & boundaryNodes,
                              globalIndex nodeOffsetGlobal );

  /**
   * @brief Construct the subregion's local to global maps, as well as other local maps (toNodes, nextWellElemId,
   *        volume, etc).
   * @param[inout] meshLevel the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] localElems set of local well elems. At this point all the well elems have been assigned
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   *
   * The function WellElementSubRegion::UpdateNodeManagerSize must have been called before this function
   */
  void constructSubRegionLocalElementMaps( MeshLevel & mesh,
                                           InternalWellGenerator const & wellGeometry,
                                           SortedArray< globalIndex > const & localElems,
                                           globalIndex nodeOffsetGlobal,
                                           globalIndex elemOffsetGlobal );

  /**
   * @brief Constructs the toElementRegionList, toElementSubRegion, toElement maps
   * @param[inout] meshLevel the mesh object (single level only)
   *
   * This function is the equivalent of NodeManager::SetElementMaps for well elements.
   * The function WellElementSubRegion::ConstructSubRegionLocalElementMaps must have been called before this function
   */
  void updateNodeManagerNodeToElementMap( MeshLevel & mesh );

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

  /// Map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// Name of the WellControls object for this well
  string m_wellControlsName;

  /// Element-to-node relation is one to one relation.
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;  // unused but needed in MeshLevel::generateAdjacencyLists

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;  // unused but needed in MeshLevel::generateAdjacencyLists

  /// Local indices of the next well element (used in solvers)
  array1d< localIndex > m_nextWellElementIndex;

  /// Indices of the next well element (to reconstruct connectivity after ghost exchange)
  array1d< localIndex > m_nextWellElementIndexGlobal;

  /// Local index of well's top segment
  localIndex m_topWellElementIndex;

  /// Perforations
  PerforationData m_perforationData;

  /// Top rank
  integer m_topRank;

  /// Radius of the well element
  array1d< real64 > m_radius;

  /// Depth of the local search to match perforation to reservoir elements
  localIndex m_searchDepth;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_WELLELEMENTSUBREGION_HPP_ */
