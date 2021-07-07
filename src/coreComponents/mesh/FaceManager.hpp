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
 * @file FaceManager.hpp
 */

#ifndef GEOSX_MESH_FACEMANAGER_HPP_
#define GEOSX_MESH_FACEMANAGER_HPP_

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "ToElementRelation.hpp"

namespace geosx
{

class NodeManager;
class ElementRegionManager;
class CellElementSubRegion;

/**
 * @class FaceManager
 * @brief The FaceManager class provides an interface to ObjectManagerBase in order to manage face data.
 *
 * The FaceManager class manages the face data using face indexed or keyed containers.
 * This means that each field is stored in a way where each container entry
 * corresponds to a face.
 */
class FaceManager : public ObjectManagerBase
{
public:

  /// FaceToNode map type
  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /// FaceToEdge map type
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /// FaceToElement map type
  using ElemMapType = FixedToManyElementRelation;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Return the name of the FaceManager in the object catalog.
   * @return string that contains the catalog name of the FaceManager
   */
  static const string catalogName()
  { return "FaceManager"; }

  /**
   * @brief Provide a virtual access to catalogName().
   * @return string that contains the catalog name of the FaceManager
   */
  virtual const string getCatalogName() const override
  { return FaceManager::catalogName(); }
  ///@}

  /**
   * @brief Get the default number of node per face in node list
   * @return the default number of node per face in node list
   *
   * @note Value forwarding is due to refactoring.
   */
  static constexpr localIndex nodeMapExtraSpacePerFace()
  { return CellBlockManagerABC::nodeMapExtraSpacePerFace(); }

  /**
   * @brief Get the default number of edge per face in edge list
   * @return the default number of edge per face in edge list
   *
   * @note Value forwarding is due to refactoring.
   */
  static constexpr localIndex edgeMapExtraSpacePerFace()
  { return CellBlockManagerABC::edgeMapExtraSpacePerFace(); }

  /**
   * @name Constructors/destructor
   */
  ///@{

  /**
   * @brief Main Constructor for FaceManager
   * @param[in] name the name of FaceManager
   * @param[in] parent the parent Group object of FaceManager
   */
  FaceManager( string const & name, Group * const parent );

  /**
   * @brief Destructor override from ObjectManager
   */
  virtual ~FaceManager() override;

  /**
   * @brief Deleted default constructor
   */
  FaceManager() = delete;

  /**
   * @brief Deleted copy constructor
   */
  FaceManager( FaceManager const & ) = delete;


  /**
   * @brief Deleted move constructor
   */
  FaceManager( FaceManager && ) = delete;

  ///@}

  /**
   * @brief Extend base class resize method resizing  m_nodeList, m_edgeList member containers.
   * @details the \p newSize of this FaceManager is the number of faces it will contain
   * @param[in] newsize new size the FaceManager.
   */
  virtual void resize( localIndex const newsize ) override;

  /**
   * @brief Build faces in filling face-to-node and face-to-element mappings.
   * @param[in] cellBlockManager cell block manager
   * @param[in] nodeManager mesh node manager
   * @param[in] elementManager element manager
   */
  void buildFaces( CellBlockManagerABC const & cellBlockManager,
                   NodeManager & nodeManager,
                   ElementRegionManager & elementManager );

  /**
   * @brief Compute faces center, area and normal.
   * @param[in] nodeManager NodeManager associated with the current DomainPartition
   */
  void computeGeometry( NodeManager const & nodeManager );

  /**
   * @brief Return the number of nodes of the faces with the greatest number of nodes.
   * @return the maximum number of nodes a face have
   */
  localIndex getMaxFaceNodes() const;

  /**
   * @brief Sort all faces nodes counterclockwisely in m_nodeList.
   * @param[in] nodeManager node manager allowing access to face nodes coordinates
   * @param[in] elemManager element manager allowing access to the cell elements
   */
  void sortAllFaceNodes( NodeManager const & nodeManager,
                         ElementRegionManager const & elemManager );

  /**
   * @brief Reorder face nodes to be labeled counter-clockwise resulting in outgoing normal.
   * @param[in] X array view of mesh nodes coordinates
   * @param[in] elementCenter coordinate of the element center
   * @param[in,out] faceNodes reordered local label list of nodes
   * @param[in] numFaceNodes number of nodes for the face
   */
  void sortFaceNodes( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                      arraySlice1d< real64 const > const elementCenter,
                      localIndex * const faceNodes,
                      localIndex const numFaceNodes );

  /**
   * @brief Flag face and nodes'face with at least one element on the boundary.
   * @param[in] nodeManager manager of mesh nodes
   */
  void setDomainBoundaryObjects( NodeManager & nodeManager );

  /**
   * @brief Flag faces on boundary or external to the DomainPartition.
   */
  void setIsExternal();

  /**
   * @name Packing methods
   */

  ///@{

  /**
   * @brief Create an array listing all excluded local face indices values.
   * @param [in,out] exclusionList Sorted array with excluded local faces indices
   */
  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  /**
   * @brief Calculate the size that a list would have if it were packed, but without actually packing it.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in] packList the list of node indices that we wish to get the size of after packing
   * @return the size of packList if it were packed
   * @note This function does not perform any packing, it just evaluates and returns the possible packed size.
   */
  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Pack an array of node indices into a buffer.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in,out] buffer buffer to pack the node index data into
   * @param[in] packList the indices of nodes that should be packed
   * @return a localIndex value representing the size of the packed data
   */
  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpack a buffer to an array of node indices.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in] buffer buffer with the packed data
   * @param [inout] packList an array of localIndex values that we wish to unpack to
   * @param [in] overwriteUpMaps boolean: true to overwrite the previous Up maps
   * @param [in] overwriteDownMaps boolean: true to overwrite the previous Down maps
   * @return a localIndex value representing the size of the unpacked list
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Call fixUpDownMaps for nodes-to-edges and nodes-to-faces maps.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param [in] clearIfUnmapped boolean: true to remove if it is not mapped
   */
  void fixUpDownMaps( bool const clearIfUnmapped );

  ///@}

  /**
   * @brief Compress FaceManager face-to-node and face-to-edge containers so that the values of
   * each array are contiguous with no extra capacity in between.
   * @note The method used here on each arrays (compress) does not free any memory.
   */
  void compressRelationMaps();

  /**
   * @brief Enforce child faces and parent faces to have opposite normals.
   * @param[in] targetIndices set of face indices for which the enforcement has to be done
   */
  virtual void enforceStateFieldConsistencyPostTopologyChange( std::set< localIndex > const & targetIndices ) override;

  /**
   * @brief Clean up the mappings from faces to element index, region, subregion on a new (updated) list of faces, in order to keep only
   * relevant mappings.
   * @param [in] receivedFaces the new list of target node indices
   * @param [in] elemRegionManager Element Region Manager
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedFaces,
                         ElementRegionManager const & elemRegionManager );


  //void SetGlobalIndexFromCompositionalObject( ObjectManagerBase const * const compositionalObject );

  /**
   * @brief Extract a face-to-nodes map with global indexed for boundary faces.
   * @param[in] nodeManager mesh nodeManager
   * @param[out] faceToNodes face-to-node map
   */
  virtual void extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                                std::vector< std::vector< globalIndex > > & faceToNodes ) override;

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   *  @struct viewKeysStruct
   *  @brief struct containing the view access keys to be bound with class data member
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @cond DO_NOT_DOCUMENT
    static constexpr char const * nodeListString() { return "nodeList"; }
    static constexpr char const * edgeListString() { return "edgeList"; }
    static constexpr char const * elementRegionListString() { return "elemRegionList"; }
    static constexpr char const * elementSubRegionListString() { return "elemSubRegionList"; }
    static constexpr char const * elementListString() { return "elemList"; }
    constexpr static char const * faceAreaString() { return "faceArea"; }
    constexpr static char const * faceCenterString() { return "faceCenter"; }
    constexpr static char const * faceNormalString() { return "faceNormal"; }

    dataRepository::ViewKey nodeList              = { nodeListString() };
    dataRepository::ViewKey edgeList              = { edgeListString() };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString() };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString() };
    dataRepository::ViewKey elementList           = { elementListString() };
    /// @endcond
  }
  /// viewKeys
  viewKeys;

  ///@}

  /**
   * @name Accessors for FaceManager fixed data
   */
  ///@{

  /**
   * @brief Get a mutable accessor to a table containing all the face area.
   * @details this table is mutable so it can be used to compute
   * or modify the face area in this FaceManager
   * @return a table containing all the face area
   */
  array1d< real64 > & faceArea()       { return m_faceArea; }

  /**
   * @brief Get an immutable accessor to a table containing all the face area.
   * @return an immutable table containing all the face area
   */
  arrayView1d< real64 const > faceArea() const { return m_faceArea; }

  /**
   * @brief Get an immutable accessor to a table containing all the face centers.
   * @return an immutable table containing all the face centers
   */
  arrayView2d< real64 const > faceCenter() const { return m_faceCenter; }

  /**
   * @brief Get a mutable accessor to a table containing all the face normals.
   * @details this table is mutable so it can be used to compute
   * or modify the face normals in this FaceManager
   * @return a table containing all the face normals
   */
  array2d< real64 > const & faceNormal() { return m_faceNormal; }

  /**
   * @brief Get an immutable accessor to a table containing all the face normals.
   * @return an immutable table containing all the face normals
   */
  arrayView2d< real64 const > faceNormal() const { return m_faceNormal; }

  /**
   * @brief Get a mutable accessor to a map containing the list of each nodes for each faces.
   * @return non-const reference to a map containing the list of each nodes for each faces
   */
  NodeMapType & nodeList()                    { return m_nodeList; }

  /**
   * @brief Get an immutable accessor to a map containing the list of each nodes for each faces
   * @return const reference to a map containing the list of each nodes for each faces
   */
  NodeMapType const & nodeList() const { return m_nodeList; }

  /**
   * @brief Get a mutable accessor to a map containing the list of each edges for each faces
   * @return non-const reference to a map containing the list of each edges for each faces
   */
  EdgeMapType & edgeList()       { return m_edgeList; }

  /**
   * @brief Get an immutable accessor to a map containing the list of each edges for each faces.
   * @return const reference to a map containing the list of each edges for each faces
   */
  EdgeMapType const & edgeList() const { return m_edgeList; }

  /**
   * @brief Get a mutable accessor to the faces-to-ElementRegion relation.
   * @return non-const reference to faces-to-ElementRegion relation
   */
  array2d< localIndex > const & elementRegionList() { return m_toElements.m_toElementRegion; }

  /**
   * @brief Get an immutable accessor to the faces-to-ElementRegion relation.
   * @return const reference to nodes-to-ElementRegion relation
   */
  arrayView2d< localIndex const > elementRegionList() const { return m_toElements.m_toElementRegion; }

  /**
   * @brief Get a mutable accessor to the faces-to-ElementSubRegion relation.
   * @return non-const reference to faces-to-ElementSubRegion relation
   */
  array2d< localIndex > const & elementSubRegionList() { return m_toElements.m_toElementSubRegion; }

  /**
   * @brief Get an immutable accessor to the faces-to-ElementSubRegion relation.
   * @return const reference to faces-to-ElementSubRegion relation
   */
  arrayView2d< localIndex const > elementSubRegionList() const { return m_toElements.m_toElementSubRegion; }

  /**
   * @brief Get a mutable accessor to the faces-to-element-index relation.
   * @return non-const reference to faces-to-element-index relation
   */
  array2d< localIndex > const & elementList() { return m_toElements.m_toElementIndex; }

  /**
   * @brief Get an imutable accessor to the faces-to-element-index relation.
   * @return const reference to faces-to-elements relation
   */
  arrayView2d< localIndex const > elementList() const { return m_toElements.m_toElementIndex; }

  /**
   * @brief Get a mutable accessor to the faces-to-element-index relation.
   * @details The returned ElemMapType gives access, for one face
   * to the element index, the element sub region, and the element region.
   * @return non-const reference to faces-to-element-index relation
   */
  ElemMapType & toElementRelation()       { return m_toElements; }

  /**
   * @brief Get an immutable accessor to the faces-to-element-index relation.
   * @details The returned ElemMapType gives access, for one face
   * to the element index, the element sub region, and the element region.
   * @return non-const reference to faces-to-element-index relation
   */
  ElemMapType const & toElementRelation() const { return m_toElements; }

  ///}@

private:

  /**
   * @brief Pack the upward and downward pointing maps.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occurred if set to TRUE.
   * @param[inout] buffer the buffer to pack data into
   * @param[in] packList the indices of faces that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

  /// face keyed map containing face-to-node relation
  NodeMapType m_nodeList;

  /// face keyed map containing face-to-edge relation
  EdgeMapType m_edgeList;

  /// face keyed map containing face-to-element relation
  ElemMapType m_toElements;

  /// map of global to local indices for nodes
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// map of global  to local  indices for edges
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToEdges;

  /// list of faces area
  array1d< real64 > m_faceArea;

  /// list of faces center
  array2d< real64 > m_faceCenter;
  /// list of faces normal
  array2d< real64 > m_faceNormal;

  /// constant expression of the maximum number of nodes per faces
  constexpr static int MAX_FACE_NODES = 9;

};

}
#endif /* FACEMANAGERT_H_ */
