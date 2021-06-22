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
 * @file EdgeManager.hpp
 */

#ifndef GEOSX_MESH_EDGEMANAGER_HPP_
#define GEOSX_MESH_EDGEMANAGER_HPP_

#include "InterObjectRelation.hpp"
#include "mesh/ObjectManagerBase.hpp"
#include "ToElementRelation.hpp"
#include "LvArray/src/tensorOps.hpp"


namespace geosx
{
class FaceManager;
class NodeManager;
class CellBlockManager;


/**
 * @class EdgeManager
 * @brief This class provides an interface to ObjectManagerBase in order to manage edge data.
 * @details An edge is bounded by two nodes. This class aims to manage all the edges of every
 * cell elements in a DomainPartition.
 */

class EdgeManager : public ObjectManagerBase
{
public:

  /// EdgeToNode map type
  using NodeMapType = FixedOneToManyRelation;

  /// EdgeToFace map type
  using FaceMapType = InterObjectRelation< ArrayOfSets< localIndex > >;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @return the string representing the edge manager name in the catalog
   */
  static const string catalogName()
  { return "EdgeManager"; }


  /**
   * @brief Getter used to access the edge manager catalog name.
   * @return the edge manager catalog name
   */
  virtual const string getCatalogName() const override
  { return EdgeManager::catalogName(); }

  ///@}

  /**
   * @brief Oversize the Face mapping by this amount for each edge (hardcoded)
   * @return extra space for each edge of the EdgetoFace map
   */
  static localIndex faceMapExtraSpacePerEdge()
  { return 4; }

  /**
   * @name Constructors/destructors
   */
  ///@{

  /**
   * @brief main constructor for EdgeManager Objects
   * @param[in] name the name of the EdgeManager object in the repository
   * @param[in] parent the parent group of the EdgeManager object being constructed
   */
  EdgeManager( string const & name,
               Group * const parent );

  /**
   * @brief default destructor
   */
  ~EdgeManager() override;

  ///@}

  /**
   * @brief Resize the EdgeManager object and all it members.
   * @details the size of the EdgeMananager is the number of edges
   * in the domain.
   * @param[in] newSize the new number of edges.
   */
  virtual void resize( localIndex const newSize ) override;

  /**
   * @brief Set the node of the domain boundary object.
   * @param[in] referenceObject the reference of the face manager.
   */
  void setDomainBoundaryObjects( FaceManager const & referenceObject );

  /**
   * @brief Set external edges.
   * @details external edges are the edges on the faces which are external
   * @param[in] faceManager the face manager to obtain external face
   */
  void setIsExternal( FaceManager const & faceManager );

  /**
   * @brief Build faces-to-edges and nodes-to-edges relation maps.
   * @param[in] nodeManager manager of all nodes in the DomainPartition
   * @param[in] faceManager manager of all faces in the DomainPartition
   */
  void buildEdges( NodeManager & nodeManager, FaceManager & faceManager );

  /**
   * @brief Build faces-to-edges and nodes-to-edges relation maps.
   * @param[in] numNodes number of nodes
   * @param[in] faceToNodeMap manager face-to-nodes map
   * @param[in] faceToEdgeMap manager face-to-edges map
   */
  void buildEdges( localIndex const numNodes,
                   ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                   ArrayOfArrays< localIndex > & faceToEdgeMap );


  /**
   * @brief Build \p globalEdgeNodes, a  vector containing all the global indices
   * of each nodes of each edges
   * @param[in] nodeManager the nodeManager object.
   * @param[out] globalEdgeNodes the globalEdgesNodes array to be built
   * [ [global_index_node_0_edge_0, global_index_node1_edge_0], [global_index_node_0_edge_1, global_index_node1_edge_1] ....]
   */
  virtual void
  extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                   std::vector< std::vector< globalIndex > > & globalEdgeNodes ) override;

  /**
   * @brief Compute the future size of a packed list.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in] packList the list of edge indices to be packed
   * @return The size of the packed list
   */
  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Pack an array of edge indices in a buffer.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in,out] buffer a buffer to pack the edge index data into
   * @param[in] packList the indices of edge to pack
   * @return The size of the packed array
   */
  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  /**
   * @brief Unpack a buffer to an array of indices that has been packed.
   * @details Packed data are meant to be communicated to other MPI ranks
   * @param[in] buffer buffer containing the packed indices
   * @param[in,out] packList the array of local index to be unpacked
   * @param[in] overwriteUpMaps boolean to state if the Up maps should be overwritten
   * @param[in] overwriteDownMaps boolean to state if the Downs maps should be overwritten
   * @return The size of the unpacked array
   */
  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       localIndex_array & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Call fixUpDownMaps of the class ObjectManagerBase for nodes-to-edges and nodes-to-faces relation maps.
   * @param[in] clearIfUnmapped boolean to indicate if the unmaped indices should be removed or not
   */
  void fixUpDownMaps( bool const clearIfUnmapped );

  /**
   * @brief Compress all nodes-to-faces relation maps
   * so that the values of each array are contiguous with no extra capacity in between.
   */
  void compressRelationMaps();

  /**
   * @brief Clean up the edges-to-faces mapping with respect to a new list of edges.
   * @param[in] receivedEdges the new list of edges indices
   * @param[in] facesToEdges map to go from faces to edges
   */
  void depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                         ArrayOfArraysView< localIndex const > const & facesToEdges );

  /**
   * @brief Build the mapping edge-to-nodes relation from the mapping global to local nodes.
   * @param[in] indices array of index of the global to local nodes
   * @param[in] nodeGlobalToLocal map of the global to local nodes
   * @param[in] faceGlobalToLocal GEOX UNUSED PARAMETER
   */
  void connectivityFromGlobalToLocal( const SortedArray< localIndex > & indices,
                                      const map< globalIndex, localIndex > & nodeGlobalToLocal,
                                      const map< globalIndex, localIndex > & faceGlobalToLocal );

  /**
   * @brief Split an edge (separate its two extremity nodes)
   * @param[in] indexToSplit Index of the edge to split
   * @param[in] parentNodeIndex index of the parent node
   * @param[in] childNodeIndex index of the child node
   * @param[in] nodesToEdges array of nodes-to-edges list
   */
  void splitEdge( const localIndex indexToSplit,
                  const localIndex parentNodeIndex,
                  const localIndex childNodeIndex[2],
                  array1d< SortedArray< localIndex > > & nodesToEdges );

  /**
   * @brief Check if edge \p edgeID contains node \p nodeID
   * @param[in] edgeID local index of the edge
   * @param[in] nodeID local index of the node
   * @return boolean : true if the node \p nodeID is part of the edge \p edgeID; false otherwise
   */
  bool hasNode( const localIndex edgeID, const localIndex nodeID ) const;

  /**
   * @brief Calculate the center of an edge given its index.
   * @tparam OUT_VECTOR type of output
   * @param edgeIndex index of the edge
   * @param X array view of nodes associated with the edge
   * @param edgeCenter output vector containing coordinate of the edge center
   */
  template< typename OUT_VECTOR >
  void calculateCenter( localIndex const edgeIndex,
                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                        OUT_VECTOR && edgeCenter ) const;

  /**
   * @brief Calculate the edge segment.
   * @tparam OUT_VECTOR type of output
   * @param edgeIndex index of the edge
   * @param X array view of the nodes associated with the NodeManager of current domain
   * @param edgeVector output vector containing edge segment
   */
  template< typename OUT_VECTOR >
  void calculateLength( localIndex const edgeIndex,
                        arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                        OUT_VECTOR && edgeVector ) const;

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{

  /**
   * @struct viewKeyStruct
   * @brief Container of keys needed to access the data of the class member
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// @cond DO_NOT_DOCUMENT
    static constexpr char const * nodeListString() { return "nodeList"; }
    static constexpr char const * faceListString() { return "faceList"; }
    static constexpr char const * elementRegionListString() { return "elemRegionList"; }
    static constexpr char const * elementSubRegionListString() { return "elemSubRegionList"; }
    static constexpr char const * elementListString() { return "elemList"; }
    static constexpr char const * edgesTofractureConnectorsEdgesString() { return "edgesToFractureConnectors"; }
    static constexpr char const * fractureConnectorEdgesToEdgesString() { return "fractureConnectorsToEdges"; }
    static constexpr char const * fractureConnectorsEdgesToFaceElementsIndexString() { return "fractureConnectorsToElementIndex"; }

    dataRepository::ViewKey nodesList             = { nodeListString() };
    dataRepository::ViewKey faceList              = { faceListString() };
    dataRepository::ViewKey elementRegionList     = { elementRegionListString() };
    dataRepository::ViewKey elementSubRegionList  = { elementSubRegionListString() };
    dataRepository::ViewKey elementList           = { elementListString() };
    /// @endcond
  }
  /// viewKeys
  viewKeys;

  ///}@

  /**
   * @brief Return the  maximum number of edges per node.
   * @return Maximum allowable number of edges connected to one node (hardcoded for now)
   */
  constexpr int maxEdgesPerNode() const { return 6; }

  /**
   * @name Getters for stored value.
   */
  ///@{

  /**
   * @brief Get a node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType & nodeList()       { return m_toNodesRelation; }

  /**
   * @brief Provide an accessor to an immutable node list.
   * @return reference to the list of edge to node relation
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

  /**
   * @brief Get a reference to a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex & nodeList( localIndex const edgeIndex, localIndex const nodeIndex )
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  /**
   * @brief Get a node list
   * @param[in] edgeIndex local index of the edge
   * @param[in] nodeIndex local index of the node
   * @return a reference of the index of the edge to node relation
   */
  localIndex nodeList( localIndex const edgeIndex, localIndex const nodeIndex ) const
  { return m_toNodesRelation( edgeIndex, nodeIndex ); }

  /**
   * @brief Get the edge-to-face relation
   * @return A FaceMapType
   */
  FaceMapType & faceList()       { return m_toFacesRelation; }

  /**
   * @brief Get an immutable accessor to the edge-to-face relation
   * @return An immutable FaceMapType
   */
  FaceMapType const & faceList() const { return m_toFacesRelation; }

  ///@}

  // TODO These should be in their own subset of edges when we add that capability.

  /// map from the edges to the fracture connectors index (edges that are fracture connectors)
  SortedArray< localIndex > m_recalculateFractureConnectorEdges;

  /// todo
  map< localIndex, localIndex > m_edgesToFractureConnectorsEdges;

  /// todo
  array1d< localIndex > m_fractureConnectorsEdgesToEdges;

  /// todo
  ArrayOfArrays< localIndex > m_fractureConnectorEdgesToFaceElements;


private:

  /// Map for the edges-to-nodes relation
  NodeMapType m_toNodesRelation;

  /// Map for the edges-to-faces relation
  FaceMapType m_toFacesRelation;

  /// Unmaped edge indices (those that are not in the nodes-to-edges relation map)
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// Unmaped edge indices (those that are not in the faces-to-edges relation map)
  map< localIndex, SortedArray< globalIndex > > m_unmappedGlobalIndicesInToFaces;

  /**
   * @brief function to pack the upward and downward pointing maps.
   * @tparam DOPACK template argument to determine whether or not to pack the buffer. If false, the buffer is not
   *                packed and the function returns the size of the packing that would have occurred if set to TRUE.
   * @param buffer the buffer to pack data into
   * @param packList the indices of nodes that should be packed.
   * @return size of data packed in terms of number of chars
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

};

template< typename OUT_VECTOR >
inline void EdgeManager::calculateCenter( localIndex const edgeIndex,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                          OUT_VECTOR && edgeCenter ) const
{
  LvArray::tensorOps::copy< 3 >( edgeCenter, X[m_toNodesRelation[edgeIndex][0]] );
  LvArray::tensorOps::add< 3 >( edgeCenter, X[m_toNodesRelation[edgeIndex][1]] );
  LvArray::tensorOps::scale< 3 >( edgeCenter, 0.5 );
}

template< typename OUT_VECTOR >
inline void EdgeManager::calculateLength( localIndex const edgeIndex,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                          OUT_VECTOR && edgeSegment ) const
{
  LvArray::tensorOps::copy< 3 >( edgeSegment, X[m_toNodesRelation[edgeIndex][1]] );
  LvArray::tensorOps::subtract< 3 >( edgeSegment, X[m_toNodesRelation[edgeIndex][0]] );
}

}
#endif /* GEOSX_MESH_EDGEMANAGER_HPP_ */
