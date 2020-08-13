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
 * @file MeshLevel.hpp
 */

#ifndef GEOSX_MESH_MESHLEVEL_HPP_
#define GEOSX_MESH_MESHLEVEL_HPP_

#include "NodeManager.hpp"
#include "EdgeManager.hpp"
#include "ElementRegionManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
class ElementRegionManager;

/**
 * @class MeshLevel
 * @brief Class facilitating the representation of a multi-level discretization of a MeshBody.
 * @details This contains the main components that compose a discretized mesh in GEOSX (nodes, faces, elements).
 *          In current practice, the code utilizes a single ``MeshLevel`` until such time as we
 *          implement a proper multi-level mesh capability.
 */
class MeshLevel : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the MeshLevel object.
   * @param[in] name the name of the MeshLevel object in the repository
   * @param[in] parent the parent group of the MeshLevel object being constructed
   */
  MeshLevel( string const & name, Group * const parent );
  virtual ~MeshLevel() override;

  /**
   * @brief Collects the nodes, edges, faces, and elements that are adjacent to a given list of nodes.
   * @param[in] seedNodeList the input nodes
   * @param[out] nodeAdjacencyList the nodes adjacent to the input nodes of seedNodeList
   * @param[out] edgeAdjacencyList the edges adjacent to the input nodes of seedNodeList
   * @param[out] faceAdjacencyList the faces adjacent to the input nodes of seedNodeList
   * @param[out] elementAdjacencyList the elements adjacent to the input nodes of seedNodeList
   * @param[in] depth the depth of the search for adjacent quantities (first-order neighbors, neighbors of neighbors, etc)
   */
  void
  GenerateAdjacencyLists(
    arrayView1d< localIndex const > const & seedNodeList,
    localIndex_array & nodeAdjacencyList,
    localIndex_array & edgeAdjacencyList,
    localIndex_array & faceAdjacencyList,
    ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > &
      elementAdjacencyList,
    integer const depth );

  virtual void
  InitializePostInitialConditions_PostSubGroups( Group * const ) override;

  /// @cond DO_NOT_DOCUMENT

  struct viewStructKeys
  {
    dataRepository::ViewKey meshLevel = { "meshLevel" };
  } viewKeys;

  struct groupStructKeys
  {
    dataRepository::GroupKey vertexManager = { "vertexManager" };
    dataRepository::GroupKey cellManager = { "cellManager" };

    static constexpr auto nodeManagerString = "nodeManager";
    static constexpr auto edgeManagerString = "edgeManager";
    static constexpr auto faceManagerString = "faceManager";

    // This key is defined in problem manager:
    static constexpr auto elemManagerString = "ElementRegions";

    dataRepository::GroupKey nodeManager = { nodeManagerString };
    dataRepository::GroupKey edgeManager = { edgeManagerString };
    dataRepository::GroupKey faceManager = { faceManagerString };
    dataRepository::GroupKey elemManager = { elemManagerString };
  } groupKeys;

  /// @endcond

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the node manager.
   * @return a pointer to the nodeManager object
   */
  NodeManager const *
  getNodeManager() const
  {
    return &m_nodeManager;
  }
  /**
   * @copydoc getNodeManager() const
   */
  NodeManager *
  getNodeManager()
  {
    return &m_nodeManager;
  }

  /**
   * @brief Get the edge manager.
   * @return a pointer to the edgeManager object
   */
  EdgeManager const *
  getEdgeManager() const
  {
    return &m_edgeManager;
  }
  /**
   * @copydoc getEdgeManager() const
   */
  EdgeManager *
  getEdgeManager()
  {
    return &m_edgeManager;
  }

  /**
   * @brief Get the face manager.
   * @return a pointer to the faceManager object
   */
  FaceManager const *
  getFaceManager() const
  {
    return &m_faceManager;
  }
  /**
   * @copydoc getFaceManager() const
   */
  FaceManager *
  getFaceManager()
  {
    return &m_faceManager;
  }

  /**
   * @brief Get the element region manager.
   * @return a pointer to the elementRegionManager object
   */
  ElementRegionManager const *
  getElemManager() const
  {
    return &m_elementManager;
  }
  /**
   * @copydoc getElemManager() const
   */
  ElementRegionManager *
  getElemManager()
  {
    return &m_elementManager;
  }

  ///@}

private:
  /// Manager for node data
  NodeManager m_nodeManager;
  /// Manager for edge data
  EdgeManager m_edgeManager;
  /// Manager for face data
  FaceManager m_faceManager;
  /// Manager for element data
  ElementRegionManager m_elementManager;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_MESHLEVEL_HPP_ */
