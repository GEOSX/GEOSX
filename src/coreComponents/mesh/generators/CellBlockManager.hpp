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
 * @file CellBlockManager.hpp
 */

#ifndef GEOSX_MESH_CELLBLOCKMANAGER_H_
#define GEOSX_MESH_CELLBLOCKMANAGER_H_

#include "mesh/generators/CellBlockManagerABC.hpp"
#include "mesh/generators/CellBlock.hpp"

namespace geosx
{

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public CellBlockManagerABC
{
public:

  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  CellBlockManager( const CellBlockManager & ) = delete;

  CellBlockManager & operator=( const CellBlockManager & ) = delete;

  ~CellBlockManager() override = default;

  virtual Group * createChild( string const & childKey, string const & childName ) override;

  /**
   * @brief Maximum number of faces allowed (in memory) per each node.
   * @return The number as an integer.
   */
  static constexpr int maxFacesPerNode()
  { return 200; }

  /**
   * @brief Extra space for node to faces mapping.
   * @return Number of extra values as an integer.
   */
  static constexpr localIndex getFaceMapOverallocation()
  { return 8; }

  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & getNodesPositions() const override;

  /**
   * @brief Returns a mutable reference to vector holding the nodes coordinates
   * @return The reference
   *
   * @note This is meant to be used as a setter. Do not resize the vector yourself.
   */
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & getNodesPositions();

  ArrayOfSets< localIndex > getNodeToEdges() const override;

  ArrayOfSets< localIndex > getNodeToFaces() const override;

  ArrayOfArrays< localIndex > getNodeToElements() const override;

  array2d< geosx::localIndex > const & getEdgeToNodes() const override;

  ArrayOfSets< geosx::localIndex > const & getEdgeToFaces() const override;

  ArrayOfArrays< localIndex > getFaceToNodes() const override;

  ArrayOfArrays< geosx::localIndex > const & getFaceToEdges() const override;

  array2d< localIndex > getFaceToElements() const override;

  array1d< globalIndex > getNodeLocalToGlobal() const override;

  /**
   * @brief Returns a mutable reference to vector holding the node to global mapping.
   * @return The reference
   *
   * @note This is meant to be used as a setter. Do not resize the vector yourself.
   */
  array1d< globalIndex > & getNodeLocalToGlobal();

  const std::map< string, SortedArray< localIndex > > & getNodeSets() const override;

  /**
   * @brief Returns a mutable reference to the node sets.
   * @return A reference to the mapping.
   *
   * The key of the map is the name of the set.
   * While the values are sorted arrays which sizes are meant to be managed by the client code.
   * This member function is meant to be used like a setter.
   */
  std::map< string, SortedArray< localIndex > > & getNodeSets();

  /**
   * @brief Defines the number of nodes and resizes some underlying arrays appropriately.
   * @param[in] numNodes The number of nodes.
   *
   * The nodes coordinates and nodes local to global mappings get resized to @p numNodes.
   */
  void setNumNodes( localIndex numNodes ); // TODO Improve doc. Is it per domain, are there duplicated nodes because of subregions?

  localIndex numNodes() const override;

  localIndex numEdges() const override;

  localIndex numFaces() const override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param numElements list of the new element numbers
   * @param regionNames list of the element region names
   * @param elementTypes list of the element types
   */
  void resize( integer_array const & numElements,
               string_array const & regionNames,
               string_array const & elementTypes );

  /**
   * @brief Trigger the computation of all the mappings.
   *
   * Call this member function to compute all the mappings.
   * Computations could be done lazily when calling getters.
   * But this is not yet implemented.
   */
  void buildMaps();

  /**
   * @brief Get cell block by name.
   * @param name Name of the cell block.
   * @return Reference to the cell block instance.
   */
  CellBlock & getCellBlock( string const & name )
  {
    return this->getGroup( viewKeyStruct::cellBlocks() ).getGroup< CellBlock >( name );
  }

  Group & getCellBlocks() override;

  /**
   * @brief Registers and returns a cell block of name @p name.
   * @param name The name of the created cell block.
   * @return A reference to the new cell block. The CellBlockManager owns this new instance.
   */
  CellBlock & registerCellBlock( string name );

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void forElementSubRegions( LAMBDA lambda )
  {
    this->getGroup( viewKeyStruct::cellBlocks() ).forSubGroups< CellBlock >( lambda );
  }

private:

  struct viewKeyStruct
  {
    /// Cell blocks key
    static constexpr char const * cellBlocks() { return "cellBlocks"; }
  };

  /**
   * @brief Returns a group containing the cell blocks as CellBlockABC instances
   * @return Reference to the Group instance.
   */
  const Group & getCellBlocks() const;

  /**
   * @brief Get cell block at index @p iCellBlock.
   * @param iCellBlock The cell block index.
   * @return Const reference to the instance.
   *
   * @note Mainly useful for iteration purposes.
   */
  const CellBlockABC & getCellBlock( localIndex iCellBlock ) const;

  /**
   * @brief Returns the number of cells blocks
   * @return Number of cell blocks
   */
  localIndex numCellBlocks() const;

  /**
   * @brief Trigger the node to edges mapping computation.
   */
  void buildNodeToEdges();

  /**
   * @brief Trigger the face to nodes, edges and elements mappings.
   */
  void buildFaceMaps();

  array2d< real64, nodes::REFERENCE_POSITION_PERM > m_nodesPositions;

  ArrayOfSets< localIndex > m_nodeToEdges;
  ArrayOfSets< localIndex > m_edgeToFaces;
  array2d< localIndex > m_edgeToNodes;
  ArrayOfArrays< localIndex >  m_faceToNodes;
  ArrayOfArrays< localIndex > m_faceToEdges;
  array2d< localIndex >  m_faceToElements;

  array1d< globalIndex >  m_nodeLocalToGlobal;

  std::map< string, SortedArray< localIndex > > m_nodeSets;

  localIndex m_numNodes;
  localIndex m_numFaces;
  localIndex m_numEdges;
};

}
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
