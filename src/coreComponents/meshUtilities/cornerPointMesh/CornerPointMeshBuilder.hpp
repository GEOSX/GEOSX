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
 * @file CornerPointMeshBuilder.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHBUILDER_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHBUILDER_HPP_

#include "dataRepository/ObjectCatalog.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshParser.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshData.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshPartition.hpp"

namespace geosx
{

namespace cornerPointMesh
{

/**
 * @class CornerPointMeshBuilder
 * @brief This class is in charge of constructing a conforming mesh from an Eclipse corner-point mesh
 */
class CornerPointMeshBuilder
{

public:

  /**
   * @brief Constructor.
   * @param name the name of the class
   */
  CornerPointMeshBuilder( string const & name );

  /// Default virtual destructor
  virtual ~CornerPointMeshBuilder() = default;

  /// Default copy constructor
  CornerPointMeshBuilder( CornerPointMeshBuilder const & ) = default;

  /// Default move constructor
  CornerPointMeshBuilder( CornerPointMeshBuilder && ) = default;

  /// Deleted copy assignment operator
  CornerPointMeshBuilder & operator=( CornerPointMeshBuilder const & ) = delete;

  /// Deleted move assignment operator
  CornerPointMeshBuilder & operator=( CornerPointMeshBuilder && ) = delete;

/// using alias for templated Catalog CornerPointMeshBuilder type
  using CatalogInterface = dataRepository::CatalogInterface< CornerPointMeshBuilder,
                                                             string const & >;

  /**
   * @brief Getter for the singleton Catalog object
   * @return a static reference to the Catalog object
   *
   */
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief Define the catalog name for this class
   * @return the catalog name
   */
  static string catalogName() { return "CornerPointMeshBuilder"; }

  /**
   * @brief Const getter for catalog name of this class
   * @return the catalog name
   */
  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Driver function of the class: reads a corner-point mesh in the Eclipse format and constructs a conforming mesh
   * @param filePath the GRDECL file that contains the description of the corner-point mesh
   */
  void buildMesh( Path const & filePath );

  // vertices

  /**
   * @brief Const getter for the (unique) vertex positions
   * @return the 2D array of vertex positions
   */
  arrayView2d< real64 const > vertexPositions() const
  { return m_vertices.m_vertexPositions.toViewConst(); }

  /**
   * @brief Const getter for the map from corner-point vertex to (unique) vertex
   * @return the map from corner-point vertex to (unique) vertex
   * @note cp stands for corner point
   */
  arrayView1d< localIndex const > cpVertexToVertex() const
  { return m_vertices.m_cpVertexToVertex.toViewConst(); }

  /**
   * @brief Const getter for the map from local corner-point vertex index to global corner-point vertex index
   * @return the map from local corner-point vertex index to global corner-point vertex index
   * @note cp stands for corner point
   */
  arrayView1d< globalIndex const > cpVertexToGlobalCPVertex() const
  { return m_vertices.m_cpVertexToGlobalCPVertex.toViewConst(); }

  /**
   * @brief Const getter for the map from local vertex index to global vertex index
   * @return the map from local vertex index to global corner-point vertex index
   */
  arrayView1d< globalIndex const > vertexToGlobalVertex() const
  { return m_vertices.m_vertexToGlobalVertex.toViewConst(); }

  // cells

  /**
   * @brief Const getter for the map from active cell inside partition to active cell
   * @return the map from active cell inside partition to active cell
   */
  arrayView1d< localIndex const > ownedActiveCellToActiveCell() const
  { return m_cells.m_ownedActiveCellToActiveCell.toViewConst(); }

  /**
   * @brief Const getter for the map from active cell to cell
   * @return the map from active cell to cell
   */
  arrayView1d< localIndex const > activeCellToCell() const
  { return m_cells.m_activeCellToCell.toViewConst(); }

  /**
   * @brief Const getter for the map from cell to its first corner-point vertex
   * @return the map from cell to its first corner-point vertex
   */
  arrayView1d< localIndex const > cellToCPVertices() const
  { return m_cells.m_cellToCPVertices.toViewConst(); }

  /**
   * @brief Const getter for the map from local active cell index to global cell index
   * @return the map from from local active cell index to global cell index
   */
  arrayView1d< globalIndex const > ownedActiveCellToGlobalCell() const
  { return m_cells.m_ownedActiveCellToGlobalCell.toViewConst(); }

  // properties

  /**
   * @brief Const getter for the porosity field (includes inactive cells!)
   * @return an array of porosity values
   */
  arrayView1d< real64 const > porosityField() const { return m_parser.poro().toViewConst(); }

  /**
   * @brief Const getter for the porosity field (includes inactive cells!)
   * @return an array of permeability values
   */
  arrayView2d< real64 const > permeabilityField() const { return m_parser.perm().toViewConst(); }

  // regions

  /**
   * @brief Const getter for the region indices
   * @return an array of arrays mapping each region to its cells
   */
  ArrayOfArraysView< localIndex const > regionId() const { return m_regions.m_regionToOwnedActiveCell.toViewConst(); }

  // MPI information

  /**
   * @brief Const getter for the topological mesh information
   * @return the topological mesh information encapsulated in an object of type CornerPointMeshData
   */
  std::set< int > const & neighborsList() const { return m_partition.neighborsList(); }


private:

  /**
   * @brief Once the file is read, process the Eclipse keywords and fill the datastructure describing the conforming mesh
   * @details The post-processing proceeds in the following steps:
   *   1) For each corner-point cell, compute the position of the eight corner-point vertices
   *   2) Filter out duplicate vertices
   *   3) Construct faces (no implemented yet)
   */
  void postProcessMesh();

  /**
   * @brief Loop over corner-point cells to compute the position of the eight corner-point vertices in each cell
   */
  void buildCornerPointCells();

  /**
   * @brief Loop over the corner-point vertices and filter out duplicates.
   */
  void filterVertices();

  /**
   * @brief Associate each region with its cells
   */
  void formRegions();

  /**
   * @brief For the non-conforming case, build faces
   */
  void buildFaces();

  /**
   * @brief For the non-conforming case, build vertical faces
   * @param[in] direction the direction considered (either X or Y)
   */
  void buildVerticalFaces( string const & direction );

  /**
   * @brief For the non-conforming case, build horizontal faces
   */
  void buildHorizontalFaces();

  /**
   * @brief Using ZCORN and COORD, compute the position of the eight corner-point vertices in cell (i,j,k)
   * @param[in] i index in the x-direction
   * @param[in] j index in the y-direction
   * @param[in] k index in the z-direction
   * @param[in] nXLocal number of cells in the x-direction on my rank
   * @param[in] nYLocal number of cells in the y-direction on my rank
   * @param[in] coord local content of the COORD keyword
   * @param[in] zcorn local content of the ZCORN keyword
   * @param[out] xPos x-position of the eight corner-point vertices
   * @param[out] yPos y-position of the eight corner-point vertices
   * @param[out] zPos z-position of the eight corner-point vertices
   * @param[out] cpVertexIsIndex vector specifying if each corner-point vertex is inside the partition
   * @return true if the cell is valid, and false otherwise
   */
  static bool computeCellCoordinates( localIndex const i,
                                      localIndex const j,
                                      localIndex const k,
                                      localIndex const nXLocal,
                                      localIndex const nYLocal,
                                      localIndex const iMinOverlap,
                                      localIndex const iMaxOverlap,
                                      localIndex const jMinOverlap,
                                      localIndex const jMaxOverlap,
                                      array1d< real64 > const & coord,
                                      array1d< real64 > const & zcorn,
                                      array1d< real64 > & xPos,
                                      array1d< real64 > & yPos,
                                      array1d< real64 > & zPos,
                                      array1d< bool > & cpVertexIsInside );

  /**
   * @brief Populate the "reverse maps" from cellToActiveCell and activeCellToOwnedActiveCell
   * @param[in] nCells number of cells in the MPI partition
   * @param[in] nActiveCells number of active cells in the MPI partition
   * @param[in] activeCellToCell
   * @param[in] ownedActiveCellToActiveCell
   * @param[out] cellToActiveCell
   * @param[out] activeCellToOwnedActiveCell
   */
  static void populateReverseCellMaps( localIndex const nCells,
                                       localIndex const nActiveCells,
                                       array1d< localIndex > const & activeCellToCell,
                                       array1d< localIndex > const & ownedActiveCellToActiveCell,
                                       array1d< localIndex > & cellToActiveCell,
                                       array1d< localIndex > & activeCellToOwnedActiveCell );

  /**
   * @brief add a horizontal face made of four vertices to the maps
   * @param[in] faceVertices the four vertices of face k
   * @param[in] iOwnedActiveCellPrev index of the owned active cell before the face (k-1)
   * @param[in] iOwnedActiveCellNext index of the owned active cell after the face (k)
   * @param[out] ownedActiveCellToFaces map from owned active cell to faces
   * @param[out] faceToVertices map from face to vertices
   */
  static void addHorizontalFace( localIndex const (&faceVertices)[ 4 ],
                                 localIndex const iOwnedActiveCellPrev,
                                 localIndex const iOwnedActiveCellNext,
                                 ArrayOfArrays< localIndex > & ownedActiveCellToFaces,
                                 ArrayOfArrays< localIndex > & faceToVertices );


  /// Object holding the mesh (and MPI partition) dimensions
  CornerPointMeshDimensions m_dims;

  /// Object storing vertex data
  CornerPointMeshVertices m_vertices;

  /// Object storing face data
  CornerPointMeshFaces m_faces;

  /// Object storing cell data
  CornerPointMeshCells m_cells;

  /// Object storing regions
  CornerPointMeshRegions m_regions;

  /// Object in charge of parsing the GRDECL file
  CornerPointMeshParser m_parser;

  /// Object holding the communication info
  CornerPointMeshPartition m_partition;

  /// Name of the mesh
  string m_meshName;

};

} // namespace CornerPointMesh

} // namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHBUILDER_HPP_
