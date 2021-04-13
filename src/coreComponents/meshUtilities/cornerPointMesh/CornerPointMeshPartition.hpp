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
 * @file CornerPointMeshPartition.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARTITION_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARTITION_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshData.hpp"

namespace geosx
{

namespace cornerPointMesh
{

/**
 * @class CornerPointMeshPartition
 * @brief This class is in charge of MPI communications during the construction of the mesh
 */
class CornerPointMeshPartition
{

public:

  /**
   * @brief Constructor.
   * @param name the name of the class
   */
  CornerPointMeshPartition( string const & name );

  /// Default virtual destructor
  virtual ~CornerPointMeshPartition() = default;

  /// Default copy constructor
  CornerPointMeshPartition( CornerPointMeshPartition const & ) = default;

  /// Default move constructor
  CornerPointMeshPartition( CornerPointMeshPartition && ) = default;

  /// Deleted copy assignment operator
  CornerPointMeshPartition & operator=( CornerPointMeshPartition const & ) = delete;

  /// Deleted move assignment operator
  CornerPointMeshPartition & operator=( CornerPointMeshPartition && ) = delete;

  /// using alias for templated Catalog CornerPointMeshBuilder type
  using CatalogInterface = dataRepository::CatalogInterface< CornerPointMeshPartition,
                                                             string const & >;

  /**
   * @brief Getter for the singleton Catalog object
   * @return a static reference to the Catalog object
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
  static string catalogName() { return "CornerPointMeshPartition"; }

  /**
   * @brief Const getter for catalog name of this class
   * @return the catalog name
   */
  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Const getter for the list of my neighbor ranks
   * @return a set containing the indices of neighbors
   */
  std::set< int > const & neighborsList() const { return m_neighborsList; }

  /**
   * @brief The main rank partitions the domain, and the communication pattern is set up
   * @param[out] dims the dimensions of the problem / the MPI partition
   * @details This involves multiple steps:
   *  1) The main rank broadcasts the full domain boundaries (nX, nY, nZ)
   *  2) The main rank partitions the domain in boxes using a simple approach
   *  3) The main rank scatters the partition boundaries
   *  4) The main rank uses the structure of the mesh to figure out the communication pattern
   */
  void setupMPIPartition( CornerPointMeshDimensions & dims );

private:

  /**
   * @brief The main rank (that has read the number of cells) broadcasts full domain dimensions (nX, nY, nZ)
   * @param[inout] dims the struct describing the local and global mesh dimensions
   */
  void broadcastDomainDimensions( CornerPointMeshDimensions & dims ) const;

  /**
   * @brief The main rank partitions the full domain and scatters the MPI partition boundaries
   * @param[inout] dims the struct describing the local and global mesh dimensions
   */
  void scatterPartitionBoundaries( CornerPointMeshDimensions & dims );


  /// Index of the main rank
  localIndex m_mainRank;

  /// Index of my rank
  localIndex m_myRank;

  /// Number of ranks
  localIndex m_size;

  /// Indices of my neighbors
  std::set< int > m_neighborsList;

  /// Name of the mesh
  string m_meshName;

};

} // namespace cornerPointMesh

} // namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARTITION_HPP_
