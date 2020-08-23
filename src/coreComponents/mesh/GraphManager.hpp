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
 * @file GraphManager.hpp
 */

#ifndef GEOSX_MESH_GRAPHMANAGER_HPP_
#define GEOSX_MESH_GRAPHMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
class SolverBase;

/**
 * @class GraphManager
 * @brief This class manages the graph objects in GEOSX (graph generated from text)
 */
class GraphManager : public dataRepository::Group
{
public:

  /**
   * @brief Constructor for the GraphManager object.
   * @param[in] name the name of the GraphManager object in the repository
   * @param[in] parent the parent group of the GraphManager object being constructed
   */
  GraphManager( std::string const & name,
               Group * const parent );

  virtual ~GraphManager() override;


  /**
   * @brief Create a new sub-graph.
   * @param[in] childKey the key of the new object in the ObjectCatalog
   * @param[in] childName the name of the new object in the collection of sub-graphes
   * @return A pointer to the Group node in the dataRepository of the new object created
   */
  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void ExpandObjectCatalogs() override;

  /**
   * @brief Generate the graphs
   */
  void GenerateGraphs();

private:

  /**
   * @brief Deleted default constructor of the GraphManager
   */
  GraphManager() = delete;

};

} /* namespace geosx */

#endif /* GEOSX_MESH_GRAPHMANAGER_HPP_ */
