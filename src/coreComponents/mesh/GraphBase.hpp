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
 * @file GraphBase.hpp
 */

#ifndef GEOSX_MESH_GRAPHBASE_HPP
#define GEOSX_MESH_GRAPHBASE_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

namespace dataRepository
{}

class NodeManager;
class DomainPartition;

/**
 *  @class GraphBase
 *  @brief The GraphBase class provides an abstract base class implementation for different graph types.
 *	   The GraphBase is the Group specialization for different type of graph handling.
 */
class GraphBase : public dataRepository::Group
{
public:

  /**
   * @brief Main constructor for Graph base class.
   * @param[in] name of the Graph object
   * @param[in] parent the parent Group pointer for the Graph object
   */
  explicit GraphBase( std::string const & name,
                              Group * const parent );

  /**
   * @brief Destructor for Graph
   */
  virtual ~GraphBase();

  /**
   * @brief Function called to construct the graphs
  */ 
  virtual void GenerateGraph() = 0;


  /**
   * @brief Return the name of the Graph in object catalog.
   * @return string that contains the catalog name of the Graph
   */
  static string CatalogName() { return "GraphBase"; }


  /// using alias for templated Catalog graph type
  using CatalogInterface = dataRepository::CatalogInterface< GraphBase, std::string const &, Group * const >;

/**
 * @brief Accessor for the singleton Catalog object
 * @return a static reference to the Catalog object
 *
 */
  static CatalogInterface::CatalogType & GetCatalog();

};
}

#endif /* GEOSX_MESH_GRAPHBASE_ */
