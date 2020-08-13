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
 * @file MeshGeneratorBase.hpp
 */

#ifndef GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP
#define GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP

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
 *  @class MeshGeneratorBase
 *  @brief The MeshGeneratorBase class provides an abstract base class implementation for different mesh types.
 *	   The MeshGeneratorBase is the Group specialization for different type of mesh handling.
 */
class MeshGeneratorBase : public dataRepository::Group
{
public:
  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the MeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  explicit MeshGeneratorBase( std::string const & name, Group * const parent );

  /**
   * @brief Destructor for MeshGenerator
   */
  virtual ~MeshGeneratorBase();

  /**
   * @brief Return the name of the MeshGenerator in object catalog.
   * @return string that contains the catalog name of the MeshGenerator
   */
  static string
  CatalogName()
  {
    return "MeshGeneratorBase";
  }

  /**
 * @brief Generate the Element regions for an input Domain.
 * @param[inout] domain the Domain object on which to generate Element regions
 */
  virtual void
  GenerateElementRegions( DomainPartition & domain ) = 0;

  /**
 * @brief Generate the mesh object the input mesh object.
 * @param[in] domain the domain partition from which to construct the mesh object
 */
  virtual void
  GenerateMesh( DomainPartition * const domain ) = 0;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) = 0;

  /**
 * @brief Get the label mapping of element vertices indexes onto node indexes for a type of element.
 * @param[in] elementType the string identifier of the element type
 * @param[in] index ndim-sptialized Element index.
 * @param[in] iEle the index of Element begin processed
 * @param[out] nodeIDInBox array to map element vertices index to node indexes
 * @param[in] size the number of node on the element
 *
 */
  virtual void
  GetElemToNodesRelationInBox( const std::string & elementType,
                               const int index[],
                               const int & iEle,
                               int nodeIDInBox[],
                               const int size ) = 0;
  /**
 * @brief Re-computing mesh tables for the input domain.
 * @param[in] domain domain point whose mesh has to be remapped
 *
 */
  virtual void
  RemapMesh( dataRepository::Group * const domain ) = 0;

  /// Integer to trigger or not mesh re-mapping at the end of GenerateMesh call
  int m_delayMeshDeformation = 0;

  /// using alias for templated Catalog meshGenerator type
  using CatalogInterface =
    dataRepository::CatalogInterface< MeshGeneratorBase, std::string const &, Group * const >;

  /**
 * @brief Accessor for the singleton Catalog object
 * @return a static reference to the Catalog object
 *
 */
  static CatalogInterface::CatalogType &
  GetCatalog();
};
}  // namespace geosx

#endif /* GEOSX_MESHUTILITIES_MESHGENERATORBASE_HPP */
