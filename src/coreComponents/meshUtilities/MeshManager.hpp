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
 * @file MeshManager.hpp
 */

#ifndef GEOSX_MESHUTILITIES_MESHMANAGER_HPP_
#define GEOSX_MESHUTILITIES_MESHMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
class SolverBase;

/**
 * @class MeshManager
 * @brief This class manages the mesh objects in GEOSX (reservoir mesh, well mesh)
 */
class MeshManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the MeshManager object.
   * @param[in] name the name of the MeshManager object in the repository
   * @param[in] parent the parent group of the MeshManager object being constructed
   */
  MeshManager( std::string const & name, Group * const parent );

  virtual ~MeshManager() override;

  /**
   * @brief Create a new sub-mesh.
   * @param[in] childKey the key of the new object in the ObjectCatalog
   * @param[in] childName the name of the new object in the collection of sub-meshes
   * @return A pointer to the Group node in the dataRepository of the new object created
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /// This function is used to expand any catalogs in the data structure
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief Generate the meshes of the physical DomainPartition.
   * @param[in] domain a pointer to the physical DomainPartition
   */
  void
  GenerateMeshes( DomainPartition * const domain );

  /**
   * @brief Generate the different mesh levels in a MeshBody of the domain.
   * @param[in] domain a pointer to the physical DomainPartition
   */
  void
  GenerateMeshLevels( DomainPartition * const domain );

private:
  /**
   * @brief Deleted default constructor of the MeshManager
   */
  MeshManager() = delete;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MESHMANAGER_HPP_ */
