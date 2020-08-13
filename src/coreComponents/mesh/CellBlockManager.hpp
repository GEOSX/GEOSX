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
 * @file CellBlockManager.hpp
 */

#ifndef GEOSX_MESH_CELLBLOCKMANAGER_H_
#define GEOSX_MESH_CELLBLOCKMANAGER_H_

#include "CellBlock.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
/// String for cellBlocks
string const cellBlocks = "cellBlocks";
}  // namespace keys
}  // namespace dataRepository

/**
 * @class CellBlockManager
 * @brief The CellBlockManager class provides an interface to ObjectManagerBase in order to manage CellBlock data.
 */
class CellBlockManager : public ObjectManagerBase
{
public:
  /**
   * @brief The function is to return the name of the CellBlockManager in the object catalog
   * @return string that contains the catalog name used to register/lookup this class in the object catalog
   */
  static string
  CatalogName()
  {
    return "CellBlockManager";
  }

  virtual const string
  getCatalogName() const override final
  {
    return CellBlockManager::CatalogName();
  }

  /**
   * @brief Constructor for CellBlockManager object.
   * @param name name of this instantiation of CellBlockManager
   * @param parent pointer to the parent Group of this instantiation of CellBlockManager
   */
  CellBlockManager( string const & name, Group * const parent );

  /**
   * @brief Destructor
   */
  virtual ~CellBlockManager() override;

  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  using Group::resize;

  /**
   * @brief Set the number of elements for a set of element regions.
   * @param numElements list of the new element numbers
   * @param regionNames list of the element region names
   * @param elementTypes list of the element types
   */
  void
  resize( integer_array const & numElements,
          string_array const & regionNames,
          string_array const & elementTypes );

  //  CellBlock & CreateRegion( string const & regionName,
  //                               string const & elementType,
  //                               integer const & numElements );

  /**
   * @brief Get element sub-region.
   * @param regionName name of the element sub-region
   * @return pointer to the element sub-region
   */
  CellBlock *
  GetRegion( string const & regionName )
  {
    return this->GetGroup( dataRepository::keys::cellBlocks )
      ->GetGroup< CellBlock >( regionName );
  }

  /**
   * @brief Launch kernel function over all the sub-regions
   * @tparam LAMBDA type of the user-provided function
   * @param lambda kernel function
   */
  template< typename LAMBDA >
  void
  forElementSubRegions( LAMBDA lambda )
  {
    Group * elementRegions = this->GetGroup( dataRepository::keys::cellBlocks );
    elementRegions->forSubGroups< CellBlock >( lambda );
  }

private:
  /**
   * @brief Copy constructor.
   */
  CellBlockManager( const CellBlockManager & );

  /**
   * @brief Copy assignment operator.
   * @return reference to this object
   */
  CellBlockManager &
  operator=( const CellBlockManager & );
};
}  // namespace geosx
#endif /* GEOSX_MESH_CELLBLOCKMANAGER_H_ */
