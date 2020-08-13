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
 * @file FunctionManager.hpp
 */

#ifndef GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "FunctionBase.hpp"

namespace geosx
{
/**
 * @class FunctionManager
 *
 * A class to manage arbitrary N-dimensional functions.
 */
class FunctionManager : public dataRepository::Group
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  FunctionManager( const std::string & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~FunctionManager() override;

  /**
   * @brief Return the function manager instance
   * @return the function manager instance
   */
  static FunctionManager &
  Instance()
  {
    static FunctionManager theFunctionManager( "Functions", nullptr );
    return theFunctionManager;
  }

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string
  CatalogName()
  {
    return "FunctionManager";
  }

  /**
   * @brief Create a new FunctionManager object as a child of this group.
   * @param functionCatalogKey the catalog key of the new FunctionManager derived type to create
   * @param functionName the name of the new FunctionManager object in the repository
   * @return the group child
   */
  virtual Group *
  CreateChild( string const & functionCatalogKey,
               string const & functionName ) override;

  /**
   * @brief This function is used to expand any catalogs in the data structure
   */
  virtual void
  ExpandObjectCatalogs() override;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_FUNCTIONMANAGER_HPP_ */
