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
 * @file schemaUtilities.hpp
 */

#ifndef GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_
#define GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_

#include "dataRepository/xmlWrapper.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include <iostream>
#include <sstream>
#include <string>

namespace geosx
{
// Forward declarations
namespace dataRepository
{
class Group;
}

namespace schemaUtilities
{
/**
 * @brief Generates the XML schema.
 *
 * @param fname schema filename
 * @param group group passed in from which the schema information is extracted
 * @param documentationType type of XML schema generated
 */
void
ConvertDocumentationToSchema( std::string const & fname,
                              dataRepository::Group * const group,
                              integer documentationType );

/**
 * @brief Generates the simple schema types.
 *
 * @param schemaRoot XML node corresponding to the root
 */
void
BuildSimpleSchemaTypes( xmlWrapper::xmlNode schemaRoot );

/**
 * @brief Recursively builds the schema from the data structure skeleton.
 *
 * @param group group passed in from which the schema information is extracted
 * @param schemaRoot XML node corresponding to the root
 * @param schemaParent XML node for the parent node
 * @param documentationType type of XML schema generated
 */
void
SchemaConstruction( dataRepository::Group * const group,
                    xmlWrapper::xmlNode schemaRoot,
                    xmlWrapper::xmlNode schemaParent,
                    integer documentationType );

}  // namespace schemaUtilities
}  // namespace geosx

#endif /* GEOSX_FILEIO_SCHEMA_SCHEMAUTILITIES_HPP_ */
