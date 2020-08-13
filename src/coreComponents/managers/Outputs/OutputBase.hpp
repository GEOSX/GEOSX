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
 * @file OutputBase.hpp
 */
#ifndef GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_
#define GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_

#include "dataRepository/Group.hpp"
#include "dataRepository/ExecutableGroup.hpp"

namespace geosx
{
/**
 * @class OutputBase
 *
 * A base class for output types.
 */
class OutputBase : public ExecutableGroup
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  explicit OutputBase( std::string const & name, Group * const parent );

  /// Destructor
  virtual ~OutputBase() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string
  CatalogName()
  {
    return "OutputBase";
  }

  /// Method for setting up output directories.
  virtual void
  SetupDirectoryStructure();

  // Catalog interface
  /// @cond DO_NOT_DOCUMENT
  using CatalogInterface =
    dataRepository::CatalogInterface< OutputBase, std::string const &, Group * const >;
  static CatalogInterface::CatalogType &
  GetCatalog();

  // Catalog view keys
  struct viewKeysStruct
  {
    static constexpr auto childDirectoryString = "childDirectory";
    static constexpr auto parallelThreadsString = "parallelThreads";
  } outputBaseViewKeys;
  /// @endcond

  /**
   * @brief Get the path of the child directory where output will be written
   * @return The directory path
   **/
  string
  childDirectory() const
  {
    return m_childDirectory;
  }

  /**
   * @brief Get the number of parallel threads to use to write plotfiles
   * @return The number of threads
   **/
  integer
  parallelThreads() const
  {
    return m_parallelThreads;
  }

protected:
  /**
   * @brief Do initialization prior to calling initialization operations
   *        on the subgroups.
   * @param group The root group
   **/
  virtual void
  InitializePreSubGroups( Group * const group ) override;

private:
  string m_childDirectory;
  integer m_parallelThreads;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_OUTPUTBASE_HPP_ */
