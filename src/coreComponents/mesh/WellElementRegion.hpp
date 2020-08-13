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
 * @file WellElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_WELLELEMENTREGION_HPP_
#define GEOSX_MESH_WELLELEMENTREGION_HPP_

#include "mesh/ElementRegionBase.hpp"
#include "meshUtilities/InternalWellGenerator.hpp"

namespace geosx
{
class MeshLevel;

/**
 * @class WellElementRegion
 * @brief This class specializes the element region for the case
 *        of a well. This class is also in charge of starting the
 *        construction of the well data structure in GenerateWell.
 */
class WellElementRegion : public ElementRegionBase
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name name of the object in the data hierarchy.
   * @param parent pointer to the parent group in the data hierarchy.
   */
  WellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~WellElementRegion() override;

  /**
   * @brief Deleted default constructor.
   */
  WellElementRegion() = delete;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static const string
  CatalogName()
  {
    return "WellElementRegion";
  }

  /**
   * @copydoc CatalogName()
   */
  virtual const string
  getCatalogName() const override final
  {
    return WellElementRegion::CatalogName();
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the name of the InternalWellGenerator object of this well.
   * @param[in] name the name of the InternalWellGenerator object
   */
  void
  SetWellGeneratorName( string const & name )
  {
    m_wellGeneratorName = name;
  }

  /**
   * @brief Get the name of the InternalWellGenerator object of this well.
   * @return the name of the InternalWellGenerator object
   */
  string const &
  GetWellGeneratorName() const
  {
    return m_wellGeneratorName;
  }

  /**
   * @brief Set the name of the WellControls object of this well.
   * @param name the name of the WellControls object
   */
  void
  SetWellControlsName( string const & name )
  {
    m_wellControlsName = name;
  }

  /**
   * @brief Get the name of the subRegion.
   * @return the name of the subRegion object
   */
  string const &
  GetSubRegionName() const
  {
    return m_subRegionName;
  }

  ///@}

  /**
   * @name Construction of the well connectivity
   */
  ///@{

  /**
   * @brief Not implemented, this task is performed in GenerateWell.
   */
  virtual void
  GenerateMesh( Group * ) override
  {}

  /**
   * @brief Build the local well elements and perforations from global well geometry.
   * @param[in] mesh the mesh object (single level only)
   * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
   * @param[in] nodeOffsetGlobal the offset of the first global well node ( = offset of last global mesh node + 1 )
   * @param[in] elemOffsetGlobal the offset of the first global well element ( = offset of last global mesh elem + 1 )
   */
  void
  GenerateWell( MeshLevel & mesh,
                InternalWellGenerator const & wellGeometry,
                globalIndex nodeOffsetGlobal,
                globalIndex elemOffsetGlobal );

  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// String key for the well control name
    static constexpr auto wellControlsString = "wellControlsName";
    /// String key for the well generator name
    static constexpr auto wellGeneratorString = "wellGeneratorName";

    /// ViewKey for the well control name
    dataRepository::ViewKey wellControlsName = { wellControlsString };
    /// ViewKey for the well generator name
    dataRepository::ViewKey wellGeneratorName = { wellGeneratorString };

  }
  /// ViewKey struct for the WellElementRegion class
  viewKeysWellElementRegion;

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ElementRegionBase::groupKeyStruct
  {
  }
  /// groupKey struct for the WellElementRegion class
  groupKeysWellElementRegion;

private:
  /// Name of the (unique) subregion
  const string m_subRegionName;

  /// Name of the well constraint
  string m_wellControlsName;

  /// Name of the well generator
  string m_wellGeneratorName;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_WELLELEMENTREGION_HPP_ */
