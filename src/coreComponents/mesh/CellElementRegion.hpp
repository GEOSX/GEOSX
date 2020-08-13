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
 * @file CellElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_CELLELEMENTREGION_HPP_
#define GEOSX_MESH_CELLELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;
class EmbeddedSurfaceGenerator;

/**
 * @class CellElementRegion
 *
 * The CellElementRegion class contains the functionality to support the concept of a CellElementRegion in the element
 * hierarchy. CellElementRegion derives from ElementRegionBase and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class CellElementRegion : public ElementRegionBase
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the name of the object in the data hierarchy.
   * @param parent a pointer to the parent group in the data hierarchy.
   */
  CellElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  CellElementRegion() = delete;

  /**
   * @brief Destructor.
   */
  virtual ~CellElementRegion() override;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string
  CatalogName()
  {
    return "CellElementRegion";
  }

  /**
   * @copydoc CatalogName()
   */
  virtual const string
  getCatalogName() const override final
  {
    return CellElementRegion::CatalogName();
  }

  ///@}

  /**
   * @name Generation of the cell element mesh region
   */
  ///@{

  /**
   * @brief Add a cellBlockRegion name to the list.
   * @param cellBlockName string containing the cell block region name.
   */
  void
  AddCellBlockName( string const & cellBlockName )
  {
    m_cellBlockNames.emplace_back( cellBlockName );
  }

  virtual void
  GenerateMesh( Group * const cellBlocks ) override;

  /**
   * @brief Generate the aggregates.
   * @param faceManager a pointer to the FaceManager
   * @param nodeManager a pointer to the NodeManager
   */
  void
  GenerateAggregates( FaceManager const * const faceManager,
                      NodeManager const * const nodeManager );

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// String key for the coarsening ratio
    static constexpr auto coarseningRatioString = "coarseningRatio";

    /// String key for the cell block names
    static constexpr auto sourceCellBlockNames = "cellBlocks";
  };

private:
  // Cell block names
  string_array m_cellBlockNames;

  // Coarsening ratio
  real64 m_coarseningRatio;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTREGION_HPP_ */
