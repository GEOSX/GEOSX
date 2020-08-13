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
 * @file FaceElementRegion.hpp
 *
 */

#ifndef GEOSX_MESH_FACEELEMENTREGION_HPP_
#define GEOSX_MESH_FACEELEMENTREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{
class EdgeManager;

/**
 * @class FaceElementRegion
 *
 * The FaceElementRegion class contains the functionality to support the concept of a FaceElementRegion in the element
 * hierarchy. FaceElementRegion derives from ElementRegion and has an entry in the ObjectManagerBase catalog.
 *
 *
 */
class FaceElementRegion : public ElementRegionBase
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
  FaceElementRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  FaceElementRegion() = delete;

  /**
   * @brief Default destructor.
   */
  virtual ~FaceElementRegion() override;

  ///@}

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief The key name for the FaceElementRegion in the object catalog.
   * @return a string containing the key name.
   */
  static const string
  CatalogName()
  {
    return "FaceElementRegion";
  }

  virtual const string
  getCatalogName() const override final
  {
    return FaceElementRegion::CatalogName();
  }

  ///@}

  /**
   * @name Generation of the face element mesh region
   */
  ///@{

  virtual void
  GenerateMesh( Group * ) override
  {}

  /**
   * @brief This function generates and adds entries to the face/fracture mesh.
   * @param time_np1 rupture time
   * @param edgeManager pointer to the EdgeManager object.
   * @param faceManager pointer to the FaceManager object.
   * @param originalFaceToEdges face-to-edge map before the rupture.
   * @param subRegionName the name of the FaceElementSubRegion to insert the new entries.
   * @param faceIndices the local indices of the new faces that define the face element.
   * @return the local index of the new FaceElement entry.
   */
  localIndex
  AddToFractureMesh(
    real64 const time_np1,
    EdgeManager * const edgeManager,
    FaceManager const * const faceManager,
    ArrayOfArraysView< localIndex const > const & originalFaceToEdges,
    string const & subRegionName,
    localIndex const faceIndices[2] );

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get default aperture value.
   * @return default aperture value
   */
  real64
  getDefaultAperture() const
  {
    return m_defaultAperture;
  }

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// Fracture set string
    static constexpr auto fractureSetString = "fractureSet";

    /// Default aperture string
    static constexpr auto defaultApertureString = "defaultAperture";

    /// Rupture time string
    constexpr static auto ruptureTimeString = "ruptureTime";
  };

protected:
  virtual void
  InitializePreSubGroups( Group * const ) override;

private:
  /// The default aperture
  real64 m_defaultAperture;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_FACEELEMENTREGION_HPP_ */
