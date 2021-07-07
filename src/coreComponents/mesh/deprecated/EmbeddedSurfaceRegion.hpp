/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file EmbeddedSurfaceRegion.hpp
 *
 */

#ifndef GEOSX_MESH_EMBEDDEDSURFACEREGION_HPP_
#define GEOSX_MESH_EMBEDDEDSURFACEREGION_HPP_

#include "ElementRegionBase.hpp"

namespace geosx
{

/**
 * @class EmbeddedSurfaceRegion
 *
 * The EmbeddedSurfaceRegion class contains the functionality to support the concept of a EmbeddedSurfaceRegion in the
 * element hierarchy. EmbeddedSurfaceRegion derives from ElementRegion and has an entry in the ObjectManagerBase
 * catalog.
 */
class EmbeddedSurfaceRegion : public ElementRegionBase
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
  EmbeddedSurfaceRegion( string const & name, Group * const parent );

  /**
   * @brief Deleted default constructor.
   */
  EmbeddedSurfaceRegion() = delete;

  /**
   * @brief Default destructor.
   */
  virtual ~EmbeddedSurfaceRegion() override;

  ///@}

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get the key name for the EmbeddedSurfaceRegion in the object catalog.
   * @return A string containing the key name.
   */
  static const string catalogName()
  { return "EmbeddedSurfaceElementRegion"; }

  virtual const string getCatalogName() const override final
  { return EmbeddedSurfaceRegion::catalogName(); }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get default aperture value.
   * @return default aperture value
   */
  real64 getDefaultAperture() const { return m_defaultAperture; }

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementRegionBase::viewKeyStruct
  {
    /// @return Fracture set string
    static constexpr char const * fractureSetString() { return "fractureSet"; }
    /// @return Default fracture aperture
    static constexpr char const * defaultApertureString() { return "defaultAperture"; }
  };

protected:
  virtual void initializePreSubGroups() override;

private:

  real64 m_defaultAperture;

};

} /* namespace geosx */

#endif /* CORECOMPONENTS_MESH_EMBEDDEDSURFACEREGION_HPP_ */
