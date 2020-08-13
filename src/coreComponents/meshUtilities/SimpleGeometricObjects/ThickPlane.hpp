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
 * @file ThickPlane.hpp
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{
/**
 * @class ThickPlane
 * @brief Class to represent a geometric thick plane  in GEOSX.
 */
class ThickPlane : public SimpleGeometricObjectBase
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
  ThickPlane( const std::string & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~ThickPlane() override;

  ///@}

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{

  /**
   * @brief Get the catalog name.
   * @return the name of this class in the catalog
   */
  static string
  CatalogName()
  {
    return "ThickPlane";
  }

  ///@}

  bool
  IsCoordInObject( const R1Tensor & coord ) const override final;

  /**
   * @name Getters
   */
  ///@{

  /**
   * @brief Get the normal to the plane.
   * @return the normal vector
   */
  R1Tensor &
  getNormal()
  {
    return m_normal;
  }

  /**
   * @copydoc getNormal()
   */
  R1Tensor const &
  getNormal() const
  {
    return m_normal;
  }

  /**
   * @brief Get the origin of the plane.
   * @return the origin of the plane
   */
  R1Tensor &
  getCenter()
  {
    return m_origin;
  }

  /**
   * @copydoc getCenter()
   */
  R1Tensor const &
  getCenter() const
  {
    return m_origin;
  }

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void
  PostProcessInput() override final;

private:
  /// Origin point (x,y,z) of the plane (basically, any point on the plane)
  R1Tensor m_origin;
  /// Normal (n_x,n_y,n_z) to the plane (will be normalized automatically)
  R1Tensor m_normal;
  /// Total thickness of the plane (with half to each side)
  real64 m_thickness;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr auto originString = "origin";
    static constexpr auto normalString = "normal";
    static constexpr auto thicknessString = "thickness";
  };

  /// @endcond
};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_THICKPLANE_HPP_
        */
