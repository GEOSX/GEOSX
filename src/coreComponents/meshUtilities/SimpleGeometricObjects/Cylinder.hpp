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
 * @file Cylinder.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
#define GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_

#include "SimpleGeometricObjectBase.hpp"

namespace geosx
{
/**
 * @class Cylinder
 * @brief Class to represent a geometric cylinder in GEOSX.
 */
class Cylinder : public SimpleGeometricObjectBase
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
  Cylinder( const std::string & name, Group * const parent );

  /**
   * @brief Default destructor.
   */
  virtual ~Cylinder() override;

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
    return "Cylinder";
  }

  ///@}

  bool
  IsCoordInObject( const R1Tensor & coord ) const override final;

private:
  /// Center point of one (upper or lower) face of the cylinder
  R1Tensor m_point1;
  /// Center point of the other face of the cylinder
  R1Tensor m_point2;
  /// Radius of the cylinder
  realT m_radius = 0.0;

  /// @cond DO_NOT_DOCUMENT

  struct viewKeyStruct
  {
    static constexpr auto point1String = "point1";
    static constexpr auto point2String = "point2";
    static constexpr auto radiusString = "radius";
  };

  /// @endcond
};
} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_SIMPLEGEOMETRICOBJECTS_CYLINDER_HPP_
        */
