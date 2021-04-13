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
 * @file CornerPointMeshGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP
#define GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/StringUtilities.hpp"
#include "common/EnumStrings.hpp"
#include "meshUtilities/MeshGeneratorBase.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshBuilder.hpp"

namespace geosx
{

/**
 *  @class CornerPointMeshGenerator
 *  @brief Class that reads and processes a corner point grid
 */
class CornerPointMeshGenerator : public MeshGeneratorBase
{
public:

  /**
   * @enum  PermeabilityUnit
   * @brief This enum class is used to specify the unit of permeability in the GRDECL file (they will be converted to m^2)
   */
  enum class PermeabilityUnit : integer
  {
    Millidarcy, //!< unit is milli-darcy
    SquareMeter,  //!< unit is square-meter
  };

  /**
   * @enum  CoordinatesUnit
   * @brief This enum class is used to specify the unit of the coordinates in the GRDECL file (they will be converted to m)
   */
  enum class CoordinatesUnit : integer
  {
    Meter, //!< unit is meter
    Foot,  //!< unit is foot
  };

/**
 * @brief Main constructor for MeshGenerator base class.
 * @param[in] name of the CornerPointMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  CornerPointMeshGenerator( const string & name,
                            Group * const parent );

  virtual ~CornerPointMeshGenerator() override;

/**
 * @brief Return the name of the CornerPointMeshGenerator in object Catalog.
 * @return string that contains the key name to CornerPointMeshGenerator in the Catalog
 */
  static string catalogName() { return "ExternalCornerPointMesh"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * filePathString() { return "file"; }
    constexpr static char const * permeabilityUnitInInputFileString() { return "permeabilityUnitInInputFile"; }
    constexpr static char const * coordinatesUnitInInputFileString() { return "coordinatesUnitInInputFile"; }
  };
/// @endcond

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey,
                               string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

protected:

  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void postProcessInput() override final;


private:

  /// Driver class for the construction of the conformal corner-point mesh
  std::unique_ptr< cornerPointMesh::CornerPointMeshBuilder > m_cpMeshBuilder;

  /// Path to the mesh file
  Path m_filePath;

  /// Permeability unit
  PermeabilityUnit m_permeabilityUnitInInputFile;

  /// Coordinates unit
  CoordinatesUnit m_coordinatesUnitInInputFile;

  /// Conversion factor for permeability
  real64 m_toSquareMeter;

  /// Conversion factor for coordinates
  real64 m_toMeter;

};

ENUM_STRINGS( CornerPointMeshGenerator::PermeabilityUnit, "Millidarcy", "SquareMeter" )
ENUM_STRINGS( CornerPointMeshGenerator::CoordinatesUnit, "Meter", "Foot" )

} // end namespace geosx

#endif /* GEOSX_MESHUTILITIES_CORNERPOINTGENERATOR_HPP */
