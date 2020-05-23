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
 * @file PAMELAMeshGenerator.hpp
 */

#pragma once

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"

//This is an include of PAMELA
#include "Mesh/Mesh.hpp"
#include "MeshDataWriters/Writer.hpp"

#include "MeshGeneratorBase.hpp"

namespace geosx
{
	
/**
 *  @class PAMELAMeshGenerator
 *  @brief The PAMELAMeshGenerator class provides a class implementation of PAMELA genrated meshes
 *
 */
class PAMELAMeshGenerator : public MeshGeneratorBase
{
public:
/**
 * @brief Main constructor for MeshGenerator base class
 * @param[in] name of the PAMELAMeshGenerator object
 * @param[in] parent the parent Group pointer for the MeshGenerator object
 */
  PAMELAMeshGenerator( const std::string & name,
                       Group * const parent );

  virtual ~PAMELAMeshGenerator() override;

/**
 * @brief Return the name of the PAMELAMeshGenerator in object Catalog
 * @return string that contains the key name to PAMELAMeshGenerator in the Catalog
 */  
  static string CatalogName() { return "PAMELAMeshGenerator"; }

///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static auto filePathString = "file";
    constexpr static auto scaleString = "scale";
    constexpr static auto fieldsToImportString = "fieldsToImport";
    constexpr static auto fieldNamesInGEOSXString = "fieldNamesInGEOSX";
    constexpr static auto reverseZString = "reverseZ";
  };
/// @endcond 

  virtual void GenerateElementRegions( DomainPartition & domain ) override;

  virtual Group * CreateChild( string const & childKey, string const & childName ) override;

  virtual void GenerateMesh( DomainPartition * const domain ) override;

  virtual void GetElemToNodesRelationInBox ( const std::string & elementType,
                                             const int index[],
                                             const int & iEle,
                                             int nodeIDInBox[],
                                             const int size ) override;

  virtual void RemapMesh ( dataRepository::Group * const domain ) override;

protected:
  void PostProcessInput() override final;

private:

  /// Unique Pointer to the Mesh in the data structure of PAMELA.
  std::unique_ptr< PAMELA::Mesh >  m_pamelaMesh;

  /// Names of the fields to be copied from PAMELA to GEOSX data structure
  string_array m_fieldsToImport;

  /// Path to the mesh file
  Path m_filePath;

  /// Scale factor that will be applied to the point coordinates
  real64 m_scale;

  /// String array of the GEOSX user decalred fields
  string_array m_fieldNamesInGEOSX;

  /// z pointing direction flag, 0 (default) is upward, 1 is downward
  int m_isZReverse;

  /// map from PAMELA enumeration element type to string
  const std::unordered_map< PAMELA::ELEMENTS::TYPE, string, PAMELA::ELEMENTS::EnumClassHash > ElementToLabel
    =
    {
    { PAMELA::ELEMENTS::TYPE::VTK_VERTEX, "VERTEX"},
    { PAMELA::ELEMENTS::TYPE::VTK_LINE, "LINE"  },
    { PAMELA::ELEMENTS::TYPE::VTK_TRIANGLE, "TRIANGLE" },
    { PAMELA::ELEMENTS::TYPE::VTK_QUAD, "QUAD" },
    { PAMELA::ELEMENTS::TYPE::VTK_TETRA, "TETRA" },
    { PAMELA::ELEMENTS::TYPE::VTK_HEXAHEDRON, "HEX" },
    { PAMELA::ELEMENTS::TYPE::VTK_WEDGE, "WEDGE" },
    { PAMELA::ELEMENTS::TYPE::VTK_PYRAMID, "PYRAMID" }
    };
};

}
