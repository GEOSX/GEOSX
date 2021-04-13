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
 * @file CornerPointMeshParser.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARSER_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARSER_HPP_

#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/ObjectCatalog.hpp"
#include "meshUtilities/cornerPointMesh/CornerPointMeshData.hpp"

namespace geosx
{

namespace cornerPointMesh
{

/**
 * @class CornerPointMeshParser
 * @brief This class is in charge of processing the GRDECL file
 */
class CornerPointMeshParser
{

public:

  /**
   * @brief Constructor.
   * @param name the name of the class
   */
  CornerPointMeshParser( string const & name );

  /**
   * @brief Destructor.
   */
  virtual ~CornerPointMeshParser() = default;

  /// Default copy constructor
  CornerPointMeshParser( CornerPointMeshParser const & ) = default;

  /// Default move constructor
  CornerPointMeshParser( CornerPointMeshParser && ) = default;

  /// Deleted copy assignment operator
  CornerPointMeshParser & operator=( CornerPointMeshParser const & ) = delete;

  /// Deleted move assignment operator
  CornerPointMeshParser & operator=( CornerPointMeshParser && ) = delete;

  /// using alias for templated Catalog CornerPointMeshParser type
  using CatalogInterface = dataRepository::CatalogInterface< CornerPointMeshParser,
                                                             string const & >;

  /**
   * @brief Getter for the singleton Catalog object
   * @return a static reference to the Catalog object
   */
  static typename CatalogInterface::CatalogType & getCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief Define the catalog name for this class
   * @return the catalog name
   */
  static string catalogName() { return "CornerPointMeshParser"; }

  /**
   * @brief Const getter for catalog name of this class
   * @return the catalog name
   */
  virtual string getCatalogName() const { return catalogName(); }

  /**
   * @brief Rank 0 reads the keyword DIMENS or SPECGRID to get the total number of cells (nX, nY, nZ)
   * @param[in] filePath the path to the file
   * @param[out] nX the number of cells in the X direction
   * @param[out] nY the number of cells in the Y direction
   * @param[out] nZ the number of cells in the Z direction
   */
  static void readNumberOfCells( Path const & filePath,
                                 localIndex & nX,
                                 localIndex & nY,
                                 localIndex & nZ );

  /**
   * @brief Each rank reads its part of the mesh
   * @details This involves reading the main keywords (COORD, ZCORN, ACTNUM), plus the property keywords
   * @param[in] filePath the path to the file
   * @param[in] dims the topogical mesh information (to get nX, nY, nZ, size of partition, etc)
   */
  void readMesh( Path const & filePath,
                 CornerPointMeshDimensions const & dims );

  /**
   * @brief Non-const getter for the local content of COORD
   * @return the local content of COORD
   */
  array1d< real64 > & coord() { return m_coord; }

  /**
   * @brief Const getter for the local content of COORD
   * @return the local content of COORD
   */
  array1d< real64 > const & coord() const { return m_coord; }

  /**
   * @brief Non-const getter for the local content of ZCORN
   * @return the local content of ZCORN
   */
  array1d< real64 > & zcorn() { return m_zcorn; }

  /**
   * @brief Const getter for the local content of ZCORN
   * @return the local content of ZCORN
   */
  array1d< real64 > const & zcorn() const { return m_zcorn; }

  /**
   * @brief Non-const getter for the local content of ACTNUM
   * @return the local content of ACTNUM
   */
  array1d< localIndex > & actnum() { return m_actnum; }

  /**
   * @brief Const getter for the local content of ACTNUM
   * @return the local content of ACTNUM
   */
  array1d< localIndex > const & actnum() const { return m_actnum; }

  /**
   * @brief Non-const getter for the local content of permeability
   * @return the local content of permeability
   */
  array2d< real64 > & perm() { return m_perm; }

  /**
   * @brief Const getter for the local content of permeability
   * @return the local content of permeability
   */
  array2d< real64 > const & perm() const { return m_perm; }

  /**
   * @brief Non-const getter for the local content of PORO
   * @return the local content of PORO
   */
  array1d< real64 > & poro() { return m_poro; }

  /**
   * @brief Const getter for the local content of PORO
   * @return the local content of PORO
   */
  array1d< real64 > const & poro() const { return m_poro; }

  /**
   * @brief Non-const getter for the local content of REGION
   * @return the local content of REGION
   */
  array1d< localIndex > & regionId() { return m_regionId; }

  /**
   * @brief Const getter for the local content of REGION
   * @return the local content of REGION
   */
  array1d< localIndex > const & regionId() const { return m_regionId; }


private:

  /**
   * @brief Utility (temporary) functions that reads all the data below a keyword
   * @details This is a carry-over from PAMELA, but this function should go if we want to load only local data
   * (i.e., data corresponding to the MPI domain owned by each rank)
   */
  static std::string extractDataBelowKeyword( std::istringstream & stringBlock );

  /**
   * @brief Read the local content of the COORD keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param dims the datastructure describing CP mesh
   */
  void readLocalCOORD( std::istringstream & meshFile,
                       CornerPointMeshDimensions const & dims );

  /**
   * @brief Read the local content of the ZCORN keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param dims the datastructure describing CP mesh
   */
  void readLocalZCORN( std::istringstream & meshFile,
                       CornerPointMeshDimensions const & dims );

  /**
   * @brief Read the local content of the ACTNUM keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @param meshFile the content of the mesh file
   * @param dims the datastructure describing CP mesh
   */
  void readLocalACTNUM( std::istringstream & meshFile,
                        CornerPointMeshDimensions const & dims );

  /**
   * @brief Read the local content of a cell-centered property keyword (i.e, data corresponding to the MPI domain owned by each rank)
   * @details This function will have to be rewritten
   * @tparam T the type of the data in prop
   * @param meshFile the content of the mesh file
   * @param dims the datastructure describing CP mesh
   * @param prop the property array to fill
   */
  template< typename T >
  void readLocalPROP( std::istringstream & meshFile,
                      CornerPointMeshDimensions const & dims,
                      array1d< T > & prop );

  /**
   * @brief If the property is absent from the file, fill the vector with the default value
   * @tparam T the type of the default value
   * @param dims the datastructure describing CP mesh
   * @param prop the property array to fill
   * @param defaultValue the default value to use here
   */
  template< typename T >
  void fillLocalPROP( CornerPointMeshDimensions const & dims,
                      array1d< T > & prop,
                      T defaultValue );



  /// local content of the COORD keyword in the GRDECL file
  array1d< real64 > m_coord;

  /// local content of the ZCORN keyword in the GRDECL file
  array1d< real64 > m_zcorn;

  /// local content of the ACTNUM keyword in the GRDECL file
  array1d< localIndex > m_actnum;

  /// local content of the PERMX keyword in the GRDECL file
  array1d< real64 > m_permx;

  /// local content of the PERMY keyword in the GRDECL file
  array1d< real64 > m_permy;

  /// local content of the PERMZ keyword in the GRDECL file
  array1d< real64 > m_permz;

  /// contenation of the local content of the permeability keyword
  array2d< real64 > m_perm;

  /// local content of the PORO keyword in the GRDECL file
  array1d< real64 > m_poro;

  /// local content of the REGION keywork in the GRDECL file
  array1d< localIndex > m_regionId;

  /// name of the mesh
  string m_meshName;

};

} // namespace cornerPointMesh

} // namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_CORNERPOINTMESHPARSER_HPP_
