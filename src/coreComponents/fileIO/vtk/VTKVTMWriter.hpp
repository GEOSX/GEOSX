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

#ifndef GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_
#define GEOSX_FILEIO_VTK_VTKVTMWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geosx
{
namespace vtk
{
/*!
 * @brief VTM Writer class.
 * @details a VTM file is the root file for one time step. It will contain
 * path to all the component of the mesh (Surfaces, Volumes, Wells etc.).
 */
class VTKVTMWriter
{
public:
  /*!
   * @brief Build the VTM Writer
   * @param[in] filePath path to the file
   */
  VTKVTMWriter( string const & filePath );

  /*!
   * @brief Triggers the file output
   */
  void
  Save() const;

  /*!
   * @brief Add a block to the VTM file
   * @details The first level of block is for the ElementRegion (\p blockName can
   * be CellElementRegion, FaceElementRegion or WellElementREgion)
   * @param[in] blockName Name of the block
   */
  void
  AddBlock( string const & blockName ) const;

  /*!
   * @brief Add a subblock to the VTM file
   * @details The second level of block is for the different Regions
   * @param[in] blockName Name of the parent block
   * @param[in] subBlockName Name of the subBlock (usually the name of the Region)
   */
  void
  AddSubBlock( string const & blockName, string const & subBlockName ) const;

  /*!
   * @brief Add data to the subblock \p subBlockName
   * @details The final level : paths to the vtu file per rank
   * @param[in] blockName Name of the parent block
   * @param[in] subBlockName Name of the subBlock (usually the name of the Region)
   * @param[in] filePath path to the vtu file containing the unstructured mesh
   * @param[in] mpiRank the mpi rank.
   */
  void
  AddDataToSubBlock( string const & blockName,
                     string const & subBlockName,
                     string const & filePath,
                     int mpiRank ) const;

private:
  /// VTM XML file
  xmlWrapper::xmlDocument m_vtmFile;

  /// Path to the XML File
  string const m_filePath;
};
}  // namespace vtk
}  // namespace geosx

#endif
