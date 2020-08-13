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

#ifndef GEOSX_FILEIO_VTK_VTKPVDWRITER_HPP_
#define GEOSX_FILEIO_VTK_VTKPVDWRITER_HPP_

#include "dataRepository/xmlWrapper.hpp"

namespace geosx
{
namespace vtk
{
/*!
 * @brief VTK PVD Writer class.
 * @details the PVD file is the root file, it contains the
 * reference to all the VTM files (one VTM file per output time-step).
 * This file is meant to be open with paraview.
 */
class VTKPVDWriter
{
public:
  /*!
   * @brief Writer for PVD file
   * @param[in] fileName the file name with the extension
   */
  VTKPVDWriter( string const & fileName );

  /*!
   * @brief Triggers the file output
   */
  void
  Save() const;

  /*!
   * @brief Add a dataset associated to a time-step
   * @param[in] time the time step
   * @param[in] filePath path to the file associated with the time-step
   */
  void
  AddData( real64 time, string const & filePath ) const;

private:
  /// PVD XML file
  xmlWrapper::xmlDocument m_pvdFile;

  /// Name of the XML File
  string const m_fileName;
};
}  // namespace vtk
}  // namespace geosx

#endif
