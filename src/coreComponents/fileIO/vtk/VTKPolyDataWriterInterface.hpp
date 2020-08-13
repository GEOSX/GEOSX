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

#ifndef GEOSX_FILEIO_VTK_VTKMULTIBLOCKWRITERINTERFACE_HPP_
#define GEOSX_FILEIO_VTK_VTKMULTIBLOCKWRITERINTERFACE_HPP_

#include "common/DataTypes.hpp"

#include "managers/DomainPartition.hpp"

#include "mesh/CellElementSubRegion.hpp"

#include "dataRepository/WrapperBase.hpp"
#include "dataRepository/Wrapper.hpp"

#include "VTKPVDWriter.hpp"
#include "VTKVTMWriter.hpp"
#include "VTKGEOSXData.hpp"

#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>

namespace geosx
{
using namespace dataRepository;
namespace vtk
{
enum struct VTKOutputMode
{
  BINARY,
  ASCII
};

/*!
 * @brief Encapsulate output methods for vtk
 */
class VTKPolyDataWriterInterface
{
public:
  /*!
   * @brief Constructor
   * @param[in] outputName folder name in which all the files will be written
   */
  VTKPolyDataWriterInterface( string const & outputName );

  /*!
   * @brief Sets the plot level
   * @details All fields have an associated plot level. If it is <= to \p plotLevel,
   * the field will be output.
   * @param[in] plotLevel the limit plotlevel
   */
  void
  SetPlotLevel( integer plotLevel )
  {
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /*!
   * @brief Set the binary mode
   * @param[in] mode output mode to be set
   */
  void
  SetOutputMode( VTKOutputMode mode )
  {
    m_outputMode = mode;
  }

  /*!
   * @brief Main method of this class. Write all the files for one time step.
   * @details This method writes a .pvd file (if a previous one was created from a precedent time step,
   * it is overwritten). The .pvd file contains relative path to every .vtm files (one vtm file per time step).
   * This method triggers also the writing of a .vtm file. A .vtm file containts relative paths to blocks
   * with the following hierarchy :
   *  - CellElementRegion
   *    - CellElementRegion1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - CellElementRegion2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - ...
   *  -WellElementRegion
   *    - Well1
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   *    - Well2
   *      - rank0
   *      - rank1
   *      - rank2
   *      - ...
   * @param[in] time the time step to be written
   * @param[in] cycle the current cycle of event
   * @param[in] domain the computation domain of this rank
   */
  void
  Write( real64 time, integer cycle, DomainPartition const & domain );

private:
  /*!
   * @brief Create a folder at the given time-step \p time
   * @details the name of the folder will be the time-step. This folder
   * will contains every files concerning the time-step \p time.
   * (aka one file per ElementRegion and per rank).
   * @param[in] time the time-step
   */
  void
  CreateTimeStepSubFolder( real64 time ) const;

  /*!
   * @brief Given a time-step \p time, returns the relative path
   * to the subfolder containing the files concerning this time-step
   * @param[in] time the time-step
   * @return the relative path to the folder of the time step
   */
  string
  GetTimeStepSubFolder( real64 time ) const;

  /*!
   * @brief Writes the files for all the CellElementRegions.
   * @details There will be one file written per CellElementRegion and per rank.
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing the CellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void
  WriteCellElementRegions( real64 time,
                           ElementRegionManager const & elemManager,
                           NodeManager const & nodeManager ) const;

  /*!
   * @brief Gets the cell connectivities as
   * a VTK object for the CellElementRegion \p er
   * @param[in] er the CellElementRegion to be written
   * @return a pair, first value is a table with the same size than the total number of element in the CellElementRegion
   * containg the type of the cells.
   * the second value is a VTK object containing the connectivity information
   * */
  std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
  GetVTKCells(
    CellElementRegion const & er ) const;

  /*!
   * @brief Gets the vertices coordinates
   * as a VTK Object for \p nodeManager
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  vtkSmartPointer< vtkPoints >
  GetVTKPoints( NodeManager const & nodeManager ) const;

  /*!
   * @brief Writes the files containing the well representation
   * @details There will be one file written per WellElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing the WellElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void
  WriteWellElementRegions( real64 time,
                           ElementRegionManager const & elemManager,
                           NodeManager const & nodeManager ) const;

  /*!
   * @brief Gets the cell connectivities and the vertices coordinates
   * as VTK objects for a specific WellElementSubRegion
   * @param[in] esr the WellElementSubRegion to be output
   * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
   * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
   * and a VTKCellArray (with the cell connectivities).
   */
  std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
  GetWell(
    WellElementSubRegion const & esr,
    NodeManager const & nodeManager ) const;

  /*!
   * @brief Writes the files containing the faces elements
   * @details There will be one file written per FaceElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing the FaceElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   */
  void
  WriteFaceElementRegions( real64 time,
                           ElementRegionManager const & elemManager,
                           NodeManager const & nodeManager ) const;

  /*!
   * @brief Gets the cell connectivities and the vertices coordinates
   * as VTK objects for a specific FaceElementSubRegion
   * @param[in] esr the FaceElementSubRegion to be output
   * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
   * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
   * and a VTKCellArray (with the cell connectivities).
   */
  std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
  GetSurface(
    FaceElementSubRegion const & esr,
    NodeManager const & nodeManager ) const;

  /*!
   * @brief Gets the cell connectivities and the vertices coordinates
   * as VTK objects for a specific EmbeddedSurafaceSubRegion
   * @param[in] esr the EmbeddedSurfaceSubRegion to be output
   * @param[in] elemManager the elemManager associated with the DomainPartition being written.
   * @param[in] nodeManager the NodeManager associated with the DomainPartition being written.
   * @param[in] edgeManager the edgeManager associated with the DomainPartition being written.
   * @return a pair containing a VTKPoints (with the information on the vertices and their coordinates)
   * and a VTKCellArray (with the cell connectivities).
   */
  std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
  GetEmbeddedSurface( EmbeddedSurfaceSubRegion const & esr,
                      ElementRegionManager const & elemManager,
                      NodeManager const & nodeManager,
                      EdgeManager const & edgeManager ) const;

  /*!
   * @brief Writes the files containing the faces elements
   * @details There will be one file written per FaceElementRegion and per rank
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing the FaceElementRegions to be output
   * @param[in] nodeManager the NodeManager containing the nodes of the domain to be output
   * @param[in] edgeManager the EdgeManager used to compute the embedded fracture intersection
   */
  void
  WriteEmbeddedSurfaceElementRegions( real64 time,
                                      ElementRegionManager const & elemManager,
                                      NodeManager const & nodeManager,
                                      EdgeManager const & edgeManager ) const;

  /*!
   * @brief Writes a VTM file for the time-step \p time.
   * @details a VTM file is a VTK Multiblock file. It contains reltive path to different files organized in blocks.
   * @param[in] time the time-step
   * @param[in] elemManager the ElementRegionManager containing all the regions to be output and refered in the VTM file
   * @param[in] vtmWrite a writer specialized for the VTM file format
   */
  void
  WriteVTMFile( real64 time,
                ElementRegionManager const & elemManager,
                VTKVTMWriter const & vtmWriter ) const;

  /*!
   * @brief Write all the fields associated to the nodes of \p nodeManager if their plotlevel is <= m_plotLevel
   * @param[in] pointdata a VTK object containing all the fields associated with the nodes
   * @param[in] nodeManager the NodeManager associated with the domain being written
   */
  void
  WriteNodeFields( vtkSmartPointer< vtkPointData > const pointdata,
                   NodeManager const & nodeManager ) const;

  /*!
   * @brief Writes all the fields associated to the elements of \p er if their plotlevel is <= m_plotLevel
   * @param[in] celldata a VTK object containing all the fields associated with the elements
   * @param[in] er ElementRegion being written
   */
  template< class SUBREGION >
  void
  WriteElementFields( vtkSmartPointer< vtkCellData > const celldata,
                      ElementRegionBase const & er ) const;

  /*!
   * @brief Writes a field from \p wrapperBase
   * @details Sets the number of components, the number of value and fill the VTK data structure using
   * a wrapper around a field.
   * @param[in] wrapperBase a wrapper around the field to be written
   * @param[in,out] data a VTK data container derived to be suitable for some GEOSX types.
   * @param[in] size the number of values in the field
   * @param[in,out] a counter that is incremented each time a value is written. This is useful
   * for CellElementSubRegion.
   */
  void
  WriteField( WrapperBase const & wrapperBase,
              vtkSmartPointer< VTKGEOSXData > data,
              localIndex size,
              localIndex & count ) const;

  /*!
   * @brief Writes an unstructured grid
   * @details The unstructured grid is the last element in the hiearchy of the output,
   * it contains the cells connectivities and the vertices coordinates as long as the
   * data fields associated with it
   * @param[in] ug a VTK SmartPointer to the VTK unstructured grid.
   * @param[in] time the current time-step
   * @param[in] name the name of the ElementRegionBase to be written
   */
  void
  WriteUnstructuredGrid( vtkSmartPointer< vtkUnstructuredGrid > ug,
                         double time,
                         string const & name ) const;

private:
  /// Folder name in which all the files will be written
  string const m_outputFolder;

  /// A writter specialized for PVD files. There is one PVD file per simulation. It is the root
  /// file containing all the paths to the VTM files.
  VTKPVDWriter m_pvd;

  /// Maximum plot level to be written.
  dataRepository::PlotLevel m_plotLevel;

  /// The previousCycle
  integer m_previousCycle;

  /// Output mode, could be ASCII or BINARAY
  VTKOutputMode m_outputMode;
};

}  // namespace vtk
}  // namespace geosx

#endif
