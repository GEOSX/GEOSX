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
 * @file VTKOutput.hpp
 */

#ifndef GEOSX_MANAGERS_OUTPUTS_VTKOUTPUT_HPP_
#define GEOSX_MANAGERS_OUTPUTS_VTKOUTPUT_HPP_

#include "OutputBase.hpp"
#include "fileIO/vtk/VTKPolyDataWriterInterface.hpp"

namespace geosx
{
/**
 * @class VTKOutput
 *
 * A class for creating vtk outputs
 */
class VTKOutput : public OutputBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group(std::string const & name, Group * const parent)

  VTKOutput( std::string const & name, Group * const parent );

  /// Destructor
  virtual ~VTKOutput() override;

  /**
   * @brief Catalog name interface
   * @return This type's catalog name
   */
  static string
  CatalogName()
  {
    return "VTK";
  }

  /**
   * @brief Writes out a set of vtk files.
   * @copydoc EventBase::Execute()
   */
  virtual void
  Execute( real64 const time_n,
           real64 const dt,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           dataRepository::Group * domain ) override;

  /**
   * @brief Write one final set of vtk files as the code exits
   * @copydoc ExecutableGroup::Cleanup()
   */
  virtual void
  Cleanup( real64 const time_n,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           dataRepository::Group * domain ) override
  {
    Execute( time_n, 0, cycleNumber, eventCounter, eventProgress, domain );
  }

  /// @cond DO_NOT_DOCUMENT
  struct viewKeysStruct : OutputBase::viewKeysStruct
  {
    static constexpr auto plotFileRoot = "plotFileRoot";
    static constexpr auto writeFEMFaces = "writeFEMFaces";
    static constexpr auto plotLevel = "plotLevel";
    static constexpr auto binaryString = "writeBinaryData";

  } vtkOutputViewKeys;
  /// @endcond

private:
  string m_plotFileRoot;
  integer m_writeFaceMesh;
  integer m_plotLevel;

  integer m_writeBinaryData;

  vtk::VTKPolyDataWriterInterface m_writer;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_OUTPUTS_VTKOUTPUT_HPP_ */
