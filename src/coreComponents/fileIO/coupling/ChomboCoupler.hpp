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

#ifndef GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_
#define GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_

#include "common/DataTypes.hpp"
#include "mesh/MeshLevel.hpp"

#include <string>

namespace geosx
{
/**
 * @class ChomboCoupler
 *
 * A class managing data exchange with CHOMBO.
 */

class ChomboCoupler
{
public:
  /**
   * @brief Main Constructor.
   *
   * @param comm communicator used in reading/writing from/to file.
   * @param outputPath filename the name of the file to write out to
   * @param inputPath filename the name of the file to read from
   * @param mesh refererence to the mesh
   */
  ChomboCoupler( MPI_Comm const comm,
                 const std::string & outputPath,
                 const std::string & inputPath,
                 MeshLevel & mesh );

  /**
   * @brief Write data to file.
   *
   * @param dt the current time step.
   */
  void
  write( double dt );

  /**
   * @brief Read data from file.
   *
   * @param usePressures if @p true, pressure are read in from file
   */
  void
  read( bool usePressures );

private:
  MPI_Comm const m_comm;
  std::string const m_outputPath;
  std::string const m_inputPath;
  std::int64_t m_face_offset;
  std::int64_t m_n_faces_written;
  std::int64_t m_node_offset;
  std::int64_t m_n_nodes_written;
  MeshLevel & m_mesh;
  int m_counter;
};

} /* namespace geosx */

#endif /* GEOSX_FILEIO_COUPLING_CHOMBOCOUPLER_HPP_ */
