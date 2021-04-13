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
 * @file OutputUtilities.hpp
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_OUTPUTUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_UTILITIES_OUTPUTUTILITIES_HPP_

#include "meshUtilities/cornerPointMesh/CornerPointMeshData.hpp"

namespace geosx
{

namespace cornerPointMesh
{

namespace outputUtilities
{

/**
 * @brief Debug function that outputs a vtk file representing the new mesh
 * @param[in] vertices the struct holding vertex information
 * @param[in] faces the struct holding face information
 * @param[in] cells the struct holding cell information
 */
void outputDebugVTKFile( CornerPointMeshVertices const & vertices,
                         CornerPointMeshFaces const & faces,
                         CornerPointMeshCells const & cells )
{
  GEOSX_UNUSED_VAR( faces ); // not needed for now

  std::ofstream myfile;
  myfile.open ( "debug.vtk" );
  myfile << "# vtk DataFile Version 3.0\n";
  myfile << "debug mesh\n";
  myfile << "ASCII\n";
  myfile << "DATASET UNSTRUCTURED_GRID\n";

  array2d< real64 > const & vertexPositions = vertices.m_vertexPositions;
  myfile << "POINTS " << vertexPositions.size() << " float\n";
  for( localIndex iVertex = 0; iVertex < vertexPositions.size( 0 ); ++iVertex )
  {
    myfile << vertexPositions( iVertex, 0 ) << " " << vertexPositions( iVertex, 1 ) << " " << vertexPositions( iVertex, 2 ) << "\n";
  }

  array1d< localIndex > const & activeCellToCell = cells.m_activeCellToCell;
  array1d< localIndex > const & cellToCPVertices = cells.m_cellToCPVertices;
  array1d< localIndex > const & cpVertexToVertex = vertices.m_cpVertexToVertex;
  localIndex const nActiveCells = activeCellToCell.size();

  myfile << "CELLS " << nActiveCells << " " << 9*nActiveCells << "\n";
  for( localIndex iActiveCell = 0; iActiveCell < nActiveCells; ++iActiveCell )
  {
    localIndex const iLocalCell = activeCellToCell( iActiveCell );
    localIndex const iFirstCPVertex = cellToCPVertices( iLocalCell );
    myfile << "8 ";

    localIndex const order[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };
    for( localIndex pos = 0; pos < 8; ++pos )
    {
      myfile << cpVertexToVertex( iFirstCPVertex + order[pos] ) << " ";
    }
    myfile << "\n";
  }

  myfile << "CELL_TYPES " << nActiveCells << std::endl;
  for( localIndex iActiveCell = 0; iActiveCell < nActiveCells; ++iActiveCell )
  {
    myfile << "12\n";
  }
  myfile.close();
}


} // namespace outputUtilities

} // namespace cornerPointMesh

} // namespace geosx

#endif //GEOSX_MESHUTILITIES_CORNERPOINTMESH_OUTPUTUTILITIES_HPP_
