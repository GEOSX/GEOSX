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
 * @file CornerPointMeshBuilder.cpp
 */

#include "CornerPointMeshBuilder.hpp"

#include "codingUtilities/Utilities.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "meshUtilities/cornerPointMesh/utilities/GeometryUtilities.hpp"
#include "meshUtilities/cornerPointMesh/utilities/OutputUtilities.hpp"

namespace geosx
{

namespace cornerPointMesh
{

using namespace geometryUtilities;
using namespace outputUtilities;

CornerPointMeshBuilder::CornerPointMeshBuilder( string const & name )
  :
  m_parser( name ),
  m_partition( name ),
  m_meshName( name )
{}

void CornerPointMeshBuilder::buildMesh( Path const & filePath )
{
  localIndex const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  if( rank == 0 )
  {
    localIndex nX = 0;
    localIndex nY = 0;
    localIndex nZ = 0;

    // read SPECGRID/DIMENS to get the number of cells in each direction
    m_parser.readNumberOfCells( filePath, nX, nY, nZ );
    m_dims.defineDomainDimensions( nX, nY, nZ );

    // TODO: read ACTNUM for the inertial partitioning
  }

  // rank 0 partitions the mesh and sends partitions to other ranks
  // rank 0 sets up the communication pattern
  m_partition.setupMPIPartition( m_dims );

  // each rank reads its part of the mesh
  m_parser.readMesh( filePath, m_dims );

  // post-process the mesh to make it conforming
  postProcessMesh();
}

void CornerPointMeshBuilder::postProcessMesh()
{
  // for each active cell, compute the position of the eight "corner-point" vertices
  // at this point, "corner-point" vertices have not been filtered yet
  buildCornerPointCells();

  // eliminate duplicates from cpVertices
  filterVertices();

  // match faces, deal with the non-conforming case
  //buildFaces();

  // construct map linking a region to its cells
  formRegions();

  // for now, do some debugging
  outputDebugVTKFile( m_vertices, m_faces, m_cells );
}

void CornerPointMeshBuilder::buildCornerPointCells()
{
  // local and global dimensions
  localIndex const nX = m_dims.nX();
  localIndex const nY = m_dims.nY();

  localIndex const iMinLocal = m_dims.iMinLocal();
  localIndex const jMinLocal = m_dims.jMinLocal();
  localIndex const kMinLocal = 0;

  localIndex const iMinOverlap = m_dims.iMinOverlap();
  localIndex const jMinOverlap = m_dims.jMinOverlap();
  localIndex const iMaxOverlap = m_dims.iMaxOverlap();
  localIndex const jMaxOverlap = m_dims.jMaxOverlap();

  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal();
  localIndex const nLocalCells = nXLocal*nYLocal*nZLocal;
  localIndex constexpr nCPVerticesPerCell = 8;
  localIndex const nCPVertices = nCPVerticesPerCell*nLocalCells;

  // corner-point mesh file information
  array1d< localIndex > const & actnum = m_parser.actnum();
  array1d< real64 > const & coord = m_parser.coord();
  array1d< real64 > const & zcorn = m_parser.zcorn();

  // vertex maps
  array2d< real64 > & cpVertexPositions = m_vertices.m_cpVertexPositions;
  array1d< globalIndex > & cpVertexToGlobalCPVertex = m_vertices.m_cpVertexToGlobalCPVertex;
  array1d< bool > & cpVertexIsInsidePartition = m_vertices.m_cpVertexIsInsidePartition;
  cpVertexPositions.resizeDimension< 0, 1 >( nCPVertices, 3 );
  cpVertexToGlobalCPVertex.resize( nCPVertices );
  cpVertexIsInsidePartition.resize( nCPVertices );

  // cell maps
  array1d< globalIndex > & ownedActiveCellToGlobalCell = m_cells.m_ownedActiveCellToGlobalCell;
  array1d< localIndex > & ownedActiveCellToActiveCell = m_cells.m_ownedActiveCellToActiveCell;
  array1d< localIndex > & activeCellToCell = m_cells.m_activeCellToCell;
  array1d< localIndex > & cellToCPVertices = m_cells.m_cellToCPVertices;
  ownedActiveCellToGlobalCell.reserve( nLocalCells );
  ownedActiveCellToActiveCell.reserve( nLocalCells );
  activeCellToCell.reserve( nLocalCells );
  cellToCPVertices.resize( nLocalCells );

  array1d< real64 > xPos( nCPVerticesPerCell );
  array1d< real64 > yPos( nCPVerticesPerCell );
  array1d< real64 > zPos( nCPVerticesPerCell );

  array1d< bool > cpVertexIsInside( nCPVerticesPerCell );

  // loop over all cells in the MPI domain, including overlaps
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        bool isValid = false;

        // compute explicit local and global indices using the ijk structure
        localIndex const iLocalCell = k*nXLocal*nYLocal + j*nXLocal + i;
        globalIndex const iGlobalCell = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + (i+iMinLocal);

        if( actnum( iLocalCell ) == 1 )
        {
          // compute the positions of the eight vertices
          isValid =
            computeCellCoordinates( i, j, k, nXLocal, nYLocal,
                                    iMinOverlap, iMaxOverlap, jMinOverlap, jMaxOverlap,
                                    coord, zcorn,
                                    xPos, yPos, zPos, cpVertexIsInside );
        }

        if( isValid )
        {
          // cell maps

          // assign local and global indices
          activeCellToCell.emplace_back( iLocalCell );
          // decide flag specifying whether the cell is in the overlap or not
          if( m_dims.columnIsInsidePartition( i, j ) )
          {
            ownedActiveCellToActiveCell.emplace_back( activeCellToCell.size()-1 );
            ownedActiveCellToGlobalCell.emplace_back( iGlobalCell );
          }

          // vertex maps

          // construct the map from local cell to first CP vertex of the cell
          localIndex const iFirstVertexLocal = nCPVerticesPerCell * iLocalCell;
          globalIndex const iFirstVertexGlobal = nCPVerticesPerCell * iGlobalCell;
          cellToCPVertices( iLocalCell ) = iFirstVertexLocal;

          // save the position of the eight vertices, and assign global CP vertex indices
          localIndex const order[nCPVerticesPerCell] = { 4, 5, 6, 7, 0, 1, 2, 3 };
          for( localIndex pos = 0; pos < nCPVerticesPerCell; ++pos )
          {
            cpVertexPositions( iFirstVertexLocal + pos, 0 ) = xPos( order[pos] );
            cpVertexPositions( iFirstVertexLocal + pos, 1 ) = yPos( order[pos] );
            cpVertexPositions( iFirstVertexLocal + pos, 2 ) = -zPos( order[pos] );
            cpVertexToGlobalCPVertex( iFirstVertexLocal + pos ) = iFirstVertexGlobal + order[pos];
            cpVertexIsInsidePartition( iFirstVertexLocal + pos ) = cpVertexIsInside( order[pos] );
          }
        }
      }
    }
  }
  populateReverseCellMaps( nLocalCells, activeCellToCell.size(),
                           activeCellToCell, ownedActiveCellToActiveCell,
                           m_cells.m_cellToActiveCell, m_cells.m_activeCellToOwnedActiveCell );
}

void CornerPointMeshBuilder::populateReverseCellMaps( localIndex const nCells,
                                                      localIndex const nActiveCells,
                                                      array1d< localIndex > const & activeCellToCell,
                                                      array1d< localIndex > const & ownedActiveCellToActiveCell,
                                                      array1d< localIndex > & cellToActiveCell,
                                                      array1d< localIndex > & activeCellToOwnedActiveCell )
{
  cellToActiveCell.resize( nCells );
  activeCellToOwnedActiveCell.resize( nActiveCells );
  cellToActiveCell.setValues< serialPolicy >( -1 );
  activeCellToOwnedActiveCell.setValues< serialPolicy >( -1 );

  for( localIndex iActiveCell = 0; iActiveCell < activeCellToCell.size(); ++iActiveCell )
  {
    cellToActiveCell( activeCellToCell( iActiveCell ) ) = iActiveCell;
  }
  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < ownedActiveCellToActiveCell.size(); ++iOwnedActiveCell )
  {
    activeCellToOwnedActiveCell( ownedActiveCellToActiveCell( iOwnedActiveCell ) ) = iOwnedActiveCell;
  }
}

bool CornerPointMeshBuilder::computeCellCoordinates( localIndex const i, localIndex const j, localIndex const k,
                                                     localIndex const nXLocal, localIndex const nYLocal,
                                                     localIndex const iMinOverlap, localIndex const iMaxOverlap,
                                                     localIndex const jMinOverlap, localIndex const jMaxOverlap,
                                                     array1d< real64 > const & coord,
                                                     array1d< real64 > const & zcorn,
                                                     array1d< real64 > & xPos,
                                                     array1d< real64 > & yPos,
                                                     array1d< real64 > & zPos,
                                                     array1d< bool > & cpVertexIsInside )
{
  localIndex const iXmLocal = k*8*nXLocal*nYLocal + j*4*nXLocal + 2*i;
  localIndex const iXpLocal = iXmLocal + 2*nXLocal;

  zPos( 0 ) = zcorn( iXmLocal );
  zPos( 1 ) = zcorn( iXmLocal + 1 );
  zPos( 2 ) = zcorn( iXpLocal );
  zPos( 3 ) = zcorn( iXpLocal + 1 );
  zPos( 4 ) = zcorn( iXmLocal + 4*nXLocal*nYLocal );
  zPos( 5 ) = zcorn( iXmLocal + 4*nXLocal*nYLocal + 1 );
  zPos( 6 ) = zcorn( iXpLocal + 4*nXLocal*nYLocal );
  zPos( 7 ) = zcorn( iXpLocal + 4*nXLocal*nYLocal + 1 );

  bool const hexaIsFlat = isZero( zPos( 0 )-zPos( 4 )
                                  + zPos( 1 )-zPos( 5 )
                                  + zPos( 2 )-zPos( 6 )
                                  + zPos( 3 )-zPos( 7 ) );

  if( hexaIsFlat )
  {
    return false;
  }

  // the following code comes from PAMELA
  // TODO: if possible, remove code duplication

  // first pillar
  localIndex const iFirstPillarLocal = j*(nXLocal+1)+i;
  localIndex const ip0 = 6*iFirstPillarLocal - 1;
  real64 const denomFirstPillar = coord( ip0 + 6 ) - coord( ip0 + 3 );
  bool const firstPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos0 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 0 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 0 ) = slopePos0 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 0 ) = slopePos0 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cpVertexIsInside( 0 ) = firstPillarIsInside;

  real64 const slopePos4 = isZero( denomFirstPillar )
    ? 1.0
    : ( zPos( 4 )-coord( ip0 + 3 ) ) / denomFirstPillar;
  xPos( 4 ) = slopePos4 * (coord( ip0 + 4 ) - coord( ip0 + 1 )) + coord( ip0 + 1 );
  yPos( 4 ) = slopePos4 * (coord( ip0 + 5 ) - coord( ip0 + 2 )) + coord( ip0 + 2 );
  cpVertexIsInside( 4 ) = firstPillarIsInside;

  // second pillar
  localIndex const iSecondPillarLocal = (nXLocal+1)*j + i+1;
  localIndex const ip1 = 6*iSecondPillarLocal - 1;
  real64 const denomSecondPillar = coord( ip1 + 6 ) - coord( ip1 + 3 );
  bool const secondPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j == 0 && jMinOverlap == 1) );

  real64 const slopePos1 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 1 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 1 ) = slopePos1 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 1 ) = slopePos1 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cpVertexIsInside( 1 ) = secondPillarIsInside;

  real64 const slopePos5 = isZero( denomSecondPillar )
    ? 1.0
    : ( zPos( 5 )-coord( ip1 + 3 ) ) / denomSecondPillar;
  xPos( 5 ) = slopePos5 * (coord( ip1 + 4 ) - coord( ip1 + 1 )) + coord( ip1 + 1 );
  yPos( 5 ) = slopePos5 * (coord( ip1 + 5 ) - coord( ip1 + 2 )) + coord( ip1 + 2 );
  cpVertexIsInside( 5 ) = secondPillarIsInside;

  // third pillar
  localIndex const iThirdPillarLocal = (nXLocal+1)*(j+1) + i;
  localIndex const ip2 = 6*iThirdPillarLocal - 1;
  real64 const denomThirdPillar = coord( ip2 + 6 ) - coord( ip2 + 3 );
  bool const thirdPillarIsInside = !( (i == 0 && iMinOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos2 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 2 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 2 ) = slopePos2 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 2 ) = slopePos2 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cpVertexIsInside( 2 ) = thirdPillarIsInside;

  real64 const slopePos6 = isZero( denomThirdPillar )
    ? 1.0
    : ( zPos( 6 )-coord( ip2 + 3 ) ) / denomThirdPillar;
  xPos( 6 ) = slopePos6 * (coord( ip2 + 4 ) - coord( ip2 + 1 )) + coord( ip2 + 1 );
  yPos( 6 ) = slopePos6 * (coord( ip2 + 5 ) - coord( ip2 + 2 )) + coord( ip2 + 2 );
  cpVertexIsInside( 6 ) = thirdPillarIsInside;

  // fourth pillar
  localIndex const iFourthPillarLocal = (nXLocal+1)*(j+1) + i+1;
  localIndex const ip3 = 6*iFourthPillarLocal - 1;
  real64 const denomFourthPillar = coord( ip3 + 6 ) - coord( ip3 + 3 );
  bool const fourthPillarIsInside = !( (i+1 == nXLocal && iMaxOverlap == 1) || (j+1 == nYLocal && jMaxOverlap == 1) );

  real64 const slopePos3 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 3 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 3 ) = slopePos3 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 3 ) = slopePos3 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cpVertexIsInside( 3 ) = fourthPillarIsInside;

  real64 const slopePos7 = isZero( denomFourthPillar )
    ? 1.0
    : ( zPos( 7 )-coord( ip3 + 3 ) ) / denomFourthPillar;
  xPos( 7 ) = slopePos7 * (coord( ip3 + 4 ) - coord( ip3 + 1 )) + coord( ip3 + 1 );
  yPos( 7 ) = slopePos7 * (coord( ip3 + 5 ) - coord( ip3 + 2 )) + coord( ip3 + 2 );
  cpVertexIsInside( 7 ) = fourthPillarIsInside;

  return true;
}


void CornerPointMeshBuilder::formRegions()
{
  array1d< localIndex > const & activeCellToCell = m_cells.m_activeCellToCell;
  array1d< localIndex > const & ownedActiveCellToActiveCell = m_cells.m_ownedActiveCellToActiveCell;
  array1d< localIndex > const & regionId = m_parser.regionId();

  // Step 1: find the largest region index
  localIndex const nOwnedActiveCells = ownedActiveCellToActiveCell.size();
  localIndex localMaxRegionId = 0;
  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < nOwnedActiveCells; ++iOwnedActiveCell )
  {
    localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
    localIndex const iCell = activeCellToCell( iActiveCell );
    localIndex const rid = regionId( iCell );
    if( rid > localMaxRegionId )
    {
      localMaxRegionId = rid;
    }
  }
  localIndex const globalMaxRegionId = MpiWrapper::max( localMaxRegionId );

  // Step 2: construct the map region to elements
  ArrayOfArrays< localIndex > & regionToOwnedActiveCell = m_regions.m_regionToOwnedActiveCell;
  localIndex const regionCapacity =
    static_cast< real64 >( nOwnedActiveCells )
    / static_cast< real64 >( globalMaxRegionId+1 );
  regionToOwnedActiveCell.resize( globalMaxRegionId+1, regionCapacity );

  for( localIndex iOwnedActiveCell = 0; iOwnedActiveCell < nOwnedActiveCells; ++iOwnedActiveCell )
  {
    localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
    localIndex const iCell = activeCellToCell( iActiveCell );

    regionToOwnedActiveCell.emplaceBack( regionId( iCell ), iOwnedActiveCell );
  }
}


void CornerPointMeshBuilder::filterVertices()
{
  array2d< real64 > const & cpVertexPositions = m_vertices.m_cpVertexPositions;
  array1d< bool > const & cpVertexIsInsidePartition = m_vertices.m_cpVertexIsInsidePartition;
  array1d< globalIndex > const & cpVertexToGlobalCPVertex = m_vertices.m_cpVertexToGlobalCPVertex;
  array1d< localIndex > & cpVertexToVertex = m_vertices.m_cpVertexToVertex;
  cpVertexToVertex.resize( cpVertexPositions.size() );

  // First step: filter the unique vertices using a set

  std::set< Vertex, CompareVertices > uniqueVerticesHelper;
  std::set< Vertex, CompareVertices >::iterator it;

  // loop over of the CP vertices (for now, including those of inactive cells)
  for( localIndex iCPVertex = 0; iCPVertex < cpVertexPositions.size( 0 ); ++iCPVertex )
  {
    // make sure cpVertices in the **exterior** pillars of the overlap are skipped
    if( cpVertexIsInsidePartition( iCPVertex ) )
    {
      Vertex v( cpVertexPositions( iCPVertex, 0 ),
                cpVertexPositions( iCPVertex, 1 ),
                cpVertexPositions( iCPVertex, 2 ) );

      // check if this vertex has already been found and inserted in the set
      it = uniqueVerticesHelper.find( v );

      // if already found, copy the local index in map and move on
      if( it != uniqueVerticesHelper.end() )
      {
        Vertex const & existingVertex = *it;
        cpVertexToVertex( iCPVertex ) = existingVertex.m_localIndex;
        if( existingVertex.m_globalIndex > cpVertexToGlobalCPVertex( iCPVertex ) )
        {
          // make sure that the global index does not depend on the order in which cpVertices are processed
          // unfortunately I cannot directly modify existingVertex.m_globalIndex (it is const)
          // TODO: find out is there is a better way to do that
          v.m_localIndex = existingVertex.m_localIndex;
          v.m_globalIndex = cpVertexToGlobalCPVertex( iCPVertex );
          uniqueVerticesHelper.erase( existingVertex );
          uniqueVerticesHelper.insert( v );
        }
      }
      // if not found yet, insert into the set
      else
      {
        v.m_localIndex = uniqueVerticesHelper.size();
        v.m_globalIndex = cpVertexToGlobalCPVertex( iCPVertex );
        cpVertexToVertex( iCPVertex ) = v.m_localIndex;
        uniqueVerticesHelper.insert( v );
      }
    }
  }

  // Second step: move the data from the set to a array2d, because we don't want the set anymore
  // TODO: check if we really need this second step
  array2d< real64 > & vertexPositions = m_vertices.m_vertexPositions;
  array1d< globalIndex > & vertexToGlobalVertex = m_vertices.m_vertexToGlobalVertex;
  vertexPositions.resizeDimension< 0, 1 >( uniqueVerticesHelper.size(), 3 );
  vertexToGlobalVertex.resize( uniqueVerticesHelper.size() );
  std::for_each( uniqueVerticesHelper.begin(), uniqueVerticesHelper.end(),
                 [&vertexPositions, &vertexToGlobalVertex]( Vertex const & v )
  {
    vertexPositions( v.m_localIndex, 0 ) = v.m_x;
    vertexPositions( v.m_localIndex, 1 ) = v.m_y;
    vertexPositions( v.m_localIndex, 2 ) = v.m_z;
    vertexToGlobalVertex( v.m_localIndex ) = v.m_globalIndex;
  } );
}


void CornerPointMeshBuilder::buildFaces()
{
  // not ready for review, in construction

  // allocate space for the maps
  localIndex const nOwnedActiveCells = m_cells.m_ownedActiveCellToActiveCell.size();
  localIndex const faceCapacity = 4;

  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;
  ownedActiveCellToFaces.resize( nOwnedActiveCells, faceCapacity );
  // TODO: reserve space for faceToNodes

  // TODO: the construction of faces involves loops in the (j,i,k) order,
  //       which is not ideal giving the layout of the arrays used here.
  //       See if it is worth permuting the arrays before this step.

  // Step 1: construct faces in the X and Y directions (between columns of cells)
  buildVerticalFaces( "X" );
  buildVerticalFaces( "Y" );

  // Step 2: construct faces in the Z direction (in a column of cells)
  buildHorizontalFaces();
}

void CornerPointMeshBuilder::buildVerticalFaces( string const & direction )
{
  // not ready for review, in construction
  GEOSX_UNUSED_VAR( direction );
}

void CornerPointMeshBuilder::buildHorizontalFaces()
{
  // not ready for review, in construction
  localIndex const nXLocal = m_dims.nXLocal();
  localIndex const nYLocal = m_dims.nYLocal();
  localIndex const nZLocal = m_dims.nZLocal();

  // cells
  array1d< localIndex > const & cellToCPVertices = m_cells.m_cellToCPVertices;
  array1d< localIndex > const & cellToActiveCell = m_cells.m_cellToActiveCell;
  array1d< localIndex > const & activeCellToOwnedActiveCell = m_cells.m_activeCellToOwnedActiveCell;
  ArrayOfArrays< localIndex > & ownedActiveCellToFaces = m_cells.m_ownedActiveCellToFaces;

  // faces
  ArrayOfArrays< localIndex > & faceToVertices = m_faces.m_faceToVertices;

  // vertices
  array1d< localIndex > const & cpVertexToVertex = m_vertices.m_cpVertexToVertex;

  // loop over owned cells only
  for( localIndex j = 0; j < nYLocal; ++j )
  {
    for( localIndex i = 0; i < nXLocal; ++i )
    {
      if( m_dims.columnIsInsidePartition( i, j ) )
      {

        // we are in column of cells, we can now loop over all the faces in this column
        localIndex const nFacesInColumn = nZLocal+1;
        for( localIndex iface = 0; iface < nFacesInColumn; ++iface )
        {
          // k-index of the cell before the face
          localIndex const kPrev = iface-1;
          // k-index of the cell after the face
          localIndex const kNext = iface;

          // local index of the previous cell
          localIndex iCellPrev = -1;
          // local index of the next cell
          localIndex iCellNext = -1;
          // local index of the previous active cell
          localIndex iActiveCellPrev = -1;
          // local index of the next active cell
          localIndex iActiveCellNext = -1;

          bool const onTopBoundary = (kPrev == -1);
          bool const onBottomBoundary = (kNext == nZLocal);
          bool const isInInterior = !onTopBoundary && !onBottomBoundary;
          if( !onTopBoundary )
          {
            iCellPrev = kPrev*nXLocal*nYLocal + j*nXLocal + i;
            iActiveCellPrev = cellToActiveCell( iCellPrev );
          }
          if( !onBottomBoundary )
          {
            iCellNext = kNext*nXLocal*nYLocal + j*nXLocal + i;
            iActiveCellNext = cellToActiveCell( iCellNext );
          }

          bool const prevIsActive = (iActiveCellPrev != -1);
          bool const nextIsActive = (iActiveCellNext != -1);

          // this is an interior cell between active cells
          if( isInInterior && prevIsActive && nextIsActive )
          {
            // get the bottom vertices of the previous face
            localIndex const iFirstVertexPrev = cellToCPVertices( iCellPrev );
            localIndex const vertices[4] = { cpVertexToVertex( iFirstVertexPrev ),
                                             cpVertexToVertex( iFirstVertexPrev + 1 ),
                                             cpVertexToVertex( iFirstVertexPrev + 2 ),
                                             cpVertexToVertex( iFirstVertexPrev + 3 ) };
            addHorizontalFace( vertices,
                               activeCellToOwnedActiveCell( iActiveCellPrev ),
                               activeCellToOwnedActiveCell( iActiveCellNext ),
                               ownedActiveCellToFaces,
                               faceToVertices );
          }
          // this is either: an interior face between an active cell and an inactive cell
          //             or: an boundary face at the bottom
          else if( prevIsActive || nextIsActive )
          {
            // get the bottom vertices of the previous face
            localIndex const iFirstVertex = prevIsActive ? cellToCPVertices( iCellPrev ) : cellToCPVertices( iCellNext );
            localIndex const iActiveCell = prevIsActive ? iActiveCellPrev : iActiveCellNext;
            localIndex const firstPos = prevIsActive ? 0 : 4;
            localIndex const vertices[4] = { cpVertexToVertex( iFirstVertex + firstPos ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 1 ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 2 ),
                                             cpVertexToVertex( iFirstVertex + firstPos + 3 ) };
            addHorizontalFace( vertices,
                               activeCellToOwnedActiveCell( iActiveCell ),
                               -1,
                               ownedActiveCellToFaces,
                               faceToVertices );
          }
        }
      }
    }
  }
}

void CornerPointMeshBuilder::addHorizontalFace( localIndex const (&faceVertices)[ 4 ],
                                                localIndex const iOwnedActiveCellPrev,
                                                localIndex const iOwnedActiveCellNext,
                                                ArrayOfArrays< localIndex > & ownedActiveCellToFaces,
                                                ArrayOfArrays< localIndex > & faceToVertices )
{
  // populate the face-to-vertices map
  localIndex const newFaceId = faceToVertices.size();
  faceToVertices.appendArray( 4 );
  faceToVertices( newFaceId, 0 ) = faceVertices[ 0 ];
  faceToVertices( newFaceId, 1 ) = faceVertices[ 1 ];
  faceToVertices( newFaceId, 2 ) = faceVertices[ 2 ];
  faceToVertices( newFaceId, 3 ) = faceVertices[ 3 ];

  // populate the cell-to-face map
  if( iOwnedActiveCellPrev != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellPrev, newFaceId );
  }
  if( iOwnedActiveCellNext != -1 )
  {
    ownedActiveCellToFaces.emplaceBack( iOwnedActiveCellNext, newFaceId );
  }
}

REGISTER_CATALOG_ENTRY( CornerPointMeshBuilder, CornerPointMeshBuilder, string const & )

} // namespace cornerPointMesh

} // namespace geosx
