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
 * @file CornerPointMeshGenerator.cpp
 */

#include "CornerPointMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace cornerPointMesh;

CornerPointMeshGenerator::CornerPointMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_permeabilityUnitInInputFile( PermeabilityUnit::Millidarcy ),
  m_coordinatesUnitInInputFile( CoordinatesUnit::Meter ),
  m_toSquareMeter( 1.0 ),
  m_toMeter( 1.0 )
{
  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "Path to the mesh file" );

  registerWrapper( viewKeyStruct::permeabilityUnitInInputFileString(), &m_permeabilityUnitInInputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( PermeabilityUnit::Millidarcy ).
    setDescription( "Flag to specify the unit of permeability in the input file. \n"
                    "Two options are available: Millidarcy (default) or SquareMeter. \n"
                    "If Millidarcy is chosen, a conversion factor is applied as GEOSX internally works with square meters." );

  registerWrapper( viewKeyStruct::coordinatesUnitInInputFileString(), &m_coordinatesUnitInInputFile ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( CoordinatesUnit::Meter ).
    setDescription( "Flag to specify the unit of the coordinates in the input file. \n"
                    "Two options are available: Meter (default) or Foot. \n"
                    "If Foot is chosen, a conversion factor is applied as GEOSX internally works with meters." );

  string const builderName = "cpMeshBuilder";
  m_cpMeshBuilder = std::make_unique< cornerPointMesh::CornerPointMeshBuilder >( builderName );
}

CornerPointMeshGenerator::~CornerPointMeshGenerator()
{}

void CornerPointMeshGenerator::postProcessInput()
{
  if( m_permeabilityUnitInInputFile == PermeabilityUnit::SquareMeter )
  {
    m_toSquareMeter = 1.0;
  }
  else if( m_permeabilityUnitInInputFile == PermeabilityUnit::Millidarcy )
  {
    m_toSquareMeter = 9.869232667160128e-16;
  }
  else
  {
    GEOSX_THROW( "This unit for permeability is not supported", InputError );
  }

  if( m_coordinatesUnitInInputFile == CoordinatesUnit::Meter )
  {
    m_toMeter = 1.0;
  }
  else if( m_coordinatesUnitInInputFile == CoordinatesUnit::Foot )
  {
    m_toMeter = 0.3048;
  }
  else
  {
    GEOSX_THROW( "This unit for coordinates is not supported", InputError );
  }


  m_cpMeshBuilder->buildMesh( m_filePath );
}

Group * CornerPointMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                               string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void CornerPointMeshGenerator::generateMesh( DomainPartition & domain )
{
  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  // Step 0: transfer the neighbor list

  // TODO: change the name, we don't use metis
  domain.getMetisNeighborList() = m_cpMeshBuilder->neighborsList();

  // we can start constructing the mesh

  // Step 1: fill vertex information

  arrayView2d< real64 const > vertexPositions = m_cpMeshBuilder->vertexPositions();
  arrayView1d< globalIndex const > vertexToGlobalVertex = m_cpMeshBuilder->vertexToGlobalVertex();
  localIndex const nVertices = vertexPositions.size( 0 );
  nodeManager.resize( nVertices );
  arrayView1d< globalIndex > const & vertexLocalToGlobal = nodeManager.localToGlobalMap();
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  Group & vertexSets = nodeManager.sets();
  SortedArray< localIndex > & allVertices =
    vertexSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();

  real64 xMin[3] = { std::numeric_limits< real64 >::max() };
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };

  for( localIndex iVertex = 0; iVertex < nVertices; ++iVertex )
  {
    X( iVertex, 0 ) = m_toMeter * vertexPositions( iVertex, 0 );
    X( iVertex, 1 ) = m_toMeter * vertexPositions( iVertex, 1 );
    X( iVertex, 2 ) = m_toMeter * vertexPositions( iVertex, 2 );
    allVertices.insert( iVertex );
    vertexLocalToGlobal( iVertex ) = vertexToGlobalVertex( iVertex );

    for( int dim = 0; dim < 3; dim++ )
    {
      if( X( iVertex, dim ) > xMax[dim] )
      {
        xMax[dim] = X( iVertex, dim );
      }
      if( X( iVertex, dim ) < xMin[dim] )
      {
        xMin[dim] = X( iVertex, dim );
      }
    }
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( xMax ) );

  // Step 2: fill cell information for each region

  ArrayOfArraysView< localIndex const > regionId = m_cpMeshBuilder->regionId();

  for( localIndex er = 0; er < regionId.size(); ++er )
  {

    CellBlock * cellBlock = &cellBlockManager.getGroup( keys::cellBlocks ).
                              registerGroup< CellBlock >( "DEFAULT_HEX_"+std::to_string( er ) );

    if( regionId.sizeOfArray( er ) == 0 )
    {
      continue;
    }

    // temporary accessors while we only support the conforming case
    // this is not what we ultimately want, which is instead:
    //  -> map from elem to faces
    //  -> map from face to nodes
    // what is below is a temporary mess (that maybe, should be hidden in CPMeshData)
    arraySlice1d< localIndex const > ownedActiveCellsInRegion = regionId[er];
    arrayView1d< localIndex const > ownedActiveCellToActiveCell = m_cpMeshBuilder->ownedActiveCellToActiveCell();
    arrayView1d< globalIndex const > ownedActiveCellToGlobalCell = m_cpMeshBuilder->ownedActiveCellToGlobalCell();

    arrayView1d< localIndex const > activeCellToCell = m_cpMeshBuilder->activeCellToCell();
    arrayView1d< localIndex const > cellToCPVertices = m_cpMeshBuilder->cellToCPVertices();
    arrayView1d< localIndex const > cpVertexToVertex = m_cpMeshBuilder->cpVertexToVertex();

    localIndex const nOwnedActiveCellsInRegion = ownedActiveCellsInRegion.size();
    cellBlock->setElementType( "C3D8" );
    cellBlock->resize( nOwnedActiveCellsInRegion );

    arrayView1d< globalIndex > cellLocalToGlobal = cellBlock->localToGlobalMap();
    auto & cellToVertex = cellBlock->nodeList(); // TODO: remove auto
    cellToVertex.resize( nOwnedActiveCellsInRegion, 8 );

    for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
    {
      localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
      // this filtering is needed because the CPMeshBuilder uses a layer of cells around each MPI partition
      // to facilitate the treatment of non-matching faces at the boundary between MPI partitions.
      // But, it should be hidden in CPMeshBuilder/CPMeshData. I will take care of this when we handle the non-conforming case.
      localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
      localIndex const iFirstCPVertex = cellToCPVertices( activeCellToCell( iActiveCell ) );

      // temporary code while we don't support the non-conforming case
      cellToVertex( iOwnedActiveCellInRegion, 0 ) = cpVertexToVertex( iFirstCPVertex );
      cellToVertex( iOwnedActiveCellInRegion, 1 ) = cpVertexToVertex( iFirstCPVertex + 1 );
      cellToVertex( iOwnedActiveCellInRegion, 2 ) = cpVertexToVertex( iFirstCPVertex + 2 );
      cellToVertex( iOwnedActiveCellInRegion, 3 ) = cpVertexToVertex( iFirstCPVertex + 3 );
      cellToVertex( iOwnedActiveCellInRegion, 4 ) = cpVertexToVertex( iFirstCPVertex + 4 );
      cellToVertex( iOwnedActiveCellInRegion, 5 ) = cpVertexToVertex( iFirstCPVertex + 5 );
      cellToVertex( iOwnedActiveCellInRegion, 6 ) = cpVertexToVertex( iFirstCPVertex + 6 );
      cellToVertex( iOwnedActiveCellInRegion, 7 ) = cpVertexToVertex( iFirstCPVertex + 7 );

      cellLocalToGlobal( iOwnedActiveCellInRegion ) = ownedActiveCellToGlobalCell( iOwnedActiveCellInRegion );
    }

    // Step 3: fill property information

    // Step 3.a: fill porosity in active cells
    arrayView1d< real64 const > porosityField = m_cpMeshBuilder->porosityField();
    if( !porosityField.empty() )
    {
      arrayView1d< real64 > referencePorosity = cellBlock->addProperty< array1d< real64 > >( "referencePorosity" ).toView();
      for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
      {
        localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
        localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
        localIndex const iCell = activeCellToCell( iActiveCell );
        referencePorosity( iOwnedActiveCellInRegion ) = porosityField( iCell );
      }
    }

    // Step 3.b: fill permeability in active cells
    arrayView2d< real64 const > permeabilityField = m_cpMeshBuilder->permeabilityField();
    if( !permeabilityField.empty() )
    {
      array2d< real64 > & permeability = cellBlock->addProperty< array2d< real64 > >( "permeability" );
      permeability.resizeDimension< 1 >( 3 );
      for( localIndex iOwnedActiveCellInRegion = 0; iOwnedActiveCellInRegion < nOwnedActiveCellsInRegion; ++iOwnedActiveCellInRegion )
      {
        localIndex const iOwnedActiveCell = ownedActiveCellsInRegion( iOwnedActiveCellInRegion );
        localIndex const iActiveCell = ownedActiveCellToActiveCell( iOwnedActiveCell );
        localIndex const iCell = activeCellToCell( iActiveCell );
        for( localIndex dim = 0; dim < 3; dim++ )
        {
          permeability( iOwnedActiveCellInRegion, dim ) =
            LvArray::math::max( 1e-19, m_toSquareMeter * permeabilityField( iCell, dim ) );
        }
      }
    }
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, CornerPointMeshGenerator, string const &, Group * const )
}
