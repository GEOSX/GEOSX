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
 * @file PAMELAMeshGenerator.cpp
 */

#include "PAMELAMeshGenerator.hpp"

#include "Elements/Element.hpp"
#include "MeshDataWriters/Variable.hpp"
#include "mesh/DomainPartition.hpp"

#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

#include "CellBlockManager.hpp"

namespace geosx
{
using namespace dataRepository;

/// TODO when we are going to port to c++17, we can remove that
string const PAMELAMeshGenerator::DecodePAMELALabels::m_separator = "_";

PAMELAMeshGenerator::PAMELAMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{

  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the mesh file" );
  registerWrapper( viewKeyStruct::fieldsToImportString(), &m_fieldsToImport ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Fields to be imported from the external mesh file" );
  registerWrapper( viewKeyStruct::fieldNamesInGEOSXString(), &m_fieldNamesInGEOSX ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the fields within GEOSX" );
  registerWrapper( viewKeyStruct::scaleString(), &m_scale ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 1. ).setDescription( "Scale the coordinates of the vertices" );
  registerWrapper( viewKeyStruct::reverseZString(), &m_isZReverse ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDefaultValue( 0 ).setDescription( "0 : Z coordinate is upward, 1 : Z coordinate is downward" );
}

PAMELAMeshGenerator::~PAMELAMeshGenerator()
{}

void PAMELAMeshGenerator::postProcessInput()
{
  m_pamelaMesh =
    std::unique_ptr< PAMELA::Mesh >
      ( PAMELA::MeshFactory::makeMesh( m_filePath ) );
  m_pamelaMesh->CreateFacesFromCells();
  m_pamelaMesh->PerformPolyhedronPartitioning( PAMELA::ELEMENTS::FAMILY::POLYGON,
                                               PAMELA::ELEMENTS::FAMILY::POLYGON );
  m_pamelaMesh->CreateLineGroupWithAdjacency(
    "TopologicalC2C",
    m_pamelaMesh->getAdjacencySet()->get_TopologicalAdjacency( PAMELA::ELEMENTS::FAMILY::POLYHEDRON, PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
                                                               PAMELA::ELEMENTS::FAMILY::POLYGON ));
}

Group * PAMELAMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void PAMELAMeshGenerator::generateMesh( DomainPartition & domain )
{
  GEOSX_LOG_RANK_0( "Writing into the GEOSX mesh data structure" );
  domain.getMetisNeighborList() = m_pamelaMesh->getNeighborList();
  Group & meshBodies = domain.getGroup( string( "MeshBodies" ) );
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  meshBody.registerGroup< MeshLevel >( string( "Level0" ) );
  CellBlockManager & cellBlockManager = domain.registerGroup< CellBlockManager >( keys::cellManager );

  // Use the PartMap of PAMELA to get the mesh
  auto const polyhedronPartMap = std::get< 0 >( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get(), 0 ));

  // Vertices are written first
  array2d< real64, nodes::REFERENCE_POSITION_PERM > & X = cellBlockManager.getNodesPositions();
  localIndex const numNodes = m_pamelaMesh->get_PointCollection()->size_all();
  cellBlockManager.setNumNodes( numNodes );

  array1d< globalIndex > & nodeLocalToGlobal = cellBlockManager.getNodeLocalToGlobal();

  auto & nodeSets = cellBlockManager.getNodeSets();
  SortedArray< localIndex > & allNodes = nodeSets["all"];

  real64 xMax[3] = { std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max() };

  double zReverseFactor = 1.;
  if( m_isZReverse )
  {
    zReverseFactor = -1.;
  }
  for( auto const & verticesIterator : *m_pamelaMesh->get_PointCollection())
  {
    localIndex const vertexLocalIndex = verticesIterator->get_localIndex();
    globalIndex const vertexGlobalIndex = verticesIterator->get_globalIndex();
    X( vertexLocalIndex, 0 ) = verticesIterator->get_coordinates().x * m_scale;
    X( vertexLocalIndex, 1 ) = verticesIterator->get_coordinates().y * m_scale;
    X( vertexLocalIndex, 2 ) = verticesIterator->get_coordinates().z * m_scale * zReverseFactor;
    allNodes.insert( vertexLocalIndex );

    nodeLocalToGlobal[vertexLocalIndex] = vertexGlobalIndex;
    for( int i = 0; i < 3; i++ )
    {
      if( X( vertexLocalIndex, i ) > xMax[i] )
      {
        xMax[i] = X( vertexLocalIndex, i );
      }
      if( X( vertexLocalIndex, i ) < xMin[i] )
      {
        xMin[i] = X( vertexLocalIndex, i );
      }
    }
  }

  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( xMax ) );

  // First loop which iterate on the regions
  array1d< globalIndex > globalIndexRegionOffset( polyhedronPartMap.size() +1 );
  for( auto const & polyhedronPart : polyhedronPartMap )
  {
    auto const regionPtr = polyhedronPart.second;
    string regionName = DecodePAMELALabels::retrieveSurfaceOrRegionName( regionPtr->Label );

    // Iterate on cell types
    for( auto const & subPart : regionPtr->SubParts )
    {
      auto const cellBlockPAMELA = subPart.second;
      PAMELA::ELEMENTS::TYPE const cellBlockType = cellBlockPAMELA->ElementType;
      string const & cellBlockName = ElementToLabel.at( cellBlockType );
      CellBlock * cellBlock = nullptr;
      if( cellBlockName == "HEX" )
      {
        localIndex const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        cellBlock = &cellBlockManager.registerCellBlock( DecodePAMELALabels::makeRegionLabel( regionName, cellBlockName ) );
        cellBlock->setElementType( "C3D8" );
        auto & cellToVertex = cellBlock->getElemToNode();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 8 );

        arrayView1d< globalIndex > const & localToGlobal = cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto const cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][5] =
            cornerList[5]->get_localIndex();
          cellToVertex[cellLocalIndex][6] =
            cornerList[7]->get_localIndex();
          cellToVertex[cellLocalIndex][7] =
            cornerList[6]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "TETRA" )
      {
        localIndex const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        cellBlock = &cellBlockManager.registerCellBlock( DecodePAMELALabels::makeRegionLabel( regionName, cellBlockName ) );
        cellBlock->setElementType( "C3D4" );
        auto & cellToVertex = cellBlock->getElemToNode();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 4 );

        arrayView1d< globalIndex > const & localToGlobal = cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto const cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[3]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "WEDGE" )
      {
        localIndex const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        cellBlock = &cellBlockManager.registerCellBlock( DecodePAMELALabels::makeRegionLabel( regionName, cellBlockName ) );
        cellBlock->setElementType( "C3D6" );
        auto & cellToVertex = cellBlock->getElemToNode();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 6 );

        arrayView1d< globalIndex > const & localToGlobal = cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto const cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][5] =
            cornerList[5]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "PYRAMID" )
      {
        localIndex const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        cellBlock = &cellBlockManager.registerCellBlock( DecodePAMELALabels::makeRegionLabel( regionName, cellBlockName ) );
        cellBlock->setElementType( "C3D5" );
        auto & cellToVertex = cellBlock->getElemToNode();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 5 );

        arrayView1d< globalIndex > const & localToGlobal = cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = (*cellItr)->get_localIndex();
          globalIndex cellGlobalIndex = (*cellItr)->get_globalIndex();
          auto const cornerList = (*cellItr)->get_vertexList();

          cellToVertex[cellLocalIndex][0] =
            cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] =
            cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] =
            cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] =
            cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][4] =
            cornerList[4]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }

      /// Import ppt
      if( cellBlock != nullptr && cellBlock->size() > 0 )
      {
        for( localIndex fieldIndex = 0; fieldIndex < m_fieldNamesInGEOSX.size(); fieldIndex++ )
        {
          auto const meshProperty = regionPtr->FindVariableByName( m_fieldsToImport[fieldIndex] );
          PAMELA::VARIABLE_DIMENSION const dimension = meshProperty->Dimension;
          if( dimension == PAMELA::VARIABLE_DIMENSION::SCALAR )
          {
            real64_array & property = cellBlock->addProperty< real64_array >( m_fieldNamesInGEOSX[fieldIndex] );
            GEOSX_ERROR_IF( property.size() != LvArray::integerConversion< localIndex >( meshProperty->size() ),
                            "Viewer size (" << property.size() << ") mismatch with property size in PAMELA ("
                                            << meshProperty->size() << ") on " <<cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < property.size(); cellIndex++ )
            {
              property[cellIndex] = meshProperty->get_data( cellIndex )[0];
            }
          }
          else if( dimension == PAMELA::VARIABLE_DIMENSION::VECTOR )
          {
            array2d< real64 > & property = cellBlock->addProperty< array2d< real64 > >( m_fieldNamesInGEOSX[fieldIndex] );
            property.resizeDimension< 1 >( 3 );
            GEOSX_ERROR_IF( property.size() != LvArray::integerConversion< localIndex >( meshProperty->size() ),
                            "Viewer size (" << property.size() << ") mismatch with property size in PAMELA ("
                                            << meshProperty->size() << ") on " << cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < cellBlock->size(); cellIndex++ )
            {
              for( int dim = 0; dim < 3; dim++ )
              {
                property( cellIndex, dim ) = meshProperty->get_data( cellIndex )[dim];
              }
            }
          }
          else
          {
            GEOSX_ERROR( "Dimension of " <<  m_fieldNamesInGEOSX[fieldIndex] << " is not supported for import in GEOSX" );
          }
        }
      }
    }
  }

  /// Import surfaces
  auto const polygonPartMap = std::get< 0 >( PAMELA::getPolygonPartMap( m_pamelaMesh.get(), 0 ));
  for( auto const & polygonPart : polygonPartMap )
  {
    auto const surfacePtr = polygonPart.second;

    string surfaceName = DecodePAMELALabels::retrieveSurfaceOrRegionName( surfacePtr->Label );
    SortedArray< localIndex > & curNodeSet = nodeSets[surfaceName];
    for( auto const & subPart : surfacePtr->SubParts )
    {
      auto const cellBlockPAMELA = subPart.second;
      PAMELA::ELEMENTS::TYPE const cellBlockType = cellBlockPAMELA->ElementType;
      string const cellBlockName = ElementToLabel.at( cellBlockType );
      if( cellBlockName == "TRIANGLE"  || cellBlockName == "QUAD" )
      {
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          auto const cornerList = (*cellItr)->get_vertexList();
          for( auto corner :cornerList )
          {
            curNodeSet.insert( corner->get_localIndex() );
          }
        }
      }
    }
  }

  cellBlockManager.buildMaps();
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, PAMELAMeshGenerator, string const &, Group * const )
}
