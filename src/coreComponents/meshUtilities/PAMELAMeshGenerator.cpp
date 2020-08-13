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
 * @file PAMELAMeshGenerator.cpp
 */

#include "PAMELAMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include <math.h>

#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

/// TODO when we are going to port to c++17, we can remove that
string const PAMELAMeshGenerator::DecodePAMELALabels::m_separator = "_";

PAMELAMeshGenerator::PAMELAMeshGenerator( string const & name, Group * const parent ) :
  MeshGeneratorBase( name, parent )
{
  registerWrapper( viewKeyStruct::filePathString, &m_filePath )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setRestartFlags( RestartFlags::NO_WRITE )
    ->setDescription( "path to the mesh file" );
  registerWrapper( viewKeyStruct::fieldsToImportString, &m_fieldsToImport )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Fields to be imported from the external mesh file" );
  registerWrapper( viewKeyStruct::fieldNamesInGEOSXString, &m_fieldNamesInGEOSX )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Name of the fields within GEOSX" );
  registerWrapper( viewKeyStruct::scaleString, &m_scale )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDefaultValue( 1. )
    ->setDescription( "Scale the coordinates of the vertices" );
  registerWrapper( viewKeyStruct::reverseZString, &m_isZReverse )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDefaultValue( 0 )
    ->setDescription(
      "0 : Z coordinate is upward, 1 : Z coordinate is downward" );
}

PAMELAMeshGenerator::~PAMELAMeshGenerator()
{}

void
PAMELAMeshGenerator::GenerateElementRegions(
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void
PAMELAMeshGenerator::PostProcessInput()
{
  m_pamelaMesh =
    std::unique_ptr< PAMELA::Mesh >( PAMELA::MeshFactory::makeMesh( m_filePath ) );
  m_pamelaMesh->CreateFacesFromCells();
  m_pamelaMesh->PerformPolyhedronPartitioning( PAMELA::ELEMENTS::FAMILY::POLYGON,
                                               PAMELA::ELEMENTS::FAMILY::POLYGON );
  m_pamelaMesh->CreateLineGroupWithAdjacency(
    "TopologicalC2C",
    m_pamelaMesh->getAdjacencySet()->get_TopologicalAdjacency(
      PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
      PAMELA::ELEMENTS::FAMILY::POLYHEDRON,
      PAMELA::ELEMENTS::FAMILY::POLYGON ) );
}

void
PAMELAMeshGenerator::RemapMesh(
  dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) )
{
  return;
}

Group *
PAMELAMeshGenerator::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                  string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void
PAMELAMeshGenerator::GenerateMesh( DomainPartition * const domain )
{
  GEOSX_LOG_RANK_0( "Writing into the GEOSX mesh data structure" );
  domain->getMetisNeighborList() = m_pamelaMesh->getNeighborList();
  Group * const meshBodies = domain->GetGroup( std::string( "MeshBodies" ) );
  MeshBody * const meshBody = meshBodies->RegisterGroup< MeshBody >( this->getName() );

  //TODO for the moment we only consider on mesh level "Level0"
  MeshLevel * const meshLevel0 =
    meshBody->RegisterGroup< MeshLevel >( std::string( "Level0" ) );
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  CellBlockManager * cellBlockManager =
    domain->GetGroup< CellBlockManager >( keys::cellManager );

  // Use the PartMap of PAMELA to get the mesh
  auto const polyhedronPartMap =
    std::get< 0 >( PAMELA::getPolyhedronPartMap( m_pamelaMesh.get(), 0 ) );

  // Vertices are written first
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();
  nodeManager->resize( m_pamelaMesh->get_PointCollection()->size_all() );

  arrayView1d< globalIndex > const & nodeLocalToGlobal =
    nodeManager->localToGlobalMap();

  Group & nodeSets = nodeManager->sets();
  SortedArray< localIndex > & allNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();

  R1Tensor xMax( std::numeric_limits< real64 >::min(),
                 std::numeric_limits< real64 >::min(),
                 std::numeric_limits< real64 >::min() );

  R1Tensor xMin( std::numeric_limits< real64 >::max(),
                 std::numeric_limits< real64 >::max(),
                 std::numeric_limits< real64 >::max() );

  double zReverseFactor = 1.;
  if( m_isZReverse )
  {
    zReverseFactor = -1.;
  }
  for( auto const & verticesIterator : *m_pamelaMesh->get_PointCollection() )
  {
    localIndex const vertexLocalIndex = verticesIterator->get_localIndex();
    globalIndex const vertexGlobalIndex = verticesIterator->get_globalIndex();
    X( vertexLocalIndex, 0 ) = verticesIterator->get_coordinates().x * m_scale;
    X( vertexLocalIndex, 1 ) = verticesIterator->get_coordinates().y * m_scale;
    X( vertexLocalIndex, 2 ) =
      verticesIterator->get_coordinates().z * m_scale * zReverseFactor;
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
  xMax -= xMin;
  meshBody->setGlobalLengthScale( std::fabs( xMax.L2_Norm() ) );

  // First loop which iterate on the regions
  array1d< globalIndex > globalIndexRegionOffset( polyhedronPartMap.size() + 1 );
  for( auto const & polyhedronPart : polyhedronPartMap )
  {
    auto const regionPtr = polyhedronPart.second;
    string regionName =
      DecodePAMELALabels::RetrieveSurfaceOrRegionName( regionPtr->Label );

    // Iterate on cell types
    for( auto const & subPart : regionPtr->SubParts )
    {
      auto const cellBlockPAMELA = subPart.second;
      auto const cellBlockType = cellBlockPAMELA->ElementType;
      auto const cellBlockName = ElementToLabel.at( cellBlockType );
      CellBlock * cellBlock = nullptr;
      if( cellBlockName == "HEX" )
      {
        auto const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )
            ->RegisterGroup< CellBlock >(
              DecodePAMELALabels::MakeRegionLabel( regionName, cellBlockName ) );
        cellBlock->SetElementType( "C3D8" );
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 8 );

        arrayView1d< globalIndex > const & localToGlobal =
          cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = ( *cellItr )->get_localIndex();
          globalIndex cellGlobalIndex = ( *cellItr )->get_globalIndex();
          auto const cornerList = ( *cellItr )->get_vertexList();

          cellToVertex[cellLocalIndex][0] = cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] = cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] = cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][3] = cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][4] = cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][5] = cornerList[5]->get_localIndex();
          cellToVertex[cellLocalIndex][6] = cornerList[7]->get_localIndex();
          cellToVertex[cellLocalIndex][7] = cornerList[6]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "TETRA" )
      {
        auto const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )
            ->RegisterGroup< CellBlock >(
              DecodePAMELALabels::MakeRegionLabel( regionName, cellBlockName ) );
        cellBlock->SetElementType( "C3D4" );
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 4 );

        arrayView1d< globalIndex > const & localToGlobal =
          cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = ( *cellItr )->get_localIndex();
          globalIndex cellGlobalIndex = ( *cellItr )->get_globalIndex();
          auto const cornerList = ( *cellItr )->get_vertexList();

          cellToVertex[cellLocalIndex][0] = cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] = cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] = cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] = cornerList[3]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "WEDGE" )
      {
        auto const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )
            ->RegisterGroup< CellBlock >(
              DecodePAMELALabels::MakeRegionLabel( regionName, cellBlockName ) );
        cellBlock->SetElementType( "C3D6" );
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 6 );

        arrayView1d< globalIndex > const & localToGlobal =
          cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = ( *cellItr )->get_localIndex();
          globalIndex cellGlobalIndex = ( *cellItr )->get_globalIndex();
          auto const cornerList = ( *cellItr )->get_vertexList();

          cellToVertex[cellLocalIndex][0] = cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] = cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] = cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] = cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][4] = cornerList[4]->get_localIndex();
          cellToVertex[cellLocalIndex][5] = cornerList[5]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }
      else if( cellBlockName == "PYRAMID" )
      {
        auto const nbCells = cellBlockPAMELA->SubCollection.size_owned();
        if( nbCells == 0 ) continue;
        cellBlock =
          cellBlockManager->GetGroup( keys::cellBlocks )
            ->RegisterGroup< CellBlock >(
              DecodePAMELALabels::MakeRegionLabel( regionName, cellBlockName ) );
        cellBlock->SetElementType( "C3D5" );
        auto & cellToVertex = cellBlock->nodeList();
        cellBlock->resize( nbCells );
        cellToVertex.resize( nbCells, 5 );

        arrayView1d< globalIndex > const & localToGlobal =
          cellBlock->localToGlobalMap();

        // Iterate on cells
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          localIndex cellLocalIndex = ( *cellItr )->get_localIndex();
          globalIndex cellGlobalIndex = ( *cellItr )->get_globalIndex();
          auto const cornerList = ( *cellItr )->get_vertexList();

          cellToVertex[cellLocalIndex][0] = cornerList[0]->get_localIndex();
          cellToVertex[cellLocalIndex][1] = cornerList[1]->get_localIndex();
          cellToVertex[cellLocalIndex][2] = cornerList[2]->get_localIndex();
          cellToVertex[cellLocalIndex][3] = cornerList[3]->get_localIndex();
          cellToVertex[cellLocalIndex][4] = cornerList[4]->get_localIndex();

          localToGlobal[cellLocalIndex] = cellGlobalIndex;
        }
      }

      /// Import ppt
      if( cellBlock != nullptr )
      {
        for( localIndex fieldIndex = 0; fieldIndex < m_fieldNamesInGEOSX.size();
             fieldIndex++ )
        {
          auto const meshProperty =
            regionPtr->FindVariableByName( m_fieldsToImport[fieldIndex] );
          auto const dimension = meshProperty->Dimension;
          if( dimension == PAMELA::VARIABLE_DIMENSION::SCALAR )
          {
            real64_array & property = cellBlock->AddProperty< real64_array >(
              m_fieldNamesInGEOSX[fieldIndex] );
            GEOSX_ERROR_IF(
              property.size() !=
                LvArray::integerConversion< localIndex >( meshProperty->size() ),
              "Viewer size ("
                << property.size() << ") mismatch with property size in PAMELA ("
                << meshProperty->size() << ") on " << cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < property.size(); cellIndex++ )
            {
              property[cellIndex] = meshProperty->get_data( cellIndex )[0];
            }
          }
          else if( dimension == PAMELA::VARIABLE_DIMENSION::VECTOR )
          {
            array1d< R1Tensor > & property =
              cellBlock->AddProperty< array1d< R1Tensor > >(
                m_fieldNamesInGEOSX[fieldIndex] );
            GEOSX_ERROR_IF(
              property.size() * 3 !=
                LvArray::integerConversion< localIndex >( meshProperty->size() ),
              "Viewer size (" << property.size() * 3
                              << ") mismatch with property size in PAMELA ("
                              << meshProperty->size() << ") on "
                              << cellBlock->getName() );
            for( int cellIndex = 0; cellIndex < cellBlock->size(); cellIndex++ )
            {
              for( int dim = 0; dim < 3; dim++ )
              {
                property[cellIndex][dim] = meshProperty->get_data( cellIndex )[dim];
              }
            }
          }
          else
          {
            GEOSX_ERROR( "Dimension of "
                         << m_fieldNamesInGEOSX[fieldIndex]
                         << " is not supported for import in GEOSX" );
          }
        }
      }
    }
  }

  /// Import surfaces
  auto const polygonPartMap =
    std::get< 0 >( PAMELA::getPolygonPartMap( m_pamelaMesh.get(), 0 ) );
  for( auto const & polygonPart : polygonPartMap )
  {
    auto const surfacePtr = polygonPart.second;

    string surfaceName =
      DecodePAMELALabels::RetrieveSurfaceOrRegionName( surfacePtr->Label );
    SortedArray< localIndex > & curNodeSet =
      nodeSets
        .registerWrapper< SortedArray< localIndex > >( std::string( surfaceName ) )
        ->reference();
    for( auto const & subPart : surfacePtr->SubParts )
    {
      auto const cellBlockPAMELA = subPart.second;
      auto const cellBlockType = cellBlockPAMELA->ElementType;
      auto const cellBlockName = ElementToLabel.at( cellBlockType );
      if( cellBlockName == "TRIANGLE" || cellBlockName == "QUAD" )
      {
        for( auto cellItr = cellBlockPAMELA->SubCollection.begin_owned();
             cellItr != cellBlockPAMELA->SubCollection.end_owned();
             cellItr++ )
        {
          auto const cornerList = ( *cellItr )->get_vertexList();
          for( auto corner : cornerList )
          {
            curNodeSet.insert( corner->get_localIndex() );
          }
        }
      }
    }
  }
}

void
PAMELAMeshGenerator::GetElemToNodesRelationInBox(
  const std::string & GEOSX_UNUSED_PARAM( elementType ),
  const int GEOSX_UNUSED_PARAM( index )[],
  const int & GEOSX_UNUSED_PARAM( iEle ),
  int GEOSX_UNUSED_PARAM( nodeIDInBox )[],
  const int GEOSX_UNUSED_PARAM( node_size ) )
{}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase,
                        PAMELAMeshGenerator,
                        std::string const &,
                        Group * const )
}  // namespace geosx
