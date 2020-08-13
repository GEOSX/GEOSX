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
 * @file VTMMeshGenerator.cpp
 */

#include "VTMMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <math.h>

#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
//#include "SimpleGeometricObjects.hpp"

#include "mesh/MeshBody.hpp"

namespace geosx
{
using namespace dataRepository;

VTMMeshGenerator::VTMMeshGenerator( string const & name, Group * const parent ) :
  MeshGeneratorBase( name, parent )
{
  /*
     for( int i=0 ; i<3 ; ++i )
     {
     m_wExtensionMin[i] = 0;
     m_wExtensionMax[i] = 0;
     m_nExtensionLayersMin[i] = 0;
     m_nExtensionLayersMax[i] = 0;
     m_commonRatioMin[i] = 1.5;
     m_commonRatioMax[i] = 1.5;
     }
   */

  registerWrapper< string >( keys::filePath )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "path to the vtm file" );
}

VTMMeshGenerator::~VTMMeshGenerator()
{
  // TODO Auto-generated destructor stub
}

//}
/**
 * @param domain
 */
void
VTMMeshGenerator::GenerateElementRegions(
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  //  lvector numElements;
  //
  //  for( array1d<string>::size_type r=0 ; r<m_regionNames.size() ; ++r )
  //  {
  //    numElements.push_back( 0 );
  //  }
  //
  //  domain.m_feElementManager->resize( numElements, m_regionNames,
  // m_elementType );
}

void
VTMMeshGenerator::PostProcessInput()
{
  m_fileName = this->getReference< string >( keys::filePath );
  m_vtmFile.Load( m_fileName, true, false );
}

void
VTMMeshGenerator::RemapMesh(
  dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) )
{}

Group *
VTMMeshGenerator::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                               string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

void
VTMMeshGenerator::GenerateMesh( DomainPartition * const domain )
{
  /// Basic mesh registration
  Group * const meshBodies = domain->GetGroup( std::string( "MeshBodies" ) );
  MeshBody * const meshBody = meshBodies->RegisterGroup< MeshBody >( this->getName() );
  MeshLevel * const meshLevel0 =
    meshBody->RegisterGroup< MeshLevel >( std::string( "Level0" ) );
  NodeManager * nodeManager = meshLevel0->getNodeManager();
  CellBlockManager * elementManager =
    domain->GetGroup< CellBlockManager >( keys::cellManager );

  /// Region registrations
  for( localIndex rankBlockIndex = 0; rankBlockIndex < m_vtmFile.NumRankBlocks();
       rankBlockIndex++ )
  {
    const auto & rankBlock = m_vtmFile.GetRankBlock( rankBlockIndex );
    for( localIndex meshBlockIndex = 0; meshBlockIndex < rankBlock.NumMeshBlocks();
         meshBlockIndex++ )
    {
      const auto & meshBlock = rankBlock.GetMeshBlock( meshBlockIndex );
      if( meshBlock.IsARegionBlock() )
      {
        const auto & mesh = meshBlock.mesh();
        /// Write nodes
        nodeManager->resize( mesh.NumVertices() );
        arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X =
          nodeManager->referencePosition();
        for( localIndex a = 0; a < mesh.NumVertices(); ++a )
        {
          for( int i = 0; i < 3; ++i )
          {
            X( a, i ) = mesh.Vertex( a )[i];
          }
        }
        /// Cell blocks registrations
        if( mesh.NumHex() > 0 )
        {
          CellBlock * cellBlock = elementManager->GetGroup( keys::cellBlocks )
                                    ->RegisterGroup< CellBlock >( "HEX" );
          cellBlock->SetElementType( "C3D8" );
          auto & cellToVertex = cellBlock->nodeList();
          cellBlock->resize( mesh.NumCells() );
          cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell( 0 ) );

          for( localIndex k = 0; k < mesh.NumCells(); ++k )
          {
            cellToVertex[k][0] = mesh.CellVertexIndex( k, 0 );
            cellToVertex[k][1] = mesh.CellVertexIndex( k, 1 );
            cellToVertex[k][2] = mesh.CellVertexIndex( k, 3 );
            cellToVertex[k][3] = mesh.CellVertexIndex( k, 2 );
            cellToVertex[k][4] = mesh.CellVertexIndex( k, 4 );
            cellToVertex[k][5] = mesh.CellVertexIndex( k, 5 );
            cellToVertex[k][6] = mesh.CellVertexIndex( k, 7 );
            cellToVertex[k][7] = mesh.CellVertexIndex( k, 6 );
          }
        }
        if( mesh.NumTetra() > 0 )
        {
          CellBlock * cellBlock = elementManager->GetGroup( keys::cellBlocks )
                                    ->RegisterGroup< CellBlock >( "TETRA" );
          cellBlock->SetElementType( "C3D4" );
          auto & cellToVertex = cellBlock->nodeList();
          cellBlock->resize( mesh.NumCells() );
          cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell( 0 ) );

          for( localIndex k = 0; k < mesh.NumCells(); ++k )
          {
            cellToVertex[k][0] = mesh.CellVertexIndex( k, 0 );
            cellToVertex[k][1] = mesh.CellVertexIndex( k, 1 );
            cellToVertex[k][2] = mesh.CellVertexIndex( k, 2 );
            cellToVertex[k][3] = mesh.CellVertexIndex( k, 3 );
          }
        }
        if( mesh.NumPrism() > 0 )
        {
          CellBlock * cellBlock = elementManager->GetGroup( keys::cellBlocks )
                                    ->RegisterGroup< CellBlock >( "WEDGE" );
          cellBlock->SetElementType( "C3D6" );
          auto & cellToVertex = cellBlock->nodeList();
          cellBlock->resize( mesh.NumCells() );
          cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell( 0 ) );

          for( localIndex k = 0; k < mesh.NumCells(); ++k )
          {
            cellToVertex[k][0] = mesh.CellVertexIndex( k, 0 );
            cellToVertex[k][1] = mesh.CellVertexIndex( k, 1 );
            cellToVertex[k][2] = mesh.CellVertexIndex( k, 2 );
            cellToVertex[k][3] = mesh.CellVertexIndex( k, 3 );
            cellToVertex[k][4] = mesh.CellVertexIndex( k, 4 );
            cellToVertex[k][5] = mesh.CellVertexIndex( k, 5 );
          }
        }
        if( mesh.NumPyr() > 0 )
        {
          CellBlock * cellBlock = elementManager->GetGroup( keys::cellBlocks )
                                    ->RegisterGroup< CellBlock >( "PYR" );
          cellBlock->SetElementType( "C3D5" );
          auto & cellToVertex = cellBlock->nodeList();
          cellBlock->resize( mesh.NumCells() );
          cellToVertex.resize( mesh.NumCells(), mesh.NumVerticesInCell( 0 ) );

          for( localIndex k = 0; k < mesh.NumCells(); ++k )
          {
            cellToVertex[k][0] = mesh.CellVertexIndex( k, 0 );
            cellToVertex[k][1] = mesh.CellVertexIndex( k, 1 );
            cellToVertex[k][2] = mesh.CellVertexIndex( k, 2 );
            cellToVertex[k][3] = mesh.CellVertexIndex( k, 3 );
            cellToVertex[k][4] = mesh.CellVertexIndex( k, 4 );
          }
        }
      }
    }
  }
}

void
VTMMeshGenerator::GetElemToNodesRelationInBox(
  const std::string & GEOSX_UNUSED_PARAM( elementType ),
  const int GEOSX_UNUSED_PARAM( index )[],
  const int & GEOSX_UNUSED_PARAM( iEle ),
  int GEOSX_UNUSED_PARAM( nodeIDInBox )[],
  const int GEOSX_UNUSED_PARAM( node_size ) )

{}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase,
                        VTMMeshGenerator,
                        std::string const &,
                        Group * const )
}  // namespace geosx
