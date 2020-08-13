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

#include "gtest/gtest.h"

#include "managers/initialization.hpp"

#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/CellElementSubRegion.hpp"

using namespace geosx;

#define STRINGIFY2( X ) #X
#define STRINGIFY( X ) STRINGIFY2( X )

#define NX 10
#define NY 11
#define NZ 12

#define MAX_COORD_X 1.0
#define MAX_COORD_Y 2.0
#define MAX_COORD_Z 3.0

constexpr localIndex numElemsInX = NX;
constexpr localIndex numElemsInY = NY;
constexpr localIndex numElemsInZ = NZ;

constexpr localIndex numNodesInX = NX + 1;
constexpr localIndex numNodesInY = NY + 1;
constexpr localIndex numNodesInZ = NZ + 1;

constexpr double dx = MAX_COORD_X / NX;
constexpr double dy = MAX_COORD_Y / NY;
constexpr double dz = MAX_COORD_Z / NZ;

constexpr localIndex node_dI = numNodesInY * numNodesInZ;
constexpr localIndex node_dJ = numNodesInZ;

constexpr localIndex elem_dI = numElemsInY * numElemsInZ;
constexpr localIndex elem_dJ = numElemsInZ;

class MeshGenerationTest : public ::testing::Test
{
protected:
  void
  SetUp() override
  {
    m_nodeManager = problemManager->getDomainPartition()
                      ->getMeshBody( 0 )
                      ->getMeshLevel( 0 )
                      ->getNodeManager();
    m_faceManager = problemManager->getDomainPartition()
                      ->getMeshBody( 0 )
                      ->getMeshLevel( 0 )
                      ->getFaceManager();
    m_edgeManager = problemManager->getDomainPartition()
                      ->getMeshBody( 0 )
                      ->getMeshLevel( 0 )
                      ->getEdgeManager();

    ElementRegionManager * const elemManager =
      problemManager->getDomainPartition()
        ->getMeshBody( 0 )
        ->getMeshLevel( 0 )
        ->getElemManager();
    GEOSX_ERROR_IF_NE_MSG( elemManager->GetRegions().size(),
                           1,
                           "Only one region should exist." );

    ElementRegionBase * const elemRegion = elemManager->GetRegion( 0 );
    GEOSX_ERROR_IF_NE_MSG( elemRegion->GetSubRegions().size(),
                           1,
                           "Only one subregion should exist." );

    m_subRegion = elemRegion->GetSubRegion< CellElementSubRegion >( 0 );
  }

  NodeManager * m_nodeManager;
  FaceManager * m_faceManager;
  EdgeManager * m_edgeManager;
  CellElementSubRegion * m_subRegion;

  static void
  SetUpTestCase()
  {
    problemManager = new ProblemManager( "Problem", nullptr );

    string const inputStream =
      "<Problem>"
      "  <Mesh>"
      "    <InternalMesh name=\"mesh1\""
      "                  elementTypes=\"{C3D8}\""
      "                  xCoords=\"{0, " STRINGIFY( MAX_COORD_X ) "}\""
                                                                  "                  yCoords=\"{0, " STRINGIFY( MAX_COORD_Y ) "}\""
                                                                                                                              "                  zCoords=\"{0, " STRINGIFY( MAX_COORD_Z ) "}\""
                                                                                                                                                                                          "                  nx=\"{" STRINGIFY( NX ) "}\""
                                                                                                                                                                                                                                     "                  ny=\"{" STRINGIFY( NY ) "}\""
                                                                                                                                                                                                                                                                                "                  nz=\"{" STRINGIFY( NZ ) "}\""
                                                                                                                                                                                                                                                                                                                           "                  cellBlockNames=\"{cb1}\"/>"
                                                                                                                                                                                                                                                                                                                           "  </Mesh>"
                                                                                                                                                                                                                                                                                                                           "  <ElementRegions>"
                                                                                                                                                                                                                                                                                                                           "    <CellElementRegion name=\"region1\" cellBlocks=\"{cb1}\" materialList=\"{dummy_material}\" />"
                                                                                                                                                                                                                                                                                                                           "  </ElementRegions>"
                                                                                                                                                                                                                                                                                                                           "</Problem>";

    xmlWrapper::xmlDocument xmlDocument;
    xmlWrapper::xmlResult xmlResult =
      xmlDocument.load_buffer( inputStream.c_str(), inputStream.size() );
    if( !xmlResult )
    {
      GEOSX_LOG_RANK_0( "XML parsed with errors!" );
      GEOSX_LOG_RANK_0( "Error description: " << xmlResult.description() );
      GEOSX_LOG_RANK_0( "Error offset: " << xmlResult.offset );
    }

    xmlWrapper::xmlNode xmlProblemNode = xmlDocument.child( "Problem" );
    problemManager->InitializePythonInterpreter();
    problemManager->ProcessInputFileRecursive( xmlProblemNode );

    // Open mesh levels
    DomainPartition * domain = problemManager->getDomainPartition();
    MeshManager * meshManager = problemManager->GetGroup< MeshManager >(
      problemManager->groupKeys.meshManager );
    meshManager->GenerateMeshLevels( domain );

    ElementRegionManager * elementManager =
      domain->getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
    xmlWrapper::xmlNode topLevelNode =
      xmlProblemNode.child( elementManager->getName().c_str() );
    elementManager->ProcessInputFileRecursive( topLevelNode );
    elementManager->PostProcessInputRecursive();

    problemManager->ProblemSetup();
  }

  /**
   * @brief Destructor.
   */
  static void
  TearDownTestCase()
  {
    delete problemManager;
    problemManager = nullptr;
  }

  static ProblemManager * problemManager;
};

ProblemManager * MeshGenerationTest::problemManager =
  nullptr;  //!< the main problemManager.

TEST_F( MeshGenerationTest, sizes )
{
  EXPECT_EQ( m_subRegion->numNodesPerElement(), 8 );
  EXPECT_EQ( m_subRegion->numFacesPerElement(), 6 );

  localIndex const numNodes = numNodesInX * numNodesInY * numNodesInZ;
  EXPECT_EQ( numNodes, m_nodeManager->size() );

  localIndex const numElements = numElemsInX * numElemsInY * numElemsInZ;
  EXPECT_EQ( numElements, m_subRegion->size() );

  localIndex const numFaces = numNodesInX * numElemsInY * numElemsInZ +
    numElemsInX * numNodesInY * numElemsInZ +
    numElemsInX * numElemsInY * numNodesInZ;
  EXPECT_EQ( numFaces, m_faceManager->size() );

  localIndex const numEdges = numNodesInX * numNodesInY * numElemsInZ +
    numNodesInX * numElemsInY * numNodesInZ +
    numElemsInX * numNodesInY * numNodesInZ;
  EXPECT_EQ( numEdges, m_edgeManager->size() );
}

TEST_F( MeshGenerationTest, nodePositions )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X =
    m_nodeManager->referencePosition();

  localIndex nodeID = 0;
  for( localIndex i = 0; i < numNodesInX; ++i )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex k = 0; k < numNodesInZ; ++k )
      {
        EXPECT_DOUBLE_EQ( X( nodeID, 0 ), i * dx );
        EXPECT_DOUBLE_EQ( X( nodeID, 1 ), j * dy );
        EXPECT_DOUBLE_EQ( X( nodeID, 2 ), k * dz );
        ++nodeID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, elementCentersAndVolumes )
{
  arrayView2d< real64 const > const & centers = m_subRegion->getElementCenter();
  arrayView1d< real64 const > const & volumes = m_subRegion->getElementVolume();

  constexpr double VOLUME = dx * dy * dz;

  localIndex elemID = 0;
  for( localIndex i = 0; i < numElemsInX; ++i )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex k = 0; k < numElemsInZ; ++k )
      {
        EXPECT_DOUBLE_EQ( centers[elemID][0], i * dx + dx / 2.0 );
        EXPECT_DOUBLE_EQ( centers[elemID][1], j * dy + dy / 2.0 );
        EXPECT_DOUBLE_EQ( centers[elemID][2], k * dz + dz / 2.0 );
        EXPECT_NEAR( volumes[elemID], VOLUME, 1e-13 );
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, elemToNodeMap )
{
  arrayView2d< localIndex const, cells::NODE_MAP_USD > const & nodeMap =
    m_subRegion->nodeList();
  GEOSX_ERROR_IF_NE( nodeMap.size( 1 ), 8 );

  localIndex elemID = 0;
  for( localIndex i = 0; i < numElemsInX; ++i )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex k = 0; k < numElemsInZ; ++k )
      {
        localIndex const firstNodeID = i * node_dI + j * node_dJ + k;

        EXPECT_EQ( firstNodeID, nodeMap( elemID, 0 ) );
        EXPECT_EQ( firstNodeID + node_dI, nodeMap( elemID, 1 ) );
        EXPECT_EQ( firstNodeID + node_dI + node_dJ, nodeMap( elemID, 3 ) );
        EXPECT_EQ( firstNodeID + node_dJ, nodeMap( elemID, 2 ) );

        EXPECT_EQ( firstNodeID + 1, nodeMap( elemID, 4 ) );
        EXPECT_EQ( firstNodeID + 1 + node_dI, nodeMap( elemID, 5 ) );
        EXPECT_EQ( firstNodeID + 1 + node_dI + node_dJ, nodeMap( elemID, 7 ) );
        EXPECT_EQ( firstNodeID + 1 + node_dJ, nodeMap( elemID, 6 ) );
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, nodeToElemMap )
{
  ArrayOfArraysView< localIndex const > const & nodeToElemMap =
    m_nodeManager->elementList().toViewConst();

  localIndex nodeID = 0;
  for( localIndex i = 0; i < numNodesInX; ++i )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex k = 0; k < numNodesInZ; ++k )
      {
        localIndex const elemID = i * elem_dI + j * elem_dJ + k;

        std::vector< localIndex > expectedElems;
        if( k < numElemsInZ )
        {
          if( i < numElemsInX && j < numElemsInY )
            expectedElems.push_back( elemID );
          if( i > 0 && j < numElemsInY )
            expectedElems.push_back( elemID - elem_dI );
          if( i > 0 && j > 0 )
            expectedElems.push_back( elemID - elem_dI - elem_dJ );
          if( i < numElemsInX && j > 0 )
            expectedElems.push_back( elemID - elem_dJ );
        }

        if( k > 0 )
        {
          if( i < numElemsInX && j < numElemsInY )
            expectedElems.push_back( elemID - 1 );
          if( i > 0 && j < numElemsInY )
            expectedElems.push_back( elemID - elem_dI - 1 );
          if( i > 0 && j > 0 )
            expectedElems.push_back( elemID - elem_dI - elem_dJ - 1 );
          if( i < numElemsInX && j > 0 )
            expectedElems.push_back( elemID - elem_dJ - 1 );
        }

        localIndex const numElems = expectedElems.size();
        ASSERT_EQ( numElems, nodeToElemMap.sizeOfArray( nodeID ) );

        localIndex const * const nodeElems = nodeToElemMap[nodeID];
        std::vector< localIndex > elems( nodeElems, nodeElems + numElems );

        std::sort( elems.begin(), elems.end() );
        std::sort( expectedElems.begin(), expectedElems.end() );

        for( localIndex a = 0; a < numElems; ++a )
        {
          EXPECT_EQ( elems[a], expectedElems[a] );
        }

        ++nodeID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, faceNodeMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap =
    m_subRegion->faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    m_faceManager->nodeList().toViewConst();
  arrayView2d< localIndex const > const & faceToElementMap =
    m_faceManager->elementList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap =
    m_nodeManager->faceList().toViewConst();

  GEOSX_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  array1d< localIndex > faceNodesFromElem( 4 );
  array1d< localIndex > faceNodesFromFace( 4 );

  localIndex elemID = 0;
  for( localIndex i = 0; i < numElemsInX; ++i )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex k = 0; k < numElemsInZ; ++k )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          m_subRegion->GetFaceNodes( elemID, f, faceNodesFromElem );
          ASSERT_EQ( faceNodesFromElem.size(), 4 );

          localIndex const faceID = elementToFaceMap( elemID, f );

          ASSERT_EQ( faceToNodeMap.sizeOfArray( faceID ), 4 );
          for( localIndex a = 0; a < 4; ++a )
          {
            faceNodesFromFace[a] = faceToNodeMap( faceID, a );
          }

          if( elemID != faceToElementMap( faceID, 0 ) )
          {
            EXPECT_EQ( elemID, faceToElementMap( faceID, 1 ) );
          }

          std::sort( faceNodesFromElem.begin(), faceNodesFromElem.end() );
          std::sort( faceNodesFromFace.begin(), faceNodesFromFace.end() );

          EXPECT_EQ( faceNodesFromElem[0], faceNodesFromFace[0] );
          EXPECT_EQ( faceNodesFromElem[1], faceNodesFromFace[1] );
          EXPECT_EQ( faceNodesFromElem[2], faceNodesFromFace[2] );
          EXPECT_EQ( faceNodesFromElem[3], faceNodesFromFace[3] );

          for( localIndex a = 0; a < 4; ++a )
          {
            EXPECT_TRUE( nodeToFaceMap.contains( faceNodesFromElem[a], faceID ) );
          }
        }
        ++elemID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, faceElementMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap = m_subRegion->faceList();
  arrayView2d< localIndex const > const & faceToElementMap =
    m_faceManager->elementList();

  GEOSX_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  localIndex const elemIDOffset[6] = { -elem_dJ, -1, -elem_dI, elem_dI, elem_dJ, 1 };

  localIndex elemID = 0;
  for( localIndex i = 0; i < numElemsInX; ++i )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex k = 0; k < numElemsInZ; ++k )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          localIndex const faceID = elementToFaceMap( elemID, f );

          if( elemID == faceToElementMap( faceID, 0 ) )
          {
            bool external = false;
            if( f == 0 && j == 0 ) external = true;
            if( f == 1 && k == 0 ) external = true;
            if( f == 2 && i == 0 ) external = true;
            if( f == 3 && i == numElemsInX - 1 ) external = true;
            if( f == 4 && j == numElemsInY - 1 ) external = true;
            if( f == 5 && k == numElemsInZ - 1 ) external = true;

            if( external )
            {
              EXPECT_EQ( -1, faceToElementMap( faceID, 1 ) );
            }
            else
            {
              EXPECT_EQ( elemID + elemIDOffset[f], faceToElementMap( faceID, 1 ) );
            }
          }
          else
          {
            EXPECT_EQ( elemID, faceToElementMap( faceID, 1 ) );
            EXPECT_EQ( elemID + elemIDOffset[f], faceToElementMap( faceID, 0 ) );
          }
        }
        ++elemID;
      }
    }
  }
}

bool
walkEdgesToFindNeighbor( localIndex const node0,
                         localIndex const node1,
                         arraySlice1d< localIndex const > const & nodeEdges,
                         arrayView2d< localIndex const > const & edgeToNodeMap )
{
  for( localIndex const edgeID : nodeEdges )
  {
    if( edgeToNodeMap[edgeID][0] == node0 )
    {
      if( edgeToNodeMap[edgeID][1] == node1 ) return true;
    }
    else
    {
      EXPECT_EQ( edgeToNodeMap[edgeID][1], node0 );
      if( edgeToNodeMap[edgeID][0] == node1 ) return true;
    }
  }

  return false;
}

TEST_F( MeshGenerationTest, edgeNodeMaps )
{
  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap =
    m_nodeManager->edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = m_edgeManager->nodeList();

  GEOSX_ERROR_IF_NE( edgeToNodeMap.size( 1 ), 2 );

  localIndex nodeID = 0;
  for( localIndex i = 0; i < numNodesInX; ++i )
  {
    for( localIndex j = 0; j < numNodesInY; ++j )
    {
      for( localIndex k = 0; k < numNodesInZ; ++k )
      {
        localIndex numEdges = 0;
        if( i != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID - node_dI,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }
        if( i != numNodesInX - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID + node_dI,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }
        if( j != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID - node_dJ,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }
        if( j != numNodesInY - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID + node_dJ,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }
        if( k != 0 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID - 1,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }
        if( k != numNodesInZ - 1 )
        {
          EXPECT_TRUE( walkEdgesToFindNeighbor( nodeID,
                                                nodeID + 1,
                                                nodeToEdgeMap[nodeID],
                                                edgeToNodeMap ) );
          ++numEdges;
        }

        EXPECT_EQ( numEdges, nodeToEdgeMap.sizeOfSet( nodeID ) );
        ++nodeID;
      }
    }
  }
}

TEST_F( MeshGenerationTest, edgeFaceMaps )
{
  arrayView2d< localIndex const > const & elementToFaceMap = m_subRegion->faceList();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    m_faceManager->nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap =
    m_faceManager->edgeList().toViewConst();
  arrayView2d< localIndex const > const & edgeToNodeMap = m_edgeManager->nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap =
    m_edgeManager->faceList().toViewConst();

  GEOSX_ERROR_IF_NE( elementToFaceMap.size( 1 ), 6 );

  localIndex elemID = 0;
  for( localIndex i = 0; i < numElemsInX; ++i )
  {
    for( localIndex j = 0; j < numElemsInY; ++j )
    {
      for( localIndex k = 0; k < numElemsInZ; ++k )
      {
        for( localIndex f = 0; f < 6; ++f )
        {
          localIndex const faceID = elementToFaceMap( elemID, f );

          ASSERT_EQ( faceToNodeMap.sizeOfArray( faceID ), 4 );
          ASSERT_EQ( faceToEdgeMap.sizeOfArray( faceID ), 4 );

          for( localIndex a = 0; a < 4; ++a )
          {
            localIndex node0 = faceToNodeMap( faceID, a );
            localIndex node1 = faceToNodeMap( faceID, ( a + 1 ) % 4 );
            if( node0 > node1 ) std::swap( node0, node1 );

            localIndex const edgeID = faceToEdgeMap( faceID, a );
            EXPECT_EQ( edgeToNodeMap( edgeID, 0 ), node0 );
            EXPECT_EQ( edgeToNodeMap( edgeID, 1 ), node1 );

            bool foundFace = false;
            for( localIndex const id : edgeToFaceMap[edgeID] )
            {
              if( id == faceID ) foundFace = true;
            }

            EXPECT_TRUE( foundFace );
          }
        }
        ++elemID;
      }
    }
  }
}

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
