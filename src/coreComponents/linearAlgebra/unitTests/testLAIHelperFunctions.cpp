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
 * @file testLAIHelperFunctions.cpp
 */

#include "gtest/gtest.h"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "managers/initialization.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/DomainPartition.hpp"
#include "meshUtilities/MeshManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

using namespace geosx;

static real64 const machinePrecision = std::numeric_limits< real64 >::epsilon();
static real64 const tolerance = machinePrecision;  //1e-10;

class LAIHelperFunctionsTest : public ::testing::Test
{
protected:
  /**
   * @brief Provide a very simple xml input, with an internally generated mesh.
   */
  static void
  SetUpTestCase()
  {
    problemManager = new ProblemManager( "Problem", nullptr );

    string const inputStream =
      "<Problem>"
      "  <Mesh>"
      "    <InternalMesh name=\"mesh1\""
      "                  elementTypes=\"{C3D8}\""
      "                  xCoords=\"{0, 1}\""
      "                  yCoords=\"{0, 1}\""
      "                  zCoords=\"{0, 1}\""
      "                  nx=\"{6}\""
      "                  ny=\"{9}\""
      "                  nz=\"{5}\""
      "                  cellBlockNames=\"{block1}\"/>"
      "  </Mesh>"
      "  <ElementRegions>"
      "    <CellElementRegion name=\"region1\" cellBlocks=\"{block1}\" "
      "materialList=\"{dummy_material}\" />"
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

    int mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
    dataRepository::Group * commandLine =
      problemManager->GetGroup< dataRepository::Group >(
        problemManager->groupKeys.commandLine );
    commandLine
      ->registerWrapper< integer >(
        problemManager->viewKeys.xPartitionsOverride.Key() )
      ->setApplyDefaultValue( mpiSize );

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

ProblemManager * LAIHelperFunctionsTest::problemManager =
  nullptr;  //!< the main problemManager.

TEST_F( LAIHelperFunctionsTest, Test_NodalVectorPermutation )
{
  DomainPartition * const domain = problemManager->getDomainPartition();
  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager * const nodeManager = meshLevel->getNodeManager();

  arrayView1d< globalIndex const > const & nodeLocalToGlobal =
    nodeManager->localToGlobalMap();

  DofManager dofManager( "test" );
  dofManager.setMesh( *domain, 0, 0 );

  string_array Region;
  Region.emplace_back( "region1" );

  dofManager.addField( "nodalVariable", DofManager::Location::Node, 3, Region );
  dofManager.addCoupling( "nodalVariable",
                          "nodalVariable",
                          DofManager::Connector::Elem );
  dofManager.reorderByRank();

  localIndex nDof = 3 * nodeManager->size();

  arrayView1d< globalIndex > const & dofNumber =
    nodeManager->getReference< globalIndex_array >(
      dofManager.getKey( "nodalVariable" ) );
  arrayView1d< integer > const & isNodeGhost = nodeManager->ghostRank();

  ParallelVector nodalVariable, expectedPermutedVector;
  nodalVariable.createWithLocalSize( nDof, MPI_COMM_GEOSX );
  nodalVariable.set( 0 );
  expectedPermutedVector.createWithLocalSize( nDof, MPI_COMM_GEOSX );
  expectedPermutedVector.set( 0 );

  nodalVariable.open();
  expectedPermutedVector.open();
  for( localIndex a = 0; a < nodeManager->size(); a++ )
  {
    if( isNodeGhost[a] < 0 )
    {
      for( localIndex d = 0; d < 3; d++ )
      {
        localIndex index = dofNumber( a ) + d;
        real64 Value = nodeLocalToGlobal[a] * 3 + d;
        nodalVariable.add( index, Value );
        index = nodeLocalToGlobal[a] * 3 + d;
        expectedPermutedVector.add( index, Value );
      }
    }
  }
  nodalVariable.close();
  expectedPermutedVector.close();

  ParallelMatrix permutationMatrix;

  LAIHelperFunctions::CreatePermutationMatrix(
    nodeManager,
    nDof,
    nDof,
    3,
    dofManager.getKey( "nodalVariable" ),
    permutationMatrix );

  ParallelVector permutedVector =
    LAIHelperFunctions::PermuteVector( nodalVariable, permutationMatrix );

  permutedVector.axpy( -1, expectedPermutedVector );

  real64 vectorNorm = permutedVector.norm1();

  EXPECT_LT( vectorNorm, tolerance );
}

TEST_F( LAIHelperFunctionsTest, Test_CellCenteredVectorPermutation )
{
  DomainPartition * const domain = problemManager->getDomainPartition();
  MeshLevel * const meshLevel = domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  ;

  DofManager dofManager( "test" );
  dofManager.setMesh( *domain, 0, 0 );

  string_array region;
  region.emplace_back( "region1" );

  dofManager.addField( "cellCentered", DofManager::Location::Elem, region );
  dofManager.addCoupling( "cellCentered",
                          "cellCentered",
                          DofManager::Connector::Face );
  dofManager.reorderByRank();

  integer nDof = dofManager.numGlobalDofs( "cellCentered" );

  ParallelVector cellCenteredVariable, expectedPermutedVector;
  cellCenteredVariable.createWithLocalSize( nDof, MPI_COMM_GEOSX );
  cellCenteredVariable.set( 0 );
  expectedPermutedVector.createWithLocalSize( nDof, MPI_COMM_GEOSX );
  expectedPermutedVector.set( 0 );

  cellCenteredVariable.open();
  expectedPermutedVector.open();
  elemManager->forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase const & elementSubRegion ) {
      localIndex const numElems = elementSubRegion.size();
      arrayView1d< globalIndex const > const & dofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >(
          dofManager.getKey( "cellCentered" ) );
      arrayView1d< integer const > const & isGhost = elementSubRegion.ghostRank();
      arrayView1d< globalIndex const > const & localToGlobal =
        elementSubRegion.localToGlobalMap();

      for( localIndex k = 0; k < numElems; ++k )
      {
        if( dofNumber[k] >= 0 && isGhost[k] < 0 )
        {
          globalIndex index = dofNumber[k];
          real64 Value = localToGlobal[k];
          cellCenteredVariable.add( index, Value );
          index = localToGlobal[k];
          Value = localToGlobal[k];
          expectedPermutedVector.add( index, Value );
        }
      }
    } );

  cellCenteredVariable.close();
  expectedPermutedVector.close();

  ParallelMatrix permutationMatrix;

  LAIHelperFunctions::CreatePermutationMatrix( elemManager,
                                               nDof,
                                               nDof,
                                               1,
                                               dofManager.getKey( "cellCentered" ),
                                               permutationMatrix );

  // permutationMatrix.print(std::cout);
  ParallelVector permutedVector =
    LAIHelperFunctions::PermuteVector( cellCenteredVariable, permutationMatrix );

  permutedVector.axpy( -1, expectedPermutedVector );

  real64 vectorNorm = permutedVector.norm1();

  EXPECT_LT( vectorNorm, tolerance );
}

/**
 * @function main
 * @brief Main function to setup the GEOSX environment, read the xml file and run all cases.
 */
int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
