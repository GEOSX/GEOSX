/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
#include "gtest/gtest.h"

#include "fileIO/xmlWrapper.hpp"

#include "tests/meshFileNames.hpp"

#include "meshUtilities/MeshManager.hpp"

#include "SetSignalHandling.hpp"

#include "stackTrace.hpp"

using namespace geosx;
using namespace geosx::dataRepository;

void TestMeshImport( string const& inputStringMesh,
                     string const& inputStringRegion,
                     string const& propertyToTest )
{
  MeshManager meshManager("mesh", nullptr);

  // Load the mesh
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child("Mesh");
  meshManager.ProcessInputFileRecursive( xmlMeshNode );
  meshManager.PostProcessInputRecursive();

  // Create the domain and generate the Mesh
  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", nullptr ) );
  meshManager.GenerateMeshes( domain.get() );

  ManagedGroup * const meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->GetGroup<MeshBody>(0);
  MeshLevel * const meshLevel = meshBody->GetGroup<MeshLevel>(0);
  NodeManager * const nodeManager = meshLevel->getNodeManager();
  ElementRegionManager * const elemManager = meshLevel->getElemManager();

  // Create the ElementRegions
  xmlDocument.load_buffer( inputStringRegion.c_str(), inputStringRegion.size() );

  xmlWrapper::xmlNode xmlRegionNode = xmlDocument.child("ElementRegions");
  elemManager->ProcessInputFileRecursive( xmlRegionNode );
  elemManager->PostProcessInputRecursive();

  ManagedGroup const * const cellBlockManager = domain->GetGroup(keys::cellManager);

  // This method will call the CopyElementSubRegionFromCellBlocks that will trigger the property transfer.
  elemManager->GenerateMesh( cellBlockManager );


  // Check if the computed center match with the imported center
  if( !propertyToTest.empty() )
  {
    auto centerProperty =  elemManager->ConstructViewAccessor<array1d<R1Tensor>, arrayView1d<R1Tensor> >( propertyToTest );
    elemManager->forElementRegions( [&]( ElementRegion * const elemRegion)->void
    {
      localIndex er = elemRegion->getIndexInParent();
      elemRegion->forElementSubRegionsIndex( [&]( localIndex const esr, auto * const elemSubRegion )
      {
        for( localIndex ei = 0; ei < elemSubRegion->size(); ei++ )
        {
          R1Tensor center = elemSubRegion->calculateElementCenter( ei, *nodeManager );
          R1Tensor centerFromProperty( centerProperty[er][esr][ei][0], 
                                       centerProperty[er][esr][ei][1],
                                       centerProperty[er][esr][ei][2] );
          center -= centerFromProperty;
          GEOS_ERROR_IF( center.L2_Norm() > meshBody->getGlobalLengthScale() * 1e-8, "Property import of centers if wrong");
        }
      });
    });
  }
}

TEST( PAMELAImport, testGMSH )
{
  MeshManager meshManager("mesh", nullptr);

  std::stringstream inputStreamMesh;
  inputStreamMesh <<
  "<?xml version=\"1.0\" ?>" <<
  "  <Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
  "  <PAMELAMeshGenerator name=\"ToyModel\" " <<
  "  fieldsToImport=\"{barycenter}\""<<
  "  fieldNamesInGEOSX=\"{barycenter}\""<<
  "  file=\"" <<gmshFilePath.c_str()<< "\"/>"<<
  "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();

  std::stringstream inputStreamRegion;
  inputStreamRegion <<
  "<?xml version=\"1.0\" ?>" <<
  "  <ElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
  "  <ElementRegion name=\"0\" cellBlocks=\"{0_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"1\" cellBlocks=\"{1_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"2\" cellBlocks=\"{2_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"3\" cellBlocks=\"{3_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"4\" cellBlocks=\"{4_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"5\" cellBlocks=\"{5_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"6\" cellBlocks=\"{6_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "  <ElementRegion name=\"7\" cellBlocks=\"{7_TETRA}\" materialList=\"{water, rock}\"/>" <<
  "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  TestMeshImport( inputStringMesh, inputStringRegion, "barycenter" );
}

TEST( PAMELAImport, testECLIPSE )
{
  MeshManager meshManager("mesh", nullptr);

  std::stringstream inputStreamMesh;
  inputStreamMesh <<
  "<?xml version=\"1.0\" ?>" <<
  "  <Mesh xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
  "  <PAMELAMeshGenerator name=\"ToyModel\" " <<
  "  fieldsToImport=\"{PERM}\""<<
  "  fieldNamesInGEOSX=\"{PERM}\""<<
  "  file=\"" << eclipseFilePath.c_str()<< "\"/>"<<
  "</Mesh>";
  const string inputStringMesh = inputStreamMesh.str();

  std::stringstream inputStreamRegion;
  inputStreamRegion <<
  "<?xml version=\"1.0\" ?>" <<
  "  <ElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
  "  <ElementRegion name=\"0\" cellBlocks=\"{0_HEX}\" materialList=\"{water rock}\"/>" <<
  "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  TestMeshImport( inputStringMesh, inputStringRegion, "" );
}

int main(int argc, char** argv)
{
  ::testing::InitGoogleTest(&argc, argv);

#ifdef GEOSX_USE_MPI

  MPI_Init(&argc,&argv);

  MPI_Comm_dup( MPI_COMM_WORLD, &MPI_COMM_GEOSX );

  logger::InitializeLogger(MPI_COMM_GEOSX);
#else
  logger::InitializeLogger():
#endif

  cxx_utilities::setSignalHandling(cxx_utilities::handler1);

  int const result = RUN_ALL_TESTS();

  logger::FinalizeLogger();

#ifdef GEOSX_USE_MPI
  MPI_Comm_free( &MPI_COMM_GEOSX );
  MPI_Finalize();
#endif

  return result;
}
