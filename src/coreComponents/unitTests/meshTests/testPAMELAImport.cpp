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


// Source includes
#include "mainInterface/initialization.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "tests/meshFileNames.hpp"
#include "mesh/MeshManager.hpp"
#include "mesh/generators/CellBlockManager.hpp"

// TPL includes
#include <gtest/gtest.h>
#include <conduit.hpp>

using namespace geosx;
using namespace geosx::dataRepository;

void TestMeshImport( string const & inputStringMesh,
                     string const & inputStringRegion,
                     string const & propertyToTest )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

  // Load the mesh
  xmlWrapper::xmlDocument xmlDocument;
  xmlDocument.load_buffer( inputStringMesh.c_str(), inputStringMesh.size() );

  xmlWrapper::xmlNode xmlMeshNode = xmlDocument.child( "Mesh" );
  meshManager.processInputFileRecursive( xmlMeshNode );
  meshManager.postProcessInputRecursive();

  // Create the domain and generate the Mesh
  auto domain = std::unique_ptr< DomainPartition >( new DomainPartition( "domain", &root ) );
  meshManager.generateMeshes( *domain );

  MeshBody & meshBody = domain->getMeshBody( 0 );
  MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
  NodeManager & nodeManager = meshLevel.getNodeManager();
  FaceManager const & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // Create the ElementRegions
  xmlDocument.load_buffer( inputStringRegion.c_str(), inputStringRegion.size() );

  xmlWrapper::xmlNode xmlRegionNode = xmlDocument.child( "ElementRegions" );
  elemManager.processInputFileRecursive( xmlRegionNode );
  elemManager.postProcessInputRecursive();

  CellBlockManager & cellBlockManager = domain->getGroup< CellBlockManager >( keys::cellManager );
  nodeManager.setNodesInformation( cellBlockManager );

  // This method will call the CopyElementSubRegionFromCellBlocks that will trigger the property transfer.
  elemManager.generateMesh( cellBlockManager );


  // Check if the computed center match with the imported center
  if( !propertyToTest.empty() )
  {
    auto centerProperty = elemManager.constructArrayViewAccessor< real64, 2 >( propertyToTest );
    elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
      [&]( localIndex const er, localIndex const esr, ElementRegionBase &, ElementSubRegionBase & elemSubRegion )
    {
      elemSubRegion.calculateElementGeometricQuantities( nodeManager, faceManager );
      arrayView2d< real64 const > const elemCenter = elemSubRegion.getElementCenter();
      for( localIndex ei = 0; ei < elemSubRegion.size(); ei++ )
      {
        real64 center[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( elemCenter[ ei ] );
        LvArray::tensorOps::subtract< 3 >( center, centerProperty[er][esr][ei] );
        GEOSX_ERROR_IF_GT_MSG( LvArray::tensorOps::l2Norm< 3 >( center ), meshBody.getGlobalLengthScale() * 1e-8, "Property import of centers if wrong" );
      }
    } );
  }
}

TEST( PAMELAImport, testGMSH )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

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
    "<ElementRegions>" <<
    "  <CellElementRegion name=\"0\" cellBlocks=\"{Overburden1_TETRA}\" materialList=\"{water, rock}\"/>" <<
    "  <CellElementRegion name=\"1\" cellBlocks=\"{Overburden2_TETRA}\" materialList=\"{water, rock}\"/>" <<
    "  <CellElementRegion name=\"2\" cellBlocks=\"{Reservoir_TETRA}\" materialList=\"{water, rock}\"/>" <<
    "  <CellElementRegion name=\"3\" cellBlocks=\"{Underburden_TETRA}\" materialList=\"{water, rock}\"/>" <<
    "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  TestMeshImport( inputStringMesh, inputStringRegion, "barycenter" );
}

TEST( PAMELAImport, testECLIPSE )
{
  conduit::Node node;
  Group root( "root", node );
  MeshManager meshManager( "mesh", &root );

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
    "  <CellElementRegions xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:noNamespaceSchemaLocation=\"geos_v0.0.xsd\">" <<
    "  <CellElementRegion name=\"0\" cellBlocks=\"{DEFAULT_HEX}\" materialList=\"{water, rock}\"/>" <<
    "</ElementRegions>";
  string inputStringRegion = inputStreamRegion.str();

  TestMeshImport( inputStringMesh, inputStringRegion, "" );
}

int main( int argc, char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
