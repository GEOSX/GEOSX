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

// Source includes
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/initialization.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;

void
RegisterAndApplyField( DomainPartition * domain,
                       string const & fieldName,
                       string const & objectPath,
                       real64 value )
{
  FieldSpecificationManager & fieldSpecificationManager =
    FieldSpecificationManager::get();

  auto fieldSpec =
    fieldSpecificationManager.RegisterGroup< FieldSpecificationBase >( fieldName );
  fieldSpec->SetFieldName( fieldName );
  fieldSpec->SetObjectPath( objectPath );
  fieldSpec->SetScale( value );
  fieldSpec->InitialCondition( true );
  fieldSpec->AddSetName( "all" );

  fieldSpecificationManager.Apply(
    0.,
    domain,
    "",
    "",
    [&]( FieldSpecificationBase const * const bc,
         string const &,
         SortedArrayView< localIndex const > const & targetSet,
         Group * const targetGroup,
         string const name ) {
      bc->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                      0.0,
                                                      targetGroup,
                                                      name );
    } );
}

TEST( FieldSpecification, Recursive )
{
  // Mesh Definitions
  localIndex nbTetReg0 = 30;
  localIndex nbHexReg0 = 60;
  localIndex nbTetReg1 = 40;
  localIndex nbHexReg1 = 50;
  auto domain =
    std::unique_ptr< DomainPartition >( new DomainPartition( "domain", nullptr ) );
  auto meshBodies = domain->getMeshBodies();
  MeshBody * const meshBody = meshBodies->RegisterGroup< MeshBody >( "body" );
  MeshLevel * const meshLevel0 =
    meshBody->RegisterGroup< MeshLevel >( std::string( "Level0" ) );

  CellBlockManager * cellBlockManager =
    domain->GetGroup< CellBlockManager >( keys::cellManager );

  CellBlock * reg0Hex = cellBlockManager->GetGroup( keys::cellBlocks )
                          ->RegisterGroup< CellBlock >( "reg0hex" );
  reg0Hex->SetElementType( "C3D8" );
  reg0Hex->resize( nbHexReg0 );
  auto & cellToVertexreg0Hex = reg0Hex->nodeList();
  cellToVertexreg0Hex.resize( nbHexReg0, 8 );

  CellBlock * reg0Tet = cellBlockManager->GetGroup( keys::cellBlocks )
                          ->RegisterGroup< CellBlock >( "reg0tet" );
  reg0Tet->SetElementType( "C3D4" );
  reg0Tet->resize( nbTetReg0 );
  auto & cellToVertexreg0Tet = reg0Tet->nodeList();
  cellToVertexreg0Tet.resize( nbTetReg0, 4 );

  CellBlock * reg1Hex = cellBlockManager->GetGroup( keys::cellBlocks )
                          ->RegisterGroup< CellBlock >( "reg1hex" );
  reg1Hex->SetElementType( "C3D8" );
  reg1Hex->resize( nbHexReg1 );
  auto & cellToVertexreg1Hex = reg1Hex->nodeList();
  cellToVertexreg1Hex.resize( nbHexReg1, 8 );

  CellBlock * reg1Tet = cellBlockManager->GetGroup( keys::cellBlocks )
                          ->RegisterGroup< CellBlock >( "reg1tet" );
  reg1Tet->SetElementType( "C3D4" );
  reg1Tet->resize( nbTetReg1 );
  auto & cellToVertexreg1Tet = reg1Tet->nodeList();
  cellToVertexreg1Tet.resize( nbTetReg1, 4 );

  ElementRegionManager * elemManager = meshLevel0->getElemManager();
  CellElementRegion * reg0 =
    elemManager->CreateChild( "CellElementRegion", "reg0" )
      ->group_cast< CellElementRegion * >();
  reg0->AddCellBlockName( reg0Hex->getName() );
  reg0->AddCellBlockName( reg0Tet->getName() );
  CellElementRegion * reg1 =
    elemManager->CreateChild( "CellElementRegion", "reg1" )
      ->group_cast< CellElementRegion * >();
  reg1->AddCellBlockName( reg1Hex->getName() );
  reg1->AddCellBlockName( reg1Tet->getName() );
  reg0->GenerateMesh( cellBlockManager->GetGroup( keys::cellBlocks ) );
  reg1->GenerateMesh( cellBlockManager->GetGroup( keys::cellBlocks ) );

  /// Field Definition
  reg0->GetSubRegion( "reg0hex" )->registerWrapper< array1d< real64 > >( "field0" );
  reg0->GetSubRegion( "reg0tet" )->registerWrapper< array1d< real64 > >( "field0" );
  reg1->GetSubRegion( "reg1tet" )->registerWrapper< array1d< real64 > >( "field0" );
  reg1->GetSubRegion( "reg1hex" )->registerWrapper< array1d< real64 > >( "field0" );

  reg0->GetSubRegion( "reg0hex" )->registerWrapper< array1d< real64 > >( "field1" );
  reg0->GetSubRegion( "reg0tet" )->registerWrapper< array1d< real64 > >( "field1" );

  reg0->GetSubRegion( "reg0hex" )->registerWrapper< array1d< real64 > >( "field2" );

  reg1->GetSubRegion( "reg1tet" )->registerWrapper< array1d< real64 > >( "field3" );

  SortedArray< localIndex > & set0hex =
    reg0->GetSubRegion( "reg0hex" )
      ->GetGroup( "sets" )
      ->registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();
  for( localIndex i = 0; i < nbHexReg0; i++ )
  {
    set0hex.insert( i );
  }

  SortedArray< localIndex > & set0tet =
    reg0->GetSubRegion( "reg0tet" )
      ->GetGroup( "sets" )
      ->registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();
  for( localIndex i = 0; i < nbTetReg0; i++ )
  {
    set0tet.insert( i );
  }

  SortedArray< localIndex > & set1hex =
    reg1->GetSubRegion( "reg1hex" )
      ->GetGroup( "sets" )
      ->registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();
  for( localIndex i = 0; i < nbHexReg1; i++ )
  {
    set1hex.insert( i );
  }

  SortedArray< localIndex > & set1tet =
    reg1->GetSubRegion( "reg1tet" )
      ->GetGroup( "sets" )
      ->registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();
  for( localIndex i = 0; i < nbTetReg1; i++ )
  {
    set1tet.insert( i );
  }

  RegisterAndApplyField( domain.get(), "field0", "ElementRegions", 1. );
  RegisterAndApplyField( domain.get(), "field1", "ElementRegions/reg0", 2. );
  RegisterAndApplyField( domain.get(),
                         "field2",
                         "ElementRegions/reg0/elementSubRegions/reg0hex",
                         3. );
  RegisterAndApplyField( domain.get(),
                         "field3",
                         "ElementRegions/reg1/elementSubRegions/reg1tet",
                         4. );

  /// Check if the values are well set

  auto field0 =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >(
      "field0" );
  auto field1 =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >(
      "field1" );
  auto field2 =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >(
      "field2" );
  auto field3 =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >(
      "field3" );
  elemManager->forElementSubRegionsComplete< ElementSubRegionBase >(
    [&]( localIndex const er,
         localIndex const esr,
         ElementRegionBase const &,
         ElementSubRegionBase const & subRegion ) {
      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei ) {
        GEOSX_ERROR_IF( field0[er][esr][ei] < 1. || field0[er][esr][ei] > 1.,
                        "Recursive fields are not set" );
      } );
    } );

  reg0->forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase & subRegion ) {
      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const ei ) {
        GEOSX_ERROR_IF( field1[0][esr][ei] < 2. || field1[0][esr][ei] > 2.,
                        "Recursive fields are not set" );
      } );
    } );

  forAll< serialPolicy >( reg0Hex->size(), [=]( localIndex const ei ) {
    GEOSX_ERROR_IF( field2[0][0][ei] < 1. || field2[0][0][ei] > 3.,
                    "Recursive fields are not set" );
  } );

  forAll< serialPolicy >( reg1Tet->size(), [=]( localIndex const ei ) {
    GEOSX_ERROR_IF( field3[1][1][ei] < 4. || field3[1][1][ei] > 4.,
                    "Recursive fields are not set" );
  } );
}

int
main( int argc, char ** argv )
{
  basicSetup( argc, argv );

  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
