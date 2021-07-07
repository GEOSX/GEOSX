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
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mainInterface/initialization.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/generators/CellBlockManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::constitutive;
using namespace geosx::dataRepository;


void RegisterAndApplyField( DomainPartition & domain,
                            string const & fieldName,
                            string const & objectPath,
                            real64 value )
{
  FieldSpecificationManager & fieldSpecificationManager = FieldSpecificationManager::getInstance();

  FieldSpecificationBase & fieldSpec = fieldSpecificationManager.registerGroup< FieldSpecificationBase >( fieldName );
  fieldSpec.setFieldName( fieldName );
  fieldSpec.setObjectPath( objectPath );
  fieldSpec.setScale( value );
  fieldSpec.initialCondition( true );
  fieldSpec.addSetName( "all" );

  fieldSpecificationManager.apply( 0., domain, "", "",
                                   [&] ( FieldSpecificationBase const & bc,
                                         string const &,
                                         SortedArrayView< localIndex const > const & targetSet,
                                         Group & targetGroup,
                                         string const name )
  {
    bc.applyFieldValue< FieldSpecificationEqual >( targetSet, 0.0, targetGroup, name );
  } );
}

TEST( FieldSpecification, Recursive )
{
  // Mesh Definitions
  localIndex const nbTetReg0 = 30;
  localIndex const nbHexReg0 = 60;
  localIndex const nbTetReg1 = 40;
  localIndex const nbHexReg1 = 50;

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( "body" );
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));

  ElementRegionManager & elemManager = meshLevel0.getElemManager();
  CellElementRegion & reg0 = dynamicCast< CellElementRegion & >( *elemManager.createChild( "CellElementRegion", "reg0" ) );
  CellElementRegion & reg1 = dynamicCast< CellElementRegion & >( *elemManager.createChild( "CellElementRegion", "reg1" ) );

  // Cell blocks should only be used to define the sub regions.
  // This scope protection is there to make them disappear from the rest of the test.
  {
    CellBlockManager & cellBlockManager = domain.registerGroup< CellBlockManager >( keys::cellManager );

    CellBlock & reg0Hex = cellBlockManager.registerCellBlock( "reg0hex" );
    reg0Hex.setElementType( "C3D8" );
    reg0Hex.resize( nbHexReg0 );

    CellBlock & reg0Tet = cellBlockManager.registerCellBlock( "reg0tet" );
    reg0Tet.setElementType( "C3D4" );
    reg0Tet.resize( nbTetReg0 );

    CellBlock & reg1Hex = cellBlockManager.registerCellBlock( "reg1hex" );
    reg1Hex.setElementType( "C3D8" );
    reg1Hex.resize( nbHexReg1 );

    CellBlock & reg1Tet = cellBlockManager.registerCellBlock( "reg1tet" );
    reg1Tet.setElementType( "C3D4" );
    reg1Tet.resize( nbTetReg1 );

    reg0.addCellBlockName( reg0Hex.getName() );
    reg0.addCellBlockName( reg0Tet.getName() );
    reg0.generateMesh( cellBlockManager.getCellBlocks() );

    reg1.addCellBlockName( reg1Hex.getName() );
    reg1.addCellBlockName( reg1Tet.getName() );
    reg1.generateMesh( cellBlockManager.getCellBlocks() );

    // The cell block manager should not be used anymore.
    domain.deregisterGroup( keys::cellManager );
  }

  /// Field Definition
  reg0.getSubRegion( "reg0hex" ).registerWrapper< array1d< real64 > >( "field0" );
  reg0.getSubRegion( "reg0tet" ).registerWrapper< array1d< real64 > >( "field0" );
  reg1.getSubRegion( "reg1tet" ).registerWrapper< array1d< real64 > >( "field0" );
  reg1.getSubRegion( "reg1hex" ).registerWrapper< array1d< real64 > >( "field0" );

  reg0.getSubRegion( "reg0hex" ).registerWrapper< array1d< real64 > >( "field1" );
  reg0.getSubRegion( "reg0tet" ).registerWrapper< array1d< real64 > >( "field1" );

  reg0.getSubRegion( "reg0hex" ).registerWrapper< array1d< real64 > >( "field2" );

  reg1.getSubRegion( "reg1tet" ).registerWrapper< array1d< real64 > >( "field3" );

  SortedArray< localIndex > & set0hex = reg0.getSubRegion( "reg0hex" )
                                          .getGroup( "sets" )
                                          .registerWrapper< SortedArray< localIndex > >( string( "all" ) )
                                          .reference();
  for( localIndex i = 0; i < nbHexReg0; i++ )
  {
    set0hex.insert( i );
  }

  SortedArray< localIndex > & set0tet = reg0.getSubRegion( "reg0tet" )
                                          .getGroup( "sets" )
                                          .registerWrapper< SortedArray< localIndex > >( string( "all" ) )
                                          .reference();
  for( localIndex i = 0; i < nbTetReg0; i++ )
  {
    set0tet.insert( i );
  }

  SortedArray< localIndex > & set1hex = reg1.getSubRegion( "reg1hex" )
                                          .getGroup( "sets" )
                                          .registerWrapper< SortedArray< localIndex > >( string( "all" ) )
                                          .reference();
  for( localIndex i = 0; i < nbHexReg1; i++ )
  {
    set1hex.insert( i );
  }

  SortedArray< localIndex > & set1tet = reg1.getSubRegion( "reg1tet" )
                                          .getGroup( "sets" )
                                          .registerWrapper< SortedArray< localIndex > >( string( "all" ) )
                                          .reference();
  for( localIndex i = 0; i < nbTetReg1; i++ )
  {
    set1tet.insert( i );
  }

  RegisterAndApplyField( domain, "field0", "ElementRegions", 1. );
  RegisterAndApplyField( domain, "field1", "ElementRegions/reg0", 2. );
  RegisterAndApplyField( domain, "field2", "ElementRegions/reg0/elementSubRegions/reg0hex", 3. );
  RegisterAndApplyField( domain, "field3", "ElementRegions/reg1/elementSubRegions/reg1tet", 4. );

  /// Check if the values are well set

  auto field0 = elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "field0" );
  auto field1 = elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "field1" );
  auto field2 = elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "field2" );
  auto field3 = elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( "field3" );
  elemManager.forElementSubRegionsComplete< ElementSubRegionBase >(
    [&] ( localIndex const er, localIndex const esr, ElementRegionBase const &, ElementSubRegionBase const & subRegion )
  {
    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      real64 const & value = field0[er][esr][ei];
      // This `value < x || value > x` is there because of compiler warning (converted to error)
      // making floating points comparisons using `==` impossible even for integers...
      GEOSX_ERROR_IF( value< 1. || value > 1., "Recursive fields are not set, value is " + std::to_string( value ) + "." );
    } );
  } );

  reg0.forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase & subRegion )
  {
    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      real64 const & value = field1[0][esr][ei];
      GEOSX_ERROR_IF( value< 2. || value > 2., "Recursive fields are not set, value is " + std::to_string( value ) + "." );
    } );
  } );

  forAll< serialPolicy >( nbHexReg0, [=] ( localIndex const ei )
  {
    real64 const & value = field2[0][0][ei];
    GEOSX_ERROR_IF( value< 3. || value > 3., "Recursive fields are not set, value is " + std::to_string( value ) + "." );
  } );

  forAll< serialPolicy >( nbTetReg1, [=] ( localIndex const ei )
  {
    real64 const & value = field3[1][1][ei];
    GEOSX_ERROR_IF( value< 4. || value > 4., "Recursive fields are not set, value is " + std::to_string( value ) + "." );
  } );
}

int main( int argc, char * * argv )
{
  GeosxState state( basicSetup( argc, argv ) );

  ::testing::InitGoogleTest( &argc, argv );
  int const result = RUN_ALL_TESTS();

  basicCleanup();

  return result;
}
