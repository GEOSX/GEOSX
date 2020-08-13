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
#include "managers/initialization.hpp"
#include "common/Logger.hpp"
#include "mesh/FaceManager.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "finiteVolume/HybridFVMInnerProduct.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::HybridFVMInnerProduct;
using namespace geosx::computationalGeometry;
using namespace geosx::testing;

struct InnerProductType
{
  static constexpr integer TPFA = 0;
  static constexpr integer QUASI_TPFA = 1;
};

void
compareTransmissibilityMatrices(
  arraySlice2d< real64 const > const & transMatrix,
  arraySlice2d< real64 const > const & transMatrixRef )
{
  for( localIndex ifaceLoc = 0; ifaceLoc < transMatrix.size( 0 ); ++ifaceLoc )
  {
    for( localIndex jfaceLoc = 0; jfaceLoc < transMatrix.size( 1 ); ++jfaceLoc )
    {
      checkRelativeError( transMatrix( ifaceLoc, jfaceLoc ),
                          transMatrixRef( ifaceLoc, jfaceLoc ),
                          1e-10 );
    }
  }
}

void
computeVolumeAndCenter(
  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & nodePosition,
  array1d< localIndex > const & toNodes,
  R1Tensor & elemCenter,
  real64 & elemVolume )
{
  localIndex const numNodes = toNodes.size();
  R1Tensor Xlocal[10];
  elemCenter( 0 ) = 0;
  elemCenter( 1 ) = 0;
  elemCenter( 2 ) = 0;
  for( localIndex a = 0; a < numNodes; ++a )
  {
    Xlocal[a][0] = nodePosition( toNodes( a ), 0 );
    Xlocal[a][1] = nodePosition( toNodes( a ), 1 );
    Xlocal[a][2] = nodePosition( toNodes( a ), 2 );
    elemCenter += Xlocal[a];
  }
  elemCenter /= numNodes;

  GEOSX_ERROR_IF( numNodes != 8 && numNodes != 4,
                  "This number of nodes is not supported in the test yet" );

  if( numNodes == 8 )
  {
    elemVolume = computationalGeometry::HexVolume( Xlocal );
  }
  else if( numNodes == 4 )
  {
    elemVolume = computationalGeometry::TetVolume( Xlocal );
  }
}

void
makeHexa( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
          FaceManager::NodeMapType & faceToNodes,
          array1d< localIndex > & elemToFaces,
          R1Tensor & elemCenter,
          real64 & elemVolume,
          R1Tensor & elemPerm,
          real64 & lengthTolerance,
          integer const ipType,
          arraySlice2d< real64 > const & transMatrixRef )
{
  localIndex constexpr numNodes = 8;
  localIndex constexpr numFaces = 6;
  localIndex constexpr numNodesPerFace = 4;

  lengthTolerance = 1.73205e-8;

  elemPerm( 0 ) = 1e-12;
  elemPerm( 1 ) = 2e-12;
  elemPerm( 2 ) = 3e-12;

  // elem-to-faces map
  elemToFaces.resize( numFaces );
  elemToFaces( 0 ) = 1;
  elemToFaces( 1 ) = 2;
  elemToFaces( 2 ) = 0;
  elemToFaces( 3 ) = 5;
  elemToFaces( 4 ) = 4;
  elemToFaces( 5 ) = 3;

  // face-to-nodes map
  faceToNodes.resize( numFaces );
  for( localIndex iface = 0; iface < numFaces; ++iface )
  {
    faceToNodes.resizeArray( iface, numNodesPerFace );
  }
  faceToNodes( 0, 0 ) = 0;
  faceToNodes( 0, 1 ) = 1;
  faceToNodes( 0, 2 ) = 3;
  faceToNodes( 0, 3 ) = 2;
  faceToNodes( 1, 0 ) = 0;
  faceToNodes( 1, 1 ) = 4;
  faceToNodes( 1, 2 ) = 5;
  faceToNodes( 1, 3 ) = 1;
  faceToNodes( 2, 0 ) = 0;
  faceToNodes( 2, 1 ) = 2;
  faceToNodes( 2, 2 ) = 6;
  faceToNodes( 2, 3 ) = 4;
  faceToNodes( 3, 0 ) = 1;
  faceToNodes( 3, 1 ) = 5;
  faceToNodes( 3, 2 ) = 7;
  faceToNodes( 3, 3 ) = 3;
  faceToNodes( 4, 0 ) = 6;
  faceToNodes( 4, 1 ) = 2;
  faceToNodes( 4, 2 ) = 3;
  faceToNodes( 4, 3 ) = 7;
  faceToNodes( 5, 0 ) = 4;
  faceToNodes( 5, 1 ) = 6;
  faceToNodes( 5, 2 ) = 7;
  faceToNodes( 5, 3 ) = 5;

  // node position
  nodePosition.resize( numNodes, 3 );
  nodePosition( 0, 0 ) = 0;
  nodePosition( 0, 1 ) = 0;
  nodePosition( 0, 2 ) = 0;
  nodePosition( 1, 0 ) = 0.75;
  nodePosition( 1, 1 ) = 0;
  nodePosition( 1, 2 ) = 1;
  nodePosition( 2, 0 ) = 0;
  nodePosition( 2, 1 ) = 1;
  nodePosition( 2, 2 ) = 0;
  nodePosition( 3, 0 ) = 0.75;
  nodePosition( 3, 1 ) = 1;
  nodePosition( 3, 2 ) = 1;
  nodePosition( 4, 0 ) = 1;
  nodePosition( 4, 1 ) = 0;
  nodePosition( 4, 2 ) = 0;
  nodePosition( 5, 0 ) = 1.75;
  nodePosition( 5, 1 ) = 0;
  nodePosition( 5, 2 ) = 1;
  nodePosition( 6, 0 ) = 1;
  nodePosition( 6, 1 ) = 1;
  nodePosition( 6, 2 ) = 0;
  nodePosition( 7, 0 ) = 1.75;
  nodePosition( 7, 1 ) = 1;
  nodePosition( 7, 2 ) = 1;

  array1d< localIndex > toNodes;
  toNodes.resize( numNodes );
  toNodes( 0 ) = 0;
  toNodes( 1 ) = 4;
  toNodes( 2 ) = 2;
  toNodes( 3 ) = 6;
  toNodes( 4 ) = 1;
  toNodes( 5 ) = 5;
  toNodes( 6 ) = 3;
  toNodes( 7 ) = 7;

  computeVolumeAndCenter( nodePosition, toNodes, elemCenter, elemVolume );

  // reference transmissibility matrix computed with MRST
  if( ipType == InnerProductType::TPFA )
  {
    transMatrixRef( 0, 0 ) = 4e-12;
    transMatrixRef( 1, 1 ) = 3.84e-12;
    transMatrixRef( 2, 2 ) = 2e-12;
    transMatrixRef( 3, 3 ) = 2e-12;
    transMatrixRef( 4, 4 ) = 4e-12;
    transMatrixRef( 5, 5 ) = 3.84e-12;
  }
  else if( ipType == InnerProductType::QUASI_TPFA )
  {
    transMatrixRef( 0, 0 ) = 4e-12;

    transMatrixRef( 1, 1 ) = 6e-12;
    transMatrixRef( 1, 2 ) = -2.250e-12;
    transMatrixRef( 1, 3 ) = 2.250e-12;

    transMatrixRef( 2, 1 ) = -2.225e-12;
    transMatrixRef( 2, 2 ) = 5.375e-12;
    transMatrixRef( 2, 5 ) = 2.225e-12;

    transMatrixRef( 3, 1 ) = 2.25e-12;
    transMatrixRef( 3, 3 ) = 5.375e-12;
    transMatrixRef( 3, 5 ) = -2.25e-12;

    transMatrixRef( 4, 4 ) = 4e-12;

    transMatrixRef( 5, 2 ) = 2.25e-12;
    transMatrixRef( 5, 5 ) = 6e-12;
    transMatrixRef( 5, 3 ) = -2.25e-12;
  }
}

void
makeTetra( array2d< real64, nodes::REFERENCE_POSITION_PERM > & nodePosition,
           FaceManager::NodeMapType & faceToNodes,
           array1d< localIndex > & elemToFaces,
           R1Tensor & elemCenter,
           real64 & elemVolume,
           R1Tensor & elemPerm,
           real64 & lengthTolerance,
           integer const ipType,
           arraySlice2d< real64 > const & transMatrixRef )
{
  localIndex constexpr numNodes = 4;
  localIndex constexpr numFaces = 4;
  localIndex constexpr numNodesPerFace = 3;

  lengthTolerance = 1.73205e-8;

  elemPerm( 0 ) = 1e-12;
  elemPerm( 1 ) = 2e-12;
  elemPerm( 2 ) = 3e-12;

  // elem-to-faces map
  elemToFaces.resize( numFaces );
  elemToFaces( 0 ) = 0;
  elemToFaces( 1 ) = 1;
  elemToFaces( 2 ) = 2;
  elemToFaces( 3 ) = 3;

  // face-to-nodes map
  faceToNodes.resize( numFaces );
  for( localIndex iface = 0; iface < numFaces; ++iface )
  {
    faceToNodes.resizeArray( iface, numNodesPerFace );
  }
  faceToNodes( 0, 0 ) = 0;
  faceToNodes( 0, 1 ) = 1;
  faceToNodes( 0, 2 ) = 2;
  faceToNodes( 1, 0 ) = 0;
  faceToNodes( 1, 1 ) = 1;
  faceToNodes( 1, 2 ) = 3;
  faceToNodes( 2, 0 ) = 1;
  faceToNodes( 2, 1 ) = 2;
  faceToNodes( 2, 2 ) = 3;
  faceToNodes( 3, 0 ) = 2;
  faceToNodes( 3, 1 ) = 0;
  faceToNodes( 3, 2 ) = 3;

  // node position
  nodePosition.resize( numNodes, 3 );
  nodePosition( 0, 0 ) = 0;
  nodePosition( 0, 1 ) = 0;
  nodePosition( 0, 2 ) = 0;
  nodePosition( 1, 0 ) = 1;
  nodePosition( 1, 1 ) = 0;
  nodePosition( 1, 2 ) = 0;
  nodePosition( 2, 0 ) = 0;
  nodePosition( 2, 1 ) = 1;
  nodePosition( 2, 2 ) = 0;
  nodePosition( 3, 0 ) = 0;
  nodePosition( 3, 1 ) = 0;
  nodePosition( 3, 2 ) = 3;

  array1d< localIndex > toNodes;
  toNodes.resize( numNodes );
  toNodes( 0 ) = 0;
  toNodes( 1 ) = 1;
  toNodes( 2 ) = 2;
  toNodes( 3 ) = 3;

  computeVolumeAndCenter( nodePosition, toNodes, elemCenter, elemVolume );

  // reference transmissibility matrix computed with MRST

  if( ipType == InnerProductType::TPFA )
  {
    transMatrixRef( 0, 0 ) = 1.952e-12;
    transMatrixRef( 1, 1 ) = 5.684e-12;
    transMatrixRef( 2, 2 ) = 9.818e-12;
    transMatrixRef( 3, 3 ) = 2.842e-12;
  }
  else if( ipType == InnerProductType::QUASI_TPFA )
  {
    transMatrixRef( 0, 0 ) = 5.25e-12;
    transMatrixRef( 0, 1 ) = 3.75e-12;
    transMatrixRef( 0, 2 ) = 2.25e-12;
    transMatrixRef( 0, 3 ) = 3.75e-12;

    transMatrixRef( 1, 0 ) = 3.75e-12;
    transMatrixRef( 1, 1 ) = 12.75e-12;
    transMatrixRef( 1, 2 ) = -5.25e-12;
    transMatrixRef( 1, 3 ) = 3.75e-12;
    transMatrixRef( 2, 0 ) = 2.25e-12;
    transMatrixRef( 2, 1 ) = -5.25e-12;
    transMatrixRef( 2, 2 ) = 18.75e-12;
    transMatrixRef( 2, 3 ) = -0.75e-12;

    transMatrixRef( 3, 0 ) = 3.75e-12;
    transMatrixRef( 3, 1 ) = 3.75e-12;
    transMatrixRef( 3, 2 ) = -0.75e-12;
    transMatrixRef( 3, 3 ) = 8.25e-12;
  }
}

TEST( testHybridFVMInnerProducts, TPFA_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  R1Tensor elemCenter;
  R1Tensor elemPerm;
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF * NF > transMatrixRef( NF, NF );

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::TPFA,
            transMatrixRef );

  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  TPFACellInnerProductKernel::Compute< NF >( nodePosition.toViewConst(),
                                             faceToNodes.toViewConst(),
                                             elemToFaces.toSliceConst(),
                                             center,
                                             perm,
                                             lengthTolerance,
                                             transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );
}

TEST( testHybridFVMInnerProducts, QTPFA_hexa )
{
  localIndex constexpr NF = 6;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  R1Tensor elemCenter;
  R1Tensor elemPerm;
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF * NF > transMatrixRef( NF, NF );

  real64 const tParam = 2;

  makeHexa( nodePosition,
            faceToNodes,
            elemToFaces,
            elemCenter,
            elemVolume,
            elemPerm,
            lengthTolerance,
            InnerProductType::QUASI_TPFA,
            transMatrixRef );

  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  QTPFACellInnerProductKernel::Compute< NF >( nodePosition.toViewConst(),
                                              faceToNodes.toViewConst(),
                                              elemToFaces.toSliceConst(),
                                              center,
                                              elemVolume,
                                              perm,
                                              tParam,
                                              lengthTolerance,
                                              transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );
}

TEST( testHybridFVMInnerProducts, TPFA_tetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  R1Tensor elemCenter;
  R1Tensor elemPerm;
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF * NF > transMatrixRef( NF, NF );

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::TPFA,
             transMatrixRef );

  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  TPFACellInnerProductKernel::Compute< NF >( nodePosition.toViewConst(),
                                             faceToNodes.toViewConst(),
                                             elemToFaces.toSliceConst(),
                                             center,
                                             perm,
                                             lengthTolerance,
                                             transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );
}

TEST( testHybridFVMInnerProducts, QTPFA_tetra )
{
  localIndex constexpr NF = 4;

  array2d< real64, nodes::REFERENCE_POSITION_PERM > nodePosition;
  FaceManager::NodeMapType faceToNodes;
  array1d< localIndex > elemToFaces;
  R1Tensor elemCenter;
  R1Tensor elemPerm;
  real64 elemVolume = 0;
  real64 lengthTolerance = 0;
  stackArray2d< real64, NF * NF > transMatrixRef( NF, NF );

  real64 const tParam = 2;

  makeTetra( nodePosition,
             faceToNodes,
             elemToFaces,
             elemCenter,
             elemVolume,
             elemPerm,
             lengthTolerance,
             InnerProductType::QUASI_TPFA,
             transMatrixRef );

  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  stackArray1d< real64, 3 > center( 3 );
  center[0] = elemCenter[0];
  center[1] = elemCenter[1];
  center[2] = elemCenter[2];
  real64 const perm[3] = { elemPerm[0], elemPerm[1], elemPerm[2] };

  QTPFACellInnerProductKernel::Compute< NF >( nodePosition.toViewConst(),
                                              faceToNodes.toViewConst(),
                                              elemToFaces.toSliceConst(),
                                              center,
                                              elemVolume,
                                              perm,
                                              tParam,
                                              lengthTolerance,
                                              transMatrix.toSlice() );

  compareTransmissibilityMatrices( transMatrix, transMatrixRef );
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
