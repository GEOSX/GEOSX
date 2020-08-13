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
 * @file MeshLevel.cpp
 */

#include "MeshLevel.hpp"

#include "ElementRegionManager.hpp"
#include "NodeManager.hpp"
//#include "EdgeManager.hpp"
#include "FaceManager.hpp"

namespace geosx
{
using namespace dataRepository;

MeshLevel::MeshLevel( string const & name, Group * const parent ) :
  Group( name, parent ),
  m_nodeManager( groupStructKeys::nodeManagerString, this ),
  m_edgeManager( groupStructKeys::edgeManagerString, this ),
  m_faceManager( groupStructKeys::faceManagerString, this ),
  m_elementManager( groupStructKeys::elemManagerString, this )
{
  RegisterGroup( groupStructKeys::nodeManagerString, &m_nodeManager );

  RegisterGroup( groupStructKeys::edgeManagerString, &m_edgeManager );

  RegisterGroup< FaceManager >( groupStructKeys::faceManagerString, &m_faceManager );
  m_faceManager.nodeList().SetRelatedObject( &m_nodeManager );

  RegisterGroup< ElementRegionManager >( groupStructKeys::elemManagerString,
                                         &m_elementManager );

  registerWrapper< integer >( viewKeys.meshLevel );
}

MeshLevel::~MeshLevel()
{}

void
MeshLevel::InitializePostInitialConditions_PostSubGroups( Group * const )
{
  m_elementManager.forElementSubRegions< FaceElementSubRegion >(
    [&]( FaceElementSubRegion & subRegion ) {
      subRegion.CalculateElementGeometricQuantities( m_nodeManager, m_faceManager );
    } );
}

void
MeshLevel::GenerateAdjacencyLists(
  arrayView1d< localIndex const > const & seedNodeList,
  localIndex_array & nodeAdjacencyList,
  localIndex_array & edgeAdjacencyList,
  localIndex_array & faceAdjacencyList,
  ElementRegionManager::ElementViewAccessor< ReferenceWrapper< localIndex_array > > &
    elementAdjacencyList,
  integer const depth )
{
  NodeManager * const nodeManager = getNodeManager();

  ArrayOfArraysView< localIndex const > const & nodeToElementRegionList =
    nodeManager->elementRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementSubRegionList =
    nodeManager->elementSubRegionList().toViewConst();

  ArrayOfArraysView< localIndex const > const & nodeToElementList =
    nodeManager->elementList().toViewConst();

  FaceManager * const faceManager = this->getFaceManager();
  ArrayOfArraysView< localIndex const > const & faceToEdges =
    faceManager->edgeList().toViewConst();

  ElementRegionManager * const elemManager = this->getElemManager();

  std::set< localIndex > nodeAdjacencySet;
  std::set< localIndex > edgeAdjacencySet;
  std::set< localIndex > faceAdjacencySet;
  std::vector< std::vector< std::set< localIndex > > > elementAdjacencySet(
    elemManager->numRegions() );

  for( localIndex a = 0; a < elemManager->numRegions(); ++a )
  {
    elementAdjacencySet[a].resize( elemManager->GetRegion( a )->numSubRegions() );
  }

  nodeAdjacencySet.insert( seedNodeList.begin(), seedNodeList.end() );

  for( integer d = 0; d < depth; ++d )
  {
    for( localIndex const nodeIndex : nodeAdjacencySet )
    {
      for( localIndex b = 0; b < nodeToElementRegionList.sizeOfArray( nodeIndex );
           ++b )
      {
        localIndex const regionIndex = nodeToElementRegionList[nodeIndex][b];
        localIndex const subRegionIndex =
          nodeToElementSubRegionList[nodeIndex][b];
        localIndex const elementIndex = nodeToElementList[nodeIndex][b];
        elementAdjacencySet[regionIndex][subRegionIndex].insert( elementIndex );
      }
    }

    for( typename dataRepository::indexType kReg = 0;
         kReg < elemManager->numRegions();
         ++kReg )
    {
      ElementRegionBase const * const elemRegion = elemManager->GetRegion( kReg );

      elemRegion->forElementSubRegionsIndex< CellElementSubRegion, WellElementSubRegion >(
        [&]( localIndex const kSubReg, auto const & subRegion ) {
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes =
            subRegion.nodeList();
          arrayView2d< localIndex const > const & elemsToFaces =
            subRegion.faceList();
          for( auto const elementIndex : elementAdjacencySet[kReg][kSubReg] )
          {
            for( localIndex a = 0; a < elemsToNodes.size( 1 ); ++a )
            {
              nodeAdjacencySet.insert( elemsToNodes[elementIndex][a] );
            }

            for( localIndex a = 0; a < elemsToFaces.size( 1 ); ++a )
            {
              faceAdjacencySet.insert( elemsToFaces[elementIndex][a] );

              localIndex const faceID = elemsToFaces[elementIndex][a];
              localIndex const numEdges = faceToEdges.sizeOfArray( faceID );
              for( localIndex b = 0; b < numEdges; ++b )
              {
                edgeAdjacencySet.insert( faceToEdges( faceID, b ) );
              }
            }
          }
        } );
    }
  }

  nodeAdjacencyList.resize(
    LvArray::integerConversion< localIndex >( nodeAdjacencySet.size() ) );
  std::copy( nodeAdjacencySet.begin(),
             nodeAdjacencySet.end(),
             nodeAdjacencyList.begin() );

  edgeAdjacencyList.resize(
    LvArray::integerConversion< localIndex >( edgeAdjacencySet.size() ) );
  std::copy( edgeAdjacencySet.begin(),
             edgeAdjacencySet.end(),
             edgeAdjacencyList.begin() );

  faceAdjacencyList.resize(
    LvArray::integerConversion< localIndex >( faceAdjacencySet.size() ) );
  std::copy( faceAdjacencySet.begin(),
             faceAdjacencySet.end(),
             faceAdjacencyList.begin() );

  for( localIndex kReg = 0; kReg < elemManager->numRegions(); ++kReg )
  {
    ElementRegionBase const * const elemRegion = elemManager->GetRegion( kReg );

    for( localIndex kSubReg = 0; kSubReg < elemRegion->numSubRegions(); ++kSubReg )
    {
      elementAdjacencyList[kReg][kSubReg].get().clear();
      elementAdjacencyList[kReg][kSubReg].get().resize(
        LvArray::integerConversion< localIndex >(
          elementAdjacencySet[kReg][kSubReg].size() ) );
      std::copy( elementAdjacencySet[kReg][kSubReg].begin(),
                 elementAdjacencySet[kReg][kSubReg].end(),
                 elementAdjacencyList[kReg][kSubReg].get().begin() );
    }
  }
}

} /* namespace geosx */
