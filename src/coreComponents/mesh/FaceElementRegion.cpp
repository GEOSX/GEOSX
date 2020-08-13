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
 * @file FaceElementRegion.cpp
 */

#include "EdgeManager.hpp"
#include "FaceElementRegion.hpp"

namespace geosx
{
using namespace dataRepository;

FaceElementRegion::FaceElementRegion( string const & name, Group * const parent ) :
  ElementRegionBase( name, parent )
{
  this->GetGroup( viewKeyStruct::elementSubRegions )
    ->RegisterGroup< FaceElementSubRegion >( "default" );

  registerWrapper( viewKeyStruct::defaultApertureString, &m_defaultAperture )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "The default aperture of for new faceElements." );
}

FaceElementRegion::~FaceElementRegion()
{}

void
FaceElementRegion::InitializePreSubGroups( Group * const )
{
  this->forElementSubRegions< FaceElementSubRegion >(
    [&]( FaceElementSubRegion & subRegion ) {
      subRegion
        .getWrapper< array1d< real64 > >(
          FaceElementSubRegion::viewKeyStruct::elementApertureString )
        ->setApplyDefaultValue( m_defaultAperture );
    } );
}

localIndex
FaceElementRegion::AddToFractureMesh(
  real64 const time_np1,
  EdgeManager * const edgeManager,
  FaceManager const * const faceManager,
  ArrayOfArraysView< localIndex const > const & originalFaceToEdgeMap,
  string const & subRegionName,
  localIndex const faceIndices[2] )
{
  localIndex rval = -1;

  SortedArray< localIndex > connectedEdges;

  arrayView2d< localIndex const > const & faceToElementRegion =
    faceManager->elementRegionList();
  arrayView2d< localIndex const > const & faceToElementSubRegion =
    faceManager->elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementIndex =
    faceManager->elementList();

  Group * elementSubRegions = this->GetGroup( viewKeyStruct::elementSubRegions );

  FaceElementSubRegion * subRegion =
    elementSubRegions->GetGroup< FaceElementSubRegion >( subRegionName );
  subRegion->resize( subRegion->size() + 1 );
  rval = subRegion->size() - 1;

  arrayView1d< real64 > const & ruptureTime =
    subRegion->getReference< real64_array >( viewKeyStruct::ruptureTimeString );

  arrayView1d< real64 > const & creationMass = subRegion->getReference< real64_array >(
    FaceElementSubRegion::viewKeyStruct::creationMassString );

  arrayView2d< real64 const > const & faceCenter = faceManager->faceCenter();
  arrayView2d< real64 > const & elemCenter = subRegion->getElementCenter();
  arrayView1d< real64 const > const & elemArea = subRegion->getElementArea();

  arrayView1d< integer > const & subRegionGhostRank = subRegion->ghostRank();

  arrayView1d< integer const > const & faceGhostRank = faceManager->ghostRank();

  FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
  FaceElementSubRegion::EdgeMapType & edgeMap = subRegion->edgeList();
  FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    faceManager->nodeList().toViewConst();

  localIndex const kfe = subRegion->size() - 1;
  ruptureTime( kfe ) = time_np1;

  LvArray::tensorOps::copy< 3 >( elemCenter[kfe], faceCenter[faceIndices[0]] );

  faceMap[kfe][0] = faceIndices[0];
  faceMap[kfe][1] = faceIndices[1];
  globalIndex const gi = faceManager->localToGlobalMap()[faceIndices[0]];
  subRegion->localToGlobalMap()[kfe] = gi;
  subRegion->updateGlobalToLocalMap( kfe );
  subRegionGhostRank[kfe] = faceGhostRank[faceIndices[0]];

  // Add the nodes that compose the new FaceElement to the nodeList
  localIndex const numNodesInFace0 = faceToNodeMap.sizeOfArray( faceIndices[0] );
  localIndex const numNodesInFace1 = faceToNodeMap.sizeOfArray( faceIndices[1] );

  //Temporarily set the map size 8 for both quadrangle and triangle faces. TODO: need to fix for arbitrary face sizes.
  nodeMap[kfe].resize( 8 );

  for( localIndex a = 0; a < numNodesInFace0; ++a )
  {
    localIndex const aa = a < 2 ? a : numNodesInFace0 - a + 1;
    localIndex const bb = aa == 0 ? aa : numNodesInFace1 - aa;

    // TODO HACK need to generalize to something other than quads
    nodeMap[kfe][a] = faceToNodeMap( faceIndices[0], aa );
    nodeMap[kfe][a + numNodesInFace0] = faceToNodeMap( faceIndices[1], bb );
  }

  if( numNodesInFace0 == 3 )
  {
    nodeMap[kfe][6] = faceToNodeMap( faceIndices[0], 2 );
    nodeMap[kfe][7] = faceToNodeMap( faceIndices[1], 2 );
  }

  // Add the edges that compose the faceElement to the edge map. This is essentially a copy of
  // the facesToEdges entry.
  localIndex const faceID = faceIndices[0];
  localIndex const numEdges = originalFaceToEdgeMap.sizeOfArray( faceID );
  edgeMap[kfe].resize( numEdges );
  for( localIndex a = 0; a < numEdges; ++a )
  {
    edgeMap[kfe][a] = originalFaceToEdgeMap( faceID, a );
    connectedEdges.insert( originalFaceToEdgeMap( faceID, a ) );
  }

  // Add the cell region/subregion/index to the faceElementToCells map

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    subRegion->m_faceElementsToCells.m_toElementRegion[kfe][ke] =
      faceToElementRegion[faceIndices[ke]][ke];
    subRegion->m_faceElementsToCells.m_toElementSubRegion[kfe][ke] =
      faceToElementSubRegion[faceIndices[ke]][ke];
    subRegion->m_faceElementsToCells.m_toElementIndex[kfe][ke] =
      faceToElementIndex[faceIndices[ke]][ke];
  }

  // Fill the connectivity between FaceElement entries. This is essentially a copy of the
  // edgesToFaces map, but with differing offsets.
  for( auto const & edge : connectedEdges )
  {
    // check to see if the edgesToFractureConnectors already have an entry
    if( edgeManager->m_edgesToFractureConnectorsEdges.count( edge ) == 0 )
    {
      // if not, then fill increase the size of the fractureConnectors to face elements map and
      // fill the fractureConnectorsToEdges map with the current edge....and the inverse map too.
      edgeManager->m_fractureConnectorEdgesToFaceElements.appendArray( 0 );
      edgeManager->m_fractureConnectorsEdgesToEdges.emplace_back( edge );
      edgeManager->m_edgesToFractureConnectorsEdges[edge] =
        edgeManager->m_fractureConnectorsEdgesToEdges.size() - 1;
    }
    // now fill the fractureConnectorsToFaceElements map. This is analogous to the edge to face map
    localIndex const connectorIndex =
      edgeManager->m_edgesToFractureConnectorsEdges[edge];
    localIndex const numCells =
      edgeManager->m_fractureConnectorEdgesToFaceElements.sizeOfArray(
        connectorIndex ) +
      1;
    edgeManager->m_fractureConnectorEdgesToFaceElements.resizeArray( connectorIndex,
                                                                     numCells );
    edgeManager->m_fractureConnectorEdgesToFaceElements[connectorIndex][numCells - 1] =
      kfe;

    // And fill the list of connectors that will need stencil modifications
    edgeManager->m_recalculateFractureConnectorEdges.insert( connectorIndex );
  }

  subRegion->CalculateElementGeometricQuantities(
    kfe,
    faceManager->faceArea(),
    faceManager->faceRotationMatrix() );

  creationMass[kfe] *= elemArea[kfe];

  // update the sets
  for( auto const & setIter : faceManager->sets().wrappers() )
  {
    SortedArrayView< localIndex const > const & faceSet =
      faceManager->sets().getReference< SortedArray< localIndex > >( setIter.first );
    SortedArray< localIndex > & faceElementSet =
      subRegion->sets()
        .registerWrapper< SortedArray< localIndex > >( setIter.first )
        ->reference();
    for( localIndex a = 0; a < faceMap.size( 0 ); ++a )
    {
      localIndex const faceIndex = faceMap[a][0];
      if( faceSet.count( faceIndex ) )
      {
        faceElementSet.insert( a );
      }
    }
  }

  return rval;
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        FaceElementRegion,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
