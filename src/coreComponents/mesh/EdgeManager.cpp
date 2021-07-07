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

/**
 * @file EdgeManager.cpp
 */

#include "EdgeManager.hpp"

#include "BufferOps.hpp"
#include "NodeManager.hpp"
#include "FaceManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "mesh/generators/CellBlockUtilities.hpp"

namespace geosx
{
using namespace dataRepository;

EdgeManager::EdgeManager( string const & name,
                          Group * const parent ):
  ObjectManagerBase( name, parent ),
  m_edgesToFractureConnectorsEdges(),
  m_fractureConnectorsEdgesToEdges(),
  m_fractureConnectorEdgesToFaceElements()
{
  this->registerWrapper( viewKeyStruct::nodeListString(), &this->m_toNodesRelation );
  this->registerWrapper( viewKeyStruct::faceListString(), &this->m_toFacesRelation );

  m_toNodesRelation.resize( 0, 2 );

  registerWrapper( viewKeyStruct::edgesTofractureConnectorsEdgesString(), &m_edgesToFractureConnectorsEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of edge local indices to the fracture connector local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorEdgesToEdgesString(), &m_fractureConnectorsEdgesToEdges ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices to edge local indices." ).
    setSizedFromParent( 0 );

  registerWrapper( viewKeyStruct::fractureConnectorsEdgesToFaceElementsIndexString(),
                   &m_fractureConnectorEdgesToFaceElements ).
    setPlotLevel( PlotLevel::NOPLOT ).
    setDescription( "A map of fracture connector local indices face element local indices" ).
    setSizedFromParent( 0 );

}

EdgeManager::~EdgeManager()
{}

void EdgeManager::resize( localIndex const newSize )
{
  // FIXME My tests say that this line could be commented out.
  //       But those tests are not exhaustive.
  m_toFacesRelation.resize( newSize, 2 * faceMapExtraSpacePerEdge() );
  ObjectManagerBase::resize( newSize );
}

void EdgeManager::buildEdges( CellBlockManagerABC const & cellBlockManager,
                              NodeManager & nodeManager,
                              FaceManager & faceManager )
{
  GEOSX_MARK_FUNCTION;

  faceManager.edgeList().setRelatedObject( *this );

  m_toNodesRelation.setRelatedObject( nodeManager );
  m_toFacesRelation.setRelatedObject( faceManager );

  resize( cellBlockManager.numEdges() );

  faceManager.edgeList().base() = cellBlockManager.getFaceToEdges();
  // FIXME I don't know what I am doing...
  auto const & edgeToFaces = cellBlockManager.getEdgeToFaces();
  m_toFacesRelation.base().reserve( edgeToFaces.size() );
  m_toFacesRelation.base().reserveValues( edgeToFaces.valueCapacity() );
  m_toFacesRelation.base() = edgeToFaces;

  m_toNodesRelation.base() = cellBlockManager.getEdgeToNodes();

  // make sets from nodesets
  auto const & nodeSets = nodeManager.sets().wrappers();
  for( int i = 0; i < nodeSets.size(); ++i )
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    createSet( setName );
  }

  // Then loop over them in parallel.
  forAll< parallelHostPolicy >( nodeSets.size(), [&]( localIndex const i ) -> void
  {
    auto const & setWrapper = nodeSets[i];
    string const & setName = setWrapper->getName();
    SortedArrayView< localIndex const > const targetSet = nodeManager.sets().getReference< SortedArray< localIndex > >( setName ).toViewConst();
    constructSetFromSetAndMap( targetSet, m_toNodesRelation, setName );
  } );

  setDomainBoundaryObjects( faceManager );
}

void EdgeManager::buildEdges( localIndex const numNodes,
                              ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                              ArrayOfArrays< localIndex > & faceToEdgeMap )
{
  localIndex const numEdges = buildEdgeMaps( numNodes, faceToNodeMap,
                                             faceToEdgeMap,
                                             m_toFacesRelation,
                                             m_toNodesRelation );
  // FIXME I've no idea of what I am doing
  m_toNodesRelation.resize( numEdges );
//  m_edgesToFractureConnectorsEdges.resize( numEdges );
  m_fractureConnectorsEdgesToEdges.resize( numEdges );
  m_fractureConnectorEdgesToFaceElements.resize( numEdges );
}


void EdgeManager::setDomainBoundaryObjects( FaceManager const & faceManager )
{
  // get the "isDomainBoundary" field from the faceManager. This should have
  // been set already!
  arrayView1d< integer const > const & isFaceOnDomainBoundary = faceManager.getDomainBoundaryIndicator();

  // get the "isDomainBoundary" field from for *this, and set it to zero
  arrayView1d< integer > const & isEdgeOnDomainBoundary = this->getDomainBoundaryIndicator();
  isEdgeOnDomainBoundary.zero();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  // loop through all faces
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isFaceOnDomainBoundary[kf] == 1 )
    {
      localIndex const numFaceEdges = faceToEdgeMap.sizeOfArray( kf );

      // loop over all nodes connected to face, and set isNodeDomainBoundary
      for( localIndex a = 0; a < numFaceEdges; ++a )
      {
        isEdgeOnDomainBoundary( faceToEdgeMap( kf, a ) ) = 1;
      }
    }
  }
}

bool EdgeManager::hasNode( const localIndex edgeID, const localIndex nodeID ) const
{
  return m_toNodesRelation( edgeID, 0 ) == nodeID || m_toNodesRelation( edgeID, 1 ) == nodeID;
}

//localIndex EdgeManager::FindEdgeFromNodeIDs(const localIndex nodeA, const localIndex nodeB, const NodeManager *
// nodeManager)
//{
//  localIndex val = std::numeric_limits<localIndex>::max();
//
//  if (nodeA == nodeB)
//    return (val);
//
//  for( SortedArray<localIndex>::const_iterator iedge=nodeManager->m_nodeToEdgeMap[nodeA].begin() ;
//       iedge!=nodeManager->m_nodeToEdgeMap[nodeA].end() ; ++iedge )
//  {
//    if (hasNode(*iedge, nodeB))
//      val = *iedge;
//  }
//  return(val);
//}

void EdgeManager::setIsExternal( FaceManager const & faceManager )
{
  // get the "isExternal" field from the faceManager...This should have been
  // set already!
  arrayView1d< integer const > const & isExternalFace = faceManager.isExternal();

  ArrayOfArraysView< localIndex const > const & faceToEdges = faceManager.edgeList().toViewConst();

  // get the "isExternal" field from for *this, and set it to zero
  m_isExternal.zero();

  // loop through all faces
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    // check to see if the face is on a domain boundary
    if( isExternalFace[kf] == 1 )
    {
      // loop over all nodes connected to face, and set isNodeDomainBoundary
      localIndex const numEdges = faceToEdges.sizeOfArray( kf );
      for( localIndex a = 0; a < numEdges; ++a )
      {
        m_isExternal[ faceToEdges( kf, a ) ] = 1;
      }
    }
  }
}


void EdgeManager::extractMapFromObjectForAssignGlobalIndexNumbers( NodeManager const & nodeManager,
                                                                   std::vector< std::vector< globalIndex > > & globalEdgeNodes )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numEdges = size();

  arrayView2d< localIndex const > const edgeNodes = this->nodeList();
  arrayView1d< integer const > const isDomainBoundary = this->getDomainBoundaryIndicator();

  globalEdgeNodes.resize( numEdges );

  forAll< parallelHostPolicy >( numEdges, [&]( localIndex const edgeID )
  {
    std::vector< globalIndex > & curEdgeGlobalNodes = globalEdgeNodes[ edgeID ];

    if( isDomainBoundary( edgeID ) )
    {
      curEdgeGlobalNodes.resize( 2 );

      for( localIndex a = 0; a < 2; ++a )
      {
        curEdgeGlobalNodes[ a ]= nodeManager.localToGlobalMap()( edgeNodes[ edgeID ][ a ] );
      }

      std::sort( curEdgeGlobalNodes.begin(), curEdgeGlobalNodes.end() );
    }
  } );
}


void EdgeManager::connectivityFromGlobalToLocal( const SortedArray< localIndex > & indices,
                                                 const map< globalIndex, localIndex > & nodeGlobalToLocal,
                                                 const map< globalIndex, localIndex > & GEOSX_UNUSED_PARAM( faceGlobalToLocal ) )
{


  for( localIndex const ke : indices )
  {
    for( localIndex a=0; a<m_toNodesRelation.size( 1 ); ++a )
    {
      const globalIndex gnode = m_toNodesRelation( ke, a );
      const localIndex lnode = stlMapLookup( nodeGlobalToLocal, gnode );
      m_toNodesRelation( ke, a ) = lnode;
    }
  }

//  array1d<SortedArray<localIndex>>* const edgesToFlowFaces =
// GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
//  if( edgesToFlowFaces != NULL )
//  {
//    for( SortedArray<localIndex>::const_iterator ke=indices.begin() ; ke!=indices.end() ; ++ke )
//    {
//      SortedArray<localIndex>& edgeToFlowFaces = (*edgesToFlowFaces)[*ke];
//      SortedArray<localIndex> newSet;
//      for( SortedArray<localIndex>::iterator faceIndex=edgeToFlowFaces.begin() ;
// faceIndex!=edgeToFlowFaces.end() ; ++faceIndex )
//      {
//
//        std::map<globalIndex,localIndex>::const_iterator MapIter =
// faceGlobalToLocal.find( static_cast<globalIndex>(*faceIndex) );
//        if( MapIter!=faceGlobalToLocal.end()  )
//        {
////          const localIndex faceLocalIndex = stlMapLookup( faceGlobalToLocal,
// static_cast<globalIndex>(*faceIndex) );
//          newSet.insert( MapIter->second );
//
//        }
//      }
//      edgeToFlowFaces = newSet;
//    }
//
//  }

}

localIndex EdgeManager::packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return packUpDownMapsPrivate< false >( junk, packList );
}

localIndex EdgeManager::packUpDownMaps( buffer_unit_type * & buffer,
                                        arrayView1d< localIndex const > const & packList ) const
{
  return packUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex EdgeManager::packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                               arrayView1d< localIndex const > const & packList ) const
{
  arrayView1d< globalIndex const > const localToGlobal = localToGlobalMap();
  arrayView1d< globalIndex const > nodeLocalToGlobal = nodeList().relatedObjectLocalToGlobal();
  arrayView1d< globalIndex const > faceLocalToGlobal = faceList().relatedObjectLocalToGlobal();

  localIndex packedSize = bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::nodeListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toNodesRelation.base().toViewConst(),
                                           m_unmappedGlobalIndicesInToNodes,
                                           packList,
                                           localToGlobal,
                                           nodeLocalToGlobal );


  packedSize += bufferOps::Pack< DOPACK >( buffer, string( viewKeyStruct::faceListString() ) );
  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           m_toFacesRelation.base().toArrayOfArraysView(),
                                           m_unmappedGlobalIndicesInToFaces,
                                           packList,
                                           localToGlobal,
                                           faceLocalToGlobal );

  return packedSize;
}



localIndex EdgeManager::unpackUpDownMaps( buffer_unit_type const * & buffer,
                                          localIndex_array & packList,
                                          bool const overwriteUpMaps,
                                          bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  // GEOSX_MARK_FUNCTION;

  localIndex unPackedSize = 0;

  string nodeListString;
  unPackedSize += bufferOps::Unpack( buffer, nodeListString );
  GEOSX_ERROR_IF_NE( nodeListString, viewKeyStruct::nodeListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toNodesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToNodes,
                                     this->globalToLocalMap(),
                                     m_toNodesRelation.relatedObjectGlobalToLocal() );

  string faceListString;
  unPackedSize += bufferOps::Unpack( buffer, faceListString );
  GEOSX_ERROR_IF_NE( faceListString, viewKeyStruct::faceListString() );

  unPackedSize += bufferOps::Unpack( buffer,
                                     m_toFacesRelation,
                                     packList,
                                     m_unmappedGlobalIndicesInToFaces,
                                     this->globalToLocalMap(),
                                     m_toFacesRelation.relatedObjectGlobalToLocal(),
                                     overwriteUpMaps );

  return unPackedSize;
}

void EdgeManager::fixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::fixUpDownMaps( m_toNodesRelation,
                                    m_unmappedGlobalIndicesInToNodes,
                                    clearIfUnmapped );

  ObjectManagerBase::fixUpDownMaps( m_toFacesRelation.base(),
                                    m_toFacesRelation.relatedObjectGlobalToLocal(),
                                    m_unmappedGlobalIndicesInToFaces,
                                    clearIfUnmapped );
}

void EdgeManager::compressRelationMaps()
{
  m_toFacesRelation.compress();
}

void EdgeManager::depopulateUpMaps( std::set< localIndex > const & receivedEdges,
                                    ArrayOfArraysView< localIndex const > const & facesToEdges )
{
  ObjectManagerBase::cleanUpMap( receivedEdges, m_toFacesRelation.toView(), facesToEdges );
}


} /// namespace geosx
