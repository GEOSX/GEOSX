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

/**
 * @file SurfaceGenerator.cpp
 */

#include "SurfaceGenerator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "MPI_Communications/SpatialPartition.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"

#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

#ifdef USE_GEOSX_PTP
#include "GEOSX_PTP/ParallelTopologyChange.hpp"
#endif

#include <set>

namespace geosx
{

void ModifiedObjectLists::clearNewFromModified()
{
  for( auto const & a : newNodes )
  {
    modifiedNodes.erase(a);
  }

  for( auto const & a : newEdges )
  {
    modifiedEdges.erase(a);
  }

  for( auto const & a : newFaces )
  {
    modifiedFaces.erase(a);
  }

}

static localIndex GetOtherFaceEdge( const map< localIndex, std::pair<localIndex, localIndex> >& localFacesToEdges,
                                    const localIndex thisFace, const localIndex thisEdge )
{
  localIndex nextEdge = LOCALINDEX_MAX;

  const std::pair<localIndex, localIndex>& faceToEdges = stlMapLookup( localFacesToEdges, thisFace );
  if( faceToEdges.first == thisEdge )
  {
    nextEdge = faceToEdges.second;
  }
  else if( faceToEdges.second == thisEdge )
  {
    nextEdge = faceToEdges.first;
  }
  else
  {
    GEOS_ERROR( "SurfaceGenerator::Couldn't find thisEdge in localFacesToEdges[thisFace]" );
  }
  return nextEdge;
}

static void CheckForAndRemoveDeadEndPath( const localIndex edgeIndex,
                                          arrayView1d<integer> const & isEdgeExternal,
                                          map< localIndex, std::set<localIndex> >& edgesToRuptureReadyFaces,
                                          map< localIndex, std::pair<localIndex, localIndex> >& localVFacesToVEdges,
                                          std::set<localIndex>& nodeToRuptureReadyFaces )
{


  localIndex thisEdge = edgeIndex;

  // if the edge is internal and the edge is only attached to one ruptured face...
  while( isEdgeExternal[thisEdge]!=1 )
  {

    //    std::set<localIndex>& edgeToRuptureReadyFaces = stlMapLookup(edgesToRuptureReadyFaces,thisEdge);
    std::set<localIndex>& edgeToRuptureReadyFaces = edgesToRuptureReadyFaces[thisEdge];

    if( edgeToRuptureReadyFaces.size()!=1 )
      break;

    // then the index for the face that is a "dead end"
    localIndex deadEndFace = *(edgeToRuptureReadyFaces.begin());


    std::pair<localIndex, localIndex>& localVFaceToVEdges = stlMapLookup( localVFacesToVEdges, deadEndFace );

    // get the edge on the other side of the "dead end" face
    localIndex nextEdge = -1;
    if( localVFaceToVEdges.first == thisEdge )
      nextEdge = localVFaceToVEdges.second;
    else if( localVFaceToVEdges.second == thisEdge )
      nextEdge = localVFaceToVEdges.first;
    else
    {
      GEOS_ERROR( "SurfaceGenerator::FindFracturePlanes: Could not find the next edge when removing dead end faces." );
    }

    // delete the face from the working arrays
    edgeToRuptureReadyFaces.erase( deadEndFace );
    edgesToRuptureReadyFaces[nextEdge].erase( deadEndFace );
    nodeToRuptureReadyFaces.erase( deadEndFace );

    // if all the faces have been deleted, then go ahead and delete the top level entry
    if( edgeToRuptureReadyFaces.empty() )
      edgesToRuptureReadyFaces.erase( thisEdge );
    if( edgesToRuptureReadyFaces[nextEdge].empty() )
      edgesToRuptureReadyFaces.erase( nextEdge );

    // now increment the "thisEdge" to point to the other edge on the face that was just deleted
    thisEdge = nextEdge;
  }

}


SurfaceGenerator::SurfaceGenerator( const std::string& name,
                                    ManagedGroup * const parent ):
  SolverBase( name, parent ),
  m_failCriterion( 1 ),
  m_separableFaceSet()
{
  this->RegisterViewWrapper( viewKeyStruct::failCriterionString,
                             &this->m_failCriterion,
                             0 );

  this->RegisterViewWrapper( viewKeyStruct::fractureRegionNameString, &m_fractureRegionName, 0 )->
    setInputFlag(dataRepository::InputFlags::OPTIONAL)->
    setApplyDefaultValue("FractureRegion");


}

SurfaceGenerator::~SurfaceGenerator()
{
  // TODO Auto-generated destructor stub
}

void SurfaceGenerator::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0);

    NodeManager * const nodeManager = meshLevel->getNodeManager();
    EdgeManager * const edgeManager = meshLevel->getEdgeManager();
    FaceManager * const faceManager = meshLevel->getFaceManager();
    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    nodeManager->RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of node.");

    nodeManager->RegisterViewWrapper<integer_array>(viewKeyStruct::degreeFromCrackString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("connectivity distance from crack.");


    edgeManager->RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of the edge.");

    faceManager->RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::parentIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("Parent index of the face.");

    faceManager->RegisterViewWrapper<localIndex_array>(ObjectManagerBase::viewKeyStruct::childIndexString)->
      setApplyDefaultValue(-1)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_1)->
      setDescription("child index of the face.");

    faceManager->RegisterViewWrapper<integer_array>(viewKeyStruct::ruptureStateString)->
      setApplyDefaultValue(0)->
      setPlotLevel(dataRepository::PlotLevel::LEVEL_0)->
      setDescription("Rupture state of the face.0=not ready for rupture. 1=ready for rupture. 2=ruptured");

//    FaceElementRegion * const faceElementRegion = elemManager->GetRegion<FaceElementRegion>(m_fractureRegionName);
//    faceElementRegion->RegisterGroup<FaceElementSubRegion>(m_fractureRegionName+"0");
  }
}

void SurfaceGenerator::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>( dataRepository::keys::domain );
  for( auto & mesh : domain->group_cast<DomainPartition *>()->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>( mesh.second )->getMeshLevel( 0 );
    NodeManager * const nodeManager = meshLevel->getNodeManager();
    FaceManager * const faceManager = meshLevel->getFaceManager();

    arrayView1d<localIndex> & parentNodeIndex =
      nodeManager->getReference<localIndex_array>( nodeManager->viewKeys.parentIndex );

    arrayView1d<localIndex> & parentFaceIndex =
      faceManager->getReference<localIndex_array>( faceManager->viewKeys.parentIndex );

    arrayView1d<localIndex> & childFaceIndex =
      faceManager->getReference<localIndex_array>( faceManager->viewKeys.childIndex );

    parentNodeIndex = -1;
    parentFaceIndex = -1;
    childFaceIndex = -1;

    m_originalNodetoFaces = nodeManager->faceList();
    m_originalNodetoEdges = nodeManager->edgeList();
    m_originalFaceToEdges = faceManager->edgeList();

    nodeManager->RegisterViewWrapper( "usedFaces", &m_usedFacesForNode, 0 );
    m_usedFacesForNode.resize( nodeManager->size() );

    m_originalFacesToElemRegion = faceManager->elementRegionList();
    m_originalFacesToElemSubRegion = faceManager->elementSubRegionList();
    m_originalFacesToElemIndex = faceManager->elementList();
  }

}

real64 SurfaceGenerator::SolverStep( real64 const & time_n,
                                     real64 const & dt,
                                     const int cycleNumber,
                                     DomainPartition * const domain )
{
  int rval = 0;
  array1d<NeighborCommunicator> & neighbors = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );

  for( auto & mesh : domain->group_cast<DomainPartition *>()->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>( mesh.second )->getMeshLevel( 0 );

    {
      NodeManager * const nodeManager = meshLevel->getNodeManager();
      EdgeManager * const edgeManager = meshLevel->getEdgeManager();
      FaceManager * const faceManager = meshLevel->getFaceManager();
      ElementRegionManager * const elemManager = meshLevel->getElemManager();
      SpatialPartition & partition = domain->getReference<SpatialPartition,PartitionBase>(dataRepository::keys::partitionManager);


        rval = SeparationDriver( meshLevel,
                                 neighbors,
                                 partition.GetColor(),
                                 partition.NumColor(),
                                 0,
                                 time_n );

    }
  }

  NumericalMethodsManager * const
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>( dataRepository::keys::numericalMethodsManager );

  FiniteVolumeManager * const
  fvManager = numericalMethodManager->GetGroup<FiniteVolumeManager>( dataRepository::keys::finiteVolumeManager );

  for( auto & mesh : domain->group_cast<DomainPartition *>()->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>( mesh.second )->getMeshLevel( 0 );

    {
      NodeManager * const nodeManager = meshLevel->getNodeManager();
      EdgeManager * const edgeManager = meshLevel->getEdgeManager();
      FaceManager * const faceManager = meshLevel->getFaceManager();
      ElementRegionManager * const elemManager = meshLevel->getElemManager();
      FaceElementRegion * const fractureRegion = elemManager->GetRegion<FaceElementRegion>(this->m_fractureRegionName);

      for( localIndex a=0 ; a<fvManager->numSubGroups() ; ++a )
      {
        FluxApproximationBase * const fluxApprox = fvManager->GetGroup<FluxApproximationBase>(a);
        if( fluxApprox!=nullptr )
        {
          fluxApprox->addToFractureStencil( *domain,
                                            this->m_fractureRegionName );
          fractureRegion->m_recalculateConnectors.clear();
        }
      }
    }
  }


  return rval;
}



int SurfaceGenerator::SeparationDriver( MeshLevel * const mesh,
                                        array1d<NeighborCommunicator> & neighbors,
                                        int const tileColor,
                                        int const numTileColors,
                                        bool const prefrac,
                                        real64 const time )
{

  NodeManager & nodeManager = *(mesh->getNodeManager());
  EdgeManager & edgeManager = *(mesh->getEdgeManager());
  FaceManager & faceManager = *(mesh->getFaceManager());
  ElementRegionManager & elementManager = *(mesh->getElemManager());


  array1d<std::set<localIndex> > nodesToRupturedFaces;
  array1d<std::set<localIndex> > edgesToRupturedFaces;

  ArrayOfArraysView<localIndex> const & nodesToElementRegion = nodeManager.elementRegionList();
  ArrayOfArraysView<localIndex> const & nodesToElementSubRegion = nodeManager.elementSubRegionList();
  ArrayOfArraysView<localIndex> const & nodesToElementList = nodeManager.elementList();


  std::map<string, string_array > fieldNames;
  fieldNames["face"].push_back(viewKeyStruct::ruptureStateString);

  MPI_iCommData icomm;
  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, icomm );
  CommunicationTools::SynchronizePackSendRecv( fieldNames, mesh, neighbors, icomm );
  CommunicationTools::SynchronizeUnpack( mesh, neighbors, icomm );



  if( !prefrac )
  {

    if( m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {
      if( m_failCriterion == 1 )
      {
//        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet,
// std::numeric_limits<real64>::max()); //this->m_failstress );
      }
      else
      {
//        faceManager.UpdateRuptureStates( elementManager, nodeManager, m_separableFaceSet, m_failstress);
      }

//      real64_array& SIFonFace = faceManager.getReference<real64_array>("SIFonFace");
//      SIFonFace = std::numeric_limits<double>::min();

//      if (m_failCriterion == 1)
//      {
//        integer_array& ruptureState = faceManager.getReference<integer_array>(viewKeyStruct::ruptureStateString);
//        for (localIndex a=0; a<faceManager.size(); ++a)
//        {
//          if (faceManager.elementList()[a][1] == -1)
//          {
//            ruptureState[a]=0;
//          }
//        }
//      }
//
//
//      IdentifyRupturedFaces( nodeManager,
//                             edgeManager,
//                             faceManager,
//                             elementManager,
//                             partition,
//                             prefrac );

    }
    else
    {
//      UpdateRuptureStates( nodeManager,
//                           edgeManager,
//                           faceManager,
//                           elementManager,
//                           nodesToRupturedFaces,
//                           edgesToRupturedFaces,
//                           prefrac );

    }
  }
  else  // In the prefrac call, we need this to get the stressNOnFace, which will be used in the initialization of
        // contacts for preexisting fractures.
  {
//    for (localIndex kf = 0; kf < faceManager.size(); ++kf)
//      faceManager.CalculateStressOnFace(elementManager, nodeManager, kf);
  }


  if( prefrac )
  {
    ModifiedObjectLists modifiedObjects;
    CalculateKinkAngles( faceManager, edgeManager, nodeManager, modifiedObjects, prefrac );
  }

  // We do this here to get the nodesToRupturedFaces etc.
  // The fail stress check inside has been disabled
  PostUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces );

  int rval = 0;
  int rank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );

  //  array1d<MaterialBaseStateDataT*>&  temp = elementManager.m_ElementRegions["PM1"].m_materialStates;

  const arrayView1d<integer>& isNodeGhost = nodeManager.GhostRank();
//  const integer_array& isSeparable = nodeManager.getReference<integer_array>("isSeparable");
//  const integer_array& layersFromDomainBoundary = nodeManager.getReference<integer_array>("LayersFromDomainBoundary");
  arrayView1d<integer>& ruptureState = faceManager.getReference<integer_array>( "ruptureState" );
  const arrayView1d<integer>& isFaceGhost = faceManager.GhostRank();


  for( int color=0 ; color<numTileColors ; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( color==tileColor )
    {
      for( localIndex a=0 ; a<nodeManager.size() ; ++a )
      {
        int didSplit = 0;
        //      const localIndex nodeID = nodeManager.GetParentIndex(a);
        if( //layersFromDomainBoundary[a]>1 &&
          (   //isSeparable[a]
            true || prefrac)&&
          isNodeGhost[a]<0 &&
          nodesToElementList.sizeOfArray(a)>1 &&
          CheckNodeSplitability( a, nodeManager, faceManager, edgeManager, prefrac ) > 0 )  //&&
        //          nodesToRupturedFaces[a].size()>0 )
        {
          didSplit += ProcessNode( a, nodeManager, edgeManager, faceManager, elementManager, nodesToRupturedFaces, edgesToRupturedFaces, elementManager,
                                   modifiedObjects, prefrac );
          if( didSplit > 0 )
          {
            rval += didSplit;
            --a;
          }
        }
      }
    }
#ifdef USE_GEOSX_PTP

    modifiedObjects.clearNewFromModified();
    /// Nodes to edges in process node is not being set on rank 2. need to check that the new node->edge map is properly communicated
    ParallelTopologyChange::SyncronizeTopologyChange( mesh,
                                                      neighbors,
                                                      modifiedObjects);
#else
    AssignNewGlobalIndicesSerial( nodeManager, modifiedObjects.newNodes );
    AssignNewGlobalIndicesSerial( edgeManager, modifiedObjects.newEdges );
    AssignNewGlobalIndicesSerial( faceManager, modifiedObjects.newFaces );

#endif



    elementManager.forElementSubRegions<FaceElementSubRegion>( [&]( auto * const subRegion )
    {
      FaceElementSubRegion::NodeMapType & nodeMap = subRegion->nodeList();
      FaceElementSubRegion::FaceMapType & faceMap = subRegion->faceList();
      OrderedVariableOneToManyRelation const & facesToNodesMap = faceManager.nodeList();

      for( localIndex kfe=0 ; kfe<subRegion->size() ; ++kfe )
      {
        arrayView1d<localIndex const> const & faceToNodesMap0 = facesToNodesMap[faceMap[kfe][0]];
        arrayView1d<localIndex const> const & faceToNodesMap1 = facesToNodesMap[faceMap[kfe][1]];
        for( localIndex a=0 ; a<faceToNodesMap0.size() ; ++a )
        {
          localIndex const aa = a < 2 ? a : faceToNodesMap0.size() - a + 1;
          localIndex const bb = aa == 0 ? aa : faceToNodesMap0.size() - aa;

          // TODO HACK need to generalize to something other than quads
          nodeMap[kfe][a]   = faceToNodesMap0[aa];
          nodeMap[kfe][a+4] = faceToNodesMap1[bb];
        }
      }
    });
//    CalculateKinkAngles(faceManager, edgeManager, nodeManager, modifiedObjects, false);
//    MarkBirthTime(faceManager, modifiedObjects, time);
  }




  return rval;
}



bool SurfaceGenerator::ProcessNode( const localIndex nodeID,
                                    NodeManager & nodeManager,
                                    EdgeManager & edgeManager,
                                    FaceManager & faceManager,
                                    ElementRegionManager & elemManager,
                                    arrayView1d<std::set<localIndex> >& nodesToRupturedFaces,
                                    arrayView1d<std::set<localIndex> >& edgesToRupturedFaces,
                                    ElementRegionManager & elementManager,
                                    ModifiedObjectLists& modifiedObjects,
                                    const bool prefrac )
{
  bool didSplit = false;
  bool fracturePlaneFlag = true;

//  while( fracturePlaneFlag )
  {
    std::set<localIndex> facialRupturePath;
    map<localIndex, int> edgeLocations;
    map<localIndex, int> faceLocations;
    map< std::pair<CellElementSubRegion*, localIndex >, int> elemLocations;


    fracturePlaneFlag = FindFracturePlanes( nodeID,
                                            nodeManager,
                                            edgeManager,
                                            faceManager,
                                            elemManager,
                                            nodesToRupturedFaces,
                                            edgesToRupturedFaces,
                                            facialRupturePath,
                                            edgeLocations,
                                            faceLocations,
                                            elemLocations );
    if( fracturePlaneFlag )
    {
      MapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

      didSplit = true;
      PerformFracture( nodeID,
                       nodeManager,
                       edgeManager,
                       faceManager,
                       elementManager,
                       modifiedObjects,
                       nodesToRupturedFaces,
                       edgesToRupturedFaces,
                       facialRupturePath,
                       edgeLocations,
                       faceLocations,
                       elemLocations );
      MapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

    }
  }

  return didSplit;
}

bool SurfaceGenerator::FindFracturePlanes( const localIndex nodeID,
                                           const NodeManager & nodeManager,
                                           const EdgeManager & edgeManager,
                                           const FaceManager & faceManager,
                                           ElementRegionManager & elemManager,
                                           const arrayView1d<std::set<localIndex> >& nodesToRupturedFaces,
                                           const arrayView1d<std::set<localIndex> >& edgesToRupturedFaces,
                                           std::set<localIndex>& separationPathFaces,
                                           map<localIndex, int>& edgeLocations,
                                           map<localIndex, int>& faceLocations,
                                           map< std::pair< CellElementSubRegion*, localIndex >, int>& elemLocations )
{

  arrayView1d<localIndex> const & parentNodeIndices =
    nodeManager.getReference<array1d<localIndex>>( nodeManager.viewKeys.parentIndex );

  localIndex const parentNodeIndex = ObjectManagerBase::GetParentRecusive( parentNodeIndices, nodeID );

  arrayView1d<localIndex> const & parentFaceIndices =
    faceManager.getReference<array1d<localIndex>>( faceManager.viewKeys.parentIndex );

  arrayView1d<localIndex> const & childFaceIndices =
    faceManager.getReference<array1d<localIndex>>( faceManager.viewKeys.childIndex );

  std::set<localIndex> const & vNodeToRupturedFaces = nodesToRupturedFaces[parentNodeIndex];

  set<localIndex> const & originalNodeToEdges = m_originalNodetoEdges[parentNodeIndex];
  set<localIndex> const & originalNodeToFaces = m_originalNodetoFaces[parentNodeIndex];

  set<localIndex> const & nodeToEdges = nodeManager.edgeList()[nodeID];
  set<localIndex> const & nodeToFaces = nodeManager.faceList()[nodeID];

  arrayView1d< array1d < localIndex > > const & facesToEdges = faceManager.edgeList();


//  array1d< ReferenceWrapper<localIndex_array> > nodeToElements
//  const std::set< std::pair<CellBlockSubRegion*,localIndex> >&
//  nodesToElements = nodeManager.m_toElementsRelation[nodeID] ;

  arraySlice1d<localIndex const> const & nodeToElementRegion = nodeManager.elementRegionList()[nodeID];
  arraySlice1d<localIndex const> const & nodeToElementSubRegion = nodeManager.elementSubRegionList()[nodeID];
  arraySlice1d<localIndex const> const & nodeToElementIndex = nodeManager.elementList()[nodeID];

  // ***** BACKWARDS COMPATIBLITY HACK
  std::set< std::pair<CellElementSubRegion*, localIndex> > nodesToElements;


  for( localIndex k=0 ; k<nodeManager.elementRegionList().sizeOfArray(nodeID) ; ++k )
  {
    nodesToElements.insert( std::make_pair( elemManager.GetRegion( nodeToElementRegion[k] )->
                                            GetSubRegion<CellElementSubRegion>( nodeToElementSubRegion[k] ),
                                            nodeToElementIndex[k] ) );
  }


  // ***** END BACKWARDS COMPATIBLITY HACK


  arrayView1d<integer> const & isEdgeExternal = edgeManager.isExternal();

  arrayView1d<localIndex_array> const & originalFaceToEdges = m_originalFaceToEdges;

//  const std::set<localIndex>& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[nodeID];

  // **** local working arrays *****************************************************************************************

  // array to hold the faces ready for rupture. It is filled with the intersection of the virtual parent faces
  // associated
  // with all faces attached to the node, and all ruptured virtual faces attached to the virtual parent node.
  std::set<localIndex> nodeToRuptureReadyFaces;
  for( auto i : nodeToFaces )
  {
    const localIndex parentFaceIndex = ( parentFaceIndices[i] == -1 ) ? i : parentFaceIndices[i];

    if( vNodeToRupturedFaces.count( parentFaceIndex ) > 0 )
    {
      nodeToRuptureReadyFaces.insert( parentFaceIndex );
    }
  }


  // local map to hold the edgesToRuptureReadyFaces
  std::map< localIndex, std::set<localIndex> > edgesToRuptureReadyFaces;
  for( auto edgeIndex : originalNodeToEdges )
  {
    if( !(edgesToRupturedFaces[edgeIndex].empty()) )
      edgesToRuptureReadyFaces[edgeIndex].insert( edgesToRupturedFaces[edgeIndex].begin(), edgesToRupturedFaces[edgeIndex].end() );
  }


  // need a map from faces to edges that are attached to the node
  std::map< localIndex, std::pair<localIndex, localIndex> > nodeLocalFacesToEdges;
  for( auto kf : originalNodeToFaces )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( auto ke : originalFaceToEdges[kf] )
    {
      if( originalNodeToEdges.count( ke ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOS_ERROR( "SurfaceGenerator::FindFracturePlanes: invalid edge." );
    }


    nodeLocalFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // ***** remove dead end paths ***************************************************************************************
  // if the edge is not external, and the size of edgesToRupturedFaces is less than 2, then the edge is a dead-end
  // as far as a rupture plane is concerned. The face associated with the edge should be removed from the working
  // list of ruptured faces.

  // loop over all the edges
  for( auto edgeIndex : originalNodeToEdges )
  {

    CheckForAndRemoveDeadEndPath( edgeIndex,
                                  isEdgeExternal,
                                  edgesToRuptureReadyFaces,
                                  nodeLocalFacesToEdges,
                                  nodeToRuptureReadyFaces );

  }

  // if there are no ruptured faces attached to the node, then we are done.
  // or if there are no faces that have not been used in a rupture path for this node...we are done.
  if( nodeToRuptureReadyFaces.empty() )//|| nodeToRuptureReadyFaces.size() == usedFaces.size() )
  {
    return false;
  }

  // ***** find separation path ****************************************************************************************

  // ***** find starting face *****
  // We need to find a starting point for the path. The path must have a face that does has not been used in a previous
  // path for this node...otherwise it is the same path as used previously.
  localIndex startingEdge = INT_MAX;
  localIndex startingFace = INT_MAX;
  bool startingEdgeExternal = false;

  for( std::set<localIndex>::const_iterator i=nodeToRuptureReadyFaces.begin() ; i!=nodeToRuptureReadyFaces.end() ; ++i )
  {
    // check to see if this face has been used to split this node as part of a previously used path
    if( m_usedFacesForNode[nodeID].count( *i )==0 )
    {
      // great! It hasn't. It's on like Donkey Kong.
      startingFace = *i;

      if( isEdgeExternal[nodeLocalFacesToEdges[startingFace].first]==1 )
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].first;
        startingEdgeExternal = true;
        break;
      }
      else if( isEdgeExternal[nodeLocalFacesToEdges[startingFace].second]==1 )
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].second;
        startingEdgeExternal = true;
        break;
      }
      else
      {
        startingEdge = nodeLocalFacesToEdges[startingFace].first;
      }
    }
  }

  // if the starting face was not set, then we don't have a rupture surface....so just quit.
  if( startingFace==INT_MAX || startingEdge==INT_MAX )
  {
    return false;
    //    GEOS_ERROR("Fracturantor3::FindFracturePlanes: couldn't set starting face/edge");
  }



  // so now the working arrays have been purged of any faces that are on a dead-end path. All remaining faces
  // are part of a separation plane...of course, there can be more than one...which is bad. We will just take the first
  // path we find, and call this function again after the selected path is processed. Since the ruptureState of a face
  // is set to 2 after it is ruptured, if we enforce that candidate paths must have a face with a ruptureState of 1,
  // then
  // everything will work out. Also since the new nodes that are created will have higher node indices than the
  // current node, they will be checked for separation prior to completion of the separation driver.



  // We now have to define the separation plane over which a node/face/edge will be split, and all elements on one side
  // of the plane get one set of objects, and all elements on the other side get the other set.



  {
    // now we start the process of setting the separation path. Begin by
    localIndex thisEdge = startingEdge;
    localIndex thisFace = startingFace;

    localIndex nextEdge = INT_MAX;
    localIndex nextFace = INT_MAX;

    //localIndex lastEdge = INT_MAX;
    //localIndex lastFace = INT_MAX;

    // the seprationPath is used to hold combinations of edge and face
    std::map<localIndex, int> facesInPath;
    std::map<localIndex, int> edgesInPath;

    int numFacesInPath = 0;
    edgesInPath[thisEdge] = numFacesInPath;
    facesInPath[thisFace] = numFacesInPath++;

    localIndex_array facePath;
    localIndex_array edgePath;

    facePath.push_back( thisFace );
    edgePath.push_back( thisEdge );

    // now walk from face->edge->face->edge etc. until we get to an external edge, or back to the startingEdge.
    // the breakFlag indicates that we have found a complete separation path
    bool breakFlag = false;
    while( !breakFlag )
    {

      // get the next edge in the path...it is on the other side of "thisFace", so assign the other edge on the face as
      // the next edge

      nextEdge = GetOtherFaceEdge( nodeLocalFacesToEdges, thisFace, thisEdge );


      // if the nextEdge has already been used in the path, and the nextEdge is not the starting edge, then we have
      // to take a step back and try a different path
      if( edgesInPath.count( nextEdge )==1 && nextEdge!=startingEdge )
      {
        // first check to see if we can use the path without the preceding
        return false;
      }

      // if we have reached an external face, or the edge is already in the path, then we are done
      if( (isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) || edgesInPath.count( nextEdge )==1 )
      {
        // check to see if nextEdge is the startingEdge. If not, then all faces must that are before the nextEdge must
        // NOT be included in the path!!!
        if( nextEdge!=startingEdge && !(isEdgeExternal[nextEdge]==1 && startingEdgeExternal ) )
        {
          std::cout<<std::endl;


          std::cout<<"  NodeID, ParentID = "<<nodeID<<", "<<parentNodeIndex<<std::endl;
          std::cout<<"  Starting Edge/Face = "<<startingEdge<<", "<<startingFace<<std::endl;
          std::cout<<"  Face Separation Path = ";
          for( localIndex_array::const_iterator kf=facePath.begin() ; kf!=facePath.end() ; ++kf )
          {
            std::cout<<*kf<<", ";
          }
          std::cout<<std::endl;

          std::cout<<"  Edge Separation Path = ";
          for( localIndex_array::const_iterator kf=edgePath.begin() ; kf!=edgePath.end() ; ++kf )
          {
            std::cout<<*kf<<", ";
          }
          std::cout<<std::endl;


          GEOS_ERROR( "crap" );
        }

        // add faces in the path to separationPathFaces
        for( std::map<localIndex, int>::const_iterator kf=facesInPath.begin() ; kf!=facesInPath.end() ; ++kf )
        {
          separationPathFaces.insert( kf->first );
        }

        // break out of the while loop
        breakFlag = true;
      }
      else
      {
        // if the previous if statement is false, then what if we have reached an external edge, but the starting edge
        // was not external?? This means that we must continue the process from the edge opposite the startingEdge on
        // the
        // startingFace....which is hard-coded as the second entry in localFacesToEdges.
        if( isEdgeExternal[nextEdge]==1 )
        {
          nextEdge = nodeLocalFacesToEdges[startingFace].second;
        }

        // I sure hope that this is true!!
        if( edgesToRuptureReadyFaces[nextEdge].size() > 1 )
        {
          // we need to pick another face attached to the "next edge"
          // increment the face and edge, and add to the separationPathFaces


          {
            // OK...so we have an iterator that points to a candidate face. We prefer to move towards a face that is
            // ruptureState 1, so that we can get as much splitting done in this event. So we will loop over all the
            // faces attached to the edge, and pick one with ruptureState==1, otherwise just pick any one.
            bool pathFound = false;

            localIndex const parentFaceIndex = parentFaceIndices[thisFace] == -1 ? thisFace : parentFaceIndices[thisFace];
            std::pair<CellElementSubRegion*, localIndex>
            thisElem0 = std::make_pair( elemManager.GetRegion( m_originalFacesToElemRegion[thisFace][0] )->
                                        GetSubRegion<CellElementSubRegion>( m_originalFacesToElemSubRegion[thisFace][0] ),
                                        m_originalFacesToElemIndex[thisFace][0] );

            std::pair<CellElementSubRegion*, localIndex>
            thisElem1 = std::make_pair( elemManager.GetRegion( m_originalFacesToElemRegion[thisFace][1] )->
                                        GetSubRegion<CellElementSubRegion>( m_originalFacesToElemSubRegion[thisFace][1] ),
                                        m_originalFacesToElemIndex[thisFace][1] );

            // nextFaceQuality is intended to keep how desirable a face is for the rupture path.
            // A value of:
            //    0 -> the face is kind of creppy
            //    1 -> the face is does not turn a corner around the elements surrounding thisFace
            //    2 -> the face has not been used in a separation path
            //    3 -> a combination of 1 and 2.
            //    4 -> other edge on the face is the startingEdge.
            //
            int nextFaceQuality = -1;

            for( std::set<localIndex>::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin() ;
                 iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end() ; ++iter_edgeToFace )
            {
              if( *iter_edgeToFace != thisFace )
              {
                pathFound = true;



                const localIndex candidateFaceIndex = *iter_edgeToFace;
                int candidateFaceQuality = 0;


                localIndex candidateEdgeIndex = GetOtherFaceEdge( nodeLocalFacesToEdges, candidateFaceIndex, nextEdge );
                if( candidateEdgeIndex == startingEdge )
                {
                  nextFace = candidateFaceIndex;
                  break;
                }

                std::pair<CellElementSubRegion*, localIndex>
                nextElem0 = std::make_pair( elemManager.GetRegion( m_originalFacesToElemRegion[candidateFaceIndex][0] )->
                                            GetSubRegion<CellElementSubRegion>( m_originalFacesToElemSubRegion[candidateFaceIndex][0] ),
                                            m_originalFacesToElemIndex[candidateFaceIndex][0] );

                std::pair<CellElementSubRegion*, localIndex>
                nextElem1 = std::make_pair( elemManager.GetRegion( m_originalFacesToElemRegion[candidateFaceIndex][1] )->
                                            GetSubRegion<CellElementSubRegion>( m_originalFacesToElemSubRegion[candidateFaceIndex][1] ),
                                            m_originalFacesToElemIndex[candidateFaceIndex][1] );

                if( thisElem0 != nextElem0 && thisElem0 != nextElem1 &&
                    thisElem1 != nextElem0 && thisElem1 != nextElem1 )
                {
                  candidateFaceQuality += 1;
                }

                if( m_usedFacesForNode[nodeID].count( candidateFaceIndex ) == 0 )
                {
                  candidateFaceQuality += 2;
                }


                if( candidateFaceQuality > nextFaceQuality )
                {
                  nextFace = candidateFaceIndex;
                  nextFaceQuality = candidateFaceQuality;
                }

                if( candidateFaceQuality == 3 )
                {
                  break;
                }
              }
            }
            if( pathFound == false )
            {
              GEOS_ERROR( "SurfaceGenerator::FindFracturePlanes: couldn't find the next face in the rupture path" );
            }
          }

          //        lastEdge = thisEdge;
          //        lastFace = thisFace;

          thisEdge = nextEdge;
          thisFace = nextFace;
          //      separationPathFaces.insert( thisFace );
          edgesInPath[thisEdge] = numFacesInPath;
          facesInPath[thisFace] = numFacesInPath++;

          facePath.push_back( thisFace );
          edgePath.push_back( thisEdge );

        }
        else
        {
          GEOS_ERROR( "SurfaceGenerator::next edge in separation path is apparently  connected to less than 2 ruptured face" );
        }

      }
    }
  }


  //***** SET LOCATIONS ************************************************************************************************



  // need a map from faces to edges that are attached to the node
  std::map< localIndex, std::pair<localIndex, localIndex> > localFacesToEdges;
  for( auto kf : nodeToFaces )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( auto ke : facesToEdges[kf] )
    {
      if( edgeManager.hasNode( ke, nodeID ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOS_ERROR( "SurfaceGenerator::FindFracturePlanes: invalid edge." );
    }


    localFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // now we want to identify the objects on either side of the separation plane. First we assign an array to indicate
  // whether a face/edge is on the fracture plane.

  for( auto kf : nodeToFaces )
  {
    // iff the face is being split NOW, the set the faceLocation = -1.
    const localIndex virtualFaceIndex = ( parentFaceIndices[kf] == -1 ) ? kf : parentFaceIndices[kf];
    if( kf == virtualFaceIndex && childFaceIndices[kf] == -1 && separationPathFaces.count( kf ) )
    {
      faceLocations[kf] = -1;
    }
    else
    {
      faceLocations[kf] = INT_MIN;
    }

  }

  for( set<localIndex>::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    edgeLocations[*ke] = INT_MIN;
  }

  for( std::set< std::pair<CellElementSubRegion*, localIndex> >::const_iterator k=nodesToElements.begin() ; k!=nodesToElements.end() ; ++k )
  {
    elemLocations[*k] = INT_MIN;
  }



  /*
     SetLocations( 0, separationPathFaces, faceManager, nodesToElements, localFacesToEdges, //nodeToEdges,
                edgeLocations, faceLocations, elemLocations );

     if( !(SetLocations( 1, separationPathFaces, faceManager, nodesToElements, localFacesToEdges, //nodeToEdges,
                      edgeLocations, faceLocations, elemLocations )) )
     {
     return false;
     }*/

  SetLocations( separationPathFaces,
                elemManager,
                faceManager,
                nodesToElements,
                localFacesToEdges,
                edgeLocations,
                faceLocations,
                elemLocations );



  bool fail = false;

  for( set<localIndex>::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
  {
    if( edgeLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }
  for( set<localIndex>::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
  {
    if( faceLocations[*ke] == INT_MIN )
    {
      fail = true;
    }
  }
  /*
     std::cout<<"  NodeID, ParentID = "<<nodeID<<", "<<nodeID<<std::endl;
     std::cout<<"  separation path = ";
     for( std::set<localIndex>::const_iterator kf=separationPathFaces.begin() ; kf!=separationPathFaces.end() ; ++kf )
     {
      std::cout<<*kf<<", ";
     }
     std::cout<<std::endl;

     std::cout<<"  Starting Edge/Face = "<<startingEdge<<", "<<startingFace<<std::endl;
     for( std::set< std::pair<CellBlockSubRegion*,localIndex> >::const_iterator k=nodesToElements.begin() ;
        k!=nodesToElements.end() ; ++k )
     {
      std::cout<<"  elemLocations["<<k->second<<"] = "<<elemLocations[*k]<<std::endl;
     }

     for( std::set<localIndex>::const_iterator ke=nodeToFaces.begin() ; ke!=nodeToFaces.end() ; ++ke )
     {
      std::cout<<"  faceLocations["<<*ke<<"] = "<<faceLocations[*ke]<<std::endl;
     }

     for( std::set<localIndex>::const_iterator ke=nodeToEdges.begin() ; ke!=nodeToEdges.end() ; ++ke )
     {
      std::cout<<"  edgeLocations["<<*ke<<"] = "<<edgeLocations[*ke]<<std::endl;
     }
   */
  if( fail )
  {

    //    GEOS_ERROR("SurfaceGenerator::FindFracturePlanes: unset element,face, or edge");
    return false;
  }
  return true;
}

bool SurfaceGenerator::SetLocations( const std::set<localIndex>& separationPathFaces,
                                     ElementRegionManager & elemManager,
                                     const FaceManager& faceManager,
                                     const std::set< std::pair<CellElementSubRegion*, localIndex> >& nodesToElements,
                                     const map< localIndex, std::pair<localIndex, localIndex> >& localFacesToEdges,
                                     map<localIndex, int>& edgeLocations,
                                     map<localIndex, int>& faceLocations,
                                     map< std::pair< CellElementSubRegion*, localIndex >, int>& elemLocations )
{
  bool rval = true;
  //  const localIndex separationFace = *(separationPathFaces.begin());

  // insert an element attached to the separation face
  //  std::pair<CellBlockSubRegion*,localIndex> elem0 = m_virtualFaces.m_FaceToElementMap[separationFace][0] ;

  std::pair<CellElementSubRegion*, localIndex> elem0 = *(nodesToElements.begin());


  SetElemLocations( 0,
                    elem0,
                    separationPathFaces,
                    elemManager,
                    faceManager,
                    nodesToElements,
                    localFacesToEdges,
                    edgeLocations,
                    faceLocations,
                    elemLocations );

  return rval;
}



bool SurfaceGenerator::SetElemLocations( const int location,
                                         const std::pair<CellElementSubRegion*, localIndex >& k,
                                         const std::set<localIndex>& separationPathFaces,
                                         ElementRegionManager & elemManager,
                                         const FaceManager & faceManager,
                                         const std::set< std::pair<CellElementSubRegion*, localIndex> >& nodesToElements,
                                         const map< localIndex, std::pair<localIndex, localIndex> >& localFacesToEdges,
                                         map<localIndex, int>& edgeLocations,
                                         map<localIndex, int>& faceLocations,
                                         map< std::pair<CellElementSubRegion*, localIndex >, int>& elemLocations )
{

  arrayView1d<localIndex> const & parentFaceIndices =
    faceManager.getReference<localIndex_array>( faceManager.viewKeys.parentIndex );

  const int otherlocation = (location==0) ? 1 : 0;

  elemLocations[k] = location;


  // loop over all faces on the element
  for( localIndex kf=0 ; kf<k.first->faceList().size( 1 ) ; ++kf )
  {

    // define the actual face index, and the virtual face index
    const localIndex faceIndex = k.first->faceList()( k.second, kf );
    const localIndex virtualFaceIndex = ( parentFaceIndices[faceIndex] == -1 ) ?
                                        faceIndex : parentFaceIndices[faceIndex];

    // see if we can find the face in the faceLocations array.
    std::map<localIndex, int>::iterator iterFace = faceLocations.find( faceIndex );
    // if we can find the face in the faceLocations array, then we must process the face, otherwise it is not
    // connected to the node, so we do nothing.
    if( iterFace != faceLocations.end() )
    {

      if( faceLocations[faceIndex]==otherlocation )
        faceLocations[faceIndex] = -1;
      else if( faceLocations[faceIndex] == INT_MIN )
        faceLocations[faceIndex] = location;

      std::map< localIndex, std::pair<localIndex, localIndex> >::const_iterator iterF2E = localFacesToEdges.find( faceIndex );

      if( iterF2E != localFacesToEdges.end() )
      {
        const localIndex edge0 = (iterF2E->second).first;
        const localIndex edge1 = (iterF2E->second).second;

        if( edgeLocations[edge0]==otherlocation )
          edgeLocations[edge0] = -1;
        else if( edgeLocations[edge0] == INT_MIN )
          edgeLocations[edge0] = location;

        if( edgeLocations[edge1]==otherlocation )
          edgeLocations[edge1] = -1;
        else if( edgeLocations[edge1] == INT_MIN )
          edgeLocations[edge1] = location;

      }



      // now we add the element that is a neighbor to the face
      // of course, this only happens if there are more than one element
      // attached to the face.
      if( m_originalFacesToElemIndex[virtualFaceIndex][1] != -1 )
      {
        localIndex const er0 = m_originalFacesToElemRegion[virtualFaceIndex][0];
        localIndex const er1 = m_originalFacesToElemRegion[virtualFaceIndex][1];

        localIndex const esr0 = m_originalFacesToElemSubRegion[virtualFaceIndex][0];
        localIndex const esr1 = m_originalFacesToElemSubRegion[virtualFaceIndex][1];


        const std::pair<CellElementSubRegion*, localIndex>
        elemIndex0 = { elemManager.GetRegion( er0 )->GetSubRegion<CellElementSubRegion>( esr0 ),
                       m_originalFacesToElemIndex[virtualFaceIndex][0] };

        const std::pair<CellElementSubRegion*, localIndex>
        elemIndex1 = { elemManager.GetRegion( er1 )->GetSubRegion<CellElementSubRegion>( esr1 ),
                       m_originalFacesToElemIndex[virtualFaceIndex][1] };

        const std::pair<CellElementSubRegion*, localIndex>& nextElem = ( elemIndex0 == k ) ? elemIndex1 : elemIndex0;
        const int nextLocation = (separationPathFaces.count( virtualFaceIndex )==0) ? location : otherlocation;

        // if the first element is the one we are on, and the element is attached
        // to the splitting node, then add the second element to the list.
        if( nodesToElements.find( nextElem )!=nodesToElements.end() )
        {
          if( elemLocations[nextElem]==INT_MIN )
          {
            SetElemLocations( nextLocation,
                              nextElem,
                              separationPathFaces,
                              elemManager,
                              faceManager,
                              nodesToElements,
                              localFacesToEdges,
                              edgeLocations,
                              faceLocations,
                              elemLocations );
          }
        }
      }
    }
  }

  return true;
}



void SurfaceGenerator::PerformFracture( const localIndex nodeID,
                                        NodeManager & nodeManager,
                                        EdgeManager & edgeManager,
                                        FaceManager & faceManager,
                                        ElementRegionManager & elementManager,
                                        ModifiedObjectLists& modifiedObjects,
                                        arrayView1d<std::set<localIndex> >& nodesToRupturedFaces,
                                        arrayView1d<std::set<localIndex> >& edgesToRupturedFaces,
                                        const std::set<localIndex>& separationPathFaces,
                                        const map<localIndex, int>& edgeLocations,
                                        const map<localIndex, int>& faceLocations,
                                        const map< std::pair<CellElementSubRegion*, localIndex >, int>& elemLocations )
{

  int rank=0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );

  arrayView1d<R1Tensor> const & X = nodeManager.referencePosition();
  arrayView1d<set<localIndex> > & nodesToEdges = nodeManager.edgeList();
  arrayView1d<set<localIndex> > & nodesToFaces = nodeManager.faceList();
  ArrayOfArrays<localIndex> & nodesToElementRegions = nodeManager.elementRegionList();
  ArrayOfArrays<localIndex> & nodesToElementSubRegions = nodeManager.elementSubRegionList();
  ArrayOfArrays<localIndex> & nodesToElementIndex = nodeManager.elementList();


  arrayView2d<localIndex> & edgesToNodes = edgeManager.nodeList();
  arrayView1d<set<localIndex> > & edgesToFaces = edgeManager.faceList();
//  arrayView1d<set<localIndex>> & edgesToElementRegions = edgeManager.elementRegionList();
//  arrayView1d<set<localIndex>> & edgesToElementSubRegions = edgeManager.elementSubRegionList();
//  arrayView1d<set<localIndex>> & edgesToElementIndex = edgeManager.elementList();

  OrderedVariableOneToManyRelation & facesToNodes = faceManager.nodeList();
  OrderedVariableOneToManyRelation & facesToEdges = faceManager.edgeList();
  arrayView2d<localIndex> & facesToElementRegions = faceManager.elementRegionList();
  arrayView2d<localIndex> & facesToElementSubRegions = faceManager.elementSubRegionList();
  arrayView2d<localIndex> & facesToElementIndex = faceManager.elementList();

  FaceElementRegion * const fractureElementRegion = elementManager.GetRegion<FaceElementRegion>("Fracture");

  arrayView1d<R1Tensor> const & faceNormals = faceManager.faceNormal();

  arrayView1d<localIndex> const & parentFaceIndices =
    faceManager.getReference<localIndex_array>( faceManager.viewKeys.parentIndex );
  arrayView1d<localIndex const> const & parentEdgeIndices =
    edgeManager.getReference<localIndex_array>( ObjectManagerBase::viewKeyStruct::parentIndexString );



//  integer_array* flowEdgeType = edgeManager.getReferencePointer<int>("flowEdgeType");
//  integer_array* flowFaceType = faceManager.getReferencePointer<int>("flowFaceType");
//  real64_array & stressNOnFace = faceManager.getReference<real64>("stressNOnFace");
//  arrayView1d<R1Tensor>& stressTOnFace = faceManager.getReference<R1Tensor>("stressTOnFace");
//  arrayView1d<set<localIndex>>* edgeToFlowFaces = edgeManager.GetUnorderedVariableOneToManyMapPointer("edgeToFlowFaces");
//arrayView1d<R1Tensor>* faceNormal = faceManager.getReferencePointer<R1Tensor>("faceNormal0");

  arrayView1d<integer> & degreeFromCrack =
    nodeManager.getReference<integer_array>( viewKeyStruct::degreeFromCrackString );

  //  const arrayView1d<localIndex_array>& childEdgeIndex = edgeManager.GetVariableOneToManyMap( "childIndices" );


  //  const localIndex nodeID = nodeManager.GetParentIndex( nodeID );



//  real64_array* delta0N = faceManager.getReferencePointer<real64>("delta0N");
//  real64_array* faceContactStiffness = faceManager.getReferencePointer<real64>("faceContactStiffness");
//  array1d<R1Tensor>* stressShear0 = faceManager.getReferencePointer<R1Tensor>("stressShear0");

  // ***** split all the objects first *****

  // Split the node into two, using the original index, and a new one.
  localIndex newNodeIndex;
  if( verboseLevel() )
  {
    GEOS_LOG_RANK("");
    std::cout<<"Splitting node "<<nodeID<<" along separation plane faces: ";
    for( std::set<localIndex>::const_iterator i=separationPathFaces.begin() ; i!=separationPathFaces.end() ; ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;
  }


  nodeManager.SplitObject( nodeID, rank, newNodeIndex );

  modifiedObjects.newNodes.insert( newNodeIndex );
  modifiedObjects.modifiedNodes.insert( nodeID );

  nodesToElementRegions.clearArray(newNodeIndex);
  nodesToElementSubRegions.clearArray(newNodeIndex);
  nodesToElementIndex.clearArray(newNodeIndex);

  nodesToEdges[newNodeIndex].clear();
  nodesToFaces[newNodeIndex].clear();

  degreeFromCrack[nodeID] = 0;
  degreeFromCrack[newNodeIndex] = 0;

  //TODO HACK...should recalculate mass
//  const real64 newMass = 0.5 * (*nodeManager.m_mass)[nodeID];
//  (*nodeManager.m_mass)[nodeID] = newMass;
//  (*nodeManager.m_mass)[newNodeIndex] = newMass;

  //TODO Either change m_usedFacesForNode to array<std::set> or add insert with iterator to SortedArray
  for (auto const val : separationPathFaces)
  {
    m_usedFacesForNode[nodeID].insert( val );
    m_usedFacesForNode[newNodeIndex].insert( val );
  }

//  set<localIndex>& usedFacesNew = nodeManager.getReference< array1d<set<localIndex>> >("usedFaces")[newNodeIndex];
//  usedFacesNew = usedFaces[nodeID];


  if( verboseLevel() )
    std::cout<<"Done splitting node "<<nodeID<<" into nodes "<<nodeID<<" and "<<newNodeIndex<<std::endl;

  // split edges
  map<localIndex, localIndex> splitEdges;
  // loop over all edges connected to the node
  for( map<localIndex, int>::const_iterator iter_edge=edgeLocations.begin() ; iter_edge!=edgeLocations.end() ; ++iter_edge )
  {
    const localIndex& parentEdgeIndex = iter_edge->first;
    const int& location = iter_edge->second;

    // if the edge is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newEdgeIndex;

      edgeManager.SplitObject( parentEdgeIndex, rank, newEdgeIndex );

      edgesToFaces[newEdgeIndex].clear();

      if( verboseLevel() )
      {
        GEOS_LOG_RANK("");
        std::cout<<"  Split edge "<<parentEdgeIndex<<" into edges "<<parentEdgeIndex<<" and "<<newEdgeIndex<<std::endl;
      }

      splitEdges[parentEdgeIndex] = newEdgeIndex;
      modifiedObjects.newEdges.insert( newEdgeIndex );
      modifiedObjects.modifiedEdges.insert( parentEdgeIndex );


      for( int a=0 ; a<2 ; ++a )
      {
        edgeManager.nodeList( newEdgeIndex, a ) = edgeManager.nodeList( parentEdgeIndex, a );
      }

    } //    if( location == -1  )
  } // for( map<localIndex,int>::const_iterator iter_edge...



  // split the faces
  arrayView1d<integer> & ruptureState = faceManager.getReference<integer_array>( "ruptureState" );
  map<localIndex, localIndex> splitFaces;


  set<localIndex> & externalFaces = faceManager.externalSet();

  // loop over all faces attached to the nodeID
  for( map<localIndex, int>::const_iterator iter_face=faceLocations.begin() ; iter_face!=faceLocations.end() ; ++iter_face )
  {
    const localIndex faceIndex = iter_face->first;
//    localIndex const parentFaceIndex = parentFaceIndices[faceIndex]==faceIndex ? faceIndex :
// parentFaceIndices[faceIndex];
    const int location = iter_face->second;
    // if the face is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newFaceIndex;

      if( faceManager.SplitObject( faceIndex, rank, newFaceIndex ) )
      {

        if( verboseLevel() )
        {
          GEOS_LOG_RANK("");
          std::cout<<"  Split face "<<faceIndex<<" into faces "<<faceIndex<<" and "<<newFaceIndex<<std::endl;
        }

        splitFaces[faceIndex] = newFaceIndex;
        modifiedObjects.newFaces.insert( newFaceIndex );
        modifiedObjects.modifiedFaces.insert( faceIndex );

        ruptureState[faceIndex] = 2;
        ruptureState[newFaceIndex] = 2;

        facesToEdges[newFaceIndex] = facesToEdges[faceIndex];
//        facesToNodes[newFaceIndex] = facesToNodes[faceIndex];
        for( localIndex a=0 ; a<facesToNodes[faceIndex].size() ; ++a )
        {
          localIndex const aa = a == 0 ? a : facesToNodes[faceIndex].size() - a;
          facesToNodes[newFaceIndex][aa] = facesToNodes[faceIndex][a];
        }
        faceNormals[newFaceIndex] *= -1;

        externalFaces.insert( newFaceIndex );
        externalFaces.insert( faceIndex );


        // Fu: All edges of the parent face should be external now.
        // We have to do the following because isExternal attribute of the tip edge is not handled by the splitter.
        for( auto edgeIndex : faceManager.edgeList()[faceIndex] )
        {
          if( edgeManager.m_isExternal[edgeIndex]==0 )
          {
            edgeManager.m_isExternal[edgeIndex] = 2;
          }
        }
        for( auto nodeIndex : faceManager.nodeList()[faceIndex] )
        {
          if( nodeManager.m_isExternal[nodeIndex] )
          {
            nodeManager.m_isExternal[nodeIndex] = 2;
          }
        }

        {
          localIndex faceIndices[2] = {faceIndex,newFaceIndex};
          modifiedObjects.newElements[ {fractureElementRegion->getIndexInParent(),0} ].
          insert( fractureElementRegion->AddToFractureMesh( &faceManager,
                                                            this->m_originalFaceToEdges,
                                                            "default",
                                                            faceIndices ) );
        }
//        externalFaceManager.SplitFace(parentFaceIndex, newFaceIndex, nodeManager);

      } // if( faceManager.SplitObject( faceIndex, newFaceIndex ) )
    } // if( location == -1 )
  } // for( map<localIndex,int>::const_iterator iter_face


  // ***** now correct all the relations between the objects *****

  /* To accomplish this annoying yet exceedingly important task, we will take a "top down"
   * approach. Note that this is a two way correction, i.e. if we are correcting
   * elementToNodes, we also correct nodesToElements. This is summarized as:
   * 1) Loop over elements attached to the split node.
   *     2a) correct all relations between the single  element and the nodes.
   *     2b) Loop over all faces on the element
   *         3a) For each face, correct the face relations with the element
   *         3b) For each face, correct the face relations with the nodes
   *         3c) Loop over all edges on the face
   *             4a) For each edge, correct the face relations
   *             4b) for each edge, correct the node relations
   *
   *  The element location will define which side of the rupture everything
   *  is on.
   *  - location 0 gets the original node,edge,face.
   *  - location 1 gets the new node,edge,face.
   */

  arrayView1d<localIndex> const & parentFaceIndex =
    faceManager.getReference<localIndex_array>( faceManager.viewKeys.parentIndex );

  arrayView1d<localIndex> const & childFaceIndex =
    faceManager.getReference<localIndex_array>( faceManager.viewKeys.childIndex );



  // 1) loop over all elements attached to the nodeID
  for( map<std::pair<CellElementSubRegion*, localIndex>, int>::const_iterator iter_elem =
         elemLocations.begin() ; iter_elem != elemLocations.end() ; ++iter_elem )
  {
    const int& location = iter_elem->second;

    if( location==1 )
    {
      const std::pair<CellElementSubRegion*, localIndex >& elem = iter_elem->first;

      CellElementSubRegion& elemSubRegion = *(elem.first);
      ElementRegion * const elemRegion = elemSubRegion.getParent()->getParent()->group_cast<ElementRegion*>();
      string const elemRegionName = elemRegion->getName();

      localIndex const regionIndex = elementManager.GetRegions().getIndex( elemRegionName );
      localIndex const subRegionIndex = elemRegion->GetSubRegions().getIndex( elemSubRegion.getName() );
      const localIndex elemIndex = elem.second;

      modifiedObjects.modifiedElements[{regionIndex,subRegionIndex}].insert( elemIndex );


      arrayView2d<localIndex> & elemsToNodes = elemSubRegion.nodeList();
      arrayView2d<localIndex> & elemsToFaces = elemSubRegion.faceList();

      if( verboseLevel() > 1 )
        std::cout<<"Element "<<elemIndex<<std::endl;

      // 2a) correct elementToNode and nodeToElement
      if( verboseLevel() > 1 )
        std::cout<<"  Looping over all nodes on element, and correcting node<->element maps:"<<std::endl;


      R1Tensor elemCenter = {0.0, 0.0, 0.0};
      {
        // loop over all nodes on element
        if( verboseLevel() > 1 )
          std::cout<<"    m_ElementToNodeMap = ( ";
        for( localIndex a=0 ; a<elemsToNodes.size( 1 ) ; ++a )
        {
          elemCenter += X[ elemsToNodes[elemIndex][a] ];
          // if the node was just split
          if( elemsToNodes[elemIndex][a] == nodeID )
          {

            if( verboseLevel() > 1 )
              std::cout<<elemsToNodes[elemIndex][a]<<"->"<<newNodeIndex<<", ";

            elemsToNodes[elemIndex][a] = newNodeIndex;

            insert( nodeManager.toElementRelation(), newNodeIndex, regionIndex, subRegionIndex, elemIndex );
            erase( nodeManager.toElementRelation(), nodeID, regionIndex, subRegionIndex, elemIndex );
          }
          else if( verboseLevel() > 1 )
            std::cout<<elemsToNodes[elemIndex][a]<<", ";
        }
        elemCenter /= elemsToNodes.size( 1 );
        if( verboseLevel() > 1 )
          std::cout<<")"<<std::endl;

        if( verboseLevel() > 1 )
        {
          for( localIndex a=0 ; a<elemsToNodes.size( 1 ) ; ++a )
          {
            if( verboseLevel() > 1 )
            {
              std::cout<<"    nodeToElemMaps["<<elemsToNodes[elemIndex][a]<<"] = ( ";
              for( localIndex k=0 ; k<nodesToElementRegions.sizeOfArray(elemsToNodes[elemIndex][a]) ; ++k )
              {
                std::cout<<"["<<nodesToElementRegions[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodesToElementSubRegions[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodesToElementIndex[elemsToNodes[elemIndex][a]][k]<<"] , ";
              }
              std::cout<<" )"<<std::endl;
            }
          }

          if( verboseLevel() > 1 )
          {
            std::cout<<"    nodeToElemMaps["<<nodeID<<"] = ( ";
            for( localIndex k=0 ; k<nodesToElementRegions.sizeOfArray(nodeID) ; ++k )
            {
              std::cout<<"["<<nodesToElementRegions[nodeID][k]<<","
                       <<nodesToElementSubRegions[nodeID][k]<<","
                       <<nodesToElementIndex[nodeID][k]<<"] , ";
            }
            std::cout<<" )"<<std::endl;
          }
        }
      }



      // 2b) loop over all faces on element.
      if( verboseLevel() > 1 )
      {
        std::cout<<"  Looping over all faces on element (parent and child):"<<std::endl;
      }

      // we need to build a list of faces that is elemToFaces FOLLOWED by any
      // parent face of those indicated in elemToFaces

      // Now we do a loop over the facelist and process all the faces
      for( int kf=0 ; kf<elemSubRegion.numFacesPerElement() ; ++kf )
      {

        // set both faceID and newFaceID to the parent face.
        localIndex const faceIndex = elemsToFaces[elemIndex][kf];
        //        map<localIndex,localIndex>::iterator iterSplitFace = splitFaces.find(faceIndex);
        bool const isNewFace = (splitFaces.count( faceIndex )>0) ? true : false;
        localIndex const newFaceIndex = isNewFace ? childFaceIndex[faceIndex] : faceIndex;


        // 3a) check to see if the face was split. If so, then we will need
        // to alter the face relation with the elements in both directions.
        if( isNewFace )
        {
          // replace the parent face with the child face in elementToFace. Now
          // faceID is the parent face, and newFaceID is the child face.
          elemsToFaces[elemIndex][kf] = childFaceIndex[faceIndex];



          // add the element to the child faceToElem
//          faceManager.m_toElements[newFaceIndex].push_back( elem );

          facesToElementRegions[newFaceIndex][0] = regionIndex;
          facesToElementSubRegions[newFaceIndex][0] = subRegionIndex;
          facesToElementIndex[newFaceIndex][0] = elemIndex;
          facesToElementRegions[newFaceIndex][1] = -1;
          facesToElementSubRegions[newFaceIndex][1] = -1;
          facesToElementIndex[newFaceIndex][1] = -1;

          // remove the element from the parent face
          if( facesToElementRegions[faceIndex][0] == regionIndex &&
              facesToElementSubRegions[faceIndex][0] == subRegionIndex &&
              facesToElementIndex[faceIndex][0] == elemIndex )
          {
            facesToElementRegions[faceIndex][0] = facesToElementRegions[faceIndex][1];
            facesToElementSubRegions[faceIndex][0] = facesToElementSubRegions[faceIndex][1];
            facesToElementIndex[faceIndex][0] = facesToElementIndex[faceIndex][1];
            facesToElementRegions[faceIndex][1] = -1;
            facesToElementSubRegions[faceIndex][1] = -1;
            facesToElementIndex[faceIndex][1] = -1;
          }
          else if( facesToElementRegions[faceIndex][1] == regionIndex &&
                   facesToElementSubRegions[faceIndex][1] == subRegionIndex &&
                   facesToElementIndex[faceIndex][1] == elemIndex )
          {
            facesToElementRegions[faceIndex][1] = -1;
            facesToElementSubRegions[faceIndex][1] = -1;
            facesToElementIndex[faceIndex][1] = -1;
          }

          if( verboseLevel() > 1 )
          {
            std::cout<<"    facesToElementRegions["<<newFaceIndex<<"][0]    = "<<facesToElementRegions[newFaceIndex][0]<<std::endl;
            std::cout<<"    facesToElementSubRegions["<<newFaceIndex<<"][0] = "<<facesToElementSubRegions[newFaceIndex][0]<<std::endl;
            std::cout<<"    facesToElementIndex["<<newFaceIndex<<"][0]      = "<<facesToElementIndex[newFaceIndex][0]<<std::endl;
            std::cout<<"    facesToElementRegions["<<newFaceIndex<<"][1]    = "<<facesToElementRegions[newFaceIndex][1]<<std::endl;
            std::cout<<"    facesToElementSubRegions["<<newFaceIndex<<"][1] = "<<facesToElementSubRegions[newFaceIndex][1]<<std::endl;
            std::cout<<"    facesToElementIndex["<<newFaceIndex<<"][1]      = "<<facesToElementIndex[newFaceIndex][1]<<std::endl;

            std::cout<<"    facesToElementRegions["<<faceIndex<<"][0]    = "<<facesToElementRegions[faceIndex][0]<<std::endl;
            std::cout<<"    facesToElementSubRegions["<<faceIndex<<"][0] = "<<facesToElementSubRegions[faceIndex][0]<<std::endl;
            std::cout<<"    facesToElementIndex["<<faceIndex<<"][0]      = "<<facesToElementIndex[faceIndex][0]<<std::endl;
            std::cout<<"    facesToElementRegions["<<faceIndex<<"][1]    = "<<facesToElementRegions[faceIndex][1]<<std::endl;
            std::cout<<"    facesToElementSubRegions["<<faceIndex<<"][1] = "<<facesToElementSubRegions[faceIndex][1]<<std::endl;
            std::cout<<"    facesToElementIndex["<<faceIndex<<"][1]      = "<<facesToElementIndex[faceIndex][1]<<std::endl;

          }

          faceManager.SortFaceNodes( X, elemCenter, facesToNodes[newFaceIndex], facesToNodes[newFaceIndex].size() );


        } // if( splitFaces.count( faceID ) > 0 )

        modifiedObjects.modifiedFaces.insert( faceIndex );



        // 3b) correct faceToNodes and nodeToFaces

        if( verboseLevel() > 1 )
        {
          localIndex const parentFace = parentFaceIndex[newFaceIndex];
          if( parentFace!=-1 )
          {
            std::cout<<"    m_FaceToNodeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToNodeMap["<<newFaceIndex<<"] = ( ";
          }
        }

        // loop over all nodes on the face.
        for( auto & nodeIndex : faceManager.nodeList()[newFaceIndex] )
        {
          if( verboseLevel() > 1 )
            std::cout<<nodeIndex;

          // if the facenode is the one that is being split
          if( nodeIndex == nodeID )
          {
            nodeIndex = newNodeIndex;

            // if it is not a new face.
            if( !isNewFace )
            {
              // remove the face from the nodeToFaceMap of the parent node.
              nodesToFaces[nodeID].erase( faceIndex );

              // add the face to the nodeToFaceMap of the new node.
              nodesToFaces[nodeIndex].insert( faceIndex );

            }
            else
            {
              // it is a new face

              // insert the newFace into the nodeToFaceMap of the newNode
              nodesToFaces[nodeIndex].insert( newFaceIndex );
            }
            if( verboseLevel() > 1 )
              std::cout<<"->"<<nodeIndex<<", ";
          }
          else // the node is not being split
          {
            nodesToFaces[nodeIndex].insert( newFaceIndex );

            if( verboseLevel() > 1 )
              std::cout<<", ";
          }

        }
        if( verboseLevel() > 1 )
          std::cout<<")"<<std::endl;



        // faceToEdges
        if( verboseLevel() > 1 )
        {
          const localIndex parentFace = parentFaceIndex[newFaceIndex];
          if( parentFace!=-1 )
          {
            std::cout<<"    m_FaceToEdgeMap["<<parentFace<<"->"<<newFaceIndex<<"] = ( ";
          }
          else
          {
            std::cout<<"    m_FaceToEdgeMap["<<newFaceIndex<<"] = ( ";
          }
        }
        // loop over all edges on face
        for( auto & edgeIndex : facesToEdges[newFaceIndex] )
        {

          // if the edge was just split
          if( splitEdges.count( edgeIndex ) > 0 )
          {
            if( faceIndex == newFaceIndex )
            {
              edgesToFaces[edgeIndex].erase( faceIndex );
            }

            edgeIndex = splitEdges[edgeIndex];
          }
          edgesToFaces[edgeIndex].insert( newFaceIndex );

          modifiedObjects.modifiedEdges.insert( edgeIndex );

          if( verboseLevel() > 1 )
            std::cout<<edgeIndex;



          //edgesToNodes
          if( verboseLevel() > 1 )
          {
            std::cout<<"(";
          }

          {
            for( localIndex a=0 ; a<edgesToNodes.size( 1 ) ; ++a )
            {
              if( edgesToNodes[edgeIndex][a] == nodeID )
              {

                if( verboseLevel() > 1 )
                  std::cout<<edgesToNodes[edgeIndex][a];

                edgesToNodes[edgeIndex][a] = newNodeIndex;
                nodesToEdges[nodeID].erase( edgeIndex );

                if( verboseLevel() > 1 )
                  std::cout<<"->"<<edgesToNodes[edgeIndex][a]<<", ";

              }
              else if( verboseLevel() > 1 )
                std::cout<<edgesToNodes[edgeIndex][a]<<", ";

              nodesToEdges[edgesToNodes[edgeIndex][a]].insert( edgeIndex );
              modifiedObjects.modifiedNodes.insert(edgesToNodes[edgeIndex][a]);
            }
            if( verboseLevel() > 1 )
              std::cout<<")";
          }
          if( verboseLevel() > 1 )
            std::cout<<", ";
        }
        if( verboseLevel() > 1 )
          std::cout<<")"<<std::endl;
      } // for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
    } // if( location==1 )
  } // for( map<std::pair<CellBlockSubRegion*, localIndex>, int>::const_iterator iter_elem = elemLocations.begin()
}


void SurfaceGenerator::MapConsistencyCheck( const localIndex nodeID,
                                            NodeManager const & nodeManager,
                                            EdgeManager const & edgeManager,
                                            FaceManager const & faceManager,
                                            ElementRegionManager const & elementManager,
                                            const map< std::pair<CellElementSubRegion*, localIndex >, int>& elemLocations )
{
  //**************************************************************************
  // THIS IS ALL JUST CONSISTENCY CHECKING
  //**************************************************************************


  arrayView1d<set<localIndex> > const & nodesToEdges = nodeManager.edgeList();
  arrayView1d<set<localIndex> > const & nodesToFaces = nodeManager.faceList();
  ArrayOfArraysView<localIndex const> const & nodesToElementRegions = nodeManager.elementRegionList();
  ArrayOfArraysView<localIndex const> const & nodesToElementSubRegions = nodeManager.elementSubRegionList();
  ArrayOfArraysView<localIndex const> const & nodesToElementIndex = nodeManager.elementList();


  arrayView2d<localIndex> const & edgesToNodes = edgeManager.nodeList();
  arrayView1d<set<localIndex> > const & edgesToFaces = edgeManager.faceList();

  OrderedVariableOneToManyRelation const & facesToNodes = faceManager.nodeList();
  OrderedVariableOneToManyRelation const & facesToEdges = faceManager.edgeList();
  arrayView2d<localIndex> const & facesToElementRegions = faceManager.elementRegionList();
  arrayView2d<localIndex> const & facesToElementSubRegions = faceManager.elementSubRegionList();
  arrayView2d<localIndex> const & facesToElementIndex = faceManager.elementList();


#if 1
  if( verboseLevel() > 2 )
  {
    std::cout<<"CONSISTENCY CHECKING OF THE MAPS"<<std::endl;

    for( map< std::pair<CellElementSubRegion*, localIndex >, int>::const_iterator iter_elem=elemLocations.begin() ;
         iter_elem!=elemLocations.end() ; ++iter_elem )
    {
      const std::pair<CellElementSubRegion*, localIndex >& elem = iter_elem->first;

      CellElementSubRegion& elemSubRegion = *(elem.first);
      const localIndex elemIndex = elem.second;

      arrayView2d<localIndex> & elemsToNodes = elemSubRegion.nodeList();
      arrayView2d<localIndex> & elemsToFaces = elemSubRegion.faceList();


      std::set<localIndex> elemNodes;


      std::cout<<"Element "<<elemIndex<<"\n";
      std::cout<<" elementToNodes = ";
      for( int a=0 ; a<8 ; ++a )
      {
        elemNodes.insert( elemsToNodes( elemIndex, a ));
        std::cout<< elemsToNodes( elemIndex, a )<<", ";
      }
      std::cout<<std::endl;

      std::cout<<" elementToFaces->edges->nodes = ";


      // Now we do a loop over the facelist and process all the faces
      for( int kf=0 ; kf<elemSubRegion.numFacesPerElement() ; ++kf )
      {
        std::set<localIndex> faceNodes;

        localIndex faceIndex  = elemsToFaces( elemIndex, kf );

        if( kf>0 )
          std::cout<<"                              = ";


        std::cout<<faceIndex<<"( ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex faceNodeID = facesToNodes[faceIndex][b];
          faceNodes.insert( faceNodeID );
          if( elemNodes.count( faceNodeID ) == 0 && kf<elemSubRegion.numFacesPerElement() )
            std::cout<<"*";
          std::cout<<faceNodeID<<",";
        }
        std::cout<<" )      ";



        std::cout<<faceIndex<<"[ ";
        for( int b=0 ; b<4 ; ++b )
        {
          localIndex edgeIndex = facesToEdges[faceIndex][b];
          std::cout<<edgeIndex<<"( ";
          for( int c=0 ; c<2 ; ++c )
          {
            localIndex edgeNodeID = edgesToNodes( edgeIndex, c );
            if( elemNodes.count( edgeNodeID ) == 0  && kf<elemSubRegion.numFacesPerElement() )
              std::cout<<"*";
            if( faceNodes.count( edgeNodeID ) == 0 )
              std::cout<<"#";
            std::cout<<edgeNodeID<<",";
          }
          std::cout<<" ), ";
        }
        std::cout<<" ] \n";

      }
      std::cout<<std::endl;

    }

  }

  if( verboseLevel() > 2 )
  {
    // nodeToEdge
    array1d<std::set<localIndex> > inverseEdgesToNodes( nodeManager.size() );

    for( localIndex ke=0 ; ke<edgeManager.size() ; ++ke )
    {
      for( localIndex b= 0 ; b<edgesToNodes.size( 1 ) ; ++b )
      {
        localIndex nodeIndex = edgesToNodes( ke, b );
        inverseEdgesToNodes[nodeIndex].insert( ke );
      }
    }
    std::cout<<"Check NodeToEdge:  nodesToEdges  inverseEdgesToNodes"<<std::endl;
    for( localIndex a=0 ; a<nodeManager.size() ; ++a )
    {
      std::cout<<"m_nodesToEdges["<<a<<"] = ( ";
      for( set<localIndex>::const_iterator iedge=nodesToEdges[a].begin() ;
           iedge!=nodesToEdges[a].end() ; ++iedge )
      {
        if( inverseEdgesToNodes[a].count( *iedge ) == 0 )
          std::cout<<"*";

        std::cout<<*iedge<<", ";
      }
      std::cout<<")    (";

      for( std::set<localIndex>::const_iterator iedge=inverseEdgesToNodes[a].begin() ;
           iedge!=inverseEdgesToNodes[a].end() ; ++iedge )
      {
        if( nodeManager.edgeList()[a].count( *iedge ) == 0 )
          std::cout<<"*";

        std::cout<<*iedge<<", ";
      }
      std::cout<<")"<<std::endl;
    }


  }

  if( verboseLevel() > 2 )
  {
    // nodeToFace
    array1d<std::set<localIndex> > inverseFacesToNodes( nodeManager.size() );
    for( localIndex kf=0 ; kf<faceManager.size() ; ++kf )
    {
      for( localIndex_array::const_iterator b=facesToNodes[kf].begin() ;
           b!=facesToNodes[kf].end() ; ++b )
      {
        inverseFacesToNodes[*b].insert( kf );
      }
    }
    std::cout<<"Check NodeToFace:  nodesToFaces  inverseFacesToNodes"<<std::endl;
    for( localIndex a=0 ; a<nodeManager.size() ; ++a )
    {
      std::cout<<"m_nodeToFaceMap["<<a<<"] = ( ";
      for( set<localIndex>::const_iterator iface=nodesToFaces[a].begin() ;
           iface!=nodesToFaces[a].end() ; ++iface )
      {
        if( inverseFacesToNodes[a].count( *iface ) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")    (";

      for( std::set<localIndex>::const_iterator iface=inverseFacesToNodes[a].begin() ;
           iface!=inverseFacesToNodes[a].end() ; ++iface )
      {
        if( nodesToFaces[a].count( *iface ) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")"<<std::endl;
    }

  }



  if( verboseLevel() > 2 )
  {


    // nodeToElement
    array1d<std::set<std::pair<CellElementSubRegion const *, localIndex > > > inverseElemsToNodes( nodeManager.size() );
    for( localIndex er=0 ; er<elementManager.numRegions() ; ++er )
    {
      ElementRegion const & elemRegion = *(elementManager.GetRegion( er ));
      for( localIndex esr=0 ; esr<elemRegion.numSubRegions() ; ++esr )
      {
        CellElementSubRegion const * const subRegion = elemRegion.GetSubRegion<CellElementSubRegion>( esr );
        if( subRegion != nullptr )
        {
          arrayView2d<localIndex> const & elemsToNodes = subRegion->nodeList();
          for( localIndex k=0 ; k<subRegion->size() ; ++k )
          {
            std::pair<CellElementSubRegion const *, localIndex > elem = std::make_pair( subRegion, k );

            for( localIndex a=0 ; a<elemsToNodes.size( 1 ) ; ++a )
            {
              inverseElemsToNodes[elemsToNodes( k, a )].insert( elem );
            }
          }
        }
      }
    }



    std::cout<<"Check NodeToElem: nodesToElems  inverseElemsToNodes "<<std::endl;


    for( localIndex a=0 ; a<nodeManager.size() ; ++a )
    {

      std::set< std::pair<CellElementSubRegion const *, localIndex> > nodeToElements;
      for( localIndex k=0 ; k<nodesToElementRegions.sizeOfArray(a) ; ++k )
      {
        if( nodesToElementRegions[a][k]!=-1 && nodesToElementSubRegions[a][k]!=-1 && nodesToElementIndex[a][k]!=-1 )
        {
          nodeToElements.insert( std::make_pair( elementManager.GetRegion( nodesToElementRegions[a][k] )->
                                                 GetSubRegion<CellElementSubRegion>( nodesToElementSubRegions[a][k] ),
                                                 nodesToElementIndex[a][k] ) );
        }
      }


      std::cout<<"m_NodeToElementMap["<<a<<"] = ( ";
      for( std::set<std::pair< CellElementSubRegion const *, localIndex > >::iterator
           ielem=nodeToElements.begin() ; ielem!=nodeToElements.end() ; ++ielem )
      {
        if( inverseElemsToNodes[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set<std::pair<CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=inverseElemsToNodes[a].begin() ;
           ielem!=inverseElemsToNodes[a].end() ; ++ielem )
      {
        if( nodeToElements.count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }


    // edgeToFace
    array1d<std::set<localIndex> > inverseFacesToEdges( edgeManager.size() );
    for( localIndex kf=0 ; kf<faceManager.size() ; ++kf )
    {
      for( localIndex_array::const_iterator b=facesToEdges[kf].begin() ;
           b!=facesToEdges[kf].end() ; ++b )
      {
        inverseFacesToEdges[*b].insert( kf );
      }
    }
    std::cout<<"Check EdgeToFace: edgesToFaces  inverseFacesToEdges "<<std::endl;
    for( localIndex ke=0 ; ke<edgeManager.size() ; ++ke )
    {
      std::cout<<"m_edgesToFaces["<<ke<<"] = ( ";
      for( set<localIndex>::const_iterator
           iface=edgeManager.faceList()[ke].begin() ;
           iface!=edgeManager.faceList()[ke].end() ; ++iface )
      {
        if( inverseFacesToEdges[ke].count( *iface ) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")    (";

      for( std::set<localIndex>::const_iterator iface=inverseFacesToEdges[ke].begin() ;
           iface!=inverseFacesToEdges[ke].end() ; ++iface )
      {
        if( edgeManager.faceList()[ke].count( *iface ) == 0 )
          std::cout<<"*";

        std::cout<<*iface<<", ";
      }
      std::cout<<")"<<std::endl;
    }

    // faceToElement
    array1d<std::set<std::pair<CellElementSubRegion const *, localIndex > > > inverseElemsToFaces( faceManager.size() );
    for( localIndex er=0 ; er<elementManager.numRegions() ; ++er )
    {
      ElementRegion const & elemRegion = *(elementManager.GetRegion( er ));
      for( localIndex esr=0 ; esr<elemRegion.numSubRegions() ; ++esr )
      {
        CellElementSubRegion const * const subRegion = elemRegion.GetSubRegion<CellElementSubRegion>( esr );
        if( subRegion != nullptr )
        {
          arrayView2d<localIndex> const & elemsToNodes = subRegion->nodeList();
          arrayView2d<localIndex> const & elemsToFaces = subRegion->faceList();

          for( localIndex k=0 ; k<subRegion->size() ; ++k )
          {
            std::pair<CellElementSubRegion const *, localIndex > elem = std::make_pair( subRegion, k );

            for( localIndex a=0 ; a<elemsToFaces.size( 1 ) ; ++a )
            {
              const localIndex faceID = elemsToFaces( k, a );
              inverseElemsToFaces[faceID].insert( elem );

  //            if( parentFaceIndex[faceID] != -1 )
  //            {
  //              inverseElemsToFaces[parentFaceIndex[faceID]].insert(elem);
  //            }
            }
          }
        }
      }
    }
    std::cout<<"Check FacesToElem: facesToElems  inverseElemsToFaces "<<std::endl;
    for( localIndex a=0 ; a<faceManager.size() ; ++a )
    {

      std::vector< std::pair<CellElementSubRegion const *, localIndex> > faceToElements;
      for( localIndex k=0 ; k<facesToElementRegions.size( 1 ) ; ++k )
      {
        // TODO This only works for a single region
        if( facesToElementRegions[a][k] != -1 )
        {
          faceToElements.push_back( std::make_pair( elementManager.GetRegion( facesToElementRegions[a][k] )->
                                                    GetSubRegion<CellElementSubRegion>( facesToElementSubRegions[a][k] ),
                                                    facesToElementIndex[a][k] ) );
        }
      }


      std::cout<<"m_FaceToElementMap["<<a<<"] = ( ";

      for( std::vector<std::pair<CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=faceToElements.begin() ;
           ielem!=faceToElements.end() ; ++ielem )
      {
        if( inverseElemsToFaces[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set<std::pair<CellElementSubRegion const *, localIndex > >::const_iterator ielem=inverseElemsToFaces[a].begin() ;
           ielem!=inverseElemsToFaces[a].end() ; ++ielem )
      {

        if( faceToElements.size() == 2 )
        {
          if( (faceToElements[0] != *ielem) && (faceToElements[1] != *ielem) )
            std::cout<<"*";
        }
        else if( faceToElements.size() )
        {
          if( (faceToElements[0] != *ielem)  )
            std::cout<<"*";
        }
        else
        {
          std::cout<<"****";
        }


        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }
  }
//  CorrectSplitNodalMass(nodeManager, nodeID, nodeManager.m_childIndices[nodeID][0]);
#endif
}



realT SurfaceGenerator::CalculateKinkAngle ( const localIndex edgeID,
                                             const NodeManager & nodeManager,
                                             EdgeManager & edgeManager,
                                             FaceManager & faceManager )
{
  localIndex_array faces;
  // realT kinkAngle;

  for( auto iface : edgeManager.faceList()[edgeID] )
  {
    if( faceManager.m_isExternal[iface] == 1 )
      faces.push_back( iface );
  }

  if( faces.size() != 2 )
  {
    return(-1.0);
  }
  else
//  {
////    // First check if the two faces are parent-child pairs
////    if (faceManager.m_parentIndex[faces[0]]==faces[1] || faceManager.m_parentIndex[faces[1]]==faces[0] )
////    {
////      return(0.0);
////    }
//
//    R1Tensor vecFace[3];
//    faceManager.InFaceVectorNormalToEdge(nodeManager, edgeManager, faces[0], edgeID, vecFace[0]);
//    faceManager.InFaceVectorNormalToEdge(nodeManager, edgeManager, faces[1], edgeID, vecFace[1]);
//    vecFace[2] = vecFace[0];
//    vecFace[2] += vecFace[1];
//    vecFace[2] /= 2.0;
//
//    kinkAngle = acos(Dot(vecFace[0],vecFace[1])*0.999999) / 3.141592653589793238462 * 180.0;
//
//    R1Tensor vecFaceNorm;
//    vecFaceNorm = faceManager.FaceNormal(nodeManager, faces[0]);
//    vecFaceNorm  += faceManager.FaceNormal(nodeManager, faces[1]);
//    vecFaceNorm /= 2.0;
//
//    if (Dot(vecFace[2], vecFaceNorm) < 0.0)
//      kinkAngle = 360.0 - kinkAngle;
//
//    return(kinkAngle);
//
//  }
    return 1e100;
}

void SurfaceGenerator::CalculateKinkAngles ( FaceManager & faceManager,
                                             EdgeManager & edgeManager,
                                             NodeManager & nodeManager,
                                             ModifiedObjectLists& modifiedObjects,
                                             const bool prefrac )
{
  arrayView1d<real64> & kinkAngle = edgeManager.getReference<real64_array>( "kinkAngle" );

  if( prefrac )
  {
    for( localIndex edgeID = 0 ; edgeID < edgeManager.size() ; ++edgeID )
    {
      kinkAngle[edgeID] = CalculateKinkAngle( edgeID, nodeManager, edgeManager, faceManager );
    }
  }
  else
  {
    for( std::set<localIndex>::const_iterator i=modifiedObjects.newEdges.begin() ; i!=modifiedObjects.newEdges.end() ; ++i )
    {
      kinkAngle[*i] = CalculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
    for( std::set<localIndex>::const_iterator i=modifiedObjects.modifiedEdges.begin() ; i!=modifiedObjects.modifiedEdges.end() ; ++i )
    {
      kinkAngle[*i] = CalculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
  }
}



void SurfaceGenerator::IdentifyRupturedFaces( NodeManager & nodeManager,
                                              EdgeManager & edgeManager,
                                              FaceManager & faceManager,
                                              ElementRegionManager & elementManager,
                                              SpatialPartition& partition,
                                              const bool prefrac )
{
  arrayView1d<integer> const & isEdgeGhost = edgeManager.GhostRank();
//  const integer_array& layersEdgeFromBoundary = edgeManager.getReference<integer_array>("LayersFromDomainBoundary");
//  localIndex_array& primaryCandidateFace = faceManager.getReference<localIndex_array>("primaryCandidateFace");
//  primaryCandidateFace = std::numeric_limits<localIndex>::max();

  int rank;
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );

  //"Heal" faces that were marked but not split.
//  if (m_failCriterion >0)
//  {
//    integer_array& ruptureState = faceManager.getReference<int>("ruptureState");
//    for (localIndex iFace = 0; iFace < faceManager.DataLengths(); ++iFace)
//    {
//      if (faceManager.m_isExternal[iFace] == 0)
//      {
//        ruptureState[iFace] = 0;
//      }
//    }
//  }
// The healing has been done in the separation driver.  Also the if criterion here should be failCriterion==1



  // We use the color map scheme because we can mark a face to be rupture ready from a partition where the face is a
  // ghost.

  // Process interior edges
  {
    ModifiedObjectLists modifiedObjects;

    for( localIndex iEdge = 0 ; iEdge != edgeManager.size() ; ++iEdge )
    {
      if( isEdgeGhost[iEdge] < 0 ) //&& layersEdgeFromBoundary[iEdge]>1 )
      {
        int edgeMode = CheckEdgeSplitability( iEdge,
                                              nodeManager,
                                              faceManager,
                                              edgeManager,
                                              prefrac );
        if( edgeMode == 0 || edgeMode == 1 ) // We need to calculate SIF
        {
          R1Tensor vecTipNorm, vecTip;
          localIndex trailFaceID = 0;
          realT SIF = CalculateEdgeSIF( iEdge, trailFaceID,
                                        nodeManager,
                                        edgeManager,
                                        faceManager,
                                        elementManager,
                                        vecTipNorm,
                                        vecTip );

          if( SIF > MinimumToughnessOnEdge( iEdge, nodeManager, edgeManager, faceManager ) * 0.5 )
          {
            MarkRuptureFaceFromEdge( iEdge, trailFaceID,
                                     nodeManager,
                                     edgeManager,
                                     faceManager,
                                     vecTipNorm,
                                     vecTip,
                                     modifiedObjects,
                                     edgeMode );
          }
        }
      }
    }
  }

  // Process near boundary edges
  {
//    for( int color=0 ; color<partition.NumColor() ; ++color )
    {
      ModifiedObjectLists modifiedObjects;
//      if( partition.Color() == color )
      {
        for( localIndex iEdge = 0 ; iEdge != edgeManager.size() ; ++iEdge )
        {

          if( isEdgeGhost[iEdge] < 0 ) //&& layersEdgeFromBoundary[iEdge]<=1 )
          {
            int edgeMode = CheckEdgeSplitability( iEdge,
                                                  nodeManager,
                                                  faceManager,
                                                  edgeManager,
                                                  prefrac );
            if( edgeMode == 0 || edgeMode == 1 ) // We need to calculate SIF
            {
              R1Tensor vecTipNorm, vecTip, vecEdge, direction;
              localIndex trailFaceID = 0;
              realT SIF = CalculateEdgeSIF( iEdge, trailFaceID,
                                            nodeManager,
                                            edgeManager,
                                            faceManager,
                                            elementManager,
                                            vecTipNorm,
                                            vecTip );

              if( SIF >  MinimumToughnessOnEdge( iEdge, nodeManager, edgeManager, faceManager ) * 0.5 ) // && edgeMode
                                                                                                        // == 1)
              {
                MarkRuptureFaceFromEdge( iEdge, trailFaceID,
                                         nodeManager,
                                         edgeManager,
                                         faceManager,
                                         vecTipNorm,
                                         vecTip,
                                         modifiedObjects,
                                         edgeMode );
              }
            }
          }
        }
//        partition.ModifyGhostsAndNeighborLists( modifiedObjects );
      }
    }
  }

}

realT SurfaceGenerator::CalculateEdgeSIF( const localIndex edgeID,
                                          localIndex& trailFaceID,
                                          NodeManager & nodeManager,
                                          EdgeManager & edgeManager,
                                          FaceManager & faceManager,
                                          ElementRegionManager & elementManager,
                                          R1Tensor& vecTipNorm,
                                          R1Tensor& vecTip )
{
#if 0
  realT rval;
  localIndex nExternalFaces = 0;
  localIndex_array faceInvolved;
  real64_array& SIF_I = edgeManager.getReference<real64_array>( "SIF_I" );
  real64_array& SIF_II = edgeManager.getReference<real64_array>( "SIF_II" );
  real64_array& SIF_III = edgeManager.getReference<real64_array>( "SIF_III" );

  array1d< set<localIndex> > const & nodesToEdges = nodeManager.edgeList();
  array1d<set<localIndex> > const & nodesToElementRegion = nodeManager.elementRegionList();
  array1d<set<localIndex> > const & nodesToElementSubRegion = nodeManager.elementSubRegionList();
  array1d<set<localIndex> > const & nodesToElementIndex = nodeManager.elementList();


  r1_array const & X = nodeManager.referencePosition();

  Array2dT<localIndex> const & edgesToNodes = edgeManager.nodeList();

  localIndex_array const & faceParentIndex = faceManager.getReference<localIndex_array>( viewKeys.faceParentIndex );
  array1d< localIndex_array > const & facesToNodes = faceManager.nodeList();
  array1d< localIndex_array > const & facesToEdges = faceManager.edgeList();
//  real64_array const * const ppFlowPressure = faceManager.getReferencePointer<FieldInfo::pressure>();
//  const real64_array& faceArea = faceManager.getReference<realT>("faceArea");

  SIF_I[edgeID] = 0.0;
  SIF_II[edgeID] = 0.0;
  SIF_III[edgeID] = 0.0;

  //Sanity check
  for( auto iface : edgeManager.faceList()[edgeID] )
  {
    if( faceManager.m_isExternal[iface] == 1 )
    {
      nExternalFaces++;
      faceInvolved.push_back( iface );
    }
  }
  if( nExternalFaces > 2 )
  {
    GEOS_ERROR( "Error! This is a singular edge, not a tip.  This should not happen!" );
  }

  // Figure out the two fracture faces connected to this edge
  localIndex faceA, faceAp;
  if( (faceParentIndex[faceInvolved[0]] == LOCALINDEX_MAX && faceParentIndex[faceInvolved[1]] == faceInvolved[0]) ||
      (faceParentIndex[faceInvolved[1]] == LOCALINDEX_MAX && faceParentIndex[faceInvolved[0]] == faceInvolved[1]) )
  {
    faceA = faceInvolved[0];
    faceAp = faceInvolved[1];
  }
  else
  {
    char msg[200];
    sprintf( msg, "Error! Edge %d has two external faces, but the parent-child relationship is wrong.", int(edgeID));
    GEOS_ERROR( msg );
  }

  trailFaceID = faceParentIndex[faceInvolved[0]];

  // We define three unit vectors
  // vecEdge: pointing from node 0 to node 1 along the tip edge
  // vecTip: pointing from the opening into the solid
  // vecTipNorm: normal of the one of the fracture faces;  vecTip X vecTipNorm should point to the direction of vecEdge
  // These new definitions are consistent with those in 2D.  See the illustration in Fractunator2D.

  real64 temp;
  R1Tensor r1Temp;
  R1Tensor faceCenterA;
  R1Tensor faceCenterAp;
  computationalGeometry::
  Centroid_3DPolygon( facesToNodes[faceA],
                      X,
                      faceCenterA,
                      vecTipNorm );

  computationalGeometry::
  Centroid_3DPolygon( facesToNodes[faceAp],
                      X,
                      faceCenterAp,
                      r1Temp );

  vecTipNorm -= r1Temp;
  vecTipNorm.Normalize();

  R1Tensor vecEdge;
  computationalGeometry::
  VectorDifference( X,
                    edgesToNodes[edgeID][0],
                    edgesToNodes[edgeID][1],
                    vecEdge );
  vecEdge.Normalize();

  vecTip.Cross( vecTipNorm, vecEdge );
  vecTip.Normalize();
  R1Tensor v0, v1;
  computationalGeometry::VectorMean<2>( X, edgesToNodes[edgeID], v0 );
  v0 -= faceCenterA;

  if( Dot( v0, vecTip ) < 0 )
    vecTip *= -1.0;
  if( Dot( Cross( vecTip, vecTipNorm ), vecEdge ) < 0 )
  {
    vecTipNorm *= -1;
    faceA = faceInvolved[1];
    faceAp = faceInvolved[0];
  }


  //Now we need to figure out if a special situation applies to this edge
  // where the fracture face is a quad and three of the nodes are still pinched
  // We use a different algorithm for this special situation.

  bool threeNodesPinched( false );
  localIndex_array openNodeID;

  if( facesToNodes[faceA].size() == 4 )  // Only quads have this problem
  {
    int numSharedNodes = 2;
    localIndex_array lNodeFaceA, lNodeFaceAp;

    lNodeFaceA = facesToNodes[faceA];
    lNodeFaceAp = facesToNodes[faceAp];


    //We remove all the shared nodes and the one remains should be the open one.
    lNodeFaceAp.erase( std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgesToNodes[edgeID][0] ));
    lNodeFaceAp.erase( std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgesToNodes[edgeID][1] ));
    lNodeFaceA.erase( std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgesToNodes[edgeID][0] ));
    lNodeFaceA.erase( std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgesToNodes[edgeID][1] ));

    for( auto j : facesToNodes[faceA] )
    {
      localIndex iNd = j;
      if( iNd != edgesToNodes[edgeID][0] && iNd != edgesToNodes[edgeID][1] )
      {
        if( std::find( facesToNodes[faceAp].begin(), facesToNodes[faceAp].end(), iNd ) != facesToNodes[faceAp].end())
        {
          numSharedNodes++;
          lNodeFaceA.erase( std::find( lNodeFaceA.begin(), lNodeFaceA.end(), iNd ));
          lNodeFaceAp.erase( std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), iNd ));
        }
      }
    }

    if( numSharedNodes == 4 )
    {
      GEOS_ERROR( "Error.  The fracture face has four shared nodes with its child.  This should not happen." );
    }
    else if( numSharedNodes == 3 )
    {
      threeNodesPinched = true;
      if( lNodeFaceA.size() != 1 || lNodeFaceAp.size() != 1 )
      {
        GEOS_ERROR( "Error. These two faces share three nodes but the number of remaining nodes is not one.  Something is wrong" );
      }
      else
      {
        openNodeID.push_back( lNodeFaceA[0] );
        openNodeID.push_back( lNodeFaceAp[0] );
      }
    }
  }

  // Now we need to identify which node on the edge is the convex point and which one is the concave corner.  The convex
  // node must share an edge with the open node.
  localIndex convexCorner( std::numeric_limits<localIndex>::max());
  if( threeNodesPinched )
  {
    localIndex iNd, jNd;
    iNd = edgesToNodes[edgeID][0];
    jNd = edgesToNodes[edgeID][1];
    for( auto const j : facesToEdges[faceA] )
    {
      localIndex edge = j;
      if((openNodeID[0] == edgesToNodes[edge][0] && iNd == edgesToNodes[edge][1]) ||
         (openNodeID[0] == edgesToNodes[edge][1] && iNd == edgesToNodes[edge][0])
         )
      {
        convexCorner = iNd;
        break;
      }
      if((openNodeID[0] == edgesToNodes[edge][0] && jNd == edgesToNodes[edge][1]) ||
         (openNodeID[0] == edgesToNodes[edge][1] && jNd == edgesToNodes[edge][0])
         )
      {
        convexCorner = jNd;
        break;
      }
    }

    if( convexCorner == std::numeric_limits<localIndex>::max())
      GEOS_ERROR( "Error.  This is a three-node-pinched edge but I cannot find the convex corner" );

  }



  // Calculate element forces acting on this edge.  Need to add nodal forces from two nodes up.
  //An element has to be within the range of this edge to be included.
  //For the threeNodesPinched case, we only use the force on the node at the convex point, not the concave point.  The
  // force at the former is ususally greater, so we just pick the great one instead of doing a geometrical check.

  localIndex nElemEachSide[2], nGhostElem;
  R1Tensor fNodeO = static_cast<R1Tensor>(0.0);
  nElemEachSide[0] = 0;
  nElemEachSide[1] = 0;
  R1Tensor xEdge;
  computationalGeometry::VectorMean<2>( X, edgesToNodes[edgeID], xEdge );

  realT GdivBeta = 0;  // Need this for opening-based SIF

  if( !threeNodesPinched )
  {
    for( localIndex a=0 ; a<edgesToNodes.size( 1 ) ; ++a ) // Loop through the two nodes
    {
      localIndex nodeID = edgesToNodes( edgeID, a );

      for( localIndex k=0 ; k<nodesToElementRegion[nodeID].size() ; ++k )
      {

        CellElementSubRegion * elementRegion = elementManager.GetRegion( nodesToElementRegion[nodeID][k] )->
                                             GetSubRegion( nodesToElementSubRegion[nodeID][k] );
        localIndex iEle = nodesToElementIndex[nodeID][k];
        R1Tensor fN, xEle;


        xEle = elementRegion->GetElementCenter( iEle, nodeManager );

        realT ndist, udist, segmentLength;
        R1Tensor ptPrj;
        GeometryUtilities::ProjectPointToLineSegment( X[nodeID],
                                                      X[edgesToNodes( edgeID, 1-a )],
                                                      xEle,
                                                      ndist, udist, segmentLength,
                                                      ptPrj );
        if( udist <= edgeManager.EdgeLength( nodeManager, edgeID ) && udist > 0.0 )
        {
          elementRegion->CalculateNodalForcesFromOneElement( nodeID, iEle, nodeManager, fN );
          GdivBeta += elementRegion->ElementGDivBeta( iEle );

          xEle -= xEdge;

          if( Dot( xEle, vecTipNorm ) > 0 )
          {
            nElemEachSide[0] += 1;
            fNodeO += fN;
          }
          else
          {
            nElemEachSide[1] +=1;
            fNodeO -= fN;
          }
        }
      } // Loop through all elements connected to this node

    }  // loop over two nodes on the edge
  }
  else
  {

    localIndex nodeID = convexCorner;
    fNodeO = 0.0;

    for( std::set< std::pair<CellElementSubRegion *, localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[nodeID].begin() ;
         k!=nodeManager.m_toElementsRelation[nodeID].end() ; ++k )
    {
      CellElementSubRegion * elementRegion = k->first;
      localIndex iEle = k->second;
      R1Tensor fN, xEle;

      xEle = elementRegion->GetElementCenter( iEle, nodeManager );

      {
        elementRegion->CalculateNodalForcesFromOneElement( nodeID, iEle, nodeManager, fN );
        GdivBeta += elementRegion->ElementGDivBeta( iEle );

        xEle -= xEdge;

        if( Dot( xEle, vecTipNorm ) > 0 )
        {
          nElemEachSide[0] += 1;
          fNodeO += fN;
        }
        else
        {
          nElemEachSide[1] +=1;
          fNodeO -= fN;
        }
      }
    } // Loop through all elements connected to this node
  }

  if( nElemEachSide[0]>=1 && nElemEachSide[1]>=1 )
    fNodeO /= 2.0;
  //We have contributions from both sides.  The two sizes are the two sides of the fracture plane.  If the fracture face
  // is on domain boundary, it's possible to have just one side.
  if( nElemEachSide[0] + nElemEachSide[1] >= 1 )
    GdivBeta /= (nElemEachSide[0] + nElemEachSide[1]);

  localIndex tipFaces[2];
  tipFaces[0] = faceA;
  tipFaces[1] = faceAp;

  // We have to subtract the nodal force at other nodes (trailing nodes) on these two open faces to take into account
  // the effects of surface traction along the fracture.
  // Finding the two trailing nodes on a hex mesh is pretty straightforward, while it is cumbersome to do in tet mesh
  // For the threeNodesPinched case, this should be the open node.
  R1Tensor fFaceA[2];
  // fFaceA is actually the average nodal force on each face.  Assuming homogeneous meshing.


  // If the two external faces connected to a trailing edge are not coplanar, then we have the risk of incomplete
  // topology.
  // In that case, we use a displacement/opening based method, not VCCT.
  bool incompleteTrailingEdgeTopology = 0;

#if 1
  for( localIndex i=0 ; i<2 ; ++i )
  {
    integer_array trailingNodes;
    localIndex trailingEdge;
    trailingNodes.clear();
    if( threeNodesPinched )
    {
      trailingNodes.push_back( openNodeID[i] );
    }
    else
    {
      localIndex faceID = tipFaces[i];
      nGhostElem = 0;
      fFaceA[i] = 0.0;

      for( auto j : facesToNodes[faceID] )
      {
        if( j != edgesToNodes( edgeID, 0 ) && j != edgesToNodes( edgeID, 1 ) ) // This is not a node along the tip edge
        {
          trailingNodes.push_back( j );
        }
      }

      if( trailingNodes.size() > 2 || trailingNodes.size() == 0 )
      {
        GEOS_ERROR( "Fatal error in finding nodes behind tip edge." );
      }
      else if( trailingNodes.size() == 1 )  // Need some work to find the other node
      {
        // First find a edge that is connected to this node and parallel to the tip edge
        realT maxCosAngle = 0.0;
        localIndex pickedTrailingEdge = std::numeric_limits<localIndex>::max();
        for( auto iedge : nodesToEdges[trailingNodes[0]] )
        {
          R1Tensor xTrailingEdge;
          edgeManager.EdgeCenter( nodeManager, *iedge, xTrailingEdge );

          realT ndist, udist, segmentLength;
          R1Tensor ptPrj;
          GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[edgesToNodes[edgeID][0]],
                                                        (*nodeManager.m_refposition)[edgesToNodes[edgeID][1]],
                                                        xTrailingEdge,
                                                        ndist, udist, segmentLength,
                                                        ptPrj );
          if( udist <= edgeManager.EdgeLength( nodeManager, edgeID ) && udist > 0.0 )
          {
            R1Tensor vEdge;
            edgeManager.EdgeVector( nodeManager, *iedge, vEdge );
            vEdge /= edgeManager.EdgeLength( nodeManager, *iedge );
            realT cosEdge = std::fabs( Dot( vEdge, vecEdge ));
            if( cosEdge > maxCosAngle )
            {
              maxCosAngle = cosEdge;
              pickedTrailingEdge = *iedge;
            }
          }
        }
        if( maxCosAngle > 0.75 )
          trailingNodes.push_back( edgesToNodes[pickedTrailingEdge][0] + edgesToNodes[pickedTrailingEdge][1] - trailingNodes[0] );
      }
    }

    if( trailingNodes.size() == 2 )
    {
      trailingEdge = edgeManager.FindEdgeFromNodeIDs( trailingNodes[0], trailingNodes[1], nodeManager );
      if( trailingEdge > edgeManager.DataLengths())
      {
        // GEOS_ERROR("Error: I have two trailing nodes but cannot find a trailing edge.");
        int rank=0;
        MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
        std::cout << "Cannot find trailing edge (edge=" << edgeID << ", rank=" << rank <<   "  )" << std::endl;
        return 0.0;
      }


      localIndex_array extFacesOnTrailingEdge;
      for( auto iface : edgesToFaces[trailingEdge] )
      {
        if( faceManager.m_isExternal[iface] == 1 )
          extFacesOnTrailingEdge.push_back( iface );
      }

      if( extFacesOnTrailingEdge.size() != 2 )
      {
        incompleteTrailingEdgeTopology = 1;
      }
      else
      {
        R1Tensor extFaceNormal[2];
        for( localIndex j = 0 ; j < 2 ; ++j )
        {
          extFaceNormal[j] = faceManager.FaceNormal( nodeManager, extFacesOnTrailingEdge[j], true );
        }

        if( std::fabs( Dot( extFaceNormal[0], extFaceNormal[1] )) < 0.9 ) //The two faces are not coplanar.
        {
          incompleteTrailingEdgeTopology = 1;
        }
      }
    }

    for( localIndex iR = 0 ; iR < trailingNodes.size() ; ++iR )
    {
      localIndex iNd = trailingNodes[iR];
      nGhostElem = 0;
      localIndex nValidElements = 0;
      R1Tensor fNode = static_cast<R1Tensor>(0.0);

      // The unbalanced force from a summation of element contributions should be the external forces, could be fluid
      // pressure or try traction.
      // In this way, we don't need to worry about the how the traction inside a fracture was applied.
      for( std::set< std::pair<CellElementSubRegion *, localIndex> >::const_iterator k=nodeManager.m_toElementsRelation[iNd].begin() ;
           k!=nodeManager.m_toElementsRelation[iNd].end() ; ++k )
      {
        CellElementSubRegion * elementRegion = k->first;
        localIndex iEle = k->second;
        integer_array& elem_is_ghost = elementRegion->getReference<FieldInfo::ghostRank>();
        R1Tensor fN, xEle;

        xEle = elementRegion->GetElementCenter( iEle, nodeManager );

        realT ndist, udist, segmentLength;
        R1Tensor ptPrj;
        GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[edgesToNodes[edgeID][0]],
                                                      (*nodeManager.m_refposition)[edgesToNodes[edgeID][1]],
                                                      xEle,
                                                      ndist, udist, segmentLength,
                                                      ptPrj );
        if(  ( udist <= edgeManager.EdgeLength( nodeManager, edgeID ) && udist > 0.0 ) || threeNodesPinched )
        {
          elementRegion->CalculateNodalForcesFromOneElement( iNd, iEle, nodeManager, fN );
          fNode += fN;
          nValidElements++;
          if( elem_is_ghost[iEle] >= 0 )
            nGhostElem +=1;
        }

      }

      fFaceA[i] += fNode;

      // I am merging this logic into one that unified the ghost situation and the intersection situation
//      if (nGhostElem == nValidElements)
//      {
//        // This is a conservative criterion. It's possible that the all elements connected to trailing edge are ghost
// but the trailing edge still have complete topological information in the local partition.
//        incompleteTrailingEdgeTopology = true;
//      }

    }
    //If we only find one node behind the tip, we do the following as a rough compensation.
    if( trailingNodes.size() == 1 && !threeNodesPinched )
      fFaceA[i] *= 2.0;


  }
#else

  R1Tensor fFaceB[2];

  if( ppFlowPressure != nullptr )
  {
    localIndex const faceID = tipFaces[0];

    const R1Tensor N[2] = { faceManager.FaceNormal( nodeManager, tipFaces[0] ),
                            faceManager.FaceNormal( nodeManager, tipFaces[1] )};
    R1Tensor Nbar = N[0];
    Nbar -= N[1];
    Nbar.Normalize();

    fFaceB[0] = Nbar;
    fFaceB[0] *= (*ppFlowPressure)[faceID] * faceArea[faceID] * 2 / facesToNodes[faceID].size();

//    std::cout<<"     facePressure["<<faceID<<"] = "<<(*ppFlowPressure)[faceID]<<std::endl;

    realT ElasticModulus = 0.0;
    Array1dT< std::pair< CellElementSubRegion *, localIndex > > toElements = faceManager.m_toElementsRelation[faceID];
    for( Array1dT< std::pair< CellElementSubRegion *, localIndex > >::iterator iterRegion=toElements.begin() ;
         iterRegion!=toElements.end() ; ++iterRegion )
    {
      ElementRegionT const & elemRegion = *(iterRegion->first);
      localIndex const element = iterRegion->second;

      const localIndex paramIndex = elemRegion.m_mat->NumParameterIndex0() > 1 ? element : 0;
      ElasticModulus += elemRegion.m_mat->ParameterData( paramIndex )->E;
    }
    ElasticModulus /= toElements.size();

    fFaceB[0] *= ElasticModulus;
    fFaceB[1] = fFaceA[0];
    fFaceB[1] *= -1;
  }


#endif

  R1Tensor tipForce;
  tipForce[0] = Dot( fNodeO, vecTipNorm ) + Dot( fFaceA[0], vecTipNorm ) / 2.0 - Dot( fFaceA[1], vecTipNorm ) /2.0;
  tipForce[1] = Dot( fNodeO, vecTip ) + Dot( fFaceA[0], vecTip ) / 2.0 - Dot( fFaceA[1], vecTip ) /2.0;
  tipForce[2] = Dot( fNodeO, vecEdge ) + Dot( fFaceA[0], vecEdge ) / 2.0 - Dot( fFaceA[1], vecEdge ) /2.0;

//  tipForce[0] = Dot(fNodeO, vecTipNorm) + ( Dot(fFaceA[0], vecTipNorm) - Dot(fFaceA[1], vecTipNorm) ) /2.0;
//  tipForce[1] = Dot(fNodeO, vecTip)     + ( Dot(fFaceA[0], vecTip)     - Dot(fFaceA[1], vecTip) ) /2.0;
//  tipForce[2] = Dot(fNodeO, vecEdge)    + ( Dot(fFaceA[0], vecEdge)    - Dot(fFaceA[1], vecEdge) ) /2.0;

//  std::cout<<"     edgeID = "<<edgeID<<std::endl;
//  std::cout<<"     Dot(fNodeO, vecTipNorm) = "<<Dot(fNodeO, vecTipNorm)<<std::endl;
//  std::cout<<"     Dot(fFaceA[0], vecTipNorm) = "<<( Dot(fFaceA[0], vecTipNorm) - Dot(fFaceA[1], vecTipNorm) )
// /2.0<<std::endl;
//  std::cout<<"     Dot(fFaceB[0], vecTipNorm) = "<<( Dot(fFaceB[0], vecTipNorm) - Dot(fFaceB[1], vecTipNorm) )
// /2.0<<std::endl;

  R1Tensor tipDisplacement, tipOpening, tipFaceDisplacement[2];

  if( !threeNodesPinched )
  {
    for( localIndex i=0 ; i<2 ; ++i )
    {
      localIndex faceID = tipFaces[i];
      tipFaceDisplacement[i] = 0.0;

      for( localIndex_array::iterator j = facesToNodes[faceID].begin() ;
           j!=facesToNodes[faceID].end() ; ++j )
      {
        localIndex iNd = *j;
        if( iNd != edgesToNodes( edgeID, 0 ) && iNd != edgesToNodes( edgeID, 1 ) )
        {
          tipFaceDisplacement[i] += (*nodeManager.m_displacement)[iNd];
        }
      }

      tipFaceDisplacement[i] /= (facesToNodes[faceID].size() - 2);
    }
    tipDisplacement = tipFaceDisplacement[1];
    tipDisplacement -= tipFaceDisplacement[0];
  }
  else
  {
    tipDisplacement = (*nodeManager.m_displacement)[openNodeID[1]];
    tipDisplacement -= (*nodeManager.m_displacement)[openNodeID[0]];


// This is Randy's old algorithm to handle the three-pinched situation.  I am keeping it because it's quite clever in
// how it matched the node pairs across the two fracture faces.
//    const localIndex numNodes = facesToNodes[tipFaces[0]].size();
//
//    for( localIndex a=0 ; a<numNodes ; ++a )
//    {
//      const localIndex aa = a == 0 ? a : numNodes - a;
//
//      const localIndex nodeID0 = facesToNodes[tipFaces[0]][a];
//      const localIndex nodeID1 = facesToNodes[tipFaces[1]][aa];
//
//      if ( nodeID0 != edgesToNodes(edgeID,0) && nodeID0 != edgesToNodes(edgeID,1) &&
//           nodeID1 != edgesToNodes(edgeID,0) && nodeID1 != edgesToNodes(edgeID,1) )
//      {
//        tipDisplacement = (*nodeManager.m_displacement)[nodeID1];
//        tipDisplacement -= (*nodeManager.m_displacement)[nodeID0];
//
//        if( tipDisplacement.L2_Norm() > maxTipDisplacement.L2_Norm() )
//        {
//          maxTipDisplacement = tipDisplacement;
//        }
//      }
//    }
  }


  tipOpening[0] = Dot( tipDisplacement, vecTipNorm );
  tipOpening[1] = Dot( tipDisplacement, vecTip );
  tipOpening[2] = Dot( tipDisplacement, vecEdge );

  //if (Dot(Cross(vecTip, vecTipNorm),vecEdge) < 0.0) tipOpening[1] *= -1.0;



  realT tipArea;
  tipArea = faceManager.SurfaceArea( nodeManager, faceA, true );
  if( facesToNodes[faceA].size() == 3 )
  {
    tipArea *= 2.0;
  }

  if( !incompleteTrailingEdgeTopology && !m_displacementBasedSIF && tipOpening[0] * tipForce[0] > 0.0 )
  {
    SIF_I[edgeID] = pow( fabs( tipForce[0] * tipOpening[0] / 2.0 / tipArea ), 0.5 );
    SIF_II[edgeID] = pow( fabs( tipForce[1] * tipOpening[1] / 2.0 / tipArea ), 0.5 );
    SIF_III[edgeID] = pow( fabs( tipForce[2] * tipOpening[2] / 2.0 / tipArea ), 0.5 );

    if( tipOpening[0] < 0 )
    {
      // We don't need this for the case of incomplete trailing edge topology.  Sign in that case should be taken care
      // of automatically because there is no sqrt involved.
      SIF_I( edgeID ) *= -1.0;
    }

  }
  else
  {
    // Opening-based SIF, based on
    // Equation 1 in Fu et al. 2012, DOI::10.1016/j.engfracmech.2012.04.010
    realT r = tipArea / edgeManager.EdgeLength( nodeManager, edgeID );
    // if (facesToNodes[faceA].size() == 3) r *= 2.0;
    // tipArea is already twice the triangle area so we don't need this.
    SIF_I[edgeID] = tipOpening[0] / 2.0 * GdivBeta / pow( r/6.28, 0.5 );
    SIF_II[edgeID] = 0.0;  // SIF is not accurate in this scenario anyway.  Let's not worry about turning.
    SIF_III[edgeID] = 0.0;
  }


  if( tipForce[1] < 0.0 )
  {
    SIF_II[edgeID] *= -1.0;
  }

//  std::cout<<"          edge, tipOpening, tipForce = "<<edgeID<<", "<<tipOpening[0]<<", "<<tipForce[0]<<std::endl;
//  std::cout<<"     SIF_I[edgeID] =  = "<<SIF_I[edgeID]<<std::endl;

  if( SIF_I[edgeID] > 0.0 )
  {
    rval = pow( SIF_I[edgeID]*SIF_I[edgeID]+SIF_II[edgeID]*SIF_II[edgeID]+SIF_III[edgeID]*SIF_III[edgeID], 0.5 );
  }
  else
  {
    rval = -1.0;
  }

  return rval;

#endif
  return 0;
}

void SurfaceGenerator::MarkRuptureFaceFromEdge ( const localIndex edgeID,
                                                 localIndex& trailFaceID,
                                                 NodeManager & nodeManager,
                                                 EdgeManager & edgeManager,
                                                 FaceManager & faceManager,
                                                 R1Tensor& vecTipNorm,
                                                 R1Tensor& vecTip,
                                                 ModifiedObjectLists& modifiedObjects,
                                                 const int edgeMode )
{
#if 0
  integer_array& ruptureStateString = faceManager.getReference<integer_array>( "ruptureState" );
  real64_array& stressNOnFace = faceManager.getReference<real64_array>( "stressNOnFace" );
  real64_array& SIFonFace = faceManager.getReference<real64_array>( "SIFonFace" );
  array1d<R1Tensor> & KIC = faceManager.getReference<r1_array>( "K_IC" );
  real64_array& SIF_I = edgeManager.getReference<real64_array>( "SIF_I" );
  real64_array& SIF_II = edgeManager.getReference<real64_array>( "SIF_II" );
  localIndex_array& primaryCandidateFace = faceManager.getReference<localIndex_array>( "primaryCandidateFace" );
  integer_array& isFaceSeparable = faceManager.getReference<integer_array>( "isSeparable" );
  integer_array* dfnIndexMap = faceManager.getReferencePointer<integer_array>( "DFN_Index" );

  integer_array eligibleFaces;
  realT lowestSIF = std::numeric_limits<realT>::max();
  realT highestSIF = std::numeric_limits<realT>::min();
  realT lowestStress = std::numeric_limits<realT>::max();
  realT lowestScore = std::numeric_limits<realT>::max();
  realT highestScore = std::numeric_limits<realT>::min();
  realT secondScore = std::numeric_limits<realT>::min();
  localIndex faceWithHighestScore = std::numeric_limits<localIndex>::max();
  localIndex faceWithSecondScore = std::numeric_limits<localIndex>::max();

  R1Tensor vecEdge, edgeCenter;
  edgeManager.EdgeVector( nodeManager, edgeID, vecEdge );
  edgeManager.EdgeCenter( nodeManager, edgeID, edgeCenter );
  vecEdge.Normalize();

  for( set<localIndex>::const_iterator a=edgeManager.m_toFacesRelation[edgeID].begin() ;
       a!=edgeManager.m_toFacesRelation[edgeID].end() ; ++a )
  {
    localIndex iface = *a;

    if( faceManager.m_toElementsRelation[iface].size() == 2  &&
        faceManager.m_isExternal[iface] != 1 &&
        CheckOrphanElement( faceManager, iface ) == 0 &&
        isFaceSeparable[iface] == 1 )
    {
      R1Tensor fc, fn, vecFace;
      faceManager.FaceCenterAndNormal( nodeManager, iface, fc, fn );
      faceManager.InFaceVectorNormalToEdge( nodeManager,
                                            edgeManager,
                                            iface, edgeID,
                                            vecFace );
      if( Dot( vecTip, vecFace ) > cos( m_maxTurnAngle ))
      {
        eligibleFaces.push_back( iface );
        realT thetaFace = acos( Dot( vecTip, vecFace )*0.999999 );  // We multiply this by 0.9999999 to avoid an
                                                                    // exception caused by acos a number slightly larger
                                                                    // than 1.

        if( Dot( Cross( vecTip, vecFace ), vecEdge ) < 0.0 )
        {
          thetaFace *= -1.0;
        }

        SIFonFace[iface] = cos( thetaFace / 2.0 ) *
                           ( SIF_I[edgeID] * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II[edgeID] * sin( thetaFace ) );

        R1Tensor direction( fc );
        direction -= edgeCenter;
        direction.Normalize();
        realT faceToughness = std::fabs( Dot( direction, KIC[iface] ));

        highestSIF = std::max( highestSIF, SIFonFace[iface]/faceToughness );
        lowestSIF = std::min( lowestSIF, SIFonFace[iface]/faceToughness );
        lowestStress = std::min( lowestStress, stressNOnFace[iface] );

      }
    }
  }


  integer_array pickedFaces;
  if( eligibleFaces.size() >=1 )
  {
    realT lengthscale = edgeManager.EdgeLength( nodeManager, edgeID );

    for( localIndex i = 0 ; i < eligibleFaces.size() ; ++i )
    {
      localIndex iface = eligibleFaces[i];
      R1Tensor fc, direction( edgeCenter );
      faceManager.FaceCenter( nodeManager, iface, fc );
      direction -= fc;
      direction.Normalize();
      realT faceToughness = std::fabs( Dot( direction, KIC[iface] ));

      realT splitabilityScore = SIFonFace[iface] - lowestSIF * faceToughness + (stressNOnFace[iface] - lowestStress) * sqrt( lengthscale );
      lowestScore = std::min( lowestScore, splitabilityScore );

      if( faceWithHighestScore == std::numeric_limits<localIndex>::max())
      {
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if( splitabilityScore > highestScore )
      {
        faceWithSecondScore = faceWithHighestScore;
        secondScore = highestScore;
        faceWithHighestScore = iface;
        highestScore = splitabilityScore;
      }
      else if( splitabilityScore > secondScore )
      {
        faceWithSecondScore = iface;
        secondScore = splitabilityScore;
      }
    }

    pickedFaces.push_back( faceWithHighestScore );

    if( eligibleFaces.size() >= 3 && (highestScore - secondScore) < 0.1 * (highestScore - lowestScore))
    {
      pickedFaces.push_back( faceWithSecondScore );
    }

  }

  for( localIndex i = 0 ; i < pickedFaces.size() ; ++i )
  {
    localIndex pickedFace = pickedFaces[i];

    if( highestSIF > 1.0 && edgeMode == 1 && i == 0 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureStateString[pickedFace] = 1;
      if( !m_dfnPrefix.empty())
        (*dfnIndexMap)[pickedFace] = (*dfnIndexMap)[trailFaceID];
      modifiedObjects.modifiedFaces.insert( pickedFace );
    }
    else if( highestSIF > 1.0 && edgeMode == 1 && i == 1 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureStateString[pickedFace] = -1;
      if( !m_dfnPrefix.empty())
        (*dfnIndexMap)[pickedFace] = (*dfnIndexMap)[trailFaceID];
      modifiedObjects.modifiedFaces.insert( pickedFace );
      primaryCandidateFace[pickedFace] = faceWithHighestScore;
    }


    // We didn't really need to do this unless the criterion above has been satisfied.
    // We are calculating this regardless the criterion for debugging purpose.
    if( m_markExtendedLayer == 1 && highestSIF > 1.0 && edgeMode == 1 )
    {
      // Next we mark the faces that are 1) connected to this face, and 2) attached to one node of the edge (implicitly
      // satisfied), and 3) almost co-plane with this face
      for( localIndex_array::iterator j = facesToEdges[pickedFace].begin() ;
           j!=facesToEdges[pickedFace].end() ; ++j )
      {
        localIndex iedge = *j;
        if( iedge != edgeID )
        {
          for( set<localIndex>::const_iterator a=edgeManager.m_toFacesRelation[iedge].begin() ;
               a!=edgeManager.m_toFacesRelation[iedge].end() ; ++a )
          {
            localIndex iface = *a;
            if( iface != pickedFace && isFaceSeparable[iface] == 1 && faceManager.m_isExternal[iface] != 1 &&
                ( faceManager.IsNodeOnFace( iface, edgesToNodes[edgeID][0] ) ||
                  faceManager.IsNodeOnFace( iface, edgesToNodes[edgeID][1] )))
            {
              R1Tensor fc, fn, vecFace, fn0, fc0, ptPrj;
              realT nDist, uDist, segmentLength;
              faceManager.FaceCenterAndNormal( nodeManager, iface, fc, fn );
              faceManager.FaceCenterAndNormal( nodeManager, pickedFace, fc0, fn0 );
              faceManager.InFaceVectorNormalToEdge( nodeManager,
                                                    edgeManager,
                                                    iface, edgeID,
                                                    vecFace );
              GeometryUtilities::ProjectPointToLineSegment( (*nodeManager.m_refposition)[edgesToNodes( edgeID, 0 )],
                                                            (*nodeManager.m_refposition)[edgesToNodes( edgeID, 1 )],
                                                            fc,
                                                            nDist, uDist, segmentLength, ptPrj );

              // thetaFace does not strictly speaking apply to this face since the tip edge is not a edge of this face.
              // We calculate it as if this face is coplane with the master face
              realT thetaFace = acos( Dot( vecTip, vecFace )*0.999999 );
              if( Dot( Cross( vecTip, vecFace ), vecEdge ) < 0.0 )
              {
                thetaFace *= -1.0;
              }

              if( Dot( vecTip, vecFace ) > cos( m_maxTurnAngle ) &&
                  uDist / segmentLength > -m_faceToEdgeProjectionTol &&
                  uDist / segmentLength < 1 + m_faceToEdgeProjectionTol &&
                  fabs( Dot( vecEdge, fn )) < m_faceToEdgeCoplaneTol &&  // this face is kind of parallel to the tip
                                                                         // edge.
                  fabs( Dot( fn0, fn )) > 1 - m_faceToFaceCoplaneTol )  // co-plane
              {
                // Calculate SIFonFace is not really necessary but we keep it here for now for debugging.
                // Marking of the extended layer is purely based on geometrical and topological criteria.
                SIFonFace[iface] = std::max( SIFonFace[iface],
                                             SIF_I[edgeID] * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II[edgeID] * sin( thetaFace ) );

                if( highestSIF > 1.0 && edgeMode == 1 )
                {
                  ruptureStateString[iface] = ruptureStateString[pickedFace];
                  if( !m_dfnPrefix.empty())
                    (*dfnIndexMap)[pickedFace] = (*dfnIndexMap)[trailFaceID];
                  modifiedObjects.modifiedFaces.insert( iface );
                  primaryCandidateFace[iface] = primaryCandidateFace[pickedFace];
                }
              }
            }
          }
        }
      }
    }
  }
#endif
}

void SurfaceGenerator::PostUpdateRuptureStates( NodeManager & nodeManager,
                                                EdgeManager & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager & elementManager,
                                                array1d<std::set<localIndex> >& nodesToRupturedFaces,
                                                array1d<std::set<localIndex> >& edgesToRupturedFaces )
{
  arrayView1d<localIndex_array> const & facesToNodes = faceManager.nodeList();
  arrayView1d<localIndex_array> const & facesToEdges = faceManager.edgeList();
  nodesToRupturedFaces.resize( nodeManager.size() );
  edgesToRupturedFaces.resize( edgeManager.size() );

  arrayView1d<integer>& faceRuptureState = faceManager.getReference<integer_array>( "ruptureState" );
  arrayView1d<localIndex const> const &
  faceParentIndex = faceManager.getReference<localIndex_array>( ObjectManagerBase::viewKeyStruct::parentIndexString );


  // assign the values of the nodeToRupturedFaces and edgeToRupturedFaces arrays.
  for( localIndex kf=0 ; kf<faceManager.size() ; ++kf )
  {
    if( faceRuptureState[kf] >0  )
    {
      localIndex const faceIndices[2] = { kf, faceParentIndex[kf] };
      int const n = faceParentIndex[kf]==-1 ? 1 : 2;
      localIndex const faceIndex = faceParentIndex[kf]==-1 ? kf : faceParentIndex[kf];

      for( int i=0 ; i<n ; ++i )
      {
        for( localIndex a=0 ; a<facesToNodes[kf].size() ; ++a )
        {
          const localIndex nodeIndex = facesToNodes[kf][a];
          nodesToRupturedFaces[nodeIndex].insert( faceIndex );
        }

        for( localIndex a=0 ; a<facesToEdges[kf].size() ; ++a )
        {
          const localIndex edgeIndex = facesToEdges[kf][a];
          edgesToRupturedFaces[edgeIndex].insert( faceIndex );
        }
      }
    }
  }
}

int SurfaceGenerator::CheckEdgeSplitability( const localIndex edgeID,
                                             NodeManager & nodeManager,
                                             FaceManager & faceManager,
                                             EdgeManager & edgeManager,
                                             const bool prefrac )
{
  //     Return value = -1, this edge won't split for sure, don't do any more work;
  //                  = 0, edge is along a tip, but the fracture connected to it is not saturated yet.  We will only
  // calculate SIF but will not perform splitting.
  //                  = 1, edge is along a tip and the adjacent fracture is saturated, more work to be done; or this is
  // a dry simulation
  //                  = 2, this is a singular edge, we need split it.
  //                  = 3, this is an eligible kink, we need to process it as a kink

  arrayView1d< set<localIndex> > const & edgesToFaces = edgeManager.faceList();

  int isSplitable = -1;

  if( edgeManager.m_isExternal[edgeID] == 0 )
  {
    return isSplitable;
  }

  // We first count the external faces connected to this edge;
  int nExternalFaces = 0;
  localIndex_array faceInvolved;
  for( auto iface : edgesToFaces[edgeID] )
  {
    if( faceManager.m_isExternal[iface] == 1 )
    {
      nExternalFaces++;
      faceInvolved.push_back( iface );
    }
  }

  if( nExternalFaces%2 == 1 )
  {
    //    char msg[200];
    //    sprintf(msg, "Error! Edge %d has an odd number of external faces.", int(edgeID));
    //    GEOS_ERROR(msg);
    //    std::cout << "Error! Edge " << int(edgeID) << " has an odd number of external faces. "
    //        << (*nodeManager.m_refposition)[edgesToNodes[edgeID][0]][0] << " ,"
    //        << (*nodeManager.m_refposition)[edgesToNodes[edgeID][0]][1] << " ,"
    //        << (*nodeManager.m_refposition)[edgesToNodes[edgeID][0]][2] << " ,";
    //    isSplitable = -1;
    return(isSplitable);
  }

  if( nExternalFaces >= 4 )
  {
    isSplitable = 2;
    return (isSplitable);
  }

  if( nExternalFaces == 2 )
  {
    isSplitable = 1;
  }

  return (isSplitable);
}

realT SurfaceGenerator::MinimumToughnessOnEdge( const localIndex edgeID,
                                                const NodeManager & nodeManager,
                                                EdgeManager & edgeManager,
                                                FaceManager & faceManager )
{
  realT val = std::numeric_limits<realT>::max();

  R1Tensor edgeCenter( 0.0 );
  R1Tensor faceCenter( 0.0 );

  arrayView1d<R1Tensor> const & X = nodeManager.referencePosition();
  edgeCenter += X[edgeManager.nodeList( edgeID, 0 )];
  edgeCenter += X[edgeManager.nodeList( edgeID, 1 )];
  edgeCenter *= 0.5;

  arrayView1d<R1Tensor> & KIC = faceManager.getReference<r1_array>( "K_IC" );
  arrayView1d<localIndex_array> const & faceToNodes = faceManager.nodeList();
  for( auto iface : edgeManager.faceList()[edgeID] )
  {
    faceCenter = 0.0;
    for( localIndex a=0 ; a<faceToNodes.size( 1 ) ; ++a )
    {
      faceCenter += X[ faceToNodes[iface][a] ];
    }
    faceCenter /= faceToNodes.size( 1 );

    R1Tensor direction( faceCenter );
    direction -= edgeCenter;
    direction.Normalize();
    val = std::min( val, fabs( Dot( KIC[iface], direction ) ) );
  }

  return val;
}

int SurfaceGenerator::CheckNodeSplitability( const localIndex nodeID,
                                             NodeManager & nodeManager,
                                             FaceManager & faceManager,
                                             EdgeManager & edgeManager,
                                             const bool prefrac )
{
  return (1);
  //
  //  if (prefrac)
  //  {
  //    return (1);
  //  }
  //  else if (m_failCriterion == 0)
  //  {
  //    return (1);
  //  }
  //  else if (nodeManager.m_isExternal[nodeID] == 1)
  //  {
  //    return (1);
  //  }
  //  else
  //  {
  //    return (0);
  //  }
}

void SurfaceGenerator::AssignNewGlobalIndicesSerial( ObjectManagerBase & object,
                                                     std::set<localIndex> const & indexList )
{
  // in serial, we can simply loop over the indexList and assign new global indices based on
  // the value of the m_maxGlobalIndex + 1;
  for( localIndex const newLocalIndex : indexList )
  {
    globalIndex const newGlobalIndex = object.m_maxGlobalIndex + 1;

    object.m_localToGlobalMap[newLocalIndex] = newGlobalIndex;
    object.m_globalToLocalMap[newGlobalIndex] = newLocalIndex;

    object.m_maxGlobalIndex = newGlobalIndex;
  }
}

void SurfaceGenerator::
AssignNewGlobalIndicesSerial( ElementRegionManager & elementManager,
                              std::map< std::pair<localIndex,localIndex>, std::set<localIndex> > const & newElems )
{
  // in serial, we can simply iterate over the entries in newElems and assign new global indices based on
  // the value of the m_maxGlobalIndex + 1 for the ElementRegionManager.

  // loop over entries of newElems, which gives elementRegion/subRegion local indices
  for( auto const & iter : newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set<localIndex> const & indexList = iter.second;

    ElementSubRegionBase * const subRegion = elementManager.GetRegion(er)->GetSubRegion(esr);

    // loop over the new elems in the subRegion
    for( localIndex const newLocalIndex : indexList )
    {
      globalIndex const newGlobalIndex = elementManager.m_maxGlobalIndex + 1;

      subRegion->m_localToGlobalMap[newLocalIndex] = newGlobalIndex;
      subRegion->m_globalToLocalMap[newGlobalIndex] = newLocalIndex;

      elementManager.m_maxGlobalIndex = newGlobalIndex;
    }
  }
}


REGISTER_CATALOG_ENTRY( SolverBase,
                        SurfaceGenerator,
                        std::string const &, dataRepository::ManagedGroup * const )

} /* namespace geosx */
