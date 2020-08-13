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

#include "WellElementSubRegion.hpp"

#include "mesh/MeshLevel.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "LvArray/src/output.hpp"

namespace geosx
{
WellElementSubRegion::WellElementSubRegion( string const & name, Group * const parent ) :
  ElementSubRegionBase( name, parent ),
  m_wellControlsName( "" ),
  m_toNodesRelation(),
  m_topWellElementIndex( -1 ),
  m_perforationData( groupKeyStruct::perforationDataString, this ),
  m_topRank( -1 ),
  m_searchDepth( 10 )
{
  registerWrapper( viewKeyStruct::wellControlsString, &m_wellControlsName );
  registerWrapper( viewKeyStruct::wellNodeListString, &m_toNodesRelation );
  registerWrapper( viewKeyStruct::nextWellElementIndexString,
                   &m_nextWellElementIndex );
  registerWrapper( viewKeyStruct::nextWellElementIndexGlobalString,
                   &m_nextWellElementIndexGlobal );
  registerWrapper( viewKeyStruct::topWellElementIndexString,
                   &m_topWellElementIndex );
  registerWrapper( viewKeyStruct::topRankString, &m_topRank );
  registerWrapper( viewKeyStruct::radiusString, &m_radius );

  RegisterGroup( groupKeyStruct::perforationDataString, &m_perforationData );

  this->setNumNodesPerElement( 2 );
  this->setNumFacesPerElement( 0 );
  m_toNodesRelation.resizeDimension< 1 >( this->numNodesPerElement() );
  m_elementTypeString = "BEAM";
}

WellElementSubRegion::~WellElementSubRegion()
{}

void
WellElementSubRegion::setupRelatedObjectsInRelations( MeshLevel const * const mesh )
{
  m_toNodesRelation.SetRelatedObject( mesh->getNodeManager() );
}

namespace
{
/**
 * @brief Now that the well elements are assigned, collect the nodes and tag the boundary nodes between ranks
          The function WellElementSubRegion::AssignUnownedElements must have been called before this function
 * @param[in] wellGeometry the InternalWellGenerator containing the global well topology
 * @param[in] localElems set of local well elems. At this point all the well elems have been assigned
 * @param[out] localNodes set of local well nodes (includes boundary nodes)
 * @param[out] boundaryNodes set of local well nodes that are at the boundary between this rank
               and another rank
 */
void
CollectLocalAndBoundaryNodes( InternalWellGenerator const & wellGeometry,
                              SortedArray< globalIndex > const & localElems,
                              SortedArray< globalIndex > & localNodes,
                              SortedArray< globalIndex > & boundaryNodes )
{
  // get the well connectivity
  arrayView1d< globalIndex const > const & nextElemIdGlobal =
    wellGeometry.GetNextElemIndex();
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal =
    wellGeometry.GetPrevElemIndices();
  arrayView2d< globalIndex const > const & elemToNodesGlobal =
    wellGeometry.GetElemToNodesMap();

  // loop over the local elements and collect the local and boundary nodes
  for( globalIndex currGlobal : localElems )
  {
    // if the element is local, its two nodes are also local
    globalIndex const inodeTopGlobal =
      elemToNodesGlobal[currGlobal][InternalWellGenerator::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal =
      elemToNodesGlobal[currGlobal][InternalWellGenerator::NodeLocation::BOTTOM];
    localNodes.insert( inodeTopGlobal );
    localNodes.insert( inodeBottomGlobal );

    localIndex const nextGlobal = LvArray::integerConversion< localIndex >(
      nextElemIdGlobal[LvArray::integerConversion< localIndex >( currGlobal )] );

    // if the next well elem is not local, add the node in between curr and next to boundaryNodes
    if( nextGlobal >= 0 && !localElems.contains( nextGlobal ) )
    {
      boundaryNodes.insert( inodeTopGlobal );
    }

    // if the prev well elem is not local, add the node in between curr and prev to boundaryNodes (relevant for
    // branches)
    for( localIndex iwelem = 0; iwelem < prevElemIdsGlobal[currGlobal].size();
         ++iwelem )
    {
      globalIndex const prevGlobal = prevElemIdsGlobal[currGlobal][iwelem];
      if( prevGlobal >= 0 && !localElems.contains( prevGlobal ) )
      {
        boundaryNodes.insert( inodeBottomGlobal );
      }
    }
  }
}

/**
 * @brief Check if "location" is contained in reservoir element ei
 * @param[in] subRegion the subRegion of reservoir element ei
 * @param[in] ei the index of the reservoir element
 * @return true if "location" is contained in reservoir element ei, false otherwise
 */
bool
IsPointInsideElement( NodeManager const * const nodeManager,
                      R1Tensor const & location,
                      CellBlock const * subRegion,
                      localIndex ei )
{
  bool isInsideElement = false;

  array1d< array1d< localIndex > > faceNodes( subRegion->numFacesPerElement() );

  // collect the faces for this element
  for( localIndex kf = 0; kf < subRegion->numFacesPerElement(); ++kf )
  {
    subRegion->GetFaceNodes( ei, kf, faceNodes[kf] );
  }

  // if the point is in the element, save the indices and stop the search
  if( computationalGeometry::IsPointInsidePolyhedron(
        nodeManager->referencePosition(),
        faceNodes,
        location ) )
  {
    isInsideElement = true;
  }
  return isInsideElement;
}

/**
 * @brief Collect the nodes of reservoir element ei
 * @param[in] subRegion the subRegion of reservoir element ei
 * @param[in] ei the index of the reservoir element
 * @param[inout] nodes the nodes that have already been visited
 */
void
CollectElementNodes( CellBlock const * subRegion,
                     localIndex ei,
                     SortedArray< localIndex > & nodes )
{
  // get all the nodes belonging to this element
  for( localIndex a = 0; a < subRegion->numNodesPerElement(); ++a )
  {
    localIndex const inode = subRegion->nodeList( ei, a );

    // if not already visited, store the newly found node
    if( !nodes.contains( inode ) )
    {
      nodes.insert( inode );
    }
  }
}

/**
 * @brief Search the reservoir elements that can be accessed from the set "nodes".
          Stop if a reservoir element containing the perforation is found.
          If not, enlarge the set "nodes"
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[inout] nodes the nodes that have already been visited
 * @param[inout] elements the reservoir elements that have already been visited
 * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
 * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 */
bool
VisitNeighborElements( MeshLevel const & mesh,
                       R1Tensor const & location,
                       SortedArray< localIndex > & nodes,
                       SortedArray< globalIndex > & elements,
                       localIndex & erMatched,
                       localIndex & esrMatched,
                       localIndex & eiMatched )
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager = mesh.getNodeManager();

  ArrayOfArraysView< localIndex const > const & toElementRegionList =
    nodeManager->elementRegionList();
  ArrayOfArraysView< localIndex const > const & toElementSubRegionList =
    nodeManager->elementSubRegionList();
  ArrayOfArraysView< localIndex const > const & toElementList =
    nodeManager->elementList();

  bool matched = false;

  // In this function, we loop over all the nodes that we have collected so far.
  // For each node, we grab all the reservoir elements that contain the node
  //     For reservoir element that has not been visited yet, we check if it contains "location"
  //           If the reservoir element contains "location" we can stop the search
  //           If the reservoir element does not contain "location", we collect its nodes

  // we will enlarge the set of nodes in the loop below
  // to do this we have to create a new set, "currNodes"
  // that contains only the nodes that have already been visited
  // the newly added nodes will be added to "nodes"
  SortedArray< localIndex > currNodes = nodes;

  // for all the nodes already visited
  for( localIndex currNode : currNodes )
  {
    // collect the elements that have not been visited yet
    for( localIndex b = 0; b < toElementRegionList.sizeOfArray( currNode ); ++b )
    {
      localIndex const er = toElementRegionList[currNode][b];
      localIndex const esr = toElementSubRegionList[currNode][b];
      localIndex const eiLocal = toElementList[currNode][b];

      CellElementRegion const * const region =
        dataRepository::Group::group_cast< CellElementRegion const * >(
          elemManager->GetRegion( er ) );
      CellBlock const * const subRegion =
        dataRepository::Group::group_cast< CellElementSubRegion const * >(
          region->GetSubRegion( esr ) );
      globalIndex const eiGlobal = subRegion->localToGlobalMap()[eiLocal];

      // if this element has not been visited yet, save it
      if( !elements.contains( eiGlobal ) )
      {
        elements.insert( eiGlobal );

        // perform the test to see if the point is in this reservoir element
        // if the point is in the resevoir element, save the indices and stop the search
        if( IsPointInsideElement( nodeManager, location, subRegion, eiLocal ) )
        {
          erMatched = er;
          esrMatched = esr;
          eiMatched = eiLocal;
          matched = true;
          break;
        }
        // otherwise add the nodes of this element to the set of new nodes to visit
        else
        {
          CollectElementNodes( subRegion, eiLocal, nodes );
        }
      }
    }

    if( matched )
    {
      break;
    }
  }

  // if not matched, insert the new nodes
  return matched;
}

/**
 * @brief Search for the reservoir element that is the *closest* from the center of well element.
          Note that this reservoir element does not necessarily contain the center of the well element.
          This "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
          contains the well element.
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[inout] erInit the region index of the reservoir element from which we start the search
 * @param[inout] esrInit the subregion index of the reservoir element from which we start the search
 * @param[inout] eiInit the element index of the reservoir element from which we start the search
 */
void
InitializeLocalSearch( MeshLevel const & mesh,
                       R1Tensor const & location,
                       localIndex & erInit,
                       localIndex & esrInit,
                       localIndex & eiInit )
{
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > resElemCenter =
    mesh.getElemManager()
      ->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >(
        ElementSubRegionBase::viewKeyStruct::elementCenterString );
  // to initialize the local search for the reservoir element that contains "location",
  // we find the reservoir element that minimizes the distance from "location" to the reservoir element center
  auto ret = minLocOverElemsInMesh(
    &mesh,
    [&]( localIndex const er, localIndex const esr, localIndex const ei ) -> real64 {
      R1Tensor v = location;
      v -= resElemCenter[er][esr][ei];
      return v.L2_Norm();
    } );

  // save the region, subregion and index of the reservoir element
  // note that this reservoir element does not necessarily contains "location"
  erInit = std::get< 0 >( ret.second );
  esrInit = std::get< 1 >( ret.second );
  eiInit = std::get< 2 >( ret.second );
}

/**
 * @brief Search for the reservoir element that contains the well element.
          To do that, loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
 * @param[in] meshLevel the mesh object (single level only)
 * @param[in] location the location of that we are trying to match with a reservoir element
 * @param[in] erInit the region index of the reservoir element from which we start the search
 * @param[in] esrInit the subregion index of the reservoir element from which we start the search
 * @param[in] eiInit the element index of the reservoir element from which we start the search
 * @param[inout] erMatched the region index of the reservoir element that contains "location", if any
 * @param[inout] esrMatched the subregion index of the reservoir element that contains "location", if any
 * @param[inout] eiMatched the element index of the reservoir element that contains "location", if any
 */
bool
SearchLocalElements( MeshLevel const & mesh,
                     R1Tensor const & location,
                     localIndex const & searchDepth,
                     localIndex const & erInit,
                     localIndex const & esrInit,
                     localIndex const & eiInit,
                     localIndex & erMatched,
                     localIndex & esrMatched,
                     localIndex & eiMatched )
{
  // search locally, starting from the location of the previous perforation
  // the assumption here is that perforations have been entered in order of depth
  bool resElemFound = false;

  CellElementRegion const * region =
    dataRepository::Group::group_cast< CellElementRegion const * >(
      mesh.getElemManager()->GetRegion( erInit ) );
  CellBlock const * subRegion =
    dataRepository::Group::group_cast< CellBlock const * >(
      region->GetSubRegion( esrInit ) );

  SortedArray< localIndex > nodes;
  SortedArray< globalIndex > elements;

  // here is how the search is done:
  //   1 - We check if "location" is within the "init" reservoir element defined by (erInit,esrMatched,eiMatched)
  //   2 - If yes, stop
  //     - If not, a) collect the nodes of the reservoir element defined by (erInit,esrMatched,eiMatched)
  //               b) use these nodes to grab the neighbors of (erInit,esrMatched,eiMatched)
  //               c) check if "location" is within the neighbors. If not, grab the neighbors of the neighbors, and so
  // on...

  // collect the nodes of the current element
  // they will be used to access the neighbors and check if they contain the perforation
  CollectElementNodes( subRegion, eiInit, nodes );

  // if no match is found, enlarge the neighborhood m_searchDepth'th times
  for( localIndex d = 0; d < searchDepth; ++d )
  {
    localIndex nNodes = nodes.size();

    // search the reservoir elements that can be accessed from the set "nodes"
    // stop if a reservoir element containing the perforation is found
    // if not, enlarge the set "nodes"
    resElemFound = VisitNeighborElements( mesh,
                                          location,
                                          nodes,
                                          elements,
                                          erMatched,
                                          esrMatched,
                                          eiMatched );
    if( resElemFound || nNodes == nodes.size() )
    {
      break;
    }
  }
  return resElemFound;
}

}  // namespace

void
WellElementSubRegion::Generate( MeshLevel & mesh,
                                InternalWellGenerator const & wellGeometry,
                                arrayView1d< integer > & elemStatusGlobal,
                                globalIndex nodeOffsetGlobal,
                                globalIndex elemOffsetGlobal )
{
  map< integer, SortedArray< globalIndex > > elemSetsByStatus;

  // convert elemStatus list into sets of indices
  for( localIndex iwelemGlobal = 0; iwelemGlobal < elemStatusGlobal.size();
       ++iwelemGlobal )
  {
    elemSetsByStatus[elemStatusGlobal[iwelemGlobal]].insert( iwelemGlobal );
  }

  // initialize the sets using the classification of well elems
  // localElems will be enlarged once boundary elements ownership is determined
  SortedArray< globalIndex > & localElems = elemSetsByStatus[WellElemStatus::LOCAL];
  SortedArray< globalIndex > & sharedElems =
    elemSetsByStatus[WellElemStatus::SHARED];
  SortedArray< globalIndex > & unownedElems =
    elemSetsByStatus[WellElemStatus::UNOWNED];

  // here we make sure that there are no shared elements
  // this is enforced in the InternalWellGenerator that currently merges two perforations
  // if they belong to the same well element. This is a temporary solution.
  // TODO: split the well elements that contain multiple perforations, so that no element is shared
  GEOSX_ERROR_IF( sharedElems.size() > 0,
                  "Well " << getName() << " contains shared well elements" );

  // In Steps 1 and 2 we determine the local objects on this rank (elems and nodes)
  // Once this is done, in Steps 3, 4, and 5, we update the nodeManager and wellElementSubRegion (size, maps)

  // 1) First assign the unowned elements to a rank
  // this is done in two steps

  // 1.a) First assign unowned elements in the reservoir based on location
  //      ie., if the center of the well element falls in the domain owned by rank k
  //      then the well element is assigned to rank k
  AssignUnownedElementsInReservoir( mesh,
                                    wellGeometry,
                                    unownedElems,
                                    localElems,
                                    elemStatusGlobal );
  // 1.b) Then we check that all the well elements have been assigned (and assigned once)
  //      This is needed because if the center of the well element falls on the boundary of
  //      a reservoir element, the assignment algorithm of 1.a) can assign the same well element
  //      to two ranks, or to no rank at all (which will break the solver).
  //      In this function we also check that the resulting well partitioning is valid, that is,
  //      we make sure that if two ranks are neighbors in the well, that are also neighbors in the
  //      reservoir mesh
  CheckPartitioningValidity( wellGeometry, localElems, elemStatusGlobal );

  SortedArray< globalIndex > localNodes;
  SortedArray< globalIndex > boundaryNodes;

  // 2) collect the local nodes and tag the boundary nodes using element info
  // now that all the elements have been assigned, we collected the local nodes
  // and tag the boundary nodes (i.e., the nodes in contact with both local and remote elems)
  CollectLocalAndBoundaryNodes( wellGeometry, localElems, localNodes, boundaryNodes );

  //DebugNodeManager( mesh );

  // 3) size update in the nodeManager
  // this is necessary to later use the node matching procedure
  // to place ghosts in DomainPartition::SetupCommunications
  UpdateNodeManagerSize( mesh,
                         wellGeometry,
                         localNodes,
                         boundaryNodes,
                         nodeOffsetGlobal );

  //DebugNodeManager( mesh );

  // 4) resize the well element subregion
  // and construct local to global, global to local, maps, etc
  ConstructSubRegionLocalElementMaps( mesh,
                                      wellGeometry,
                                      localElems,
                                      nodeOffsetGlobal,
                                      elemOffsetGlobal );

  //DebugWellElementSubRegions( elemStatusGlobal, elemOffsetGlobal );

  // 5) node-to-elem map update in the nodeManager
  // This map will be used by MeshLevel::GenerateAdjacencyLists
  // this assumes that the elemToNodes maps has been filled at Step 5)
  UpdateNodeManagerNodeToElementMap( mesh );
}

void
WellElementSubRegion::AssignUnownedElementsInReservoir(
  MeshLevel & mesh,
  InternalWellGenerator const & wellGeometry,
  SortedArray< globalIndex > const & unownedElems,
  SortedArray< globalIndex > & localElems,
  arrayView1d< integer > & elemStatusGlobal ) const
{
  // get the well and reservoir element coordinates
  arrayView1d< R1Tensor const > const & wellElemCoordsGlobal =
    wellGeometry.GetElemCoords();

  // assign the well elements based on location wrt the reservoir elements
  // if the center of the well element falls in the domain owned by rank k
  // then the well element is assigned to rank k
  for( globalIndex currGlobal : unownedElems )
  {
    R1Tensor const & location = wellElemCoordsGlobal[currGlobal];

    // this will contain the indices of the reservoir element
    // in which the center of the well element is located
    localIndex erMatched = -1;
    localIndex esrMatched = -1;
    localIndex eiMatched = -1;

    // this will contain the indices of the reservoir element
    // from which we are going to start the search
    localIndex erInit = -1;
    localIndex esrInit = -1;
    localIndex eiInit = -1;

    // Step 1: first, we search for the reservoir element that is the *closest* from the center of well element
    //         note that this reservoir element does not necessarily contain the center of the well element
    //         this "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
    //         contains the well element
    InitializeLocalSearch( mesh, location, erInit, esrInit, eiInit );

    // Step 2: then, search for the reservoir element that contains the well element
    //         to do that, we loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
    bool resElemFound = SearchLocalElements( mesh,
                                             location,
                                             m_searchDepth,
                                             erInit,
                                             esrInit,
                                             eiInit,
                                             erMatched,
                                             esrMatched,
                                             eiMatched );

    // if the element was found
    if( resElemFound )
    {
      // the well element is in the reservoir element (erMatched,esrMatched,eiMatched), so tag it as local
      localElems.insert( currGlobal );
      elemStatusGlobal[currGlobal] = WellElemStatus::LOCAL;
    }
  }
}

void
WellElementSubRegion::CheckPartitioningValidity(
  InternalWellGenerator const & wellGeometry,
  SortedArray< globalIndex > & localElems,
  arrayView1d< integer > & elemStatusGlobal ) const
{
  arrayView1d< arrayView1d< globalIndex const > const > const & prevElemIdsGlobal =
    wellGeometry.GetPrevElemIndices();

  // we are going to make sure that the partitioning is good,
  // well element per well element, starting from the bottom of the well
  for( globalIndex iwelemGlobal = wellGeometry.GetNumElements() - 1;
       iwelemGlobal >= 0;
       --iwelemGlobal )
  {
    // communicate the status of this element
    array1d< integer > thisElemStatusGlobal;
    MpiWrapper::allGather( elemStatusGlobal[iwelemGlobal], thisElemStatusGlobal );
    // group the ranks by well element status
    map< integer, SortedArray< globalIndex > > rankSetsByStatus;
    for( globalIndex irank = 0; irank < thisElemStatusGlobal.size(); ++irank )
    {
      rankSetsByStatus[thisElemStatusGlobal[irank]].insert( irank );
    }
    globalIndex const numLocalRanks =
      rankSetsByStatus[WellElemStatus::LOCAL].size();

    // in this case, this element has not been assigned
    //    => we assign it to the rank that owns
    //       the well element below iwelemGlobal (prevGlobal, already assigned and checked)
    if( numLocalRanks == 0 )
    {
      globalIndex const numBranches = prevElemIdsGlobal[iwelemGlobal].size();
      globalIndex const prevGlobal =
        prevElemIdsGlobal[iwelemGlobal][numBranches - 1];

      GEOSX_ERROR_IF( prevGlobal <= iwelemGlobal || prevGlobal < 0,
                      "Invalid partitioning in well " << getName() );

      if( elemStatusGlobal[prevGlobal] == WellElemStatus::LOCAL )
      {
        localElems.insert( iwelemGlobal );
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::LOCAL;
      }
      else
      {
        elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
      }
    }
    // in this case, everything is fine,
    // we just update the elemStatusGlobal array for all ranks
    else if( numLocalRanks == 1 )
    {
      for( globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL] )
      {
        if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != iownerRank )
        {
          elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
        }
      }
    }
    // in this case, this element has been assigned to more than rank
    //    => the smallest rank keeps it
    else  // (numLocalRanks > 1)
    {
      localIndex rankCount = 0;
      for( globalIndex iownerRank : rankSetsByStatus[WellElemStatus::LOCAL] )
      {
        if( rankCount == 0 )
        {
          // update the elemStatusGlobal array for all ranks
          if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != iownerRank )
          {
            elemStatusGlobal[iwelemGlobal] = WellElemStatus::REMOTE;
          }
        }
        else  // (rankCount > 0)
        {
          // remove the duplicate elements
          if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) == iownerRank )
          {
            localElems.remove( iwelemGlobal );
          }
        }
        rankCount++;
      }
    }

    // TODO: check neighbor rank
  }
}

void
WellElementSubRegion::UpdateNodeManagerSize(
  MeshLevel & mesh,
  InternalWellGenerator const & wellGeometry,
  SortedArray< globalIndex > const & localNodes,
  SortedArray< globalIndex > const & boundaryNodes,
  globalIndex nodeOffsetGlobal )
{
  // get the node manager to compute the total number of mesh nodes
  NodeManager * const nodeManager = mesh.getNodeManager();
  localIndex const numWellNodesLocal = localNodes.size();
  localIndex const oldNumNodesLocal = nodeManager->size();

  // resize nodeManager to account for the new well nodes and update the properties
  nodeManager->resize( oldNumNodesLocal + numWellNodesLocal );

  array1d< integer > & isDomainBoundary = nodeManager->getReference< integer_array >(
    m_ObjectManagerBaseViewKeys.domainBoundaryIndicator );

  arrayView1d< globalIndex > const & nodeLocalToGlobal =
    nodeManager->localToGlobalMap();

  arrayView1d< R1Tensor const > const & nodeCoordsGlobal =
    wellGeometry.GetNodeCoords();

  // local *well* index
  localIndex iwellNodeLocal = 0;
  // loop over global *well* indices

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();
  for( globalIndex iwellNodeGlobal : localNodes )
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal;

    // update node manager maps and position
    nodeLocalToGlobal[inodeLocal] =
      nodeOffsetGlobal + iwellNodeGlobal;  // global *nodeManager* index
    for( localIndex i = 0; i < 3; ++i )
    {
      X( inodeLocal, i ) = nodeCoordsGlobal[iwellNodeGlobal][i];
    }

    // mark the boundary nodes for ghosting in DomainPartition::SetupCommunications
    if( boundaryNodes.contains( iwellNodeGlobal ) )
    {
      isDomainBoundary[inodeLocal] = 1;
    }

    iwellNodeLocal++;
  }

  // now with update the relevant node indices in nodeManager->globalToLocalMap
  // this is to avoid a call to nodeManager->ConstructGlobalToLocalMap everytime we add a well
  for( iwellNodeLocal = 0; iwellNodeLocal < numWellNodesLocal; ++iwellNodeLocal )
  {
    // local *nodeManager* index
    localIndex const inodeLocal = oldNumNodesLocal + iwellNodeLocal;
    nodeManager->updateGlobalToLocalMap( inodeLocal );
  }
}

void
WellElementSubRegion::ConstructSubRegionLocalElementMaps(
  MeshLevel & mesh,
  InternalWellGenerator const & wellGeometry,
  SortedArray< globalIndex > const & localElems,
  globalIndex nodeOffsetGlobal,
  globalIndex elemOffsetGlobal )
{
  // get the well geometry
  arrayView1d< globalIndex const > const & nextElemIdGlobal =
    wellGeometry.GetNextElemIndex();
  arrayView1d< R1Tensor const > const & elemCoordsGlobal =
    wellGeometry.GetElemCoords();
  arrayView2d< globalIndex const > const & elemToNodesGlobal =
    wellGeometry.GetElemToNodesMap();
  arrayView1d< real64 const > const & elemVolumeGlobal =
    wellGeometry.GetElemVolume();

  NodeManager const * const nodeManager = mesh.getNodeManager();

  resize( localElems.size() );

  // create local elem numbering

  // local well elem ordering
  localIndex iwelemLocal = 0;
  // loop over global well elem indices
  for( globalIndex iwelemGlobal : localElems )
  {
    // create a global *elemManager* index
    m_localToGlobalMap[iwelemLocal++] = elemOffsetGlobal + iwelemGlobal;
  }
  ConstructGlobalToLocalMap();

  // recreate local wellbore tree by connecting locally relevant elems
  for( iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    globalIndex const ielemGlobal =
      m_localToGlobalMap[iwelemLocal];  // global index in elemManager ordering
    globalIndex const iwelemGlobal =
      ielemGlobal - elemOffsetGlobal;  // global index in well ordering
    globalIndex const iwelemNextGlobal =
      nextElemIdGlobal[iwelemGlobal];  // global index in well ordering
    globalIndex const ielemNextGlobal = elemOffsetGlobal +
      iwelemNextGlobal;  // global index in elemManager ordering

    if( iwelemNextGlobal < 0 )
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = -1;  // wellhead
      m_nextWellElementIndex[iwelemLocal] = -1;
      m_topWellElementIndex = iwelemLocal;
    }
    else
    {
      m_nextWellElementIndexGlobal[iwelemLocal] = ielemNextGlobal;  // wellhead

      if( globalToLocalMap().count( ielemNextGlobal ) > 0 )
      {
        m_nextWellElementIndex[iwelemLocal] = globalToLocalMap( ielemNextGlobal );
      }
      else
      {
        m_nextWellElementIndex[iwelemLocal] = -2;  // remote elem
      }
    }

    // TODO Change to LvArray::tensorOps::copy
    m_elementCenter[iwelemLocal][0] = elemCoordsGlobal[iwelemGlobal][0];
    m_elementCenter[iwelemLocal][1] = elemCoordsGlobal[iwelemGlobal][1];
    m_elementCenter[iwelemLocal][2] = elemCoordsGlobal[iwelemGlobal][2];

    m_elementVolume[iwelemLocal] = elemVolumeGlobal[iwelemGlobal];
    m_radius[iwelemLocal] = wellGeometry.GetElementRadius();

    // update local well elem to node map (note: nodes are in nodeManager ordering)

    // first get the global node indices in nodeManager ordering
    globalIndex const inodeTopGlobal = nodeOffsetGlobal +
      elemToNodesGlobal[iwelemGlobal][InternalWellGenerator::NodeLocation::TOP];
    globalIndex const inodeBottomGlobal = nodeOffsetGlobal +
      elemToNodesGlobal[iwelemGlobal][InternalWellGenerator::NodeLocation::BOTTOM];

    // then get the local node indices in nodeManager ordering
    localIndex const inodeTopLocal =
      nodeManager->globalToLocalMap( inodeTopGlobal );
    localIndex const inodeBottomLocal =
      nodeManager->globalToLocalMap( inodeBottomGlobal );

    m_toNodesRelation[iwelemLocal][InternalWellGenerator::NodeLocation::TOP] =
      inodeTopLocal;
    m_toNodesRelation[iwelemLocal][InternalWellGenerator::NodeLocation::BOTTOM] =
      inodeBottomLocal;
  }
}

void
WellElementSubRegion::UpdateNodeManagerNodeToElementMap( MeshLevel & mesh )
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager * const nodeManager = mesh.getNodeManager();

  // at this point, NodeManager::SetElementMaps has already been called for the mesh nodes
  // we have to update the following maps for the well nodes
  ArrayOfArrays< localIndex > & toElementRegionList =
    nodeManager->elementRegionList();
  ArrayOfArrays< localIndex > & toElementSubRegionList =
    nodeManager->elementSubRegionList();
  ArrayOfArrays< localIndex > & toElementList = nodeManager->elementList();

  // we get the region and subregion indices in the elemManager
  WellElementRegion const * const elemRegion =
    this->getParent()->getParent()->group_cast< WellElementRegion * >();
  string const elemRegionName = elemRegion->getName();

  localIndex const iregion = elemManager->GetRegions().getIndex( elemRegionName );
  localIndex const isubRegion = elemRegion->GetSubRegions().getIndex( getName() );

  // for each (new) well element
  for( localIndex iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    for( localIndex a = 0; a < numNodesPerElement(); ++a )
    {
      // get the local node index (in nodeManager ordering) using the elem-to-nodes maps constructed above
      localIndex const inodeLocal = m_toNodesRelation[iwelemLocal][a];

      // update the reverse map from well node to well element
      // this is needed to generate the adjacency list in communication setup phase
      toElementRegionList.emplaceBack( inodeLocal, iregion );
      toElementSubRegionList.emplaceBack( inodeLocal, isubRegion );
      toElementList.emplaceBack( inodeLocal, iwelemLocal );
    }
  }
  setupRelatedObjectsInRelations( &mesh );
}

void
WellElementSubRegion::ConnectPerforationsToMeshElements(
  MeshLevel & mesh,
  InternalWellGenerator const & wellGeometry )
{
  arrayView1d< R1Tensor const > const & perfCoordsGlobal =
    wellGeometry.GetPerfCoords();
  arrayView1d< real64 const > const & perfWellTransmissibilityGlobal =
    wellGeometry.GetPerfTransmissibility();

  m_perforationData.resize( perfCoordsGlobal.size() );
  localIndex iperfLocal = 0;

  // loop over all the perforations
  for( globalIndex iperfGlobal = 0; iperfGlobal < perfCoordsGlobal.size();
       ++iperfGlobal )
  {
    R1Tensor const & location = perfCoordsGlobal[iperfGlobal];

    localIndex erMatched = -1;
    localIndex esrMatched = -1;
    localIndex eiMatched = -1;

    localIndex erInit = -1;
    localIndex esrInit = -1;
    localIndex eiInit = -1;

    // for each perforation, we have to find the reservoir element that contains the perforation

    if( iperfLocal > 0 )
    {
      // get the info of the element matched with the previous perforation
      // this will be used next to search around this reservoir element
      erInit =
        m_perforationData.GetMeshElements().m_toElementRegion[iperfLocal - 1];
      esrInit =
        m_perforationData.GetMeshElements().m_toElementSubRegion[iperfLocal - 1];
      eiInit =
        m_perforationData.GetMeshElements().m_toElementIndex[iperfLocal - 1];
    }
    else
    {
      // Step 1: first, we search for the reservoir element that is the *closest* from the center of well element
      //         note that this reservoir element does not necessarily contain the center of the well element
      //         this "init" reservoir element will be used in SearchLocalElements to find the reservoir element that
      //         contains the well element
      InitializeLocalSearch( mesh, location, erInit, esrInit, eiInit );
    }

    // Step 2: then, search for the reservoir element that contains the well element
    //         to do that, we loop over the reservoir elements that are in the neighborhood of (erInit,esrInit,eiInit)
    bool resElemFound = SearchLocalElements( mesh,
                                             location,
                                             m_searchDepth,
                                             erInit,
                                             esrInit,
                                             eiInit,
                                             erMatched,
                                             esrMatched,
                                             eiMatched );

    // if the element was found
    if( resElemFound )
    {
      // set the indices for the matched reservoir element
      m_perforationData.GetMeshElements().m_toElementRegion[iperfLocal] =
        erMatched;
      m_perforationData.GetMeshElements().m_toElementSubRegion[iperfLocal] =
        esrMatched;
      m_perforationData.GetMeshElements().m_toElementIndex[iperfLocal] =
        eiMatched;

      // construct the local wellTransmissibility and location maps
      m_perforationData.GetWellTransmissibility()[iperfLocal] =
        perfWellTransmissibilityGlobal[iperfGlobal];
      m_perforationData.GetLocation()[iperfLocal] = location;

      // increment the local to global map
      m_perforationData.localToGlobalMap()[iperfLocal++] = iperfGlobal;
    }
  }

  // set the size based on the number of perforations matched with local reservoir elements
  m_perforationData.resize( iperfLocal );
  m_perforationData.ConstructGlobalToLocalMap();
}

void
WellElementSubRegion::ReconstructLocalConnectivity()
{
  // here we reconstruct the array m_nextWellElementIndexGlobal
  // this is needed after the addition of ghost well elements

  for( localIndex iwelemLocal = 0; iwelemLocal < size(); ++iwelemLocal )
  {
    globalIndex const nextGlobal = m_nextWellElementIndexGlobal[iwelemLocal];

    if( nextGlobal < 0 )  // well head
    {
      m_nextWellElementIndex[iwelemLocal] = -1;
      m_topWellElementIndex =
        iwelemLocal;  // reset this is case top element was added as ghost
    }
    else if( globalToLocalMap().count( nextGlobal ) == 0 )  // next is remote
    {
      m_nextWellElementIndex[iwelemLocal] = -2;
    }
    else  // local
    {
      m_nextWellElementIndex[iwelemLocal] = this->globalToLocalMap( nextGlobal );
    }
  }
}

bool
WellElementSubRegion::IsLocallyOwned() const
{
  return m_topRank == MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
}

void
WellElementSubRegion::ViewPackingExclusionList(
  SortedArray< localIndex > & exclusionList ) const
{
  ObjectManagerBase::ViewPackingExclusionList( exclusionList );
  exclusionList.insert( this->getWrapperIndex( viewKeyStruct::nodeListString ) );
}

localIndex
WellElementSubRegion::PackUpDownMapsSize(
  arrayView1d< localIndex const > const & packList ) const
{
  buffer_unit_type * junk = nullptr;
  return PackUpDownMapsPrivate< false >( junk, packList );
}

localIndex
WellElementSubRegion::PackUpDownMaps(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  return PackUpDownMapsPrivate< true >( buffer, packList );
}

template< bool DOPACK >
localIndex
WellElementSubRegion::PackUpDownMapsPrivate(
  buffer_unit_type *& buffer,
  arrayView1d< localIndex const > const & packList ) const
{
  localIndex packedSize = 0;

  packedSize += bufferOps::Pack< DOPACK >( buffer,
                                           nodeList().Base().toViewConst(),
                                           m_unmappedGlobalIndicesInNodelist,
                                           packList,
                                           this->localToGlobalMap(),
                                           nodeList().RelatedObjectLocalToGlobal() );

  return packedSize;
}

localIndex
WellElementSubRegion::UnpackUpDownMaps(
  buffer_unit_type const *& buffer,
  localIndex_array & packList,
  bool const GEOSX_UNUSED_PARAM( overwriteUpMaps ),
  bool const GEOSX_UNUSED_PARAM( overwriteDownMaps ) )
{
  localIndex unPackedSize = 0;

  unPackedSize += bufferOps::Unpack( buffer,
                                     nodeList().Base().toView(),
                                     packList,
                                     m_unmappedGlobalIndicesInNodelist,
                                     this->globalToLocalMap(),
                                     nodeList().RelatedObjectGlobalToLocal() );

  return unPackedSize;
}

void
WellElementSubRegion::FixUpDownMaps( bool const clearIfUnmapped )
{
  ObjectManagerBase::FixUpDownMaps( nodeList(),
                                    m_unmappedGlobalIndicesInNodelist,
                                    clearIfUnmapped );
}

void
WellElementSubRegion::DebugNodeManager( MeshLevel const & mesh ) const
{
  NodeManager const * const nodeManager = mesh.getNodeManager();
  // arrayView2d< real64 const > const & X = nodeManager->referencePosition();

  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1 )
  {
    return;
  }

  GEOSX_LOG_RANK( "++++++++++++++++++++++++++" );
  GEOSX_LOG_RANK( "WellElementSubRegion = " << getName() );
  GEOSX_LOG_RANK( "Number of local well elements = " << size() );
  GEOSX_LOG_RANK( "Number of local node elements = " << nodeManager->size() );

  if( nodeManager->size() > 0 )
  {
    return;
  }

  arrayView1d< globalIndex const > const & nodeLocalToGlobal =
    nodeManager->localToGlobalMap();
  for( localIndex inodeLocal = 0; inodeLocal < nodeManager->size(); ++inodeLocal )
  {
    std::cout << "nodeManager->localToGlobalMap[" << inodeLocal
              << "] = " << nodeLocalToGlobal[inodeLocal] << std::endl;
  }
}

void
WellElementSubRegion::DebugWellElementSubRegions(
  arrayView1d< integer const > const & elemStatusGlobal,
  globalIndex elemOffsetGlobal ) const
{
  if( size() == 0 )
  {
    return;
  }

  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) < 1 )
  {
    return;
  }

  GEOSX_LOG_RANK( "++++++++++++++++++++++++++" );
  GEOSX_LOG_RANK( "WellElementSubRegion = " << getName() );
  GEOSX_LOG_RANK( "Number of local well elements = " << size() );

  for( localIndex iwelem = 0; iwelem < size(); ++iwelem )
  {
    GEOSX_LOG_RANK_VAR( iwelem );
    GEOSX_LOG_RANK_VAR( m_nextWellElementIndex[iwelem] );
    GEOSX_LOG_RANK_VAR( m_nextWellElementIndexGlobal[iwelem] );
    GEOSX_LOG_RANK_VAR( m_elementCenter[iwelem] );
    GEOSX_LOG_RANK_VAR( m_localToGlobalMap[iwelem] );
    GEOSX_LOG_RANK_VAR(
      elemStatusGlobal[m_localToGlobalMap[iwelem] - elemOffsetGlobal] );
    GEOSX_LOG_RANK_VAR( m_topWellElementIndex );
  }
}

void
WellElementSubRegion::DebugWellElementSubRegionsAfterSetupCommunications() const
{
  if( size() == 0 )
  {
    return;
  }

  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1 )
  {
    return;
  }

  GEOSX_LOG_RANK( "++++++++++++++++++++++++++" );
  GEOSX_LOG_RANK( "WellElementSubRegion = " << getName() );
  GEOSX_LOG_RANK( "Number of local well elements = " << size() );
  GEOSX_LOG_RANK( "Number of ghost well elements = " << this->GetNumberOfGhosts() );

  for( localIndex iwelem = 0; iwelem < size(); ++iwelem )
  {
    GEOSX_LOG_RANK_VAR( iwelem );
    GEOSX_LOG_RANK_VAR( m_nextWellElementIndex[iwelem] );
    GEOSX_LOG_RANK_VAR( m_nextWellElementIndexGlobal[iwelem] );
    GEOSX_LOG_RANK_VAR( m_elementCenter[iwelem] );
    GEOSX_LOG_RANK_VAR( m_localToGlobalMap[iwelem] );
    GEOSX_LOG_RANK_VAR( m_ghostRank[iwelem] );
    GEOSX_LOG_RANK_VAR( m_topWellElementIndex );
  }
}

}  // namespace geosx
