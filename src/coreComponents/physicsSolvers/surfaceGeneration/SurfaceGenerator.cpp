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
 * @file SurfaceGenerator.cpp
 */

#include "SurfaceGenerator.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"


#ifdef USE_GEOSX_PTP
#include "physicsSolvers/GEOSX_PTP/ParallelTopologyChange.hpp"
#endif

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

constexpr real64 SurfaceGenerator::m_nonRuptureTime;

void ModifiedObjectLists::clearNewFromModified()
{
  for( localIndex const a : newNodes )
  {
    modifiedNodes.erase( a );
  }

  for( localIndex const a : newEdges )
  {
    modifiedEdges.erase( a );
  }

  for( localIndex const a : newFaces )
  {
    modifiedFaces.erase( a );
  }
}

void ModifiedObjectLists::insert( ModifiedObjectLists const & modifiedObjects )
{
  newNodes.insert( modifiedObjects.newNodes.begin(),
                   modifiedObjects.newNodes.end() );
  modifiedNodes.insert( modifiedObjects.modifiedNodes.begin(),
                        modifiedObjects.modifiedNodes.end() );

  newEdges.insert( modifiedObjects.newEdges.begin(),
                   modifiedObjects.newEdges.end() );
  modifiedEdges.insert( modifiedObjects.modifiedEdges.begin(),
                        modifiedObjects.modifiedEdges.end() );

  newFaces.insert( modifiedObjects.newFaces.begin(),
                   modifiedObjects.newFaces.end() );
  modifiedFaces.insert( modifiedObjects.modifiedFaces.begin(),
                        modifiedObjects.modifiedFaces.end() );

  for( auto & iter : modifiedObjects.newElements )
  {
    std::pair< localIndex, localIndex > const & key = iter.first;
    std::set< localIndex > const & values = iter.second;
    newElements[key].insert( values.begin(), values.end() );
  }

  for( auto & iter : modifiedObjects.modifiedElements )
  {
    std::pair< localIndex, localIndex > const & key = iter.first;
    std::set< localIndex > const & values = iter.second;
    modifiedElements[key].insert( values.begin(), values.end() );
  }

}

static localIndex GetOtherFaceEdge( const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                                    const localIndex thisFace, const localIndex thisEdge )
{
  localIndex nextEdge = LOCALINDEX_MAX;

  const std::pair< localIndex, localIndex > & faceToEdges = stlMapLookup( localFacesToEdges, thisFace );
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
    GEOSX_ERROR( "SurfaceGenerator::Couldn't find thisEdge in localFacesToEdges[thisFace]" );
  }
  return nextEdge;
}

static void CheckForAndRemoveDeadEndPath( const localIndex edgeIndex,
                                          arrayView1d< integer const > const & isEdgeExternal,
                                          map< localIndex, std::set< localIndex > > & edgesToRuptureReadyFaces,
                                          map< localIndex, std::pair< localIndex, localIndex > > & localVFacesToVEdges,
                                          std::set< localIndex > & nodeToRuptureReadyFaces )
{


  localIndex thisEdge = edgeIndex;

  // if the edge is internal and the edge is only attached to one ruptured face...
  while( isEdgeExternal[thisEdge]!=1 )
  {

    //    std::set<localIndex>& edgeToRuptureReadyFaces = stlMapLookup(edgesToRuptureReadyFaces,thisEdge);
    std::set< localIndex > & edgeToRuptureReadyFaces = edgesToRuptureReadyFaces[thisEdge];

    if( edgeToRuptureReadyFaces.size()!=1 )
      break;

    // then the index for the face that is a "dead end"
    localIndex deadEndFace = *(edgeToRuptureReadyFaces.begin());


    std::pair< localIndex, localIndex > & localVFaceToVEdges = stlMapLookup( localVFacesToVEdges, deadEndFace );

    // get the edge on the other side of the "dead end" face
    localIndex nextEdge = -1;
    if( localVFaceToVEdges.first == thisEdge )
      nextEdge = localVFaceToVEdges.second;
    else if( localVFaceToVEdges.second == thisEdge )
      nextEdge = localVFaceToVEdges.first;
    else
    {
      GEOSX_ERROR( "SurfaceGenerator::FindFracturePlanes: Could not find the next edge when removing dead end faces." );
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


SurfaceGenerator::SurfaceGenerator( const string & name,
                                    Group * const parent ):
  SolverBase( name, parent ),
  m_failCriterion( 1 ),
//  m_maxTurnAngle(91.0),
  m_nodeBasedSIF( 0 ),
  m_rockToughness( 1.0e99 ),
  m_mpiCommOrder( 0 )
{
  this->registerWrapper( viewKeyStruct::failCriterionString(), &this->m_failCriterion );

  registerWrapper( viewKeyStruct::solidMaterialNameString(), &m_solidMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid material used in solid mechanic solver" );

  registerWrapper( viewKeyStruct::rockToughnessString(), &m_rockToughness ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Rock toughness of the solid material" );

  registerWrapper( viewKeyStruct::nodeBasedSIFString(), &m_nodeBasedSIF ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Rock toughness of the solid material" );

  registerWrapper( viewKeyStruct::mpiCommOrderString(), &m_mpiCommOrder ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to enable MPI consistent communication ordering" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( "Fracture" );

  registerWrapper( viewKeyStruct::tipNodesString(), &m_tipNodes ).
    setDescription( "Set containing all the nodes at the fracture tip" );

  registerWrapper( viewKeyStruct::tipEdgesString(), &m_tipEdges ).
    setDescription( "Set containing all the tip edges" );

  registerWrapper( viewKeyStruct::tipFacesString(), &m_tipFaces ).
    setDescription( "Set containing all the tip faces" );

  registerWrapper( viewKeyStruct::trailingFacesString(), &m_trailingFaces ).
    setDescription( "Set containing all the trailing faces" );


  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

SurfaceGenerator::~SurfaceGenerator()
{
  // TODO Auto-generated destructor stub
}

void SurfaceGenerator::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::K_IC_00,
                                       extrinsicMeshData::K_IC_01,
                                       extrinsicMeshData::K_IC_02,
                                       extrinsicMeshData::K_IC_10,
                                       extrinsicMeshData::K_IC_11,
                                       extrinsicMeshData::K_IC_12,
                                       extrinsicMeshData::K_IC_20,
                                       extrinsicMeshData::K_IC_21,
                                       extrinsicMeshData::K_IC_22 >( this->getName() );
    } );

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      subRegion.registerExtrinsicData< extrinsicMeshData::K_IC_00,
                                       extrinsicMeshData::K_IC_01,
                                       extrinsicMeshData::K_IC_02,
                                       extrinsicMeshData::K_IC_10,
                                       extrinsicMeshData::K_IC_11,
                                       extrinsicMeshData::K_IC_12,
                                       extrinsicMeshData::K_IC_20,
                                       extrinsicMeshData::K_IC_21,
                                       extrinsicMeshData::K_IC_22,
                                       extrinsicMeshData::RuptureTime,
                                       extrinsicMeshData::RuptureRate >( this->getName() );
    } );

    NodeManager & nodeManager = meshLevel.getNodeManager();
    EdgeManager & edgeManager = meshLevel.getEdgeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();

    nodeManager.registerExtrinsicData< extrinsicMeshData::ParentIndex,
                                       extrinsicMeshData::ChildIndex,
                                       extrinsicMeshData::DegreeFromCrack,
                                       extrinsicMeshData::DegreeFromCrackTip,
                                       extrinsicMeshData::SIFNode,
                                       extrinsicMeshData::RuptureTime >( this->getName() );

    edgeManager.registerExtrinsicData< extrinsicMeshData::ParentIndex,
                                       extrinsicMeshData::ChildIndex,
                                       extrinsicMeshData::SIF_I,
                                       extrinsicMeshData::SIF_II,
                                       extrinsicMeshData::SIF_III >( this->getName() );

    faceManager.registerExtrinsicData< extrinsicMeshData::ParentIndex,
                                       extrinsicMeshData::ChildIndex,
                                       extrinsicMeshData::RuptureState,
                                       extrinsicMeshData::RuptureTime,
                                       extrinsicMeshData::SIFonFace,
                                       extrinsicMeshData::K_IC,
                                       extrinsicMeshData::PrimaryCandidateFace,
                                       extrinsicMeshData::IsFaceSeparable,
                                       extrinsicMeshData::DegreeFromCrackTip >( this->getName() );

    // TODO: handle this automatically in registerExtrinsicData()
    faceManager.getExtrinsicData< extrinsicMeshData::K_IC >().resizeDimension< 1 >( 3 );
  } );
}

void SurfaceGenerator::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );//this->getGroupByPath<DomainPartition>("/Problem/domain");
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    NodeManager & nodeManager = meshLevel.getNodeManager();
    FaceManager & faceManager = meshLevel.getFaceManager();

    arrayView1d< localIndex > const & parentNodeIndex = nodeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

    arrayView1d< localIndex > const & parentFaceIndex = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

    arrayView1d< localIndex > const & childFaceIndex = faceManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();

    parentNodeIndex.setValues< serialPolicy >( -1 );
    parentFaceIndex.setValues< serialPolicy >( -1 );
    childFaceIndex.setValues< serialPolicy >( -1 );

    m_originalNodetoFaces = nodeManager.faceList();
    m_originalNodetoEdges = nodeManager.edgeList();
    m_originalFaceToEdges = faceManager.edgeList();

    nodeManager.registerWrapper( "usedFaces", &m_usedFacesForNode );
    m_usedFacesForNode.resize( nodeManager.size() );

    localIndex const numFaces = faceManager.size();
    m_originalFacesToElemRegion.resize( numFaces, 2 );
    m_originalFacesToElemSubRegion.resize( numFaces, 2 );
    m_originalFacesToElemIndex.resize( numFaces, 2 );

    for( localIndex faceID = 0; faceID < numFaces; ++faceID )
    {
      for( localIndex side = 0; side < 2; ++side )
      {
        m_originalFacesToElemRegion( faceID, side ) = faceManager.elementRegionList()( faceID, side );
        m_originalFacesToElemSubRegion( faceID, side ) = faceManager.elementSubRegionList()( faceID, side );
        m_originalFacesToElemIndex( faceID, side ) = faceManager.elementList()( faceID, side );
      }
    }
  }

  domain.getMeshBodies().forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
    FaceManager & faceManager = meshLevel.getFaceManager();
    ElementRegionManager & elementManager = meshLevel.getElemManager();
    arrayView2d< real64 const > const & faceNormals = faceManager.faceNormal();

    //TODO: roughness to KIC should be made a material constitutive relationship.
    arrayView2d< real64 > const & KIC = faceManager.getExtrinsicData< extrinsicMeshData::K_IC >();

    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      if( m_rockToughness >= 0 )
      {
        KIC[kf][0] = m_rockToughness;
        KIC[kf][1] = m_rockToughness;
        KIC[kf][2] = m_rockToughness;
      }
      else
      {
        arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
        arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
        arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

        for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
        {
          localIndex const er = faceToRegionMap[kf][k];
          localIndex const esr = faceToSubRegionMap[kf][k];
          localIndex const ei = faceToElementMap[kf][k];

          if( er != -1 &&  esr != -1 && ei != -1 )
          {
            CellElementSubRegion & elementSubRegion = elementManager.getRegion( faceToRegionMap[kf][k] ).
                                                        getSubRegion< CellElementSubRegion >( faceToSubRegionMap[kf][k] );
            localIndex iEle = faceToElementMap[kf][k];

            arrayView1d< real64 const > const K_IC_00 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_00 >();
            arrayView1d< real64 const > const K_IC_01 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_01 >();
            arrayView1d< real64 const > const K_IC_02 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_02 >();
            arrayView1d< real64 const > const K_IC_10 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_10 >();
            arrayView1d< real64 const > const K_IC_11 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_11 >();
            arrayView1d< real64 const > const K_IC_12 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_12 >();
            arrayView1d< real64 const > const K_IC_20 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_20 >();
            arrayView1d< real64 const > const K_IC_21 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_21 >();
            arrayView1d< real64 const > const K_IC_22 = elementSubRegion.getExtrinsicData< extrinsicMeshData::K_IC_22 >();

            real64 k0[3];
            k0[0] = K_IC_00[iEle]*faceNormals[kf][0] + K_IC_10[iEle]*faceNormals[kf][1] + K_IC_20[iEle]*faceNormals[kf][2];
            k0[1] = K_IC_01[iEle]*faceNormals[kf][0] + K_IC_11[iEle]*faceNormals[kf][1] + K_IC_21[iEle]*faceNormals[kf][2];
            k0[2] = K_IC_02[iEle]*faceNormals[kf][0] + K_IC_12[iEle]*faceNormals[kf][1] + K_IC_22[iEle]*faceNormals[kf][2];

            KIC[kf][0] = std::min( std::fabs( k0[0] ), std::fabs( KIC[kf][0] ) );
            KIC[kf][1] = std::min( std::fabs( k0[1] ), std::fabs( KIC[kf][1] ) );
            KIC[kf][2] = std::min( std::fabs( k0[2] ), std::fabs( KIC[kf][2] ) );
          }
        }
      }
    }
  } );
}


void SurfaceGenerator::postRestartInitialization()
{
  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );//this->getGroupByPath<DomainPartition>("/Problem/domain");

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  // repopulate the fracture stencil
  domain.getMeshBodies().forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    EdgeManager & edgeManager = meshLevel.getEdgeManager();
    ElementRegionManager & elemManager = meshLevel.getElemManager();
    SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
    FaceElementSubRegion & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( 0 );

    for( localIndex fce = 0; fce < edgeManager.m_fractureConnectorEdgesToFaceElements.size(); ++fce )
    {
      edgeManager.m_recalculateFractureConnectorEdges.insert( fce );
    }

    for( localIndex fe = 0; fe < fractureSubRegion.size(); ++fe )
    {
      fractureSubRegion.m_newFaceElements.insert( fe );
    }

    for( localIndex a = 0; a < fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addToFractureStencil( meshLevel, this->m_fractureRegionName, false );
        edgeManager.m_recalculateFractureConnectorEdges.clear();
        fractureSubRegion.m_newFaceElements.clear();
      }
    }
  } );
}


real64 SurfaceGenerator::solverStep( real64 const & time_n,
                                     real64 const & dt,
                                     const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                     DomainPartition & domain )
{
  int rval = 0;

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    {
      SpatialPartition & partition = dynamicCast< SpatialPartition & >( domain.getReference< PartitionBase >( dataRepository::keys::partitionManager ) );

      rval = separationDriver( domain,
                               meshLevel,
                               domain.getNeighbors(),
                               partition.getColor(),
                               partition.numColor(),
                               0,
                               time_n + dt );
    }
  }

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    {
      ElementRegionManager & elemManager = meshLevel.getElemManager();
      EdgeManager & edgeManager = meshLevel.getEdgeManager();
      SurfaceElementRegion & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );

      for( localIndex a=0; a<fvManager.numSubGroups(); ++a )
      {
        FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
        if( fluxApprox!=nullptr )
        {
          fluxApprox->addToFractureStencil( meshLevel, this->m_fractureRegionName, true );
          edgeManager.m_recalculateFractureConnectorEdges.clear();
          fractureRegion.getSubRegion< FaceElementSubRegion >( 0 ).m_newFaceElements.clear();
        }
      }
    }
  }

  return rval;
}



int SurfaceGenerator::separationDriver( DomainPartition & domain,
                                        MeshLevel & mesh,
                                        std::vector< NeighborCommunicator > & neighbors,
                                        int const tileColor,
                                        int const numTileColors,
                                        bool const prefrac,
                                        real64 const time_np1 )
{
  GEOSX_MARK_FUNCTION;

  m_faceElemsRupturedThisSolve.clear();
  NodeManager & nodeManager = mesh.getNodeManager();
  EdgeManager & edgeManager = mesh.getEdgeManager();
  FaceManager & faceManager = mesh.getFaceManager();
  ElementRegionManager & elementManager = mesh.getElemManager();

  std::vector< std::set< localIndex > > nodesToRupturedFaces;
  std::vector< std::set< localIndex > > edgesToRupturedFaces;

  ArrayOfArrays< localIndex > const & nodeToElementMap = nodeManager.elementList();

  map< string, string_array > fieldNames;
  fieldNames["face"].emplace_back( string( extrinsicMeshData::RuptureState::key() ) );
  fieldNames["node"].emplace_back( string( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternalString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

  elementManager.forElementSubRegions< CellElementSubRegion >( [] ( auto & elemSubRegion )
  {
    elemSubRegion.moveSets( LvArray::MemorySpace::CPU );
  } );
  faceManager.moveSets( LvArray::MemorySpace::CPU );
  edgeManager.moveSets( LvArray::MemorySpace::CPU );
  nodeManager.moveSets( LvArray::MemorySpace::CPU );

  if( !prefrac )
  {
    if( m_failCriterion >0 )  // Stress intensity factor based criterion and mixed criterion.
    {

      identifyRupturedFaces( domain,
                             nodeManager,
                             edgeManager,
                             faceManager,
                             elementManager,
                             prefrac );

    }
  }


  if( prefrac )
  {
    ModifiedObjectLists modifiedObjects;
    calculateKinkAngles( faceManager, edgeManager, nodeManager, modifiedObjects, prefrac );
  }

  // We do this here to get the nodesToRupturedFaces etc.
  // The fail stress check inside has been disabled
  postUpdateRuptureStates( nodeManager,
                           edgeManager,
                           faceManager,
                           elementManager,
                           nodesToRupturedFaces,
                           edgesToRupturedFaces );

  int rval = 0;
  //  array1d<MaterialBaseStateDataT*>&  temp = elementManager.m_ElementRegions["PM1"].m_materialStates;

  array1d< integer > const & isNodeGhost = nodeManager.ghostRank();

  for( int color=0; color<numTileColors; ++color )
  {
    ModifiedObjectLists modifiedObjects;
    if( color==tileColor )
    {
      for( localIndex a=0; a<nodeManager.size(); ++a )
      {
        int didSplit = 0;
        if( isNodeGhost[a]<0 &&
            nodeToElementMap.sizeOfArray( a )>1 )
        {
          didSplit += processNode( a,
                                   time_np1,
                                   nodeManager,
                                   edgeManager,
                                   faceManager,
                                   elementManager,
                                   nodesToRupturedFaces,
                                   edgesToRupturedFaces,
                                   elementManager,
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

    // 1) Assign new global indices to the new objects
    CommunicationTools::getInstance().assignNewGlobalIndices( nodeManager, modifiedObjects.newNodes );
    CommunicationTools::getInstance().assignNewGlobalIndices( edgeManager, modifiedObjects.newEdges );
    CommunicationTools::getInstance().assignNewGlobalIndices( faceManager, modifiedObjects.newFaces );
//    CommunicationTools::getInstance().AssignNewGlobalIndices( elementManager, modifiedObjects.newElements );

    ModifiedObjectLists receivedObjects;

    /// Nodes to edges in process node is not being set on rank 2. need to check that the new node->edge map is properly
    /// communicated
    ParallelTopologyChange::synchronizeTopologyChange( &mesh,
                                                       neighbors,
                                                       modifiedObjects,
                                                       receivedObjects,
                                                       m_mpiCommOrder );

    synchronizeTipSets( faceManager,
                        edgeManager,
                        nodeManager,
                        receivedObjects );


#else

    GEOSX_UNUSED_VAR( neighbors );
    assignNewGlobalIndicesSerial( nodeManager, modifiedObjects.newNodes );
    assignNewGlobalIndicesSerial( edgeManager, modifiedObjects.newEdges );
    assignNewGlobalIndicesSerial( faceManager, modifiedObjects.newFaces );

#endif

    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    elementManager.forElementSubRegionsComplete< FaceElementSubRegion >( [&]( localIndex const er,
                                                                              localIndex const esr,
                                                                              ElementRegionBase &,
                                                                              FaceElementSubRegion & subRegion )
    {
      std::set< localIndex > & newFaceElems = modifiedObjects.newElements[{er, esr}];
      for( localIndex const newFaceElemIndex : newFaceElems )
      {
        subRegion.m_newFaceElements.insert( newFaceElemIndex );
      }
    } );


    elementManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      FaceElementSubRegion::NodeMapType & nodeMap = subRegion.nodeList();
      FaceElementSubRegion::FaceMapType & faceMap = subRegion.faceList();

      for( localIndex kfe=0; kfe<subRegion.size(); ++kfe )
      {
        nodeMap.resizeArray( kfe, 8 );

        localIndex const numNodesInFace = faceToNodeMap.sizeOfArray( faceMap[ kfe ][ 0 ] );
        for( localIndex a = 0; a < numNodesInFace; ++a )
        {
          localIndex const aa = a < 2 ? a : numNodesInFace - a + 1;
          localIndex const bb = aa == 0 ? aa : numNodesInFace - aa;

          // TODO HACK need to generalize to something other than quads
          //wu40: I temporarily make it work for tet mesh. Need further check with Randy.
          nodeMap[ kfe ][ a ]   = faceToNodeMap( faceMap[ kfe ][ 0 ], aa );
          nodeMap[ kfe ][ a + numNodesInFace ] = faceToNodeMap( faceMap[ kfe ][ 1 ], bb );
        }

        if( numNodesInFace == 3 )
        {
          nodeMap[kfe][6] = faceToNodeMap( faceMap[ kfe ][ 0 ], 2 );
          nodeMap[kfe][7] = faceToNodeMap( faceMap[ kfe ][ 1 ], 2 );
        }
      }
    } );
  }


  real64 ruptureRate = calculateRuptureRate( elementManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName ), edgeManager );

  GEOSX_LOG_LEVEL_RANK_0( 3, "rupture rate is " << ruptureRate );
  if( ruptureRate > 0 )
    m_nextDt = ruptureRate < 1e99 ? m_cflFactor / ruptureRate : 1e99;


//  if( rval>0 )
  {
    elementManager.forElementSubRegions< CellElementSubRegion >( [] ( auto & elemSubRegion )
    {
      elemSubRegion.nodeList().registerTouch( LvArray::MemorySpace::CPU );
      elemSubRegion.edgeList().registerTouch( LvArray::MemorySpace::CPU );
      elemSubRegion.faceList().registerTouch( LvArray::MemorySpace::CPU );
    } );


    faceManager.nodeList().toView().registerTouch( LvArray::MemorySpace::CPU );
//    faceManager.edgeList().registerTouch( LvArray::MemorySpace::CPU );
    faceManager.elementList().registerTouch( LvArray::MemorySpace::CPU );
    faceManager.elementRegionList().registerTouch( LvArray::MemorySpace::CPU );
    faceManager.elementSubRegionList().registerTouch( LvArray::MemorySpace::CPU );

    edgeManager.nodeList().registerTouch( LvArray::MemorySpace::CPU );

//    nodeManager.edgeList().registerTouch( LvArray::MemorySpace::CPU );
//    nodeManager.faceList()().registerTouch( LvArray::MemorySpace::CPU );
//    nodeManager.elementList().registerTouch( LvArray::MemorySpace::CPU );
//    nodeManager.elementRegionList().registerTouch( LvArray::MemorySpace::CPU );
//    nodeManager.elementSubRegionList().registerTouch( LvArray::MemorySpace::CPU );

  }

  return rval;
}

void SurfaceGenerator::synchronizeTipSets ( FaceManager & faceManager,
                                            EdgeManager & edgeManager,
                                            NodeManager & nodeManager,
                                            ModifiedObjectLists & receivedObjects )
{
  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  for( localIndex const nodeIndex : receivedObjects.newNodes )
  {
    localIndex const parentNodeIndex = parentNodeIndices[nodeIndex];

    GEOSX_ERROR_IF( parentNodeIndex == -1, "parentNodeIndex should not be -1" );

    m_tipNodes.remove( parentNodeIndex );
  }

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();
  arrayView1d< integer > const & edgeIsExternal = edgeManager.isExternal();
  arrayView1d< integer > const & nodeIsExternal = nodeManager.isExternal();


  arrayView1d< localIndex const > const &
  parentEdgeIndices = edgeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  arrayView1d< localIndex const > const &
  childEdgeIndices = edgeManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();


  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  for( localIndex const edgeIndex : receivedObjects.newEdges )
  {
    localIndex const parentEdgeIndex = parentEdgeIndices[edgeIndex];

    GEOSX_ERROR_IF( parentEdgeIndex == -1, "parentEdgeIndex should not be -1" );

    m_tipEdges.remove( parentEdgeIndex );
    for( localIndex const faceIndex : edgeToFaceMap[ parentEdgeIndex ] )
    {
      bool trailingFace = false;
      if( m_trailingFaces.contains( faceIndex ))
      {
        for( localIndex const faceLocalEdgeIndex : faceToEdgeMap[ faceIndex ] )
        {
          if( m_tipEdges.contains( faceLocalEdgeIndex ))
          {
            trailingFace = true;
          }
        }

        if( trailingFace == false )
        {
          m_trailingFaces.remove( faceIndex );
        }
      }
    }
  }

  arrayView1d< integer const > const & isFaceSeparable = faceManager.getExtrinsicData< extrinsicMeshData::IsFaceSeparable >();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

  arrayView1d< localIndex const > const & childNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();
  arrayView1d< localIndex > const & parentFaceIndices = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  for( localIndex const faceIndex : receivedObjects.newFaces )
  {
    localIndex const parentFaceIndex = parentFaceIndices[faceIndex];
    GEOSX_ERROR_IF( parentFaceIndex == -1, "parentFaceIndex should not be -1" );

    m_trailingFaces.insert( parentFaceIndex );
    m_tipFaces.remove( parentFaceIndex );

    for( localIndex const edgeIndex : faceManager.edgeList()[ parentFaceIndex ] )
    {
      if( parentEdgeIndices[edgeIndex]==-1 && childEdgeIndices[edgeIndex]==-1 )
      {
        m_tipEdges.insert( edgeIndex );

        for( localIndex const iface: edgeManager.faceList()[ edgeIndex ] )
        {
          if( faceToElementMap.size( 1 ) == 2  &&
              faceIsExternal[iface] < 1 &&
              isFaceSeparable[iface] == 1 )
          {
            m_tipFaces.insert( iface );
          }
        }
      }
      if( edgeIsExternal[edgeIndex]==0 )
      {
        edgeIsExternal[edgeIndex] = 2;
      }
    }
    for( localIndex const nodeIndex : faceManager.nodeList()[ parentFaceIndex ] )
    {
      if( parentNodeIndices[nodeIndex]==-1 && childNodeIndices[nodeIndex]==-1 )
      {
        m_tipNodes.insert( nodeIndex );
      }
      if( nodeIsExternal[nodeIndex] )
      {
        nodeIsExternal[nodeIndex] = 2;
      }
    }
  }
}


//void SurfaceGenerator::setDegreeFromCrackTip( NodeManager & nodeManager,
//                                              FaceManager & faceManager )
//{
//
//  arrayView1d<integer> &
//  nodeDegreeFromCrackTip = nodeManager.getReference<integer_array>( viewKeyStruct::degreeFromCrackTipString() );
//
//  arrayView1d<integer> &
//  faceDegreeFromCrackTip = faceManager.getReference<integer_array>( viewKeyStruct::degreeFromCrackTipString() );
//
//  ArrayOfArraysView< localIndex const > const & facesToNodes = faceManager.nodeList();
//
//  arrayView1d<integer const > const & ruptureState = faceManager.getReference<integer_array>( "ruptureState" );
//
//  faceDegreeFromCrackTip = 100000;
//
//  for( localIndex kf=0 ; kf<faceManager.size() ; ++kf )
//  {
//    if( ruptureState(kf) >=2 )
//    {
//      for( localIndex a=0 ; a<facesToNodes.sizeOfArray(kf) ; ++a )
//      {
//        localIndex const nodeIndex = facesToNodes(kf,a);
//        if( )
//      }
//    }
//  }
//}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::processNode( const localIndex nodeID,
                                    real64 const time_np1,
                                    NodeManager & nodeManager,
                                    EdgeManager & edgeManager,
                                    FaceManager & faceManager,
                                    ElementRegionManager & elemManager,
                                    std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                                    std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                                    ElementRegionManager & elementManager,
                                    ModifiedObjectLists & modifiedObjects,
                                    const bool GEOSX_UNUSED_PARAM( prefrac ) )
{
  bool didSplit = false;
  bool fracturePlaneFlag = true;

  {
    std::set< localIndex > facialRupturePath;
    map< localIndex, int > edgeLocations;
    map< localIndex, int > faceLocations;
    map< std::pair< CellElementSubRegion *, localIndex >, int > elemLocations;


    fracturePlaneFlag = findFracturePlanes( nodeID,
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
      mapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

      didSplit = true;
      performFracture( nodeID,
                       time_np1,
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
      mapConsistencyCheck( nodeID, nodeManager, edgeManager, faceManager, elementManager, elemLocations );

    }
  }

  return didSplit;
}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::findFracturePlanes( const localIndex nodeID,
                                           const NodeManager & nodeManager,
                                           const EdgeManager & edgeManager,
                                           const FaceManager & faceManager,
                                           ElementRegionManager & elemManager,
                                           const std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                                           const std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                                           std::set< localIndex > & separationPathFaces,
                                           map< localIndex, int > & edgeLocations,
                                           map< localIndex, int > & faceLocations,
                                           map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations )
{
  arrayView1d< localIndex const > const & parentNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  localIndex const parentNodeIndex = ObjectManagerBase::getParentRecusive( parentNodeIndices, nodeID );

  arrayView1d< localIndex const > const & parentFaceIndices = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();
  arrayView1d< localIndex const > const & childFaceIndices = faceManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();

  std::set< localIndex > const & vNodeToRupturedFaces = nodesToRupturedFaces[parentNodeIndex];

  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

//  array1d< ReferenceWrapper<localIndex_array> > nodeToElements
//  const std::set< std::pair<CellBlockSubRegion*,localIndex> >&
//  nodeToElementMaps = nodeManager.m_toElementsRelation[nodeID] ;

  arraySlice1d< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList()[nodeID];
  arraySlice1d< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList()[nodeID];
  arraySlice1d< localIndex const > const & nodeToElementMap = nodeManager.elementList()[nodeID];

  // ***** BACKWARDS COMPATIBLITY HACK
  std::set< std::pair< CellElementSubRegion *, localIndex > > nodeToElementMaps;


  for( localIndex k=0; k<nodeManager.elementRegionList().sizeOfArray( nodeID ); ++k )
  {
    nodeToElementMaps.emplace( &elemManager.getRegion( nodeToRegionMap[k] ).
                                 getSubRegion< CellElementSubRegion >( nodeToSubRegionMap[k] ),
                               nodeToElementMap[k] );
  }


  // ***** END BACKWARDS COMPATIBLITY HACK


  arrayView1d< integer const > const & isEdgeExternal = edgeManager.isExternal();

//  const std::set<localIndex>& usedFaces = nodeManager.GetUnorderedVariableOneToManyMap("usedFaces")[nodeID];

  // **** local working arrays *****************************************************************************************

  // array to hold the faces ready for rupture. It is filled with the intersection of the virtual parent faces
  // associated
  // with all faces attached to the node, and all ruptured virtual faces attached to the virtual parent node.
  std::set< localIndex > nodeToRuptureReadyFaces;
  for( localIndex const i : nodeToFaceMap[ nodeID ] )
  {
    const localIndex parentFaceIndex = ( parentFaceIndices[i] == -1 ) ? i : parentFaceIndices[i];

    if( vNodeToRupturedFaces.count( parentFaceIndex ) > 0 )
    {
      nodeToRuptureReadyFaces.insert( parentFaceIndex );
    }
  }


  // local map to hold the edgesToRuptureReadyFaces
  map< localIndex, std::set< localIndex > > edgesToRuptureReadyFaces;
  for( localIndex const edgeIndex : m_originalNodetoEdges[ parentNodeIndex ] )
  {
    if( !(edgesToRupturedFaces[edgeIndex].empty()) )
      edgesToRuptureReadyFaces[edgeIndex].insert( edgesToRupturedFaces[edgeIndex].begin(), edgesToRupturedFaces[edgeIndex].end() );
  }


  // need a map from faces to edges that are attached to the node
  map< localIndex, std::pair< localIndex, localIndex > > nodeLocalFacesToEdges;
  for( localIndex const kf : m_originalNodetoFaces[ parentNodeIndex ] )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( localIndex const ke : m_originalFaceToEdges[ kf ] )
    {
      if( m_originalNodetoEdges.contains( parentNodeIndex, ke ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOSX_ERROR( "SurfaceGenerator::FindFracturePlanes: invalid edge." );
    }


    nodeLocalFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // ***** remove dead end paths ***************************************************************************************
  // if the edge is not external, and the size of edgesToRupturedFaces is less than 2, then the edge is a dead-end
  // as far as a rupture plane is concerned. The face associated with the edge should be removed from the working
  // list of ruptured faces.

  // loop over all the edges
  for( localIndex const edgeIndex : m_originalNodetoEdges[ parentNodeIndex ] )
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

  for( std::set< localIndex >::const_iterator i=nodeToRuptureReadyFaces.begin(); i!=nodeToRuptureReadyFaces.end(); ++i )
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
    //    GEOSX_ERROR("Fracturantor3::FindFracturePlanes: couldn't set starting face/edge");
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
    map< localIndex, int > facesInPath;
    map< localIndex, int > edgesInPath;

    int numFacesInPath = 0;
    edgesInPath[thisEdge] = numFacesInPath;
    facesInPath[thisFace] = numFacesInPath++;

    localIndex_array facePath;
    localIndex_array edgePath;

    facePath.emplace_back( thisFace );
    edgePath.emplace_back( thisEdge );

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
          std::cout<<"  Face Separation Path = " << facePath << std::endl;
          std::cout<<"  Edge Separation Path = " << facePath << std::endl;

          GEOSX_ERROR( "crap" );
        }

        // add faces in the path to separationPathFaces
        for( map< localIndex, int >::const_iterator kf=facesInPath.begin(); kf!=facesInPath.end(); ++kf )
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

            std::pair< CellElementSubRegion *, localIndex >
            thisElem0 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[thisFace][0] ).
                                          getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[thisFace][0] ),
                                        m_originalFacesToElemIndex[thisFace][0] );

            std::pair< CellElementSubRegion *, localIndex >
            thisElem1 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[thisFace][1] ).
                                          getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[thisFace][1] ),
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

            for( std::set< localIndex >::const_iterator iter_edgeToFace = edgesToRuptureReadyFaces[nextEdge].begin();
                 iter_edgeToFace!=edgesToRuptureReadyFaces[nextEdge].end(); ++iter_edgeToFace )
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

                std::pair< CellElementSubRegion *, localIndex >
                nextElem0 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[candidateFaceIndex][0] ).
                                              getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[candidateFaceIndex][0] ),
                                            m_originalFacesToElemIndex[candidateFaceIndex][0] );

                std::pair< CellElementSubRegion *, localIndex >
                nextElem1 = std::make_pair( &elemManager.getRegion( m_originalFacesToElemRegion[candidateFaceIndex][1] ).
                                              getSubRegion< CellElementSubRegion >( m_originalFacesToElemSubRegion[candidateFaceIndex][1] ),
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
              GEOSX_ERROR( "SurfaceGenerator::FindFracturePlanes: couldn't find the next face in the rupture path" );
            }
          }

          //        lastEdge = thisEdge;
          //        lastFace = thisFace;

          thisEdge = nextEdge;
          thisFace = nextFace;
          //      separationPathFaces.insert( thisFace );
          edgesInPath[thisEdge] = numFacesInPath;
          facesInPath[thisFace] = numFacesInPath++;

          facePath.emplace_back( thisFace );
          edgePath.emplace_back( thisEdge );

        }
        else
        {
          GEOSX_ERROR( "SurfaceGenerator::next edge in separation path is apparently  connected to less than 2 ruptured face" );
        }

      }
    }
  }


  //***** SET LOCATIONS ************************************************************************************************



  // need a map from faces to edges that are attached to the node
  map< localIndex, std::pair< localIndex, localIndex > > localFacesToEdges;
  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
  {
    localIndex edge[2] = { INT_MAX, INT_MAX };
    int count = 0;
    for( auto ke : faceToEdgeMap[ kf ] )
    {
      if( edgeManager.hasNode( ke, nodeID ) )
      {
        edge[count++] = ke;
      }
    }

    if( edge[0] == INT_MAX || edge[1] == INT_MAX )
    {
      GEOSX_ERROR( "SurfaceGenerator::FindFracturePlanes: invalid edge." );
    }


    localFacesToEdges[kf] = std::make_pair( edge[0], edge[1] );

  }


  // now we want to identify the objects on either side of the separation plane. First we assign an array to indicate
  // whether a face/edge is on the fracture plane.

  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
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
  for( localIndex const edgeID : nodeToEdgeMap[ nodeID ] )
  {
    edgeLocations[edgeID] = INT_MIN;
  }

  for( std::set< std::pair< CellElementSubRegion *, localIndex > >::const_iterator k=nodeToElementMaps.begin(); k!=nodeToElementMaps.end(); ++k )
  {
    elemLocations[*k] = INT_MIN;
  }



  /*
     SetLocations( 0, separationPathFaces, faceManager, nodeToElementMaps, localFacesToEdges, //nodeToEdges,
                edgeLocations, faceLocations, elemLocations );

     if( !(SetLocations( 1, separationPathFaces, faceManager, nodeToElementMaps, localFacesToEdges, //nodeToEdges,
                      edgeLocations, faceLocations, elemLocations )) )
     {
     return false;
     }*/

  setLocations( separationPathFaces,
                elemManager,
                faceManager,
                nodeToElementMaps,
                localFacesToEdges,
                edgeLocations,
                faceLocations,
                elemLocations );



  bool fail = false;

  for( localIndex const edgeID : nodeToEdgeMap[ nodeID ] )
  {
    if( edgeLocations[edgeID] == INT_MIN )
    {
      fail = true;
    }
  }
  for( localIndex const kf : nodeToFaceMap[ nodeID ] )
  {
    if( faceLocations[kf] == INT_MIN )
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
     for( std::set< std::pair<CellBlockSubRegion*,localIndex> >::const_iterator k=nodeToElementMaps.begin() ;
        k!=nodeToElementMaps.end() ; ++k )
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

    //    GEOSX_ERROR("SurfaceGenerator::FindFracturePlanes: unset element,face, or edge");
    return false;
  }
  return true;
}

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::setLocations( const std::set< localIndex > & separationPathFaces,
                                     ElementRegionManager & elemManager,
                                     const FaceManager & faceManager,
                                     const std::set< std::pair< CellElementSubRegion *, localIndex > > & nodeToElementMaps,
                                     const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                                     map< localIndex, int > & edgeLocations,
                                     map< localIndex, int > & faceLocations,
                                     map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations )
{
  bool rval = true;
  //  const localIndex separationFace = *(separationPathFaces.begin());

  // insert an element attached to the separation face
  //  std::pair<CellBlockSubRegion*,localIndex> elem0 = m_virtualFaces.m_FaceToElementMap[separationFace][0] ;

  std::pair< CellElementSubRegion *, localIndex > elem0 = *(nodeToElementMaps.begin());


  setElemLocations( 0,
                    elem0,
                    separationPathFaces,
                    elemManager,
                    faceManager,
                    nodeToElementMaps,
                    localFacesToEdges,
                    edgeLocations,
                    faceLocations,
                    elemLocations );

  return rval;
}


//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
bool SurfaceGenerator::setElemLocations( const int location,
                                         const std::pair< CellElementSubRegion *, localIndex > & k,
                                         const std::set< localIndex > & separationPathFaces,
                                         ElementRegionManager & elemManager,
                                         const FaceManager & faceManager,
                                         const std::set< std::pair< CellElementSubRegion *, localIndex > > & nodeToElementMaps,
                                         const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                                         map< localIndex, int > & edgeLocations,
                                         map< localIndex, int > & faceLocations,
                                         map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations )
{
  arrayView1d< localIndex const > const & parentFaceIndices = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  const int otherlocation = (location==0) ? 1 : 0;

  elemLocations[k] = location;


  // loop over all faces on the element
  for( localIndex kf=0; kf<k.first->faceList().size( 1 ); ++kf )
  {

    // define the actual face index, and the virtual face index
    const localIndex faceIndex = k.first->faceList()( k.second, kf );
    const localIndex virtualFaceIndex = ( parentFaceIndices[faceIndex] == -1 ) ?
                                        faceIndex : parentFaceIndices[faceIndex];

    // see if we can find the face in the faceLocations array.
    map< localIndex, int >::iterator iterFace = faceLocations.find( faceIndex );
    // if we can find the face in the faceLocations array, then we must process the face, otherwise it is not
    // connected to the node, so we do nothing.
    if( iterFace != faceLocations.end() )
    {

      if( faceLocations[faceIndex]==otherlocation )
        faceLocations[faceIndex] = -1;
      else if( faceLocations[faceIndex] == INT_MIN )
        faceLocations[faceIndex] = location;

      map< localIndex, std::pair< localIndex, localIndex > >::const_iterator iterF2E = localFacesToEdges.find( faceIndex );

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


        std::pair< CellElementSubRegion *, localIndex > const
        elemIndex0 = { &elemManager.getRegion( er0 ).getSubRegion< CellElementSubRegion >( esr0 ),
                       m_originalFacesToElemIndex[virtualFaceIndex][0] };

        std::pair< CellElementSubRegion *, localIndex > const
        elemIndex1 = { &elemManager.getRegion( er1 ).getSubRegion< CellElementSubRegion >( esr1 ),
                       m_originalFacesToElemIndex[virtualFaceIndex][1] };

        std::pair< CellElementSubRegion *, localIndex > const & nextElem = ( elemIndex0 == k ) ? elemIndex1 : elemIndex0;
        int const nextLocation = (separationPathFaces.count( virtualFaceIndex )==0) ? location : otherlocation;

        // if the first element is the one we are on, and the element is attached
        // to the splitting node, then add the second element to the list.
        if( nodeToElementMaps.find( nextElem )!=nodeToElementMaps.end() )
        {
          if( elemLocations[nextElem]==INT_MIN )
          {
            setElemLocations( nextLocation,
                              nextElem,
                              separationPathFaces,
                              elemManager,
                              faceManager,
                              nodeToElementMaps,
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


//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************
void SurfaceGenerator::performFracture( const localIndex nodeID,
                                        real64 const time_np1,
                                        NodeManager & nodeManager,
                                        EdgeManager & edgeManager,
                                        FaceManager & faceManager,
                                        ElementRegionManager & elementManager,
                                        ModifiedObjectLists & modifiedObjects,
                                        std::vector< std::set< localIndex > > & GEOSX_UNUSED_PARAM( nodesToRupturedFaces ),
                                        std::vector< std::set< localIndex > > & GEOSX_UNUSED_PARAM( edgesToRupturedFaces ),
                                        const std::set< localIndex > & separationPathFaces,
                                        const map< localIndex, int > & edgeLocations,
                                        const map< localIndex, int > & faceLocations,
                                        const map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations )
{
  int const rank = MpiWrapper::commRank( MPI_COMM_WORLD );

  array2d< real64, nodes::REFERENCE_POSITION_PERM > const & X = nodeManager.referencePosition();
  ArrayOfSets< localIndex > & nodeToEdgeMap = nodeManager.edgeList();
  ArrayOfSets< localIndex > & nodeToFaceMap = nodeManager.faceList();
  ArrayOfArrays< localIndex > & nodeToRegionMap = nodeManager.elementRegionList();
  ArrayOfArrays< localIndex > & nodeToSubRegionMap = nodeManager.elementSubRegionList();
  ArrayOfArrays< localIndex > & nodeToElementMap = nodeManager.elementList();

  array2d< localIndex > & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSets< localIndex > & edgeToFaceMap = edgeManager.faceList();

  ArrayOfArrays< localIndex > & faceToNodeMap = faceManager.nodeList();
  ArrayOfArrays< localIndex > & faceToEdgeMap = faceManager.edgeList();
  array2d< localIndex > const & faceToRegionMap = faceManager.elementRegionList();
  array2d< localIndex > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  array2d< localIndex > const & faceToElementMap = faceManager.elementList();

  array1d< integer > const & faceIsExternal = faceManager.isExternal();
  array1d< integer > const & edgeIsExternal = edgeManager.isExternal();
  array1d< integer > const & nodeIsExternal = nodeManager.isExternal();

  SurfaceElementRegion & fractureElementRegion = elementManager.getRegion< SurfaceElementRegion >( "Fracture" );
  array1d< integer > const & isFaceSeparable = faceManager.getExtrinsicData< extrinsicMeshData::IsFaceSeparable >();

  array2d< real64 > const & faceNormals = faceManager.faceNormal();

  array1d< localIndex > const & parentEdgeIndices = edgeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();
  array1d< localIndex > const & childEdgeIndices = edgeManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();
  array1d< localIndex > const & parentNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();
  array1d< localIndex > const & childNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();

  array1d< integer > const & degreeFromCrack = nodeManager.getExtrinsicData< extrinsicMeshData::DegreeFromCrack >();
  array1d< integer > const & nodeDegreeFromCrackTip = nodeManager.getExtrinsicData< extrinsicMeshData::DegreeFromCrackTip >();
  array1d< integer > const & faceDegreeFromCrackTip = faceManager.getExtrinsicData< extrinsicMeshData::DegreeFromCrackTip >();

  array1d< real64 > const & nodeRuptureTime = nodeManager.getExtrinsicData< extrinsicMeshData::RuptureTime >();
  array1d< real64 > const & faceRuptureTime = faceManager.getExtrinsicData< extrinsicMeshData::RuptureTime >();

  // ***** split all the objects first *****

  // Split the node into two, using the original index, and a new one.
  localIndex newNodeIndex;
  if( getLogLevel() )
  {
    GEOSX_LOG_RANK( "" );
    std::cout<<"Splitting node "<<nodeID<<" along separation plane faces: ";
    for( std::set< localIndex >::const_iterator i=separationPathFaces.begin(); i!=separationPathFaces.end(); ++i )
    {
      std::cout<<*i<<", ";
    }
    std::cout<<std::endl;
  }


  nodeManager.splitObject( nodeID, rank, newNodeIndex );

  modifiedObjects.newNodes.insert( newNodeIndex );
  modifiedObjects.modifiedNodes.insert( nodeID );

  nodeToRegionMap.clearArray( newNodeIndex );
  nodeToSubRegionMap.clearArray( newNodeIndex );
  nodeToElementMap.clearArray( newNodeIndex );

  nodeToEdgeMap.clearSet( newNodeIndex );
  nodeToFaceMap.clearSet( newNodeIndex );

  degreeFromCrack[nodeID] = 0;
  degreeFromCrack[newNodeIndex] = 0;
  m_tipNodes.remove( nodeID );
  nodeDegreeFromCrackTip( nodeID ) = 1;
  nodeRuptureTime( nodeID ) = time_np1;
  nodeRuptureTime( newNodeIndex ) = time_np1;

  //TODO HACK...should recalculate mass
//  const real64 newMass = 0.5 * (*nodeManager.m_mass)[nodeID];
//  (*nodeManager.m_mass)[nodeID] = newMass;
//  (*nodeManager.m_mass)[newNodeIndex] = newMass;

  //TODO Either change m_usedFacesForNode to array<std::set> or add insert with iterator to SortedArray
  for( auto const val : separationPathFaces )
  {
    m_usedFacesForNode[nodeID].insert( val );
    m_usedFacesForNode[newNodeIndex].insert( val );
  }

//  SortedArray<localIndex>& usedFacesNew = nodeManager.getReference< array1d<SortedArray<localIndex>>
// >("usedFaces")[newNodeIndex];
//  usedFacesNew = usedFaces[nodeID];


  if( getLogLevel() )
    std::cout<<"Done splitting node "<<nodeID<<" into nodes "<<nodeID<<" and "<<newNodeIndex<<std::endl;

  // split edges
  map< localIndex, localIndex > splitEdges;
  // loop over all edges connected to the node
  for( map< localIndex, int >::const_iterator iter_edge=edgeLocations.begin(); iter_edge!=edgeLocations.end(); ++iter_edge )
  {
    const localIndex & parentEdgeIndex = iter_edge->first;
    const int & location = iter_edge->second;

    // if the edge is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newEdgeIndex;

      edgeManager.splitObject( parentEdgeIndex, rank, newEdgeIndex );

      m_tipEdges.remove( parentEdgeIndex );

      edgeToFaceMap.clearSet( newEdgeIndex );

      if( getLogLevel() )
      {
        GEOSX_LOG_RANK( "" );
        std::cout<<"  Split edge "<<parentEdgeIndex<<" into edges "<<parentEdgeIndex<<" and "<<newEdgeIndex<<std::endl;
      }

      splitEdges[parentEdgeIndex] = newEdgeIndex;
      modifiedObjects.newEdges.insert( newEdgeIndex );
      modifiedObjects.modifiedEdges.insert( parentEdgeIndex );

      for( localIndex const faceIndex : edgeToFaceMap[ parentEdgeIndex ] )
      {
        bool trailingFace = false;
        if( m_trailingFaces.contains( faceIndex ))
        {
          for( localIndex const edgeIndex : faceToEdgeMap[ faceIndex ] )
          {
            if( m_tipEdges.contains( edgeIndex ))
            {
              trailingFace = true;
            }
          }

          if( trailingFace == false )
          {
            m_trailingFaces.remove( faceIndex );
          }
        }
      }

      for( int a=0; a<2; ++a )
      {
        edgeManager.nodeList( newEdgeIndex, a ) = edgeManager.nodeList( parentEdgeIndex, a );
      }

    } //    if( location == -1  )
  } // for( map<localIndex,int>::const_iterator iter_edge...


  // split the faces
  array1d< integer > const & ruptureState = faceManager.getExtrinsicData< extrinsicMeshData::RuptureState >();
  map< localIndex, localIndex > splitFaces;


  SortedArray< localIndex > & externalFaces = faceManager.externalSet();

  // loop over all faces attached to the nodeID
  for( map< localIndex, int >::const_iterator iter_face=faceLocations.begin(); iter_face!=faceLocations.end(); ++iter_face )
  {
    const localIndex faceIndex = iter_face->first;
//    localIndex const parentFaceIndex = parentFaceIndices[faceIndex]==faceIndex ? faceIndex :
// parentFaceIndices[faceIndex];
    const int location = iter_face->second;
    // if the face is on the separation plane, then split it
    if( location == -1 )
    {
      localIndex newFaceIndex;

      if( faceManager.splitObject( faceIndex, rank, newFaceIndex ) )
      {

        if( getLogLevel() )
        {
          GEOSX_LOG_RANK( "" );
          std::cout<<"  Split face "<<faceIndex<<" into faces "<<faceIndex<<" and "<<newFaceIndex<<std::endl;
        }

        splitFaces[faceIndex] = newFaceIndex;
        modifiedObjects.newFaces.insert( newFaceIndex );
        modifiedObjects.modifiedFaces.insert( faceIndex );

        ruptureState[faceIndex] = 2;
        ruptureState[newFaceIndex] = 2;

        faceRuptureTime( faceIndex ) = time_np1;
        faceRuptureTime( newFaceIndex ) = time_np1;


        m_trailingFaces.insert( faceIndex );
        m_tipFaces.remove( faceIndex );
        faceDegreeFromCrackTip( faceIndex ) = 0;
        faceDegreeFromCrackTip( newFaceIndex ) = 0;

        localIndex const numFaceEdges = faceToEdgeMap.sizeOfArray( faceIndex );
        faceToEdgeMap.resizeArray( newFaceIndex, numFaceEdges );
        for( localIndex a = 0; a < numFaceEdges; ++a )
        {
          faceToEdgeMap( newFaceIndex, a ) = faceToEdgeMap( faceIndex, a );
        }

        localIndex const numFaceNodes = faceToNodeMap.sizeOfArray( faceIndex );
        faceToNodeMap.resizeArray( newFaceIndex, numFaceNodes );
        for( localIndex a=0; a<numFaceNodes; ++a )
        {
          localIndex const aa = a == 0 ? a : numFaceNodes - a;
          faceToNodeMap( newFaceIndex, aa ) = faceToNodeMap( faceIndex, a );
        }
        LvArray::tensorOps::scale< 3 >( faceNormals[ newFaceIndex ], -1 );

        externalFaces.insert( newFaceIndex );
        externalFaces.insert( faceIndex );


        // Fu: All edges of the parent face should be external now.
        // We have to do the following because isExternal attribute of the tip edge is not handled by the splitter.
        for( localIndex const edgeIndex : faceManager.edgeList()[ faceIndex ] )
        {
          if( parentEdgeIndices[edgeIndex]==-1 && childEdgeIndices[edgeIndex]==-1 )
          {
            m_tipEdges.insert( edgeIndex );

            for( localIndex const iface: edgeManager.faceList()[ edgeIndex ] )
            {
              if( faceToElementMap.size( 1 ) == 2  &&
                  faceIsExternal[iface] < 1 &&
                  checkOrphanElement( elementManager, faceManager, iface ) == 0 &&
                  isFaceSeparable[iface] == 1
//                  && fabs(LvArray::tensorOps::AiBi< 3 >(faceNormals[faceIndex], faceNormals[iface])) > cos( m_maxTurnAngle )
                  )
              {
                m_tipFaces.insert( iface );
              }
            }
          }
          if( edgeIsExternal[edgeIndex]==0 )
          {
            edgeIsExternal[edgeIndex] = 2;
          }
        }
        for( localIndex const nodeIndex : faceToNodeMap[ faceIndex ] )
        {
          if( parentNodeIndices[nodeIndex]==-1 && childNodeIndices[nodeIndex]==-1 )
          {
            m_tipNodes.insert( nodeIndex );
            nodeDegreeFromCrackTip( nodeIndex ) = 0;
          }
          if( nodeIsExternal[nodeIndex] )
          {
            nodeIsExternal[nodeIndex] = 2;
          }
        }

        {
          localIndex faceIndices[2] = {faceIndex, newFaceIndex};
          localIndex const
          newFaceElement = fractureElementRegion.addToFractureMesh( time_np1,
                                                                    &edgeManager,
                                                                    &faceManager,
                                                                    this->m_originalFaceToEdges.toViewConst(),
                                                                    "faceElementSubRegion",
                                                                    faceIndices );
          m_faceElemsRupturedThisSolve.insert( newFaceElement );
          modifiedObjects.newElements[ {fractureElementRegion.getIndexInParent(), 0} ].insert( newFaceElement );
        }
//        externalFaceManager.SplitFace(parentFaceIndex, newFaceIndex, nodeManager);

      } // if( faceManager.SplitObject( faceIndex, newFaceIndex ) )
    } // if( location == -1 )
  } // for( map<localIndex,int>::const_iterator iter_face


  // ***** now correct all the relations between the objects *****

  /* To accomplish this annoying yet exceedingly important task, we will take a "top down"
   * approach. Note that this is a two way correction, i.e. if we are correcting
   * elementToNodes, we also correct nodeToElementMaps. This is summarized as:
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

  array1d< localIndex > const & parentFaceIndex = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();
  array1d< localIndex > const & childFaceIndex = faceManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();

  // 1) loop over all elements attached to the nodeID
  for( map< std::pair< CellElementSubRegion *, localIndex >, int >::const_iterator iter_elem =
         elemLocations.begin(); iter_elem != elemLocations.end(); ++iter_elem )
  {
    const int & location = iter_elem->second;

    if( location==1 )
    {
      const std::pair< CellElementSubRegion *, localIndex > & elem = iter_elem->first;

      CellElementSubRegion & elemSubRegion = *(elem.first);
      ElementRegionBase & elemRegion = dynamicCast< ElementRegionBase & >( elemSubRegion.getParent().getParent() );
      string const & elemRegionName = elemRegion.getName();

      localIndex const regionIndex = elementManager.getRegions().getIndex( elemRegionName );
      localIndex const subRegionIndex = elemRegion.getSubRegions().getIndex( elemSubRegion.getName() );
      const localIndex elemIndex = elem.second;

      modifiedObjects.modifiedElements[{regionIndex, subRegionIndex}].insert( elemIndex );


      array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elemsToNodes = elemSubRegion.nodeList();
      array2d< localIndex > const & elemsToFaces = elemSubRegion.faceList();

      if( getLogLevel() > 1 )
        std::cout<<"Element "<<elemIndex<<std::endl;

      // 2a) correct elementToNode and nodeToElement
      if( getLogLevel() > 1 )
        std::cout<<"  Looping over all nodes on element, and correcting node<->element maps:"<<std::endl;


      real64 elemCenter[3] = {0.0, 0.0, 0.0};
      {
        // loop over all nodes on element
        if( getLogLevel() > 1 )
          std::cout<<"    m_ElementToNodeMap = ( ";
        for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
        {
          LvArray::tensorOps::add< 3 >( elemCenter, X[ elemsToNodes[elemIndex][a] ] );
          // if the node was just split
          if( elemsToNodes[elemIndex][a] == nodeID )
          {

            if( getLogLevel() > 1 )
              std::cout<<elemsToNodes[elemIndex][a]<<"->"<<newNodeIndex<<", ";

            elemsToNodes[elemIndex][a] = newNodeIndex;

            insert( nodeManager.toElementRelation(), newNodeIndex, regionIndex, subRegionIndex, elemIndex );
            erase( nodeManager.toElementRelation(), nodeID, regionIndex, subRegionIndex, elemIndex );
          }
          else if( getLogLevel() > 1 )
            std::cout<<elemsToNodes[elemIndex][a]<<", ";
        }
        LvArray::tensorOps::scale< 3 >( elemCenter, 1.0 / elemsToNodes.size( 1 ) );
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;

        if( getLogLevel() > 1 )
        {
          for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
          {
            if( getLogLevel() > 1 )
            {
              std::cout<<"    nodeToElemMaps["<<elemsToNodes[elemIndex][a]<<"] = ( ";
              for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( elemsToNodes[elemIndex][a] ); ++k )
              {
                std::cout<<"["<<nodeToRegionMap[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodeToSubRegionMap[elemsToNodes[elemIndex][a]][k]<<","
                         <<nodeToElementMap[elemsToNodes[elemIndex][a]][k]<<"] , ";
              }
              std::cout<<" )"<<std::endl;
            }
          }

          if( getLogLevel() > 1 )
          {
            std::cout<<"    nodeToElemMaps["<<nodeID<<"] = ( ";
            for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
            {
              std::cout<<"["<<nodeToRegionMap[nodeID][k]<<","
                       <<nodeToSubRegionMap[nodeID][k]<<","
                       <<nodeToElementMap[nodeID][k]<<"] , ";
            }
            std::cout<<" )"<<std::endl;
          }
        }
      }



      // 2b) loop over all faces on element.
      if( getLogLevel() > 1 )
      {
        std::cout<<"  Looping over all faces on element (parent and child):"<<std::endl;
      }

      // we need to build a list of faces that is elemToFaces FOLLOWED by any
      // parent face of those indicated in elemToFaces

      // Now we do a loop over the facelist and process all the faces
      for( int kf=0; kf<elemSubRegion.numFacesPerElement(); ++kf )
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
//          faceManager.m_toElements[newFaceIndex].emplace_back( elem );

          faceToRegionMap[newFaceIndex][0] = regionIndex;
          faceToSubRegionMap[newFaceIndex][0] = subRegionIndex;
          faceToElementMap[newFaceIndex][0] = elemIndex;
          faceToRegionMap[newFaceIndex][1] = -1;
          faceToSubRegionMap[newFaceIndex][1] = -1;
          faceToElementMap[newFaceIndex][1] = -1;

          // remove the element from the parent face
          if( faceToRegionMap[faceIndex][0] == regionIndex &&
              faceToSubRegionMap[faceIndex][0] == subRegionIndex &&
              faceToElementMap[faceIndex][0] == elemIndex )
          {
            faceToRegionMap[faceIndex][0] = faceToRegionMap[faceIndex][1];
            faceToSubRegionMap[faceIndex][0] = faceToSubRegionMap[faceIndex][1];
            faceToElementMap[faceIndex][0] = faceToElementMap[faceIndex][1];
            faceToRegionMap[faceIndex][1] = -1;
            faceToSubRegionMap[faceIndex][1] = -1;
            faceToElementMap[faceIndex][1] = -1;
          }
          else if( faceToRegionMap[faceIndex][1] == regionIndex &&
                   faceToSubRegionMap[faceIndex][1] == subRegionIndex &&
                   faceToElementMap[faceIndex][1] == elemIndex )
          {
            faceToRegionMap[faceIndex][1] = -1;
            faceToSubRegionMap[faceIndex][1] = -1;
            faceToElementMap[faceIndex][1] = -1;
          }

          if( getLogLevel() > 1 )
          {
            GEOSX_LOG( "    faceToRegionMap["<<newFaceIndex<<"][0]    = "<<faceToRegionMap[newFaceIndex][0] );
            GEOSX_LOG( "    faceToSubRegionMap["<<newFaceIndex<<"][0] = "<<faceToSubRegionMap[newFaceIndex][0] );
            GEOSX_LOG( "    faceToElementMap["<<newFaceIndex<<"][0]      = "<<faceToElementMap[newFaceIndex][0] );
            GEOSX_LOG( "    faceToRegionMap["<<newFaceIndex<<"][1]    = "<<faceToRegionMap[newFaceIndex][1] );
            GEOSX_LOG( "    faceToSubRegionMap["<<newFaceIndex<<"][1] = "<<faceToSubRegionMap[newFaceIndex][1] );
            GEOSX_LOG( "    faceToElementMap["<<newFaceIndex<<"][1]      = "<<faceToElementMap[newFaceIndex][1] );

            GEOSX_LOG( "    faceToRegionMap["<<faceIndex<<"][0]    = "<<faceToRegionMap[faceIndex][0] );
            GEOSX_LOG( "    faceToSubRegionMap["<<faceIndex<<"][0] = "<<faceToSubRegionMap[faceIndex][0] );
            GEOSX_LOG( "    faceToElementMap["<<faceIndex<<"][0]      = "<<faceToElementMap[faceIndex][0] );
            GEOSX_LOG( "    faceToRegionMap["<<faceIndex<<"][1]    = "<<faceToRegionMap[faceIndex][1] );
            GEOSX_LOG( "    faceToSubRegionMap["<<faceIndex<<"][1] = "<<faceToSubRegionMap[faceIndex][1] );
            GEOSX_LOG( "    faceToElementMap["<<faceIndex<<"][1]      = "<<faceToElementMap[faceIndex][1] );

          }

          for( int i = 0; i < 2; i++ )
          {
            localIndex iFace = i == 0 ? faceIndex : newFaceIndex;

            localIndex elementIndex = faceToElementMap[iFace][0];
            CellElementSubRegion & elementSubRegion = elementManager.getRegion( faceToRegionMap[iFace][0] ).
                                                        getSubRegion< CellElementSubRegion >( faceToSubRegionMap[iFace][0] );
            arrayView2d< real64 const > const subRegionElemCenter = elementSubRegion.getElementCenter();

            faceManager.sortFaceNodes( X, subRegionElemCenter[ elementIndex ], faceToNodeMap[ iFace ], faceToNodeMap.sizeOfArray( iFace ) );

            //Face normal need to be updated here
            real64 fCenter[ 3 ];
            computationalGeometry::Centroid_3DPolygon( faceToNodeMap[ iFace ],
                                                       X,
                                                       fCenter,
                                                       faceNormals[ iFace ] );
          }

        } // if( splitFaces.count( faceID ) > 0 )

        modifiedObjects.modifiedFaces.insert( faceIndex );



        // 3b) correct faceToNodes and nodeToFaces

        if( getLogLevel() > 1 )
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
        for( localIndex & nodeIndex : faceToNodeMap[ newFaceIndex ] )
        {
          if( getLogLevel() > 1 )
            std::cout<<nodeIndex;

          // if the facenode is the one that is being split
          if( nodeIndex == nodeID )
          {
            nodeIndex = newNodeIndex;

            // if it is not a new face.
            if( !isNewFace )
            {
              // remove the face from the nodeToFaceMap of the parent node.
              nodeToFaceMap.removeFromSet( nodeID, faceIndex );

              // add the face to the nodeToFaceMap of the new node.
              nodeToFaceMap.insertIntoSet( nodeIndex, faceIndex );
            }
            else
            {
              // it is a new face

              // insert the newFace into the nodeToFaceMap of the newNode
              nodeToFaceMap.insertIntoSet( nodeIndex, newFaceIndex );
            }
            if( getLogLevel() > 1 )
              std::cout<<"->"<<nodeIndex<<", ";
          }
          else // the node is not being split
          {
            nodeToFaceMap.insertIntoSet( nodeIndex, newFaceIndex );

            if( getLogLevel() > 1 )
              std::cout<<", ";
          }

        }
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;



        // faceToEdges
        if( getLogLevel() > 1 )
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
        for( localIndex & edgeIndex : faceToEdgeMap[ newFaceIndex ] )
        {

          // if the edge was just split
          if( splitEdges.count( edgeIndex ) > 0 )
          {
            if( faceIndex == newFaceIndex )
            {
              edgeToFaceMap.removeFromSet( edgeIndex, faceIndex );
            }

            edgeIndex = splitEdges[edgeIndex];
          }
          edgeToFaceMap.insertIntoSet( edgeIndex, newFaceIndex );

          modifiedObjects.modifiedEdges.insert( edgeIndex );

          if( getLogLevel() > 1 )
            std::cout<<edgeIndex;



          //edgeToNodeMap
          if( getLogLevel() > 1 )
          {
            std::cout<<"(";
          }

          {
            for( localIndex a=0; a<edgeToNodeMap.size( 1 ); ++a )
            {
              if( edgeToNodeMap[edgeIndex][a] == nodeID )
              {

                if( getLogLevel() > 1 )
                  std::cout<<edgeToNodeMap[edgeIndex][a];

                edgeToNodeMap[edgeIndex][a] = newNodeIndex;
                nodeToEdgeMap.removeFromSet( nodeID, edgeIndex );

                if( getLogLevel() > 1 )
                  std::cout<<"->"<<edgeToNodeMap[edgeIndex][a]<<", ";

              }
              else if( getLogLevel() > 1 )
                std::cout<<edgeToNodeMap[edgeIndex][a]<<", ";

              nodeToEdgeMap.insertIntoSet( edgeToNodeMap[edgeIndex][a], edgeIndex );
              modifiedObjects.modifiedNodes.insert( edgeToNodeMap[edgeIndex][a] );
            }
            if( getLogLevel() > 1 )
              std::cout<<")";
          }
          if( getLogLevel() > 1 )
            std::cout<<", ";
        }
        if( getLogLevel() > 1 )
          std::cout<<")"<<std::endl;
      } // for( int kf=0 ; kf<elemRegion.m_numFacesPerElement ; ++kf )
    } // if( location==1 )
  } // for( map<std::pair<CellBlockSubRegion*, localIndex>, int>::const_iterator iter_elem = elemLocations.begin()
}


void SurfaceGenerator::mapConsistencyCheck( localIndex const GEOSX_UNUSED_PARAM( nodeID ),
                                            NodeManager const & nodeManager,
                                            EdgeManager const & edgeManager,
                                            FaceManager const & faceManager,
                                            ElementRegionManager const & elementManager,
                                            map< std::pair< CellElementSubRegion *, localIndex >, int > const & elemLocations )
{
  //**************************************************************************
  // THIS IS ALL JUST CONSISTENCY CHECKING
  //**************************************************************************


  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList();


  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();


#if 1
  if( getLogLevel() > 2 )
  {
    std::cout<<"CONSISTENCY CHECKING OF THE MAPS"<<std::endl;

    for( map< std::pair< CellElementSubRegion *, localIndex >, int >::const_iterator iter_elem=elemLocations.begin();
         iter_elem!=elemLocations.end(); ++iter_elem )
    {
      const std::pair< CellElementSubRegion *, localIndex > & elem = iter_elem->first;

      CellElementSubRegion & elemSubRegion = *(elem.first);
      const localIndex elemIndex = elem.second;

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elemSubRegion.nodeList();
      arrayView2d< localIndex > & elemsToFaces = elemSubRegion.faceList();


      std::set< localIndex > elemNodes;


      GEOSX_LOG( "Element " << elemIndex );
      std::cout << " elementToNodes = ";
      for( int a=0; a<8; ++a )
      {
        elemNodes.insert( elemsToNodes( elemIndex, a ));
        std::cout << elemsToNodes( elemIndex, a )<<", ";
      }
      std::cout << std::endl;

      std::cout << " elementToFaces->edges->nodes = ";


      // Now we do a loop over the facelist and process all the faces
      for( int kf=0; kf<elemSubRegion.numFacesPerElement(); ++kf )
      {
        std::set< localIndex > faceNodes;

        localIndex faceIndex  = elemsToFaces( elemIndex, kf );

        std::cout << "                              = ";
        std::cout << faceIndex << "( ";
        for( int b=0; b<4; ++b )
        {
          localIndex faceNodeID = faceToNodeMap( faceIndex, b );
          faceNodes.insert( faceNodeID );
          if( elemNodes.count( faceNodeID ) == 0 )
          {
            std::cout << "*";
          }
          std::cout << faceNodeID << ",";
        }
        std::cout << " )      ";



        std::cout << faceIndex << "[ ";
        for( int b=0; b<4; ++b )
        {
          localIndex edgeIndex = faceToEdgeMap( faceIndex, b );
          std::cout << edgeIndex << "( ";
          for( int c=0; c<2; ++c )
          {
            localIndex edgeNodeID = edgeToNodeMap( edgeIndex, c );
            if( elemNodes.count( edgeNodeID ) == 0  && kf<elemSubRegion.numFacesPerElement() )
            {
              std::cout << "*";
            }
            if( faceNodes.count( edgeNodeID ) == 0 )
            {
              std::cout << "#";
            }
            std::cout << edgeNodeID << ",";
          }
          std::cout << " ), ";
        }
        std::cout << " ] \n";

      }
      std::cout << std::endl;

    }

  }

  if( getLogLevel() > 2 )
  {
    // nodeToEdge
    std::vector< std::set< localIndex > > inverseEdgesToNodes( nodeManager.size() );

    for( localIndex ke=0; ke<edgeManager.size(); ++ke )
    {
      for( localIndex b= 0; b<edgeToNodeMap.size( 1 ); ++b )
      {
        localIndex nodeIndex = edgeToNodeMap( ke, b );
        inverseEdgesToNodes[nodeIndex].insert( ke );
      }
    }
    std::cout << "Check NodeToEdge:  nodeToEdgeMap  inverseEdgesToNodes" << std::endl;
    for( localIndex a=0; a<nodeManager.size(); ++a )
    {
      std::cout << "nodeToEdgeMap[" << a << "] = ( ";

      for( localIndex const edgeID : nodeToEdgeMap[ a ] )
      {
        if( inverseEdgesToNodes[a].count( edgeID ) == 0 )
        {
          std::cout << "*";
        }
        std::cout << edgeID << ", ";
      }

      std::cout<<")    (";

      for( localIndex const edgeID : inverseEdgesToNodes[a] )
      {
        if( !nodeToEdgeMap.contains( a, edgeID ) )
          std::cout << "*";
        std::cout << edgeID <<", ";
      }
      std::cout<< ")" <<std::endl;
    }
  }

  if( getLogLevel() > 2 )
  {
    // nodeToFace
    std::vector< std::set< localIndex > > inverseFacesToNodes( nodeManager.size() );
    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      for( localIndex const b : faceToNodeMap[ kf ] )
      {
        inverseFacesToNodes[b].insert( kf );
      }
    }
    std::cout << "Check NodeToFace:  nodeToFaceMap  inverseFacesToNodes" << std::endl;
    for( localIndex a=0; a<nodeManager.size(); ++a )
    {
      std::cout << "m_nodeToFaceMap[ "<< a << "] = ( ";
      for( localIndex const & faceID : nodeToFaceMap[ a ] )
      {
        if( inverseFacesToNodes[a].count( faceID ) == 0 )
          std::cout << "*";
        std::cout << faceID << ", ";
      }
      std::cout<<")    (";

      for( localIndex const edgeID : inverseFacesToNodes[a] )
      {
        if( !nodeToFaceMap.contains( a, edgeID ) )
          std::cout << "*";
        std::cout << edgeID << ", ";
      }
      std::cout<<")"<<std::endl;
    }
  }



  if( getLogLevel() > 2 )
  {


    // nodeToElement
    std::vector< std::set< std::pair< CellElementSubRegion const *, localIndex > > > inverseElemsToNodes( nodeManager.size() );
    elementManager.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = subRegion.nodeList();
      for( localIndex k=0; k<subRegion.size(); ++k )
      {
        for( localIndex a=0; a<elemsToNodes.size( 1 ); ++a )
        {
          inverseElemsToNodes[elemsToNodes( k, a )].emplace( &subRegion, k );
        }
      }
    } );

    std::cout<<"Check NodeToElem: nodesToElems  inverseElemsToNodes "<<std::endl;


    for( localIndex a=0; a<nodeManager.size(); ++a )
    {

      std::set< std::pair< CellElementSubRegion const *, localIndex > > nodeToElements;
      for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( a ); ++k )
      {
        if( nodeToRegionMap[a][k]!=-1 && nodeToSubRegionMap[a][k]!=-1 && nodeToElementMap[a][k]!=-1 )
        {
          nodeToElements.emplace( &elementManager.getRegion( nodeToRegionMap( a, k ) ).
                                    getSubRegion< CellElementSubRegion >( nodeToSubRegionMap( a, k ) ),
                                  nodeToElementMap( a, k ) );
        }
      }


      std::cout<<"m_NodeToElementMap["<<a<<"] = ( ";
      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::iterator
           ielem=nodeToElements.begin(); ielem!=nodeToElements.end(); ++ielem )
      {
        if( inverseElemsToNodes[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=inverseElemsToNodes[a].begin();
           ielem!=inverseElemsToNodes[a].end(); ++ielem )
      {
        if( nodeToElements.count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")"<<std::endl;
    }


    // edgeToFace
    std::vector< std::set< localIndex > > inverseFacesToEdges( edgeManager.size() );
    for( localIndex kf=0; kf<faceManager.size(); ++kf )
    {
      for( localIndex const b : faceToEdgeMap[ kf ] )
      {
        inverseFacesToEdges[ b ].insert( kf );
      }
    }
    std::cout<<"Check EdgeToFace: edgeToFaceMap  inverseFacesToEdges "<<std::endl;
    for( localIndex ke=0; ke<edgeManager.size(); ++ke )
    {
      std::cout<<"m_edgeToFaceMap["<<ke<<"] = ( ";
      for( localIndex const faceID : edgeManager.faceList()[ ke ] )
      {
        if( inverseFacesToEdges[ke].count( faceID ) == 0 )
          std::cout << "*";
        std::cout<<faceID<<", ";
      }
      std::cout<<")    (";

      for( std::set< localIndex >::const_iterator iface=inverseFacesToEdges[ke].begin();
           iface!=inverseFacesToEdges[ke].end(); ++iface )
      {
        if( !edgeManager.faceList().contains( ke, *iface ) )
          std::cout<<"*";
        std::cout<< *iface <<", ";
      }
      std::cout<<")"<<std::endl;
    }

    // faceToElement
    std::vector< std::set< std::pair< CellElementSubRegion const *, localIndex > > > inverseElemsToFaces( faceManager.size() );
    elementManager.forElementSubRegions< CellElementSubRegion >( [&] ( CellElementSubRegion const & subRegion )
    {
      arrayView2d< localIndex > const & elemsToFaces = subRegion.faceList();

      for( localIndex k=0; k<subRegion.size(); ++k )
      {
        for( localIndex a=0; a<elemsToFaces.size( 1 ); ++a )
        {
          const localIndex faceID = elemsToFaces( k, a );
          inverseElemsToFaces[ faceID ].emplace( &subRegion, k );

          //            if( parentFaceIndex[faceID] != -1 )
          //            {
          //              inverseElemsToFaces[parentFaceIndex[faceID]].insert(elem);
          //            }
        }
      }
    } );

    std::cout<<"Check FacesToElem: facesToElems  inverseElemsToFaces "<<std::endl;
    for( localIndex a=0; a<faceManager.size(); ++a )
    {

      std::vector< std::pair< CellElementSubRegion const *, localIndex > > faceToElements;
      for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
      {
        // TODO This only works for a single region
        if( faceToRegionMap( a, k ) != -1 )
        {
          faceToElements.emplace_back( &elementManager.getRegion( faceToRegionMap( a, k ) ).
                                         getSubRegion< CellElementSubRegion >( faceToSubRegionMap( a, k ) ),
                                       faceToElementMap( a, k ) );
        }
      }


      std::cout<<"m_FaceToElementMap["<<a<<"] = ( ";

      for( std::vector< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator
           ielem=faceToElements.begin();
           ielem!=faceToElements.end(); ++ielem )
      {
        if( inverseElemsToFaces[a].count( *ielem ) == 0 )
          std::cout<<"*";

        std::cout<<ielem->second<<", ";
      }
      std::cout<<")    (";

      for( std::set< std::pair< CellElementSubRegion const *, localIndex > >::const_iterator ielem=inverseElemsToFaces[a].begin();
           ielem!=inverseElemsToFaces[a].end(); ++ielem )
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



real64 SurfaceGenerator::calculateKinkAngle ( const localIndex edgeID,
                                              const NodeManager & GEOSX_UNUSED_PARAM( nodeManager ),
                                              EdgeManager & edgeManager,
                                              FaceManager & faceManager )
{
  // TODO: This method should be re-implemented.
  localIndex_array faces;
  // real64 kinkAngle;

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  for( localIndex const iface : edgeManager.faceList()[ edgeID ] )
  {
    if( faceIsExternal[iface] == 1 )
      faces.emplace_back( iface );
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
//    kinkAngle = acos(LvArray::tensorOps::AiBi< 3 >(vecFace[0],vecFace[1])*0.999999) / 3.141592653589793238462 * 180.0;
//
//    R1Tensor vecFaceNorm;
//    vecFaceNorm = faceManager.FaceNormal(nodeManager, faces[0]);
//    vecFaceNorm  += faceManager.FaceNormal(nodeManager, faces[1]);
//    vecFaceNorm /= 2.0;
//
//    if (LvArray::tensorOps::AiBi< 3 >(vecFace[2], vecFaceNorm) < 0.0)
//      kinkAngle = 360.0 - kinkAngle;
//
//    return(kinkAngle);
//
//  }
    return 1e100;
}

void SurfaceGenerator::calculateKinkAngles ( FaceManager & faceManager,
                                             EdgeManager & edgeManager,
                                             NodeManager & nodeManager,
                                             ModifiedObjectLists & modifiedObjects,
                                             const bool prefrac )
{
  arrayView1d< real64 > & kinkAngle = edgeManager.getReference< real64_array >( "kinkAngle" );

  if( prefrac )
  {
    for( localIndex edgeID = 0; edgeID < edgeManager.size(); ++edgeID )
    {
      kinkAngle[edgeID] = calculateKinkAngle( edgeID, nodeManager, edgeManager, faceManager );
    }
  }
  else
  {
    for( std::set< localIndex >::const_iterator i=modifiedObjects.newEdges.begin(); i!=modifiedObjects.newEdges.end(); ++i )
    {
      kinkAngle[*i] = calculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
    for( std::set< localIndex >::const_iterator i=modifiedObjects.modifiedEdges.begin(); i!=modifiedObjects.modifiedEdges.end(); ++i )
    {
      kinkAngle[*i] = calculateKinkAngle( *i, nodeManager, edgeManager, faceManager );
    }
  }
}



void SurfaceGenerator::identifyRupturedFaces( DomainPartition & domain,
                                              NodeManager & nodeManager,
                                              EdgeManager & edgeManager,
                                              FaceManager & faceManager,
                                              ElementRegionManager & elementManager,
                                              const bool prefrac )
{
  arrayView1d< integer > const & isEdgeGhost = edgeManager.ghostRank();
  arrayView1d< real64 const > const & SIFNode = nodeManager.getExtrinsicData< extrinsicMeshData::SIFNode >();


  // We use the color map scheme because we can mark a face to be rupture ready from a partition where the face is a
  // ghost.

  if( !m_nodeBasedSIF )
  {
    {
//      for( int color=0 ; color<partition.NumColor() ; ++color )
      {
        ModifiedObjectLists modifiedObjects;
//        if( partition.Color() == color )
        {
          for( localIndex iEdge = 0; iEdge != edgeManager.size(); ++iEdge )
          {

            if( isEdgeGhost[iEdge] < 0 )
            {
              int edgeMode = checkEdgeSplitability( iEdge,
                                                    nodeManager,
                                                    faceManager,
                                                    edgeManager,
                                                    prefrac );
              if( edgeMode == 0 || edgeMode == 1 ) // We need to calculate SIF
              {
                real64 vecTipNorm[3], vecTip[3];
                localIndex trailFaceID = 0;
                real64 const SIF = calculateEdgeSif( domain, iEdge, trailFaceID,
                                                     nodeManager,
                                                     edgeManager,
                                                     faceManager,
                                                     elementManager,
                                                     vecTipNorm,
                                                     vecTip );

                if( SIF >  minimumToughnessOnEdge( iEdge, nodeManager, edgeManager, faceManager ) * 0.5 ) // && edgeMode
                                                                                                          // == 1)
                {
                  markRuptureFaceFromEdge( iEdge, trailFaceID,
                                           nodeManager,
                                           edgeManager,
                                           faceManager,
                                           elementManager,
                                           vecTipNorm,
                                           vecTip,
                                           modifiedObjects,
                                           edgeMode );
                }
              }
            }
          }
        }
      }
    }
  }
  else
  {
    ModifiedObjectLists modifiedObjects;

    calculateNodeAndFaceSif( domain, nodeManager, edgeManager, faceManager, elementManager );

    for( auto nodeIndex: m_tipNodes )
    {
      if( SIFNode[nodeIndex] > minimumToughnessOnNode( nodeIndex, nodeManager, edgeManager, faceManager ))
      {
        markRuptureFaceFromNode( nodeIndex,
                                 nodeManager,
                                 edgeManager,
                                 faceManager,
                                 elementManager,
                                 modifiedObjects );
      }
    }
  }



}

void SurfaceGenerator::calculateNodeAndFaceSif( DomainPartition & domain,
                                                NodeManager & nodeManager,
                                                EdgeManager & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager & elementManager )
{
  arrayView1d< real64 > const & SIFNode = nodeManager.getExtrinsicData< extrinsicMeshData::SIFNode >();
  arrayView1d< real64 > const & SIFonFace = faceManager.getExtrinsicData< extrinsicMeshData::SIFonFace >();

  std::vector< std::vector< real64 > > SIFNode_All, SIFonFace_All;
  std::vector< real64 > SIFOnEdge;
  SIFNode_All.resize( nodeManager.size() );
  SIFonFace_All.resize( faceManager.size() );
  SIFOnEdge.resize( edgeManager.size() );


  SIFNode.setValues< parallelHostPolicy >( 0 );
  SIFonFace.setValues< parallelHostPolicy >( 0 );

  arrayView2d< real64 const > const &
  fext = nodeManager.getReference< array2d< real64 > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternalString() );
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & displacement = nodeManager.totalDisplacement();
  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  const arrayView1d< integer > & isNodeGhost = nodeManager.ghostRank();

  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  arrayView1d< localIndex const > const & childFaceIndices = faceManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();
  arrayView1d< localIndex const > const & childNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ChildIndex >();
  arrayView1d< localIndex > const & parentNodeIndices = nodeManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  ConstitutiveBase const & solid = cm.getConstitutiveRelation< ConstitutiveBase >( m_solidMaterialNames[0] );
  m_solidMaterialFullIndex = solid.getIndexInParent();

  ConstitutiveManager & constitutiveManager =
    domain.getGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "shearModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "bulkModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const, solid::STRESS_USD > > const
  stress = elementManager.constructFullMaterialViewAccessor< array3d< real64, solid::STRESS_PERMUTATION >,
                                                             arrayView3d< real64 const, solid::STRESS_USD > >( SolidBase::viewKeyStruct::stressString(),
                                                                                                               constitutiveManager );


  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > const
  dNdX = elementManager.constructViewAccessor< array4d< real64 >, arrayView4d< real64 const > >( keys::dNdX );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const
  detJ = elementManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( keys::detJ );



  nodeManager.totalDisplacement().move( LvArray::MemorySpace::CPU, false );
  elementManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
  {
    for( localIndex mat=0; mat<m_solidMaterialNames.size(); ++mat )
    {
      subRegion.getConstitutiveModel( m_solidMaterialNames[mat] ).
        getReference< array3d< real64, solid::STRESS_PERMUTATION > >( SolidBase::viewKeyStruct::stressString() ).move( LvArray::MemorySpace::CPU,
                                                                                                                       false );
    }
  } );
  displacement.move( LvArray::MemorySpace::CPU, false );


  for( localIndex const trailingFaceIndex : m_trailingFaces )
//  RAJA::forall< parallelHostPolicy >( RAJA::TypedRangeSegment< localIndex >( 0, m_trailingFaces.size() ),
//                                      [=] GEOSX_HOST_DEVICE ( localIndex const trailingFacesCounter )
  {
//    localIndex const trailingFaceIndex = m_trailingFaces[ trailingFacesCounter ];

    real64 const faceNormalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[trailingFaceIndex] );
    //TODO: check if a ghost face still has the correct attributes such as normal vector, face center, face index.
    localIndex_array unpinchedNodeID;
    localIndex_array pinchedNodeID;
    localIndex_array tipEdgesID;

    for( localIndex const nodeIndex : faceToNodeMap[ trailingFaceIndex ] )
    {
      if( m_tipNodes.contains( nodeIndex ))
      {
        pinchedNodeID.emplace_back( nodeIndex );
      }
      else
      {
        unpinchedNodeID.emplace_back( nodeIndex );
      }
    }

    for( localIndex const edgeIndex : faceToEdgeMap[ trailingFaceIndex ] )
    {
      if( m_tipEdges.contains( edgeIndex ))
      {
        tipEdgesID.emplace_back( edgeIndex );
      }
    }

    if( tipEdgesID.size() >= 1 )
    {
      for( localIndex const nodeIndex : pinchedNodeID )
      {
        if( isNodeGhost[nodeIndex] < 0 )
        {
          real64 nodeDisconnectForce[3] = { 0 };
          real64 const nodePosition[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[nodeIndex] );
          localIndex tralingNodeID = std::numeric_limits< localIndex >::max();
          localIndex nElemEachSide[2];
          nElemEachSide[0] = 0;
          nElemEachSide[1] = 0;

          for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeIndex ); ++k )
          {
            localIndex const er  = nodeToRegionMap[nodeIndex][k];
            localIndex const esr = nodeToSubRegionMap[nodeIndex][k];
            localIndex const ei  = nodeToElementMap[nodeIndex][k];

            CellElementSubRegion & elementSubRegion = elementManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr );

            arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementsToNodes = elementSubRegion.nodeList();
            arrayView2d< real64 const > const & elementCenter = elementSubRegion.getElementCenter().toViewConst();
            real64 K = bulkModulus[er][esr][m_solidMaterialFullIndex][ei];
            real64 G = shearModulus[er][esr][m_solidMaterialFullIndex][ei];
            real64 YoungModulus = 9 * K * G / ( 3 * K + G );
            real64 poissonRatio = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );

            localIndex const numQuadraturePoints = detJ[er][esr].size( 1 );

            for( localIndex n=0; n<elementsToNodes.size( 1 ); ++n )
            {
              if( elementsToNodes( ei, n ) == nodeIndex )
              {
                real64 temp[ 3 ] = {0};
                real64 xEle[ 3 ]  = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( elementCenter[ei] );

                SolidMechanicsLagrangianFEMKernels::ExplicitKernel::
                  calculateSingleNodalForce( ei,
                                             n,
                                             numQuadraturePoints,
                                             dNdX[er][esr],
                                             detJ[er][esr],
                                             stress[er][esr][m_solidMaterialFullIndex],
                                             temp );

                //wu40: the nodal force need to be weighted by Young's modulus and possion's ratio.
                LvArray::tensorOps::scale< 3 >( temp, YoungModulus );
                LvArray::tensorOps::scale< 3 >( temp, 1.0 / (1 - poissonRatio * poissonRatio) );

                LvArray::tensorOps::subtract< 3 >( xEle, nodePosition );
                if( LvArray::tensorOps::AiBi< 3 >( xEle, faceNormalVector ) > 0 ) //TODO: check the sign.
                {
                  nElemEachSide[0] += 1;
                  LvArray::tensorOps::add< 3 >( nodeDisconnectForce, temp );
                }
                else
                {
                  nElemEachSide[1] +=1;
                  LvArray::tensorOps::subtract< 3 >( nodeDisconnectForce, temp );
                }
              }
            }
          }

          if( nElemEachSide[0]>=1 && nElemEachSide[1]>=1 )
          {
            LvArray::tensorOps::scale< 3 >( nodeDisconnectForce, 0.5 );
          }

          //Find the trailing node according to the node index and face index
          if( unpinchedNodeID.size() == 0 ) //Tet mesh under three nodes pinched scenario. Need to find the other
                                            // trailing face that containing the trailing node.
          {
            for( localIndex const edgeIndex: faceToEdgeMap[ trailingFaceIndex ] )
            {
              for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
              {
                if( faceIndex != trailingFaceIndex && m_tipFaces.contains( faceIndex ))
                {
                  for( localIndex const iNode: faceToNodeMap[ faceIndex ] )
                  {
                    if( !m_tipNodes.contains( iNode ))
                    {
                      tralingNodeID = iNode;
                    }
                  }
                }
              }
            }

            if( tralingNodeID == std::numeric_limits< localIndex >::max())
            {
              GEOSX_ERROR( "Error. The triangular trailing face has three tip nodes but cannot find the other trailing face containing the trailing node." );
            }
          }
          else if( unpinchedNodeID.size() == 1 )
          {
            tralingNodeID = unpinchedNodeID[0];
          }
          else if( unpinchedNodeID.size() == 2 )
          {
            for( localIndex const edgeIndex : nodeToEdgeMap[ nodeIndex ] )
            {
              auto const faceToEdgeMapIterator = faceToEdgeMap[ trailingFaceIndex ];
              if( std::find( faceToEdgeMapIterator.begin(), faceToEdgeMapIterator.end(), edgeIndex ) != faceToEdgeMapIterator.end() &&
                  !m_tipEdges.contains( edgeIndex ) )
              {
                tralingNodeID = edgeToNodeMap[edgeIndex][0] == nodeIndex ? edgeToNodeMap[edgeIndex][1] : edgeToNodeMap[edgeIndex][0];
              }
            }
          }

          //Calculate SIF for the node.
          real64 tipNodeSIF;
          real64 tipNodeForce[3];
          real64 trailingNodeDisp[3];
          localIndex theOtherTrailingNodeID;

          if( childNodeIndices[tralingNodeID] == -1 )
          {
            theOtherTrailingNodeID = parentNodeIndices[tralingNodeID];
          }
          else
          {
            theOtherTrailingNodeID = childNodeIndices[tralingNodeID];
          }

          LvArray::tensorOps::copy< 3 >( trailingNodeDisp, displacement[theOtherTrailingNodeID] );
          LvArray::tensorOps::subtract< 3 >( trailingNodeDisp, displacement[tralingNodeID] );

          //Calculate average young's modulus and poisson ratio for fext.
          real64 fExternal[2][3];
          for( localIndex i=0; i<2; ++i )
          {
            real64 averageYoungModulus( 0 ), averagePoissonRatio( 0 );
            localIndex nodeID = i == 0 ? tralingNodeID : theOtherTrailingNodeID;
            for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
            {
              localIndex const er  = nodeToRegionMap[nodeIndex][k];
              localIndex const esr = nodeToSubRegionMap[nodeIndex][k];
              localIndex const ei  = nodeToElementMap[nodeIndex][k];

              real64 K = bulkModulus[er][esr][m_solidMaterialFullIndex][ei];
              real64 G = shearModulus[er][esr][m_solidMaterialFullIndex][ei];
              averageYoungModulus += 9 * K * G / ( 3 * K + G );
              averagePoissonRatio += ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
            }

            averageYoungModulus /= nodeToRegionMap.sizeOfArray( nodeID );
            averagePoissonRatio /= nodeToRegionMap.sizeOfArray( nodeID );

            LvArray::tensorOps::copy< 3 >( fExternal[i], fext[nodeID] );
            LvArray::tensorOps::scale< 3 >( fExternal[i], averageYoungModulus / (1 - averagePoissonRatio * averagePoissonRatio) );
          }

          //TODO: The sign of fext here is opposite to the sign of fFaceA in function "CalculateEdgeSIF".
          tipNodeForce[0] = nodeDisconnectForce[0] - ( fExternal[0][0] - fExternal[1][0] ) / 2.0;
          tipNodeForce[1] = nodeDisconnectForce[1] - ( fExternal[0][1] - fExternal[1][1] ) / 2.0;
          tipNodeForce[2] = nodeDisconnectForce[2] - ( fExternal[0][2] - fExternal[1][2] ) / 2.0;

//          tipNodeForce[0] = nodeDisconnectForce[0];
//          tipNodeForce[1] = nodeDisconnectForce[1];
//          tipNodeForce[2] = nodeDisconnectForce[2];

          real64 tipArea;
          tipArea = faceArea( trailingFaceIndex );
          if( faceToNodeMap.sizeOfArray( trailingFaceIndex ) == 3 )
          {
            tipArea *= 2.0;
          }

          tipNodeSIF = pow( (fabs( tipNodeForce[0] * trailingNodeDisp[0] / 2.0 / tipArea ) + fabs( tipNodeForce[1] * trailingNodeDisp[1] / 2.0 / tipArea )
                             + fabs( tipNodeForce[2] * trailingNodeDisp[2] / 2.0 / tipArea )), 0.5 );

          SIFNode_All[nodeIndex].emplace_back( tipNodeSIF );


          //Calculate SIF on tip faces connected to this trailing face and the tip node.
          for( localIndex const edgeIndex: tipEdgesID )
          {
            if( edgeToNodeMap[edgeIndex][0] == nodeIndex || edgeToNodeMap[edgeIndex][1] == nodeIndex )
            {
              real64 SIF_I = 0, SIF_II = 0, /*SIF_III,*/ SIF_Face;
              real64 vecTipNorm[3], vecTip[3], tipForce[3], tipOpening[3];

              LvArray::tensorOps::copy< 3 >( vecTipNorm, faceNormal[trailingFaceIndex] );
              LvArray::tensorOps::subtract< 3 >( vecTipNorm, faceNormal[childFaceIndices[trailingFaceIndex]] );
              LvArray::tensorOps::normalize< 3 >( vecTipNorm );

              real64 vecEdge[3];
              edgeManager.calculateLength( edgeIndex, X, vecEdge );
              LvArray::tensorOps::normalize< 3 >( vecEdge );

              LvArray::tensorOps::crossProduct( vecTip, vecTipNorm, vecEdge );
              LvArray::tensorOps::normalize< 3 >( vecTip );
              real64 v0[3];
              edgeManager.calculateCenter( edgeIndex, X, v0 );
              LvArray::tensorOps::subtract< 3 >( v0, faceCenter[ trailingFaceIndex ] );

              if( LvArray::tensorOps::AiBi< 3 >( v0, vecTip ) < 0 )
                LvArray::tensorOps::scale< 3 >( vecTip, -1.0 );

              tipForce[0] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTipNorm ) -
                            ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecTipNorm ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecTipNorm ) ) / 2.0;
              tipForce[1] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTip ) -
                            ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecTip ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecTip ) ) / 2.0;
              tipForce[2] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecEdge ) -
                            ( LvArray::tensorOps::AiBi< 3 >( fExternal[0], vecEdge ) - LvArray::tensorOps::AiBi< 3 >( fExternal[1], vecEdge ) ) / 2.0;

//              tipForce[0] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTipNorm );
//              tipForce[1] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecTip );
//              tipForce[2] = LvArray::tensorOps::AiBi< 3 >( nodeDisconnectForce, vecEdge );

              tipOpening[0] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecTipNorm );
              tipOpening[1] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecTip );
              tipOpening[2] = LvArray::tensorOps::AiBi< 3 >( trailingNodeDisp, vecEdge );

//              if( tipForce[0] > 0.0 )
              {
                SIF_I = pow( fabs( tipForce[0] * tipOpening[0] / 2.0 / tipArea ), 0.5 );
                SIF_II = pow( fabs( tipForce[1] * tipOpening[1] / 2.0 / tipArea ), 0.5 );
//              SIF_III = pow( fabs( tipForce[2] * tipOpening[2] / 2.0 / tipArea ), 0.5 );
              }

              if( tipOpening[0] < 0 )
              {
                SIF_I *= -1.0;
              }

              if( tipForce[1] < 0.0 )
              {
                SIF_II *= -1.0;
              }

              for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
              {
                if( m_tipFaces.contains( faceIndex ))
                {
                  real64 vecFace[ 3 ];
                  real64 fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( faceCenter[faceIndex] );

                  //Get the vector in the face and normal to the edge.
                  real64 udist;

                  real64 x0_x1[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeIndex][0]] );
                  real64 x0_fc[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );

                  LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeIndex][1]] );
                  LvArray::tensorOps::normalize< 3 >( x0_x1 );
                  LvArray::tensorOps::subtract< 3 >( x0_fc, X[edgeToNodeMap[edgeIndex][1]] );
                  udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_fc );

                  real64 ptPrj[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( x0_x1 );
                  LvArray::tensorOps::scale< 3 >( ptPrj, udist );
                  LvArray::tensorOps::add< 3 >( ptPrj, X[edgeToNodeMap[edgeIndex][1]] );
                  LvArray::tensorOps::copy< 3 >( vecFace, fc );
                  LvArray::tensorOps::subtract< 3 >( vecFace, ptPrj );
                  LvArray::tensorOps::normalize< 3 >( vecFace );

//                  if( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace ) > cos( m_maxTurnAngle ))
                  {
                    // We multiply this by 0.9999999 to avoid an exception caused by acos a number slightly larger than
                    // 1.
                    real64 thetaFace = acos( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace )*0.999999 );

                    real64 tipCrossFace[ 3 ];
                    LvArray::tensorOps::crossProduct( tipCrossFace, vecTip, vecEdge );

                    if( LvArray::tensorOps::AiBi< 3 >( tipCrossFace, vecEdge ) < 0.0 )
                    {
                      thetaFace *= -1.0;
                    }

                    SIF_Face = cos( thetaFace / 2.0 ) *
                               ( SIF_I * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II * sin( thetaFace ) );

                    SIFonFace_All[faceIndex].emplace_back( SIF_Face );
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //wu40: the tip node may be included in multiple trailing faces and SIF of the node/face will be calculated multiple
  // times. We chose the smaller node SIF and the larger face SIF.
  for( localIndex const nodeIndex : m_tipNodes )
  {
    if( isNodeGhost[nodeIndex] < 0 )
    {
      SIFNode[nodeIndex] = *min_element( SIFNode_All[nodeIndex].begin(), SIFNode_All[nodeIndex].end());

      for( localIndex const edgeIndex: m_tipEdges )
      {
        if( edgeToNodeMap[edgeIndex][0] == nodeIndex || edgeToNodeMap[edgeIndex][1] == nodeIndex )
        {
          for( localIndex const faceIndex: edgeToFaceMap[ edgeIndex ] )
          {
            if( m_tipFaces.contains( faceIndex ))
            {
              SIFonFace[faceIndex] = *max_element( SIFonFace_All[faceIndex].begin(), SIFonFace_All[faceIndex].end());
            }
          }
        }
      }
    }
  }
}

real64 SurfaceGenerator::calculateEdgeSif( DomainPartition & domain,
                                           const localIndex edgeID,
                                           localIndex & trailFaceID,
                                           NodeManager & nodeManager,
                                           EdgeManager & edgeManager,
                                           FaceManager & faceManager,
                                           ElementRegionManager & elementManager,
                                           real64 ( & vecTipNorm )[ 3 ],
                                           real64 ( & vecTip )[ 3 ] )
{
  real64 rval;
  localIndex_array faceInvolved;
  arrayView1d< real64 > const & SIF_I = edgeManager.getExtrinsicData< extrinsicMeshData::SIF_I >();
  arrayView1d< real64 > const & SIF_II = edgeManager.getExtrinsicData< extrinsicMeshData::SIF_II >();
  arrayView1d< real64 > const & SIF_III = edgeManager.getExtrinsicData< extrinsicMeshData::SIF_III >();

  ArrayOfSetsView< localIndex const > const & nodeToEdgeMap = nodeManager.edgeList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & displacement = nodeManager.totalDisplacement();

  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  arrayView1d< localIndex > const & faceParentIndex = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  SIF_I[edgeID] = 0.0;
  SIF_II[edgeID] = 0.0;
  SIF_III[edgeID] = 0.0;

  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceIsExternal[iface] >= 1 )
    {
      faceInvolved.emplace_back( iface );
    }
  }

  // Figure out the two fracture faces connected to this edge
  localIndex faceA( 0 ), faceAp( 0 );
  if( (faceParentIndex[faceInvolved[0]] == -1 && faceParentIndex[faceInvolved[1]] == faceInvolved[0]) ||
      (faceParentIndex[faceInvolved[1]] == -1 && faceParentIndex[faceInvolved[0]] == faceInvolved[1]) )
  {
    faceA = faceInvolved[0];
    faceAp = faceInvolved[1];
  }
  else
  {
    char msg[200];
    sprintf( msg, "Error! Edge %d has two external faces, but the parent-child relationship is wrong.", int(edgeID));
    GEOSX_ERROR( msg );
  }

  trailFaceID = faceParentIndex[faceInvolved[0]]==-1 ? faceInvolved[0] : faceParentIndex[faceInvolved[0]];


  // We define three unit vectors
  // vecEdge: pointing from node 0 to node 1 along the tip edge
  // vecTip: pointing from the opening into the solid
  // vecTipNorm: normal of the one of the fracture faces;  vecTip X vecTipNorm should point to the direction of vecEdge

  LvArray::tensorOps::copy< 3 >( vecTipNorm, faceNormal[faceA] );
  LvArray::tensorOps::subtract< 3 >( vecTipNorm, faceNormal[faceAp] );
  LvArray::tensorOps::normalize< 3 >( vecTipNorm );

  //TODO: wu40: There is a function for EdgeVector in EdgeManager.cpp but has been commented.
  real64 vecEdge[3];
  edgeManager.calculateLength( edgeID, X, vecEdge );
  real64 const edgeLength = LvArray::tensorOps::l2Norm< 3 >( vecEdge );

  LvArray::tensorOps::crossProduct( vecTip, vecTipNorm, vecEdge );
  LvArray::tensorOps::normalize< 3 >( vecTip );
  real64 v0[3];
  edgeManager.calculateCenter( edgeID, X, v0 );
  LvArray::tensorOps::subtract< 3 >( v0, faceCenter[faceA] );

  if( LvArray::tensorOps::AiBi< 3 >( v0, vecTip ) < 0 )
    LvArray::tensorOps::scale< 3 >( vecTip, -1.0 );

  real64 tipCrossTipNorm[ 3 ];
  LvArray::tensorOps::crossProduct( tipCrossTipNorm, vecTip, vecTipNorm );
  if( LvArray::tensorOps::AiBi< 3 >( tipCrossTipNorm, vecEdge ) < 0 )
  {
    LvArray::tensorOps::scale< 3 >( vecTipNorm, -1 );
    faceA = faceInvolved[1];
    faceAp = faceInvolved[0];
  }


  //Now we need to figure out if a special situation applies to this edge
  // where the fracture face is a quad and three of the nodes are still pinched
  // We use a different algorithm for this special situation.

  bool threeNodesPinched( false );
  localIndex_array openNodeID;

  if( faceToNodeMap.sizeOfArray( faceA ) == 4 )  // Only quads have this problem
  {
    int numSharedNodes = 2;
    localIndex_array lNodeFaceA, lNodeFaceAp;

    lNodeFaceA.insert( 0, faceToNodeMap[ faceA ].begin(), faceToNodeMap[ faceA ].end() );
    lNodeFaceAp.insert( 0, faceToNodeMap[ faceAp ].begin(), faceToNodeMap[ faceAp ].end() );

    //We remove all the shared nodes and the one remains should be the open one.
    lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgeToNodeMap[edgeID][0] ))));
    lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), edgeToNodeMap[edgeID][1] ))));
    lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgeToNodeMap[edgeID][0] ))));
    lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), edgeToNodeMap[edgeID][1] ))));

    for( localIndex const j : faceToNodeMap[ faceA ] )
    {
      localIndex iNd = j;
      if( iNd != edgeToNodeMap[edgeID][0] && iNd != edgeToNodeMap[edgeID][1] )
      {
        auto faceToNodeMapIterator = faceToNodeMap[ faceAp ];
        if( std::find( faceToNodeMapIterator.begin(), faceToNodeMapIterator.end(), iNd ) != faceToNodeMapIterator.end())
        {
          numSharedNodes++;
          lNodeFaceA.erase( std::distance( lNodeFaceA.begin(), (std::find( lNodeFaceA.begin(), lNodeFaceA.end(), iNd ))));
          lNodeFaceAp.erase( std::distance( lNodeFaceAp.begin(), (std::find( lNodeFaceAp.begin(), lNodeFaceAp.end(), iNd ))));
        }
      }
    }

    if( numSharedNodes == 4 )
    {
      GEOSX_ERROR( "Error.  The fracture face has four shared nodes with its child.  This should not happen." );
    }
    else if( numSharedNodes == 3 )
    {
      threeNodesPinched = true;

      //wu40: I think the following check is not necessary.
      if( lNodeFaceA.size() != 1 || lNodeFaceAp.size() != 1 )
      {
        GEOSX_ERROR( "Error. These two faces share three nodes but the number of remaining nodes is not one.  Something is wrong" );
      }
      else
      {
        openNodeID.emplace_back( lNodeFaceA[0] );
        openNodeID.emplace_back( lNodeFaceAp[0] );
      }
    }
  }


  // Now we need to identify which node on the edge is the convex point and which one is the concave corner.  The convex
  // node must share an edge with the open node.
  localIndex convexCorner( std::numeric_limits< localIndex >::max());
  if( threeNodesPinched )
  {
    localIndex iNd, jNd;
    iNd = edgeToNodeMap[edgeID][0];
    jNd = edgeToNodeMap[edgeID][1];
    for( localIndex const j : faceToEdgeMap[ faceA ] )
    {
      localIndex edge = j;
      if((openNodeID[0] == edgeToNodeMap[edge][0] && iNd == edgeToNodeMap[edge][1]) ||
         (openNodeID[0] == edgeToNodeMap[edge][1] && iNd == edgeToNodeMap[edge][0])
         )
      {
        convexCorner = iNd;
        break;
      }
      if((openNodeID[0] == edgeToNodeMap[edge][0] && jNd == edgeToNodeMap[edge][1]) ||
         (openNodeID[0] == edgeToNodeMap[edge][1] && jNd == edgeToNodeMap[edge][0])
         )
      {
        convexCorner = jNd;
        break;
      }
    }

    if( convexCorner == std::numeric_limits< localIndex >::max())
      GEOSX_ERROR( "Error.  This is a three-node-pinched edge but I cannot find the convex corner" );

  }


  // Calculate element forces acting on this edge, i.e., f_disconnect.  Need to add nodal forces from two nodes up.
  //An element has to be within the range of this edge to be included.
  //For the threeNodesPinched case, we only use the force on the node at the convex point, not the concave point.  The
  // force at the former is usually greater, so we just pick the great one instead of doing a geometrical check.
  real64 fNodeO[3] = { 0.0, 0.0, 0.0 };
  real64 GdivBeta = 0.0;  // Need this for opening-based SIF

  localIndex_array nodeIndices;

  if( !threeNodesPinched )
  {
    for( localIndex a=0; a<edgeToNodeMap.size( 1 ); ++a )
    {
      nodeIndices.emplace_back( edgeToNodeMap( edgeID, a ) );
    }
  }
  else
  {
    nodeIndices.emplace_back( convexCorner );
  }

  calculateElementForcesOnEdge ( domain, edgeID, edgeLength, nodeIndices,
                                 nodeManager, edgeManager, elementManager, vecTipNorm, fNodeO, GdivBeta, threeNodesPinched, false );


  localIndex tipFaces[2];
  tipFaces[0] = faceA;
  tipFaces[1] = faceAp;

  // Now calculate f_u. We have to subtract the nodal force at other nodes (trailing nodes) on these two open faces to
  // take into account
  // the effects of surface traction along the fracture.
  // Finding the two trailing nodes on a hex mesh is pretty straightforward, while it is cumbersome to do in tet mesh
  // For the threeNodesPinched case, this should be the open node.
  real64 fFaceA[2][3];

  // If the two external faces connected to a trailing edge are not coplanar, then we have the risk of incomplete
  // topology.
  // In that case, we use a displacement/opening based method, not VCCT.
  bool incompleteTrailingEdgeTopology = false;

  for( localIndex i=0; i<2; ++i )
  {
    localIndex_array trailingNodes;
    trailingNodes.clear();
    if( threeNodesPinched )
    {
      trailingNodes.emplace_back( openNodeID[i] );
    }
    else
    {
      localIndex faceID = tipFaces[i];
      LvArray::tensorOps::fill< 3 >( fFaceA[i], 0.0 );

      for( localIndex const j : faceToNodeMap[ faceID ] )
      {
        if( j != edgeToNodeMap( edgeID, 0 ) && j != edgeToNodeMap( edgeID, 1 ) ) // This is not a node along the tip
                                                                                 // edge
        {
          trailingNodes.emplace_back( j );
        }
      }

      if( trailingNodes.size() > 2 || trailingNodes.size() == 0 )
      {
        GEOSX_ERROR( "Fatal error in finding nodes behind tip edge." );
      }
      else if( trailingNodes.size() == 1 )  // Need some work to find the other node
      {
        // First find an edge that is connected to this node and parallel to the tip edge
        real64 maxCosAngle = 0.0;
        localIndex pickedTrailingEdge = std::numeric_limits< localIndex >::max();
        for( localIndex const iedge : nodeToEdgeMap[ trailingNodes[ 0 ] ] )
        {
          real64 xTrailingEdge[3];
          edgeManager.calculateCenter( iedge, X, xTrailingEdge );

          real64 udist;
          real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
          real64 x0_xTrailingEdge[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xTrailingEdge );
          LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap( edgeID, 1 )] );
          LvArray::tensorOps::normalize< 3 >( x0_x1 );
          LvArray::tensorOps::subtract< 3 >( x0_xTrailingEdge, X[edgeToNodeMap( edgeID, 1 )] );
          udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_xTrailingEdge );

          if( udist <= edgeLength && udist > 0.0 )
          {
            real64 vEdge[3];
            edgeManager.calculateLength( iedge, X, vEdge );
            LvArray::tensorOps::normalize< 3 >( vEdge );

            real64 cosEdge = std::fabs( LvArray::tensorOps::AiBi< 3 >( vEdge, vecEdge ));
            if( cosEdge > maxCosAngle )
            {
              maxCosAngle = cosEdge;
              pickedTrailingEdge = iedge;
            }
          }
        }
        if( maxCosAngle > 0.75 )
          trailingNodes.emplace_back( edgeToNodeMap[pickedTrailingEdge][0] + edgeToNodeMap[pickedTrailingEdge][1] - trailingNodes[0] );
      }
    }

    localIndex trailingEdge;
    trailingEdge = std::numeric_limits< localIndex >::max();

    if( trailingNodes.size() == 2 )
    {
      //wu40: TODO: This check is from GEOS. I think this may not be necessary. Check with Randy and PC.
      if( trailingNodes[0] != trailingNodes[1] )
      {
        for( localIndex const iedge : nodeToEdgeMap[ trailingNodes[ 0 ] ] )
        {
          if( edgeToNodeMap[iedge][0] == trailingNodes[1] || edgeToNodeMap[iedge][1] == trailingNodes[1] )
          {
            trailingEdge = iedge;
          }
        }
      }

      if( trailingEdge > edgeManager.size())
      {
        int const rank = MpiWrapper::commRank( MPI_COMM_WORLD );
        std::cout << "Cannot find trailing edge (edge=" << edgeID << ", rank=" << rank <<   "  )" << std::endl;
        return 0.0;
      }

      localIndex_array extFacesOnTrailingEdge;
      for( localIndex const iface : edgeToFaceMap[ trailingEdge ] )
      {
        if( faceIsExternal[iface] >= 1 )
          extFacesOnTrailingEdge.emplace_back( iface );
      }

      if( extFacesOnTrailingEdge.size() != 2 )
      {
        incompleteTrailingEdgeTopology = true;
      }
      else
      {
        real64 extFaceNormal[2][3];
        for( localIndex j = 0; j < 2; ++j )
        {
          LvArray::tensorOps::copy< 3 >( extFaceNormal[j], faceNormal[extFacesOnTrailingEdge[j]] );
        }

        if( std::fabs( LvArray::tensorOps::AiBi< 3 >( extFaceNormal[0], extFaceNormal[1] )) < 0.9 ) //The two faces are not coplanar.
        {
          incompleteTrailingEdgeTopology = true;
        }
      }
    }

    calculateElementForcesOnEdge ( domain, edgeID, edgeLength, trailingNodes,
                                   nodeManager, edgeManager, elementManager, vecTipNorm, fFaceA[i], GdivBeta, threeNodesPinched, true );

  }


  real64 tipForce[3];
  tipForce[0] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecTipNorm ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecTipNorm ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecTipNorm ) / 2.0;
  tipForce[1] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecTip ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecTip ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecTip ) / 2.0;
  tipForce[2] = LvArray::tensorOps::AiBi< 3 >( fNodeO, vecEdge ) + LvArray::tensorOps::AiBi< 3 >( fFaceA[0], vecEdge ) / 2.0 - LvArray::tensorOps::AiBi< 3 >( fFaceA[1], vecEdge ) / 2.0;

  real64 tipDisplacement[3], tipOpening[3], tipFaceDisplacement[2][3];

  if( !threeNodesPinched )
  {
    for( localIndex i=0; i<2; ++i )
    {
      localIndex faceID = tipFaces[i];
      LvArray::tensorOps::fill< 3 >( tipFaceDisplacement[i], 0.0 );

      for( localIndex const j : faceToNodeMap[ faceID ] )
      {
        if( j != edgeToNodeMap( edgeID, 0 ) && j != edgeToNodeMap( edgeID, 1 ))
        {
          LvArray::tensorOps::add< 3 >( tipFaceDisplacement[i], displacement[j] );
        }
      }

      LvArray::tensorOps::scale< 3 >( tipFaceDisplacement[i], 1.0 / (faceToNodeMap.sizeOfArray( faceID ) - 2) );
    }
    LvArray::tensorOps::copy< 3 >( tipDisplacement, tipFaceDisplacement[1] );
    LvArray::tensorOps::subtract< 3 >( tipDisplacement, tipFaceDisplacement[0] );
  }
  else
  {
    LvArray::tensorOps::copy< 3 >( tipDisplacement, displacement[openNodeID[1]] );
    LvArray::tensorOps::subtract< 3 >( tipDisplacement, displacement[openNodeID[0]] );
  }

  tipOpening[0] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecTipNorm );
  tipOpening[1] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecTip );
  tipOpening[2] = LvArray::tensorOps::AiBi< 3 >( tipDisplacement, vecEdge );

  real64 tipArea;
  tipArea = faceArea( faceA );
  if( faceToNodeMap.sizeOfArray( faceA ) == 3 )
  {
    tipArea *= 2.0;
  }

  if( !incompleteTrailingEdgeTopology && tipOpening[0] * tipForce[0] > 0.0 )
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
    real64 r = tipArea / edgeLength;
    SIF_I[edgeID] = tipOpening[0] / 2.0 * GdivBeta / pow( r/6.28, 0.5 );
    SIF_II[edgeID] = 0.0;  // SIF is not accurate in this scenario anyway.  Let's not worry about turning.
    SIF_III[edgeID] = 0.0;
  }

  if( tipForce[1] < 0.0 )
  {
    SIF_II[edgeID] *= -1.0;
  }

  if( SIF_I[edgeID] > 0.0 )
  {
    rval = pow( SIF_I[edgeID]*SIF_I[edgeID]+SIF_II[edgeID]*SIF_II[edgeID]+SIF_III[edgeID]*SIF_III[edgeID], 0.5 );
  }
  else
  {
    rval = -1.0;
  }

//  std::cout << "EdgeID: " << edgeID << " SIF: " << rval << std::endl;

  return rval;
}


int SurfaceGenerator::calculateElementForcesOnEdge( DomainPartition & domain,
                                                    const localIndex edgeID,
                                                    real64 edgeLength,
                                                    localIndex_array & nodeIndices,
                                                    NodeManager & nodeManager,
                                                    EdgeManager & edgeManager,
                                                    ElementRegionManager & elementManager,
                                                    real64 ( & vecTipNorm )[3],
                                                    real64 ( & fNode )[3],
                                                    real64 & GdivBeta,
                                                    bool threeNodesPinched,
                                                    bool calculatef_u )
{
  ArrayOfArraysView< localIndex const > const & nodeToRegionMap = nodeManager.elementRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToSubRegionMap = nodeManager.elementSubRegionList().toViewConst();
  ArrayOfArraysView< localIndex const > const & nodeToElementMap = nodeManager.elementList().toViewConst();

  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  ConstitutiveManager const & cm = domain.getConstitutiveManager();
  ConstitutiveBase const & solid = cm.getConstitutiveRelation< ConstitutiveBase >( m_solidMaterialNames[0] );
  m_solidMaterialFullIndex = solid.getIndexInParent();

  ConstitutiveManager & constitutiveManager =
    domain.getGroup< ConstitutiveManager >( keys::ConstitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "shearModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elementManager.constructFullMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( "bulkModulus", constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< arrayView3d< real64 const, solid::STRESS_USD > > const
  stress = elementManager.constructFullMaterialViewAccessor< array3d< real64, solid::STRESS_PERMUTATION >,
                                                             arrayView3d< real64 const, solid::STRESS_USD > >( SolidBase::viewKeyStruct::stressString(),
                                                                                                               constitutiveManager );

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > const
  dNdX = elementManager.constructViewAccessor< array4d< real64 >, arrayView4d< real64 const > >( keys::dNdX );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const
  detJ = elementManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( keys::detJ );

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > const elemCenter =
    elementManager.constructViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementCenterString() );

  localIndex nElemEachSide[2];
  nElemEachSide[0] = 0;
  nElemEachSide[1] = 0;

  real64 xEdge[3] = { 0.0 };

  if( !calculatef_u )
  {
    edgeManager.calculateCenter( edgeID, X, xEdge );
  }

  for( localIndex i=0; i < nodeIndices.size(); ++i )
  {
    localIndex nodeID = nodeIndices( i );
//    localIndex_array temp11;
//    for (int ii = 0; ii < nodeToElementMap.sizeOfArray(nodeID); ii++)
//    {
//      temp11.emplace_back(nodeToElementMap[nodeID][ii]);
//    }

    for( localIndex k=0; k<nodeToRegionMap.sizeOfArray( nodeID ); ++k )
    {
      localIndex const er  = nodeToRegionMap[nodeID][k];
      localIndex const esr = nodeToSubRegionMap[nodeID][k];
      localIndex const ei  = nodeToElementMap[nodeID][k];

      CellElementSubRegion const & elementSubRegion = elementManager.getRegion( er ).getSubRegion< CellElementSubRegion >( esr );

      real64 xEle[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( elemCenter[er][esr][ei] );

      real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
      LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeID][1]] );
      LvArray::tensorOps::normalize< 3 >( x0_x1 );

      real64 x0_xEle[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( xEle );
      LvArray::tensorOps::subtract< 3 >( x0_xEle, X[edgeToNodeMap[edgeID][1]] );
      real64 const udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_xEle );

      localIndex const numQuadraturePoints = detJ[er][esr].size( 1 );

      if(( udist <= edgeLength && udist > 0.0 ) || threeNodesPinched )
      {
        real64 K = bulkModulus[er][esr][m_solidMaterialFullIndex][ei];
        real64 G = shearModulus[er][esr][m_solidMaterialFullIndex][ei];
        real64 YoungModulus = 9 * K * G / ( 3 * K + G );
        real64 poissonRatio = ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );

        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elementsToNodes = elementSubRegion.nodeList();
        for( localIndex n=0; n<elementsToNodes.size( 1 ); ++n )
        {
          if( elementsToNodes( ei, n ) == nodeID )
          {
            real64 temp[3]{};
            LvArray::tensorOps::copy< 3 >( xEle, elemCenter[er][esr][ei] ); //For C3D6 element type, elementsToNodes map may include
            // repeated indices and the following may run multiple
            // times for the same element.

            //wu40: the nodal force need to be weighted by Young's modulus and possion's ratio.
            SolidMechanicsLagrangianFEMKernels::ExplicitKernel::
              calculateSingleNodalForce( ei,
                                         n,
                                         numQuadraturePoints,
                                         dNdX[er][esr],
                                         detJ[er][esr],
                                         stress[er][esr][m_solidMaterialFullIndex],
                                         temp );

            LvArray::tensorOps::scale< 3 >( temp, YoungModulus );
            LvArray::tensorOps::scale< 3 >( temp, 1.0 / (1 - poissonRatio * poissonRatio) );

            if( !calculatef_u )
            {
              LvArray::tensorOps::subtract< 3 >( xEle, xEdge );
              if( LvArray::tensorOps::AiBi< 3 >( xEle, vecTipNorm ) > 0 )
              {
                nElemEachSide[0] += 1;
                LvArray::tensorOps::add< 3 >( fNode, temp );

                //wu40: for debug purpose
//                std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//                std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//                std::cout << "Add to total nodal force (fdisc): " << fNode[0] << ", " << fNode[1] << ", " << fNode[2]
// << std::endl;
              }
              else
              {
                nElemEachSide[1] +=1;
                LvArray::tensorOps::subtract< 3 >( fNode, temp );

                //wu40: for debug purpose
//                std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//                std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//                std::cout << "Minus from total nodal force (fdisc): " << fNode[0] << ", " << fNode[1] << ", " <<
// fNode[2] << std::endl;
              }
            }
            else
            {
              LvArray::tensorOps::add< 3 >( fNode, temp );

              //wu40: for debug purpose
//              std::cout << "ElementID: " << iEle << ", NodeID: " << nodeID << std::endl;
//              std::cout << "Nodal force: " << temp[0] << ", " << temp[1] << ", " << temp[2] << std::endl;
//              std::cout << "Add to total nodal force (fext): " << fNode[0] << ", " << fNode[1] << ", " << fNode[2] <<
// std::endl;
            }
          }
        }

        if( !calculatef_u )
        {
          GdivBeta += G /2/(1-poissonRatio);
        }
      }

    }

    //If we only find one node behind the tip for the non-threeNodesPinched scenario, we do the following as a rough
    // compensation for f_u.
    if( calculatef_u && nodeIndices.size() == 1 && !threeNodesPinched )
    {
      LvArray::tensorOps::scale< 3 >( fNode, 2.0 );
    }
  }

  if( !calculatef_u )
  {
    if( nElemEachSide[0]>=1 && nElemEachSide[1]>=1 )
      LvArray::tensorOps::scale< 3 >( fNode, 0.5 );
    //We have contributions from both sides. The two sizes are the two sides of the fracture plane.  If the fracture
    // face
    // is on domain boundary, it's possible to have just one side.
    if( nElemEachSide[0] + nElemEachSide[1] >= 1 )
      GdivBeta /= (nElemEachSide[0] + nElemEachSide[1]);
  }

  return 0;
}

int SurfaceGenerator::checkOrphanElement ( ElementRegionManager & elementManager,
                                           FaceManager & faceManager,
                                           localIndex iFace )
{
  arrayView2d< localIndex const > const & faceToRegionMap = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & faceToSubRegionMap = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();

  arrayView1d< integer const > const & ruptureState = faceManager.getExtrinsicData< extrinsicMeshData::RuptureState >();

  int flagOrphan = 0;
  for( localIndex k=0; k<faceToRegionMap.size( 1 ); ++k )
  {
    localIndex const er = faceToRegionMap[iFace][k];
    localIndex const esr = faceToSubRegionMap[iFace][k];
    localIndex const ei = faceToElementMap[iFace][k];
    if( er != -1 &&  esr != -1 && ei != -1 )
    {
      CellElementSubRegion & elementSubRegion = elementManager.getRegion( faceToRegionMap[iFace][k] ).
                                                  getSubRegion< CellElementSubRegion >( faceToSubRegionMap[iFace][k] );


      int nRuptureFace = 0;
      arrayView2d< localIndex > & elementsToFaces = elementSubRegion.faceList();
      for( localIndex a=0; a < elementsToFaces.size( 1 ); ++a )
      {
        localIndex jFace = elementsToFaces[ei][a];
        if( (ruptureState[jFace] == 1 || faceIsExternal[jFace] >= 1) && jFace != iFace )
        {
          nRuptureFace +=1;
        }
      }

      if( nRuptureFace == elementsToFaces.size( 1 ) - 1 )
      {
        flagOrphan = 1;
      }
    }
  }
  return flagOrphan;

}

void SurfaceGenerator::markRuptureFaceFromNode ( const localIndex nodeIndex,
                                                 NodeManager & nodeManager,
                                                 EdgeManager & edgeManager,
                                                 FaceManager & faceManager,
                                                 ElementRegionManager & GEOSX_UNUSED_PARAM( elementManager ),
                                                 ModifiedObjectLists & modifiedObjects )
{
  arrayView1d< integer > const & ruptureState = faceManager.getExtrinsicData< extrinsicMeshData::RuptureState >();
  arrayView1d< real64 const > const & SIFonFace = faceManager.getExtrinsicData< extrinsicMeshData::SIFonFace >();
  arrayView2d< real64 const > const & KIC = faceManager.getExtrinsicData< extrinsicMeshData::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  localIndex_array eligibleFaces;
  real64_array faceSIFToToughnessRatio;
  real64 lowestSIF = std::numeric_limits< real64 >::max();
  real64 highestSIF = std::numeric_limits< real64 >::min();

  for( localIndex const faceIndex : nodeToFaceMap[ nodeIndex ] )
  {
    if( m_tipFaces.contains( faceIndex ))
    {
      real64 faceToughness;
      real64 fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[faceIndex] );

      eligibleFaces.emplace_back( faceIndex );

      for( localIndex const edgeIndex : faceToEdgeMap[ faceIndex ] )
      {
        if( m_tipEdges.contains( edgeIndex ))
        {
          real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
          real64 edgeCenter[3];
          edgeManager.calculateCenter( edgeIndex, X, edgeCenter );
          LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
          LvArray::tensorOps::normalize< 3 >( direction );
          faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[faceIndex] ));

          faceSIFToToughnessRatio.emplace_back( SIFonFace[faceIndex]/faceToughness );
          highestSIF = std::max( highestSIF, SIFonFace[faceIndex]/faceToughness );
          lowestSIF = std::min( lowestSIF, SIFonFace[faceIndex]/faceToughness );
        }
      }
    }
  }

  for( localIndex i = 0; i < eligibleFaces.size(); ++i )
  {
    localIndex pickedFace = eligibleFaces[i];

    if( highestSIF > 1.0 &&
        ((eligibleFaces.size() < 3) || (eligibleFaces.size() >= 3 && (highestSIF - faceSIFToToughnessRatio[i]) <= 0.2 * (highestSIF - lowestSIF))))
    {
      ruptureState[pickedFace] = 1;
      modifiedObjects.modifiedFaces.insert( pickedFace );

      // Next we mark the faces that are 1) connected to this face, and 2) attached to one node of the edge (implicitly
      // satisfied), and 3) almost co-plane with this face
//      if( m_markExtendedLayer == 1)
//      {
//        for( auto iedge : faceToEdgeMap[pickedFace] )
//        {
//          for( auto iface : edgeToFaceMap[iedge] )
//          {
//            if( iface != pickedFace && isFaceSeparable[iface] == 1 && faceManager.isExternal()[iface] < 1 &&
//                fabs(LvArray::tensorOps::AiBi< 3 >(faceNormals[pickedFace], faceNormals[iface])) > cos( m_maxTurnAngle ) &&
//                ((faceToNodeMap[iface].size() == 3) || (faceToNodeMap[iface].size() == 4 &&
//                    (std::find(faceToNodeMap[iface].begin(), faceToNodeMap[iface].end(), nodeIndex) !=
// faceToNodeMap[iface].end()))))
//            {
//              //wu40: Under tet mesh scenario, the face next to the pickedFace should also be marked but it may not
// necessarily connect to the tip node.
//              bool ruptureFace = true;
//              for (auto edgeIndex: faceToEdgeMap[pickedFace])
//              {
//                if (m_tipEdges.contains(edgeIndex))
//                {
//                  R1Tensor fc;
//                  real64 uDist, segmentLength;
//
//                  fc = faceCenter[iface];
//
//                  R1Tensor x0_x1(X[edgeToNodeMap[edgeIndex][0]]), x0_fc(fc);
//                  x0_x1 -= X[edgeToNodeMap[edgeIndex][1]];
//                  segmentLength = x0_x1.Normalize();
//                  x0_fc -= X[edgeToNodeMap[edgeIndex][1]];
//                  uDist = LvArray::tensorOps::AiBi< 3 >(x0_x1, x0_fc);
//
//                  if (uDist / segmentLength < -m_faceToEdgeProjectionTol || uDist / segmentLength > 1 +
// m_faceToEdgeProjectionTol)
//                  {
//                    ruptureFace = false;
//                  }
//                }
//              }
//
//              if (ruptureFace)
//              {
//                ruptureState[iface] = 1;
//                modifiedObjects.modifiedFaces.insert( iface );
//              }
//            }
//          }
//        }
//      }
    }
  }
}

void SurfaceGenerator::markRuptureFaceFromEdge ( localIndex const edgeID,
                                                 localIndex & GEOSX_UNUSED_PARAM( trailFaceID ),
                                                 NodeManager & nodeManager,
                                                 EdgeManager & edgeManager,
                                                 FaceManager & faceManager,
                                                 ElementRegionManager & elementManager,
                                                 real64 ( &GEOSX_UNUSED_PARAM( vecTipNorm ) )[ 3 ],
                                                 real64 ( & vecTip )[ 3 ],
                                                 ModifiedObjectLists & modifiedObjects,
                                                 int const edgeMode )
{
  arrayView1d< integer > const & ruptureState = faceManager.getExtrinsicData< extrinsicMeshData::RuptureState >();
  arrayView1d< real64 > const & SIFonFace = faceManager.getExtrinsicData< extrinsicMeshData::SIFonFace >();
  arrayView2d< real64 const > const & KIC = faceManager.getExtrinsicData< extrinsicMeshData::K_IC >();
  arrayView1d< real64 > const & SIF_I = edgeManager.getExtrinsicData< extrinsicMeshData::SIF_I >();
  arrayView1d< real64 > const & SIF_II = edgeManager.getExtrinsicData< extrinsicMeshData::SIF_II >();
  arrayView1d< localIndex > const & primaryCandidateFace = faceManager.getExtrinsicData< extrinsicMeshData::PrimaryCandidateFace >();
  arrayView1d< integer const > const & isFaceSeparable = faceManager.getExtrinsicData< extrinsicMeshData::IsFaceSeparable >();
//  integer_array* dfnIndexMap = faceManager.getReferencePointer<integer_array>( "DFN_Index" );

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();


  arrayView2d< localIndex > const & edgeToNodeMap = edgeManager.nodeList();
  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  arrayView2d< localIndex const > const & faceToElementMap = faceManager.elementList();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  localIndex_array eligibleFaces;
  real64 lowestSIF = std::numeric_limits< real64 >::max();
  real64 highestSIF = std::numeric_limits< real64 >::min();
  real64 lowestScore = std::numeric_limits< real64 >::max();
  real64 highestScore = std::numeric_limits< real64 >::min();
  real64 secondScore = std::numeric_limits< real64 >::min();
  localIndex faceWithHighestScore = std::numeric_limits< localIndex >::max();
  localIndex faceWithSecondScore = std::numeric_limits< localIndex >::max();

  real64 vecEdge[3], edgeCenter[3];
  edgeManager.calculateLength( edgeID, X, vecEdge );
  edgeManager.calculateCenter( edgeID, X, edgeCenter );


  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceToElementMap.size( 1 ) == 2  &&
        faceIsExternal[iface] < 1 &&
        checkOrphanElement( elementManager, faceManager, iface ) == 0 &&
        isFaceSeparable[iface] == 1 )
    {
      real64 const fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[iface] );

      //Get the vector in the face and normal to the edge.
      //wu40: there is a function in GEOS for this calculation. Maybe it's worth to have a function in GEOSX too.
      real64 x0_x1[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[edgeToNodeMap[edgeID][0]] );
      LvArray::tensorOps::subtract< 3 >( x0_x1, X[edgeToNodeMap[edgeID][1]] );
      LvArray::tensorOps::normalize< 3 >( x0_x1 );

      real64 x0_fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
      LvArray::tensorOps::subtract< 3 >( x0_fc, X[edgeToNodeMap[edgeID][1]] );
      real64 const udist = LvArray::tensorOps::AiBi< 3 >( x0_x1, x0_fc );

      real64 ptPrj[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( x0_x1 );
      LvArray::tensorOps::scale< 3 >( ptPrj, udist );
      LvArray::tensorOps::add< 3 >( ptPrj, X[edgeToNodeMap[edgeID][1]] );

      real64 vecFace[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
      LvArray::tensorOps::subtract< 3 >( vecFace, ptPrj );
      LvArray::tensorOps::normalize< 3 >( vecFace );

//      if( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace ) > cos( m_maxTurnAngle ))
      {
        eligibleFaces.emplace_back( iface );
        real64 thetaFace = acos( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace )*0.999999 );  // We multiply this by 0.9999999 to avoid an
        // exception caused by acos a number slightly larger
        // than 1.

        real64 tipCrossFace[ 3 ];
        LvArray::tensorOps::crossProduct( tipCrossFace, vecTip, vecFace );
        if( LvArray::tensorOps::AiBi< 3 >( tipCrossFace, vecEdge ) < 0.0 )
        {
          thetaFace *= -1.0;
        }

        SIFonFace[iface] = cos( thetaFace / 2.0 ) *
                           ( SIF_I[edgeID] * cos( thetaFace / 2.0 ) * cos( thetaFace / 2.0 ) - 1.5 * SIF_II[edgeID] * sin( thetaFace ) );

        real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
        LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
        LvArray::tensorOps::normalize< 3 >( direction );
        real64 faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[iface] ));

        highestSIF = std::max( highestSIF, SIFonFace[iface]/faceToughness );
        lowestSIF = std::min( lowestSIF, SIFonFace[iface]/faceToughness );

      }
    }
  }

  localIndex_array pickedFaces;
  if( eligibleFaces.size() >=1 )
  {
    for( localIndex i = 0; i < eligibleFaces.size(); ++i )
    {
      localIndex iface = eligibleFaces[i];
      real64 direction[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( edgeCenter );
      LvArray::tensorOps::subtract< 3 >( direction, faceCenter[iface] );
      LvArray::tensorOps::normalize< 3 >( direction );
      real64 faceToughness = std::fabs( LvArray::tensorOps::AiBi< 3 >( direction, KIC[iface] ));

      real64 splitabilityScore = SIFonFace[iface] - lowestSIF * faceToughness;
      lowestScore = std::min( lowestScore, splitabilityScore );

      if( faceWithHighestScore == std::numeric_limits< localIndex >::max())
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

    pickedFaces.emplace_back( faceWithHighestScore );

    if( eligibleFaces.size() >= 3 && (highestScore - secondScore) < 0.1 * (highestScore - lowestScore))
    {
      pickedFaces.emplace_back( faceWithSecondScore );
    }

  }


  for( localIndex i = 0; i < pickedFaces.size(); ++i )
  {
    localIndex pickedFace = pickedFaces[i];

    if( highestSIF > 1.0 && edgeMode == 1 && i == 0 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureState[pickedFace] = 1;
//      if( !m_dfnPrefix.empty())
//        (*dfnIndexMap)[pickedFace] = (*dfnIndexMap)[trailFaceID];
      modifiedObjects.modifiedFaces.insert( pickedFace );
    }
    else if( highestSIF > 1.0 && edgeMode == 1 && i == 1 && isFaceSeparable[pickedFace] == 1 )
    {
      ruptureState[pickedFace] = -1;
//      if( !m_dfnPrefix.empty())
//        (*dfnIndexMap)[pickedFace] = (*dfnIndexMap)[trailFaceID];
      modifiedObjects.modifiedFaces.insert( pickedFace );
      primaryCandidateFace[pickedFace] = faceWithHighestScore;
    }


    // We didn't really need to do this unless the criterion above has been satisfied.
    // We are calculating this regardless the criterion for debugging purpose.
//    if( m_markExtendedLayer == 1 && highestSIF > 1.0 && edgeMode == 1 )
//    {
//      // Next we mark the faces that are 1) connected to this face, and 2) attached to one node of the edge
// (implicitly
//      // satisfied), and 3) almost co-plane with this face
//      for( auto iedge : faceToEdgeMap[pickedFace] )
//      {
//        if( iedge != edgeID )
//        {
//          for( auto iface : edgeToFaceMap[iedge] )
//          {
//            if( iface != pickedFace && isFaceSeparable[iface] == 1 && faceManager.isExternal()[iface] < 1 &&
//                ( std::find(faceToNodeMap[iface].begin(), faceToNodeMap[iface].end(), edgeToNodeMap[edgeID][0]) !=
// faceToNodeMap[iface].end() ||
//                  std::find(faceToNodeMap[iface].begin(), faceToNodeMap[iface].end(), edgeToNodeMap[edgeID][1]) !=
// faceToNodeMap[iface].end()))
//            {
//              R1Tensor fc, fn, vecFace, fn0, ptPrj;
//              real64 uDist, segmentLength;
//
//              fc = faceCenter[iface];
//              fn = faceNormal[iface];
//              fn0 = faceNormal[pickedFace];
//
//              R1Tensor x0_x1(X[edgeToNodeMap[edgeID][0]]), x0_fc(fc);
//              x0_x1 -= X[edgeToNodeMap[edgeID][1]];
//              segmentLength = x0_x1.Normalize();
//              x0_fc -= X[edgeToNodeMap[edgeID][1]];
//              uDist = LvArray::tensorOps::AiBi< 3 >(x0_x1, x0_fc);
//
//              ptPrj = x0_x1;
//              ptPrj *= uDist;
//              ptPrj += X[edgeToNodeMap[edgeID][1]];
//              vecFace = fc;
//              vecFace -= ptPrj;
//              vecFace.Normalize();
//
//              // thetaFace does not strictly speaking apply to this face since the tip edge is not a edge of this
// face.
//              // We calculate it as if this face is coplane with the master face
//              real64 thetaFace = acos( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace )*0.999999 );
//              if( LvArray::tensorOps::AiBi< 3 >( Cross( vecTip, vecFace ), vecEdge ) < 0.0 )
//              {
//                thetaFace *= -1.0;
//              }
//
//              if( LvArray::tensorOps::AiBi< 3 >( vecTip, vecFace ) > cos( m_maxTurnAngle ) &&
//                  uDist / segmentLength > -m_faceToEdgeProjectionTol &&
//                  uDist / segmentLength < 1 + m_faceToEdgeProjectionTol &&
//                  fabs( LvArray::tensorOps::AiBi< 3 >( vecEdge, fn )) < m_faceToEdgeCoplaneTol &&  // this face is kind of parallel to the
// tip
//                                                                         // edge.
//                  fabs( LvArray::tensorOps::AiBi< 3 >( fn0, fn )) > 1 - m_faceToFaceCoplaneTol )  // co-plane
//              {
//                if( highestSIF > 1.0 && edgeMode == 1 )
//                {
//                  ruptureState[iface] = ruptureState[pickedFace];
//                  modifiedObjects.modifiedFaces.insert( iface );
//                  primaryCandidateFace[iface] = primaryCandidateFace[pickedFace];
//                }
//              }
//            }
//          }
//        }
//      }
//    }
  }
}

void SurfaceGenerator::postUpdateRuptureStates( NodeManager & nodeManager,
                                                EdgeManager & edgeManager,
                                                FaceManager & faceManager,
                                                ElementRegionManager & GEOSX_UNUSED_PARAM( elementManager ),
                                                std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                                                std::vector< std::set< localIndex > > & edgesToRupturedFaces )
{
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  nodesToRupturedFaces.resize( nodeManager.size() );
  edgesToRupturedFaces.resize( edgeManager.size() );

  arrayView1d< integer const > const & faceRuptureState = faceManager.getExtrinsicData< extrinsicMeshData::RuptureState >();
  arrayView1d< localIndex const > const & faceParentIndex = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  // assign the values of the nodeToRupturedFaces and edgeToRupturedFaces arrays.
  for( localIndex kf=0; kf<faceManager.size(); ++kf )
  {
    if( faceRuptureState[kf] >0 )
    {
      int const n = faceParentIndex[kf]==-1 ? 1 : 2;
      localIndex const faceIndex = faceParentIndex[kf]==-1 ? kf : faceParentIndex[kf];

      for( int i=0; i<n; ++i )
      {
        for( localIndex a=0; a<faceToNodeMap.sizeOfArray( kf ); ++a )
        {
          const localIndex nodeIndex = faceToNodeMap( kf, a );
          nodesToRupturedFaces[nodeIndex].insert( faceIndex );
        }

        for( localIndex a=0; a<faceToEdgeMap.sizeOfArray( kf ); ++a )
        {
          const localIndex edgeIndex = faceToEdgeMap( kf, a );
          edgesToRupturedFaces[edgeIndex].insert( faceIndex );
        }
      }
    }
  }
}

int SurfaceGenerator::checkEdgeSplitability( localIndex const edgeID,
                                             NodeManager & GEOSX_UNUSED_PARAM( nodeManager ),
                                             FaceManager & faceManager,
                                             EdgeManager & edgeManager,
                                             bool const GEOSX_UNUSED_PARAM( prefrac ) )
{
  //     Return value = -1, this edge won't split for sure, don't do any more work;
  //                  = 0, edge is along a tip, but the fracture connected to it is not saturated yet.  We will only
  // calculate SIF but will not perform splitting.
  //                  = 1, edge is along a tip and the adjacent fracture is saturated, more work to be done; or this is
  // a dry simulation
  //                  = 2, this is a singular edge, we need split it.
  //                  = 3, this is an eligible kink, we need to process it as a kink

  ArrayOfSetsView< localIndex const > const & edgeToFaceMap = edgeManager.faceList().toViewConst();
  arrayView1d< localIndex > const & faceParentIndex = faceManager.getExtrinsicData< extrinsicMeshData::ParentIndex >();

  arrayView1d< integer const > const & faceIsExternal = faceManager.isExternal();
  arrayView1d< integer const > const & edgeIsExternal = edgeManager.isExternal();

  int isSplitable = -1;

  if( edgeIsExternal[edgeID] == 0 )
  {
    return isSplitable;
  }

  // We first count the external faces connected to this edge;
  int nExternalFaces = 0;
  localIndex_array faceInvolved;
  for( localIndex const iface : edgeToFaceMap[ edgeID ] )
  {
    if( faceIsExternal[iface] >= 1 )
    {
      nExternalFaces++;
      faceInvolved.emplace_back( iface );
    }
  }

  if( nExternalFaces%2 == 1 )
  {
    //    char msg[200];
    //    sprintf(msg, "Error! Edge %d has an odd number of external faces.", int(edgeID));
    //    GEOSX_ERROR(msg);
    //    std::cout << "Error! Edge " << int(edgeID) << " has an odd number of external faces. "
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][0] << " ,"
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][1] << " ,"
    //        << (*nodeManager.m_refposition)[edgeToNodeMap[edgeID][0]][2] << " ,";
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
    localIndex parentFace = -1;
    if( faceParentIndex[faceInvolved[0]] == -1 && faceParentIndex[faceInvolved[1]] == faceInvolved[0] )
    {
      parentFace = faceInvolved[0];
    }
    else if( faceParentIndex[faceInvolved[1]] == -1 && faceParentIndex[faceInvolved[0]] == faceInvolved[1] )
    {
      parentFace = faceInvolved[1];
    }

    if( parentFace == -1 )
    {
      isSplitable = -1;
    }
    else
    {
      isSplitable = 1;
    }

  }

  return (isSplitable);
}

real64 SurfaceGenerator::minimumToughnessOnEdge( const localIndex edgeID,
                                                 const NodeManager & nodeManager,
                                                 EdgeManager & edgeManager,
                                                 FaceManager & faceManager )
{
  real64 val = std::numeric_limits< real64 >::max();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  real64 edgeCenter[3] = { 0.0 };
  edgeManager.calculateCenter( edgeID, X, edgeCenter );

  arrayView2d< real64 const > const & KIC = faceManager.getExtrinsicData< extrinsicMeshData::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();
  for( localIndex const iface : edgeManager.faceList()[ edgeID ] )
  {
    localIndex const numFaceNodes = faceToNodes.sizeOfArray( iface );
    real64 faceCenter[3] = { 0.0 };
    for( localIndex a=0; a<numFaceNodes; ++a )
    {
      LvArray::tensorOps::add< 3 >( faceCenter, X[ faceToNodes[iface][a] ] );
    }
    LvArray::tensorOps::scale< 3 >( faceCenter, 1.0 / numFaceNodes );

    real64 direction[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter );
    LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
    LvArray::tensorOps::normalize< 3 >( direction );
    val = std::min( val, fabs( LvArray::tensorOps::AiBi< 3 >( KIC[iface], direction ) ) );
  }

  return val;
}

real64 SurfaceGenerator::minimumToughnessOnNode( const localIndex nodeID,
                                                 const NodeManager & nodeManager,
                                                 EdgeManager & edgeManager,
                                                 FaceManager & faceManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  ArrayOfSetsView< localIndex const > const & nodeToFaceMap = nodeManager.faceList().toViewConst();

  arrayView2d< real64 const > const & KIC = faceManager.getExtrinsicData< extrinsicMeshData::K_IC >();
  ArrayOfArraysView< localIndex const > const & faceToEdgeMap = faceManager.edgeList().toViewConst();
  arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

  real64 val = std::numeric_limits< real64 >::max();

  for( localIndex const faceIndex: nodeToFaceMap[ nodeID ] )
  {
    if( m_tipFaces.contains( faceIndex ))
    {
      real64 const fc[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceCenter[faceIndex] );

      for( localIndex const edgeIndex: faceToEdgeMap[ faceIndex ] )
      {
        if( m_tipEdges.contains( edgeIndex ))
        {
          real64 direction[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fc );
          real64 edgeCenter[ 3 ];
          edgeManager.calculateCenter( edgeIndex, X, edgeCenter );
          LvArray::tensorOps::subtract< 3 >( direction, edgeCenter );
          LvArray::tensorOps::normalize< 3 >( direction );

          val = std::min( val, fabs( LvArray::tensorOps::AiBi< 3 >( KIC[faceIndex], direction ) ) );
        }
      }
    }
  }

  return val;
}

void SurfaceGenerator::assignNewGlobalIndicesSerial( ObjectManagerBase & object,
                                                     std::set< localIndex > const & indexList )
{
  // in serial, we can simply loop over the indexList and assign new global indices based on
  // the value of the maxGlobalIndex() + 1;
  arrayView1d< globalIndex > const & localToGlobal = object.localToGlobalMap();
  for( localIndex const newLocalIndex : indexList )
  {
    globalIndex const newGlobalIndex = object.maxGlobalIndex() + 1;

    localToGlobal[newLocalIndex] = newGlobalIndex;
    object.updateGlobalToLocalMap( newLocalIndex );
  }

  object.setMaxGlobalIndex();
}

void SurfaceGenerator::
  assignNewGlobalIndicesSerial( ElementRegionManager & elementManager,
                                map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & newElems )
{
  // in serial, we can simply iterate over the entries in newElems and assign new global indices based on
  // the value of the maxGlobalIndex() + 1 for the ElementRegionManager.

  // loop over entries of newElems, which gives elementRegion/subRegion local indices
  for( auto const & iter : newElems )
  {
    localIndex const er = iter.first.first;
    localIndex const esr = iter.first.second;
    std::set< localIndex > const & indexList = iter.second;

    ElementSubRegionBase & subRegion = elementManager.getRegion( er ).getSubRegion( esr );
    arrayView1d< globalIndex > const & localToGlobal = subRegion.localToGlobalMap();

    // loop over the new elems in the subRegion
    for( localIndex const newLocalIndex : indexList )
    {
      globalIndex const newGlobalIndex = elementManager.maxGlobalIndex() + 1;

      localToGlobal[newLocalIndex] = newGlobalIndex;
      subRegion.updateGlobalToLocalMap( newLocalIndex );
    }
  }

  elementManager.setMaxGlobalIndex();
}

real64
SurfaceGenerator::calculateRuptureRate( SurfaceElementRegion & faceElementRegion,
                                        EdgeManager const & edgeManager )
{
  real64 maxRuptureRate = 0;
  FaceElementSubRegion & subRegion = faceElementRegion.getSubRegion< FaceElementSubRegion >( 0 );

  ArrayOfArraysView< localIndex const > const &
  fractureConnectorEdgesToFaceElements = edgeManager.m_fractureConnectorEdgesToFaceElements.toViewConst();

  arrayView1d< real64 > const &
  ruptureTime = subRegion.getExtrinsicData< extrinsicMeshData::RuptureTime >();

  arrayView1d< real64 > const &
  ruptureRate = subRegion.getExtrinsicData< extrinsicMeshData::RuptureRate >();

  arrayView2d< real64 const > const & elemCenter = subRegion.getElementCenter();

  for( localIndex kfc=0; kfc<fractureConnectorEdgesToFaceElements.size(); ++kfc )
  {
    for( localIndex kfe0=0; kfe0<fractureConnectorEdgesToFaceElements.sizeOfArray( kfc ); ++kfe0 )
    {
      localIndex const faceElem0 = fractureConnectorEdgesToFaceElements( kfc, kfe0 );
      for( localIndex kfe1=kfe0+1; kfe1<fractureConnectorEdgesToFaceElements.sizeOfArray( kfc ); ++kfe1 )
      {
        localIndex const faceElem1 = fractureConnectorEdgesToFaceElements( kfc, kfe1 );
        if( !( m_faceElemsRupturedThisSolve.count( faceElem0 ) && m_faceElemsRupturedThisSolve.count( faceElem1 ) ) )
        {
          real64 const deltaRuptureTime = fabs( ruptureTime( faceElem0 ) - ruptureTime( faceElem1 ));
          if( deltaRuptureTime > 1.0e-14 * (ruptureTime( faceElem0 ) + ruptureTime( faceElem1 )) )
          {
            real64 distance[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3 ( elemCenter[ faceElem0 ] );
            LvArray::tensorOps::subtract< 3 >( distance, elemCenter[ faceElem1 ] );
            real64 const pairwiseRuptureRate = 1.0 / deltaRuptureTime;
            if( m_faceElemsRupturedThisSolve.count( faceElem0 ) )
            {
              ruptureRate( faceElem0 ) = std::min( ruptureRate( faceElem0 ), pairwiseRuptureRate );
            }
            else if( m_faceElemsRupturedThisSolve.count( faceElem1 ) )
            {
              ruptureRate( faceElem1 ) = std::min( ruptureRate( faceElem1 ), pairwiseRuptureRate );
            }
          }
        }
      }
    }
  }


  for( localIndex faceElemIndex : m_faceElemsRupturedThisSolve )
  {
    maxRuptureRate = std::max( maxRuptureRate, ruptureRate( faceElemIndex ) );
  }

  real64 globalMaxRuptureRate;
  MpiWrapper::allReduce( &maxRuptureRate,
                         &globalMaxRuptureRate,
                         1,
                         MPI_MAX,
                         MPI_COMM_GEOSX );

  return globalMaxRuptureRate;
}



REGISTER_CATALOG_ENTRY( SolverBase,
                        SurfaceGenerator,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
