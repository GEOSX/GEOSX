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
 * @file EmbeddedSurfaceGenerator.cpp
 */

#include "EmbeddedSurfaceGenerator.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/ExtrinsicMeshData.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp"
#include "meshUtilities/SimpleGeometricObjects/GeometricObjectManager.hpp"
#include "meshUtilities/SimpleGeometricObjects/BoundedPlane.hpp"

#ifdef USE_GEOSX_PTP
#include "physicsSolvers/GEOSX_PTP/ParallelTopologyChange.hpp"
#endif
#include <set>

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

EmbeddedSurfaceGenerator::EmbeddedSurfaceGenerator( const string & name,
                                                    Group * const parent ):
  SolverBase( name, parent )
{
  registerWrapper( viewKeyStruct::solidMaterialNameString(), &m_solidMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid material used in solid mechanic solver" );

  registerWrapper( viewKeyStruct::fractureRegionNameString(), &m_fractureRegionName ).
    setInputFlag( dataRepository::InputFlags::OPTIONAL ).
    setApplyDefaultValue( "FractureRegion" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

EmbeddedSurfaceGenerator::~EmbeddedSurfaceGenerator()
{}

void EmbeddedSurfaceGenerator::initializePostSubGroups()
{
  /*
   * Here we generate embedded elements for fractures (or faults) that already exist in the domain and
   * were specified in the input file.
   */

  // Get domain
  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();

  // Get geometric object manager
  GeometricObjectManager & geometricObjManager = getGlobalState().getProblemManager().getGroup< GeometricObjectManager >( "Geometry" );

  // Get meshLevel
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // Get managers
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager & nodeManager = meshLevel.getNodeManager();
  EdgeManager & edgeManager = meshLevel.getEdgeManager();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord = nodeManager.referencePosition();

  // Get EmbeddedSurfaceSubRegions
  SurfaceElementRegion & embeddedSurfaceRegion = elemManager.getRegion< SurfaceElementRegion >( this->m_fractureRegionName );
  EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion = embeddedSurfaceRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  localIndex localNumberOfSurfaceElems         = 0;

  // Loop over all the fracture planes
  geometricObjManager.forSubGroups< BoundedPlane >( [&]( BoundedPlane & fracture )
  {
    /* 1. Find out if an element is cut by the fracture or not.
     * Loop over all the elements and for each one of them loop over the nodes and compute the
     * dot product between the distance between the plane center and the node and the normal
     * vector defining the plane. If two scalar products have different signs the plane cuts the
     * cell. If a nodes gives a 0 dot product it has to be neglected or the method won't work.
     */
    real64 const planeCenter[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getCenter() );
    real64 const normalVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( fracture.getNormal() );
    // Initialize variables
    globalIndex nodeIndex;
    integer isPositive, isNegative;
    real64 distVec[ 3 ];

    elemManager.forElementSubRegionsComplete< CellElementSubRegion >(
      [&]( localIndex const er, localIndex const esr, ElementRegionBase &, CellElementSubRegion & subRegion )
    {
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const cellToNodes = subRegion.nodeList();
      FixedOneToManyRelation const & cellToEdges = subRegion.edgeList();

      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      for( localIndex cellIndex = 0; cellIndex < subRegion.size(); cellIndex++ )
      {
        if( ghostRank[cellIndex] < 0 )
        {
          isPositive = 0;
          isNegative = 0;
          for( localIndex kn = 0; kn < subRegion.numNodesPerElement(); kn++ )
          {
            nodeIndex = cellToNodes[cellIndex][kn];
            LvArray::tensorOps::copy< 3 >( distVec, nodesCoord[nodeIndex] );
            LvArray::tensorOps::subtract< 3 >( distVec, planeCenter );
            // check if the dot product is zero
            if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) > 0 )
            {
              isPositive = 1;
            }
            else if( LvArray::tensorOps::AiBi< 3 >( distVec, normalVector ) < 0 )
            {
              isNegative = 1;
            }
          } // end loop over nodes
          if( isPositive * isNegative == 1 )
          {
            bool added = embeddedSurfaceSubRegion.addNewEmbeddedSurface( cellIndex,
                                                                         esr,
                                                                         er,
                                                                         nodeManager,
                                                                         edgeManager,
                                                                         cellToEdges,
                                                                         &fracture );

            if( added )
            {
              GEOSX_LOG_LEVEL_RANK_0( 2, "Element " << cellIndex << " is fractured" );


              // Add the information to the CellElementSubRegion
              subRegion.addFracturedElement( cellIndex, localNumberOfSurfaceElems );
              localNumberOfSurfaceElems++;
            }
          }
        }
      } // end loop over cells
    } );// end loop over subregions
  } );// end loop over thick planes


  // Set the ghostRank form the parent cell
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & cellElemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  embeddedSurfaceSubRegion.inheritGhostRank( cellElemGhostRank );

  setGlobalIndices( elemManager, embeddedSurfaceSubRegion );

  addEmbeddedElementsToSets( elemManager, embeddedSurfaceSubRegion );

  // Populate EdgeManager for embedded surfaces.
  EdgeManager & embSurfEdgeManager = meshLevel.getEmbdSurfEdgeManager();

  EmbeddedSurfaceSubRegion::EdgeMapType & embSurfToEdgeMap = embeddedSurfaceSubRegion.edgeList();
  EmbeddedSurfaceSubRegion::NodeMapType & embSurfToNodeMap = embeddedSurfaceSubRegion.nodeList();

  localIndex numOfPoints = nodeManager.embSurfNodesPosition().size( 0 );

  embSurfEdgeManager.buildEdges( numOfPoints, embSurfToNodeMap.toViewConst(), embSurfToEdgeMap );
}

void EmbeddedSurfaceGenerator::initializePostInitialConditionsPreSubGroups()
{
  // I don't think there is  much to do here.
}

real64 EmbeddedSurfaceGenerator::solverStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                             real64 const & GEOSX_UNUSED_PARAM( dt ),
                                             const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                             DomainPartition & domain )
{
  real64 rval = 0;
  /*
   * This should be the method that generates new fracture elements based on the propagation criterion of choice.
   */
  // Add the embedded elements to the fracture stencil.
  addToFractureStencil( domain );

  return rval;
}

void EmbeddedSurfaceGenerator::addToFractureStencil( DomainPartition & domain )
{
  // Add embedded elements to the fracture Stencil
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    for( localIndex a=0; a<fvManager.numSubGroups(); ++a )
    {
      FluxApproximationBase * const fluxApprox = fvManager.getGroupPointer< FluxApproximationBase >( a );
      if( fluxApprox!=nullptr )
      {
        fluxApprox->addEDFracToFractureStencil( meshLevel, this->m_fractureRegionName );
      }
    }
  }

}

void EmbeddedSurfaceGenerator::setGlobalIndices( ElementRegionManager & elemManager,
                                                 EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // Add new globalIndices
  int const thisRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const commSize = MpiWrapper::commSize( MPI_COMM_GEOSX );

  localIndex_array numberOfSurfaceElemsPerRank( commSize );
  localIndex_array globalIndexOffset( commSize );
  MpiWrapper::allGather( embeddedSurfaceSubRegion.size(), numberOfSurfaceElemsPerRank );

  globalIndexOffset[0] = 0; // offSet for the globalIndex
  localIndex totalNumberOfSurfaceElements = numberOfSurfaceElemsPerRank[ 0 ];  // Sum across all ranks
  for( int rank = 1; rank < commSize; ++rank )
  {
    globalIndexOffset[rank] = globalIndexOffset[rank - 1] + numberOfSurfaceElemsPerRank[rank - 1];
    totalNumberOfSurfaceElements += numberOfSurfaceElemsPerRank[rank];
  }

  GEOSX_LOG_LEVEL_RANK_0( 1, "Number of embedded surface elements: " << totalNumberOfSurfaceElements );

  arrayView1d< globalIndex > const & localToGlobal = embeddedSurfaceSubRegion.localToGlobalMap();

  for( localIndex ei = 0; ei < embeddedSurfaceSubRegion.size(); ei++ )
  {
    localToGlobal( ei ) = ei + globalIndexOffset[ thisRank ] + elemManager.maxGlobalIndex() + 1;
    embeddedSurfaceSubRegion.updateGlobalToLocalMap( ei );
  }

  embeddedSurfaceSubRegion.setMaxGlobalIndex();

  elemManager.setMaxGlobalIndex();
}

void EmbeddedSurfaceGenerator::addEmbeddedElementsToSets( ElementRegionManager const & elemManager,
                                                          EmbeddedSurfaceSubRegion & embeddedSurfaceSubRegion )
{
  // We want to create the sets for the embeddedSurfaceSubRegion which was empty when they
  // were created for the other subRegions. This way, for example, if the parent cell belongs to the set "source"
  // the embedded element will belong to the same set.

  dataRepository::Group & setGroupEmbSurf =
    embeddedSurfaceSubRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

  elemManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion const & subRegion )
  {
    dataRepository::Group const & setGroupCell = subRegion.getGroup( ObjectManagerBase::groupKeyStruct::setsString() );

    SortedArrayView< localIndex const > const fracturedElements = subRegion.fracturedElementsList();

    ArrayOfArraysView< localIndex const > const cellToEmbSurf = subRegion.embeddedSurfacesList().toViewConst();

    for( localIndex ei = 0; ei < fracturedElements.size(); ei++ )
    {
      localIndex const cellIndex    = fracturedElements[ei];
      localIndex const embSurfIndex = cellToEmbSurf[cellIndex][0];
      setGroupCell.forWrappers< SortedArray< localIndex > >( [&]( auto const & wrapper )
      {
        SortedArrayView< const localIndex > const & targetSetCell = wrapper.reference();
        SortedArray< localIndex > & targetSetEmbSurf =
          setGroupEmbSurf.getWrapper< SortedArray< localIndex > >( wrapper.getName() ).reference();

        if( targetSetCell.contains( cellIndex ) )
        {
          targetSetEmbSurf.insert( embSurfIndex );
        }
      } );
    }
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        EmbeddedSurfaceGenerator,
                        string const &, dataRepository::Group * const )

} /* namespace geosx */
