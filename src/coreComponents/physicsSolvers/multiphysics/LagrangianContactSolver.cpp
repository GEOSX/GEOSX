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
 * @file LagrangianContactSolver.cpp
 *
 */

#include "LagrangianContactSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/solvers/PreconditionerJacobi.hpp"
#include "linearAlgebra/solvers/PreconditionerBlockJacobi.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "math/interpolation/Interpolation.hpp"
#include "finiteElement/elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "finiteElement/elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace interpolation;
using namespace finiteElement;

constexpr integer LagrangianContactSolver::FractureState::STICK;
constexpr integer LagrangianContactSolver::FractureState::SLIP;
constexpr integer LagrangianContactSolver::FractureState::NEW_SLIP;
constexpr integer LagrangianContactSolver::FractureState::OPEN;

LagrangianContactSolver::LagrangianContactSolver( const string & name,
                                                  Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_solidSolver( nullptr ),
  m_stabilizationName(),
  m_contactRelationName(),
  m_activeSetMaxIter()
{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::stabilizationNameString(), &m_stabilizationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the stabilization to use in the lagrangian contact solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the constitutive law used for fracture elements" );

  registerWrapper( viewKeyStruct::activeSetMaxIterString(), &m_activeSetMaxIter ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum number of iteration for the active set strategy in the lagrangian contact solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::lagrangianContactMechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void LagrangianContactSolver::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
    {
      region.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
      {
        subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::rotationMatrixString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the rotation matrices on the fracture." ).
          reference().resizeDimension< 1, 2 >( 3, 3 );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::tractionString() ).
          setApplyDefaultValue( 0.0 ).
          setPlotLevel( PlotLevel::LEVEL_0 ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the tractions on the fracture." ).
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaTractionString() ).
          setApplyDefaultValue( 0.0 ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the traction increments on the fracture." ).
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::fractureStateString() ).
          setPlotLevel( PlotLevel::LEVEL_0 ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the fracture state." );
        initializeFractureState( meshLevel, viewKeyStruct::fractureStateString() );

        subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::previousFractureStateString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the fracture state." );
        initializeFractureState( meshLevel, viewKeyStruct::previousFractureStateString() );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::localJumpString() ).
          setApplyDefaultValue( 0.0 ).
          setPlotLevel( PlotLevel::LEVEL_0 ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the local jump on the fracture at the current time step." ).
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::previousLocalJumpString() ).
          setApplyDefaultValue( 0.0 ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the local jump on the fracture at the previous time step." ).
          reference().resizeDimension< 1 >( 3 );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the normal traction tolerance." );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the normal displacement tolerance." );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::slidingToleranceString() ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName()).
          setDescription( "An array that holds the sliding tolerance." );

        // Needed just because SurfaceGenerator initialize the field "pressure" (NEEDED!!!)
        // It is used in "TwoPointFluxApproximation.cpp", called by "SurfaceGenerator.cpp"
        subRegion.registerWrapper< real64_array >( "pressure" ).
          setPlotLevel( PlotLevel::NOPLOT ).
          setRegisteringObjects( this->getName());
      } );
    } );
  } );
}

void LagrangianContactSolver::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  ConstitutiveBase const & contactRelation = cm.getConstitutiveRelation< ConstitutiveBase >( m_contactRelationName );
  m_contactRelationFullIndex = contactRelation.getIndexInParent();
}

void LagrangianContactSolver::setupSystem( DomainPartition & domain,
                                           DofManager & dofManager,
                                           CRSMatrix< real64, globalIndex > & localMatrix,
                                           array1d< real64 > & localRhs,
                                           array1d< real64 > & localSolution,
                                           bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
  {
    createPreconditioner( domain );
  }
}

void LagrangianContactSolver::implicitStepSetup( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  computeRotationMatrices( domain );
  computeTolerances( domain );
  computeFaceDisplacementJump( domain );

  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void LagrangianContactSolver::implicitStepComplete( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView2d< real64 > const &
      deltaTraction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString() );
      arrayView2d< real64 const > const &
      localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString() );
      arrayView2d< real64 > const &
      previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString() );
      arrayView1d< integer const > const &
      fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );
      arrayView1d< integer > const &
      previousFractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::previousFractureStateString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          deltaTraction[kfe][i] = 0.0;
          previousLocalJump[kfe][i] = localJump[kfe][i];
        }
        previousFractureState[kfe] = fractureState[kfe];
      } );
    }
  } );

  // Need a synchronization of deltaTraction as will be used in AssembleStabilization
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaTractionString() ) );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );

  GEOSX_LOG_LEVEL_RANK_0( 1, " ***** ImplicitStepComplete *****" );
}

void LagrangianContactSolver::postProcessInput()
{
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  SolverBase::postProcessInput();
}

void LagrangianContactSolver::initializePostInitialConditionsPreSubGroups()
{}

LagrangianContactSolver::~LagrangianContactSolver()
{
  // TODO Auto-generated destructor stub
}

void LagrangianContactSolver::computeTolerances( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion;
  arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion;
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex;

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString(),
                                                                                                 m_solidSolver->targetRegionNames(),
                                                                                                 m_solidSolver->solidMaterialNames() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString(),
                                                                                                 m_solidSolver->targetRegionNames(),
                                                                                                 m_solidSolver->solidMaterialNames() );

  using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
    elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
  ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & faceArea = subRegion.getElementArea().toViewConst();
      arrayView3d< real64 const > const &
      faceRotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      arrayView1d< real64 > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
      arrayView1d< real64 > const & normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );
      arrayView1d< real64 > const & slidingTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          real64 const area = faceArea[kfe];
          // approximation of the stiffness along coordinate directions
          // ( first, second ) index -> ( element index, direction )
          // 1. T -> top (index 0), B -> bottom (index 1)
          // 2. the coordinate direction (x, y, z)
          real64 stiffDiagApprox[ 2 ][ 3 ];
          real64 averageYoungModulus = 0.0;
          real64 averageConstrainedModulus = 0.0;
          real64 averageBoxSize0 = 0.0;

          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex const faceIndex = elemsToFaces[kfe][i];
            localIndex const er = faceToElemRegion[faceIndex][0];
            localIndex const esr = faceToElemSubRegion[faceIndex][0];
            localIndex const ei = faceToElemIndex[faceIndex][0];

            real64 const volume = elemVolume[er][esr][ei];

            // Get the "element to node" map for the specific region/subregion
            NodeMapViewType const & cellElemsToNodes = elemToNodeView[er][esr];
            localIndex const numNodesPerElem = cellElemsToNodes.size( 1 );

            // Compute the box size
            real64 maxSize[3];
            real64 minSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            }
            for( localIndex a = 1; a < numNodesPerElem; ++a )
            {
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = fmax( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                minSize[j] = fmin( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              }
            }

            real64 boxSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );
            real64 const M = K + 4.0 / 3.0 * G;

            // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
            for( localIndex j = 0; j < 3; ++j )
            {
              stiffDiagApprox[ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 4.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
            }

            averageYoungModulus += 0.5*E;
            averageConstrainedModulus += 0.5*M;
            averageBoxSize0 += 0.5*boxSize[0];
          }

          // Average the stiffness and compute the inverse
          real64 invStiffApprox[ 3 ][ 3 ] = { { 0 } };
          for( localIndex j = 0; j < 3; ++j )
          {
            invStiffApprox[ j ][ j ] = ( stiffDiagApprox[ 0 ][ j ] + stiffDiagApprox[ 1 ][ j ] ) / ( stiffDiagApprox[ 0 ][ j ] * stiffDiagApprox[ 1 ][ j ] );
          }

          // Rotate in the local reference system, computing R^T * (invK) * R
          real64 temp[ 3 ][ 3 ];
          LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, faceRotationMatrix[ kfe ], invStiffApprox );
          real64 rotatedInvStiffApprox[ 3 ][ 3 ];
          LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, faceRotationMatrix[ kfe ] );
          LvArray::tensorOps::scale< 3, 3 >( rotatedInvStiffApprox, area );

          // Finally, compute tolerances for the given fracture element
          normalDisplacementTolerance[kfe] = rotatedInvStiffApprox[ 0 ][ 0 ] * averageYoungModulus / 2.e+7;
          slidingTolerance[kfe] = sqrt( rotatedInvStiffApprox[ 1 ][ 1 ] * rotatedInvStiffApprox[ 1 ][ 1 ] +
                                        rotatedInvStiffApprox[ 2 ][ 2 ] * rotatedInvStiffApprox[ 2 ][ 2 ] ) * averageYoungModulus / 2.e+7;
          normalTractionTolerance[kfe] = 1.0 / 2.0 * averageConstrainedModulus / averageBoxSize0 * normalDisplacementTolerance[kfe];
        }
      } );
    }
  } );
}

void LagrangianContactSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_solidSolver->resetStateToBeginningOfStep( domain );

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView2d< real64 > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString() );
      arrayView2d< real64 > const & deltaTraction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString() );
      arrayView2d< real64 > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString() );
      arrayView2d< real64 const > const & previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString() );

      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );
      arrayView1d< integer const > const & previousFractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::previousFractureStateString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        for( localIndex i = 0; i < 3; ++i )
        {
          traction[kfe][i] -= deltaTraction[kfe][i];
          deltaTraction[kfe][i] = 0.0;

          localJump[kfe][i] = previousLocalJump[kfe][i];
        }
        fractureState[kfe] = previousFractureState[kfe];
      } );
    }
  } );
}

real64 LagrangianContactSolver::solverStep( real64 const & time_n,
                                            real64 const & dt,
                                            int const cycleNumber,
                                            DomainPartition & domain )
{
  real64 dtReturn = dt;

  implicitStepSetup( time_n,
                     dt,
                     domain );

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_localRhs,
               m_localSolution );

  // currently the only method is implicit time integration
  dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. Typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void LagrangianContactSolver::computeFaceDisplacementJump( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = meshLevel.getNodeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager.totalDisplacement();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView2d< real64 > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString() );
      arrayView1d< real64 const > const & area = subRegion.getElementArea().toViewConst();

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        // Contact constraints
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );

        array1d< real64 > nodalArea0, nodalArea1;
        computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][0], nodalArea0 );
        computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][1], nodalArea1 );

        real64 globalJumpTemp[ 3 ] = { 0 };
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( localIndex i = 0; i < 3; ++i )
          {
            globalJumpTemp[ i ] +=
              ( -u[faceToNodeMap( elemsToFaces[kfe][0], a )][i] * nodalArea0[a]
                +u[faceToNodeMap( elemsToFaces[kfe][1], a )][i] * nodalArea1[a] ) / area[kfe];
          }
        }

        real64 localJumpTemp[ 3 ];
        LvArray::tensorOps::Ri_eq_AjiBj< 3, 3 >( localJumpTemp, rotationMatrix[ kfe ], globalJumpTemp );
        LvArray::tensorOps::copy< 3 >( localJump[ kfe ], localJumpTemp );
      } );
    }
  } );

  return;
}

real64 LagrangianContactSolver::explicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & dt,
                                              const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                              DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_MARK_FUNCTION;
  GEOSX_ERROR( "ExplicitStep non available for LagrangianContactSolver!" );
  return dt;
}

real64 LagrangianContactSolver::nonlinearImplicitStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       integer const cycleNumber,
                                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged = m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numdtAttempts;

  // a flag to denote whether we have converged
  bool isNewtonConverged = false;

  bool isActiveSetConverged = false;

  // outer loop attempts to apply full timestep, and manages timestep cut if required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      resetStateToBeginningOfStep( domain );
      globalIndex numStick, numSlip, numOpen;
      computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
    }

    bool useElasticStep = !isFractureAllInStickCondition( domain );

    integer & activeSetIter = m_activeSetIter;
    for( activeSetIter = 0; activeSetIter < m_activeSetMaxIter; ++activeSetIter )
    {
      // *******************************
      // Newton loop: begin
      // *******************************
      isNewtonConverged = false;
      // keep residual from previous iteration in case we need to do a line search
      real64 lastResidual = 1e99;
      integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
      real64 scaleFactor = 1.0;

      // main Newton loop
      bool computeResidual = true;
      for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
      {
        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          char output[55] = {0};
          sprintf( output, "    Attempt: %2d, ActiveSetIter: %2d ; NewtonIter: %2d ; ",
                   dtAttempt, activeSetIter, newtonIter );
          std::cout<<output<<std::endl;
        }

        // zero out matrix/rhs before assembly
        m_localMatrix.zero();
        m_localRhs.zero();

        // call assemble to fill the matrix and the rhs
        assembleSystem( time_n,
                        stepDt,
                        domain,
                        m_dofManager,
                        m_localMatrix.toViewConstSizes(),
                        m_localRhs );

        // apply boundary conditions to system
        applyBoundaryConditions( time_n,
                                 stepDt,
                                 domain,
                                 m_dofManager,
                                 m_localMatrix.toViewConstSizes(),
                                 m_localRhs );

        // TODO: maybe add scale function here?
        // Scale()

        real64 residualNorm;
        // get residual norm
        if( computeResidual )
        {
          residualNorm = calculateResidualNorm( domain, m_dofManager, m_localRhs );
        }
        else
        {
          residualNorm = lastResidual;
        }

        if( getLogLevel() >= 1 && logger::internal::rank==0 )
        {
          if( newtonIter!=0 )
          {
            char output[46] = {0};
            sprintf( output,
                     "Last LinSolve(iter,tol) = (%4d, %4.2e) ; ",
                     m_linearSolverResult.numIterations,
                     m_linearSolverResult.residualReduction );
            std::cout<<output;
          }
          std::cout<<std::endl;
        }

        // if the residual norm is less than the Newton tolerance we denote that we have
        // converged and break from the Newton loop immediately.
        if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
        {
          isNewtonConverged = true;
          break;
        }

        // if using adaptive Krylov tolerance scheme, update tolerance.
        LinearSolverParameters::Krylov & krylovParams = m_linearSolverParameters.get().krylov;
        if( krylovParams.useAdaptiveTol )
        {
          krylovParams.relTolerance = eisenstatWalker( residualNorm, lastResidual, krylovParams.weakestTol );
        }

        // Compose parallel LA matrix/rhs out of local LA matrix/rhs
        m_matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
        m_rhs.create( m_localRhs, MPI_COMM_GEOSX );
        m_solution.createWithLocalSize( m_matrix.numLocalCols(), MPI_COMM_GEOSX );

        // Output the linear system matrix/rhs for debugging purposes
        debugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

        // Solve the linear system
        solveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

        // Output the linear system solution for debugging purposes
        debugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

        // Copy solution from parallel vector back to local
        // TODO: This step will not be needed when we teach LA vectors to wrap our pointers
        m_solution.extract( m_localSolution );

        scaleFactor = scalingForSystemSolution( domain, m_dofManager, m_localSolution );

        // do line search in case residual has increased
        if( m_nonlinearSolverParameters.m_lineSearchAction != NonlinearSolverParameters::LineSearchAction::None && newtonIter > 0 )
        {
          bool lineSearchSuccess = lineSearch( time_n,
                                               stepDt,
                                               cycleNumber,
                                               domain,
                                               m_dofManager,
                                               m_localMatrix.toViewConstSizes(),
                                               m_localRhs,
                                               m_localSolution,
                                               scaleFactor,
                                               residualNorm );

          if( !lineSearchSuccess )
          {
            if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Attempt )
            {
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Accepting iteration." );
            }
            else if( m_nonlinearSolverParameters.m_lineSearchAction == NonlinearSolverParameters::LineSearchAction::Require )
            {
              // if line search failed, then break out of the main Newton loop. Timestep will be cut.
              GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search failed to produce reduced residual. Exiting Newton Loop." );
              break;
            }
          }
          // Residual norm already computed in line search and stored in "residualNorm"
          computeResidual = false;
        }
        else
        {
          // apply the system solution to the fields/variables
          applySystemSolution( m_dofManager, m_localSolution, scaleFactor, domain );
          // Need to compute the residual norm
          computeResidual = true;
        }

        if( !checkSystemSolution( domain, m_dofManager, m_localSolution, scaleFactor ) )
        {
          // TODO try chopping (similar to line search)
          GEOSX_LOG_RANK_0( "    Solution check failed. Newton loop terminated." );
          break;
        }

        lastResidual = residualNorm;
      }
      // *******************************
      // Newton loop: end
      // *******************************

      // *******************************
      // Active set check: begin
      // *******************************
      bool const isPreviousFractureStateValid = updateFractureState( domain );
      GEOSX_LOG_LEVEL_RANK_0( 1, "active set flag: " << std::boolalpha << isPreviousFractureStateValid );

      if( getLogLevel() >= 1 )
      {
        globalIndex numStick, numSlip, numOpen;
        computeFractureStateStatistics( domain, numStick, numSlip, numOpen, true );
      }
      // *******************************
      // Active set check: end
      // *******************************

      GEOSX_LOG_LEVEL_RANK_0( 1, "isPreviousFractureStateValid: " << std::boolalpha << isPreviousFractureStateValid <<
                              " | isNewtonConverged: " << isNewtonConverged << " | useElasticStep: " << useElasticStep );
      if( isNewtonConverged )
      {
        isActiveSetConverged = isPreviousFractureStateValid;
        if( isActiveSetConverged )
        {
          break;
        }
      }
      else if( useElasticStep )
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "Trying with an elastic step" );
        useElasticStep = false;
        resetStateToBeginningOfStep( domain );
        setFractureStateForElasticStep( domain );
      }
      else
      {
        GEOSX_LOG_LEVEL_RANK_0( 1, "Newton did not converge in active set loop" );
        break;
      }
    }
    if( !isNewtonConverged )
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      stepDt *= dtCutFactor;
      GEOSX_LOG_LEVEL_RANK_0 ( 1, "New dt = " <<  stepDt );
    }
    if( isActiveSetConverged )
    {
      break;
    }
  }

  if( !isNewtonConverged )
  {
    GEOSX_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOSX_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOSX_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  if( !isActiveSetConverged )
  {
    GEOSX_ERROR( "Active set did not reached a solution. Terminating..." );
  }
  else
  {
    GEOSX_LOG_RANK_0( "Number of active set iterations: " << m_activeSetIter );
  }

  // return the achieved timestep
  return stepDt;
}

bool LagrangianContactSolver::lineSearch( real64 const & time_n,
                                          real64 const & dt,
                                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs,
                                          arrayView1d< real64 const > const & localSolution,
                                          real64 const scaleFactor,
                                          real64 & lastResidual )
{
  bool lineSearchSuccess = true;

  integer const maxNumberLineSearchCuts = m_nonlinearSolverParameters.m_lineSearchMaxCuts;

  real64 const sigma1 = 0.5;
  real64 const alpha = 1.e-4;

  real64 localScaleFactor = scaleFactor;
  real64 lamm = scaleFactor;
  real64 lamc = localScaleFactor;
  integer lineSearchIteration = 0;

  // get residual norm
  real64 residualNorm0 = lastResidual;

  applySystemSolution( dofManager, localSolution, scaleFactor, domain );

  // re-assemble system
  localMatrix.zero();
  localRhs.zero();
  assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // apply boundary conditions to system
  applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // get residual norm
  real64 residualNormT = calculateResidualNorm( domain, dofManager, localRhs );

  real64 ff0 = residualNorm0*residualNorm0;
  real64 ffT = residualNormT*residualNormT;
  real64 ffm = ffT;
  real64 cumulativeScale = scaleFactor;

  while( residualNormT >= (1.0 - alpha*localScaleFactor)*residualNorm0 )
  {
    real64 const previousLocalScaleFactor = localScaleFactor;
    // Apply the three point parabolic model
    if( lineSearchIteration == 0 )
    {
      localScaleFactor *= sigma1;
    }
    else
    {
      localScaleFactor = ParabolicInterpolationThreePoints( lamc, lamm, ff0, ffT, ffm );
    }

    // Update x; keep the books on lambda
    real64 const deltaLocalScaleFactor = ( localScaleFactor - previousLocalScaleFactor );
    cumulativeScale += deltaLocalScaleFactor;

    if( !checkSystemSolution( domain, dofManager, localSolution, deltaLocalScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "        Line search " << lineSearchIteration << ", solution check failed" );
      continue;
    }

    applySystemSolution( dofManager, localSolution, deltaLocalScaleFactor, domain );
    lamm = lamc;
    lamc = localScaleFactor;

    // Keep the books on the function norms
    // re-assemble system
    // TODO: add a flag to avoid a completely useless Jacobian computation: rhs is enough
    localMatrix.zero();
    localRhs.zero();
    assembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

    // apply boundary conditions to system
    applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout<<output;
    }

    // get residual norm
    residualNormT = calculateResidualNorm( domain, dofManager, localRhs );
    ffm = ffT;
    ffT = residualNormT*residualNormT;
    lineSearchIteration += 1;

    if( lineSearchIteration > maxNumberLineSearchCuts )
    {
      lineSearchSuccess = false;
      break;
    }
  }

  lastResidual = residualNormT;

  return lineSearchSuccess;
}

void LagrangianContactSolver::setupDofs( DomainPartition const & domain,
                                         DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );

  // restrict coupling to fracture regions only
  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addField( viewKeyStruct::tractionString(),
                       DofManager::Location::Elem,
                       3,
                       fractureRegions );
  dofManager.addCoupling( viewKeyStruct::tractionString(),
                          viewKeyStruct::tractionString(),
                          DofManager::Connector::Face,
                          fractureRegions );
  dofManager.addCoupling( keys::TotalDisplacement,
                          viewKeyStruct::tractionString(),
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void LagrangianContactSolver::assembleSystem( real64 const time,
                                              real64 const dt,
                                              DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  synchronizeFractureState( domain );

  m_solidSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  assembleForceResidualDerivativeWrtTraction( domain, dofManager, localMatrix, localRhs );
  assembleTractionResidualDerivativeWrtDisplacementAndTraction( domain, dofManager, localMatrix, localRhs );
  assembleStabilization( domain, dofManager, localMatrix, localRhs );
}

void LagrangianContactSolver::applyBoundaryConditions( real64 const time,
                                                       real64 const dt,
                                                       DomainPartition & domain,
                                                       DofManager const & dofManager,
                                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                       arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->applyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );
}

real64 LagrangianContactSolver::calculateResidualNorm( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();

  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( keys::TotalDisplacement ) );

  string const & dofKey = dofManager.getKey( viewKeyStruct::tractionString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< integer const > const & elemGhostRank = nodeManager.ghostRank();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum0( 0.0 );
  forAll< parallelDevicePolicy<> >( nodeManager.size(),
                                    [localRhs, localSum0, dispDofNumber, rankOffset, elemGhostRank] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    if( elemGhostRank[k] < 0 )
    {
      localIndex const localRow = LvArray::integerConversion< localIndex >( dispDofNumber[k] - rankOffset );
      for( localIndex dim = 0; dim < 3; ++dim )
      {
        localSum0 += localRhs[localRow + dim] * localRhs[localRow + dim];
      }
    }
  } );
  real64 const momentumR2 = localSum0.get();

  real64 contactR2 = 0.0;

  forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const, FaceElementSubRegion const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

    RAJA::ReduceSum< parallelHostReduce, real64 > localSum( 0.0 );
    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const k )
    {
      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
        for( localIndex dim = 0; dim < 3; ++dim )
        {
          localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
        }
      }
    } );
    contactR2 += localSum.get();
  } );

  real64 localR2[2] = { momentumR2, contactR2 };
  real64 globalResidualNorm[3]{};

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  array1d< real64 > globalR2( 2 * size );
  globalR2.zero();

  // Everything is done on rank 0
  MpiWrapper::gather( localR2,
                      2,
                      globalR2.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    globalResidualNorm[0] = 0.0;
    globalResidualNorm[1] = 0.0;
    for( int r=0; r<size; ++r )
    {
      // sum across all ranks
      globalResidualNorm[0] += globalR2[2 * r + 0];
      globalResidualNorm[1] += globalR2[2 * r + 1];
    }
    globalResidualNorm[2] = globalResidualNorm[0] + globalResidualNorm[1];
    globalResidualNorm[0] = sqrt( globalResidualNorm[0] );
    globalResidualNorm[1] = sqrt( globalResidualNorm[1] );
    globalResidualNorm[2] = sqrt( globalResidualNorm[2] );
  }

  MpiWrapper::bcast( globalResidualNorm, 3, 0, MPI_COMM_GEOSX );

  if( m_nonlinearSolverParameters.m_numNewtonIterations == 0 )
  {
    m_initialResidual[0] = globalResidualNorm[0];
    m_initialResidual[1] = globalResidualNorm[1];
    m_initialResidual[2] = globalResidualNorm[2];
    globalResidualNorm[0] = 1.0;
    globalResidualNorm[1] = 1.0;
    globalResidualNorm[2] = 1.0;
  }
  else
  {
    globalResidualNorm[0] /= (m_initialResidual[0]+1.0);
    globalResidualNorm[1] /= (m_initialResidual[1]+1.0);
    // Add 0 just to match Matlab code results
    globalResidualNorm[2] /= (m_initialResidual[2]+1.0);
  }

  char output[94] = {0};
  sprintf( output,
           "( Rdisplacement, Rtraction, Rtotal ) = ( %15.6e, %15.6e, %15.6e );",
           globalResidualNorm[0],
           globalResidualNorm[1],
           globalResidualNorm[2] );
  GEOSX_LOG_LEVEL_RANK_0( 1, output );

  return globalResidualNorm[2];
}

void LagrangianContactSolver::createPreconditioner( DomainPartition const & domain )
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    // TODO: move among inputs (xml)
    string const leadingBlockApproximation = "blockJacobi";

    LinearSolverParameters mechParams = m_solidSolver->getLinearSolverParameters();
    // Because of boundary conditions
    mechParams.isSymmetric = false;

    std::unique_ptr< BlockPreconditioner< LAInterface > > precond;
    std::unique_ptr< PreconditionerBase< LAInterface > > tracPrecond;

    if( leadingBlockApproximation == "jacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockDiagonal,
                                                                        BlockScalingOption::UserProvided );
      // Using GEOSX implementation of Jacobi preconditioner
      // tracPrecond = std::make_unique< PreconditionerJacobi< LAInterface > >();

      // Using LAI implementation of Jacobi preconditioner
      LinearSolverParameters tracParams;
      tracParams.preconditionerType = LinearSolverParameters::PreconditionerType::jacobi;
      tracPrecond = LAInterface::createPreconditioner( tracParams );
    }
    else if( leadingBlockApproximation == "blockJacobi" )
    {
      precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::LowerUpperTriangular,
                                                                        SchurComplementOption::FirstBlockUserDefined,
                                                                        BlockScalingOption::UserProvided );
      tracPrecond = std::make_unique< PreconditionerBlockJacobi< LAInterface > >( mechParams.dofsPerNode );
    }
    else
    {
      GEOSX_ERROR( "LagrangianContactSolver::CreatePreconditioner leadingBlockApproximation option " << leadingBlockApproximation << " not supported" );
    }

    // Preconditioner for the leading block: tracPrecond
    precond->setupBlock( 0,
                         { { viewKeyStruct::tractionString(), { 3, true } } },
                         std::move( tracPrecond ) );

    if( mechParams.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
    {
      if( m_solidSolver->getRigidBodyModes().empty() )
      {
        MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
        LAIHelperFunctions::computeRigidBodyModes( mesh,
                                                   m_dofManager,
                                                   { keys::TotalDisplacement },
                                                   m_solidSolver->getRigidBodyModes() );
      }
    }

    // Preconditioner for the Schur complement: mechPrecond
    std::unique_ptr< PreconditionerBase< LAInterface > > mechPrecond = LAInterface::createPreconditioner( mechParams, m_solidSolver->getRigidBodyModes() );
    precond->setupBlock( 1,
                         { { keys::TotalDisplacement, { 3, true } } },
                         std::move( mechPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void LagrangianContactSolver::computeRotationMatrices( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  ElementRegionManager & elemManager = mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      arrayView3d< real64 > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const kfe )
      {
        stackArray1d< real64, 3 > Nbar( 3 );
        Nbar[ 0 ] = faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0];
        Nbar[ 1 ] = faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1];
        Nbar[ 2 ] = faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2];
        LvArray::tensorOps::normalize< 3 >( Nbar );

        computationalGeometry::RotationMatrix_3D( Nbar.toSliceConst(), rotationMatrix[kfe] );
      } );
    }
  } );
}

void LagrangianContactSolver::computeFaceNodalArea( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
                                                    ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                                                    localIndex const kf0,
                                                    array1d< real64 > & nodalArea ) const
{
  // I've tried to access the finiteElement::dispatch3D with
  // finiteElement::FiniteElementBase const &
  // fe = fractureSubRegion->getReference< finiteElement::FiniteElementBase >( surfaceGenerator->getDiscretizationName() );
  // but it's either empty (unknown discretization) or for 3D only (e.g., hexahedra)
  GEOSX_MARK_FUNCTION;

  localIndex const TriangularPermutation[3] = { 0, 1, 2 };
  localIndex const QuadrilateralPermutation[4] = { 0, 1, 3, 2 };

  localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

  nodalArea.resize( numNodesPerFace );
  for( localIndex a = 0; a < numNodesPerFace; ++a )
  {
    nodalArea[a] = 0.0;
  }
  localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
  if( numNodesPerFace == 3 )
  {
    real64 xLocal[3][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[3];
    for( localIndex q=0; q<H1_TriangleFace_Lagrange1_Gauss1::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, xLocal );
      H1_TriangleFace_Lagrange1_Gauss1::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        nodalArea[a] += detJ * N[permutation[a]];
      }
    }
  }
  else if( numNodesPerFace == 4 )
  {
    real64 xLocal[4][3];
    for( localIndex a = 0; a < numNodesPerFace; ++a )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        xLocal[a][j] = nodePosition[faceToNodeMap( kf0, permutation[a] )][j];
      }
    }
    real64 N[4];
    for( localIndex q=0; q<H1_QuadrilateralFace_Lagrange1_GaussLegendre2::numQuadraturePoints; ++q )
    {
      real64 const detJ = H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal );
      H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
      for( localIndex a = 0; a < numNodesPerFace; ++a )
      {
        nodalArea[a] += detJ * N[permutation[a]];
      }
    }
  }
  else
  {
    GEOSX_ERROR( "LagrangianContactSolver: face with " << numNodesPerFace << " nodes. Only triangles and quadrilaterals are supported." );
  }
}

void LagrangianContactSolver::
  assembleForceResidualDerivativeWrtTraction( DomainPartition & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const & tracDofKey = dofManager.getKey( viewKeyStruct::tractionString() );
  string const & dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString() );
      arrayView3d< real64 const > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      constexpr localIndex TriangularPermutation[3] = { 0, 1, 2 };
      constexpr localIndex QuadrilateralPermutation[4] = { 0, 1, 3, 2 };

      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
        localIndex const numQuadraturePointsPerElem = numNodesPerFace==3 ? 1 : 4;

        globalIndex rowDOF[12];
        real64 nodeRHS[12];
        stackArray2d< real64, 3*4*3 > dRdT( 3*numNodesPerFace, 3 );
        globalIndex colDOF[3];
        for( localIndex i = 0; i < 3; ++i )
        {
          colDOF[i] = tracDofNumber[kfe] + i;
        }

        localIndex const * const permutation = ( numNodesPerFace == 3 ) ? TriangularPermutation : QuadrilateralPermutation;
        real64 xLocal[2][4][3];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];
          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            for( localIndex j = 0; j < 3; ++j )
            {
              xLocal[kf][a][j] = nodePosition[ faceToNodeMap( faceIndex, permutation[a] ) ][j];
            }
          }
        }

        real64 N[4];

        for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
        {
          if( numNodesPerFace==3 )
          {
            using NT = real64[3];
            H1_TriangleFace_Lagrange1_Gauss1::calcN( q, reinterpret_cast< NT & >(N) );
          }
          else if( numNodesPerFace==4 )
          {
            H1_QuadrilateralFace_Lagrange1_GaussLegendre2::calcN( q, N );
          }

          constexpr int normalSign[2] = { 1, -1 };
          for( localIndex kf = 0; kf < 2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];
            using xLocalTriangle = real64[3][3];
            real64 const detJxW = numNodesPerFace==3 ?
                                  H1_TriangleFace_Lagrange1_Gauss1::transformedQuadratureWeight( q, reinterpret_cast< xLocalTriangle & >( xLocal[kf] ) ) :
                                  H1_QuadrilateralFace_Lagrange1_GaussLegendre2::transformedQuadratureWeight( q, xLocal[kf] );

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              real64 const NaDetJxQ = N[permutation[a]] * detJxW;
              real64 const localNodalForce[ 3 ] = { traction( kfe, 0 ) * NaDetJxQ,
                                                    traction( kfe, 1 ) * NaDetJxQ,
                                                    traction( kfe, 2 ) * NaDetJxQ };
              real64 globalNodalForce[ 3 ];
              LvArray::tensorOps::Ri_eq_AijBj< 3, 3 >( globalNodalForce, rotationMatrix[ kfe ], localNodalForce );

              for( localIndex i = 0; i < 3; ++i )
              {
                rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
                // Opposite sign w.r.t. to formulation presented in
                // Algebraically Stabilized Lagrange Multiplier Method for Frictional Contact Mechanics with
                // Hydraulically Active Fractures
                // Franceschini, A., Castelletto, N., White, J. A., Tchelepi, H. A.
                // Computer Methods in Applied Mechanics and Engineering (2020) 368, 113161
                // doi: 10.1016/j.cma.2020.113161
                nodeRHS[3*a+i] = +globalNodalForce[i] * normalSign[ kf ];

                // Opposite sign w.r.t. to the same formulation as above
                dRdT( 3*a+i, 0 ) = rotationMatrix( kfe, i, 0 ) * normalSign[ kf ] * NaDetJxQ;
                dRdT( 3*a+i, 1 ) = rotationMatrix( kfe, i, 1 ) * normalSign[ kf ] * NaDetJxQ;
                dRdT( 3*a+i, 2 ) = rotationMatrix( kfe, i, 2 ) * normalSign[ kf ] * NaDetJxQ;
              }
            }

            for( localIndex idof = 0; idof < numNodesPerFace * 3; ++idof )
            {
              localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

              if( localRow >= 0 && localRow < localMatrix.numRows() )
              {
                localMatrix.addToRow< parallelHostAtomic >( localRow,
                                                            colDOF,
                                                            dRdT[idof].dataIfContiguous(),
                                                            3 );
                RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow], nodeRHS[idof] );
              }
            }
          }
        }
      } );
    }
  } );
}

void LagrangianContactSolver::
  assembleTractionResidualDerivativeWrtDisplacementAndTraction( DomainPartition const & domain,
                                                                DofManager const & dofManager,
                                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();
  ContactRelationBase const & contactRelation = constitutiveManager.getGroup< ContactRelationBase const >( m_contactRelationName );

  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const & tracDofKey = dofManager.getKey( viewKeyStruct::tractionString() );
  string const & dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< globalIndex const > const & tracDofNumber = subRegion.getReference< globalIndex_array >( tracDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView3d< real64 const > const &
      rotationMatrix = subRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString() );
      arrayView1d< integer const > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );
      arrayView2d< real64 const > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString() );
      arrayView2d< real64 const > const & previousLocalJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::previousLocalJumpString() );
      arrayView1d< real64 const > const & slidingTolerance = subRegion.getReference< array1d< real64 > >( viewKeyStruct::slidingToleranceString() );

      forAll< parallelHostPolicy >( subRegion.size(), [=, &contactRelation] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[kfe][0] );
          globalIndex nodeDOF[24];
          globalIndex elemDOF[3];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[i] = tracDofNumber[kfe] + i;
          }

          real64 elemRHS[3] = {0.0, 0.0, 0.0};
          real64 const Ja = area[kfe];

          stackArray2d< real64, 2 * 3 * 4 * 3 > dRdU( 3, 2 * 3 * numNodesPerFace );
          stackArray2d< real64, 3 * 3 > dRdT( 3, 3 );

          switch( fractureState[kfe] )
          {
            case FractureState::STICK:
              {
                for( localIndex i = 0; i < 3; ++i )
                {
                  if( i == 0 )
                  {
                    elemRHS[i] = +Ja * localJump[kfe][i];
                  }
                  else
                  {
                    elemRHS[i] = +Ja * ( localJump[kfe][i] - previousLocalJump[kfe][i] );
                  }
                }

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  array1d< real64 > nodalArea;
                  computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] + i;

                      dRdU( 0, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                      dRdU( 1, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 1 ) * pow( -1, kf );
                      dRdU( 2, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 2 ) * pow( -1, kf );
                    }
                  }
                }
                break;
              }
            case FractureState::SLIP:
            case FractureState::NEW_SLIP:
              {
                elemRHS[0] = +Ja * localJump[kfe][0];

                for( localIndex kf = 0; kf < 2; ++kf )
                {
                  // Compute local area contribution for each node
                  array1d< real64 > nodalArea;
                  computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                  for( localIndex a = 0; a < numNodesPerFace; ++a )
                  {
                    for( localIndex i = 0; i < 3; ++i )
                    {
                      nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[kfe][kf], a )] +
                                                                      LvArray::integerConversion< globalIndex >( i );
                      dRdU( 0, kf * 3 * numNodesPerFace + 3 * a + i ) = -nodalArea[a] * rotationMatrix( kfe, i, 0 ) * pow( -1, kf );
                    }
                  }
                }

                real64 const limitTau = contactRelation.limitTangentialTractionNorm( traction[kfe][0] );
                real64 sliding[ 2 ] = { localJump[kfe][1] - previousLocalJump[kfe][1], localJump[kfe][2] - previousLocalJump[kfe][2] };
                real64 slidingNorm = sqrt( sliding[ 0 ]*sliding[ 0 ] + sliding[ 1 ]*sliding[ 1 ] );

//                GEOSX_LOG_LEVEL_BY_RANK( 3, "element: " << kfe << " sliding: " << sliding[0] << " " << sliding[1] );

                if( !( ( m_nonlinearSolverParameters.m_numNewtonIterations == 0 ) && ( fractureState[kfe] == FractureState::NEW_SLIP ) )
                    && slidingNorm > slidingTolerance[kfe] )
                {
                  for( localIndex i = 1; i < 3; ++i )
                  {
                    elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * sliding[ i-1 ] / slidingNorm );
                  }

                  // A symmetric 2x2 matrix.
                  real64 dUdgT[ 3 ];
                  dUdgT[ 0 ] = (slidingNorm * slidingNorm - sliding[ 0 ] * sliding[ 0 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 1 ] = (slidingNorm * slidingNorm - sliding[ 1 ] * sliding[ 1 ]) * limitTau / std::pow( slidingNorm, 3 );
                  dUdgT[ 2 ] = -sliding[ 0 ] * sliding[ 1 ] * limitTau / std::pow( slidingNorm, 3 );

                  for( localIndex kf = 0; kf < 2; ++kf )
                  {
                    // Compute local area contribution for each node
                    array1d< real64 > nodalArea;
                    computeFaceNodalArea( nodePosition, faceToNodeMap, elemsToFaces[kfe][kf], nodalArea );

                    for( localIndex a = 0; a < numNodesPerFace; ++a )
                    {
                      for( localIndex i = 0; i < 3; ++i )
                      {
                        real64 const localRowB[ 2 ] = { rotationMatrix( kfe, i, 1 ), rotationMatrix( kfe, i, 2 ) };
                        real64 localRowE[ 2 ];
                        LvArray::tensorOps::Ri_eq_symAijBj< 2 >( localRowE, dUdgT, localRowB );

                        dRdU( 1, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 0 ] * pow( -1, kf );
                        dRdU( 2, kf * 3 * numNodesPerFace + 3 * a + i ) = nodalArea[a] * localRowE[ 1 ] * pow( -1, kf );
                      }
                    }
                  }
                  for( localIndex i = 1; i < 3; ++i )
                  {
                    dRdT( i, 0 ) = Ja * contactRelation.dLimitTangentialTractionNorm_dNormalTraction( traction[kfe][0] ) * sliding[ i-1 ] / slidingNorm;
                    dRdT( i, i ) = Ja;
                  }
                }
                else
                {
                  real64 vaux[ 2 ] = { traction[kfe][1], traction[kfe][2] };
                  real64 vauxNorm = sqrt( vaux[ 0 ]*vaux[ 0 ] + vaux[ 1 ]*vaux[ 1 ] );
                  if( vauxNorm > 0.0 )
                  {
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = +Ja * ( traction[kfe][i] - limitTau * vaux[ i-1 ] / vauxNorm );
                    }
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                  else
                  {
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      elemRHS[i] = 0.0;
                    }
                    for( localIndex i = 1; i < 3; ++i )
                    {
                      dRdT( i, i ) = Ja;
                    }
                  }
                }
                break;
              }
            case FractureState::OPEN:
              {
//                GEOSX_LOG_LEVEL_BY_RANK( 3, "element: " << kfe << " opening: " << localJump[kfe][0] );

                for( localIndex i = 0; i < 3; ++i )
                {
                  elemRHS[i] = +Ja * traction[kfe][i];
                }

                for( localIndex i = 0; i < 3; ++i )
                {
                  dRdT( i, i ) = Ja;
                }
                break;
              }
          }

          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[0] - rankOffset );

          for( localIndex idof = 0; idof < 3; ++idof )
          {
            localRhs[localRow + idof] += elemRHS[idof];

            if( fractureState[kfe] != FractureState::OPEN )
            {
              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow + idof,
                                                                        nodeDOF,
                                                                        dRdU[idof].dataIfContiguous(),
                                                                        2 * 3 * numNodesPerFace );
            }

            if( fractureState[kfe] != FractureState::STICK )
            {
              localMatrix.addToRow< serialAtomic >( localRow + idof,
                                                    elemDOF,
                                                    dRdT[idof].dataIfContiguous(),
                                                    3 );
            }
          }
        }
      } );
    }
  } );
}

void LagrangianContactSolver::assembleStabilization( DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const & tracDofKey = dofManager.getKey( viewKeyStruct::tractionString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // Get the finite volume method used to compute the stabilization
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & stabilizationMethod = fvManager.getFluxApproximation( m_stabilizationName );

  // Get the "face to element" map (valid for the entire mesh)
  FaceManager::ElemMapType const & faceToElem = faceManager.toElementRelation();
  arrayView2d< localIndex const > const & faceToElemRegion = faceToElem.m_toElementRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemSubRegion = faceToElem.m_toElementSubRegion.toViewConst();
  arrayView2d< localIndex const > const & faceToElemIndex = faceToElem.m_toElementIndex.toViewConst();

  // Form the SurfaceGenerator, get the fracture name and use it to retrieve the faceMap (from fracture element to face)
  SurfaceGenerator const &
  surfaceGenerator = this->getParent().getGroup< SurfaceGenerator >( "SurfaceGen" );
  SurfaceElementRegion const & fractureRegion = elemManager.getRegion< SurfaceElementRegion >( surfaceGenerator.getFractureRegionName() );
  FaceElementSubRegion const & fractureSubRegion = fractureRegion.getSubRegion< FaceElementSubRegion >( "faceElementSubRegion" );
  GEOSX_ERROR_IF( !fractureSubRegion.hasWrapper( m_tractionKey ), "The fracture subregion must contain traction field." );
  arrayView2d< localIndex const > const faceMap = fractureSubRegion.faceList();
  GEOSX_ERROR_IF( faceMap.size( 1 ) != 2, "A fracture face has to be shared by two cells." );

  // Get the state of fracture elements
  arrayView1d< integer const > const & fractureState =
    fractureSubRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );

  // Get the tractions and stabilization contribution to the local jump
  arrayView2d< real64 const > const & traction =
    fractureSubRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString() );
  arrayView2d< real64 const > const & deltaTraction =
    fractureSubRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaTractionString() );

  // Get the volume for all elements
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const elemVolume =
    elemManager.constructViewAccessor< real64_array, arrayView1d< real64 const > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  // Get the coordinates for all nodes
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // Get area and rotation matrix for all faces
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();
  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView3d< real64 const > const &
  faceRotationMatrix = fractureSubRegion.getReference< array3d< real64 > >( viewKeyStruct::rotationMatrixString() );

  // Bulk modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const bulkModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::bulkModulusString(),
                                                                                                 m_solidSolver->targetRegionNames(),
                                                                                                 m_solidSolver->solidMaterialNames() );
  // Shear modulus accessor
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const shearModulus =
    elemManager.constructMaterialViewAccessor< array1d< real64 >, arrayView1d< real64 const > >( ElasticIsotropic::viewKeyStruct::shearModulusString(),
                                                                                                 m_solidSolver->targetRegionNames(),
                                                                                                 m_solidSolver->solidMaterialNames() );

  using NodeMapViewType = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  ElementRegionManager::ElementViewAccessor< NodeMapViewType > const elemToNode =
    elemManager.constructViewAccessor< CellElementSubRegion::NodeMapType, NodeMapViewType >( ElementSubRegionBase::viewKeyStruct::nodeListString() );
  ElementRegionManager::ElementViewConst< NodeMapViewType > const elemToNodeView = elemToNode.toNestedViewConst();

  arrayView1d< globalIndex const > const & tracDofNumber = fractureSubRegion.getReference< globalIndex_array >( tracDofKey );

  stabilizationMethod.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

    forAll< serialPolicy >( stencil.size(), [=] ( localIndex const iconn )
    {
      localIndex const numFluxElems = sei.sizeOfArray( iconn );

      // A fracture connector has to be an edge shared by two faces
      if( numFluxElems == 2 )
      {
        // Find shared edge (pair of nodes)
        array1d< real64 > Nbar0( 3 ), Nbar1( 3 );
        Nbar0[ 0 ] = faceRotationMatrix[ sei[iconn][0] ][0][0];
        Nbar0[ 1 ] = faceRotationMatrix[ sei[iconn][0] ][1][0];
        Nbar0[ 2 ] = faceRotationMatrix[ sei[iconn][0] ][2][0];
        Nbar1[ 0 ] = faceRotationMatrix[ sei[iconn][1] ][0][0];
        Nbar1[ 1 ] = faceRotationMatrix[ sei[iconn][1] ][1][0];
        Nbar1[ 2 ] = faceRotationMatrix[ sei[iconn][1] ][2][0];

        real64 normalProduct = LvArray::tensorOps::AiBi< 3 >( Nbar0, Nbar1 );

        localIndex const id1 = ( normalProduct > 0.0 ) ? 0 : 1;

        localIndex const numNodesPerFace0 = faceToNodeMap.sizeOfArray( faceMap[sei[iconn][0]][0] );
        array1d< localIndex > nodes0( numNodesPerFace0 );
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          nodes0[i] = faceToNodeMap( faceMap[sei[iconn][0]][0], i );
        }
        localIndex const numNodesPerFace1 = faceToNodeMap.sizeOfArray( faceMap[sei[iconn][1]][0] );
        array1d< localIndex > nodes1( numNodesPerFace1 );
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          nodes1[i] = faceToNodeMap( faceMap[sei[iconn][1]][id1], i );
        }
        std::sort( nodes0.begin(), nodes0.end() );
        std::sort( nodes1.begin(), nodes1.end() );
        array1d< localIndex > edge( std::max( numNodesPerFace0, numNodesPerFace1 ) );
        edge.setValues< serialPolicy >( -1 );
        std::set_intersection( nodes0.begin(), nodes0.end(), nodes1.begin(), nodes1.end(), edge.begin() );
        localIndex realNodes = 0;
        for( localIndex i = 0; i < edge.size(); ++i )
        {
          if( edge[i] > -1 )
          {
            realNodes++;
          }
        }
        GEOSX_ERROR_IF( realNodes != 2, "An edge shared by two fracture elements must have 2 nodes." );
        edge.resize( realNodes );

        // Compute nodal area factor
        localIndex node0index0 = -1;
        localIndex node1index0 = -1;
        for( localIndex i = 0; i < numNodesPerFace0; ++i )
        {
          if( edge[0] == faceToNodeMap( faceMap[sei[iconn][0]][0], i ) )
          {
            node0index0 = i;
          }
          if( edge[1] == faceToNodeMap( faceMap[sei[iconn][0]][0], i ) )
          {
            node1index0 = i;
          }
        }
        localIndex node0index1 = -1;
        localIndex node1index1 = -1;
        for( localIndex i = 0; i < numNodesPerFace1; ++i )
        {
          if( edge[0] == faceToNodeMap( faceMap[sei[iconn][1]][id1], i ) )
          {
            node0index1 = i;
          }
          if( edge[1] == faceToNodeMap( faceMap[sei[iconn][1]][id1], i ) )
          {
            node1index1 = i;
          }
        }
        array1d< real64 > nodalArea0, nodalArea1;
        computeFaceNodalArea( nodePosition, faceToNodeMap, faceMap[sei[iconn][0]][0], nodalArea0 );
        computeFaceNodalArea( nodePosition, faceToNodeMap, faceMap[sei[iconn][1]][id1], nodalArea1 );
        real64 const areafac = nodalArea0[node0index0] * nodalArea1[node0index1] + nodalArea0[node1index0] * nodalArea1[node1index1];

        // first index: face, second index: element (T/B), third index: dof (x, y, z)
        real64 stiffDiagApprox[ 2 ][ 2 ][ 3 ];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          // Get fracture, face and region/subregion/element indices (for elements on both sides)
          localIndex const fractureIndex = sei[iconn][kf];

          for( localIndex i = 0; i < 2; ++i )
          {
            localIndex const faceIndex = ( kf == 0 || id1 == 0 ) ? faceMap[fractureIndex][i] : faceMap[fractureIndex][1-i];
            localIndex const ke = faceToElemIndex[faceIndex][0] >= 0 ? 0 : 1;

            localIndex const er  = faceToElemRegion[faceIndex][ke];
            localIndex const esr = faceToElemSubRegion[faceIndex][ke];
            localIndex const ei  = faceToElemIndex[faceIndex][ke];

            real64 const volume = elemVolume[er][esr][ei];

            // Get the "element to node" map for the specific region/subregion
            NodeMapViewType const & cellElemsToNodes = elemToNodeView[er][esr];
            localIndex const numNodesPerElem = cellElemsToNodes.size( 1 );

            // Compute the box size
            real64 maxSize[3];
            real64 minSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              maxSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
              minSize[j] = nodePosition[cellElemsToNodes[ei][0]][j];
            }
            for( localIndex a = 1; a < numNodesPerElem; ++a )
            {
              for( localIndex j = 0; j < 3; ++j )
              {
                maxSize[j] = fmax( maxSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
                minSize[j] = fmin( minSize[j], nodePosition[cellElemsToNodes[ei][a]][j] );
              }
            }
            real64 boxSize[3];
            for( localIndex j = 0; j < 3; ++j )
            {
              boxSize[j] = maxSize[j] - minSize[j];
            }

            // Get linear elastic isotropic constitutive parameters for the element
            real64 const K = bulkModulus[er][esr][ei];
            real64 const G = shearModulus[er][esr][ei];
            real64 const E = 9.0 * K * G / ( 3.0 * K + G );
            real64 const nu = ( 3.0 * K - 2.0 * G ) / ( 2.0 * ( 3.0 * K + G ) );

            // Combine E and nu to obtain a stiffness approximation (like it was an hexahedron)
            for( localIndex j = 0; j < 3; ++j )
            {
              stiffDiagApprox[ kf ][ i ][ j ] = E / ( ( 1.0 + nu )*( 1.0 - 2.0*nu ) ) * 2.0 / 9.0 * ( 2.0 - 3.0 * nu ) * volume / ( boxSize[j]*boxSize[j] );
            }
          }
        }
        real64 invTotStiffApprox[ 3 ][ 3 ] = { { 0 } };
        for( localIndex i = 0; i < 3; ++i )
        {
          // K(i,i)^-1 = Ka(i,i)^-1 + Kb(i,i)^-1
          // T -> top (index 0), B -> bottom (index 1)
          // Ka(i,i) = KT(i,i) + KB(i,i)
          // Kb(i,i) = KT(i,i) + KB(i,i)
          invTotStiffApprox[ i ][ i ] = 1.0 / ( stiffDiagApprox[ 0 ][ 0 ][ i ] + stiffDiagApprox[ 1 ][ 0 ][ i ] )
                                        + 1.0 / ( stiffDiagApprox[ 0 ][ 1 ][ i ] + stiffDiagApprox[ 1 ][ 1 ][ i ] );
        }

        array2d< real64 > avgRotationMatrix( 3, 3 );

        // To be able to compute an average rotation matrix, normal has to point in the same direction.
        if( normalProduct < 0.0 )
        {
          LvArray::tensorOps::scale< 3 >( Nbar1, -1.0 );
          normalProduct *= -1.0;
        }
        // If the surfaces are co-planar, then use the first rotation matrix
        if( std::abs( normalProduct - 1.0 ) < 1.e+2*machinePrecision )
        {
          LvArray::tensorOps::copy< 3, 3 >( avgRotationMatrix, faceRotationMatrix[ sei[iconn][0] ] );
        }
        // otherwise, compute the average rotation matrix
        else
        {
          array1d< real64 > avgNbar( 3 );
          avgNbar[ 0 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[0] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[0];
          avgNbar[ 1 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[1] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[1];
          avgNbar[ 2 ] = faceArea[faceMap[ sei[iconn][0] ][0]] * Nbar0[2] + faceArea[faceMap[ sei[iconn][1] ][0]] * Nbar1[2];
          LvArray::tensorOps::normalize< 3 >( avgNbar );

          computationalGeometry::RotationMatrix_3D( avgNbar.toSliceConst(), avgRotationMatrix );
        }

        // Compute R^T * (invK) * R
        real64 temp[ 3 ][ 3 ];
        real64 rotatedInvStiffApprox[ 3 ][ 3 ];
        LvArray::tensorOps::Rij_eq_AkiBkj< 3, 3, 3 >( temp, avgRotationMatrix, invTotStiffApprox );
        LvArray::tensorOps::Rij_eq_AikBkj< 3, 3, 3 >( rotatedInvStiffApprox, temp, avgRotationMatrix );

        // Add nodal area contribution
        stackArray2d< real64, 3*3 > totalInvStiffApprox( 3, 3 );
        for( localIndex i = 0; i < 3; ++i )
        {
          for( localIndex j = 0; j < 3; ++j )
          {
            totalInvStiffApprox( i, j ) = -rotatedInvStiffApprox[ i ][ j ] * areafac;
          }
        }

        // Get DOF numbering
        localIndex fractureIndex[2];
        localIndex nDof[2];
        globalIndex elemDOF[2][3];
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          fractureIndex[kf] = sei[iconn][kf];
          for( localIndex i = 0; i < 3; ++i )
          {
            elemDOF[kf][i] = tracDofNumber[fractureIndex[kf]] + i;
          }
          nDof[kf] = 0;
          switch( fractureState[fractureIndex[kf]] )
          {
            case ( FractureState::STICK ):
              {
                nDof[kf] = 3;
                break;
              }
            case ( FractureState::NEW_SLIP ):
            case ( FractureState::SLIP ):
              {
                nDof[kf] = 1;
                break;
              }
            case ( FractureState::OPEN ):
              {
                nDof[kf] = 0;
                break;
              }
          }
        }

        // Define local "transmissibility" matrices
        stackArray2d< real64, 3*3 > totalInvStiffApprox00( nDof[0], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox01( nDof[0], nDof[1] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox10( nDof[1], nDof[0] );
        stackArray2d< real64, 3*3 > totalInvStiffApprox11( nDof[1], nDof[1] );
        for( localIndex i = 0; i < nDof[0]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox00( i, j ) = totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox01( i, j ) = -totalInvStiffApprox( i, j );
          }
        }

        for( localIndex i = 0; i < nDof[1]; ++i )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            totalInvStiffApprox10( i, j ) = -totalInvStiffApprox( i, j );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            totalInvStiffApprox11( i, j ) = totalInvStiffApprox( i, j );
          }
        }

        // Compute rhs
        stackArray1d< real64, 3 > rhs0( 3 );
        if( nDof[0] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox00( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs0( 0 ) += totalInvStiffApprox01( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[0]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox00( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs0( i ) += totalInvStiffApprox01( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        stackArray1d< real64, 3 > rhs1( 3 );
        if( nDof[1] > 0 )
        {
          for( localIndex j = 0; j < nDof[0]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox10( 0, j ) * ( traction[fractureIndex[0]][j] );
          }
          for( localIndex j = 0; j < nDof[1]; ++j )
          {
            rhs1( 0 ) += totalInvStiffApprox11( 0, j ) * ( traction[fractureIndex[1]][j] );
          }
          for( localIndex i = 1; i < nDof[1]; ++i )
          {
            for( localIndex j = 0; j < nDof[0]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox10( i, j ) * ( deltaTraction[fractureIndex[0]][j] );
            }
            for( localIndex j = 0; j < nDof[1]; ++j )
            {
              rhs1( i ) += totalInvStiffApprox11( i, j ) * ( deltaTraction[fractureIndex[1]][j] );
            }
          }
        }

        // Global matrix and rhs assembly
        for( localIndex kf = 0; kf < 2; ++kf )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( elemDOF[kf][0] - rankOffset );

          stackArray2d< real64, 3*3 > const & totalInvStiffApproxDiag = ( kf == 0 ) ? totalInvStiffApprox00 : totalInvStiffApprox11;
          stackArray2d< real64, 3*3 > const & totalInvStiffApproxOffDiag = ( kf == 0 ) ? totalInvStiffApprox01 : totalInvStiffApprox10;
          stackArray1d< real64, 3 > const & rhs = ( kf == 0 ) ? rhs0 : rhs1;

          // Only assemble contribution if "row" fracture element is local
          // TODO: use parallel atomics
          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            for( localIndex idof = 0; idof < nDof[kf]; ++idof )
            {
              // (i,i)-block
              localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                              elemDOF[kf],
                                                                              totalInvStiffApproxDiag[idof].dataIfContiguous(),
                                                                              nDof[kf] );
              // (i,j)-block
              if( nDof[1-kf] > 0 )
              {
                localMatrix.addToRowBinarySearchUnsorted< parallelHostAtomic >( localRow + idof,
                                                                                elemDOF[1 - kf],
                                                                                totalInvStiffApproxOffDiag[idof].dataIfContiguous(),
                                                                                nDof[1 - kf] );
              }

              // residual
              RAJA::atomicAdd( parallelHostAtomic{}, &localRhs[localRow + idof], rhs[idof] );
            }
          }
        }
      }
    } );
  } );
}

void LagrangianContactSolver::applySystemSolution( DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localSolution,
                                                   real64 const scalingFactor,
                                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );

  dofManager.addVectorToField( localSolution, viewKeyStruct::tractionString(), viewKeyStruct::deltaTractionString(), -scalingFactor );
  dofManager.addVectorToField( localSolution, viewKeyStruct::tractionString(), viewKeyStruct::tractionString(), -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::tractionString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaTractionString() ) );
  // This is used locally only, synchronized just for output reasons
  fieldNames["elems"].emplace_back( string( viewKeyStruct::localJumpString() ) );
  // fractureStateString is synchronized in UpdateFractureState
  // previousFractureStateString and previousLocalJumpString used locally only

  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );

  computeFaceDisplacementJump( domain );
}

void LagrangianContactSolver::initializeFractureState( MeshLevel & mesh,
                                                       string const & fieldName ) const
{
  GEOSX_MARK_FUNCTION;
  ElementRegionManager & elemManager = mesh.getElemManager();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( fieldName );
      fractureState.setValues< parallelHostPolicy >( FractureState::STICK );
    }
  } );
}

void LagrangianContactSolver::setFractureStateForElasticStep( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( fractureState[kfe] != FractureState::OPEN )
        {
          fractureState[kfe] = FractureState::STICK;
        }
      } );
    }
  } );
}

bool LagrangianContactSolver::updateFractureState( DomainPartition & domain ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();
  ContactRelationBase const & contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_contactRelationName );

  bool checkActiveSet = true;

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >( viewKeyStruct::tractionString() );
      arrayView2d< real64 const > const & localJump = subRegion.getReference< array2d< real64 > >( viewKeyStruct::localJumpString() );
      arrayView1d< integer > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );

      arrayView1d< real64 const > const & normalTractionTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalTractionToleranceString() );
      arrayView1d< real64 const > const & normalDisplacementTolerance =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::normalDisplacementToleranceString() );

      RAJA::ReduceMin< parallelHostReduce, integer > checkActiveSetSub( 1 );

      forAll< parallelHostPolicy >( subRegion.size(), [=, &contactRelation] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          integer const originalFractureState = fractureState[kfe];
          if( originalFractureState == FractureState::OPEN )
          {
            if( localJump[kfe][0] > -normalDisplacementTolerance[kfe] )
            {
              fractureState[kfe] = FractureState::OPEN;
            }
            else
            {
              fractureState[kfe] = FractureState::STICK;
            }
          }
          else if( traction[kfe][0] > normalTractionTolerance[kfe] )
          {
            fractureState[kfe] = FractureState::OPEN;
          }
          else
          {
            real64 currentTau = sqrt( traction[kfe][1]*traction[kfe][1] + traction[kfe][2]*traction[kfe][2] );
            real64 const limitTau = contactRelation.limitTangentialTractionNorm( traction[kfe][0] );
            if( originalFractureState == FractureState::STICK && currentTau >= limitTau )
            {
              currentTau *= (1.0 - m_slidingCheckTolerance);
            }
            else if( originalFractureState != FractureState::STICK && currentTau <= limitTau )
            {
              currentTau *= (1.0 + m_slidingCheckTolerance);
            }
            if( currentTau > limitTau )
            {
              if( originalFractureState == FractureState::STICK )
              {
                fractureState[kfe] = FractureState::NEW_SLIP;
              }
              else
              {
                fractureState[kfe] = FractureState::SLIP;
              }
            }
            else
            {
              fractureState[kfe] = FractureState::STICK;
            }
          }

          if( originalFractureState != fractureState[kfe] )
          {
//            GEOSX_LOG_LEVEL_BY_RANK( 3, "element " << kfe << " traction: " << traction[kfe]
//                                                   << " previous state <"
//                                                   << FractureStateToString( originalFractureState )
//                                                   << "> current state <"
//                                                   << FractureStateToString( fractureState[kfe] )
//                                                   << ">" );
          }
          checkActiveSetSub.min( compareFractureStates( originalFractureState, fractureState[kfe] ) );
        }
      } );

      checkActiveSet &= checkActiveSetSub.get();
    }
  } );

  // Need to synchronize the fracture state due to the use will be made of in AssemblyStabilization
  synchronizeFractureState( domain );

  // Compute if globally the fracture state has changed
  bool globalCheckActiveSet;
  MpiWrapper::allReduce( &checkActiveSet,
                         &globalCheckActiveSet,
                         1,
                         MPI_LAND,
                         MPI_COMM_GEOSX );

  return globalCheckActiveSet;
}

void LagrangianContactSolver::synchronizeFractureState( DomainPartition & domain ) const
{
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::fractureStateString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );
}

bool LagrangianContactSolver::isFractureAllInStickCondition( DomainPartition const & domain ) const
{
  globalIndex numStick, numSlip, numOpen;
  computeFractureStateStatistics( domain, numStick, numSlip, numOpen, false );
  return ( ( numSlip + numOpen ) == 0 );
}

void LagrangianContactSolver::computeFractureStateStatistics( DomainPartition const & domain,
                                                              globalIndex & numStick,
                                                              globalIndex & numSlip,
                                                              globalIndex & numOpen,
                                                              bool printAll ) const
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  array1d< localIndex > localCounter( 3 );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {
    if( subRegion.hasWrapper( m_tractionKey ) )
    {
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
      arrayView1d< integer const > const & fractureState = subRegion.getReference< array1d< integer > >( viewKeyStruct::fractureStateString() );
//      arrayView2d< real64 const > const & traction = subRegion.getReference< array2d< real64 > >(
// viewKeyStruct::tractionString );

      RAJA::ReduceSum< parallelHostReduce, localIndex > stickCount( 0 ), slipCount( 0 ), openCount( 0 );
      forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        if( ghostRank[kfe] < 0 )
        {
          switch( fractureState[kfe] )
          {
            case FractureState::STICK:
              {
                stickCount += 1;
                break;
              }
            case FractureState::NEW_SLIP:
            case FractureState::SLIP:
              {
                slipCount += 1;
                break;
              }
            case FractureState::OPEN:
              {
                openCount += 1;
                break;
              }
          }
          if( printAll )
          {
//            GEOSX_LOG_LEVEL_BY_RANK( 3, "element " << kfe << " traction: " << traction[kfe]
//                                                   << " state <"
//                                                   << FractureStateToString( fractureState[kfe] )
//                                                   << ">" );
          }
        }
      } );

      localCounter[0] += stickCount.get();
      localCounter[1] += slipCount.get();
      localCounter[2] += openCount.get();
    }
  } );

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );

  array1d< globalIndex > globalCounter( 3*size );

  // Everything is done on rank 0
  MpiWrapper::gather( localCounter.data(),
                      3,
                      globalCounter.data(),
                      3,
                      0,
                      MPI_COMM_GEOSX );

  array1d< globalIndex > totalCounter( 3 );

  if( rank == 0 )
  {
    for( int r = 0; r < size; ++r )
    {
      // sum across all ranks
      totalCounter[0] += globalCounter[3*r];
      totalCounter[1] += globalCounter[3*r+1];
      totalCounter[2] += globalCounter[3*r+2];
    }
  }

  MpiWrapper::bcast( totalCounter.data(), 3, 0, MPI_COMM_GEOSX );

  numStick = totalCounter[0];
  numSlip  = totalCounter[1];
  numOpen  = totalCounter[2];

  char output[108] = {0};
  sprintf( output,
           " Number of element for each fracture state:"
#if defined(GEOSX_USE_HYPRE_CUDA) && defined(GEOSX_LA_INTERFACE_HYPRE)
           " stick: %12i | slip:  %12i | open:  %12i",
#else
           " stick: %12lli | slip:  %12lli | open:  %12lli",
#endif
           numStick,
           numSlip,
           numOpen );
  GEOSX_LOG_RANK_0( output );
}

void LagrangianContactSolver::solveSystem( DofManager const & dofManager,
                                           ParallelMatrix & matrix,
                                           ParallelVector & rhs,
                                           ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  if( getLogLevel() > 3 )
  {
    matrix.write( "matrix.mtx", LAIOutputFormat::MATRIX_MARKET );
    rhs.write( "rhs.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() > 3 )
  {
    solution.write( "sol.mtx", LAIOutputFormat::MATRIX_MARKET );
  }

  // int rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  // if( rank == 0 )
  // {
  //   string str;
  //   std::getline( std::cin, str );
  //   if( str.length() > 0 )
  //   {
  //     GEOSX_ERROR( "STOP" );
  //   }
  // }
  // MpiWrapper::Barrier( MPI_COMM_GEOSX );
}

void LagrangianContactSolver::setNextDt( real64 const & currentDt,
                                         real64 & nextDt )
{
  nextDt = currentDt;
}

REGISTER_CATALOG_ENTRY( SolverBase, LagrangianContactSolver, string const &, Group * const )
} /* namespace geosx */
