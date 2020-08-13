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
 * @file SolidMechanicsLagrangianFEM.cpp
 */

#include "SolidMechanicsLagrangianFEM.hpp"
#include "SolidMechanicsPoroElasticKernel.hpp"
#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"
#include "SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp"
#include "SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"
#include "SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

SolidMechanicsLagrangianFEM::SolidMechanicsLagrangianFEM( const std::string & name,
                                                          Group * const parent ) :
  SolverBase( name, parent ),
  m_newmarkGamma( 0.5 ),
  m_newmarkBeta( 0.25 ),
  m_massDamping( 0.0 ),
  m_stiffnessDamping( 0.0 ),
  m_timeIntegrationOption( TimeIntegrationOption::ExplicitDynamic ),
  m_useVelocityEstimateForQS( 0 ),
  m_maxForce( 0.0 ),
  m_maxNumResolves( 10 ),
  m_strainTheory( 0 ),
  //  m_elemsAttachedToSendOrReceiveNodes(),
  //  m_elemsNotAttachedToSendOrReceiveNodes(),
  m_sendOrReceiveNodes(),
  m_nonSendOrReceiveNodes(),
  m_iComm(),
  m_effectiveStress( 0 )
{
  m_sendOrReceiveNodes.setName(
    "SolidMechanicsLagrangianFEM::m_sendOrReceiveNodes" );
  m_nonSendOrReceiveNodes.setName(
    "SolidMechanicsLagrangianFEM::m_nonSendOrReceiveNodes" );

  registerWrapper( viewKeyStruct::newmarkGammaString, &m_newmarkGamma )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic "
      "time integration option" );

  registerWrapper( viewKeyStruct::newmarkBetaString, &m_newmarkBeta )
    ->setApplyDefaultValue( 0.25 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time "
      "integration option. "
      "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you "
      "are doing." );

  registerWrapper( viewKeyStruct::massDampingString, &m_massDamping )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Value of mass based damping coefficient. " );

  registerWrapper( viewKeyStruct::stiffnessDampingString, &m_stiffnessDamping )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Value of stiffness based damping coefficient. " );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString,
                   &m_timeIntegrationOption )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setApplyDefaultValue( m_timeIntegrationOption )
    ->setDescription(
      "Time integration method. Options are: \n QuasiStatic \n ImplicitDynamic "
      "\n ExplicitDynamic" );

  registerWrapper( viewKeyStruct::useVelocityEstimateForQSString,
                   &m_useVelocityEstimateForQS )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Flag to indicate the use of the incremental displacement from the "
      "previous step as an "
      "initial estimate for the incremental displacement of the current step." );

  registerWrapper( viewKeyStruct::maxNumResolvesString, &m_maxNumResolves )
    ->setApplyDefaultValue( 10 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Value to indicate how many resolves may be executed after some other "
      "event is executed. "
      "For example, if a SurfaceGenerator is specified, it will be executed "
      "after the mechanics solve. "
      "However if a new surface is generated, then the mechanics solve must be "
      "executed again due to the "
      "change in topology." );

  registerWrapper( viewKeyStruct::strainTheoryString, &m_strainTheory )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Indicates whether or not to use "
      "`Infinitesimal Strain Theory "
      "<https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
      "`Finite Strain Theory "
      "<https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs "
      "are:\n"
      " 0 - Infinitesimal Strain \n"
      " 1 - Finite Strain" );

  registerWrapper( viewKeyStruct::solidMaterialNamesString, &m_solidMaterialNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "The name of the material that should be used in the constitutive "
      "updates" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )
    ->setApplyDefaultValue( viewKeyStruct::noContactRelationNameString )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxForce, &m_maxForce )
    ->setInputFlag( InputFlags::FALSE )
    ->setDescription( "The maximum force contribution in the problem domain." );

  registerWrapper( viewKeyStruct::effectiveStress, &m_effectiveStress )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Apply fluid pressure to produce effective stress when integrating "
      "stress." );
}

void
SolidMechanicsLagrangianFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  CheckModelNames( m_solidMaterialNames, viewKeyStruct::solidMaterialNamesString );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.isSymmetric = true;
  linParams.dofsPerNode = 3;
  linParams.amg.separateComponents = true;
}

SolidMechanicsLagrangianFEM::~SolidMechanicsLagrangianFEM()
{
  // TODO Auto-generated destructor stub
}

void
SolidMechanicsLagrangianFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes =
      mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getNodeManager();

    nodes
      ->registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >(
        keys::TotalDisplacement )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the total displacements on the nodes." )
      ->reference()
      .resizeDimension< 1 >( 3 );

    nodes
      ->registerWrapper< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >(
        keys::IncrementalDisplacement )
      ->setPlotLevel( PlotLevel::LEVEL_3 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the incremental displacements for the current "
        "time step on the nodes." )
      ->reference()
      .resizeDimension< 1 >( 3 );

    nodes
      ->registerWrapper< array2d< real64, nodes::VELOCITY_PERM > >( keys::Velocity )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription( "An array that holds the current velocity on the nodes." )
      ->reference()
      .resizeDimension< 1 >( 3 );

    nodes
      ->registerWrapper< array2d< real64, nodes::ACCELERATION_PERM > >(
        keys::Acceleration )
      ->setPlotLevel( PlotLevel::LEVEL_1 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the current acceleration on the nodes. This array "
        "also is used "
        "to hold the summation of nodal forces resulting from the governing "
        "equations." )
      ->reference()
      .resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array2d< real64 > >( viewKeyStruct::forceExternal )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the external forces on the nodes. This includes "
        "any boundary"
        " conditions as well as coupling forces such as hydraulic forces." )
      ->reference()
      .resizeDimension< 1 >( 3 );

    nodes->registerWrapper< array1d< real64 > >( keys::Mass )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription( "An array that holds the mass on the nodes." );

    nodes->registerWrapper< array1d< R1Tensor > >( viewKeyStruct::vTildeString )
      ->setPlotLevel( PlotLevel::NOPLOT )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the velocity predictors on the nodes." );

    nodes->registerWrapper< array1d< R1Tensor > >( viewKeyStruct::uhatTildeString )
      ->setPlotLevel( PlotLevel::NOPLOT )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the incremental displacement predictors on the "
        "nodes." );

    nodes->registerWrapper< array1d< R1Tensor > >( viewKeyStruct::contactForceString )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription( "An array that holds the contact force." );

    ElementRegionManager * const elementRegionManager =
      mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getElemManager();
    elementRegionManager->forElementSubRegions< CellElementSubRegion >(
      [&]( CellElementSubRegion & subRegion ) {
        subRegion
          .registerWrapper< array3d< real64, solid::STRESS_PERMUTATION > >(
            viewKeyStruct::stress_n )
          ->setPlotLevel( PlotLevel::NOPLOT )
          ->setRestartFlags( RestartFlags::NO_WRITE )
          ->setRegisteringObjects( this->getName() )
          ->setDescription(
            "Array to hold the beginning of step stress for implicit problem "
            "rewinds" )
          ->reference()
          .resizeDimension< 2 >( 6 );

        subRegion
          .registerWrapper< SortedArray< localIndex > >(
            viewKeyStruct::elemsAttachedToSendOrReceiveNodes )
          ->setPlotLevel( PlotLevel::NOPLOT )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< SortedArray< localIndex > >(
            viewKeyStruct::elemsNotAttachedToSendOrReceiveNodes )
          ->setPlotLevel( PlotLevel::NOPLOT )
          ->setRestartFlags( RestartFlags::NO_WRITE );
      } );
  }
}

void
SolidMechanicsLagrangianFEM::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );

  // Validate solid models in target regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel =
      *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< SolidBase >( *meshLevel.getElemManager(),
                                       m_solidMaterialNames );
  }

  NumericalMethodsManager const & numericalMethodManager =
    domain->getNumericalMethodManager();

  FiniteElementDiscretizationManager const & feDiscretizationManager =
    numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const * feDiscretization =
    feDiscretizationManager.GetGroup< FiniteElementDiscretization >(
      m_discretizationName );
  GEOSX_ERROR_IF(
    feDiscretization == nullptr,
    getName() << ": FE discretization not found: " << m_discretizationName );
}

template< typename... PARAMS >
real64
SolidMechanicsLagrangianFEM::explicitKernelDispatch( PARAMS &&... params )
{
  real64 rval = 0;
  if( m_strainTheory == 0 )
  {
    rval = finiteElement::regionBasedKernelApplication<
      parallelDevicePolicy< 32 >,
      constitutive::SolidBase,
      CellElementSubRegion,
      SolidMechanicsLagrangianFEMKernels::ExplicitSmallStrain >(
      std::forward< PARAMS >( params )... );
  }
  else if( m_strainTheory == 1 )
  {
    rval = finiteElement::regionBasedKernelApplication<
      parallelDevicePolicy< 32 >,
      constitutive::SolidBase,
      CellElementSubRegion,
      SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrain >(
      std::forward< PARAMS >( params )... );
  }
  else
  {
    GEOSX_ERROR(
      "Invalid option for strain theory (0 = infinitesimal strain, 1 = finite "
      "strain" );
  }
  return rval;
}

void
SolidMechanicsLagrangianFEM::updateIntrinsicNodalData(
  DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodes = *mesh.getNodeManager();
  //FaceManager * const faceManager = mesh->getFaceManager();

  ElementRegionManager const & elementRegionManager = *mesh.getElemManager();

  arrayView1d< real64 > & mass = nodes.getReference< array1d< real64 > >( keys::Mass );
  mass.setValues< serialPolicy >( 0.0 );

  arrayView1d< integer const > const & nodeGhostRank = nodes.ghostRank();

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > rho =
    elementRegionManager
      .ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 const > >(
        "density",
        targetRegionNames(),
        solidMaterialNames() );

  forTargetRegionsComplete(
    mesh,
    [&]( localIndex const, localIndex const er, ElementRegionBase const & elemRegion ) {
      elemRegion.forElementSubRegionsIndex< CellElementSubRegion >(
        [&]( localIndex const esr, CellElementSubRegion const & elementSubRegion ) {
          arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
          arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes =
            elementSubRegion.nodeList();

          finiteElement::FiniteElementBase const & fe =
            elementSubRegion.getReference< finiteElement::FiniteElementBase >(
              this->getDiscretizationName() );
          finiteElement::dispatch3D( fe, [&]( auto const finiteElement ) {
            using FE_TYPE = TYPEOFREF( finiteElement );

            constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
            constexpr localIndex numQuadraturePointsPerElem =
              FE_TYPE::numQuadraturePoints;

            real64 N[numNodesPerElem];
            for( localIndex k = 0; k < elemsToNodes.size( 0 ); ++k )
            {
              real64 elemMass = 0;
              for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
              {
                elemMass += rho[er][esr][k][q] * detJ[k][q];
                FE_TYPE::shapeFunctionValues( q, N );

                for( localIndex a = 0; a < numNodesPerElem; ++a )
                {
                  mass[elemsToNodes[k][a]] +=
                    rho[er][esr][k][q] * detJ[k][q] * N[a];
                }
              }

              for( localIndex a = 0; a < elementSubRegion.numNodesPerElement();
                   ++a )
              {
                if( nodeGhostRank[elemsToNodes[k][a]] >= -1 )
                {
                  m_sendOrReceiveNodes.insert( elemsToNodes[k][a] );
                }
                else
                {
                  m_nonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
                }
              }
            }
          } );
        } );
    } );
}

void
SolidMechanicsLagrangianFEM::InitializePostInitialConditions_PreSubGroups(
  Group * const problemManager )
{
  DomainPartition * domain =
    problemManager->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  NodeManager & nodes = *mesh.getNodeManager();
  //FaceManager * const faceManager = mesh.getFaceManager();

  ElementRegionManager & elementRegionManager = *mesh.getElemManager();

  arrayView1d< real64 > & mass = nodes.getReference< array1d< real64 > >( keys::Mass );

  arrayView1d< integer const > const & nodeGhostRank = nodes.ghostRank();

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > rho =
    elementRegionManager
      .ConstructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 const > >(
        "density",
        targetRegionNames(),
        solidMaterialNames() );

  forTargetRegionsComplete( mesh, [&]( localIndex const, localIndex const er, ElementRegionBase & elemRegion ) {
    elemRegion.forElementSubRegionsIndex<
      CellElementSubRegion >( [&]( localIndex const esr,
                                   CellElementSubRegion & elementSubRegion ) {
      SortedArray< localIndex > & elemsAttachedToSendOrReceiveNodes =
        getElemsAttachedToSendOrReceiveNodes( elementSubRegion );
      SortedArray< localIndex > & elemsNotAttachedToSendOrReceiveNodes =
        getElemsNotAttachedToSendOrReceiveNodes( elementSubRegion );

      elemsAttachedToSendOrReceiveNodes.setName(
        "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes[" +
        std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

      elemsNotAttachedToSendOrReceiveNodes.setName(
        "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes[" +
        std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes =
        elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const & fe =
        elementSubRegion.getReference< finiteElement::FiniteElementBase >(
          getDiscretizationName() );
      finiteElement::dispatch3D( fe, [&]( auto const finiteElement ) {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem =
          FE_TYPE::numQuadraturePoints;

        real64 N[numNodesPerElem];
        for( localIndex k = 0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 elemMass = 0;
          for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
          {
            elemMass += rho[er][esr][k][q] * detJ[k][q];
            FE_TYPE::shapeFunctionValues( q, N );

            for( localIndex a = 0; a < numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[er][esr][k][q] * detJ[k][q] * N[a];
            }
          }

          bool isAttachedToGhostNode = false;
          for( localIndex a = 0; a < elementSubRegion.numNodesPerElement(); ++a )
          {
            if( nodeGhostRank[elemsToNodes[k][a]] >= -1 )
            {
              isAttachedToGhostNode = true;
              m_sendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
            else
            {
              m_nonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
          }

          if( isAttachedToGhostNode )
          {
            elemsAttachedToSendOrReceiveNodes.insert( k );
          }
          else
          {
            elemsNotAttachedToSendOrReceiveNodes.insert( k );
          }
        }
      } );
    } );
  } );
}

real64
SolidMechanicsLagrangianFEM::SolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         const int cycleNumber,
                                         DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =
    this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );

    if( surfaceGenerator != nullptr )
    {
      surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain );
    }
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    int const maxNumResolves = m_maxNumResolves;
    int locallyFractured = 0;
    int globallyFractured = 0;
    ImplicitStepSetup( time_n, dt, domain );
    for( int solveIter = 0; solveIter < maxNumResolves; ++solveIter )
    {
      SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

      if( solveIter > 0 )
      {
        ResetStressToBeginningOfStep( domain );
      }

      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

      //      updateStress( domain );
      if( surfaceGenerator != nullptr )
      {
        if( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 )
        {
          locallyFractured = 1;
        }
        MpiWrapper::allReduce( &locallyFractured,
                               &globallyFractured,
                               1,
                               MPI_MAX,
                               MPI_COMM_GEOSX );
      }
      if( globallyFractured == 0 )
      {
        break;
      }
      else
      {
        GEOSX_LOG_RANK_0( "Fracture Occurred. Resolve" );
      }
    }
    ImplicitStepComplete( time_n, dt, domain );
  }

  return dtReturn;
}

real64
SolidMechanicsLagrangianFEM::ExplicitStep(
  real64 const & time_n,
  real64 const & dt,
  const int GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

#define USE_PHYSICS_LOOP

  // updateIntrinsicNodalData(domain);

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodes = *mesh.getNodeManager();

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  arrayView1d< real64 const > const & mass =
    nodes.getReference< array1d< real64 > >( keys::Mass );
  arrayView2d< real64, nodes::VELOCITY_USD > const & vel = nodes.velocity();

  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u =
    nodes.totalDisplacement();
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat =
    nodes.incrementalDisplacement();
  arrayView2d< real64, nodes::ACCELERATION_USD > const & acc = nodes.acceleration();

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::Velocity );
  fieldNames["node"].emplace_back( keys::Acceleration );

  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames,
                                                    &mesh,
                                                    domain.getNeighbors(),
                                                    m_iComm,
                                                    true );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n,
                                                             &domain,
                                                             "nodeManager",
                                                             keys::Acceleration );

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt / 2 );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n,
                                                             &domain,
                                                             "nodeManager",
                                                             keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

  fsManager.ApplyFieldValue(
    time_n + dt,
    &domain,
    "nodeManager",
    NodeManager::viewKeyStruct::totalDisplacementString,
    [&]( FieldSpecificationBase const * const bc,
         SortedArrayView< localIndex const > const & targetSet ) {
      integer const component = bc->GetComponent();
      forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
                                              [=] GEOSX_DEVICE( localIndex const i ) {
                                                localIndex const a = targetSet[i];
                                                vel( a, component ) = u( a, component );
                                              } );
    },
    [&]( FieldSpecificationBase const * const bc,
         SortedArrayView< localIndex const > const & targetSet ) {
      integer const component = bc->GetComponent();
      forAll< parallelDevicePolicy< 1024 > >(
        targetSet.size(),
        [=] GEOSX_DEVICE( localIndex const i ) {
          localIndex const a = targetSet[i];
          uhat( a, component ) = u( a, component ) - vel( a, component );
          vel( a, component ) = uhat( a, component ) / dt;
        } );
    } );

  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  explicitKernelDispatch(
    mesh,
    targetRegionNames(),
    this->getDiscretizationName(),
    m_solidMaterialNames,
    dt,
    string( viewKeyStruct::elemsAttachedToSendOrReceiveNodes ) );

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate(
    acc,
    mass,
    vel,
    dt / 2,
    m_sendOrReceiveNodes.toViewConst() );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n,
                                                             &domain,
                                                             "nodeManager",
                                                             keys::Velocity );

  CommunicationTools::SynchronizePackSendRecv( fieldNames,
                                               &mesh,
                                               domain.getNeighbors(),
                                               m_iComm,
                                               true );

  explicitKernelDispatch(
    mesh,
    targetRegionNames(),
    this->getDiscretizationName(),
    m_solidMaterialNames,
    dt,
    string( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodes ) );

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate(
    acc,
    mass,
    vel,
    dt / 2,
    m_nonSendOrReceiveNodes.toViewConst() );

  fsManager.ApplyFieldValue< parallelDevicePolicy< 1024 > >( time_n,
                                                             &domain,
                                                             "nodeManager",
                                                             keys::Velocity );

  CommunicationTools::SynchronizeUnpack( &mesh, domain.getNeighbors(), m_iComm, true );

  return dt;
}

void
SolidMechanicsLagrangianFEM::ApplyDisplacementBC_implicit(
  real64 const time,
  DofManager const & dofManager,
  DomainPartition & domain,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const fieldName ) {
                     bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                                         parallelDevicePolicy< 32 > >(
                       targetSet,
                       time,
                       targetGroup,
                       fieldName,
                       dofKey,
                       dofManager.rankOffset(),
                       localMatrix,
                       localRhs );
                   } );
}

void
SolidMechanicsLagrangianFEM::CRSApplyTractionBC(
  real64 const time,
  DofManager const & dofManager,
  DomainPartition & domain,
  arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  FunctionManager const & functionManager = FunctionManager::Instance();

  FaceManager const & faceManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getFaceManager();
  NodeManager const & nodeManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();

  arrayView1d< real64 const > const & faceArea =
    faceManager.getReference< real64_array >( "faceArea" );
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    faceManager.nodeList().toViewConst();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & blockLocalDofNumber =
    nodeManager.getReference< globalIndex_array >( dofKey );
  globalIndex const dofRankOffset = dofManager.rankOffset();

  fsManager.Apply(
    time,
    &domain,
    "faceManager",
    string( "Traction" ),
    [&]( FieldSpecificationBase const * const bc,
         string const &,
         SortedArrayView< localIndex const > const & targetSet,
         Group * const,
         string const & ) {
      string const & functionName = bc->GetFunctionName();

      globalIndex_array nodeDOF;
      real64_array nodeRHS;
      integer const component = bc->GetComponent();

      if( functionName.empty() ||
          functionManager.getGroupReference< FunctionBase >( functionName )
              .isFunctionOfTime() == 2 )
      {
        real64 value = bc->GetScale();
        if( !functionName.empty() )
        {
          FunctionBase const & function =
            functionManager.getGroupReference< FunctionBase >( functionName );
          value *= function.Evaluate( &time );
        }

        forAll< parallelDevicePolicy< 32 > >(
          targetSet.size(),
          [=] GEOSX_HOST_DEVICE( localIndex const i ) {
            localIndex const kf = targetSet[i];
            localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
            for( localIndex a = 0; a < numNodes; ++a )
            {
              localIndex const dof = blockLocalDofNumber[faceToNodeMap( kf, a )] +
                component - dofRankOffset;
              if( dof < 0 || dof >= localRhs.size() ) continue;
              RAJA::atomicAdd< parallelDeviceAtomic >(
                &localRhs[dof],
                value * faceArea[kf] / numNodes );
            }
          } );
      }
      else
      {
        FunctionBase const & function =
          functionManager.getGroupReference< FunctionBase >( functionName );
        array1d< real64 > resultsArray( targetSet.size() );
        resultsArray.setName(
          "SolidMechanicsLagrangianFEM::TractionBC function results" );
        function.Evaluate( &faceManager, time, targetSet, resultsArray );
        arrayView1d< real64 const > const & results = resultsArray.toView();

        forAll< parallelDevicePolicy< 32 > >(
          targetSet.size(),
          [=] GEOSX_HOST_DEVICE( localIndex const i ) {
            localIndex const kf = targetSet[i];
            localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
            for( localIndex a = 0; a < numNodes; ++a )
            {
              localIndex const dof = blockLocalDofNumber[faceToNodeMap( kf, a )] +
                component - dofRankOffset;
              if( dof < 0 || dof >= localRhs.size() ) continue;
              RAJA::atomicAdd< parallelDeviceAtomic >(
                &localRhs[dof],
                results[kf] * faceArea[kf] / numNodes );
            }
          } );
      }
    } );
}

void
SolidMechanicsLagrangianFEM::ApplyChomboPressure(
  DofManager const & dofManager,
  DomainPartition & domain,
  arrayView1d< real64 > const & localRhs )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *mesh.getFaceManager();
  NodeManager & nodeManager = *mesh.getNodeManager();

  arrayView1d< real64 const > const & faceArea = faceManager.faceArea();
  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    faceManager.nodeList().toViewConst();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dofNumber =
    nodeManager.getReference< globalIndex_array >( dofKey );
  arrayView1d< real64 const > const & facePressure =
    faceManager.getReference< array1d< real64 > >( "ChomboPressure" );

  forAll< serialPolicy >( faceManager.size(), [=]( localIndex const kf ) {
    int const numNodes =
      LvArray::integerConversion< int >( faceToNodeMap.sizeOfArray( kf ) );
    for( int a = 0; a < numNodes; ++a )
    {
      for( int component = 0; component < 3; ++component )
      {
        localIndex const dof = dofNumber[faceToNodeMap( kf, a )] + component;
        real64 const value = -facePressure[kf] * faceNormal( kf, component ) *
          faceArea[kf] / numNodes;
        localRhs[dof] += value;
      }
    }
  } );
}

void
SolidMechanicsLagrangianFEM::ImplicitStepSetup(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();

  arrayView2d< real64 const, nodes::VELOCITY_USD > const & v_n =
    nodeManager.velocity();
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat =
    nodeManager.incrementalDisplacement();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & disp =
    nodeManager.totalDisplacement();

  localIndex const numNodes = nodeManager.size();

  if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d< real64 const, nodes::ACCELERATION_USD > const & a_n =
      nodeManager.acceleration();
    arrayView1d< R1Tensor > const & vtilde =
      nodeManager.getReference< array1d< R1Tensor > >( solidMechanicsViewKeys.vTilde );
    arrayView1d< R1Tensor > const & uhatTilde =
      nodeManager.getReference< array1d< R1Tensor > >(
        solidMechanicsViewKeys.uhatTilde );

    real64 const newmarkGamma =
      this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
    real64 const newmarkBeta =
      this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

    forAll< parallelDevicePolicy< 32 > >(
      numNodes,
      [=] GEOSX_HOST_DEVICE( localIndex const a ) {
        for( int i = 0; i < 3; ++i )
        {
          vtilde[a][i] = v_n( a, i ) + ( 1.0 - newmarkGamma ) * a_n( a, i ) * dt;
          uhatTilde[a][i] =
            ( v_n( a, i ) + 0.5 * ( 1.0 - 2.0 * newmarkBeta ) * a_n( a, i ) * dt ) * dt;
          uhat( a, i ) = uhatTilde[a][i];
          disp( a, i ) += uhatTilde[a][i];
        }
      } );
  }
  else if( this->m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    if( m_useVelocityEstimateForQS == 1 )
    {
      forAll< parallelDevicePolicy< 32 > >( numNodes,
                                            [=] GEOSX_HOST_DEVICE( localIndex const a ) {
                                              for( int i = 0; i < 3; ++i )
                                              {
                                                uhat( a, i ) = v_n( a, i ) * dt;
                                                disp( a, i ) += uhat( a, i );
                                              }
                                            } );
    }
    else
    {
      forAll< parallelDevicePolicy< 32 > >( numNodes,
                                            [=] GEOSX_HOST_DEVICE( localIndex const a ) {
                                              for( int i = 0; i < 3; ++i )
                                              {
                                                uhat( a, i ) = 0.0;
                                              }
                                            } );
    }
  }

  ElementRegionManager * const elementRegionManager = mesh.getElemManager();
  ConstitutiveManager * const constitutiveManager =
    domain.getConstitutiveManager();
  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
    constitutiveRelations =
      elementRegionManager->ConstructFullConstitutiveAccessor< ConstitutiveBase >(
        constitutiveManager );

  forTargetSubRegions< CellElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, CellElementSubRegion & subRegion ) {
      SolidBase const & constitutiveRelation =
        GetConstitutiveModel< SolidBase >( subRegion,
                                           m_solidMaterialNames[targetIndex] );

      arrayView3d< real64 const, solid::STRESS_USD > const & stress =
        constitutiveRelation.getStress();

      array3d< real64, solid::STRESS_PERMUTATION > & stress_n =
        subRegion.getReference< array3d< real64, solid::STRESS_PERMUTATION > >(
          viewKeyStruct::stress_n );
      // TODO: eliminate
      stress_n.resize( stress.size( 0 ), stress.size( 1 ), 6 );

      arrayView3d< real64, solid::STRESS_USD > const & vstress_n = stress_n.toView();

      forAll< parallelDevicePolicy<> >( stress.size( 0 ),
                                        [=] GEOSX_HOST_DEVICE( localIndex const k ) {
                                          for( localIndex a = 0; a < stress.size( 1 );
                                               ++a )
                                          {
                                            for( localIndex i = 0; i < 6; ++i )
                                            {
                                              vstress_n( k, a, i ) = stress( k, a, i );
                                            }
                                          }
                                        } );
    } );
}

void
SolidMechanicsLagrangianFEM::ImplicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();
  localIndex const numNodes = nodeManager.size();

  arrayView2d< real64, nodes::VELOCITY_USD > const & v_n = nodeManager.velocity();
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & uhat =
    nodeManager.incrementalDisplacement();

  if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d< real64, nodes::ACCELERATION_USD > const & a_n =
      nodeManager.acceleration();
    arrayView1d< R1Tensor const > const & vtilde =
      nodeManager.getReference< r1_array >( solidMechanicsViewKeys.vTilde );
    arrayView1d< R1Tensor const > const & uhatTilde =
      nodeManager.getReference< r1_array >( solidMechanicsViewKeys.uhatTilde );
    real64 const newmarkGamma =
      this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
    real64 const newmarkBeta =
      this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

    RAJA::forall< parallelDevicePolicy<> >(
      RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
      [=] GEOSX_HOST_DEVICE( localIndex const a ) {
        for( int i = 0; i < 3; ++i )
        {
          a_n( a, i ) =
            1.0 / ( newmarkBeta * dt * dt ) * ( uhat( a, i ) - uhatTilde[a][i] );
          v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n( a, i ) * dt;
        }
      } );
  }
  else if( this->m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic &&
           dt > 0.0 )
  {
    RAJA::forall< parallelDevicePolicy<> >(
      RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
      [=] GEOSX_HOST_DEVICE( localIndex const a ) {
        for( int i = 0; i < 3; ++i )
        {
          v_n[a][i] = uhat( a, i ) / dt;
        }
      } );
  }
}

void
SolidMechanicsLagrangianFEM::SetupDofs(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  dofManager.addField( keys::TotalDisplacement, DofManager::Location::Node, 3 );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );
}

void
SolidMechanicsLagrangianFEM::SetupSystem(
  DomainPartition & domain,
  DofManager & dofManager,
  CRSMatrix< real64, globalIndex > & localMatrix,
  array1d< real64 > & localRhs,
  array1d< real64 > & localSolution,
  bool const setSparisty )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::SetupSystem( domain,
                           dofManager,
                           localMatrix,
                           localRhs,
                           localSolution,
                           setSparisty );

  MeshLevel & mesh =
    *( domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 ) );
  NodeManager const & nodeManager = *( mesh.getNodeManager() );
  arrayView1d< globalIndex const > const & dofNumber =
    nodeManager.getReference< globalIndex_array >(
      dofManager.getKey( keys::TotalDisplacement ) );

  SparsityPattern< globalIndex > pattern( dofManager.numLocalDofs(),
                                          dofManager.numGlobalDofs() );
  array1d< localIndex > rowSizes( dofManager.numLocalDofs() );

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString )
  {
    ElementRegionManager const & elemManager = *mesh.getElemManager();
    array1d< string > allFaceElementRegions;
    elemManager.forElementRegions< FaceElementRegion >(
      [&]( FaceElementRegion const & elemRegion ) {
        allFaceElementRegions.emplace_back( elemRegion.getName() );
      } );

    finiteElement::fillSparsity< serialPolicy,
                                 FaceElementSubRegion,
                                 SolidMechanicsLagrangianFEMKernels::QuasiStatic >(
      mesh,
      allFaceElementRegions,
      this->getDiscretizationName(),
      dofNumber,
      dofManager.rankOffset(),
      pattern,
      rowSizes );
  }
  finiteElement::fillSparsity< serialPolicy,
                               CellElementSubRegion,
                               SolidMechanicsLagrangianFEMKernels::QuasiStatic >(
    mesh,
    targetRegionNames(),
    this->getDiscretizationName(),
    dofNumber,
    dofManager.rankOffset(),
    pattern,
    rowSizes );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
}

void
SolidMechanicsLagrangianFEM::AssembleSystem(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  localMatrix.setValues< parallelDevicePolicy< 32 > >( 0 );
  localRhs.setValues< parallelDevicePolicy< 32 > >( 0 );

  if( m_effectiveStress == 1 )
  {
    GEOSX_UNUSED_VAR( dt );
    AssemblyLaunch< constitutive::PoroElasticBase,
                    SolidMechanicsLagrangianFEMKernels::QuasiStaticPoroElastic >(
      domain,
      dofManager,
      localMatrix,
      localRhs );
  }
  else
  {
    if( m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
    {
      GEOSX_UNUSED_VAR( dt );
      AssemblyLaunch< constitutive::SolidBase,
                      SolidMechanicsLagrangianFEMKernels::QuasiStatic >( domain,
                                                                         dofManager,
                                                                         localMatrix,
                                                                         localRhs );
    }
    else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
    {
      AssemblyLaunch< constitutive::SolidBase,
                      SolidMechanicsLagrangianFEMKernels::ImplicitNewmark >(
        domain,
        dofManager,
        localMatrix,
        localRhs,
        m_newmarkGamma,
        m_newmarkBeta,
        m_massDamping,
        m_stiffnessDamping,
        dt );
    }
  }

  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK_0( "After SolidMechanicsLagrangianFEM::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    //    std::cout<< matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    //    std::cout<< rhs;
  }
}

void
SolidMechanicsLagrangianFEM::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager & faceManager = *mesh.getFaceManager();
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  fsManager.Apply( time_n + dt,
                   &domain,
                   "nodeManager",
                   keys::Force,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) ) {
                     bc->ApplyBoundaryConditionToSystem< FieldSpecificationAdd,
                                                         parallelDevicePolicy< 32 > >(
                       targetSet,
                       time_n + dt,
                       targetGroup,
                       keys::TotalDisplacement,  // TODO fix use of dummy
                                                 // name
                       dofKey,
                       dofManager.rankOffset(),
                       localMatrix,
                       localRhs );
                   } );

  CRSApplyTractionBC( time_n + dt, dofManager, domain, localRhs );

  if( faceManager.hasWrapper( "ChomboPressure" ) )
  {
    fsManager.ApplyFieldValue( time_n, &domain, "faceManager", "ChomboPressure" );
    ApplyChomboPressure( dofManager, domain, m_localRhs );
  }

  ApplyDisplacementBC_implicit( time_n + dt,
                                dofManager,
                                domain,
                                localMatrix,
                                localRhs );
}

real64
SolidMechanicsLagrangianFEM::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();

  arrayView1d< globalIndex const > const & dofNumber =
    nodeManager.getReference< array1d< globalIndex > >(
      dofManager.getKey( keys::TotalDisplacement ) );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< integer const > const & ghostRank = nodeManager.ghostRank();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  forAll< parallelDevicePolicy<> >(
    nodeManager.size(),
    [localRhs, localSum, dofNumber, rankOffset, ghostRank] GEOSX_HOST_DEVICE(
      localIndex const k ) {
      if( ghostRank[k] < 0 )
      {
        localIndex const localRow =
          LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
        for( localIndex dim = 0; dim < 3; ++dim )
        {
          localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
        }
      }
    } );

  real64 const localResidualNorm[2] = { localSum.get(), this->m_maxForce };

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = { 0, 0 };

  int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  array1d< real64 > globalValues( size * 2 );

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm, 2, globalValues.data(), 2, 0, MPI_COMM_GEOSX );

  if( rank == 0 )
  {
    for( int r = 0; r < size; ++r )
    {
      // sum/max across all ranks
      globalResidualNorm[0] += globalValues[r * 2];
      globalResidualNorm[1] =
        std::max( globalResidualNorm[1], globalValues[r * 2 + 1] );
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );

  real64 const residual = sqrt( globalResidualNorm[0] ) /
    ( globalResidualNorm[1] + 1 );  // the + 1 is for the first
                                    // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    char output[200] = { 0 };
    sprintf( output, "( RSolid ) = (%4.2e) ; ", residual );
    std::cout << output;
  }

  return residual;
}

void
SolidMechanicsLagrangianFEM::ApplySystemSolution(
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor,
  DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               keys::TotalDisplacement,
                               keys::IncrementalDisplacement,
                               -scalingFactor );

  dofManager.addVectorToField( localSolution,
                               keys::TotalDisplacement,
                               keys::TotalDisplacement,
                               -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::IncrementalDisplacement );
  fieldNames["node"].emplace_back( keys::TotalDisplacement );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );
}

void
SolidMechanicsLagrangianFEM::SolveSystem( DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution )
{
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void
SolidMechanicsLagrangianFEM::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();

  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & incdisp =
    nodeManager.incrementalDisplacement();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & disp =
    nodeManager.totalDisplacement();

  // TODO need to finish this rewind
  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(),
                                        [=] GEOSX_HOST_DEVICE( localIndex const a ) {
                                          for( localIndex i = 0; i < 3; ++i )
                                          {
                                            disp( a, i ) -= incdisp( a, i );
                                            incdisp( a, i ) = 0.0;
                                          }
                                        } );

  ResetStressToBeginningOfStep( domain );
}

void
SolidMechanicsLagrangianFEM::ResetStressToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, CellElementSubRegion & subRegion ) {
      SolidBase & constitutiveRelation =
        GetConstitutiveModel< SolidBase >( subRegion,
                                           m_solidMaterialNames[targetIndex] );

      arrayView3d< real64, solid::STRESS_USD > const & stress =
        constitutiveRelation.getStress();

      arrayView3d< real64 const, solid::STRESS_USD > const & stress_n =
        subRegion.getReference< array3d< real64, solid::STRESS_PERMUTATION > >(
          viewKeyStruct::stress_n );

      forAll< parallelDevicePolicy<> >( stress.size( 0 ),
                                        [=] GEOSX_HOST_DEVICE( localIndex const k ) {
                                          for( localIndex a = 0; a < stress.size( 1 );
                                               ++a )
                                          {
                                            for( localIndex i = 0; i < 6; ++i )
                                            {
                                              stress( k, a, i ) = stress_n( k, a, i );
                                            }
                                          }
                                        } );
    } );
}

void
SolidMechanicsLagrangianFEM::ApplyContactConstraint(
  DofManager const & dofManager,
  DomainPartition & domain,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString )
  {
    MeshLevel * const mesh =
      domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
    FaceManager const * const faceManager = mesh->getFaceManager();
    NodeManager * const nodeManager = mesh->getNodeManager();
    ElementRegionManager * const elemManager = mesh->getElemManager();

    ConstitutiveManager const * const constitutiveManager =
      domain.GetGroup< ConstitutiveManager >( keys::ConstitutiveManager );

    ContactRelationBase const * const contactRelation =
      constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

    real64 const contactStiffness = contactRelation->stiffness();

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
      nodeManager->totalDisplacement();
    arrayView1d< R1Tensor > const & fc =
      nodeManager->getReference< array1d< R1Tensor > >(
        viewKeyStruct::contactForceString );
    fc.setValues< serialPolicy >( { 0, 0, 0 } );

    arrayView2d< real64 const > const & faceNormal = faceManager->faceNormal();
    ArrayOfArraysView< localIndex const > const & facesToNodes =
      faceManager->nodeList().toViewConst();

    string const dofKey = dofManager.getKey( keys::TotalDisplacement );
    arrayView1d< globalIndex > const & nodeDofNumber =
      nodeManager->getReference< globalIndex_array >( dofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // TODO: this bound may need to change
    constexpr localIndex maxNodexPerFace = 4;
    constexpr localIndex maxDofPerElem = maxNodexPerFace * 3 * 2;

    elemManager->forElementSubRegions< FaceElementSubRegion >(
      [&]( FaceElementSubRegion & subRegion ) {
        arrayView1d< real64 > const & area = subRegion.getElementArea();
        arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

        // TODO: use parallel policy?
        forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe ) {
          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
          Nbar -= faceNormal[elemsToFaces[kfe][1]];
          Nbar.Normalize();

          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const kf1 = elemsToFaces[kfe][1];
          localIndex const numNodesPerFace = facesToNodes.sizeOfArray( kf0 );
          real64 const Ja = area[kfe] / numNodesPerFace;

          stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace * 3 * 2 );
          stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace * 3 * 2 );
          stackArray2d< real64, maxDofPerElem * maxDofPerElem > dRdP(
            numNodesPerFace * 3 * 2,
            numNodesPerFace * 3 * 2 );

          for( localIndex a = 0; a < numNodesPerFace; ++a )
          {
            R1Tensor penaltyForce = Nbar;
            localIndex const node0 = facesToNodes[kf0][a];
            localIndex const node1 =
              facesToNodes[kf1][a == 0 ? a : numNodesPerFace - a];
            R1Tensor gap = u[node1];
            gap -= u[node0];
            real64 const gapNormal = Dot( gap, Nbar );

            for( int i = 0; i < 3; ++i )
            {
              rowDOF[3 * a + i] = nodeDofNumber[node0] + i;
              rowDOF[3 * ( numNodesPerFace + a ) + i] = nodeDofNumber[node1] + i;
            }

            if( gapNormal < 0 )
            {
              penaltyForce *= -contactStiffness * gapNormal * Ja;
              for( int i = 0; i < 3; ++i )
              {
                fc[node0] -= penaltyForce;
                fc[node1] += penaltyForce;
                nodeRHS[3 * a + i] -= penaltyForce[i];
                nodeRHS[3 * ( numNodesPerFace + a ) + i] += penaltyForce[i];

                dRdP( 3 * a + i, 3 * a + i ) -=
                  contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3 * a + i, 3 * ( numNodesPerFace + a ) + i ) +=
                  contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3 * ( numNodesPerFace + a ) + i, 3 * a + i ) +=
                  contactStiffness * Ja * Nbar[i] * Nbar[i];
                dRdP( 3 * ( numNodesPerFace + a ) + i,
                      3 * ( numNodesPerFace + a ) + i ) -=
                  contactStiffness * Ja * Nbar[i] * Nbar[i];
              }
            }
          }

          for( localIndex idof = 0; idof < numNodesPerFace * 3 * 2; ++idof )
          {
            localIndex const localRow =
              LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              localMatrix.addToRowBinarySearchUnsorted< serialAtomic >(
                localRow,
                rowDOF.data(),
                dRdP[idof].dataIfContiguous(),
                numNodesPerFace * 3 * 2 );
              RAJA::atomicAdd( serialAtomic {}, &localRhs[localRow], nodeRHS[idof] );
            }
          }
        } );
      } );
  }
}

real64
SolidMechanicsLagrangianFEM::ScalingForSystemSolution(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;
  real64 scalingFactor = 1.0;

  GEOSX_UNUSED_VAR( domain )
  GEOSX_UNUSED_VAR( dofManager )
  GEOSX_UNUSED_VAR( localSolution )

  //  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  //  FaceManager const * const faceManager = mesh->getFaceManager();
  //  NodeManager const * const nodeManager = mesh->getNodeManager();
  //  ElementRegionManager const * const elemManager = mesh->getElemManager();
  //
  //
  //  arrayView1d<R1Tensor const> const & u = nodeManager->getReference< array1d<R1Tensor> >( keys::TotalDisplacement );
  //
  //  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  //  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();
  //
  //  string const dofKey = dofManager.getKey( keys::TotalDisplacement );
  //  arrayView1d<globalIndex> const & nodeDofNumber = nodeManager->getReference<globalIndex_array>( dofKey );
  //
  //  globalIndex const rankOffset = dofManager.rankOffset();
  //
  //  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion const * const subRegion )->void
  //  {
  //      arrayView1d<integer const> const & ghostRank = subRegion->ghostRank();
  //      arrayView1d<real64 const > const & area = subRegion->getElementArea();
  //      arrayView2d< localIndex const > const & elemsToFaces = subRegion->faceList();
  //
  //      RAJA::ReduceMin<RAJA::seq_reduce, real64> newScaleFactor(1.0);
  //      forAll<serialPolicy>( subRegion->size(), [=] ( localIndex const kfe )
  //      {
  //
  //        if( ghostRank[kfe] < 0 )
  //        {
  //          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
  //          Nbar -= faceNormal[elemsToFaces[kfe][1]];
  //          Nbar.Normalize();
  //
  //          localIndex const kf0 = elemsToFaces[kfe][0];
  //          localIndex const kf1 = elemsToFaces[kfe][1];
  //          localIndex const numNodesPerFace=facesToNodes[kf0].size();
  //          localIndex const * const nodelist0 = facesToNodes[kf0];
  //          localIndex const * const nodelist1 = facesToNodes[kf1];
  //
  //          for( localIndex a=0 ; a<numNodesPerFace ; ++a )
  //          {
  //            R1Tensor penaltyForce = Nbar;
  //            localIndex const node0 = facesToNodes[kf0][a];
  //            localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
  //
  //            localIndex const lid0 = nodeDofNumber[node0] - rankOffset;
  //            localIndex const lid1 = nodeDofNumber[node1] - rankOffset;
  //            GEOSX_ASSERT( lid0 >= 0 && lid1 >=0 ); // since vectors are partitioned same as the mesh
  //
  //            R1Tensor gap = u[node1];
  //            gap -= u[node0];
  //            real64 const gapNormal = Dot(gap,Nbar);
  //
  //            R1Tensor deltaGap ;
  //            for( int i=0 ; i<3 ; ++i )
  //            {
  //              deltaGap[i] = localSolution[lid1] - localSolution[lid0];
  //            }
  //            real64 const deltaGapNormal = Dot( deltaGap, Nbar );
  //
  //            if( ( gapNormal + deltaGapNormal ) * gapNormal < 0 )
  //            {
  //              newScaleFactor.min( - gapNormal / deltaGapNormal * 1.05 );
  //              newScaleFactor.min( 1.0 );
  //              std::cout<< "gapNormal, deltaGapNormal, scaleFactor, newGap = "<<
  //                  gapNormal<<", "<<deltaGapNormal<<", "<<newScaleFactor.get()<<",
  // "<<gapNormal+newScaleFactor.get()*deltaGapNormal<<std::endl;
  //            }
  //          }
  //        }
  //      });
  //      scalingFactor = newScaleFactor.get();
  //  });

  return scalingFactor;
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        SolidMechanicsLagrangianFEM,
                        string const &,
                        dataRepository::Group * const )
}  // namespace geosx
