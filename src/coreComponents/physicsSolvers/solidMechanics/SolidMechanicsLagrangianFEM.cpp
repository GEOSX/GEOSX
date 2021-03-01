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
#include "managers/GeosxState.hpp"
#include "managers/ProblemManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SolidMechanicsLagrangianFEM::SolidMechanicsLagrangianFEM( const string & name,
                                                          Group * const parent ):
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
  m_targetNodes(),
  m_iComm(),
  m_effectiveStress( 0 )
{
  m_sendOrReceiveNodes.setName( "SolidMechanicsLagrangianFEM::m_sendOrReceiveNodes" );
  m_nonSendOrReceiveNodes.setName( "SolidMechanicsLagrangianFEM::m_nonSendOrReceiveNodes" );
  m_targetNodes.setName( "SolidMechanicsLagrangianFEM::m_targetNodes" );

  registerWrapper( viewKeyStruct::newmarkGammaString(), &m_newmarkGamma ).
    setApplyDefaultValue( 0.5 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic time integration option" );

  registerWrapper( viewKeyStruct::newmarkBetaString(), &m_newmarkBeta ).
    setApplyDefaultValue( 0.25 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time integration option. "
                    "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing." );

  registerWrapper( viewKeyStruct::massDampingString(), &m_massDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of mass based damping coefficient. " );

  registerWrapper( viewKeyStruct::stiffnessDampingString(), &m_stiffnessDamping ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value of stiffness based damping coefficient. " );

  registerWrapper( viewKeyStruct::timeIntegrationOptionString(), &m_timeIntegrationOption ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( m_timeIntegrationOption ).
    setDescription( "Time integration method. Options are:\n* " + EnumStrings< TimeIntegrationOption >::concat( "\n* " ) );

  registerWrapper( viewKeyStruct::useVelocityEstimateForQSString(), &m_useVelocityEstimateForQS ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag to indicate the use of the incremental displacement from the previous step as an "
                    "initial estimate for the incremental displacement of the current step." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed after some other event is executed. "
                    "For example, if a SurfaceGenerator is specified, it will be executed after the mechanics solve. "
                    "However if a new surface is generated, then the mechanics solve must be executed again due to the "
                    "change in topology." );

  registerWrapper( viewKeyStruct::strainTheoryString(), &m_strainTheory ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Indicates whether or not to use "
                    "`Infinitesimal Strain Theory <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
                    "`Finite Strain Theory <https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs are:\n"
                    " 0 - Infinitesimal Strain \n"
                    " 1 - Finite Strain" );

  registerWrapper( viewKeyStruct::solidMaterialNamesString(), &m_solidMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setApplyDefaultValue( viewKeyStruct::noContactRelationNameString() ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxForceString(), &m_maxForce ).
    setInputFlag( InputFlags::FALSE ).
    setDescription( "The maximum force contribution in the problem domain." );

  registerWrapper( viewKeyStruct::effectiveStressString(), &m_effectiveStress ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Apply fluid pressure to produce effective stress when integrating stress." );

}

void SolidMechanicsLagrangianFEM::postProcessInput()
{
  SolverBase::postProcessInput();

  checkModelNames( m_solidMaterialNames, viewKeyStruct::solidMaterialNamesString() );

  LinearSolverParameters & linParams = m_linearSolverParameters.get();
  linParams.isSymmetric = true;
  linParams.dofsPerNode = 3;
  linParams.amg.separateComponents = true;
}

SolidMechanicsLagrangianFEM::~SolidMechanicsLagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsLagrangianFEM::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodes = meshBody.getMeshLevel( 0 ).getNodeManager();

    nodes.registerWrapper< array2d< real64, nodes::TOTAL_DISPLACEMENT_PERM > >( keys::TotalDisplacement ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the total displacements on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::INCR_DISPLACEMENT_PERM > >( keys::IncrementalDisplacement ).
      setPlotLevel( PlotLevel::LEVEL_3 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the incremental displacements for the current time step on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::VELOCITY_PERM > >( keys::Velocity ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the current velocity on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64, nodes::ACCELERATION_PERM > >( keys::Acceleration ).
      setPlotLevel( PlotLevel::LEVEL_1 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the current acceleration on the nodes. This array also is used "
                      "to hold the summation of nodal forces resulting from the governing equations." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::forceExternalString() ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the external forces on the nodes. This includes any boundary"
                      " conditions as well as coupling forces such as hydraulic forces." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array1d< real64 > >( keys::Mass ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the mass on the nodes." );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::vTildeString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the velocity predictors on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::uhatTildeString() ).
      setPlotLevel( PlotLevel::NOPLOT ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the incremental displacement predictors on the nodes." ).
      reference().resizeDimension< 1 >( 3 );

    nodes.registerWrapper< array2d< real64 > >( viewKeyStruct::contactForceString() ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setRegisteringObjects( this->getName()).
      setDescription( "An array that holds the contact force." ).
      reference().resizeDimension< 1 >( 3 );

    ElementRegionManager &
    elementRegionManager = meshBody.getMeshLevel( 0 ).getElemManager();
    elementRegionManager.forElementSubRegions< CellElementSubRegion >( [&]( CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );

      subRegion.registerWrapper< SortedArray< localIndex > >( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ).
        setPlotLevel( PlotLevel::NOPLOT ).
        setRestartFlags( RestartFlags::NO_WRITE );
    } );

  } );
}


void SolidMechanicsLagrangianFEM::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();

  // Validate solid models in target regions
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    validateModelMapping< SolidBase >( meshLevel.getElemManager(), m_solidMaterialNames );
  }

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteElementDiscretizationManager const &
  feDiscretizationManager = numericalMethodManager.getFiniteElementDiscretizationManager();

  FiniteElementDiscretization const &
  feDiscretization = feDiscretizationManager.getGroup< FiniteElementDiscretization >( m_discretizationName );
  GEOSX_UNUSED_VAR( feDiscretization );
}



template< typename ... PARAMS >
real64 SolidMechanicsLagrangianFEM::explicitKernelDispatch( MeshLevel & mesh,
                                                            arrayView1d< string const > const & targetRegions,
                                                            string const & finiteElementName,
                                                            arrayView1d< string const > const & constitutiveNames,
                                                            PARAMS && ... params )
{
  GEOSX_MARK_FUNCTION;
  real64 rval = 0;
  if( m_strainTheory==0 )
  {
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion,
                                           SolidMechanicsLagrangianFEMKernels::ExplicitSmallStrain >( mesh,
                                                                                                      targetRegions,
                                                                                                      finiteElementName,
                                                                                                      constitutiveNames,
                                                                                                      "geosx::SolidMechanicsLagrangianFEMKernels::ExplicitSmallStrain",
                                                                                                      std::forward< PARAMS >( params )... );
  }
  else if( m_strainTheory==1 )
  {
    rval = finiteElement::
             regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                           constitutive::SolidBase,
                                           CellElementSubRegion,
                                           SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrain >( mesh,
                                                                                                       targetRegions,
                                                                                                       finiteElementName,
                                                                                                       constitutiveNames,
                                                                                                       "geosx::SolidMechanicsLagrangianFEMKernels::ExplicitFiniteStrain",
                                                                                                       std::forward< PARAMS >( params )... );
  }
  else
  {
    GEOSX_ERROR( "Invalid option for strain theory (0 = infinitesimal strain, 1 = finite strain" );
  }
  return rval;
}


void SolidMechanicsLagrangianFEM::updateIntrinsicNodalData( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain->getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodes = mesh.getNodeManager();

  ElementRegionManager const & elementRegionManager = mesh.getElemManager();

  arrayView1d< real64 > & mass = nodes.getReference< array1d< real64 > >( keys::Mass );
  mass.setValues< serialPolicy >( 0.0 );

  arrayView1d< integer const > const & nodeGhostRank = nodes.ghostRank();

  // to fill m_sendOrReceiveNodes and m_nonSendOrReceiveNodes, we first insert
  // the nodes one-by-one in the following std::sets. Then, when all the nodes
  // have been collected, we do a batch insertion into m_sendOrReceiveNodes and
  // m_nonSendOrReceiveNodes
  std::set< localIndex > tmpSendOrReceiveNodes;
  std::set< localIndex > tmpNonSendOrReceiveNodes;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
  rho = elementRegionManager.constructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( "density",
                                                                                                              targetRegionNames(),
                                                                                                              solidMaterialNames() );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const er,
                                       ElementRegionBase const & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr,
                                                                       CellElementSubRegion const & elementSubRegion )
    {
      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( this->getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

        real64 N[numNodesPerElem];
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 elemMass = 0;
          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            elemMass += rho[er][esr][k][q] * detJ[k][q];
            FE_TYPE::calcN( q, N );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[er][esr][k][q] * detJ[k][q] * N[a];
            }
          }


          for( localIndex a=0; a<elementSubRegion.numNodesPerElement(); ++a )
          {
            if( nodeGhostRank[elemsToNodes[k][a]] >= -1 )
            {
              tmpSendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
            else
            {
              tmpNonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
          }
        }
      } );
    } );
  } );
  m_sendOrReceiveNodes.insert( tmpSendOrReceiveNodes.begin(),
                               tmpSendOrReceiveNodes.end() );
  m_nonSendOrReceiveNodes.insert( tmpNonSendOrReceiveNodes.begin(),
                                  tmpNonSendOrReceiveNodes.end() );
  m_targetNodes = m_sendOrReceiveNodes;
  m_targetNodes.insert( m_nonSendOrReceiveNodes.begin(),
                        m_nonSendOrReceiveNodes.end() );
}

void SolidMechanicsLagrangianFEM::initializePostInitialConditionsPreSubGroups()
{
  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodes = mesh.getNodeManager();

  ElementRegionManager & elementRegionManager = mesh.getElemManager();

  arrayView1d< real64 > & mass = nodes.getReference< array1d< real64 > >( keys::Mass );

  arrayView1d< integer const > const & nodeGhostRank = nodes.ghostRank();

  // to fill m_sendOrReceiveNodes and m_nonSendOrReceiveNodes, we first insert
  // the nodes one-by-one in the following std::sets. Then, when all the nodes
  // have been collected, we do a batch insertion into m_sendOrReceiveNodes and
  // m_nonSendOrReceiveNodes
  std::set< localIndex > tmpSendOrReceiveNodes;
  std::set< localIndex > tmpNonSendOrReceiveNodes;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
  rho = elementRegionManager.constructMaterialViewAccessor< array2d< real64 >, arrayView2d< real64 const > >( "density",
                                                                                                              targetRegionNames(),
                                                                                                              solidMaterialNames() );

  forTargetRegionsComplete( mesh, [&]( localIndex const,
                                       localIndex const er,
                                       ElementRegionBase & elemRegion )
  {
    elemRegion.forElementSubRegionsIndex< CellElementSubRegion >( [&]( localIndex const esr, CellElementSubRegion & elementSubRegion )
    {
      SortedArray< localIndex > & elemsAttachedToSendOrReceiveNodes = getElemsAttachedToSendOrReceiveNodes( elementSubRegion );
      SortedArray< localIndex > & elemsNotAttachedToSendOrReceiveNodes = getElemsNotAttachedToSendOrReceiveNodes( elementSubRegion );

      std::set< localIndex > tmpElemsAttachedToSendOrReceiveNodes;
      std::set< localIndex > tmpElemsNotAttachedToSendOrReceiveNodes;

      elemsAttachedToSendOrReceiveNodes.setName(
        "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes["
        + std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

      elemsNotAttachedToSendOrReceiveNodes.setName(
        "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes["
        + std::to_string( er ) + "][" + std::to_string( esr ) + "]" );

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();
      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

      finiteElement::FiniteElementBase const &
      fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( getDiscretizationName() );
      finiteElement::dispatch3D( fe,
                                 [&] ( auto const finiteElement )
      {
        using FE_TYPE = TYPEOFREF( finiteElement );

        constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
        constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

        real64 N[numNodesPerElem];
        for( localIndex k=0; k < elemsToNodes.size( 0 ); ++k )
        {
          real64 elemMass = 0;
          for( localIndex q=0; q<numQuadraturePointsPerElem; ++q )
          {
            elemMass += rho[er][esr][k][q] * detJ[k][q];
            FE_TYPE::calcN( q, N );

            for( localIndex a=0; a< numNodesPerElem; ++a )
            {
              mass[elemsToNodes[k][a]] += rho[er][esr][k][q] * detJ[k][q] * N[a];
            }
          }

          bool isAttachedToGhostNode = false;
          for( localIndex a=0; a<elementSubRegion.numNodesPerElement(); ++a )
          {
            if( nodeGhostRank[elemsToNodes[k][a]] >= -1 )
            {
              isAttachedToGhostNode = true;
              tmpSendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
            else
            {
              tmpNonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
            }
          }

          if( isAttachedToGhostNode )
          {
            tmpElemsAttachedToSendOrReceiveNodes.insert( k );
          }
          else
          {
            tmpElemsNotAttachedToSendOrReceiveNodes.insert( k );
          }
        }
      } );
      elemsAttachedToSendOrReceiveNodes.insert( tmpElemsAttachedToSendOrReceiveNodes.begin(),
                                                tmpElemsAttachedToSendOrReceiveNodes.end() );
      elemsNotAttachedToSendOrReceiveNodes.insert( tmpElemsNotAttachedToSendOrReceiveNodes.begin(),
                                                   tmpElemsNotAttachedToSendOrReceiveNodes.end() );

    } );
  } );

  m_sendOrReceiveNodes.insert( tmpSendOrReceiveNodes.begin(),
                               tmpSendOrReceiveNodes.end() );
  m_nonSendOrReceiveNodes.insert( tmpNonSendOrReceiveNodes.begin(),
                                  tmpNonSendOrReceiveNodes.end() );
  m_targetNodes = m_sendOrReceiveNodes;
  m_targetNodes.insert( m_nonSendOrReceiveNodes.begin(),
                        m_nonSendOrReceiveNodes.end() );
}



real64 SolidMechanicsLagrangianFEM::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                const int cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator = this->getParent().getGroupPointer< SolverBase >( "SurfaceGen" );

  if( m_timeIntegrationOption == TimeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = explicitStep( time_n, dt, cycleNumber, domain );

    if( surfaceGenerator!=nullptr )
    {
      surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain );
    }
  }
  else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    int const maxNumResolves = m_maxNumResolves;
    int locallyFractured = 0;
    int globallyFractured = 0;
    implicitStepSetup( time_n, dt, domain );
    for( int solveIter=0; solveIter<maxNumResolves; ++solveIter )
    {
      setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

      dtReturn = nonlinearImplicitStep( time_n,
                                        dt,
                                        cycleNumber,
                                        domain );

      if( surfaceGenerator!=nullptr )
      {
        if( surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain ) > 0 )
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
    implicitStepComplete( time_n, dt, domain );
  }

  return dtReturn;
}

real64 SolidMechanicsLagrangianFEM::explicitStep( real64 const & time_n,
                                                  real64 const & dt,
                                                  const int GEOSX_UNUSED_PARAM( cycleNumber ),
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  #define USE_PHYSICS_LOOP

  // updateIntrinsicNodalData(domain);

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodes = mesh.getNodeManager();

  // save previous constitutive state data in preparation for next timestep
  forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          CellElementSubRegion & subRegion )
  {
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, m_solidMaterialNames[targetIndex] );
    constitutiveRelation.saveConvergedState();
  } );

  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();

  arrayView1d< real64 const > const & mass = nodes.getReference< array1d< real64 > >( keys::Mass );
  arrayView2d< real64, nodes::VELOCITY_USD > const & vel = nodes.velocity();

  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodes.totalDisplacement();
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat = nodes.incrementalDisplacement();
  arrayView2d< real64, nodes::ACCELERATION_USD > const & acc = nodes.acceleration();

  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::Velocity );
  fieldNames["node"].emplace_back( keys::Acceleration );

  getGlobalState().getCommunicationTools().synchronizePackSendRecvSizes( fieldNames, mesh, domain.getNeighbors(), m_iComm, true );

  fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Acceleration );

  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt/2 );

  fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

  fsManager.applyFieldValue( time_n + dt,
                             domain, "nodeManager",
                             NodeManager::viewKeyStruct::totalDisplacementString(),
                             [&]( FieldSpecificationBase const & bc,
                                  SortedArrayView< localIndex const > const & targetSet )
  {
    integer const component = bc.getComponent();
    forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
                                            [=] GEOSX_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      vel( a, component ) = u( a, component );
    } );
  },
                             [&]( FieldSpecificationBase const & bc,
                                  SortedArrayView< localIndex const > const & targetSet )
  {
    integer const component = bc.getComponent();
    forAll< parallelDevicePolicy< 1024 > >( targetSet.size(),
                                            [=] GEOSX_DEVICE ( localIndex const i )
    {
      localIndex const a = targetSet[ i ];
      uhat( a, component ) = u( a, component ) - vel( a, component );
      vel( a, component )  = uhat( a, component ) / dt;
    } );
  } );

  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  explicitKernelDispatch( mesh,
                          targetRegionNames(),
                          this->getDiscretizationName(),
                          m_solidMaterialNames,
                          dt,
                          string( viewKeyStruct::elemsAttachedToSendOrReceiveNodesString() ) );

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes.toViewConst() );

  fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  getGlobalState().getCommunicationTools().synchronizePackSendRecv( fieldNames, mesh, domain.getNeighbors(), m_iComm, true );

  explicitKernelDispatch( mesh,
                          targetRegionNames(),
                          this->getDiscretizationName(),
                          m_solidMaterialNames,
                          dt,
                          string( viewKeyStruct::elemsNotAttachedToSendOrReceiveNodesString() ) );

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes.toViewConst() );

  fsManager.applyFieldValue< parallelDevicePolicy< 1024 > >( time_n, domain, "nodeManager", keys::Velocity );

  getGlobalState().getCommunicationTools().synchronizeUnpack( mesh, domain.getNeighbors(), m_iComm, true );

  return dt;
}



void SolidMechanicsLagrangianFEM::applyDisplacementBCImplicit( real64 const time,
                                                               DofManager const & dofManager,
                                                               DomainPartition & domain,
                                                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                               arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  FieldSpecificationManager const & fsManager = getGlobalState().getFieldSpecificationManager();

  fsManager.apply( time,
                   domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const fieldName )
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time,
                                                                     targetGroup,
                                                                     fieldName,
                                                                     dofKey,
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );
}

void SolidMechanicsLagrangianFEM::crsApplyTractionBC( real64 const time,
                                                      DofManager const & dofManager,
                                                      DomainPartition & domain,
                                                      arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();
  FunctionManager const & functionManager = getGlobalState().getFunctionManager();

  FaceManager const & faceManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getFaceManager();
  NodeManager const & nodeManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getNodeManager();

  arrayView1d< real64 const > const faceArea  = faceManager.getReference< real64_array >( "faceArea" );
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const blockLocalDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
  globalIndex const dofRankOffset = dofManager.rankOffset();

  fsManager.apply( time,
                   domain,
                   "faceManager",
                   string( "Traction" ),
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group &,
                        string const & )
  {
    string const & functionName = bc.getFunctionName();

    globalIndex_array nodeDOF;
    real64_array nodeRHS;
    integer const component = bc.getComponent();


    if( functionName.empty() || functionManager.getGroup< FunctionBase >( functionName ).isFunctionOfTime() == 2 )
    {
      real64 value = bc.getScale();
      if( !functionName.empty() )
      {
        FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
        value *= function.evaluate( &time );
      }

      forAll< parallelDevicePolicy< 32 > >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex const kf = targetSet[ i ];
        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a<numNodes; ++a )
        {
          localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component - dofRankOffset;
          if( dof < 0 || dof >= localRhs.size() )
            continue;
          RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], value * faceArea[kf] / numNodes );
        }
      } );
    }
    else
    {
      FunctionBase const & function = functionManager.getGroup< FunctionBase >( functionName );
      array1d< real64 > resultsArray( targetSet.size() );
      resultsArray.setName( "SolidMechanicsLagrangianFEM::TractionBC function results" );
      function.evaluate( faceManager, time, targetSet, resultsArray );
      arrayView1d< real64 const > const results = resultsArray;

      forAll< parallelDevicePolicy< 32 > >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex const kf = targetSet[ i ];
        localIndex const numNodes = faceToNodeMap.sizeOfArray( kf );
        for( localIndex a=0; a<numNodes; ++a )
        {
          localIndex const dof = blockLocalDofNumber[ faceToNodeMap( kf, a ) ] + component - dofRankOffset;
          if( dof < 0 || dof >= localRhs.size() )
            continue;
          RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], results[ kf ] * faceArea[kf] / numNodes );
        }
      } );
    }
  } );
}

void SolidMechanicsLagrangianFEM::applyChomboPressure( DofManager const & dofManager,
                                                       DomainPartition & domain,
                                                       arrayView1d< real64 > const & localRhs )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 const > const faceArea  = faceManager.faceArea();
  arrayView2d< real64 const > const faceNormal  = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const dofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
  arrayView1d< real64 const > const facePressure = faceManager.getReference< array1d< real64 > >( "ChomboPressure" );

  forAll< serialPolicy >( faceManager.size(), [=] ( localIndex const kf )
  {
    int const numNodes = LvArray::integerConversion< int >( faceToNodeMap.sizeOfArray( kf ));
    for( int a=0; a<numNodes; ++a )
    {
      localIndex const dof = dofNumber[ faceToNodeMap( kf, a ) ];
      if( dof < 0 || dof >= localRhs.size() )
        continue;

      for( int component=0; component<3; ++component )
      {
        real64 const value = -facePressure[ kf ] * faceNormal( kf, component ) * faceArea[kf] / numNodes;
        localRhs[ dof + component ] += value;
      }
    }
  } );
}



void
SolidMechanicsLagrangianFEM::
  implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                     real64 const & dt,
                     DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const, nodes::VELOCITY_USD > const v_n = nodeManager.velocity();
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const uhat = nodeManager.incrementalDisplacement();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const disp = nodeManager.totalDisplacement();

  localIndex const numNodes = nodeManager.size();

  if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d< real64 const, nodes::ACCELERATION_USD > const a_n = nodeManager.acceleration();
    arrayView1d< R1Tensor > const vtilde = nodeManager.getReference< array1d< R1Tensor > >( solidMechanicsViewKeys.vTilde );
    arrayView1d< R1Tensor > const uhatTilde = nodeManager.getReference< array1d< R1Tensor > >( solidMechanicsViewKeys.uhatTilde );

    real64 const newmarkGamma = this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
    real64 const newmarkBeta = this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

    forAll< parallelDevicePolicy< 32 > >( numNodes, [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      for( int i=0; i<3; ++i )
      {
        vtilde[a][i] = v_n( a, i ) + (1.0-newmarkGamma) * a_n( a, i ) * dt;
        uhatTilde[a][i] = ( v_n( a, i ) + 0.5 * ( 1.0 - 2.0*newmarkBeta ) * a_n( a, i ) * dt ) *dt;
        uhat( a, i ) = uhatTilde[a][i];
        disp( a, i ) += uhatTilde[a][i];
      }
    } );
  }
  else if( this->m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
  {
    if( m_useVelocityEstimateForQS==1 )
    {
      forAll< parallelDevicePolicy< 32 > >( numNodes, [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        for( int i=0; i<3; ++i )
        {
          uhat( a, i ) = v_n( a, i ) * dt;
          disp( a, i ) += uhat( a, i );
        }
      } );
    }
    else
    {
      forAll< parallelDevicePolicy< 32 > >( numNodes, [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        for( int i=0; i<3; ++i )
        {
          uhat( a, i ) = 0.0;
        }
      } );
    }
  }

  ElementRegionManager & elementRegionManager = mesh.getElemManager();
  ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();
  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase >
  constitutiveRelations = elementRegionManager.constructFullConstitutiveAccessor< ConstitutiveBase >( constitutiveManager );

  forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          CellElementSubRegion & subRegion )
  {
    SolidBase const & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, m_solidMaterialNames[targetIndex] );
    constitutiveRelation.saveConvergedState();
  } );

}

void SolidMechanicsLagrangianFEM::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();
  localIndex const numNodes = nodeManager.size();

  arrayView2d< real64, nodes::VELOCITY_USD > const v_n = nodeManager.velocity();
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const uhat  = nodeManager.incrementalDisplacement();

  if( this->m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
  {
    arrayView2d< real64, nodes::ACCELERATION_USD > const a_n = nodeManager.acceleration();
    arrayView1d< R1Tensor const > const vtilde    = nodeManager.getReference< r1_array >( solidMechanicsViewKeys.vTilde );
    arrayView1d< R1Tensor const > const uhatTilde = nodeManager.getReference< r1_array >( solidMechanicsViewKeys.uhatTilde );
    real64 const newmarkGamma = this->getReference< real64 >( solidMechanicsViewKeys.newmarkGamma );
    real64 const newmarkBeta = this->getReference< real64 >( solidMechanicsViewKeys.newmarkBeta );

    RAJA::forall< parallelDevicePolicy<> >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                            [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      for( int i=0; i<3; ++i )
      {
        a_n( a, i ) = 1.0 / ( newmarkBeta * dt*dt) * ( uhat( a, i ) - uhatTilde[a][i] );
        v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n( a, i ) * dt;
      }
    } );
  }
  else if( this->m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic && dt > 0.0 )
  {
    RAJA::forall< parallelDevicePolicy<> >( RAJA::TypedRangeSegment< localIndex >( 0, numNodes ),
                                            [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      for( int i=0; i<3; ++i )
      {
        v_n[a][i] = uhat( a, i ) / dt;
      }
    } );
  }

  // save (converged) constitutive state data
  forTargetSubRegions< CellElementSubRegion >( mesh, [&]( localIndex const targetIndex,
                                                          CellElementSubRegion & subRegion )
  {
    SolidBase & constitutiveRelation = getConstitutiveModel< SolidBase >( subRegion, m_solidMaterialNames[targetIndex] );
    constitutiveRelation.saveConvergedState();
  } );

}

void SolidMechanicsLagrangianFEM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                             DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  dofManager.addField( keys::TotalDisplacement,
                       DofManager::Location::Node,
                       3,
                       targetRegionNames() );

  dofManager.addCoupling( keys::TotalDisplacement,
                          keys::TotalDisplacement,
                          DofManager::Connector::Elem );
}


void SolidMechanicsLagrangianFEM::setupSystem( DomainPartition & domain,
                                               DofManager & dofManager,
                                               CRSMatrix< real64, globalIndex > & localMatrix,
                                               array1d< real64 > & localRhs,
                                               array1d< real64 > & localSolution,
                                               bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();
  arrayView1d< globalIndex const > const
  dofNumber = nodeManager.getReference< globalIndex_array >( dofManager.getKey( keys::TotalDisplacement ) );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3*1.2 );

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString() )
  {
    ElementRegionManager const & elemManager = mesh.getElemManager();
    array1d< string > allFaceElementRegions;
    elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elemRegion )
    {
      allFaceElementRegions.emplace_back( elemRegion.getName() );
    } );

    finiteElement::
      fillSparsity< FaceElementSubRegion,
                    SolidMechanicsLagrangianFEMKernels::QuasiStatic >( mesh,
                                                                       allFaceElementRegions,
                                                                       this->getDiscretizationName(),
                                                                       dofNumber,
                                                                       dofManager.rankOffset(),
                                                                       sparsityPattern );

  }
  finiteElement::
    fillSparsity< CellElementSubRegion,
                  SolidMechanicsLagrangianFEMKernels::QuasiStatic >( mesh,
                                                                     targetRegionNames(),
                                                                     this->getDiscretizationName(),
                                                                     dofNumber,
                                                                     dofManager.rankOffset(),
                                                                     sparsityPattern );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );


}

void SolidMechanicsLagrangianFEM::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  localMatrix.setValues< parallelDevicePolicy< 32 > >( 0 );
  localRhs.setValues< parallelDevicePolicy< 32 > >( 0 );

  if( m_effectiveStress==1 )
  {
    GEOSX_UNUSED_VAR( dt );
    assemblyLaunch< constitutive::PoroElasticBase,
                    SolidMechanicsLagrangianFEMKernels::QuasiStaticPoroElastic >( domain,
                                                                                  dofManager,
                                                                                  localMatrix,
                                                                                  localRhs,
                                                                                  "geosx::SolidMechanicsLagrangianFEMKernels::QuasiStaticPoroElastic" );

  }
  else
  {
    if( m_timeIntegrationOption == TimeIntegrationOption::QuasiStatic )
    {
      GEOSX_UNUSED_VAR( dt );
      assemblyLaunch< constitutive::SolidBase,
                      SolidMechanicsLagrangianFEMKernels::QuasiStatic >( domain,
                                                                         dofManager,
                                                                         localMatrix,
                                                                         localRhs,
                                                                         "geosx::SolidMechanicsLagrangianFEMKernels::QuasiStatic" );
    }
    else if( m_timeIntegrationOption == TimeIntegrationOption::ImplicitDynamic )
    {
      assemblyLaunch< constitutive::SolidBase,
                      SolidMechanicsLagrangianFEMKernels::ImplicitNewmark >( domain,
                                                                             dofManager,
                                                                             localMatrix,
                                                                             localRhs,
                                                                             "geosx::SolidMechanicsLagrangianFEMKernels::ImplicitNewmark",
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
    //std::cout<< localMatrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout<< localRhs;
  }
}

void
SolidMechanicsLagrangianFEM::
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager & faceManager = mesh.getFaceManager();
  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );

  fsManager.apply( time_n + dt,
                   domain,
                   "nodeManager",
                   keys::Force,
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationAdd,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time_n + dt,
                                                                     targetGroup,
                                                                     keys::TotalDisplacement,    // TODO fix use of
                                                                                                 // dummy
                                                                                                 // name
                                                                     dofKey,
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );

  crsApplyTractionBC( time_n + dt, dofManager, domain, localRhs );

  if( faceManager.hasWrapper( "ChomboPressure" ) )
  {
    fsManager.applyFieldValue( time_n, domain, "faceManager", "ChomboPressure" );
    applyChomboPressure( dofManager, domain, m_localRhs );
  }

  applyDisplacementBCImplicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

real64
SolidMechanicsLagrangianFEM::
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();

  arrayView1d< globalIndex const > const dofNumber =
    nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( keys::TotalDisplacement ) );
  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< integer const > const ghostRank = nodeManager.ghostRank();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  SortedArrayView< localIndex const > const targetNodes = m_targetNodes.toViewConst();

  forAll< parallelDevicePolicy<> >( targetNodes.size(),
                                    [localRhs, localSum, dofNumber, rankOffset, ghostRank, targetNodes] GEOSX_HOST_DEVICE ( localIndex const k )
  {
    localIndex const nodeIndex = targetNodes[k];
    if( ghostRank[nodeIndex] < 0 )
    {
      localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[nodeIndex] - rankOffset );

      for( localIndex dim = 0; dim < 3; ++dim )
      {
        localSum += localRhs[localRow + dim] * localRhs[localRow + dim];
      }
    }
  } );

  real64 const localResidualNorm[2] = { localSum.get(), this->m_maxForce };

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = {0, 0};

  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  array1d< real64 > globalValues( size * 2 );

  // Everything is done on rank 0
  MpiWrapper::gather( localResidualNorm,
                      2,
                      globalValues.data(),
                      2,
                      0,
                      MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0; r<size; ++r )
    {
      // sum/max across all ranks
      globalResidualNorm[0] += globalValues[r*2];
      globalResidualNorm[1] = std::max( globalResidualNorm[1], globalValues[r*2+1] );
    }
  }

  MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );


  real64 const residual = sqrt( globalResidualNorm[0] )/(globalResidualNorm[1]+1);  // the + 1 is for the first
                                                                                    // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output,
             "( RSolid ) = (%4.2e) ; ",
             residual );
    std::cout<<output;
  }


  return residual;
}



void
SolidMechanicsLagrangianFEM::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
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

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames,
                                                              domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                              domain.getNeighbors(),
                                                              true );
}

void SolidMechanicsLagrangianFEM::solveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void SolidMechanicsLagrangianFEM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & incdisp  = nodeManager.incrementalDisplacement();
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & disp = nodeManager.totalDisplacement();

  // TODO need to finish this rewind
  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    for( localIndex i = 0; i < 3; ++i )
    {
      disp( a, i ) -= incdisp( a, i );
      incdisp( a, i ) = 0.0;
    }
  } );
}


void SolidMechanicsLagrangianFEM::applyContactConstraint( DofManager const & dofManager,
                                                          DomainPartition & domain,
                                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_contactRelationName != viewKeyStruct::noContactRelationNameString() )
  {
    MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
    FaceManager const & faceManager = mesh.getFaceManager();
    NodeManager & nodeManager = mesh.getNodeManager();
    ElementRegionManager & elemManager = mesh.getElemManager();


    ConstitutiveManager const &
    constitutiveManager = domain.getGroup< ConstitutiveManager >( keys::ConstitutiveManager );

    ContactRelationBase const &
    contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_contactRelationName );

    real64 const contactStiffness = contactRelation.stiffness();

    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u = nodeManager.totalDisplacement();
    arrayView2d< real64 > const fc = nodeManager.getReference< array2d< real64 > >( viewKeyStruct::contactForceString() );
    fc.setValues< serialPolicy >( 0 );

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
    ArrayOfArraysView< localIndex const > const facesToNodes = faceManager.nodeList().toViewConst();

    string const dofKey = dofManager.getKey( keys::TotalDisplacement );
    arrayView1d< globalIndex > const nodeDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
    globalIndex const rankOffset = dofManager.rankOffset();

    // TODO: this bound may need to change
    constexpr localIndex maxNodexPerFace = 4;
    constexpr localIndex maxDofPerElem = maxNodexPerFace * 3 * 2;

    elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 > const area = subRegion.getElementArea();
      arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

      // TODO: use parallel policy?
      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        real64 Nbar[ 3 ] = { faceNormal[elemsToFaces[kfe][0]][0] - faceNormal[elemsToFaces[kfe][1]][0],
                             faceNormal[elemsToFaces[kfe][0]][1] - faceNormal[elemsToFaces[kfe][1]][1],
                             faceNormal[elemsToFaces[kfe][0]][2] - faceNormal[elemsToFaces[kfe][1]][2] };

        LvArray::tensorOps::normalize< 3 >( Nbar );


        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const kf1 = elemsToFaces[kfe][1];
        localIndex const numNodesPerFace=facesToNodes.sizeOfArray( kf0 );
        real64 const Ja = area[kfe] / numNodesPerFace;

        stackArray1d< globalIndex, maxDofPerElem > rowDOF( numNodesPerFace*3*2 );
        stackArray1d< real64, maxDofPerElem > nodeRHS( numNodesPerFace*3*2 );
        stackArray2d< real64, maxDofPerElem *maxDofPerElem > dRdP( numNodesPerFace*3*2, numNodesPerFace*3*2 );

        for( localIndex a=0; a<numNodesPerFace; ++a )
        {
          real64 penaltyForce[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( Nbar );
          localIndex const node0 = facesToNodes[kf0][a];
          localIndex const node1 = facesToNodes[kf1][ a==0 ? a : numNodesPerFace-a ];
          real64 gap[ 3 ] = LVARRAY_TENSOROPS_INIT_LOCAL_3( u[node1] );
          LvArray::tensorOps::subtract< 3 >( gap, u[node0] );
          real64 const gapNormal = LvArray::tensorOps::AiBi< 3 >( gap, Nbar );

          for( int i=0; i<3; ++i )
          {
            rowDOF[3*a+i]                     = nodeDofNumber[node0]+i;
            rowDOF[3*(numNodesPerFace + a)+i] = nodeDofNumber[node1]+i;
          }

          if( gapNormal < 0 )
          {
            LvArray::tensorOps::scale< 3 >( penaltyForce, -contactStiffness * gapNormal * Ja );
            for( int i=0; i<3; ++i )
            {
              LvArray::tensorOps::subtract< 3 >( fc[node0], penaltyForce );
              LvArray::tensorOps::add< 3 >( fc[node1], penaltyForce );
              nodeRHS[3*a+i]                     -= penaltyForce[i];
              nodeRHS[3*(numNodesPerFace + a)+i] += penaltyForce[i];

              dRdP( 3*a+i, 3*a+i )                                         -= contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*a+i, 3*(numNodesPerFace + a)+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*(numNodesPerFace + a)+i, 3*a+i )                     += contactStiffness * Ja * Nbar[i] * Nbar[i];
              dRdP( 3*(numNodesPerFace + a)+i, 3*(numNodesPerFace + a)+i ) -= contactStiffness * Ja * Nbar[i] * Nbar[i];
            }
          }
        }

        for( localIndex idof = 0; idof < numNodesPerFace*3*2; ++idof )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[idof] - rankOffset );

          if( localRow >= 0 && localRow < localMatrix.numRows() )
          {
            localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( localRow,
                                                                      rowDOF.data(),
                                                                      dRdP[idof].dataIfContiguous(),
                                                                      numNodesPerFace*3*2 );
            RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow], nodeRHS[idof] );
          }
        }
      } );
    } );
  }
}

real64
SolidMechanicsLagrangianFEM::scalingForSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( domain )
  GEOSX_UNUSED_VAR( dofManager )
  GEOSX_UNUSED_VAR( localSolution )

  return 1.0;
}

REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianFEM, string const &, dataRepository::Group * const )
}
