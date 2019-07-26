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
 * @file SolidMechanicsLagrangianFEM.cpp
 */

#include "SolidMechanicsLagrangianFEM.hpp"

#include <vector>
#include <math.h>
#include <sys/time.h>

#include "common/TimingMacros.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/ElementLibrary/FiniteElement.h"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "codingUtilities/Utilities.hpp"

#include "managers/DomainPartition.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "../../rajaInterface/GEOS_RAJA_Interface.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"


//#define verbose 0 //Need to move this somewhere else

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;


inline void LinearElasticIsotropic_Kernel(R2SymTensor & Dadt, R2SymTensor & TotalStress, R2Tensor & Rot,
                                          localIndex i, real64 bulkModulus, real64 shearModulus,
                                          arrayView1d<real64> & meanStress,
                                          arrayView1d<R2SymTensor> & devStress)
{
  real64 volumeStrain = Dadt.Trace();
  TotalStress = Dadt; 

  meanStress[i] += volumeStrain * bulkModulus;
  TotalStress *= 2.0 * shearModulus;
  
  devStress[i] += TotalStress;
  
  TotalStress.QijAjkQlk(devStress[i],Rot);
  devStress[i] = TotalStress;
  
  TotalStress.PlusIdentity(meanStress[i]);
  
}

SolidMechanicsLagrangianFEM::SolidMechanicsLagrangianFEM( const std::string& name,
                                                          ManagedGroup * const parent ):
  SolverBase( name, parent ),
  m_newmarkGamma(0.5),
  m_newmarkBeta(0.25),
  m_massDamping(0.0),
  m_stiffnessDamping(0.0),
  m_timeIntegrationOptionString(),
  m_timeIntegrationOption(timeIntegrationOption::ExplicitDynamic),
  m_useVelocityEstimateForQS(0),
  m_maxForce(0.0),
  m_maxNumResolves(10),
  m_strainTheory(0),
  m_solidMaterialName(""),
  m_solidMaterialFullIndex(0),
  m_elemsAttachedToSendOrReceiveNodes(),
  m_elemsNotAttachedToSendOrReceiveNodes(),
  m_sendOrReceiveNodes(),
  m_nonSendOrReceiveNodes(),
  m_iComm()
{
  // To generate the schema, multiple solvers of that use this command are constructed
  // Doing this can cause an error in the block setup, so move it to InitializePreSubGroups
  // getLinearSystemRepository()->SetBlockID( BlockIDs::displacementBlock, this->getName() );
  m_sendOrReceiveNodes.setUserCallBack("SolidMechanicsLagrangianFEM::m_sendOrReceiveNodes");
  m_nonSendOrReceiveNodes.setUserCallBack("SolidMechanicsLagrangianFEM::m_nonSendOrReceiveNodes");

  RegisterViewWrapper(viewKeyStruct::newmarkGammaString, &m_newmarkGamma, false )->
    setApplyDefaultValue(0.5)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of :math:`\\gamma` in the Newmark Method for Implicit Dynamic time integration option");

  RegisterViewWrapper(viewKeyStruct::newmarkBetaString, &m_newmarkBeta, false )->
    setApplyDefaultValue(0.25)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of :math:`\\beta` in the Newmark Method for Implicit Dynamic time integration option. "
                   "This should be pow(newmarkGamma+0.5,2.0)/4.0 unless you know what you are doing.");

  RegisterViewWrapper(viewKeyStruct::massDampingString, &m_massDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of mass based damping coefficient. ");

  RegisterViewWrapper(viewKeyStruct::stiffnessDampingString, &m_stiffnessDamping, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Value of stiffness based damping coefficient. ");

  RegisterViewWrapper(viewKeyStruct::timeIntegrationOptionString, &m_timeIntegrationOption, false )->
    setInputFlag(InputFlags::FALSE)->
    setDescription("Time integration enum class value.");

  RegisterViewWrapper(viewKeyStruct::timeIntegrationOptionStringString, &m_timeIntegrationOptionString, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Time integration method. Options are: \n QuasiStatic \n ImplicitDynamic \n ExplicitDynamic");

  RegisterViewWrapper(viewKeyStruct::useVelocityEstimateForQSString, &m_useVelocityEstimateForQS, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Flag to indicate the use of the incremental displacement from the previous step as an "
                    "initial estimate for the incremental displacement of the current step.");

  RegisterViewWrapper(viewKeyStruct::maxNumResolvesString, &m_maxNumResolves, false )->
    setApplyDefaultValue(10)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Value to indicate how many resolves may be executed after some other event is executed. "
                    "For example, if a SurfaceGenerator is specified, it will be executed after the mechanics solve. "
                    "However if a new surface is generated, then the mechanics solve must be executed again due to the "
                    "change in topology.");

  RegisterViewWrapper(viewKeyStruct::strainTheoryString, &m_strainTheory, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription( "Indicates whether or not to use "
                    "`Infinitesimal Strain Theory <https://en.wikipedia.org/wiki/Infinitesimal_strain_theory>`_, or "
                    "`Finite Strain Theory <https://en.wikipedia.org/wiki/Finite_strain_theory>`_. Valid Inputs are:\n"
                    " 0 - Infinitesimal Strain \n"
                    " 1 - Finite Strain");

  RegisterViewWrapper(viewKeyStruct::solidMaterialNameString, &m_solidMaterialName, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription( "The name of the material that should be used in the constitutive updates");


}

void SolidMechanicsLagrangianFEM::PostProcessInput()
{
  if( !m_timeIntegrationOptionString.empty() )
  {
    SetTimeIntegrationOption( m_timeIntegrationOptionString );
  }
}

SolidMechanicsLagrangianFEM::~SolidMechanicsLagrangianFEM()
{
  // TODO Auto-generated destructor stub
}


void SolidMechanicsLagrangianFEM::RegisterDataOnMesh( ManagedGroup * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes = mesh.second->group_cast<MeshBody*>()->getMeshLevel(0)->getNodeManager();


    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::TotalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the total displacements on the nodes.");

    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::IncrementalDisplacement )->
      setPlotLevel(PlotLevel::LEVEL_3)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the incremental displacements for the current time step on the nodes.");

    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::Velocity )->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current velocity on the nodes.");

    nodes->RegisterViewWrapper<array1d<R1Tensor> >( keys::Acceleration )->setPlotLevel(PlotLevel::LEVEL_1)->
      setPlotLevel(PlotLevel::LEVEL_0)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the current acceleration on the nodes. This array also is used "
                      "to hold the summation of nodal forces resulting from the governing equations.");

    nodes->RegisterViewWrapper<array1d<real64> >( keys::Mass )->setPlotLevel(PlotLevel::LEVEL_0)->
        setPlotLevel(PlotLevel::LEVEL_0)->
        setRegisteringObjects(this->getName())->
        setDescription( "An array that holds the mass on the nodes.");

    nodes->RegisterViewWrapper<array1d<R1Tensor> >( viewKeyStruct::vTildeString )->
      setPlotLevel(PlotLevel::NOPLOT)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the velocity predictors on the nodes.");

    nodes->RegisterViewWrapper<array1d<R1Tensor> >( viewKeyStruct::uhatTildeString )->
      setPlotLevel(PlotLevel::NOPLOT)->
      setRegisteringObjects(this->getName())->
      setDescription( "An array that holds the incremental displacement predictors on the nodes.");

    nodes->RegisterViewWrapper<array1d<globalIndex> >( viewKeyStruct::globalDofNumberString )->setPlotLevel(PlotLevel::LEVEL_1);

  }
}


void SolidMechanicsLagrangianFEM::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  SolverBase::InitializePreSubGroups(rootGroup);

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::displacementBlock, this->getName() );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);
  ConstitutiveManager const * const cm = domain->getConstitutiveManager();

  ConstitutiveBase const * const solid  = cm->GetConstitutiveRelation<ConstitutiveBase>( m_solidMaterialName );
  GEOS_ERROR_IF( solid == nullptr, "constitutive model " + m_solidMaterialName + " not found" );
  m_solidMaterialFullIndex = solid->getIndexInParent();

}

void SolidMechanicsLagrangianFEM::updateIntrinsicNodalData( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();

  ElementRegionManager const * const elementRegionManager = mesh->getElemManager();
  ConstitutiveManager const * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<real64> & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  mass = 0.0;

  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  FiniteElementDiscretization const *
  feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > rho =
    elementRegionManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("density", constitutiveManager);

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(er);

    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      for( localIndex k=0 ; k < elemsToNodes.size(0) ; ++k )
      {

        // TODO this integration needs to be be carried out properly.
        real64 elemMass = 0;
        for( localIndex q=0 ; q<fe->n_quadrature_points() ; ++q )
        {
          elemMass += rho[er][esr][m_solidMaterialFullIndex][k][q] * detJ[k][q];
        }
        for( localIndex a=0 ; a< elemsToNodes.size(1) ; ++a )
        {
          mass[elemsToNodes[k][a]] += elemMass/elemsToNodes.size(1);
        }

        for( localIndex a=0 ; a<elementSubRegion->numNodesPerElement() ; ++a )
        {
          if( nodes->GhostRank()[elemsToNodes[k][a]] >= -1 )
          {
            m_sendOrReceiveNodes.insert( elemsToNodes[k][a] );
          }
          else
          {
            m_nonSendOrReceiveNodes.insert( elemsToNodes[k][a] );
          }
        }
      }
    });
  }
}

void SolidMechanicsLagrangianFEM::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const problemManager )
{
  DomainPartition * domain = problemManager->GetGroup<DomainPartition>(keys::domain);
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  NodeManager * const nodes = mesh->getNodeManager();
  FaceManager * const faceManager = mesh->getFaceManager();


  ElementRegionManager * elementRegionManager = mesh->getElemManager();
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  arrayView1d<real64> & mass = nodes->getReference<array1d<real64>>(keys::Mass);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > rho =
    elementRegionManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("density", constitutiveManager);


  NumericalMethodsManager const *
  numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteElementDiscretizationManager const *
  feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  FiniteElementDiscretization const *
  feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

  m_elemsAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );
  m_elemsNotAttachedToSendOrReceiveNodes.resize( elementRegionManager->numRegions() );

  for( localIndex er=0 ; er<elementRegionManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elementRegionManager->GetRegion(er);
    m_elemsAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );
    m_elemsNotAttachedToSendOrReceiveNodes[er].resize( elemRegion->numSubRegions() );

    elemRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      m_elemsAttachedToSendOrReceiveNodes[er][esr].setUserCallBack(
        "SolidMechanicsLagrangianFEM::m_elemsAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      m_elemsNotAttachedToSendOrReceiveNodes[er][esr].setUserCallBack(
        "SolidMechanicsLagrangianFEM::m_elemsNotAttachedToSendOrReceiveNodes["
        + std::to_string(er) + "][" + std::to_string(esr) + "]" );

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);
      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      for( localIndex k=0 ; k < elemsToNodes.size(0) ; ++k )
      {

        // TODO this integration needs to be be carried out properly.
        real64 elemMass = 0;
        for( localIndex q=0 ; q<fe->n_quadrature_points() ; ++q )
        {
          elemMass += rho[er][esr][m_solidMaterialFullIndex][k][q] * detJ[k][q];
        }
        for( localIndex a=0 ; a< elemsToNodes.size(1) ; ++a )
        {
          mass[elemsToNodes[k][a]] += elemMass/elemsToNodes.size(1);
        }

        bool isAttachedToGhostNode = false;
        for( localIndex a=0 ; a<elementSubRegion->numNodesPerElement() ; ++a )
        {
          if( nodes->GhostRank()[elemsToNodes[k][a]] >= -1 )
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
          m_elemsAttachedToSendOrReceiveNodes[er][esr].insert(k);
        }
        else
        {
          m_elemsNotAttachedToSendOrReceiveNodes[er][esr].insert(k);
        }
      }
    });
  }

#if HAVE_TRIBOLCOUPLING
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back("Velocity");

  CommunicationTools::SynchronizeFields( fieldNames,
                                         mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );
#endif
}

real64 SolidMechanicsLagrangianFEM::SolverStep( real64 const& time_n,
                                             real64 const& dt,
                                             const int cycleNumber,
                                             DomainPartition * domain )
{
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =  this->getParent()->GetGroup<SolverBase>("SurfaceGen");

  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitDynamic )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, ManagedGroup::group_cast<DomainPartition*>(domain) );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic ||
           m_timeIntegrationOption == timeIntegrationOption::QuasiStatic )
  {
    int const maxNumResolves = m_maxNumResolves;
    int locallyFractured = 0;
    int globallyFractured = 0;
    for( int solveIter=0 ; solveIter<maxNumResolves ; ++solveIter )
    {
      ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain->group_cast<DomainPartition*>(),
                                        getLinearSystemRepository() );
      if( surfaceGenerator!=nullptr )
      {
        if( !( surfaceGenerator->SolverStep( time_n, dt, cycleNumber, domain ) > 0 ) )
        {
          locallyFractured = 1;
        }
        MPI_Allreduce( &locallyFractured,
                       &globallyFractured,
                       1,
                       MPI_INT,
                       MPI_MAX,
                       MPI_COMM_GEOSX);
      }
      if( globallyFractured == 0 )
      {
        break;
      }
    }
    ImplicitStepComplete( time_n, dt,  domain );
  }

  return dtReturn;
}

real64 SolidMechanicsLagrangianFEM::ExplicitStep( real64 const& time_n,
                                                  real64 const& dt,
                                                  const int cycleNumber,
                                                  DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // updateIntrinsicNodalData(domain);

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodes = mesh->getNodeManager();
  ElementRegionManager * elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);
  ConstitutiveManager * constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);

  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
  localIndex const numNodes = nodes->size();

  arrayView1d<R1Tensor> const & X = nodes->getReference<array1d<R1Tensor>>(nodes->viewKeys.referencePosition);
  arrayView1d<real64> const & mass = nodes->getReference<array1d<real64>>(keys::Mass);
  arrayView1d<R1Tensor> & vel = nodes->getReference<array1d<R1Tensor>>(keys::Velocity);

  arrayView1d<R1Tensor> & u = nodes->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);
  arrayView1d<R1Tensor> & uhat = nodes->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
  arrayView1d<R1Tensor> & acc = nodes->getReference<array1d<R1Tensor>>(keys::Acceleration);

  array1d<NeighborCommunicator> & neighbors = domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors );
  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back("Velocity");

  CommunicationTools::SynchronizePackSendRecvSizes( fieldNames, mesh, neighbors, m_iComm );

  fsManager->ApplyFieldValue( time_n, domain, "nodeManager", keys::Acceleration );

#if HAVE_TRIBOLCOUPLING
  // rewind for leapfrog
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, -dt/2, numNodes );
#endif


#if HAVE_TRIBOLCOUPLING
  TribolCoupling::ApplyTribolForces(domain, time_n, dt, cycleNumber) ;
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt, numNodes );
#else
  //3: v^{n+1/2} = v^{n} + a^{n} dt/2
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, vel, dt/2, numNodes );
#endif

  fsManager->ApplyFieldValue( time_n, domain, "nodeManager", keys::Velocity );

  //4. x^{n+1} = x^{n} + v^{n+{1}/{2}} dt (x is displacement)
  SolidMechanicsLagrangianFEMKernels::displacementUpdate( vel, uhat, u, dt );

  fsManager->ApplyFieldValue( time_n + dt, domain, "nodeManager", keys::TotalDisplacement,
    [&]( FieldSpecificationBase const * const bc, set<localIndex> const & targetSet )->void
    {
      integer const component = bc->GetComponent();
      for( auto const a : targetSet )
      {
        vel[a][component] = u[a][component];
      }
    },
    [&]( FieldSpecificationBase const * const bc, set<localIndex> const & targetSet )->void
    {
      integer const component = bc->GetComponent();
      for( auto const a : targetSet )
      {
        uhat[a][component] = u[a][component] - vel[a][component];
        vel[a][component]  = uhat[a][component] / dt;
      }
    }
  );

  ElementRegionManager::MaterialViewAccessor< arrayView2d<real64> > meanStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<real64>, arrayView2d<real64> >("MeanStress", constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< arrayView2d<R2SymTensor> > const devStress =
    elemManager->ConstructFullMaterialViewAccessor< array2d<R2SymTensor>, arrayView2d<R2SymTensor> >("DeviatorStress", constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase> constitutiveRelations =
    elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  //Step 5. Calculate deformation input to constitutive model and update state to
  // Q^{n+1}
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);
    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor > const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      ExplicitElementKernelLaunch( numNodesPerElement,
                                   numQuadraturePoints,
                                   constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                   this->m_elemsAttachedToSendOrReceiveNodes[er][esr],
                                   elemsToNodes,
                                   dNdX,
                                   detJ,
                                   u,
                                   vel,
                                   acc,
                                   meanStress[er][esr][m_solidMaterialFullIndex],
                                   devStress[er][esr][m_solidMaterialFullIndex],
                                   dt );

    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_sendOrReceiveNodes );

  fsManager->ApplyFieldValue( time_n, domain, "nodeManager", keys::Velocity );

  CommunicationTools::SynchronizePackSendRecv( fieldNames, mesh, neighbors, m_iComm );

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const * feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      arrayView3d< R1Tensor > const & dNdX = elementSubRegion->getReference< array3d< R1Tensor > >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();

      localIndex const numNodesPerElement = elemsToNodes.size(1);

      localIndex const numQuadraturePoints = feDiscretization->m_finiteElement->n_quadrature_points();

      ExplicitElementKernelLaunch( numNodesPerElement,
                                   numQuadraturePoints,
                                   constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                   this->m_elemsNotAttachedToSendOrReceiveNodes[er][esr],
                                   elemsToNodes,
                                   dNdX,
                                   detJ,
                                   u,
                                   vel,
                                   acc,
                                   meanStress[er][esr][m_solidMaterialFullIndex],
                                   devStress[er][esr][m_solidMaterialFullIndex],
                                   dt );
    }); //Element Region

  } //Element Manager

  // apply this over a set
  SolidMechanicsLagrangianFEMKernels::velocityUpdate( acc, mass, vel, dt / 2, m_nonSendOrReceiveNodes );

  fsManager->ApplyFieldValue( time_n, domain, "nodeManager", keys::Velocity );

  CommunicationTools::SynchronizeUnpack( mesh, neighbors, m_iComm );

  return dt;
}



void SolidMechanicsLagrangianFEM::ApplyDisplacementBC_implicit( real64 const time,
                                                                 DomainPartition & domain,
                                                                 EpetraBlockSystem & blockSystem )
{

  FieldSpecificationManager const * const fsManager = FieldSpecificationManager::get();

  fsManager->Apply( time,
                     &domain,
                     "nodeManager",
                     keys::TotalDisplacement,
                     [&]( FieldSpecificationBase const * const bc,
                     string const &,
                     set<localIndex> const & targetSet,
                     ManagedGroup * const targetGroup,
                     string const fieldName )->void
    {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationEqual>( targetSet,
                                                                 false,
                                                                 time,
                                                                 targetGroup,
                                                                 fieldName,
                                                                 viewKeyStruct::globalDofNumberString,
                                                                 3,
                                                                 &blockSystem,
                                                                 BlockIDs::displacementBlock );
  });
}


void SolidMechanicsLagrangianFEM::ApplyTractionBC( DomainPartition * const domain,
                                                    real64 const time,
                                                    systemSolverInterface::EpetraBlockSystem & blockSystem )
{
  FieldSpecificationManager * const fsManager = FieldSpecificationManager::get();
  NewFunctionManager * const functionManager = NewFunctionManager::Instance();

  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  real64_array const & faceArea  = faceManager->getReference<real64_array>("faceArea");
  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();

  arrayView1d<globalIndex> const & blockLocalDofNumber =
    nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.globalDofNumber);

  Epetra_FEVector * const rhs = blockSystem.GetResidualVector( BlockIDs::displacementBlock );

  fsManager->Apply( time,
                    domain,
                    "faceManager",
                    string("Traction"),
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName ) -> void
  {
    string const & functionName = bc->getReference<string>( FieldSpecificationBase::viewKeyStruct::functionNameString);

    globalIndex_array nodeDOF;
    real64_array nodeRHS;
    integer const component = bc->GetComponent();

    if( functionName.empty() )
    {
      for( auto kf : targetSet )
      {
        localIndex const numNodes = facesToNodes[kf].size();
        nodeDOF.resize( numNodes );
        nodeRHS.resize( numNodes );
        for( localIndex a=0 ; a<numNodes ; ++a )
        {
          nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
          nodeRHS[a] = bc->GetScale() * faceArea[kf] / numNodes;
        }
        rhs->SumIntoGlobalValues( integer_conversion<int>(numNodes), nodeDOF.data(), nodeRHS.data() );
      }
    }
    else
    {
      FunctionBase const * const function  = functionManager->GetGroup<FunctionBase>(functionName);
      assert( function!=nullptr);

        if( function->isFunctionOfTime()==2 )
        {
          real64 value = bc->GetScale() * function->Evaluate( &time );
          for( auto kf : targetSet )
          {
            localIndex const numNodes = facesToNodes[kf].size();
            nodeDOF.resize( numNodes );
            nodeRHS.resize( numNodes );
            for( localIndex a=0 ; a<numNodes ; ++a )
            {
              nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = value * faceArea[kf] / numNodes;
            }
            rhs->SumIntoGlobalValues( integer_conversion<int>(nodeDOF.size()), nodeDOF.data(), nodeRHS.data() );
          }
        }
        else
        {
          real64_array result;
          result.resize( targetSet.size() );
          function->Evaluate( faceManager, time, targetSet, result );

          for( auto kf : targetSet )
          {
            localIndex const numNodes = facesToNodes[kf].size();
            nodeDOF.resize( numNodes );
            nodeRHS.resize( numNodes );
            for( localIndex a=0 ; a<numNodes ; ++a )
            {
              nodeDOF[a] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
              nodeRHS[a] = result[kf] * faceArea[kf] / numNodes;
            }
            rhs->SumIntoGlobalValues( integer_conversion<int>(nodeDOF.size()), nodeDOF.data(), nodeRHS.data() );
          }
      }
    }
  });
}

void SolidMechanicsLagrangianFEM::ApplyChomboPressure( DomainPartition * const domain,
                                                       systemSolverInterface::EpetraBlockSystem & blockSystem )
{
  FaceManager * const faceManager = domain->getMeshBody(0)->getMeshLevel(0)->getFaceManager();
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  arrayView1d<real64 const> const & faceArea  = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceNormal  = faceManager->faceNormal();
  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();

  arrayView1d<globalIndex> const &
  blockLocalDofNumber =  nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.globalDofNumber);

  Epetra_FEVector * const rhs = blockSystem.GetResidualVector( BlockIDs::displacementBlock );
  arrayView1d<real64 const> const & facePressure = faceManager->getReference< array1d<real64> >("ChomboPressure");

  for( localIndex kf=0 ; kf<faceManager->size() ; ++kf )
  {
    globalIndex nodeDOF[20];
    real64 nodeRHS[20];

    int const numNodes = integer_conversion<int>(facesToNodes[kf].size());
    for( int a=0 ; a<numNodes ; ++a )
    {
      for( int component=0 ; component<3 ; ++component )
      {
        nodeDOF[3*a+component] = 3*blockLocalDofNumber[facesToNodes[kf][a]]+component;
        nodeRHS[3*a+component] = - facePressure[kf] * faceNormal[kf][component] * faceArea[kf] / numNodes;
      }
    }
    rhs->SumIntoGlobalValues( numNodes*3, nodeDOF, nodeRHS );
  }

}

void
SolidMechanicsLagrangianFEM::
ImplicitStepSetup( real64 const& time_n, real64 const& dt, DomainPartition * const domain,
                   systemSolverInterface::EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    arrayView1d<R1Tensor> const& v_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Velocity);
    arrayView1d<R1Tensor> const& a_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Acceleration);
    arrayView1d<R1Tensor>& vtilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.vTilde);
    arrayView1d<R1Tensor>& uhatTilde   = nodeManager->getReference<array1d<R1Tensor>>(solidMechanicsViewKeys.uhatTilde);

    arrayView1d<R1Tensor>& uhat  = nodeManager->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
    arrayView1d<R1Tensor>& disp = nodeManager->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);

    localIndex const numNodes = nodeManager->size();
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

    for( localIndex a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        vtilde[a][i] = v_n[a][i] + (1.0-newmarkGamma) * a_n[a][i] * dt;
        uhatTilde[a][i] = ( v_n[a][i] + 0.5 * ( 1.0 - 2.0*newmarkBeta ) * a_n[a][i] * dt ) *dt;
        uhat[a][i] = uhatTilde[a][i];
        disp[a][i] += uhatTilde[a][i];
      }
    }
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic  )
  {

    arrayView1d<R1Tensor>& uhat = nodeManager->getReference<array1d<R1Tensor>>(keys::IncrementalDisplacement);
    integer const useVelocityEstimateForQS = this->getReference<integer>(solidMechanicsViewKeys.useVelocityEstimateForQS);
    localIndex const numNodes = nodeManager->size();

    if( useVelocityEstimateForQS==1 )
    {
      arrayView1d<R1Tensor> const& v_n = nodeManager->getReference<array1d<R1Tensor>>(keys::Velocity);
      arrayView1d<R1Tensor>& disp = nodeManager->getReference<array1d<R1Tensor>>(keys::TotalDisplacement);

      for( localIndex a = 0 ; a < numNodes ; ++a )
      {
        for( int i=0 ; i<3 ; ++i )
        {
          uhat[a][i] = v_n[a][i] * dt;
          disp[a][i] += uhat[a][i];
        }
      }
    }
    else
    {
      for( localIndex a = 0 ; a < numNodes ; ++a )
      {
        uhat[a] = 0.0;
      }
    }
  }

  SetupSystem( domain, blockSystem );
}

void SolidMechanicsLagrangianFEM::ImplicitStepComplete( real64 const & time_n,
                                                             real64 const & dt,
                                                             DomainPartition * const domain)
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  localIndex const numNodes = nodeManager->size();

  r1_array& v_n = nodeManager->getReference<r1_array>(keys::Velocity);
  r1_array& uhat  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  if( this->m_timeIntegrationOption == timeIntegrationOption::ImplicitDynamic )
  {
    r1_array& a_n = nodeManager->getReference<r1_array>(keys::Acceleration);
    r1_array& vtilde    = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.vTilde);
    r1_array& uhatTilde = nodeManager->getReference<r1_array>(solidMechanicsViewKeys.uhatTilde);
    real64 const newmarkGamma = this->getReference<real64>(solidMechanicsViewKeys.newmarkGamma);
    real64 const newmarkBeta = this->getReference<real64>(solidMechanicsViewKeys.newmarkBeta);

    for( auto a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        //        real64 a_np1 = 4.0 / (dt*dt) * ( uhat[a][i] - uhatTilde[a][i]
        // );
        a_n[a][i] = 1.0 / ( newmarkBeta * dt*dt) * ( uhat[a][i] - uhatTilde[a][i] );
        v_n[a][i] = vtilde[a][i] + newmarkGamma * a_n[a][i] * dt;
      }
    }
  }
  else if( this->m_timeIntegrationOption == timeIntegrationOption::QuasiStatic && dt > 0.0)
  {
    for( auto a = 0 ; a < numNodes ; ++a )
    {
      for( int i=0 ; i<3 ; ++i )
      {
        v_n[a][i] = uhat[a][i] / dt;
      }
    }
  }
}

void SolidMechanicsLagrangianFEM::SetNumRowsAndTrilinosIndices( ManagedGroup * const nodeManager,
                                                                 localIndex & numLocalRows,
                                                                 globalIndex & numGlobalRows,
                                                                 localIndex_array& localIndices,
                                                                 localIndex offset )
{
  int n_mpi_processes;
  MPI_Comm_size( MPI_COMM_GEOSX, &n_mpi_processes );

  int this_mpi_process = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &this_mpi_process );

  std::vector<int> gather(n_mpi_processes);

  int intNumLocalRows = static_cast<int>(numLocalRows);
  m_linearSolverWrapper.m_epetraComm.GatherAll( &intNumLocalRows,
                                                &gather.front(),
                                                1 );
  numLocalRows = intNumLocalRows;

  localIndex first_local_row = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<n_mpi_processes ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<this_mpi_process)
      first_local_row += gather[p];
  }

  // create trilinos dof indexing

  globalIndex_array& trilinos_index = nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.globalDofNumber);
  integer_array& is_ghost       = nodeManager->getReference<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString);

  trilinos_index = -1;

  integer local_count = 0;
  for(integer r=0 ; r<trilinos_index.size() ; ++r )
  {
    if(is_ghost[r] < 0)
    {
      trilinos_index[r] = first_local_row+local_count+offset;
      local_count++;
    }
    else
    {
      trilinos_index[r] = -INT_MAX;
    }
  }

  assert(local_count == numLocalRows );

//  partition.SynchronizeFields(m_syncedFields, CommRegistry::lagrangeSolver02);

}


void SolidMechanicsLagrangianFEM :: SetupSystem ( DomainPartition * const domain,
                                                   EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  localIndex dim = 3;
  localIndex n_ghost_rows  = nodeManager->GetNumberOfGhosts();
  localIndex n_local_rows  = nodeManager->size()-n_ghost_rows;
  globalIndex n_global_rows = 0;

  localIndex_array displacementIndices;
  SetNumRowsAndTrilinosIndices( nodeManager, n_local_rows, n_global_rows, displacementIndices, 0 );

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(viewKeyStruct::globalDofNumberString);

  CommunicationTools::SynchronizeFields( fieldNames, mesh,
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // create epetra map
  Epetra_Map * const rowMap = blockSystem->SetRowMap( BlockIDs::displacementBlock,
                                                      std::make_unique<Epetra_Map>( dim*n_global_rows, 
                                                                                    dim*n_local_rows,
                                                                                    0,
                                                                                    m_linearSolverWrapper.m_epetraComm ) );

  Epetra_FECrsGraph * const sparsity = blockSystem->SetSparsity( BlockIDs::displacementBlock,
                                                                 BlockIDs::displacementBlock,
                                                                 std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );

  SetSparsityPattern( domain, sparsity );

  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  blockSystem->SetMatrix( BlockIDs::displacementBlock,
                          BlockIDs::displacementBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  blockSystem->SetSolutionVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  blockSystem->SetResidualVector( BlockIDs::displacementBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );
}

void SolidMechanicsLagrangianFEM::SetSparsityPattern( DomainPartition const * const domain,
                                                       Epetra_FECrsGraph * const sparsity )
{
  int dim=3;
  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup const * const nodeManager = mesh->getNodeManager();

  arrayView1d<globalIndex> const & trilinos_index = nodeManager->getReference<array1d<globalIndex>>(solidMechanicsViewKeys.globalDofNumber);

  ElementRegionManager const * const elemManager = mesh->getElemManager();

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elementRegion = elemManager->GetRegion(er);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr, CellElementSubRegion const * const elementSubRegion )
    {
      localIndex const numElems = elementSubRegion->size();
      arrayView2d<localIndex> const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      globalIndex_array elementLocalDofIndex(dim * numNodesPerElement);

      for( localIndex k=0 ; k<numElems ; ++k )
      {
        for( localIndex a=0 ; a<numNodesPerElement ; ++a )
        {
          for(localIndex i=0 ; i<numNodesPerElement ; ++i)
          {
            for( int d=0 ; d<dim ; ++d )
            {
              elementLocalDofIndex[i * dim + d] = dim * trilinos_index[elemsToNodes[k][i]] + d;
            }
          }

          sparsity->InsertGlobalIndices(static_cast<int>(elementLocalDofIndex.size()),
                                        elementLocalDofIndex.data(),
                                        static_cast<int>(elementLocalDofIndex.size()),
                                        elementLocalDofIndex.data());
        }
      }
    });
  }
}


void SolidMechanicsLagrangianFEM::AssembleSystem ( DomainPartition * const  domain,
                                                EpetraBlockSystem * const blockSystem,
                                                real64 const time_n,
                                                real64 const dt )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ManagedGroup * const nodeManager = mesh->getNodeManager();
  ConstitutiveManager  * const constitutiveManager = domain->GetGroup<ConstitutiveManager >(keys::ConstitutiveManager);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  NumericalMethodsManager const * numericalMethodManager = domain->getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);
  FiniteElementDiscretizationManager const * feDiscretizationManager = numericalMethodManager->GetGroup<FiniteElementDiscretizationManager>(keys::finiteElementDiscretizations);

  ElementRegionManager::MaterialViewAccessor<real64> const biotCoefficient =
    elemManager->ConstructFullMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const fluidPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("pressure");

  ElementRegionManager::ElementViewAccessor<arrayView1d<real64>> const dPres =
    elemManager->ConstructViewAccessor<array1d<real64>, arrayView1d<real64>>("deltaPressure");

  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                              BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );
  Epetra_FEVector * const solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  matrix->Scale(0.0);
  rhs->Scale(0.0);

  r1_array const& disp = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  r1_array const& uhat = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);
  r1_array const& vel  = nodeManager->getReference<r1_array>(keys::Velocity);

  r1_array const uhattilde;
  r1_array const vtilde;

  globalIndex_array const &
  trilinos_index = nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.globalDofNumber);


  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->ConstructFullConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< real64 > const
  density = elemManager->ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                                  constitutiveManager );

  // begin region loop
  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion * const elementRegion = elemManager->GetRegion(er);

    FiniteElementDiscretization const *
    feDiscretization = feDiscretizationManager->GetGroup<FiniteElementDiscretization>(m_discretizationName);

    elementRegion->forElementSubRegionsIndex<CellElementSubRegion>([&]( localIndex const esr,
                                                                        CellElementSubRegion const * const elementSubRegion )
    {
      array3d<R1Tensor> const &
      dNdX = elementSubRegion->getReference< array3d<R1Tensor> >(keys::dNdX);

      arrayView2d<real64> const & detJ = elementSubRegion->getReference< array2d<real64> >(keys::detJ);

      arrayView2d< localIndex > const & elemsToNodes = elementSubRegion->nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size(1);

      std::unique_ptr<FiniteElementBase>
      fe = feDiscretization->getFiniteElement( elementSubRegion->GetElementTypeString() );

      // space for element matrix and rhs
      m_maxForce = ImplicitElementKernelLaunch( numNodesPerElement,
                                                fe->n_quadrature_points(),
                                                constitutiveRelations[er][esr][m_solidMaterialFullIndex],
                                                elementSubRegion->size(),
                                                dt,
                                                dNdX,
                                                detJ,
                                                fe.get(),
                                                elementSubRegion->m_ghostRank,
                                                elemsToNodes,
                                                trilinos_index,
                                                disp,
                                                uhat,
                                                vtilde,
                                                uhattilde,
                                                density[er][esr],
                                                fluidPres[er][esr],
                                                dPres[er][esr],
                                                biotCoefficient[er][esr],
                                                m_timeIntegrationOption,
                                                this->m_stiffnessDamping,
                                                this->m_massDamping,
                                                this->m_newmarkBeta,
                                                this->m_newmarkGamma,
                                                matrix,
                                                rhs );

    });
  }


  // Global assemble
  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }
}

void
SolidMechanicsLagrangianFEM::
ApplyBoundaryConditions( DomainPartition * const domain,
                         systemSolverInterface::EpetraBlockSystem * const blockSystem,
                         real64 const time_n,
                         real64 const dt )
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  FaceManager * const faceManager = mesh->getFaceManager();
  NodeManager * const nodeManager = mesh->getNodeManager();
  ElementRegionManager * const elemManager = mesh->getElemManager();

  FieldSpecificationManager * fsManager = FieldSpecificationManager::get();
//  fsManager->ApplyBoundaryCondition( this, &SolidMechanics_LagrangianFEM::ForceBC,
//                                     nodeManager, keys::Force, time_n + dt, *blockSystem );

  fsManager->Apply( time_n+dt,
                    domain,
                    "nodeManager",
                    keys::Force,
                    [&]( FieldSpecificationBase const * const bc,
                    string const &,
                    set<localIndex> const & targetSet,
                    ManagedGroup * const targetGroup,
                    string const fieldName )->void
  {
    bc->ApplyBoundaryConditionToSystem<FieldSpecificationAdd>( targetSet,
                                                               false,
                                                               time_n+dt,
                                                               targetGroup,
                                                               keys::TotalDisplacement, // TODO fix use of dummy name for
                                                               viewKeyStruct::globalDofNumberString,
                                                               3,
                                                               blockSystem,
                                                               BlockIDs::displacementBlock );
  });

  ApplyTractionBC( domain,
                   time_n+dt,
                   *blockSystem );

  ApplyDisplacementBC_implicit( time_n + dt, *domain, *blockSystem );

  if( faceManager->hasView("ChomboPressure") )
  {
    fsManager->ApplyFieldValue( time_n, domain, "faceManager", "ChomboPressure" );
    ApplyChomboPressure( domain, *blockSystem );
  }

  {
  arrayView1d<real64 const>   const & faceArea   = faceManager->faceArea();
  arrayView1d<R1Tensor const> const & faceNormal = faceManager->faceNormal();
  array1d<localIndex_array> const & facesToNodes = faceManager->nodeList();

  arrayView1d<globalIndex> const &
  blockLocalDofNumber =  nodeManager->getReference<globalIndex_array>(solidMechanicsViewKeys.globalDofNumber);
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  elemManager->forElementSubRegions<FaceElementSubRegion>([&]( FaceElementSubRegion * const subRegion )->void
  {
    if( subRegion->hasView("pressure") )
    {
      arrayView1d<real64 const> const & fluidPressure = subRegion->getReference<array1d<real64> >("pressure");
      arrayView1d<real64 const> const & deltaFluidPressure = subRegion->getReference<array1d<real64> >("deltaPressure");

      FaceElementSubRegion::FaceMapType const & faceMap = subRegion->faceList();

      forall_in_range<serialPolicy>( 0,
                                   subRegion->size(),
                                   GEOSX_LAMBDA ( localIndex const kfe )
      {

        R1Tensor Nbar = faceNormal[faceMap[kfe][0]];
        Nbar -= faceNormal[faceMap[kfe][1]];
        Nbar.Normalize();
        std::cout<<"kfe, faceIndices = "<<kfe<<", ("<<faceMap[kfe][0]<<", "<<faceMap[kfe][1]<<")"<<std::endl;

        std::cout<<"Nbar = "<<Nbar<<std::endl;

        globalIndex nodeDOF[20];
        real64 nodeRHS[20];

        localIndex const numNodes = facesToNodes[faceMap[kfe][0]].size();
        real64 nodalForce = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * faceArea[faceMap[kfe][0]] / numNodes;
        std::cout<<"nodalForce = "<<nodalForce<<std::endl;

        for( localIndex kf=0 ; kf<2 ; ++kf )
        {
          localIndex const faceIndex = faceMap[kfe][kf];
          std::cout<<"faceIndex = "<<faceIndex<<std::endl;

          for( localIndex a=0 ; a<numNodes ; ++a )
          {

            std::cout<<facesToNodes[faceIndex][a]<<", ";
            for( int component=0 ; component<3 ; ++component )
            {
              nodeDOF[3*a+component] = 3*blockLocalDofNumber[facesToNodes[faceIndex][a]]+component;
              nodeRHS[3*a+component] = - nodalForce * pow(-1,kf) * Nbar[component];
            }
          }
          std::cout<<std::endl;

          rhs->SumIntoGlobalValues( integer_conversion<int>(numNodes*3), nodeDOF, nodeRHS );
        }
      });
    }
  });
  }

  Epetra_FECrsMatrix * const matrix = blockSystem->GetMatrix( BlockIDs::displacementBlock,
                                                                    BlockIDs::displacementBlock );
  Epetra_FEVector * const rhs = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  if( verboseLevel() >= 2 )
  {
    matrix->Print(std::cout);
    rhs->Print(std::cout);
  }

  matrix->GlobalAssemble(true);
  rhs->GlobalAssemble();


}

real64
SolidMechanicsLagrangianFEM::
CalculateResidualNorm(systemSolverInterface::EpetraBlockSystem const *const blockSystem, DomainPartition *const domain)
{

  Epetra_FEVector const * const
  residual = blockSystem->GetResidualVector( BlockIDs::displacementBlock );

  real64 localResidual[2] = {0.0, this->m_maxForce};
//  residual->Norm2(&scalarResidual);

  real64 * residualData = nullptr;
  int length;
  residual->ExtractView(&residualData,&length);
  for( localIndex i=0 ; i<length ; ++i )
  {
    localResidual[0] += residualData[i]*residualData[i];
  }


  real64 globalResidualNorm[2] = {0,0};
//  MPI_Allreduce (&localResidual,&globalResidualNorm,1,MPI_DOUBLE,MPI_SUM ,MPI_COMM_GEOSX);


  int rank, size;
  MPI_Comm_rank(MPI_COMM_GEOSX, &rank);
  MPI_Comm_size(MPI_COMM_GEOSX, &size);
  array1d<real64> globalValues( size * 2 );
  globalValues = 0;
  MPI_Gather( localResidual,
              2,
              MPI_DOUBLE,
              globalValues.data(),
              2,
              MPI_DOUBLE,
              0,
              MPI_COMM_GEOSX );

  if( rank==0 )
  {
    for( int r=0 ; r<size ; ++r )
    {
      globalResidualNorm[0] += globalValues[r*2];

      if( globalResidualNorm[1] < globalValues[r*2+1] )
      {
        globalResidualNorm[1] = globalValues[r*2+1];
      }
    }
  }

  MPI_Bcast( globalResidualNorm, 2, MPI_DOUBLE, 0, MPI_COMM_GEOSX );



  return sqrt(globalResidualNorm[0])/(globalResidualNorm[1]+1);

}



void SolidMechanicsLagrangianFEM::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                        real64 const scalingFactor,
                                                        DomainPartition * const domain )
{
  NodeManager * const nodeManager = domain->getMeshBody(0)->getMeshLevel(0)->getNodeManager();

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::displacementBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  int solutionLength;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&solutionLength);


  globalIndex_array const& trilinos_index = nodeManager->getReference< globalIndex_array >(solidMechanicsViewKeys.globalDofNumber);

  r1_array& X        = nodeManager->getReference<r1_array>(nodeManager->viewKeys.referencePosition);
  r1_array& disp     = nodeManager->getReference<r1_array>(keys::TotalDisplacement);
  r1_array& incdisp  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  localIndex const numNodes = nodeManager->size();
  realT maxpos = 0.0;
  realT maxdisp = 0.0;

  integer const dim = 3;

  for(integer r=0 ; r<numNodes ; ++r)
  {
    {
      for( int d=0 ; d<dim ; ++d )
      {
        int lid = rowMap->LID( static_cast<int>(dim*trilinos_index[r]) + d );

        if( lid >=0 )
        {
          incdisp[r][d] -= scalingFactor*local_solution[lid];
          disp[r][d] -= scalingFactor*local_solution[lid];
          maxpos = std::max( maxpos, fabs(X[r][d]+disp[r][d]) );
          maxdisp = std::max( maxdisp, fabs(disp[r][d]) );
        }
      }
    }
  }
//  m_maxDofVal = maxpos;
//  std::cout<<"Maximum DeltaDisplacement, Position = "<<maxinc<<",
// "<<maxpos<<", "<<maxinc/maxpos<<std::endl;

  std::map<string, string_array > fieldNames;
  fieldNames["node"].push_back(keys::IncrementalDisplacement);
  fieldNames["node"].push_back(keys::TotalDisplacement);

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain->getMeshBody(0)->getMeshLevel(0),
                                         domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );



}

void SolidMechanicsLagrangianFEM::SolveSystem( EpetraBlockSystem * const blockSystem,
                                        SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::displacementBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::displacementBlock );
//  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                 params,
                                                 BlockIDs::displacementBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

void SolidMechanicsLagrangianFEM::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  NodeManager * const nodeManager = mesh->getNodeManager();

  r1_array& incdisp  = nodeManager->getReference<r1_array>(keys::IncrementalDisplacement);

  // TODO need to finish this rewind
  forall_in_range(0, nodeManager->size(), GEOSX_LAMBDA (localIndex a) mutable
  {
    incdisp[a] = 0.0;
  });
}




REGISTER_CATALOG_ENTRY( SolverBase, SolidMechanicsLagrangianFEM, string const &, dataRepository::ManagedGroup * const )
}
