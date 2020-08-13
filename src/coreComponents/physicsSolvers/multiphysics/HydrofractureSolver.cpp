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
 * @file HydrofractureSolver.cpp
 *
 */

#include "HydrofractureSolver.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const std::string & name,
                                          Group * const parent ) :
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString( "FIM" ),
  m_couplingTypeOption(),
  m_solidSolver( nullptr ),
  m_flowSolver( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString,
                   &m_couplingTypeOptionString )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Coupling option: (FIM, SIM_FixedStress)" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString, &m_maxNumResolves )
    ->setApplyDefaultValue( 10 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Value to indicate how many resolves may be executed to perform surface "
      "generation after the execution of flow and mechanics solver. " );

  m_numResolves[0] = 0;
}

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
void
HydrofractureSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel =
      Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    elemManager->forElementRegions< FaceElementRegion >(
      [&]( FaceElementRegion * const region ) {
        region->forElementSubRegions< FaceElementSubRegion >(
          [&]( FaceElementSubRegion * const subRegion ) {
            subRegion
              ->registerWrapper< array1d< real64 > >(
                viewKeyStruct::separationCoeff0String )
              ->setRestartFlags( RestartFlags::NO_WRITE );
            subRegion
              ->registerWrapper< array1d< real64 > >(
                viewKeyStruct::apertureAtFailureString )
              ->setApplyDefaultValue( -1.0 )
              ->setPlotLevel( PlotLevel::LEVEL_0 );

            subRegion
              ->registerWrapper< array1d< real64 > >(
                FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString )
              ->setRestartFlags( RestartFlags::NO_WRITE );
          } );
      } );
  }
}
#endif

void
HydrofractureSolver::ImplicitStepSetup( real64 const & time_n,
                                        real64 const & dt,
                                        DomainPartition & domain )
{
  UpdateDeformationForCoupling( domain );
  m_solidSolver->ImplicitStepSetup( time_n, dt, domain );
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  mesh.getElemManager()->forElementRegions< FaceElementRegion >(
    [&]( FaceElementRegion & faceElemRegion ) {
      faceElemRegion.forElementSubRegions< FaceElementSubRegion >(
        [&]( FaceElementSubRegion & subRegion ) {
          arrayView1d< real64 > const & separationCoeff0 =
            subRegion.getReference< array1d< real64 > >(
              viewKeyStruct::separationCoeff0String );
          arrayView1d< real64 const > const & separationCoeff =
            subRegion.getSeparationCoefficient();
          for( localIndex k = 0; k < separationCoeff0.size(); ++k )
          {
            separationCoeff0[k] = separationCoeff[k];
          }
        } );
    } );
#endif
}

void
HydrofractureSolver::ImplicitStepComplete( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition & domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

void
HydrofractureSolver::PostProcessInput()
{
  string ctOption =
    this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  if( ctOption == "SIM_FixedStress" )
  {
    this->m_couplingTypeOption = couplingTypeOption::SIM_FixedStress;
  }
  else if( ctOption == "FIM" )
  {
    this->m_couplingTypeOption = couplingTypeOption::FIM;
  }
  else
  {
    GEOSX_ERROR( "invalid coupling type option: " + ctOption );
  }

  m_solidSolver =
    this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr,
                  this->getName()
                    << ": invalid solid solver name: " << m_solidSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr,
                  this->getName()
                    << ": invalid flow solver name: " << m_flowSolverName );
}

void
HydrofractureSolver::InitializePostInitialConditions_PreSubGroups(
  Group * const GEOSX_UNUSED_PARAM( problemManager ) )
{}

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
}

void
HydrofractureSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_solidSolver->ResetStateToBeginningOfStep( domain );
}

real64
HydrofractureSolver::SolverStep( real64 const & time_n,
                                 real64 const & dt,
                                 int const cycleNumber,
                                 DomainPartition & domain )
{
  real64 dtReturn = dt;

  SolverBase * const surfaceGenerator =
    this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );

  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    dtReturn = SplitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FIM )
  {
    ImplicitStepSetup( time_n, dt, domain );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;
    for( solveIter = 0; solveIter < maxIter; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

      if( solveIter > 0 )
      {
        m_solidSolver->ResetStressToBeginningOfStep( domain );
      }

      // currently the only method is implicit time integration
      dtReturn = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

      //      m_solidSolver->updateStress( domain );

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
        std::map< string, string_array > fieldNames;
        fieldNames["node"].emplace_back( keys::IncrementalDisplacement );
        fieldNames["node"].emplace_back( keys::TotalDisplacement );
        fieldNames["elems"].emplace_back(
          string( FlowSolverBase::viewKeyStruct::pressureString ) );
        fieldNames["elems"].emplace_back( "elementAperture" );

        CommunicationTools::SynchronizeFields(
          fieldNames,
          domain.getMeshBody( 0 )->getMeshLevel( 0 ),
          domain.getNeighbors() );

        this->UpdateDeformationForCoupling( domain );

        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0(
            "++ Fracture propagation. Re-entering Newton Solve." );
        }
        m_flowSolver->ResetViews( *domain.getMeshBody( 0 )->getMeshLevel( 0 ) );
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    ImplicitStepComplete( time_n, dtReturn, domain );
    m_numResolves[1] = solveIter;
  }

  return dtReturn;
}

void
HydrofractureSolver::UpdateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel * const meshLevel = domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  NodeManager const * const nodeManager = meshLevel->getNodeManager();
  FaceManager * const faceManager = meshLevel->getFaceManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
    nodeManager->totalDisplacement();
  arrayView2d< real64 const > const & faceNormal = faceManager->faceNormal();
  // arrayView1d<real64 const> const & faceArea = faceManager->faceArea();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    faceManager->nodeList().toViewConst();

  ConstitutiveManager const * const constitutiveManager =
    domain.getConstitutiveManager();

  ContactRelationBase const * const contactRelation =
    constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  elemManager->forElementSubRegions< FaceElementSubRegion >(
    [&]( FaceElementSubRegion & subRegion ) {
      arrayView1d< real64 > const & aperture = subRegion.getElementAperture();
      arrayView1d< real64 > const & effectiveAperture =
        subRegion.getReference< array1d< real64 > >(
          FlowSolverBase::viewKeyStruct::effectiveApertureString );
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 > const & deltaVolume =
        subRegion.getReference< array1d< real64 > >(
          FlowSolverBase::viewKeyStruct::deltaVolumeString );
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
      arrayView1d< real64 const > const & apertureF =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::apertureAtFailureString );

      arrayView1d< real64 > const & separationCoeff =
        subRegion.getSeparationCoefficient();

      arrayView1d< real64 > const & dSeparationCoeff_dAper =
        subRegion.getReference< array1d< real64 > >(
          FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
      arrayView1d< real64 const > const & separationCoeff0 =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::separationCoeff0String );
#endif

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe ) {
        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const kf1 = elemsToFaces[kfe][1];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
        real64 temp[3] = { 0 };
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          LvArray::tensorOps::add< 3 >( temp, u[faceToNodeMap( kf0, a )] );
          LvArray::tensorOps::subtract< 3 >( temp, u[faceToNodeMap( kf1, a )] );
        }

        // TODO this needs a proper contact based strategy for aperture
        aperture[kfe] =
          -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[kf0] ) / numNodesPerFace;

        effectiveAperture[kfe] =
          contactRelation->effectiveAperture( aperture[kfe] );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
        real64 const s = aperture[kfe] / apertureF[kfe];
        if( separationCoeff0[kfe] < 1.0 && s > separationCoeff0[kfe] )
        {
          if( s >= 1.0 )
          {
            separationCoeff[kfe] = 1.0;
            dSeparationCoeff_dAper[kfe] = 0.0;
          }
          else
          {
            separationCoeff[kfe] = s;
            dSeparationCoeff_dAper[kfe] = 1.0 / apertureF[kfe];
          }
        }
#endif
        deltaVolume[kfe] = effectiveAperture[kfe] * area[kfe] - volume[kfe];
      } );

      //#if defined(USE_CUDA)
      //    deltaVolume.move( LvArray::MemorySpace::GPU );
      //    aperture.move( LvArray::MemorySpace::GPU );
      //    effectiveAperture.move( LvArray::MemorySpace::GPU );
      //#endif
    } );
}

real64
HydrofractureSolver::SplitOperatorStep(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & dt,
  integer const GEOSX_UNUSED_PARAM( cycleNumber ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "Not implemented" );
  real64 dtReturn = dt;
  //  real64 dtReturnTemporary = dtReturn;
  //
  //  m_flowSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  //  m_solidSolver->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  //  this->ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );
  //
  //
  //
  //  fluidSolver.ImplicitStepSetup( time_n, dt, domain,
  //                                 fluidSolver.getDofManager(),
  //                                 fluidSolver.getSystemMatrix(),
  //                                 fluidSolver.getSystemRhs(),
  //                                 fluidSolver.getSystemSolution() );
  //
  //  solidSolver.ImplicitStepSetup( time_n, dt, domain,
  //                                 solidSolver.getDofManager(),
  //                                 solidSolver.getSystemMatrix(),
  //                                 solidSolver.getSystemRhs(),
  //                                 solidSolver.getSystemSolution() );
  //
  //  this->UpdateDeformationForCoupling(domain);
  //
  //  int iter = 0;
  //  while (iter < solverParams->maxIterNewton() )
  //  {
  //    if (iter == 0)
  //    {
  //      // reset the states of all child solvers if any of them has been reset
  //      m_flowSolver->ResetStateToBeginningOfStep( domain );
  //      m_solidSolver->ResetStateToBeginningOfStep( domain );
  //      ResetStateToBeginningOfStep( domain );
  //    }
  //    LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );
  //
  //    // call assemble to fill the matrix and the rhs
  //    m_flowSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
  //
  //    // apply boundary conditions to system
  //    m_flowSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
  //
  //    // call the default linear solver on the system
  //    m_flowSolver->SolveSystem( getLinearSystemRepository(),
  //                 getSystemSolverParameters() );
  //
  //    // apply the system solution to the fields/variables
  //    m_flowSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
  //
  //    if (dtReturnTemporary < dtReturn)
  //    {
  //      iter = 0;
  //      dtReturn = dtReturnTemporary;
  //      continue;
  //    }
  //
  ////    if (m_fluidSolver->getSystemSolverParameters()->numNewtonIterations() == 0 && iter > 0 && getLogLevel() >= 1)
  ////    {
  ////      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
  ////      break;
  ////    }
  //
  //    if (getLogLevel() >= 1)
  //    {
  //      GEOSX_LOG_RANK_0( "\tIteration: " << iter+1  << ", MechanicsSolver: " );
  //    }
  //
  //    // call assemble to fill the matrix and the rhs
  //    m_solidSolver->AssembleSystem( domain, getLinearSystemRepository(), time_n+dt, dt );
  //
  //
  //    ApplyFractureFluidCoupling( domain, *getLinearSystemRepository() );
  //
  //    // apply boundary conditions to system
  //    m_solidSolver->ApplyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
  //
  //    // call the default linear solver on the system
  //    m_solidSolver->SolveSystem( getLinearSystemRepository(),
  //                 getSystemSolverParameters() );
  //
  //    // apply the system solution to the fields/variables
  //    m_solidSolver->ApplySystemSolution( getLinearSystemRepository(), 1.0, domain );
  //
  //    if( m_flowSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
  //        m_solidSolver->CalculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
  //    {
  //      GEOSX_LOG_RANK_0( "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
  //      break;
  //    }
  //
  //    if (dtReturnTemporary < dtReturn)
  //    {
  //      iter = 0;
  //      dtReturn = dtReturnTemporary;
  //      continue;
  //    }
  ////    if (m_solidSolver->getSystemSolverParameters()->numNewtonIterations() > 0)
  //    {
  //      this->UpdateDeformationForCoupling(domain);
  ////      m_fluidSolver->UpdateState(domain);
  //    }
  //    ++iter;
  //  }
  //
  //  this->ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

real64
HydrofractureSolver::ExplicitStep( real64 const & time_n,
                                   real64 const & dt,
                                   const int cycleNumber,
                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ExplicitStep( time_n, dt, cycleNumber, domain );
  m_flowSolver->SolverStep( time_n, dt, cycleNumber, domain );

  return dt;
}

void
HydrofractureSolver::SetupDofs( DomainPartition const & domain,
                                DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  // restrict coupling to fracture regions only (as done originally in SetupSystem)
  ElementRegionManager const & elemManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< FaceElementRegion >(
    [&]( FaceElementRegion const & elementRegion ) {
      fractureRegions.emplace_back( elementRegion.getName() );
    } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem,
                          fractureRegions );
}

void
HydrofractureSolver::SetupSystem(
  DomainPartition & domain,
  DofManager & dofManager,
  CRSMatrix< real64, globalIndex > & GEOSX_UNUSED_PARAM( localMatrix ),
  array1d< real64 > & GEOSX_UNUSED_PARAM( localRhs ),
  array1d< real64 > & GEOSX_UNUSED_PARAM( localSolution ),
  bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  m_flowSolver->ResetViews( mesh );

  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getLocalMatrix(),
                              m_solidSolver->getLocalRhs(),
                              m_solidSolver->getLocalSolution() );

  m_flowSolver->SetupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getLocalRhs(),
                             m_flowSolver->getLocalSolution(),
                             setSparsity );

  // setup coupled DofManager
  m_dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );

  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.
  m_matrix01.createWithLocalSize( m_solidSolver->getLocalMatrix().numRows(),
                                  m_flowSolver->getLocalMatrix().numRows(),
                                  9,
                                  MPI_COMM_GEOSX );
  m_matrix10.createWithLocalSize( m_flowSolver->getLocalMatrix().numRows(),
                                  m_solidSolver->getLocalMatrix().numRows(),
                                  24,
                                  MPI_COMM_GEOSX );

#if 0
  dofManager.setSparsityPattern( m_matrix01, keys::TotalDisplacement, FlowSolverBase::viewKeyStruct::pressureString );
  dofManager.setSparsityPattern( m_matrix10, FlowSolverBase::viewKeyStruct::pressureString, keys::TotalDisplacement );
#else

  NodeManager const & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const presDofKey = m_flowSolver->getDofManager().getKey(
    FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey =
    m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< globalIndex_array >( dispDofKey );

  m_matrix01.open();
  m_matrix10.open();

  elemManager.forElementSubRegions< FaceElementSubRegion >(
    [&]( FaceElementSubRegion const & elementSubRegion ) {
      localIndex const numElems = elementSubRegion.size();
      array1d< array1d< localIndex > > const & elemsToNodes =
        elementSubRegion.nodeList();
      arrayView1d< globalIndex const > const & faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k = 0; k < numElems; ++k )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[k];
        localIndex const numNodesPerElement = elemsToNodes[k].size();
        array1d< globalIndex > activeDisplacementDOF( 3 * numNodesPerElement );
        array1d< real64 > values( 3 * numNodesPerElement );
        values.setValues< serialPolicy >( 1 );

        for( localIndex a = 0; a < numNodesPerElement; ++a )
        {
          for( int d = 0; d < 3; ++d )
          {
            activeDisplacementDOF[a * 3 + d] =
              dispDofNumber[elemsToNodes[k][a]] + d;
          }
        }

        m_matrix01.insert( activeDisplacementDOF.data(),
                           &activeFlowDOF,
                           values.data(),
                           activeDisplacementDOF.size(),
                           1 );

        m_matrix10.insert( &activeFlowDOF,
                           activeDisplacementDOF.data(),
                           values.data(),
                           1,
                           activeDisplacementDOF.size() );
      }
    } );

  NumericalMethodsManager const & numericalMethodManager =
    domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager =
    numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox =
    fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >(
    mesh,
    [&]( FaceElementStencil const & stencil ) {
      for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
      {
        localIndex const numFluxElems = stencil.stencilSize( iconn );
        typename FaceElementStencil::IndexContainerViewConstType const & seri =
          stencil.getElementRegionIndices();
        typename FaceElementStencil::IndexContainerViewConstType const & sesri =
          stencil.getElementSubRegionIndices();
        typename FaceElementStencil::IndexContainerViewConstType const & sei =
          stencil.getElementIndices();

        FaceElementSubRegion const & elementSubRegion =
          *elemManager.GetRegion( seri[iconn][0] )
             ->GetSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

        array1d< array1d< localIndex > > const & elemsToNodes =
          elementSubRegion.nodeList();

        arrayView1d< globalIndex const > const & faceElementDofNumber =
          elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

        for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
        {
          globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];

          for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
          {
            localIndex const numNodesPerElement =
              elemsToNodes[sei[iconn][k1]].size();
            array1d< globalIndex > activeDisplacementDOF( 3 * numNodesPerElement );
            array1d< real64 > values( 3 * numNodesPerElement );
            values.setValues< serialPolicy >( 1 );

            for( localIndex a = 0; a < numNodesPerElement; ++a )
            {
              for( int d = 0; d < 3; ++d )
              {
                activeDisplacementDOF[a * 3 + d] =
                  dispDofNumber[elemsToNodes[sei[iconn][k1]][a]] + d;
              }
            }

            m_matrix10.insert( &activeFlowDOF,
                               activeDisplacementDOF.data(),
                               values.data(),
                               1,
                               activeDisplacementDOF.size() );
          }
        }
      }  //);
    } );

  m_matrix01.close();
  m_matrix10.close();
#endif
}

void
HydrofractureSolver::AssembleSystem(
  real64 const time,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->getLocalMatrix().setValues< parallelDevicePolicy<> >( 0.0 );
  m_solidSolver->getLocalRhs().setValues< parallelDevicePolicy<> >( 0.0 );

  m_solidSolver->AssembleSystem(
    time,
    dt,
    domain,
    m_solidSolver->getDofManager(),
    m_solidSolver->getLocalMatrix().toViewConstSizes(),
    m_solidSolver->getLocalRhs().toView() );

  m_flowSolver->getLocalMatrix().setValues< parallelDevicePolicy<> >( 0.0 );
  m_flowSolver->getLocalRhs().setValues< parallelDevicePolicy<> >( 0.0 );

  //  CRSMatrixView< real64 const, globalIndex const > const & matrix11 = m_flowSolver->getLocalMatrix().toViewConst();
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "matrix11" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  matrix11.move( LvArray::MemorySpace::CPU );
  //  std::cout<<matrix11<<std::endl;
  //
  //  arrayView1d< real64 const > const & rhs1 = m_flowSolver->getLocalRhs().toView();
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "rhs1" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  rhs1.move( LvArray::MemorySpace::CPU );
  //  std::cout<<rhs1<<std::endl;

  m_flowSolver->ResetViews( *( domain.getMeshBody( 0 )->getMeshLevel( 0 ) ) );

  m_flowSolver->AssembleSystem( time,
                                dt,
                                domain,
                                m_flowSolver->getDofManager(),
                                m_flowSolver->getLocalMatrix().toViewConstSizes(),
                                m_flowSolver->getLocalRhs().toView() );

  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "matrix11" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  matrix11.move( LvArray::MemorySpace::CPU );
  //  std::cout<<matrix11<<std::endl;
  //
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "rhs1" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  rhs1.move( LvArray::MemorySpace::CPU );
  //  std::cout<<rhs1<<std::endl;

  m_matrix01.zero();
  AssembleForceResidualDerivativeWrtPressure( domain,
                                              &m_matrix01,
                                              m_solidSolver->getLocalRhs() );

  m_matrix10.zero();
  AssembleFluidMassResidualDerivativeWrtDisplacement( domain, &m_matrix10 );
}

void
HydrofractureSolver::ApplyBoundaryConditions(
  real64 const time,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;

  //  arrayView1d< real64 const > const & rhs0 = m_solidSolver->getLocalRhs().toView();
  //  arrayView1d< real64 const > const & rhs1 = m_flowSolver->getLocalRhs().toView();
  //  CRSMatrixView< real64 const, globalIndex const > const & matrix11 = m_flowSolver->getLocalMatrix().toViewConst();
  //
  //
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "matrix10" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  m_matrix10.print( std::cout );
  //
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "matrix11" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  matrix11.move( LvArray::MemorySpace::CPU );
  //  std::cout<<matrix11<<std::endl;
  //
  //  GEOSX_LOG_RANK_0("***********************************************************");
  //  GEOSX_LOG_RANK_0("residual0");
  //  GEOSX_LOG_RANK_0("***********************************************************");
  //  std::cout<<rhs0<<std::endl;
  //
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  GEOSX_LOG_RANK_0( "residual1" );
  //  GEOSX_LOG_RANK_0( "***********************************************************" );
  //  std::cout<<rhs1<<std::endl;

  m_solidSolver->ApplyBoundaryConditions(
    time,
    dt,
    domain,
    m_solidSolver->getDofManager(),
    m_solidSolver->getLocalMatrix().toViewConstSizes(),
    m_solidSolver->getLocalRhs().toView() );

  MeshLevel * const mesh = domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();
  string const dispDofKey =
    m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager->getReference< globalIndex_array >( dispDofKey );
  arrayView1d< integer const > const & nodeGhostRank = nodeManager->ghostRank();

  m_matrix01.open();
  fsManager.Apply( time + dt,
                   &domain,
                   "nodeManager",
                   keys::TotalDisplacement,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const,
                        string const & ) {
                     SortedArray< localIndex > localSet;
                     for( auto const & a : targetSet )
                     {
                       if( nodeGhostRank[a] < 0 )
                       {
                         localSet.insert( a );
                       }
                     }
                     bc->ZeroSystemRowsForBoundaryCondition< LAInterface >(
                       localSet.toViewConst(),
                       dispDofNumber,
                       m_matrix01 );
                   } );
  m_matrix01.close();

  m_flowSolver->ApplyBoundaryConditions(
    time,
    dt,
    domain,
    m_flowSolver->getDofManager(),
    m_flowSolver->getLocalMatrix().toViewConstSizes(),
    m_flowSolver->getLocalRhs().toView() );

  string const presDofKey = m_flowSolver->getDofManager().getKey(
    FlowSolverBase::viewKeyStruct::pressureString );

  m_matrix10.open();
  fsManager.Apply( time + dt,
                   &domain,
                   "ElementRegions",
                   FlowSolverBase::viewKeyStruct::pressureString,
                   [&]( FieldSpecificationBase const * const fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group * subRegion,
                        string const & ) {
                     arrayView1d< globalIndex const > const & dofNumber =
                       subRegion->getReference< array1d< globalIndex > >( presDofKey );
                     arrayView1d< integer const > const & ghostRank =
                       subRegion->group_cast< ObjectManagerBase * >()->ghostRank();

                     SortedArray< localIndex > localSet;
                     for( auto const & a : lset )
                     {
                       if( ghostRank[a] < 0 )
                       {
                         localSet.insert( a );
                       }
                     }

                     fs->ZeroSystemRowsForBoundaryCondition< LAInterface >(
                       localSet.toViewConst(),
                       dofNumber,
                       m_matrix10 );
                   } );
  m_matrix10.close();

  // debugging info.  can be trimmed once everything is working.
  if( getLogLevel() == 2 )
  {
    // Before outputting anything generate permuation matrix and permute.
    //    ElementRegionManager * const elemManager = mesh->getElemManager();

    //    LAIHelperFunctions::CreatePermutationMatrix(nodeManager,
    //                                                m_solidSolver->getSystemMatrix().numGlobalRows(),
    //                                                m_solidSolver->getSystemMatrix().numGlobalCols(),
    //                                                3,
    //                                                m_solidSolver->getDofManager().getKey( keys::TotalDisplacement ),
    //                                                m_permutationMatrix0);
    //
    //    LAIHelperFunctions::CreatePermutationMatrix(elemManager,
    //                                                m_flowSolver->getSystemMatrix().numGlobalRows(),
    //                                                m_flowSolver->getSystemMatrix().numGlobalCols(),
    //                                                1,
    //                                                m_flowSolver->getDofManager().getKey(
    // FlowSolverBase::viewKeyStruct::pressureString ),
    //                                                m_permutationMatrix1);

    m_solidSolver->getSystemMatrix().create(
      m_solidSolver->getLocalMatrix().toViewConst(),
      MPI_COMM_GEOSX );
    m_solidSolver->getSystemRhs().create(
      m_solidSolver->getLocalRhs().toViewConst(),
      MPI_COMM_GEOSX );
    m_flowSolver->getSystemMatrix().create(
      m_flowSolver->getLocalMatrix().toViewConst(),
      MPI_COMM_GEOSX );
    m_flowSolver->getSystemRhs().create( m_flowSolver->getLocalRhs().toViewConst(),
                                         MPI_COMM_GEOSX );

    //    GEOSX_LOG_RANK_0("***********************************************************");
    //    GEOSX_LOG_RANK_0("matrix00");
    //    GEOSX_LOG_RANK_0("***********************************************************");
    //    LAIHelperFunctions::PrintPermutedMatrix(m_solidSolver->getSystemMatrix(), m_permutationMatrix0, std::cout);
    //    m_solidSolver->getSystemMatrix().print(std::cout);
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix01" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    //    LAIHelperFunctions::PrintPermutedMatrix(m_matrix01, m_permutationMatrix0, m_permutationMatrix1, std::cout);
    m_matrix01.print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix10" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    //    LAIHelperFunctions::PrintPermutedMatrix(m_matrix10, m_permutationMatrix1, m_permutationMatrix0, std::cout);
    m_matrix10.print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "matrix11" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    //    LAIHelperFunctions::PrintPermutedMatrix(m_flowSolver->getSystemMatrix(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemMatrix().print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "residual0" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    //    LAIHelperFunctions::PrintPermutedVector(m_solidSolver->getSystemRhs(), m_permutationMatrix0, std::cout);
    m_solidSolver->getSystemRhs().print( std::cout );
    MpiWrapper::Barrier();

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "residual1" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    //    LAIHelperFunctions::PrintPermutedVector(m_flowSolver->getSystemRhs(), m_permutationMatrix1, std::cout);
    m_flowSolver->getSystemRhs().print( std::cout );
    MpiWrapper::Barrier();
  }

  if( getLogLevel() >= 10 )
  {
    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

    {
      string filename = "matrix00_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemMatrix().write( filename,
                                              LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix00: written to " << filename );
    }
    {
      string filename = "matrix01_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_matrix01.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix01: written to " << filename );
    }
    {
      string filename = "matrix10_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_matrix10.write( filename, LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix10: written to " << filename );
    }
    {
      string filename = "matrix11_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemMatrix().write( filename,
                                             LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "matrix11: written to " << filename );
    }
    {
      string filename = "residual0_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_solidSolver->getSystemRhs().write( filename,
                                           LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual0: written to " << filename );
    }
    {
      string filename = "residual1_" + std::to_string( time ) + "_" +
        std::to_string( newtonIter ) + ".mtx";
      m_flowSolver->getSystemRhs().write( filename,
                                          LAIOutputFormat::MATRIX_MARKET );
      GEOSX_LOG_RANK_0( "residual1: written to " << filename );
    }
  }
}

real64
HydrofractureSolver::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;

  real64 const fluidResidual =
    m_flowSolver->CalculateResidualNorm( domain,
                                         m_flowSolver->getDofManager(),
                                         m_flowSolver->getLocalRhs() );

  real64 const solidResidual =
    m_solidSolver->CalculateResidualNorm( domain,
                                          m_solidSolver->getDofManager(),
                                          m_solidSolver->getLocalRhs() );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    char output[200] = { 0 };
    sprintf( output,
             "    ( Rfluid, Rsolid ) = ( %4.2e, %4.2e )",
             fluidResidual,
             solidResidual );
    std::cout << output << std::endl;
  }

  return std::sqrt( fluidResidual * fluidResidual + solidResidual * solidResidual );
}

void
HydrofractureSolver::AssembleForceResidualDerivativeWrtPressure(
  DomainPartition & domain,
  ParallelMatrix * const matrix01,
  arrayView1d< real64 > const & rhs0 )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FaceManager const & faceManager = *mesh.getFaceManager();
  NodeManager & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap =
    faceManager.nodeList().toViewConst();

  arrayView2d< real64 > const & fext = nodeManager.getReference< array2d< real64 > >(
    SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternal );
  fext.setValues< serialPolicy >( 0 );

  string const presDofKey = m_flowSolver->getDofManager().getKey(
    FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey =
    m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  globalIndex const dispRankOffset = m_solidSolver->getDofManager().rankOffset();
  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< globalIndex_array >( dispDofKey );

  matrix01->open();

  elemManager.forElementSubRegions< FaceElementSubRegion >(
    [&]( FaceElementSubRegion const & subRegion ) {
      arrayView1d< globalIndex const > const & faceElementDofNumber =
        subRegion.getReference< array1d< globalIndex > >( presDofKey );

      if( subRegion.hasWrapper( "pressure" ) )
      {
        arrayView1d< real64 const > const & fluidPressure =
          subRegion.getReference< array1d< real64 > >( "pressure" );
        arrayView1d< real64 const > const & deltaFluidPressure =
          subRegion.getReference< array1d< real64 > >( "deltaPressure" );
        arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();
        arrayView1d< real64 const > const & area = subRegion.getElementArea();
        arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

        forAll< serialPolicy >( subRegion.size(), [=]( localIndex const kfe ) {
          R1Tensor Nbar = faceNormal[elemsToFaces[kfe][0]];
          Nbar -= faceNormal[elemsToFaces[kfe][1]];
          Nbar.Normalize();

          localIndex const kf0 = elemsToFaces[kfe][0];
          localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

          globalIndex rowDOF[24];
          real64 nodeRHS[24];
          stackArray2d< real64, 12 * 12 > dRdP( numNodesPerFace * 3, 1 );
          globalIndex colDOF = faceElementDofNumber[kfe];

          real64 const Ja = area[kfe] / numNodesPerFace;

          real64 nodalForceMag =
            ( fluidPressure[kfe] + deltaFluidPressure[kfe] ) * Ja;
          R1Tensor nodalForce( Nbar );
          nodalForce *= nodalForceMag;

          for( localIndex kf = 0; kf < 2; ++kf )
          {
            localIndex const faceIndex = elemsToFaces[kfe][kf];

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              for( int i = 0; i < 3; ++i )
              {
                rowDOF[3 * a + i] =
                  dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
                nodeRHS[3 * a + i] = -nodalForce[i] * pow( -1, kf );
                fext[faceToNodeMap( faceIndex, a )][i] +=
                  -nodalForce[i] * pow( -1, kf );

                dRdP( 3 * a + i, 0 ) = -Ja * Nbar[i] * pow( -1, kf );
              }
            }

            for( localIndex a = 0; a < numNodesPerFace; ++a )
            {
              localIndex const localRow = LvArray::integerConversion< localIndex >(
                rowDOF[3 * a] - dispRankOffset );
              if( localRow >= 0 && localRow < rhs0.size() )
              {
                for( int i = 0; i < 3; ++i )
                {
                  // TODO: use parallel atomic when loop is parallel
                  RAJA::atomicAdd( serialAtomic {},
                                   &rhs0[localRow + i],
                                   nodeRHS[3 * a + i] );
                }
              }
            }

            if( ghostRank[kfe] < 0 )
            {
              matrix01->add( rowDOF, &colDOF, dRdP.data(), numNodesPerFace * 3, 1 );
            }
          }
        } );
      }
    } );

  matrix01->close();
}

void
HydrofractureSolver::AssembleFluidMassResidualDerivativeWrtDisplacement(
  DomainPartition const & domain,
  ParallelMatrix * const matrix10 )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();
  NodeManager const & nodeManager = *mesh.getNodeManager();
  ConstitutiveManager const & constitutiveManager =
    *domain.getConstitutiveManager();

  string const presDofKey = m_flowSolver->getDofManager().getKey(
    FlowSolverBase::viewKeyStruct::pressureString );
  string const dispDofKey =
    m_solidSolver->getDofManager().getKey( keys::TotalDisplacement );

  CRSMatrixView< real64 const, localIndex const > const & dFluxResidual_dAperture =
    m_flowSolver->getDerivativeFluxResidual_dAperture().toViewConst();

  ContactRelationBase const * const contactRelation =
    constitutiveManager.GetGroup< ContactRelationBase >( m_contactRelationName );

  matrix10->open();

  forTargetSubRegionsComplete< FaceElementSubRegion >(
    mesh,
    [&]( localIndex const,
         localIndex const,
         localIndex const,
         ElementRegionBase const & region,
         FaceElementSubRegion const & subRegion ) {
      string const & fluidName =
        m_flowSolver
          ->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
      SingleFluidBase const & fluid =
        GetConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< globalIndex const > const & presDofNumber =
        subRegion.getReference< array1d< globalIndex > >( presDofKey );
      arrayView1d< globalIndex const > const & dispDofNumber =
        nodeManager.getReference< array1d< globalIndex > >( dispDofKey );

      arrayView2d< real64 const > const & dens = fluid.density();

      arrayView1d< real64 const > const & aperture = subRegion.getElementAperture();
      arrayView1d< real64 const > const & area = subRegion.getElementArea();

      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();
      ArrayOfArraysView< localIndex const > const & faceToNodeMap =
        faceManager.nodeList().toViewConst();

      arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();

      //    arrayView1d< real64 const > const & separationCoeff = subRegion.getSeparationCoefficient();
      //    arrayView1d<real64 const> const & dseparationCoeff_dAper  =
      //      subRegion.getReference<array1d<real64>>(FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString);

      forAll< serialPolicy >( subRegion.size(), [=]( localIndex ei ) {
        //if (elemGhostRank[ei] < 0)
        {
          globalIndex const elemDOF = presDofNumber[ei];
          localIndex const numNodesPerFace =
            faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );
          real64 const dAccumulationResidualdAperture = dens[ei][0] * area[ei];
          //* ( separationCoeff[ei] + aperture[ei] * dseparationCoeff_dAper[ei] );

          globalIndex nodeDOF[8 * 3];

          R1Tensor Nbar = faceNormal[elemsToFaces[ei][0]];
          Nbar -= faceNormal[elemsToFaces[ei][1]];
          Nbar.Normalize();

          stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );

          // Accumulation derivative
          if( elemGhostRank[ei] < 0 )
          {
            for( localIndex kf = 0; kf < 2; ++kf )
            {
              for( localIndex a = 0; a < numNodesPerFace; ++a )
              {
                for( int i = 0; i < 3; ++i )
                {
                  nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] =
                    dispDofNumber[faceToNodeMap( elemsToFaces[ei][kf], a )] + i;
                  real64 const dGap_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                  real64 const dAper_dU =
                    contactRelation->dEffectiveAperture_dAperture( aperture[ei] ) *
                    dGap_dU;
                  dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) =
                    dAccumulationResidualdAperture * dAper_dU;
                }
              }
            }
            matrix10->add( elemDOF, nodeDOF, dRdU.data(), 2 * numNodesPerFace * 3 );
          }

          // flux derivative
          localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
          arraySlice1d< localIndex const > const & columns =
            dFluxResidual_dAperture.getColumns( ei );
          arraySlice1d< real64 const > const & values =
            dFluxResidual_dAperture.getEntries( ei );

          for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
          {
            real64 dRdAper = values[kfe2];
            localIndex const ei2 = columns[kfe2];

            for( localIndex kf = 0; kf < 2; ++kf )
            {
              for( localIndex a = 0; a < numNodesPerFace; ++a )
              {
                for( int i = 0; i < 3; ++i )
                {
                  nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] =
                    dispDofNumber[faceToNodeMap( elemsToFaces[ei2][kf], a )] + i;
                  real64 const dGap_dU = -pow( -1, kf ) * Nbar[i] / numNodesPerFace;
                  real64 const dAper_dU =
                    contactRelation->dEffectiveAperture_dAperture( aperture[ei2] ) *
                    dGap_dU;
                  dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dAper_dU;
                }
              }
            }
            matrix10->add( elemDOF, nodeDOF, dRdU.data(), 2 * numNodesPerFace * 3 );
          }
        }
      } );
    } );

  matrix10->close();
}

void
HydrofractureSolver::ApplySystemSolution(
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localSolution ),
  real64 const scalingFactor,
  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->ApplySystemSolution( m_solidSolver->getDofManager(),
                                      m_solidSolver->getLocalSolution(),
                                      scalingFactor,
                                      domain );
  m_flowSolver->ApplySystemSolution( m_flowSolver->getDofManager(),
                                     m_flowSolver->getLocalSolution(),
                                     -scalingFactor,
                                     domain );

  UpdateDeformationForCoupling( domain );
}

}  // namespace geosx

#ifdef GEOSX_LA_INTERFACE_TRILINOS

  // For some reason this needs to be defined when using cuda
  #if defined( __CUDACC__ )
    #define KOKKOS_ENABLE_SERIAL_ATOMICS
  #endif

  #include "Epetra_FEVector.h"
  #include "Epetra_FECrsMatrix.h"
  #include "EpetraExt_MatrixMatrix.h"
  #include "Thyra_OperatorVectorClientSupport.hpp"
  #include "Thyra_AztecOOLinearOpWithSolveFactory.hpp"
  #include "Thyra_AztecOOLinearOpWithSolve.hpp"
  #include "Thyra_EpetraThyraWrappers.hpp"
  #include "Thyra_EpetraLinearOp.hpp"
  #include "Thyra_EpetraLinearOpBase.hpp"
  #include "Thyra_LinearOpBase.hpp"
  #include "Thyra_LinearOpWithSolveBase.hpp"
  #include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
  #include "Thyra_DefaultBlockedLinearOp.hpp"
  #include "Thyra_DefaultIdentityLinearOp.hpp"
  #include "Thyra_DefaultZeroLinearOp.hpp"
  #include "Thyra_DefaultLinearOpSource.hpp"
  #include "Thyra_DefaultPreconditioner.hpp"
  #include "Thyra_EpetraThyraWrappers.hpp"
  #include "Thyra_PreconditionerFactoryHelpers.hpp"
  #include "Thyra_VectorStdOps.hpp"
  #include "Thyra_PreconditionerFactoryHelpers.hpp"
  #include "Thyra_DefaultInverseLinearOp.hpp"
  #include "Thyra_PreconditionerFactoryBase.hpp"
  #include "Thyra_get_Epetra_Operator.hpp"
  #include "Thyra_MLPreconditionerFactory.hpp"
  #include "Teuchos_ParameterList.hpp"
  #include "Teuchos_RCP.hpp"
  #include "Teuchos_Time.hpp"
  #include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#endif

namespace geosx
{
void
HydrofractureSolver::SolveSystem(
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  ParallelMatrix &,
  ParallelVector &,
  ParallelVector & )
{
  GEOSX_MARK_FUNCTION;

#if defined( GEOSX_LA_INTERFACE_TRILINOS )
  /*
     globalIndex numU = m_solidSolver->getSystemRhs().globalSize();
     globalIndex numP = m_flowSolver->getSystemRhs().globalSize();
     GEOSX_LOG_RANK_0("size = " << numU << " + " << numP);
   */

  integer const newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;

  using namespace Teuchos;
  using namespace Thyra;

  Teuchos::Time clock( "solveClock" );

  GEOSX_MARK_BEGIN( Setup );
  Epetra_FECrsMatrix * p_matrix[2][2];
  Epetra_FEVector * p_rhs[2];
  Epetra_FEVector * p_solution[2];

  m_solidSolver->getSystemRhs().create( m_solidSolver->getLocalRhs(),
                                        MPI_COMM_GEOSX );
  m_flowSolver->getSystemRhs().create( m_flowSolver->getLocalRhs(),
                                       MPI_COMM_GEOSX );

  p_rhs[0] = &m_solidSolver->getSystemRhs().unwrapped();
  p_rhs[1] = &m_flowSolver->getSystemRhs().unwrapped();

  m_solidSolver->getSystemSolution().create( m_solidSolver->getLocalSolution(),
                                             MPI_COMM_GEOSX );
  m_flowSolver->getSystemSolution().create( m_flowSolver->getLocalSolution(),
                                            MPI_COMM_GEOSX );

  p_solution[0] = &m_solidSolver->getSystemSolution().unwrapped();
  p_solution[1] = &m_flowSolver->getSystemSolution().unwrapped();

  m_solidSolver->getSystemMatrix().create(
    m_solidSolver->getLocalMatrix().toViewConst(),
    MPI_COMM_GEOSX );
  m_flowSolver->getSystemMatrix().create(
    m_flowSolver->getLocalMatrix().toViewConst(),
    MPI_COMM_GEOSX );

  p_matrix[0][0] = &m_solidSolver->getSystemMatrix().unwrapped();
  p_matrix[0][1] = &m_matrix01.unwrapped();
  p_matrix[1][0] = &m_matrix10.unwrapped();
  p_matrix[1][1] = &m_flowSolver->getSystemMatrix().unwrapped();

  // scale and symmetrize

  m_densityScaling = 1e-3;
  m_pressureScaling = 1e9;

  p_matrix[0][1]->Scale( m_pressureScaling );
  p_matrix[1][0]->Scale( m_pressureScaling * m_densityScaling );
  p_matrix[1][1]->Scale( m_pressureScaling * m_pressureScaling * m_densityScaling );
  p_rhs[1]->Scale( m_pressureScaling * m_densityScaling );

  // SCHEME CHOICES
  //
  // there are several flags to control solver behavior.
  // these should be compared in a scaling study.
  //
  // -- whether to use a block diagonal or a
  //    block triangular preconditioner.
  // -- whether to use BiCGstab or GMRES for the
  //    krylov solver.  GMRES is generally more robust,
  //    BiCGstab sometimes shows better parallel performance.
  //    false is probably better.

  LinearSolverParameters const & linParams = m_linearSolverParameters.get();

  const bool use_diagonal_prec = true;
  const bool use_bicgstab = ( linParams.solverType == "bicgstab" );

  // set initial guess to zero

  p_solution[0]->PutScalar( 0.0 );
  p_solution[1]->PutScalar( 0.0 );

  // create separate displacement component matrix

  clock.start( true );
  if( newtonIter == 0 )
  {
    m_blockDiagUU.reset( new ParallelMatrix() );
    LAIHelperFunctions::SeparateComponentFilter( m_solidSolver->getSystemMatrix(),
                                                 *m_blockDiagUU,
                                                 3 );
  }

  // create schur complement approximation matrix

  Epetra_CrsMatrix * schurApproxPP =
    NULL;  // confirm we delete this at end of function!
  {
    Epetra_Vector diag( p_matrix[0][0]->RowMap() );
    Epetra_Vector diagInv( p_matrix[0][0]->RowMap() );

    p_matrix[0][0]->ExtractDiagonalCopy( diag );
    diagInv.Reciprocal( diag );

    Epetra_FECrsMatrix DB( *p_matrix[0][1] );
    DB.LeftScale( diagInv );
    DB.FillComplete();

    Epetra_FECrsMatrix BtDB( Epetra_DataAccess::Copy, p_matrix[1][1]->RowMap(), 1 );
    EpetraExt::MatrixMatrix::Multiply( *p_matrix[1][0], false, DB, false, BtDB );
    EpetraExt::MatrixMatrix::Add( BtDB,
                                  false,
                                  -1.0,
                                  *p_matrix[1][1],
                                  false,
                                  1.0,
                                  schurApproxPP );

    schurApproxPP->FillComplete();
  }
  double auxTime = clock.stop();
  GEOSX_MARK_END( Setup );

  // we want to use thyra to wrap epetra operators and vectors
  // for individual blocks.  this is an ugly conversion, but
  // it is basically just window dressing.
  //
  // note the use of Teuchos::RCP reference counted pointers.
  // The general syntax is usually one of:
  //
  //   RCP<T> Tptr = rcp(new T)
  //   RCP<T> Tptr = nonMemberConstructor();
  //   RCP<T> Tptr (t_ptr,false)
  //
  // where "false" implies the RCP does not own the object and
  // should not attempt to delete it when finished.

  GEOSX_MARK_BEGIN( THYRA_SETUP );

  RCP< const Thyra::LinearOpBase< double > > matrix_block[2][2];
  RCP< Thyra::MultiVectorBase< double > > lhs_block[2];
  RCP< Thyra::MultiVectorBase< double > > rhs_block[2];

  for( unsigned i = 0; i < 2; ++i )
    for( unsigned j = 0; j < 2; ++j )
    {
      RCP< Epetra_Operator > mmm( &*p_matrix[i][j], false );
      matrix_block[i][j] = Thyra::epetraLinearOp( mmm );
    }

  RCP< Epetra_Operator > bbb( &m_blockDiagUU->unwrapped(), false );
  RCP< Epetra_Operator > ppp( schurApproxPP, false );

  RCP< const Thyra::LinearOpBase< double > > blockDiagOp = Thyra::epetraLinearOp( bbb );
  RCP< const Thyra::LinearOpBase< double > > schurOp = Thyra::epetraLinearOp( ppp );

  for( unsigned i = 0; i < 2; ++i )
  {
    RCP< Epetra_MultiVector > lll( &*p_solution[i], false );
    RCP< Epetra_MultiVector > rrr( &*p_rhs[i], false );

    lhs_block[i] = Thyra::create_MultiVector( lll, matrix_block[i][i]->domain() );
    rhs_block[i] = Thyra::create_MultiVector( rrr, matrix_block[i][i]->range() );
  }

  // now use thyra to create an operator representing
  // the full block 2x2 system

  RCP< const Thyra::LinearOpBase< double > > matrix =
    Thyra::block2x2( matrix_block[0][0],
                     matrix_block[0][1],
                     matrix_block[1][0],
                     matrix_block[1][1] );

  // creating a representation of the blocked
  // rhs and lhs is a little uglier.

  RCP< Thyra::ProductMultiVectorBase< double > > rhs;
  {
    Teuchos::Array< RCP< Thyra::MultiVectorBase< double > > > mva;
    Teuchos::Array< RCP< const Thyra::VectorSpaceBase< double > > > mvs;

    for( unsigned i = 0; i < 2; ++i )
    {
      mva.push_back( rhs_block[i] );
      mvs.push_back( rhs_block[i]->range() );
    }

    RCP< const Thyra::DefaultProductVectorSpace< double > > vs =
      Thyra::productVectorSpace< double >( mvs );

    rhs = Thyra::defaultProductMultiVector< double >( vs, mva );
  }

  RCP< Thyra::ProductMultiVectorBase< double > > lhs;

  {
    Teuchos::Array< RCP< Thyra::MultiVectorBase< double > > > mva;
    Teuchos::Array< RCP< const Thyra::VectorSpaceBase< double > > > mvs;

    for( unsigned i = 0; i < 2; ++i )
    {
      mva.push_back( lhs_block[i] );
      mvs.push_back( lhs_block[i]->range() );
    }

    RCP< const Thyra::DefaultProductVectorSpace< double > > vs =
      Thyra::productVectorSpace< double >( mvs );

    lhs = Thyra::defaultProductMultiVector< double >( vs, mva );
  }

  GEOSX_MARK_END( THYRA_SETUP );

  // for the preconditioner, we need two approximate inverses,
  // we store both "sub operators" in a 1x2 array:

  RCP< const Thyra::LinearOpBase< double > > sub_op[2];

  clock.start( true );
  GEOSX_MARK_BEGIN( PRECONDITIONER );

  for( unsigned i = 0; i < 2; ++i )  // loop over diagonal blocks
  {
    RCP< Teuchos::ParameterList > list =
      rcp( new Teuchos::ParameterList( "precond_list" ), true );

    if( linParams.preconditionerType == "amg" )
    {
      list->set( "Preconditioner Type", "ML" );
      list->sublist( "Preconditioner Types" )
        .sublist( "ML" )
        .set( "Base Method Defaults", "SA" );
      list->sublist( "Preconditioner Types" )
        .sublist( "ML" )
        .sublist( "ML Settings" )
        .set( "PDE equations", ( i == 0 ? 3 : 1 ) );
      list->sublist( "Preconditioner Types" )
        .sublist( "ML" )
        .sublist( "ML Settings" )
        .set( "ML output", 0 );
      list->sublist( "Preconditioner Types" )
        .sublist( "ML" )
        .sublist( "ML Settings" )
        .set( "aggregation: type", "Uncoupled" );
      list->sublist( "Preconditioner Types" )
        .sublist( "ML" )
        .sublist( "ML Settings" )
        .set( "aggregation: threshold", 1e-3 );

      if( i == 0 )  // smoother for mechanics block
      {
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "smoother: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "smoother: sweeps", 3 );
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "coarse: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "coarse: sweeps", 3 );
      }
      else  // smoother for flow block
      {
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "smoother: type", "Chebyshev" );
        list->sublist( "Preconditioner Types" )
          .sublist( "ML" )
          .sublist( "ML Settings" )
          .set( "smoother: sweeps", 3 );
      }
    }
    else  // use ILU for both blocks
    {
      list->set( "Preconditioner Type", "Ifpack" );
      list->sublist( "Preconditioner Types" )
        .sublist( "Ifpack" )
        .set( "Prec Type", "ILU" );
    }

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( list );

    RCP< const Thyra::PreconditionerFactoryBase< double > > strategy =
      createPreconditioningStrategy( builder );
    RCP< Thyra::PreconditionerBase< double > > tmp;

    if( i == 0 )
      tmp = prec( *strategy, blockDiagOp );
    else
      tmp = prec( *strategy, schurOp );
    //tmp = prec(*strategy,matrix_block[i][i]);

    sub_op[i] = tmp->getUnspecifiedPrecOp();
  }

  // create zero operators for off diagonal blocks

  RCP< const Thyra::LinearOpBase< double > > zero_01 =
    rcp( new Thyra::DefaultZeroLinearOp< double >( matrix_block[0][0]->range(),
                                                   matrix_block[1][1]->domain() ) );

  RCP< const Thyra::LinearOpBase< double > > zero_10 =
    rcp( new Thyra::DefaultZeroLinearOp< double >( matrix_block[1][1]->range(),
                                                   matrix_block[0][0]->domain() ) );

  // now build the block preconditioner

  RCP< const Thyra::LinearOpBase< double > > preconditioner;

  if( use_diagonal_prec )
  {
    preconditioner = Thyra::block2x2( sub_op[0], zero_01, zero_10, sub_op[1] );
  }
  else
  {
    RCP< const Thyra::LinearOpBase< double > > eye_00 = Teuchos::rcp(
      new Thyra::DefaultIdentityLinearOp< double >( matrix_block[0][0]->range() ) );

    RCP< const Thyra::LinearOpBase< double > > eye_11 = Teuchos::rcp(
      new Thyra::DefaultIdentityLinearOp< double >( matrix_block[1][1]->range() ) );

    RCP< const Thyra::LinearOpBase< double > > mAinvB1, mB2Ainv;

    mAinvB1 = Thyra::scale( -1.0, Thyra::multiply( sub_op[0], matrix_block[0][1] ) );
    mB2Ainv = Thyra::scale( -1.0, Thyra::multiply( matrix_block[1][0], sub_op[0] ) );

    RCP< const Thyra::LinearOpBase< double > > Linv, Dinv, Uinv, Eye;

    Linv = Thyra::block2x2( eye_00, zero_01, mB2Ainv, eye_11 );
    Dinv = Thyra::block2x2( sub_op[0], zero_01, zero_10, sub_op[1] );
    Uinv = Thyra::block2x2( eye_00, mAinvB1, zero_10, eye_11 );

    //preconditioner = Thyra::multiply(Uinv,Dinv);
    //preconditioner = Thyra::multiply(Dinv,Linv);
    preconditioner = Thyra::multiply( Uinv, Dinv, Linv );
  }

  GEOSX_MARK_END( PRECONDITIONER );
  double setupTime = clock.stop();

  // define solver strategy for blocked system. this is
  // similar but slightly different from the sub operator
  // construction, since now we have a user defined preconditioner

  {
    RCP< Teuchos::ParameterList > list = rcp( new Teuchos::ParameterList( "list" ) );

    list->set( "Linear Solver Type", "AztecOO" );
    list->set( "Preconditioner Type", "None" );  // will use user-defined P
    list->sublist( "Linear Solver Types" )
      .sublist( "AztecOO" )
      .sublist( "Forward Solve" )
      .set( "Max Iterations", linParams.krylov.maxIterations );
    list->sublist( "Linear Solver Types" )
      .sublist( "AztecOO" )
      .sublist( "Forward Solve" )
      .set( "Tolerance", linParams.krylov.relTolerance );

    if( use_bicgstab )
      list->sublist( "Linear Solver Types" )
        .sublist( "AztecOO" )
        .sublist( "Forward Solve" )
        .sublist( "AztecOO Settings" )
        .set( "Aztec Solver", "BiCGStab" );
    else
      list->sublist( "Linear Solver Types" )
        .sublist( "AztecOO" )
        .sublist( "Forward Solve" )
        .sublist( "AztecOO Settings" )
        .set( "Aztec Solver", "GMRES" );

    if( linParams.logLevel > 1 )
      list->sublist( "Linear Solver Types" )
        .sublist( "AztecOO" )
        .sublist( "Forward Solve" )
        .sublist( "AztecOO Settings" )
        .set( "Output Frequency", 1 );

    Stratimikos::DefaultLinearSolverBuilder builder;
    builder.setParameterList( list );

    RCP< const Thyra::LinearOpWithSolveFactoryBase< double > > strategy =
      createLinearSolveStrategy( builder );
    RCP< Thyra::LinearOpWithSolveBase< double > > solver = strategy->createOp();

    Thyra::initializePreconditionedOp< double >(
      *strategy,
      matrix,
      Thyra::rightPrec< double >( preconditioner ),
      solver.ptr() );

    clock.start( true );
    GEOSX_MARK_BEGIN( SOLVER );

    // !!!! Actual Solve !!!!
    Thyra::SolveStatus< double > status =
      solver->solve( Thyra::NOTRANS, *rhs, lhs.ptr() );

    GEOSX_MARK_END( SOLVER );
    double solveTime = clock.stop();

    /* TODO: replace with SolverBase status output */

    integer numKrylovIter = status.extraParameters->get< int >( "Iteration Count" );
    if( getLogLevel() >= 2 )
    {
      GEOSX_LOG_RANK_0( "\t\tLinear Solver | Iter = "
                        << numKrylovIter << " | TargetReduction "
                        << linParams.krylov.relTolerance << " | AuxTime "
                        << auxTime << " | SetupTime " << setupTime
                        << " | SolveTime " << solveTime );
    }

    p_solution[1]->Scale( m_pressureScaling );
    p_rhs[1]->Scale( 1 / ( m_pressureScaling * m_densityScaling ) );
  }

  delete schurApproxPP;

  //TODO: remove all this once everything is working
  if( getLogLevel() == 2 )
  {
    /*
       ParallelVector permutedSol;
       ParallelVector const & solution = m_solidSolver->getSystemSolution();
       permutedSol.createWithLocalSize(m_solidSolver->getSystemMatrix().numLocalRows(), MPI_COMM_GEOSX);
       m_permutationMatrix0.multiply(solution, permutedSol);
       permutedSol.close();
     */

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "solution0" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    p_solution[0]->Print( std::cout );
    std::cout << std::endl;
    MPI_Barrier( MPI_COMM_GEOSX );

    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    GEOSX_LOG_RANK_0( "solution1" );
    GEOSX_LOG_RANK_0(
      "***********************************************************" );
    p_solution[1]->Print( std::cout );
  }

  m_solidSolver->getSystemSolution().extract( m_solidSolver->getLocalSolution() );
  m_flowSolver->getSystemSolution().extract( m_flowSolver->getLocalSolution() );
#else
  GEOSX_ERROR( "Only implemented for trilinos." );
#endif
}

real64
HydrofractureSolver::ScalingForSystemSolution(
  DomainPartition const & domain,
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localSolution ) )
{
  return m_solidSolver->ScalingForSystemSolution(
    domain,
    m_solidSolver->getDofManager(),
    m_solidSolver->getLocalSolution() );
}

void
HydrofractureSolver::SetNextDt( real64 const & currentDt, real64 & nextDt )
{
  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->SetNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    SolverBase * const surfaceGenerator =
      this->getParent()->GetGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator->GetTimestepRequest() < 1e99
      ? surfaceGenerator->GetTimestepRequest()
      : currentDt;
  }
  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is " << nextDt );
}

void
HydrofractureSolver::initializeNewFaceElements( DomainPartition const & )
{
  //  m_flowSolver->
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        HydrofractureSolver,
                        std::string const &,
                        Group * const )
} /* namespace geosx */
