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
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/SurfaceElementRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "physicsSolvers/surfaceGeneration/SurfaceGenerator.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

HydrofractureSolver::HydrofractureSolver( const string & name,
                                          Group * const parent ):
  PoroelasticSolver( name, parent ),
  m_contactRelationName(),
  m_surfaceGeneratorName(),
  m_surfaceGenerator( nullptr ),
  m_maxNumResolves( 10 )
{
  registerWrapper( viewKeyStruct::surfaceGeneratorNameString(), &m_surfaceGeneratorName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the surface generator to use in the hydrofracture solver" );

  registerWrapper( viewKeyStruct::contactRelationNameString(), &m_contactRelationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of contact relation to enforce constraints on fracture boundary." );

  registerWrapper( viewKeyStruct::maxNumResolvesString(), &m_maxNumResolves ).
    setApplyDefaultValue( 10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Value to indicate how many resolves may be executed to perform surface generation after the execution of flow and mechanics solver. " );

  m_numResolves[0] = 0;

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hydrofracture;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
void HydrofractureSolver::RegisterDataOnMesh( dataRepository::Group & MeshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = *meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion * const region )
    {
      region->forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion * const subRegion )
      {
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::separationCoeff0String() ).
          setRestartFlags( RestartFlags::NO_WRITE );
        subRegion->registerWrapper< array1d< real64 > >( viewKeyStruct::apertureAtFailureString() ).
          setApplyDefaultValue( -1.0 ).
          setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion->registerWrapper< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
      } );
    } );
  } );
}
#endif

HydrofractureSolver::~HydrofractureSolver()
{
  // TODO Auto-generated destructor stub
}

void HydrofractureSolver::implicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  updateDeformationForCoupling( domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
  m_flowSolver->implicitStepSetup( time_n, dt, domain );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  mesh.getElemManager().forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion & faceElemRegion )
  {
    faceElemRegion.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
    {
      arrayView1d< real64 > const &
      separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String() );
      arrayView1d< real64 const > const &
      separationCoeff = subRegion.getSeparationCoefficient();
      for( localIndex k=0; k<separationCoeff0.size(); ++k )
      {
        separationCoeff0[k] = separationCoeff[k];
      }
    } );
  } );
#endif

}

void HydrofractureSolver::postProcessInput()
{
  PoroelasticSolver::postProcessInput();
  m_surfaceGenerator = &this->getParent().getGroup< SurfaceGenerator >( m_surfaceGeneratorName );
}

void HydrofractureSolver::initializePostInitialConditionsPreSubGroups()
{}



void HydrofractureSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );
}

real64 HydrofractureSolver::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
{
  real64 dtReturn = dt;

  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    dtReturn = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {

    implicitStepSetup( time_n, dt, domain );

    int const maxIter = m_maxNumResolves + 1;
    m_numResolves[1] = m_numResolves[0];
    int solveIter;
    for( solveIter=0; solveIter<maxIter; ++solveIter )
    {
      int locallyFractured = 0;
      int globallyFractured = 0;

      setupSystem( domain,
                   m_dofManager,
                   m_localMatrix,
                   m_localRhs,
                   m_localSolution );

      // currently the only method is implicit time integration
      dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

      if( m_surfaceGenerator!=nullptr )  // is this check really necessary ?
      {
        if( m_surfaceGenerator->solverStep( time_n, dt, cycleNumber, domain ) > 0 )
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
        fieldNames["elems"].emplace_back( string( FlowSolverBase::viewKeyStruct::pressureString() ) );
        fieldNames["elems"].emplace_back( "elementAperture" );

        CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                             domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                             domain.getNeighbors(),
                                                             false );

        this->updateDeformationForCoupling( domain );

        if( getLogLevel() >= 1 )
        {
          GEOSX_LOG_RANK_0( "++ Fracture propagation. Re-entering Newton Solve." );
        }
        m_flowSolver->resetViews( domain.getMeshBody( 0 ).getMeshLevel( 0 ) );
      }
    }

    // final step for completion of timestep. typically secondary variable updates and cleanup.
    implicitStepComplete( time_n, dtReturn, domain );
    m_numResolves[1] = solveIter;
  }

  return dtReturn;
}

void HydrofractureSolver::updateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();
  NodeManager const & nodeManager = meshLevel.getNodeManager();
  FaceManager & faceManager = meshLevel.getFaceManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const u = nodeManager.totalDisplacement();
  arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();
  // arrayView1d<real64 const> const faceArea = faceManager.faceArea();
  ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();

  ContactRelationBase const & contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_contactRelationName );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion & subRegion )
  {
    arrayView1d< real64 > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 > const effectiveAperture = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString() );
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();
    arrayView1d< real64 > const deltaVolume = subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString() );
    arrayView1d< real64 const > const area = subRegion.getElementArea();
    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    arrayView1d< real64 const > const &
    apertureF = subRegion.getReference< array1d< real64 > >( viewKeyStruct::apertureAtFailureString() );

    arrayView1d< real64 > const &
    separationCoeff = subRegion.getSeparationCoefficient();

    arrayView1d< real64 > const &
    dSeparationCoeff_dAper = subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
    arrayView1d< real64 const > const &
    separationCoeff0 = subRegion.getReference< array1d< real64 > >( viewKeyStruct::separationCoeff0String() );
#endif

    forAll< serialPolicy >( subRegion.size(), [=, &contactRelation] ( localIndex const kfe )
    {
      localIndex const kf0 = elemsToFaces[kfe][0];
      localIndex const kf1 = elemsToFaces[kfe][1];
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );
      real64 temp[ 3 ] = { 0 };
      for( localIndex a=0; a<numNodesPerFace; ++a )
      {
        LvArray::tensorOps::add< 3 >( temp, u[ faceToNodeMap( kf0, a ) ] );
        LvArray::tensorOps::subtract< 3 >( temp, u[ faceToNodeMap( kf1, a ) ] );
      }

      // TODO this needs a proper contact based strategy for aperture
      aperture[kfe] = -LvArray::tensorOps::AiBi< 3 >( temp, faceNormal[ kf0 ] ) / numNodesPerFace;

      effectiveAperture[kfe] = contactRelation.effectiveAperture( aperture[kfe] );


#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
      real64 const s = aperture[kfe] / apertureF[kfe];
      if( separationCoeff0[kfe]<1.0 && s>separationCoeff0[kfe] )
      {
        if( s >= 1.0 )
        {
          separationCoeff[kfe] = 1.0;
          dSeparationCoeff_dAper[kfe] = 0.0;
        }
        else
        {
          separationCoeff[kfe] = s;
          dSeparationCoeff_dAper[kfe] = 1.0/apertureF[kfe];
        }
      }
#endif
      deltaVolume[kfe] = effectiveAperture[kfe] * area[kfe] - volume[kfe];
    } );

//#if defined(USE_CUDA)
//    deltaVolume.move( LvArray::MemorySpace::cuda );
//    aperture.move( LvArray::MemorySpace::cuda );
//    effectiveAperture.move( LvArray::MemorySpace::cuda );
//#endif
  } );
}

real64 HydrofractureSolver::splitOperatorStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
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
//    m_flowSolver->applyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_flowSolver->SolveSystem( getLinearSystemRepository(),
//                 getLinearSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_flowSolver->applySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if (dtReturnTemporary < dtReturn)
//    {
//      iter = 0;
//      dtReturn = dtReturnTemporary;
//      continue;
//    }
//
////    if (m_fluidSolver->getLinearSolverParameters()->numNewtonIterations() == 0 && iter > 0 && getLogLevel() >= 1)
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
//    m_solidSolver->applyBoundaryConditions( domain, getLinearSystemRepository(), time_n, dt );
//
//    // call the default linear solver on the system
//    m_solidSolver->SolveSystem( getLinearSystemRepository(),
//                 getLinearSolverParameters() );
//
//    // apply the system solution to the fields/variables
//    m_solidSolver->applySystemSolution( getLinearSystemRepository(), 1.0, domain );
//
//    if( m_flowSolver->calculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() &&
//        m_solidSolver->calculateResidualNorm( getLinearSystemRepository(), domain ) < solverParams->newtonTol() )
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
////    if (m_solidSolver->getLinearSolverParameters()->numNewtonIterations() > 0)
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

real64 HydrofractureSolver::explicitStep( real64 const & time_n,
                                          real64 const & dt,
                                          const int cycleNumber,
                                          DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->explicitStep( time_n, dt, cycleNumber, domain );
  m_flowSolver->solverStep( time_n, dt, cycleNumber, domain );

  return dt;
}


void HydrofractureSolver::setupDofs( DomainPartition const & domain,
                                     DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;

  // TODO: do we still need to avoid poroelasticity in the rock matrix? If not this can be done by
  // calling PoroelasticSolver::setupDofs( domain, dofManager );

  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  // restrict coupling to fracture regions only (as done originally in setupSystem)
  ElementRegionManager const & elemManager = domain.getMeshBody( 0 ).getMeshLevel( 0 ).getElemManager();
  string_array fractureRegions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & elementRegion )
  {
    fractureRegions.emplace_back( elementRegion.getName() );
  } );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem,
                          fractureRegions );

}

void HydrofractureSolver::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       array1d< real64 > & localRhs,
                                       array1d< real64 > & localSolution,
                                       bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  m_flowSolver->resetViews( mesh );

  dofManager.setMesh( mesh );

  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > patternOriginal;
  dofManager.setSparsityPattern( patternOriginal );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternOriginal.numRows() );
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternOriginal.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling
  addFluxApertureCouplingNNZ( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternOriginal.numRows(),
                                                         patternOriginal.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternOriginal.numRows(); ++localRow )
  {
    globalIndex const * cols = patternOriginal.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternOriginal.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addFluxApertureCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );

  localRhs.resize( numLocalRows );
  localSolution.resize( numLocalRows );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );

  setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
}

void HydrofractureSolver::addFluxApertureCouplingNNZ( DomainPartition & domain,
                                                      DofManager & dofManager,
                                                      arrayView1d< localIndex > const & rowLengths ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

      arrayView1d< globalIndex const > const faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];
        globalIndex const rowNumber = activeFlowDOF - rankOffset;

        if( rowNumber >= 0 && rowNumber < rowLengths.size() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
            // so we only add the coupling with the nodal displacements of the neighbours.
            if( k1 != k0 )
            {
              localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();
              rowLengths[rowNumber] += 3*numNodesPerElement;
            }
          }
        }
      }
    }  //);
  } );

}

void HydrofractureSolver::addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                                                  DofManager & dofManager,
                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = dofManager.rankOffset();

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      FaceElementSubRegion const & elementSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< FaceElementSubRegion >( sesri[iconn][0] );

      ArrayOfArraysView< localIndex const > const elemsToNodes = elementSubRegion.nodeList().toViewConst();

      arrayView1d< globalIndex const > const faceElementDofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( presDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = faceElementDofNumber[sei[iconn][k0]];

        globalIndex const rowIndex = activeFlowDOF - rankOffset;

        if( rowIndex >= 0 && rowIndex < pattern.numRows() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the nodal displacements of the cell itself has already been added by the dofManager
            // so we only add the coupling with the nodal displacements of the neighbours.
            if( k1 != k0 )
            {
              localIndex const numNodesPerElement = elemsToNodes[sei[iconn][k1]].size();

              for( localIndex a=0; a<numNodesPerElement; ++a )
              {
                for( int d=0; d<3; ++d )
                {
                  globalIndex const colIndex = dispDofNumber[elemsToNodes[sei[iconn][k1]][a]] + d;
                  pattern.insertNonZero( rowIndex, colIndex );
                }
              }
            }
          }
        }
      }
    }   //);
  } );
}

void HydrofractureSolver::assembleSystem( real64 const time,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->assembleSystem( time,
                                 dt,
                                 domain,
                                 dofManager,
                                 localMatrix,
                                 localRhs );

  m_flowSolver->resetViews( domain.getMeshBody( 0 ).getMeshLevel( 0 ) );

  m_flowSolver->assembleAccumulationTerms< parallelDevicePolicy<> >( domain,
                                                                     dofManager,
                                                                     localMatrix,
                                                                     localRhs );
  m_flowSolver->assembleFluxTerms( time,
                                   dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );


  assembleForceResidualDerivativeWrtPressure( domain, localMatrix, localRhs );

  assembleFluidMassResidualDerivativeWrtDisplacement( domain, localMatrix );
}

void HydrofractureSolver::applyBoundaryConditions( real64 const time,
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

  m_flowSolver->applyBoundaryConditions( time,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}


void
HydrofractureSolver::
  assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();

  arrayView2d< real64 const > const & faceNormal = faceManager.faceNormal();
  ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

  arrayView2d< real64 > const &
  fext = nodeManager.getReference< array2d< real64 > >( SolidMechanicsLagrangianFEM::viewKeyStruct::forceExternalString() );
  fext.zero();

  string const presDofKey = m_dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dispDofKey );

  elemManager.forElementSubRegions< FaceElementSubRegion >( [&]( FaceElementSubRegion const & subRegion )
  {

    arrayView1d< globalIndex const > const &
    pressureDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );

    if( subRegion.hasWrapper( "pressure" ) )
    {
      arrayView1d< real64 const > const & fluidPressure = subRegion.getReference< array1d< real64 > >( "pressure" );
      arrayView1d< real64 const > const & deltaFluidPressure = subRegion.getReference< array1d< real64 > >( "deltaPressure" );
      arrayView1d< real64 const > const & area = subRegion.getElementArea();
      arrayView2d< localIndex const > const & elemsToFaces = subRegion.faceList();

      forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const kfe )
      {
        constexpr int kfSign[2] = { -1, 1 };

        real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[kfe][0]] );
        LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[kfe][1]] );
        LvArray::tensorOps::normalize< 3 >( Nbar );

        localIndex const kf0 = elemsToFaces[kfe][0];
        localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( kf0 );

        // TODO make if work for any element type.
        globalIndex rowDOF[24];   // Here it assumes 8 nodes?
        real64 nodeRHS[24];   // Here it assumes 8 nodes?
        stackArray2d< real64, 12*12 > dRdP( numNodesPerFace*3, 1 );
        globalIndex colDOF = pressureDofNumber[kfe];

        real64 const Ja = area[kfe] / numNodesPerFace;

        real64 nodalForceMag = ( fluidPressure[kfe]+deltaFluidPressure[kfe] ) * Ja;
        real64 nodalForce[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( Nbar );
        LvArray::tensorOps::scale< 3 >( nodalForce, nodalForceMag );

        for( localIndex kf=0; kf<2; ++kf )
        {
          localIndex const faceIndex = elemsToFaces[kfe][kf];


          for( localIndex a=0; a<numNodesPerFace; ++a )
          {

            for( int i=0; i<3; ++i )
            {
              rowDOF[3*a+i] = dispDofNumber[faceToNodeMap( faceIndex, a )] + i;
              nodeRHS[3*a+i] = nodalForce[i] * kfSign[kf];
              fext[faceToNodeMap( faceIndex, a )][i] += nodalForce[i] * kfSign[kf];

              dRdP( 3*a+i, 0 ) = Ja * Nbar[i] * kfSign[kf];
            }
          }

          for( localIndex a=0; a<numNodesPerFace; ++a )
          {
            localIndex const localRow = LvArray::integerConversion< localIndex >( rowDOF[3*a] - rankOffset );
            if( localRow >= 0 && localRow < localMatrix.numRows() )
            {
              for( int i=0; i<3; ++i )
              {
                // TODO: use parallel atomic when loop is parallel
                RAJA::atomicAdd( serialAtomic{}, &localRhs[localRow + i], nodeRHS[3*a+i] );
                localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + i,
                                                                                  &colDOF,
                                                                                  &dRdP[3*a+i][0],
                                                                                  1 );
              }
            }
          }

        }
      } );
    }
  } );
}

void
HydrofractureSolver::
  assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager const & faceManager = mesh.getFaceManager();
  NodeManager const & nodeManager = mesh.getNodeManager();
  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();

  string const presDofKey = m_dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );
  string const dispDofKey = m_dofManager.getKey( keys::TotalDisplacement );

  globalIndex const rankOffset = m_dofManager.rankOffset();

  CRSMatrixView< real64 const, localIndex const > const
  dFluxResidual_dAperture = getDerivativeFluxResidual_dAperture().toViewConst();

  ContactRelationBase const &
  contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_contactRelationName );

  forTargetSubRegionsComplete< FaceElementSubRegion >( mesh,
                                                       [&]( localIndex const,
                                                            localIndex const,
                                                            localIndex const,
                                                            ElementRegionBase const & region,
                                                            FaceElementSubRegion const & subRegion )
  {
    string const & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

    arrayView1d< globalIndex const > const presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< globalIndex const > const dispDofNumber = nodeManager.getReference< array1d< globalIndex > >( dispDofKey );

    arrayView2d< real64 const > const dens = fluid.density();

    arrayView1d< real64 const > const aperture = subRegion.getElementAperture();
    arrayView1d< real64 const > const area = subRegion.getElementArea();

    arrayView2d< localIndex const > const elemsToFaces = subRegion.faceList();
    ArrayOfArraysView< localIndex const > const faceToNodeMap = faceManager.nodeList().toViewConst();

    arrayView2d< real64 const > const faceNormal = faceManager.faceNormal();

    forAll< serialPolicy >( subRegion.size(), [=, &contactRelation]( localIndex ei )
    {
      constexpr int kfSign[2] = { -1, 1 };

      globalIndex const elemDOF = presDofNumber[ei];
      // row number associated to the pressure dof
      globalIndex rowNumber = elemDOF - rankOffset;
      localIndex const numNodesPerFace = faceToNodeMap.sizeOfArray( elemsToFaces[ei][0] );
      real64 const dAccumulationResidualdAperture = dens[ei][0] * area[ei];

      globalIndex nodeDOF[8 * 3];

      real64 Nbar[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( faceNormal[elemsToFaces[ei][0]] );
      LvArray::tensorOps::subtract< 3 >( Nbar, faceNormal[elemsToFaces[ei][1]] );
      LvArray::tensorOps::normalize< 3 >( Nbar );

      stackArray1d< real64, 24 > dRdU( 2 * numNodesPerFace * 3 );

      // Accumulation derivative
      for( localIndex kf = 0; kf < 2; ++kf )
      {
        for( localIndex a = 0; a < numNodesPerFace; ++a )
        {
          for( int i = 0; i < 3; ++i )
          {
            nodeDOF[kf * 3 * numNodesPerFace + 3 * a + i] = dispDofNumber[faceToNodeMap( elemsToFaces[ei][kf], a )] + i;
            real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;
            real64 const dAper_dU = contactRelation.dEffectiveAperture_dAperture( aperture[ei] ) * dGap_dU;
            dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dAccumulationResidualdAperture * dAper_dU;
          }
        }
      }

      if( rowNumber >= 0  && rowNumber < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                          nodeDOF,
                                                                          dRdU.data(),
                                                                          2 * numNodesPerFace * 3 );
      }

      // flux derivative
      localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( ei );
      arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( ei );

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
              real64 const dGap_dU = kfSign[kf] * Nbar[i] / numNodesPerFace;
              real64 const
              dAper_dU = contactRelation.dEffectiveAperture_dAperture( aperture[ei2] ) * dGap_dU;
              dRdU( kf * 3 * numNodesPerFace + 3 * a + i ) = dRdAper * dAper_dU;
            }
          }
        }

        if( rowNumber >= 0 && rowNumber < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                            nodeDOF,
                                                                            dRdU.data(),
                                                                            2 * numNodesPerFace * 3 );
        }
      }
    } );
  } );
}

void
HydrofractureSolver::
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  PoroelasticSolver::applySystemSolution( dofManager, localSolution, scalingFactor, domain );

  updateDeformationForCoupling( domain );
}


real64
HydrofractureSolver::scalingForSystemSolution( DomainPartition const & domain,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution )
{
  return m_solidSolver->scalingForSystemSolution( domain,
                                                  dofManager,
                                                  localSolution );
}

void HydrofractureSolver::setNextDt( real64 const & currentDt,
                                     real64 & nextDt )
{

  if( m_numResolves[0] == 0 && m_numResolves[1] == 0 )
  {
    this->setNextDtBasedOnNewtonIter( currentDt, nextDt );
  }
  else
  {
    SolverBase & surfaceGenerator = this->getParent().getGroup< SolverBase >( "SurfaceGen" );
    nextDt = surfaceGenerator.GetTimestepRequest() < 1e99 ? surfaceGenerator.GetTimestepRequest() : currentDt;
  }

  GEOSX_LOG_LEVEL_RANK_0( 3, this->getName() << ": nextDt request is "  << nextDt );
}

void HydrofractureSolver::initializeNewFaceElements( DomainPartition const & )
{
//  m_flowSolver->
}

void HydrofractureSolver::setUpDflux_dApertureMatrix( DomainPartition & domain,
                                                      DofManager const & dofManager,
                                                      CRSMatrix< real64, globalIndex > & localMatrix )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::unique_ptr< CRSMatrix< real64, localIndex > > &
  derivativeFluxResidual_dAperture = this->getRefDerivativeFluxResidual_dAperture();

  {
    localIndex numRows = 0;
    this->template forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                                           auto const & elementSubRegion )
    {
      numRows += elementSubRegion.size();
    } );

    derivativeFluxResidual_dAperture = std::make_unique< CRSMatrix< real64, localIndex > >( numRows, numRows );
    derivativeFluxResidual_dAperture->setName( this->getName() + "/derivativeFluxResidual_dAperture" );

    derivativeFluxResidual_dAperture->reserveNonZeros( localMatrix.numNonZeros() );
    localIndex maxRowSize = -1;
    for( localIndex row = 0; row < localMatrix.numRows(); ++row )
    {
      localIndex const rowSize = localMatrix.numNonZeros( row );
      maxRowSize = maxRowSize > rowSize ? maxRowSize : rowSize;
    }
    // TODO This is way too much. The With the full system rowSize is not a good estimate for this.
    for( localIndex row = localMatrix.numRows(); row < numRows; ++row )
    {
      derivativeFluxResidual_dAperture->reserveNonZeros( row, maxRowSize );
    }
  }

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( this->getDiscretization() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn = 0; iconn < stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      for( localIndex k0 = 0; k0 < numFluxElems; ++k0 )
      {
        for( localIndex k1 = 0; k1 < numFluxElems; ++k1 )
        {
          derivativeFluxResidual_dAperture->insertNonZero( sei[iconn][k0], sei[iconn][k1], 0.0 );
        }
      }
    }
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, HydrofractureSolver, string const &, Group * const )
} /* namespace geosx */
