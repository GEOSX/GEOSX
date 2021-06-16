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
 * @file PoroelasticSolver.cpp
 *
 */


#include "PoroelasticSolver.hpp"

#include "../solidMechanics/SolidMechanicsPoroElasticKernel.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/permeability/permeabilitySelector.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "finiteElement/FiniteElementDispatch.hpp"
#include "SinglePhasePoroelasticKernel.hpp"
#include "constitutive/permeability/PoroMechanicsPermeabilityKernel.hpp"


namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace PermeabilityKernels;

PoroelasticSolver::PoroelasticSolver( const string & name,
                                      Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOption( CouplingTypeOption::FIM )

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString(), &m_couplingTypeOption ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Coupling method. Valid options:\n* " + EnumStrings< CouplingTypeOption >::concat( "\n* " ) );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoroelastic;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void PoroelasticSolver::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    ElementRegionManager & elemManager = meshBody.getMeshLevel( 0 ).getElemManager();

    elemManager.forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >( [&]( ElementSubRegionBase & elementSubRegion )
    {
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString() ).
        setDescription( "Total Mean Stress" );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() ).
        setDescription( "Total Mean Stress" );
      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPerm_dDisplacementString() ).
        setDescription( "Derivative of the permeability w.r.t displacement" );
    } );
  } );
}

void PoroelasticSolver::setupDofs( DomainPartition const & domain,
                                   DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );
}

void PoroelasticSolver::setupSystem( DomainPartition & domain,
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
    createPreconditioner();
  }
}

void PoroelasticSolver::implicitStepSetup( real64 const & time_n,
                                           real64 const & dt,
                                           DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );

  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&] ( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const & totalMeanStress =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString() );
      arrayView1d< real64 > const & oldTotalMeanStress =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        oldTotalMeanStress[ei] = totalMeanStress[ei];
      } );
    } );
  }
}

void PoroelasticSolver::implicitStepComplete( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}

void PoroelasticSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  m_solidSolver->setEffectiveStress( 1 );

  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver,
    // otherwise it will never converge.
    m_flowSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
    m_solidSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
  }
}

void PoroelasticSolver::initializePostInitialConditionsPreSubGroups()
{
  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    m_flowSolver->setPoroElasticCoupling();
    // Calculate initial total mean stress
    updateDeformationForCoupling( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
  }
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

void PoroelasticSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&] ( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & oldTotalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );
    arrayView1d< real64 > const & totalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      totalMeanStress[ei] = oldTotalMeanStress[ei];
    } );
  } );
}

real64 PoroelasticSolver::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      int const cycleNumber,
                                      DomainPartition & domain )
{
  real64 dt_return = dt;
  if( m_couplingTypeOption == CouplingTypeOption::SIM_FixedStress )
  {
    dt_return = splitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == CouplingTypeOption::FIM )
  {
    setupSystem( domain,
                 m_dofManager,
                 m_localMatrix,
                 m_localRhs,
                 m_localSolution );

    implicitStepSetup( time_n, dt, domain );

    dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    implicitStepComplete( time_n, dt_return, domain );
  }
  return dt_return;
}

void PoroelasticSolver::updateDeformationForCoupling( DomainPartition & domain )
{

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u = nodeManager.totalDisplacement();

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                  localIndex const,
                                                                  localIndex const,
                                                                  ElementRegionBase & elemRegion,
                                                                  CellElementSubRegion & elementSubRegion )
  {
    string const & solidName = m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex( elemRegion.getName() )];
    SolidBase const & solid = getConstitutiveModel< SolidBase >( elementSubRegion, solidName );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();

    arrayView1d< real64 > const &
    totalMeanStress = elementSubRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    arrayView1d< real64 > const &
    oldTotalMeanStress = elementSubRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );

    arrayView1d< real64 const > const &
    pres = elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );

    arrayView1d< real64 const > const &
    dPres = elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );

    arrayView1d< real64 > const &
    poro = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityString() );

    arrayView1d< real64 const > const &
    poroOld = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::porosityOldString() );

    arrayView1d< real64 const > const &
    volume = elementSubRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString() );

    arrayView1d< real64 > const &
    dVol = elementSubRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaVolumeString() );

    arrayView1d< real64 const > const & bulkModulus = solid.getReference< array1d< real64 > >( ElasticIsotropic::viewKeyStruct::bulkModulusString() );

    real64 const biotCoefficient = solid.getReference< real64 >( "BiotCoefficient" );

    arrayView3d< real64 const, solid::STRESS_USD > const & stress = solid.getStress();


    localIndex const numNodesPerElement = elemsToNodes.size( 1 );
    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_solidSolver->getDiscretizationName() );
    localIndex const numQuadraturePoints = fe.getNumQuadraturePoints();

    forAll< parallelDevicePolicy< 32 > >( elementSubRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      real64 effectiveMeanStress = 0.0;
      for( localIndex q=0; q<numQuadraturePoints; ++q )
      {
        effectiveMeanStress += ( stress( ei, q, 0 ) + stress( ei, q, 1 ) + stress( ei, q, 2 ) );
      }
      effectiveMeanStress /= ( 3 * numQuadraturePoints );

      totalMeanStress[ei] = effectiveMeanStress - biotCoefficient * (pres[ei] + dPres[ei]);

      poro[ei] = poroOld[ei] + (biotCoefficient - poroOld[ei]) / bulkModulus[ei]
                 * (totalMeanStress[ei] - oldTotalMeanStress[ei] + dPres[ei]);

      // update element volume
      real64 Xlocal[ElementRegionManager::maxNumNodesPerElem][3];
      for( localIndex a = 0; a < numNodesPerElement; ++a )
      {
        LvArray::tensorOps::copy< 3 >( Xlocal[a], X[elemsToNodes[ei][a]] );
        LvArray::tensorOps::add< 3 >( Xlocal[a], u[elemsToNodes[ei][a]] );
      }

      dVol[ei] = computationalGeometry::HexVolume( Xlocal ) - volume[ei];
    } );
  } );
}

void PoroelasticSolver::assembleSystem( real64 const time_n,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );

  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

//  m_solidSolver->resetStressToBeginningOfStep( domain );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  PoroelasticKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                              pDofKey,
                                                              dofManager.rankOffset(),
                                                              localMatrix,
                                                              localRhs,
                                                              gravityVectorData,
                                                              m_flowSolver->fluidModelNames() );

  // Cell-based contributions
  m_solidSolver->getMaxForce() =
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::SolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            m_solidSolver->solidMaterialNames(),
                                                            kernelFactory );

  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              " " );

}

void PoroelasticSolver::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  m_solidSolver->applyBoundaryConditions( time_n, dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64 PoroelasticSolver::calculateResidualNorm( DomainPartition const & domain,
                                                 DofManager const & dofManager,
                                                 arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_solidSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rsolid, Rfluid ) = ( %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void PoroelasticSolver::createPreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == LinearSolverParameters::PreconditionerType::block )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >( BlockShapeOption::UpperTriangular,
                                                                           SchurComplementOption::RowsumDiagonalProbing,
                                                                           BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner( m_solidSolver->getLinearSolverParameters() );
    precond->setupBlock( 0,
                         { { keys::TotalDisplacement, { 3, true } } },
                         std::make_unique< SeparateComponentPreconditioner< LAInterface > >( 3, std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner( m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { SinglePhaseBase::viewKeyStruct::pressureString(), { 1, true } } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void PoroelasticSolver::solveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void PoroelasticSolver::applySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void PoroelasticSolver::updateState( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();

  this->template forTargetSubRegions< CellElementSubRegion >( mesh, [&] ( localIndex const targetIndex,
                                                                          auto & subRegion )
  {
    updatePermeability( nodeManager, subRegion, targetIndex );
    m_flowSolver->updateFluidState( subRegion, targetIndex );
  } );
}

void PoroelasticSolver::updatePermeability( NodeManager const & nodeManager,
                                            CellElementSubRegion & subRegion,
                                            localIndex const targetIndex ) const
{
  PermeabilityBase & perm =
    getConstitutiveModel< PermeabilityBase >( subRegion, m_flowSolver->permeabilityModelNames()[targetIndex] );

  string const & solidName = m_solidSolver->solidMaterialNames()[targetIndex];

  constitutive::constitutiveUpdatePassThru( perm, [&] ( auto & castedPerm )
  {
    finiteElement::FiniteElementBase &
    subRegionFE = subRegion.template getReference< finiteElement::FiniteElementBase >( m_solidSolver->getDiscretizationName() );

    finiteElement::dispatch3D( subRegionFE,
                               [&solidName,
                                &subRegion,
                                &castedPerm,
                                &nodeManager] ( auto const finiteElement )
    {

      using FE_TYPE = TYPEOFREF( finiteElement );


      using KERNEL_TYPE = PoroMechanicsPermeabilityKernel< CellElementSubRegion,
                                                           FE_TYPE >;

      KERNEL_TYPE permKernel( subRegion,
                              finiteElement,
                              nodeManager );

      typename TYPEOFREF( castedPerm ) ::KernelWrapper permWrapper = castedPerm.createKernelWrapper();

      arrayView3d< real64 > const & dPerm_dDisplacement =
        subRegion.template getReference< array3d< real64 > >( viewKeyStruct::dPerm_dDisplacementString() );

      arrayView1d< real64 const > const & pressure =
        subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );

      arrayView1d< real64 const > const & deltaPressure =
        subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );

      SolidBase & solidModel =
        getConstitutiveModel< SolidBase >( subRegion,
                                           solidName );

      arrayView2d< real64 const > const & porosity = solidModel.getPorosity();
      array2d< real64 > dPorosity_dVolStrain;
      dPorosity_dVolStrain.setValues< serialPolicy >( 0.0 );

      KERNEL_TYPE::template launch< parallelDevicePolicy<>, KERNEL_TYPE >( subRegion.size(),
                                                                           permKernel,
                                                                           permWrapper,
                                                                           pressure,
                                                                           deltaPressure,
                                                                           porosity,
                                                                           dPorosity_dVolStrain.toViewConst(),
                                                                           dPerm_dDisplacement );
    } );
  } );

}

real64 PoroelasticSolver::splitOperatorStep( real64 const & time_n,
                                             real64 const & dt,
                                             integer const cycleNumber,
                                             DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  m_flowSolver->setupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getLocalRhs(),
                             m_flowSolver->getLocalSolution() );

  m_solidSolver->setupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getLocalMatrix(),
                              m_solidSolver->getLocalRhs(),
                              m_solidSolver->getLocalSolution() );

  implicitStepSetup( time_n, dt, domain );

  int iter = 0;
  while( iter < m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all child solvers if any of them has been reset
      resetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = m_flowSolver->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( m_flowSolver->getNonlinearSolverParameters().m_numNewtonIterations == 0 && iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", MechanicsSolver: " );

    //m_solidSolver->resetStressToBeginningOfStep( domain );
    dtReturnTemporary = m_solidSolver->nonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    if( m_solidSolver->getNonlinearSolverParameters().m_numNewtonIterations > 0 )
    {
      updateDeformationForCoupling( domain );
    }
    ++iter;
  }

  implicitStepComplete( time_n, dt, domain );

  return dtReturn;
}


REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolver, string const &, Group * const )

} /* namespace geosx */
