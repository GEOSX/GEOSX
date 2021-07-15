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
 * @file SinglePhasePoromechanicsSolver.cpp
 *
 */


#include "SinglePhasePoromechanicsSolver.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PorousSolid.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
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

#include "SinglePhasePoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhasePoromechanicsSolver::SinglePhasePoromechanicsSolver( const string & name,
                                                                Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName()

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poromechanics solver" );

  registerWrapper( viewKeyStruct::porousMaterialNamesString(), &m_porousMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void SinglePhasePoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                                DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );
}

void SinglePhasePoromechanicsSolver::setupSystem( DomainPartition & domain,
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

void SinglePhasePoromechanicsSolver::implicitStepSetup( real64 const & time_n,
                                                        real64 const & dt,
                                                        DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolver::implicitStepComplete( real64 const & time_n,
                                                           real64 const & dt,
                                                           DomainPartition & domain )
{
  m_solidSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    ConstitutiveBase const & porousMaterial = getConstitutiveModel< ConstitutiveBase >( subRegion, porousMaterialNames()[targetIndex] );
    porousMaterial.saveConvergedState();
  } );
}

void SinglePhasePoromechanicsSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

void SinglePhasePoromechanicsSolver::initializePostInitialConditionsPreSubGroups()
{
  if( m_flowSolver->getLinearSolverParameters().mgr.strategy == LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM )
  {
    m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::hybridSinglePhasePoromechanics;
  }
}

SinglePhasePoromechanicsSolver::~SinglePhasePoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void SinglePhasePoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );
}

real64 SinglePhasePoromechanicsSolver::solverStep( real64 const & time_n,
                                                   real64 const & dt,
                                                   int const cycleNumber,
                                                   DomainPartition & domain )
{
  real64 dt_return = dt;

  setupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_localRhs,
               m_localSolution );

  implicitStepSetup( time_n, dt, domain );

  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void SinglePhasePoromechanicsSolver::assembleSystem( real64 const time_n,
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

  PoromechanicsKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
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
                                    constitutive::PorousSolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            this->getDiscretizationName(),
                                                            porousMaterialNames(),
                                                            kernelFactory );

  // Face-based contributions
  m_flowSolver->assembleFluxTerms( time_n, dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );

}

void SinglePhasePoromechanicsSolver::assembleCouplingTerms( DomainPartition const & domain,
                                                            DofManager const & dofManager,
                                                            CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                            arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();

  string const uDofKey = dofManager.getKey( keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & uDofNumber = nodeManager.getReference< globalIndex_array >( uDofKey );

  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incr_disp = nodeManager.incrementalDisplacement();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  // begin subregion loop
  forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                  localIndex const,
                                                                  localIndex const,
                                                                  ElementRegionBase const & region,
                                                                  CellElementSubRegion const & elementSubRegion )
  {
    string const & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( elementSubRegion, fluidName );

    string const & solidName = m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex( region.getName() )];
    SolidBase const & solid = getConstitutiveModel< SolidBase >( elementSubRegion, solidName );

    arrayView4d< real64 const > const & dNdX = elementSubRegion.dNdX();

    arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

    arrayView1d< globalIndex const > const & pDofNumber = elementSubRegion.getReference< globalIndex_array >( pDofKey );

    arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes = elementSubRegion.nodeList();
    localIndex const numNodesPerElement = elemsToNodes.size( 1 );

    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_solidSolver->getDiscretizationName() );
    localIndex const numQuadraturePoints = fe.getNumQuadraturePoints();

    real64 const biotCoefficient = solid.getReference< real64 >( "BiotCoefficient" );

    arrayView2d< real64 const > const & density = fluid.density();

    int dim = 3;
    localIndex constexpr maxNumUDof = 24;   // TODO: assuming linear HEX at most for the moment
    localIndex constexpr maxNumPDof = 1;   // TODO: assuming piecewise constant (P0) only for the moment
    localIndex const nUDof = dim * numNodesPerElement;
    localIndex const nPDof = m_flowSolver->numDofPerCell();
    GEOSX_ERROR_IF_GT( nPDof, maxNumPDof );

    forAll< parallelDevicePolicy< 32 > >( elementSubRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      stackArray2d< real64, maxNumUDof * maxNumPDof > dRsdP( nUDof, nPDof );
      stackArray2d< real64, maxNumUDof * maxNumPDof > dRfdU( nPDof, nUDof );
      stackArray1d< real64, maxNumPDof > Rf( nPDof );

      for( integer q = 0; q < numQuadraturePoints; ++q )
      {
        const real64 detJq = detJ[k][q];

        for( integer a = 0; a < numNodesPerElement; ++a )
        {

          dRsdP( a * dim + 0, 0 ) += biotCoefficient * dNdX[k][q][a][0] * detJq;
          dRsdP( a * dim + 1, 0 ) += biotCoefficient * dNdX[k][q][a][1] * detJq;
          dRsdP( a * dim + 2, 0 ) += biotCoefficient * dNdX[k][q][a][2] * detJq;
          dRfdU( 0, a * dim + 0 ) += density[k][0] * biotCoefficient * dNdX[k][q][a][0] * detJq;
          dRfdU( 0, a * dim + 1 ) += density[k][0] * biotCoefficient * dNdX[k][q][a][1] * detJq;
          dRfdU( 0, a * dim + 2 ) += density[k][0] * biotCoefficient * dNdX[k][q][a][2] * detJq;

          localIndex localNodeIndex = elemsToNodes[k][a];

          real64 Rf_tmp = dNdX[k][q][a][0] * incr_disp[localNodeIndex][0]
                          + dNdX[k][q][a][1] * incr_disp[localNodeIndex][1]
                          + dNdX[k][q][a][2] * incr_disp[localNodeIndex][2];
          Rf_tmp *= density[k][0] * biotCoefficient * detJq;
          Rf[0] += Rf_tmp;
        }
      }

      stackArray1d< globalIndex, maxNumUDof > elementULocalDofIndex( nUDof );
      stackArray1d< globalIndex, maxNumPDof > elementPLocalDofIndex( nPDof );

      // Get dof local to global mapping
      for( localIndex a = 0; a < numNodesPerElement; ++a )
      {
        for( int i = 0; i < dim; ++i )
        {
          elementULocalDofIndex[a * dim + i] = uDofNumber[elemsToNodes[k][a]] + i;
        }
      }
      for( localIndex i = 0; i < nPDof; ++i )
      {
        elementPLocalDofIndex[i] = pDofNumber[k] + i;
      }

      for( localIndex i = 0; i < nUDof; ++i )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( elementULocalDofIndex[ i ] - rankOffset );
        if( dof < 0 || dof >= localMatrix.numRows() )
          continue;
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                          elementPLocalDofIndex.data(),
                                                                          dRsdP[i].dataIfContiguous(),
                                                                          nPDof );
      }
      for( localIndex i = 0; i < nPDof; ++i )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( elementPLocalDofIndex[ i ] - rankOffset );
        if( dof < 0 || dof >= localMatrix.numRows() )
          continue;
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                          elementULocalDofIndex.data(),
                                                                          dRfdU[i].dataIfContiguous(),
                                                                          nUDof );

        RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[ dof ], Rf[i] );
      }
    } );
  } );
}

void SinglePhasePoromechanicsSolver::applyBoundaryConditions( real64 const time_n,
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

real64 SinglePhasePoromechanicsSolver::calculateResidualNorm( DomainPartition const & domain,
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

void SinglePhasePoromechanicsSolver::createPreconditioner()
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

void SinglePhasePoromechanicsSolver::solveSystem( DofManager const & dofManager,
                                                  ParallelMatrix & matrix,
                                                  ParallelVector & rhs,
                                                  ParallelVector & solution )
{
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void SinglePhasePoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                          arrayView1d< real64 const > const & localSolution,
                                                          real64 const scalingFactor,
                                                          DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
