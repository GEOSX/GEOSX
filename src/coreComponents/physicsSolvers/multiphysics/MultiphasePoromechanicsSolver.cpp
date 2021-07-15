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
 * @file MultiphasePoroelasticSolver.cpp
 *
 */

#include "MultiphasePoromechanicsSolver.hpp"

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
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

#include "physicsSolvers/multiphysics/MultiphasePoromechanicsKernel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

MultiphasePoromechanicsSolver::MultiphasePoromechanicsSolver( const string & name,
                                                              Group * const parent ):
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName()

{
  registerWrapper( viewKeyStruct::solidSolverNameString(), &m_solidSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::porousMaterialNamesString(), &m_porousMaterialNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "The name of the material that should be used in the constitutive updates" );

  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::multiphasePoromechanics;
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName = keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void MultiphasePoromechanicsSolver::setupDofs( DomainPartition const & domain,
                                               DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString(),
                          DofManager::Connector::Elem );
}

void MultiphasePoromechanicsSolver::setupSystem( DomainPartition & domain,
                                                 DofManager & dofManager,
                                                 CRSMatrix< real64, globalIndex > & localMatrix,
                                                 array1d< real64 > & localRhs,
                                                 array1d< real64 > & localSolution,
                                                 bool const setSparsity )
{
//  if( m_precond )
//  {
//    m_precond->clear();
//  }

  // setup monolithic coupled system
  SolverBase::setupSystem( domain, dofManager, localMatrix, localRhs, localSolution, setSparsity );

//  if( !m_precond && m_linearSolverParameters.get().solverType != LinearSolverParameters::SolverType::direct )
//  {
//    createPreconditioner();
//  }
}

void MultiphasePoromechanicsSolver::implicitStepSetup( real64 const & time_n,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_solidSolver->implicitStepSetup( time_n, dt, domain );
}

void MultiphasePoromechanicsSolver::implicitStepComplete( real64 const & time_n,
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

void MultiphasePoromechanicsSolver::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< CompositionalMultiphaseBase >( m_flowSolverName );
  m_solidSolver = &this->getParent().getGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );
}

MultiphasePoromechanicsSolver::~MultiphasePoromechanicsSolver()
{
  // TODO Auto-generated destructor stub
}

void MultiphasePoromechanicsSolver::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_solidSolver->resetStateToBeginningOfStep( domain );
}

real64 MultiphasePoromechanicsSolver::solverStep( real64 const & time_n,
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

void MultiphasePoromechanicsSolver::assembleSystem( real64 const time_n,
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( time_n );

  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();

  string const displacementDofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & displacementDofNumber = nodeManager.getReference< globalIndex_array >( displacementDofKey );

  string const flowDofKey = dofManager.getKey( CompositionalMultiphaseBase::viewKeyStruct::elemDofFieldString() );

  localIndex const numComponents = m_flowSolver->numFluidComponents();
  localIndex const numPhases = m_flowSolver->numFluidPhases();

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  PoromechanicsKernels::MultiphaseKernelFactory kernelFactory( displacementDofNumber,
                                                               flowDofKey,
                                                               dofManager.rankOffset(),
                                                               gravityVectorData,
                                                               numComponents,
                                                               numPhases,
                                                               m_flowSolver->fluidModelNames(),
                                                               localMatrix,
                                                               localRhs );

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
  m_flowSolver->assembleFluxTerms( dt,
                                   domain,
                                   dofManager,
                                   localMatrix,
                                   localRhs );
}

void MultiphasePoromechanicsSolver::applyBoundaryConditions( real64 const time_n,
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

real64 MultiphasePoromechanicsSolver::calculateResidualNorm( DomainPartition const & domain,
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

void MultiphasePoromechanicsSolver::solveSystem( DofManager const & dofManager,
                                                 ParallelMatrix & matrix,
                                                 ParallelVector & rhs,
                                                 ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void MultiphasePoromechanicsSolver::applySystemSolution( DofManager const & dofManager,
                                                         arrayView1d< real64 const > const & localSolution,
                                                         real64 const scalingFactor,
                                                         DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

REGISTER_CATALOG_ENTRY( SolverBase, MultiphasePoromechanicsSolver, string const &, Group * const )

} /* namespace geosx */
