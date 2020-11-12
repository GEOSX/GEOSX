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
 * @file FlowReactiveTransportSolver.cpp
 *
 */


#include "FlowReactiveTransportSolver.hpp"

#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/ReactiveTransport.hpp"
#include "physicsSolvers/fluidFlow/GeochemicalModel.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowReactiveTransportSolver::FlowReactiveTransportSolver( const std::string & name, Group * const parent ):
  SolverBase( name, parent ),
  m_flowSolverName(),
  m_reactiveTransportSolverName(),
  m_geochemicalModelName(),
  m_flowSolver{},
  m_reactiveTransportSolver{},
  m_geochemicalModel{}
{
  registerWrapper( viewKeyStruct::geochemicalModelNameString, &m_geochemicalModelName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the geochemical model to use in the flowReactiveTransport solver" );

  registerWrapper( viewKeyStruct::reactiveTransportSolverNameString, &m_reactiveTransportSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the reactive transport solver to use in the flowReactiveTransport solver" );

  registerWrapper( viewKeyStruct::flowSolverNameString, &m_flowSolverName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of the flow solver to use in the flowReactiveTransport solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString )->
    setInputFlag( InputFlags::FALSE );
}

void FlowReactiveTransportSolver::RegisterDataOnMesh( dataRepository::Group * const )
{}

void FlowReactiveTransportSolver::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_geochemicalModel = this->getParent()->GetGroup< GeochemicalModel >( m_geochemicalModelName );
  GEOSX_ERROR_IF( m_geochemicalModel == nullptr, "Invalid geochemical model name: " << m_geochemicalModelName );

  m_reactiveTransportSolver = this->getParent()->GetGroup< ReactiveTransport >( m_reactiveTransportSolverName );
  GEOSX_ERROR_IF( m_reactiveTransportSolver == nullptr, "Invalid reactive transport solver name: " << m_reactiveTransportSolverName );

  m_flowSolver = this->getParent()->GetGroup< FlowSolverBase >( m_flowSolverName );
  GEOSX_ERROR_IF( m_flowSolver == nullptr, "Invalid flow solver name: " << m_flowSolverName );
}

FlowReactiveTransportSolver::~FlowReactiveTransportSolver() = default;

void FlowReactiveTransportSolver::PreStepUpdate( real64 const & time_n,
                                                 real64 const & dt,
                                                 int const cycleNumber,
                                                 DomainPartition & domain )
{


  if( time_n <= 0.0 )
  {
    if( 0 )
      m_geochemicalModel->SolverStep( time_n, dt, cycleNumber, domain );

    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();
    boundaryConditionManager.ApplyInitialConditions( &domain );

  }

  m_flowSolver->SetupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getLocalRhs(),
                             m_flowSolver->getLocalSolution() );


  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );

  m_reactiveTransportSolver->SetupSystem( domain,
                                          m_reactiveTransportSolver->getDofManager(),
                                          m_reactiveTransportSolver->getLocalMatrix(),
                                          m_reactiveTransportSolver->getLocalRhs(),
                                          m_reactiveTransportSolver->getLocalSolution() );

  m_reactiveTransportSolver->ImplicitStepSetup( time_n, dt, domain );

}

void FlowReactiveTransportSolver::PostStepUpdate( real64 const & time_n,
                                                  real64 const & dt,
                                                  int const cycleNumber,
                                                  DomainPartition & domain )
{
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
  m_reactiveTransportSolver->ImplicitStepComplete( time_n, dt, domain );
  if( 1 )
    m_geochemicalModel->SolverStep( time_n, dt, cycleNumber, domain );

}

void FlowReactiveTransportSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{

  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_reactiveTransportSolver->ResetStateToBeginningOfStep( domain );

}


real64 FlowReactiveTransportSolver::SolverStep( real64 const & time_n,
                                                real64 const & dt,
                                                int const cycleNumber,
                                                DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  PreStepUpdate( time_n, dt, cycleNumber, domain );

  int iter = 0;
  while( iter < this->m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all slave solvers if any of them has been reset
      ResetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", FlowSolver: " );

    dtReturnTemporary = m_flowSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    /*
       NonlinearSolverParameters const & fluidNonLinearParams = m_flowSolver->getNonlinearSolverParameters();
       if( fluidNonLinearParams.m_numNewtonIterations <= this->m_nonlinearSolverParameters.m_minIterNewton && iter > 0 )
       {
       GEOSX_LOG_LEVEL_RANK_0( 1, "***** The iterative coupling has converged in " << iter  << " iterations! *****\n" );
       break;
       }
     */

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter+1  << ", Reactive Transport Solver: " );

    dtReturnTemporary = m_reactiveTransportSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    ++iter;
  }

  PostStepUpdate( time_n, dtReturn, cycleNumber, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase, FlowReactiveTransportSolver, std::string const &, Group * const )

} /* namespace geosx */
