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
 * @file LaplaceFEM.cpp
 *
 */

// Source includes
#include "LaplaceFEM.hpp"
#include "LaplaceFEMKernels.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/TimingMacros.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteElement/FiniteElementDiscretization.hpp"
#include "finiteElement/FiniteElementDiscretizationManager.hpp"
#include "finiteElement/Kinematics.h"
#include "managers/NumericalMethodsManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{}
}  // namespace dataRepository

using namespace dataRepository;

//START_SPHINX_INCLUDE_01
LaplaceFEM::LaplaceFEM( const std::string & name, Group * const parent ) :
  SolverBase( name, parent ),
  m_fieldName( "primaryField" ),
  m_timeIntegrationOption( timeIntegrationOption::ImplicitTransient )
{
  registerWrapper< string >( laplaceFEMViewKeys.timeIntegrationOption.Key() )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "option for default time integration method" );

  registerWrapper< string >( laplaceFEMViewKeys.fieldVarName.Key(), &m_fieldName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "name of field variable" );
}
//END_SPHINX_INCLUDE_01

LaplaceFEM::~LaplaceFEM()
{
  // TODO Auto-generated destructor stub
}

//START_SPHINX_INCLUDE_02
void
LaplaceFEM::RegisterDataOnMesh( Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    NodeManager * const nodes =
      mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getNodeManager();

    nodes->registerWrapper< real64_array >( m_fieldName )
      ->setApplyDefaultValue( 0.0 )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setDescription( "Primary field variable" );
  }
}
//END_SPHINX_INCLUDE_02

//START_SPHINX_INCLUDE_03
void
LaplaceFEM::PostProcessInput()
{
  SolverBase::PostProcessInput();

  string tiOption =
    this->getReference< string >( laplaceFEMViewKeys.timeIntegrationOption );

  if( tiOption == "SteadyState" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::SteadyState;
  }
  else if( tiOption == "ImplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ImplicitTransient;
  }
  else if( tiOption == "ExplicitTransient" )
  {
    this->m_timeIntegrationOption = timeIntegrationOption::ExplicitTransient;
  }
  else
  {
    GEOSX_ERROR( "invalid time integration option" );
  }
}
//END_SPHINX_INCLUDE_03

real64
LaplaceFEM::SolverStep( real64 const & time_n,
                        real64 const & dt,
                        const int cycleNumber,
                        DomainPartition & domain )
{
  real64 dtReturn = dt;
  if( m_timeIntegrationOption == timeIntegrationOption::ExplicitTransient )
  {
    dtReturn = ExplicitStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_timeIntegrationOption == timeIntegrationOption::ImplicitTransient ||
           m_timeIntegrationOption == timeIntegrationOption::SteadyState )
  {
    dtReturn = this->LinearImplicitStep( time_n, dt, cycleNumber, domain );
  }
  return dtReturn;
}

void
LaplaceFEM::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                               real64 const & GEOSX_UNUSED_PARAM( dt ),
                               DomainPartition & domain )
{
  // Computation of the sparsity pattern
  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void
LaplaceFEM::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void
LaplaceFEM::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                       DofManager & dofManager ) const
{
  dofManager.addField( m_fieldName, DofManager::Location::Node );

  dofManager.addCoupling( m_fieldName, m_fieldName, DofManager::Connector::Elem );
}

void
LaplaceFEM::SetupSystem( DomainPartition & domain,
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

  MeshLevel * const mesh =
    domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  NodeManager const * const nodeManager = mesh->getNodeManager();
  arrayView1d< globalIndex const > const & dofIndex =
    nodeManager->getReference< globalIndex_array >( dofManager.getKey( m_fieldName ) );

  SparsityPattern< globalIndex > pattern( dofManager.numLocalDofs(),
                                          dofManager.numGlobalDofs() );
  array1d< localIndex > rowSizes( dofManager.numLocalDofs() );

  finiteElement::fillSparsity< serialPolicy, CellElementSubRegion, LaplaceFEMKernel >(
    *mesh,
    targetRegionNames(),
    this->getDiscretizationName(),
    dofIndex,
    dofManager.rankOffset(),
    pattern,
    rowSizes );

  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
}

//START_SPHINX_INCLUDE_04
void
LaplaceFEM::AssembleSystem(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  MeshLevel * const mesh =
    domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  NodeManager & nodeManager = *( mesh->getNodeManager() );

  arrayView1d< globalIndex const > const & dofIndex =
    nodeManager.getReference< array1d< globalIndex > >( dofManager.getKey( m_fieldName ) );

  finiteElement::regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                               constitutive::NullModel,
                                               CellElementSubRegion,
                                               LaplaceFEMKernel >(
    *mesh,
    targetRegionNames(),
    this->getDiscretizationName(),
    array1d< string >(),
    dofIndex,
    dofManager.rankOffset(),
    localMatrix,
    localRhs,
    m_fieldName );

  //END_SPHINX_INCLUDE_04

  if( getLogLevel() == 2 )
  {
    GEOSX_LOG_RANK_0( "After LaplaceFEM::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    //    std::cout << matrix;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    //    std::cout << rhs;
  }

  //  if( getLogLevel() >= 3 )
  //  {
  //    integer newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  //>>>>>>> develop
  //
  //      finiteElementLaunchDispatch< LaplaceFEMKernels::ImplicitKernel >( numNodesPerElement,
  //                                                                        n_q_points,
  //                                                                        dNdX,
  //                                                                        detJ,
  //                                                                        elemNodes,
  //                                                                        dofIndex,
  //                                                                        dofManager.rankOffset(),
  //                                                                        localMatrix );
  //
  //    } );
  //  }

  //END_SPHINX_INCLUDE_04
}

void
LaplaceFEM::ApplySystemSolution( DofManager const & dofManager,
                                 arrayView1d< real64 const > const & localSolution,
                                 real64 const scalingFactor,
                                 DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               m_fieldName,
                               m_fieldName,
                               scalingFactor );

  // Synchronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( m_fieldName );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );
}

void
LaplaceFEM::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  ApplyDirichletBC_implicit( time_n + dt, dofManager, domain, localMatrix, localRhs );
}

void
LaplaceFEM::SolveSystem( DofManager const & dofManager,
                         ParallelMatrix & matrix,
                         ParallelVector & rhs,
                         ParallelVector & solution )
{
  rhs.scale( -1.0 );  // TODO decide if we want this here
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void
LaplaceFEM::ApplyDirichletBC_implicit(
  real64 const time,
  DofManager const & dofManager,
  DomainPartition & domain,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::get();

  fsManager.Apply( time,
                   &domain,
                   "nodeManager",
                   m_fieldName,
                   [&]( FieldSpecificationBase const * const bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group * const targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) ) {
                     bc->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                                         parallelDevicePolicy< 32 > >(
                       targetSet,
                       time,
                       targetGroup,
                       m_fieldName,
                       dofManager.getKey( m_fieldName ),
                       dofManager.rankOffset(),
                       localMatrix,
                       localRhs );
                   } );
}

void
LaplaceFEM::ResetStateToBeginningOfStep(
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

//START_SPHINX_INCLUDE_00
REGISTER_CATALOG_ENTRY( SolverBase, LaplaceFEM, std::string const &, Group * const )
//END_SPHINX_INCLUDE_00
}  // namespace geosx
