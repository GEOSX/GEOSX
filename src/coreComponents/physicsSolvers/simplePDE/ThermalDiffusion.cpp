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

#include "ThermalDiffusion.hpp"
#include "ThermalDiffusionKernels.hpp"

namespace geosx
{
using namespace dataRepository;

ThermalDiffusion::ThermalDiffusion( string const & name,
                                    Group * const parent ):
  SolverBase( name, parent ),
  m_thermalDiffusion( 1.0 )
{
  registerWrapper( viewKeyStruct::thermalDiffusionString(), &m_thermalDiffusion ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Thermal diffusion coefficient" );
}

void ThermalDiffusion::registerDataOnMesh( Group & meshBodies )
{
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    NodeManager & nodes = meshBody.getMeshLevel( 0 ).getNodeManager();

    nodes.registerWrapper< array1d< real64 > >( keys::Temperature ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_0 ).
      setDescription( "Temperature field" );

    nodes.registerWrapper< array1d< real64 > >( viewKeyStruct::newDeltaTemperatureString() ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_1 ).
      setDescription( "Incremental temperature field" );

    nodes.registerWrapper< array1d< real64 > >( keys::IncrementalTemperature ).
      setApplyDefaultValue( 0.0 ).
      setPlotLevel( PlotLevel::LEVEL_1 ).
      setDescription( "Old incremental temperature field" );
  } );
}

real64 ThermalDiffusion::solverStep( real64 const & time_n,
                                     real64 const & dt,
                                     const integer cycleNumber,
                                     DomainPartition & domain )
{
  this->setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
  return this->nonlinearImplicitStep( time_n, dt, cycleNumber, domain );
}

real64 ThermalDiffusion::calculateResidualNorm( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                                DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  return 0;
}

void ThermalDiffusion::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                  DofManager & dofManager ) const
{
  dofManager.addField( keys::Temperature,
                       DofManager::Location::Node,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( keys::Temperature,
                          keys::Temperature,
                          DofManager::Connector::Elem );
}

void ThermalDiffusion::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time ),
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NodeManager & nodeManager = mesh.getNodeManager();
  string const dofKey = dofManager.getKey( keys::Temperature );
  arrayView1d< globalIndex const > const &
  dofIndex =  nodeManager.getReference< array1d< globalIndex > >( dofKey );


  ThermalDiffusionKernelFactory kernelFactory( dofIndex,
                                               dofManager.rankOffset(),
                                               localMatrix,
                                               localRhs,
                                               m_thermalDiffusion,
                                               dt );

  finiteElement::
    regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                  constitutive::NullModel,
                                  CellElementSubRegion >( mesh,
                                                          targetRegionNames(),
                                                          this->getDiscretizationName(),
                                                          arrayView1d< string const >(),
                                                          kernelFactory );
}

void ThermalDiffusion::setupSystem( DomainPartition & domain,
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
  string const dofKey = dofManager.getKey( keys::Temperature );
  arrayView1d< globalIndex const > const &
  dofIndex = nodeManager.getReference< globalIndex_array >( dofKey );

  SparsityPattern< globalIndex > sparsityPattern( dofManager.numLocalDofs(),
                                                  dofManager.numGlobalDofs(),
                                                  8*8*3 );

  finiteElement::fillSparsity< CellElementSubRegion,
                               ThermalDiffusionKernel >( mesh,
                                                         targetRegionNames(),
                                                         this->getDiscretizationName(),
                                                         dofIndex,
                                                         dofManager.rankOffset(),
                                                         sparsityPattern );

  sparsityPattern.compress();
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( sparsityPattern ) );
}

void ThermalDiffusion::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                          real64 const & GEOSX_UNUSED_PARAM( dt ),
                                          DomainPartition & domain )
{
  setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
}

void ThermalDiffusion::applySystemSolution( DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localSolution,
                                            real64 const scalingFactor,
                                            DomainPartition & domain )
{
  dofManager.addVectorToField( localSolution,
                               keys::Temperature,
                               keys::Temperature,
                               scalingFactor );

  dofManager.addVectorToField( localSolution,
                               keys::Temperature,
                               viewKeyStruct::newDeltaTemperatureString(),
                               scalingFactor );

  // Synchronize ghost nodes
  std::map< string, string_array > fieldNames;
  fieldNames["node"].emplace_back( keys::Temperature );
  fieldNames["node"].emplace_back( viewKeyStruct::newDeltaTemperatureString() );

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames,
                                                              domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                              domain.getNeighbors(),
                                                              true );
}


void ThermalDiffusion::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                             real64 const & GEOSX_UNUSED_PARAM( dt ),
                                             DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ThermalDiffusion::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager & nodeManager = mesh.getNodeManager();

  arrayView1d< real64 > const &
  newDeltaTemperature  = nodeManager.template getReference< array1d< real64 > >( viewKeyStruct::newDeltaTemperatureString() );

  arrayView1d< real64 > const &
  oldDeltaTemperature  = nodeManager.template getReference< array1d< real64 > >( keys::IncrementalTemperature );

  forAll< parallelDevicePolicy< 32 > >( nodeManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    oldDeltaTemperature( a ) = newDeltaTemperature( a );
    newDeltaTemperature( a ) = 0.0;
  } );
}

void ThermalDiffusion::solveSystem( DofManager const & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  rhs.scale( -1.0 );
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void ThermalDiffusion::applyBoundaryConditions( real64 const time_n,
                                                real64 const dt,
                                                DomainPartition & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs )
{
  FieldSpecificationManager const & fsManager = FieldSpecificationManager::getInstance();

  fsManager.apply( time_n + dt,
                   domain,
                   "nodeManager",
                   keys::Temperature,
                   [&]( FieldSpecificationBase const & bc,
                        string const &,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & targetGroup,
                        string const & GEOSX_UNUSED_PARAM( fieldName ) )
  {
    bc.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy< 32 > >( targetSet,
                                                                     time_n + dt,
                                                                     targetGroup,
                                                                     keys::Temperature,
                                                                     dofManager.getKey( keys::Temperature ),
                                                                     dofManager.rankOffset(),
                                                                     localMatrix,
                                                                     localRhs );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, ThermalDiffusion, string const &, Group * const )

} /* namespace geosx */
