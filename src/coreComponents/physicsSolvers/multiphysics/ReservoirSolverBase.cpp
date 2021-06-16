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
 * @file ReservoirSolverBase.cpp
 *
 */

#include "ReservoirSolverBase.hpp"

#include "common/TimingMacros.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/wells/WellSolverBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

ReservoirSolverBase::ReservoirSolverBase( const string & name,
                                          Group * const parent ):
  SolverBase( name, parent ),
  m_flowSolverName(),
  m_wellSolverName()
{
  registerWrapper( viewKeyStruct::flowSolverNameString(), &m_flowSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the flow solver to use in the reservoir-well system solver" );

  registerWrapper( viewKeyStruct::wellSolverNameString(), &m_wellSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the well solver to use in the reservoir-well system solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );

}

ReservoirSolverBase::~ReservoirSolverBase()
{}

void ReservoirSolverBase::postProcessInput()
{
  SolverBase::postProcessInput();

  m_flowSolver = &this->getParent().getGroup< FlowSolverBase >( m_flowSolverName );
  m_wellSolver = &this->getParent().getGroup< WellSolverBase >( m_wellSolverName );

  m_wellSolver->setFlowSolverName( m_flowSolverName );
  m_flowSolver->setReservoirWellsCoupling();
}

void ReservoirSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups( );

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  // loop over the wells
  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion & subRegion )
  {
    // get the string to access the permeability
    string const permeabilityKey = PermeabilityBase::viewKeyStruct::permeabilityString();

    PerforationData * const perforationData = subRegion.getPerforationData();

    // compute the Peaceman index (if not read from XML)
    perforationData->computeWellTransmissibility( meshLevel, subRegion,
                                                  permeabilityKey,
                                                  m_flowSolver->targetRegionNames(),
                                                  m_flowSolver->permeabilityModelNames() );
  } );

  // bind the stored reservoir views to the current domain
  resetViews( domain );
}


real64 ReservoirSolverBase::solverStep( real64 const & time_n,
                                        real64 const & dt,
                                        int const cycleNumber,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return = dt;

  // setup the coupled linear system
  setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

  // setup reservoir and well systems
  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // complete time step in reservoir and well systems
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void ReservoirSolverBase::setupDofs( DomainPartition const & domain,
                                     DofManager & dofManager ) const
{
  m_flowSolver->setupDofs( domain, dofManager );
  m_wellSolver->setupDofs( domain, dofManager );
  // TODO: add coupling when dofManager can support perforation connectors
}

void ReservoirSolverBase::addCouplingNumNonzeros( DomainPartition & domain,
                                                  DofManager & dofManager,
                                                  arrayView1d< localIndex > const & rowLengths ) const
{
  localIndex const resNDOF = m_wellSolver->numDofPerResElement();
  localIndex const wellNDOF = m_wellSolver->numDofPerWellElement();

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
  string const resDofKey = dofManager.getKey( m_wellSolver->resElementDofName() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resElemDofNumber =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > const & resElemGhostRank =
    elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

  globalIndex const rankOffset = dofManager.rankOffset();
  elemManager.forElementSubRegions< WellElementSubRegion >( [&]( WellElementSubRegion const & subRegion )
  {
    PerforationData const * const perforationData = subRegion.getPerforationData();

    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion = perforationData->getMeshElements().m_toElementRegion;
    arrayView1d< localIndex const > const & resElementSubRegion = perforationData->getMeshElements().m_toElementSubRegion;
    arrayView1d< localIndex const > const & resElementIndex = perforationData->getMeshElements().m_toElementIndex;

    // Loop over perforations and increase row lengths for reservoir and well elements accordingly
    forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
    {
      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      if( resElemGhostRank[er][esr][ei] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( resElemDofNumber[er][esr][ei] - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( rowLengths.size(), localRow + resNDOF );

        for( localIndex idof = 0; idof < resNDOF; ++idof )
        {
          rowLengths[localRow + idof] += wellNDOF;
        }
      }

      if( wellElemGhostRank[iwelem] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( wellElemDofNumber[iwelem] - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( rowLengths.size(), localRow + wellNDOF );

        for( localIndex idof = 0; idof < wellNDOF; ++idof )
        {
          rowLengths[localRow + idof] += resNDOF;
        }
      }
    } );
  } );
}

void ReservoirSolverBase::setupSystem( DomainPartition & domain,
                                       DofManager & dofManager,
                                       CRSMatrix< real64, globalIndex > & localMatrix,
                                       array1d< real64 > & localRhs,
                                       array1d< real64 > & localSolution,
                                       bool const )
{
  GEOSX_MARK_FUNCTION;

  dofManager.setMesh( domain.getMeshBody( 0 ).getMeshLevel( 0 ) );

  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without reservoir-well coupling
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling on perforations
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localRhs.resize( localMatrix.numRows() );
  localSolution.resize( localMatrix.numRows() );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );
}


void ReservoirSolverBase::implicitStepSetup( real64 const & time_n,
                                             real64 const & dt,
                                             DomainPartition & domain )
{
  // setup the individual solvers
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_wellSolver->implicitStepSetup( time_n, dt, domain );
}


void ReservoirSolverBase::assembleSystem( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  // assemble J_RR (excluding perforation rates)
  m_flowSolver->assembleSystem( time_n, dt,
                                domain,
                                dofManager,
                                localMatrix,
                                localRhs );

  /*
   * This redundant call to UpdateStateAll is here to make sure that we compute the
   * perforation rates AFTER the reservoir phase compositions have been moved to device.
   *
   * An issue with ElementViewAccessors is that if the outer arrays are already on device,
   * but an inner array gets touched and updated on host, capturing outer arrays in a device kernel
   * DOES NOT call move() on the inner array (see implementation of NewChaiBuffer::moveNested()).
   * Here we force the move by launching a dummy kernel.
   *
   * If the perforation rates are computed BEFORE the reservoir phase compositions have been
   * moved to device, the calculation is wrong. the problem should go away when fluid updates
   * are executed on device.
   */
  m_wellSolver->updateState( domain );

  // assemble J_WW (excluding perforation rates)
  m_wellSolver->assembleSystem( time_n, dt,
                                domain,
                                dofManager,
                                localMatrix,
                                localRhs );

  // assemble perforation rates in J_WR, J_RW, J_RR and J_WW
  assembleCouplingTerms( time_n, dt,
                         domain,
                         dofManager,
                         localMatrix,
                         localRhs );
}


void ReservoirSolverBase::applyBoundaryConditions( real64 const time_n,
                                                   real64 const dt,
                                                   DomainPartition & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  m_flowSolver->applyBoundaryConditions( time_n,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
  // no boundary conditions for wells
}

real64 ReservoirSolverBase::calculateResidualNorm( DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   arrayView1d< real64 const > const & localRhs )
{
  // compute norm of reservoir equations residuals
  real64 const reservoirResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );
  // compute norm of well equations residuals
  real64 const wellResidualNorm      = m_wellSolver->calculateResidualNorm( domain, dofManager, localRhs );

  return sqrt( reservoirResidualNorm * reservoirResidualNorm + wellResidualNorm * wellResidualNorm );
}

void ReservoirSolverBase::solveSystem( DofManager const & dofManager,
                                       ParallelMatrix & matrix,
                                       ParallelVector & rhs,
                                       ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();
  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

bool ReservoirSolverBase::checkSystemSolution( DomainPartition const & domain,
                                               DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor )
{
  bool const validReservoirSolution = m_flowSolver->checkSystemSolution( domain, dofManager, localSolution, scalingFactor );
  bool const validWellSolution      = m_wellSolver->checkSystemSolution( domain, dofManager, localSolution, scalingFactor );

  return ( validReservoirSolution && validWellSolution );
}

void ReservoirSolverBase::applySystemSolution( DofManager const & dofManager,
                                               arrayView1d< real64 const > const & localSolution,
                                               real64 const scalingFactor,
                                               DomainPartition & domain )
{
  // update the reservoir variables
  m_flowSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update the well variables
  m_wellSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
}

void ReservoirSolverBase::updateState( DomainPartition & domain )
{
  m_flowSolver->updateState( domain );
  m_wellSolver->updateState( domain );
}

void ReservoirSolverBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // reset reservoir variables
  m_flowSolver->resetStateToBeginningOfStep( domain );
  // reset well variables
  m_wellSolver->resetStateToBeginningOfStep( domain );
}

void ReservoirSolverBase::implicitStepComplete( real64 const & time_n,
                                                real64 const & dt,
                                                DomainPartition & domain )
{
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
  m_wellSolver->implicitStepComplete( time_n, dt, domain );
}

void ReservoirSolverBase::resetViews( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

real64 ReservoirSolverBase::scalingForSystemSolution( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localSolution )
{
  real64 const flowScalingFactor = m_flowSolver->scalingForSystemSolution( domain, dofManager, localSolution );
  real64 const wellScalingFactor = m_wellSolver->scalingForSystemSolution( domain, dofManager, localSolution );

  GEOSX_LOG_LEVEL_RANK_0( 2, "Scaling factor for the reservoir: " << flowScalingFactor
                                                                  << "; for the well(s): " << wellScalingFactor );

  return LvArray::math::min( flowScalingFactor, wellScalingFactor );
}


} /* namespace geosx */
