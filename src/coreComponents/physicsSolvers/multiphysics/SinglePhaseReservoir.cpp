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
 * @file SinglePhaseReservoir.cpp
 *
 */

#include "SinglePhaseReservoir.hpp"

#include "common/TimingMacros.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWell.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVM.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhaseReservoir::SinglePhaseReservoir( const string & name,
                                            Group * const parent ):
  ReservoirSolverBase( name, parent )
{}

SinglePhaseReservoir::~SinglePhaseReservoir()
{}

void SinglePhaseReservoir::setupSystem( DomainPartition & domain,
                                        DofManager & dofManager,
                                        CRSMatrix< real64, globalIndex > & localMatrix,
                                        array1d< real64 > & localRhs,
                                        array1d< real64 > & localSolution,
                                        bool const setSparsity )
{
  ReservoirSolverBase::setupSystem( domain,
                                    dofManager,
                                    localMatrix,
                                    localRhs,
                                    localSolution,
                                    setSparsity );

//  // we need to set the dR_dAper CRS matrix in SinglePhaseFVM to handle the presence of fractures
//  if( dynamicCast< SinglePhaseFVM< SinglePhaseBase > * >( m_flowSolver ) )
//  {
//    SinglePhaseFVM< SinglePhaseBase > * fvmSolver = dynamicCast< SinglePhaseFVM< SinglePhaseBase > * >( m_flowSolver );
//    fvmSolver->setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
//  }
}


void SinglePhaseReservoir::addCouplingSparsityPattern( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       SparsityPatternView< globalIndex > const & pattern ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  // TODO: remove this and just call SolverBase::setupSystem when DofManager can handle the coupling

  // Populate off-diagonal sparsity between well and reservoir

  string const resDofKey  = dofManager.getKey( m_wellSolver->resElementDofName() );
  string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );

  localIndex const wellNDOF = m_wellSolver->numDofPerWellElement();

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const & resDofNumber =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    PerforationData const * const perforationData = subRegion.getPerforationData();

    // get the well degrees of freedom and ghosting info
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get the well element indices corresponding to each perforation
    arrayView1d< localIndex const > const & perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

    // Insert the entries corresponding to reservoir-well perforations
    // This will fill J_WR, and J_RW
    forAll< serialPolicy >( perforationData->size(), [=] ( localIndex const iperf )
    {
      // working arrays
      globalIndex eqnRowIndexRes = 0;
      stackArray1d< globalIndex, 2 > eqnRowIndicesWell( wellNDOF );
      globalIndex dofColIndexRes = 0;
      stackArray1d< globalIndex, 2 > dofColIndicesWell( wellNDOF );

      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];
      localIndex const iwelem = perfWellElemIndex[iperf];

      eqnRowIndexRes = resDofNumber[er][esr][ei] - rankOffset;
      dofColIndexRes = resDofNumber[er][esr][ei];

      for( localIndex idof = 0; idof < wellNDOF; ++idof )
      {
        eqnRowIndicesWell[idof] = wellElemDofNumber[iwelem] + idof - rankOffset;
        dofColIndicesWell[idof] = wellElemDofNumber[iwelem] + idof;
      }

      if( eqnRowIndexRes >= 0 && eqnRowIndexRes < pattern.numRows() )
      {
        for( localIndex j = 0; j < dofColIndicesWell.size(); ++j )
        {
          pattern.insertNonZero( eqnRowIndexRes, dofColIndicesWell[j] );
        }
      }

      for( localIndex i = 0; i < eqnRowIndicesWell.size(); ++i )
      {
        if( eqnRowIndicesWell[i] >= 0 && eqnRowIndicesWell[i] < pattern.numRows() )
        {
          pattern.insertNonZero( eqnRowIndicesWell[i], dofColIndexRes );
        }
      }
    } );
  } );
}

void SinglePhaseReservoir::assembleCouplingTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  string const resDofKey = dofManager.getKey( m_wellSolver->resElementDofName() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > const resDofNumberAccessor =
    elemManager.constructArrayViewAccessor< globalIndex, 1 >( resDofKey );
  ElementRegionManager::ElementViewConst< arrayView1d< globalIndex const > > const resDofNumber =
    resDofNumberAccessor.toNestedViewConst();
  globalIndex const rankOffset = dofManager.rankOffset();

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    PerforationData const * const perforationData = subRegion.getPerforationData();

    // get the degrees of freedom
    string const wellDofKey = dofManager.getKey( m_wellSolver->wellElementDofName() );
    arrayView1d< globalIndex const > const wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );

    // get well variables on perforations
    arrayView1d< real64 const > const perfRate =
      perforationData->getReference< array1d< real64 > >( SinglePhaseWell::viewKeyStruct::perforationRateString() );
    arrayView2d< real64 const > const dPerfRate_dPres =
      perforationData->getReference< array2d< real64 > >( SinglePhaseWell::viewKeyStruct::dPerforationRate_dPresString() );

    arrayView1d< localIndex const > const perfWellElemIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const resElementRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const resElementSubRegion =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const resElementIndex =
      perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

    // loop over the perforations and add the rates to the residual and jacobian
    forAll< parallelDevicePolicy<> >( perforationData->size(), [=] GEOSX_HOST_DEVICE ( localIndex const iperf )
    {
      // local working variables and arrays
      localIndex eqnRowIndices[ 2 ] = { -1 };
      globalIndex dofColIndices[ 2 ] = { -1 };

      real64 localPerf[ 2 ] = { 0.0 };
      real64 localPerfJacobian[ 2 ][ 2 ] = {{ 0.0 }};

      // get the reservoir (sub)region and element indices
      localIndex const er = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei = resElementIndex[iperf];

      // get the well element index for this perforation
      localIndex const iwelem = perfWellElemIndex[iperf];
      globalIndex const elemOffset = wellElemDofNumber[iwelem];

      // row index on reservoir side
      eqnRowIndices[WellSolverBase::SubRegionTag::RES] = resDofNumber[er][esr][ei] - rankOffset;
      // column index on reservoir side
      dofColIndices[WellSolverBase::SubRegionTag::RES] = resDofNumber[er][esr][ei];

      // row index on well side
      eqnRowIndices[WellSolverBase::SubRegionTag::WELL] = LvArray::integerConversion< localIndex >( elemOffset - rankOffset )
                                                          + SinglePhaseWell::RowOffset::MASSBAL;
      // column index on well side
      dofColIndices[WellSolverBase::SubRegionTag::WELL] = elemOffset
                                                          + SinglePhaseWell::ColOffset::DPRES;

      // populate local flux vector and derivatives
      localPerf[WellSolverBase::SubRegionTag::RES] = dt * perfRate[iperf];
      localPerf[WellSolverBase::SubRegionTag::WELL] = -localPerf[WellSolverBase::SubRegionTag::RES];

      for( localIndex ke = 0; ke < 2; ++ke )
      {
        localPerfJacobian[WellSolverBase::SubRegionTag::RES][ke] = dt * dPerfRate_dPres[iperf][ke];
        localPerfJacobian[WellSolverBase::SubRegionTag::WELL][ke] = -localPerfJacobian[WellSolverBase::SubRegionTag::RES][ke];
      }

      for( localIndex i = 0; i < 2; ++i )
      {
        if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnRowIndices[i],
                                                                            &dofColIndices[0],
                                                                            &localPerfJacobian[0][0] + 2 * i,
                                                                            2 );
          atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnRowIndices[i]], localPerf[i] );
        }
      }
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseReservoir, string const &, Group * const )

} /* namespace geosx */
