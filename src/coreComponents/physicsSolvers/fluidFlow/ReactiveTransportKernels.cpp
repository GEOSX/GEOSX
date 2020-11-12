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
 * @file ReactiveTransportKernels.cpp
 */

#include "ReactiveTransportKernels.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

namespace ReactiveTransportKernels
{

GEOSX_HOST_DEVICE
void
AccumulationKernel::
  Compute( localIndex const NC,
           arraySlice1d< real64 const > const &,
           arraySlice1d< real64 const > const & dComponentConc,
           real64 const volume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
{

  // component mass conservation
  for( localIndex c = 0; c < NC; ++c )
  {

    localAccum[c] = dComponentConc[c] * volume;
    localAccumJacobian[c][c] = volume;

  }
}

void
AccumulationKernel::
  Launch( localIndex const size,
          localIndex const NC,
          localIndex const NDOF,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView2d< real64 const > const & componentConc,
          arrayView2d< real64 const > const & dComponentConc,
          arrayView1d< real64 const > const & porosity,
          arrayView1d< real64 const > const & volume,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] < 0 )
    {
      localIndex constexpr MAX_NC = ReactiveTransport::MAX_NUM_COMPONENTS;
      stackArray1d< globalIndex, MAX_NC > localAccumDOF( NDOF );
      stackArray1d< real64, MAX_NC > localAccum( NDOF );
      stackArray2d< real64, MAX_NC * MAX_NC > localAccumJacobian( NDOF, NDOF );

      real64 effectiveVolume = volume[ei] * porosity[ei] /dt;

      Compute( NC,
               componentConc[ei],
               dComponentConc[ei],
               effectiveVolume,
               localAccum,
               localAccumJacobian );

      globalIndex const elemDOF = dofNumber[ei];

      for( localIndex idof = 0; idof < NDOF; ++idof )
      {
        localAccumDOF[idof] = elemDOF + idof;
      }

      localIndex const localRow = dofNumber[ei] - rankOffset;
      for( localIndex idof = 0; idof < NDOF; ++idof )
      {
        localRhs[localRow + idof] += localAccum[idof];
        localMatrix.addToRow< serialAtomic >( localRow + idof,
                                              localAccumDOF.data(),
                                              localAccumJacobian[idof].dataIfContiguous(),
                                              NDOF );
      }
    }
  } );
}

GEOSX_HOST_DEVICE
void
FluxKernel::
  Compute( localIndex const numElems,
           localIndex const NC,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView2d< real64 const > > const & componentConc,
           ElementViewConst< arrayView2d< real64 const > > const & dComponentConc,
           real64 const viscosity,
           arraySlice1d< real64 > const & localFlux,
           arraySlice2d< real64 > const & localFluxJacobian )
{

  localIndex constexpr maxNumFluxElems = CellElementStencilTPFA::MAX_STENCIL_SIZE;
  constexpr localIndex maxNumComponents = ReactiveTransport::MAX_NUM_COMPONENTS;

  stackArray2d< real64, maxNumFluxElems * maxNumComponents > concentration( numElems, NC );

  real64 potDif = 0.0;

  for( localIndex ke = 0; ke < numElems; ++ke )
  {
    localIndex const er  = seri[ke];
    localIndex const esr = sesri[ke];
    localIndex const ei  = sei[ke];

    real64 const weight = stencilWeights[ke];
    real64 const pressure = pres[er][esr][ei] + dPres[er][esr][ei];
    real64 const pot = weight * pressure;

    for( localIndex c = 0; c < NC; ++c )
      concentration[ke][c] = componentConc[er][esr][ei][c] + dComponentConc[er][esr][ei][c];

    potDif += pot;

  }

  localIndex id = potDif >= 0 ? 0 : 1;
  real64 const fluxVal = potDif / viscosity;
  /*
     real64 fluxVal;
     if(fabs(potDif) > 0.0)
     fluxVal = potDif / fabs(potDif);
     else
     fluxVal = 0.0*viscosity;
   */

  for( localIndex c = 0; c < NC; ++c )
  {

    localFlux[c] = fluxVal * concentration[id][c];
    localFlux[NC + c] = -fluxVal * concentration[id][c];

    localFluxJacobian[c][id * NC + c] = fluxVal;
    localFluxJacobian[NC + c][id * NC + c] = -fluxVal;

  }

}


template<>
void FluxKernel::
  Launch< CellElementStencilTPFA >( CellElementStencilTPFA const & stencil,
                                    localIndex const numDofPerCell,
                                    real64 const GEOSX_UNUSED_PARAM( dt ),
                                    globalIndex const rankOffset,
                                    ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
                                    ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                    ElementViewConst< arrayView1d< real64 const > > const & pres,
                                    ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                    ElementViewConst< arrayView2d< real64 const > > const & componentConc,
                                    ElementViewConst< arrayView2d< real64 const > > const & dComponentConc,
                                    real64 const viscosity,
                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                    arrayView1d< real64 > const & localRhs )
{

  constexpr localIndex maxNumFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex numFluxElems = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  constexpr localIndex stencilSize  = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  constexpr localIndex maxDOF = maxNumFluxElems * constitutive::ReactiveFluidBase::MAX_NUM_SPECIES;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename CellElementStencilTPFA::WeightContainerViewConstType const & weights = stencil.getWeights();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    // working arrays

    localIndex const DOF = numFluxElems * numDofPerCell;

    stackArray1d< globalIndex, maxDOF > dofColIndices( DOF );
    stackArray1d< real64, maxDOF > localFlux( DOF );
    stackArray2d< real64, maxDOF * maxDOF > localFluxJacobian( DOF, DOF );

    Compute( stencilSize,
             numDofPerCell,
             seri[iconn],
             sesri[iconn],
             sei[iconn],
             weights[iconn],
             pres,
             dPres,
             componentConc,
             dComponentConc,
             viscosity,
             localFlux,
             localFluxJacobian );

    for( localIndex i = 0; i < stencilSize; ++i )
    {
      for( localIndex j = 0; j < numDofPerCell; ++j )
      {
        dofColIndices[i * numDofPerCell + j] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] + j;
      }
    }

    for( localIndex i = 0; i < numFluxElems; ++i )
    {
      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( localMatrix.numRows(), localRow + numDofPerCell );

        for( localIndex idof = 0; idof < numDofPerCell; ++idof )
        {
          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + idof], localFlux[i * numDofPerCell + idof] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + idof,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i * numDofPerCell + idof].dataIfContiguous(),
                                                                            stencilSize * numDofPerCell );
        }
      }
    }
  } );

}

template<>
void FluxKernel::
  Launch< FaceElementStencil >( FaceElementStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                localIndex const GEOSX_UNUSED_PARAM( numDofPerCell ),
                                real64 const GEOSX_UNUSED_PARAM( dt ),
                                globalIndex const GEOSX_UNUSED_PARAM( rankOffset ),
                                ElementViewConst< arrayView1d< globalIndex const > > const & GEOSX_UNUSED_PARAM( dofNumber ),
                                ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( ghostRank ),
                                ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                                ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( componentConc ),
                                ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponenConc ),
                                real64 const GEOSX_UNUSED_PARAM( viscosity ),
                                CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{

  GEOSX_ERROR( "Not implemented" );

}

} // namespace ReactiveTransportKernels

} // namespace geosx
