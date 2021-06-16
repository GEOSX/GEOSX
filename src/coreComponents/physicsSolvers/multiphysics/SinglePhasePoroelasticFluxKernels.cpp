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
 * @file SinglePhasePoroelasticFluxKernels.cpp
 */

#include "SinglePhasePoroelasticFluxKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

namespace SinglePhasePoroelasticFluxKernels
{

using namespace FluxKernelsHelper;

/******************************** EmbeddedSurfaceFluxKernel ********************************/

template<>
void EmbeddedSurfaceFluxKernel::
  launch< CellElementStencilTPFAWrapper >( CellElementStencilTPFAWrapper const & stencilWrapper,
                                           real64 const dt,
                                           globalIndex const rankOffset,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                           ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                                           ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                           ElementViewConst< arrayView2d< real64 const > > const & dens,
                                           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & mob,
                                           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                           ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                           R1Tensor const & gravityVector,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( jumpDofNumber );
  GEOSX_UNUSED_VAR( dPerm_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             transTMultiplier,
                                             gravityVector,
                                             localMatrix,
                                             localRhs );
}

template<>
void EmbeddedSurfaceFluxKernel::
  launch< EmbeddedSurfaceToCellStencilWrapper >( EmbeddedSurfaceToCellStencilWrapper const & stencilWrapper,
                                                 real64 const dt,
                                                 globalIndex const rankOffset,
                                                 ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                                 ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                                 ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                                 ElementViewConst< arrayView1d< real64 const > > const & pres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                                 ElementViewConst< arrayView2d< real64 const > > const & dens,
                                                 ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                                 ElementViewConst< arrayView1d< real64 const > > const & mob,
                                                 ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                                 ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                                 ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                                 ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                                 ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                                 R1Tensor const & gravityVector,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( jumpDofNumber );
  GEOSX_UNUSED_VAR( dPerm_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             transTMultiplier,
                                             gravityVector,
                                             localMatrix,
                                             localRhs );
}


template<>
void EmbeddedSurfaceFluxKernel::
  launch< SurfaceElementStencilWrapper >( SurfaceElementStencilWrapper const & stencilWrapper,
                                          real64 const dt,
                                          globalIndex const rankOffset,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & jumpDofNumber,
                                          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                          ElementViewConst< arrayView1d< real64 const > > const & pres,
                                          ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                          ElementViewConst< arrayView2d< real64 const > > const & dens,
                                          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & mob,
                                          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM ( transTMultiplier ),
                                          R1Tensor const & GEOSX_UNUSED_PARAM ( gravityVector ),
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex maxNumFluxElems = SurfaceElementStencilWrapper::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = SurfaceElementStencilWrapper::MAX_STENCIL_SIZE;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();


  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize * 2;// pressure and normal jump

    // working arrays
    stackArray1d< globalIndex, maxNumFluxElems > dofColIndices( numDofs );
    stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
    stackArray2d< real64, maxNumFluxElems * maxStencilSize > localFluxJacobian( numFluxElems, numDofs );

    // compute transmissibility
    real64 transmissibility[2], dTrans_dPres[2], dTrans_dAper[2];
    stencilWrapper.computeTransmissibility( iconn,
                                            permeability,
                                            dPerm_dPres,
                                            dPerm_dAper,
                                            transmissibility,
                                            dTrans_dPres,
                                            dTrans_dAper );

    compute( stencilSize,
             seri[iconn],
             sesri[iconn],
             sei[iconn],
             transmissibility,
             dTrans_dPres,
             dTrans_dAper,
             pres,
             dPres,
             gravCoef,
             dens,
             dDens_dPres,
             mob,
             dMob_dPres,
             dt,
             localFlux,
             localFluxJacobian );

    // extract DOF numbers
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex localDofIndex = 2 * i;
      dofColIndices[ localDofIndex ]     = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
      dofColIndices[ localDofIndex + 1 ] = jumpDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
    }

    for( localIndex i = 0; i < numFluxElems; ++i )
    {

      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                          dofColIndices.data(),
                                                                          localFluxJacobian[i].dataIfContiguous(),
                                                                          numDofs );

      }
    }
  } );
}

void EmbeddedSurfaceFluxKernel::
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[2],
           real64 const (&dTrans_dPres)[2],
           real64 const (&dTrans_dAper)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
{
  GEOSX_UNUSED_VAR( numFluxElems );

  real64 fluxVal = 0.0;
  real64 dFlux_dTrans[2] = {0.0, 0.0};
  real64 dFlux_dP[2] = {0.0, 0.0};

  computeSinglePhaseFlux( seri, sesri, sei,
                          transmissibility,
                          dTrans_dPres,
                          pres,
                          dPres,
                          gravCoef,
                          dens,
                          dDens_dPres,
                          mob,
                          dMob_dPres,
                          fluxVal,
                          dFlux_dP,
                          dFlux_dTrans );

  real64 dFlux_dAper[2] = {0.0, 0.0};
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    dFlux_dAper[ke] = dFlux_dTrans[ke] * dTrans_dAper[ke];
  }

  // populate local flux vector and derivatives
  flux[0] =  dt * fluxVal;
  flux[1] = -dt * fluxVal;

  for( localIndex ke = 0; ke < 2; ++ke )
  {
    localIndex const dofIndex = 2*ke;

    fluxJacobian[0][dofIndex]   =  dt * dFlux_dP[ke];
    fluxJacobian[0][dofIndex+1] =  dt * dFlux_dAper[ke];
    fluxJacobian[1][dofIndex]   = -dt * dFlux_dP[ke];
    fluxJacobian[1][dofIndex+1] = -dt * dFlux_dAper[ke];
  }
}

/******************************** FaceElementFluxKernel ********************************/

template<>
void FaceElementFluxKernel::
  launch< CellElementStencilTPFAWrapper >( CellElementStencilTPFAWrapper const & stencilWrapper,
                                           real64 const dt,
                                           globalIndex const rankOffset,
                                           ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                           ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                           ElementViewConst< arrayView1d< real64 const > > const & pres,
                                           ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                           ElementViewConst< arrayView2d< real64 const > > const & dens,
                                           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                           ElementViewConst< arrayView1d< real64 const > > const & mob,
                                           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                           ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                           ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                           R1Tensor const & gravityVector,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs,
                                           CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_UNUSED_VAR( dPerm_dAper );
  GEOSX_UNUSED_VAR( dR_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             transTMultiplier,
                                             gravityVector,
                                             localMatrix,
                                             localRhs );
}

template<>
void FaceElementFluxKernel::
  launch< FaceElementToCellStencilWrapper >( FaceElementToCellStencilWrapper const & stencilWrapper,
                                             real64 const dt,
                                             globalIndex const rankOffset,
                                             ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                             ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                             ElementViewConst< arrayView1d< real64 const > > const & pres,
                                             ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                             ElementViewConst< arrayView2d< real64 const > > const & dens,
                                             ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                             ElementViewConst< arrayView1d< real64 const > > const & mob,
                                             ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                             ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                             ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                             ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                             ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                             R1Tensor const & gravityVector,
                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                             arrayView1d< real64 > const & localRhs,
                                             CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_UNUSED_VAR( dPerm_dAper );
  GEOSX_UNUSED_VAR( dR_dAper );

  SinglePhaseFVMKernels::FluxKernel::launch( stencilWrapper,
                                             dt,
                                             rankOffset,
                                             pressureDofNumber,
                                             ghostRank,
                                             pres,
                                             dPres,
                                             gravCoef,
                                             dens,
                                             dDens_dPres,
                                             mob,
                                             dMob_dPres,
                                             permeability,
                                             dPerm_dPres,
                                             transTMultiplier,
                                             gravityVector,
                                             localMatrix,
                                             localRhs );
}

template<>
void FaceElementFluxKernel::
  launch< SurfaceElementStencilWrapper >( SurfaceElementStencilWrapper const & stencilWrapper,
                                          real64 const dt,
                                          globalIndex const rankOffset,
                                          ElementViewConst< arrayView1d< globalIndex const > > const & pressureDofNumber,
                                          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                          ElementViewConst< arrayView1d< real64 const > > const & pres,
                                          ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                          ElementViewConst< arrayView2d< real64 const > > const & dens,
                                          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                          ElementViewConst< arrayView1d< real64 const > > const & mob,
                                          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & permeability,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
                                          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dAper,
                                          ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                          R1Tensor const & gravityVector,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs,
                                          CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  constexpr localIndex maxNumFluxElems = SurfaceElementStencilWrapper::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize  = SurfaceElementStencilWrapper::MAX_STENCIL_SIZE;

  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
  typename SurfaceElementStencilWrapper::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencilWrapper.getCellCenterToEdgeCenters();

  static constexpr real64 TINY = 1e-10;

  forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
    localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );
    localIndex const numDofs = stencilSize; // pressures

    // working arrays
    stackArray1d< globalIndex, maxNumFluxElems > dofColIndices( numDofs );
    stackArray1d< localIndex, maxNumFluxElems > localColIndices( numFluxElems );

    stackArray1d< real64, maxNumFluxElems > localFlux( numFluxElems );
    stackArray2d< real64, maxNumFluxElems * maxStencilSize > localFluxJacobian( numFluxElems, numDofs );

    // need to store this for later use in determining the dFlux_dU terms when using better permeabilty approximations.
    stackArray2d< real64, maxNumFluxElems *maxStencilSize > dFlux_dAper( numFluxElems, stencilSize );

    // compute transmissibility
    real64 transmissiblity[2], dTrans_dPres[2], dTrans_dAper[2];
    stencilWrapper.computeTransmissibility( iconn,
                                            permeability,
                                            dPerm_dPres,
                                            dPerm_dAper,
                                            transmissiblity,
                                            dTrans_dPres,
                                            dTrans_dAper );

    // TODO does this need to be here??
    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];
    for( localIndex k = 0; k < numFluxElems; ++k )
    {

      localIndex const ei = sei[iconn][k];

      if( fabs( LvArray::tensorOps::AiBi< 3 >( cellCenterToEdgeCenters[iconn][k], gravityVector ) ) > TINY )
      {
        transmissiblity[k] *= transTMultiplier[er][esr][ei][1];
      }
      else
      {
        transmissiblity[k] *= transTMultiplier[er][esr][ei][0];
      }

    }

    compute( stencilSize,
             seri[iconn],
             sesri[iconn],
             sei[iconn],
             transmissiblity,
             dTrans_dPres,
             dTrans_dAper,
             pres,
             dPres,
             gravCoef,
             dens,
             dDens_dPres,
             mob,
             dMob_dPres,
             dt,
             localFlux,
             localFluxJacobian,
             dFlux_dAper );

    // extract DOF numbers
    for( localIndex i = 0; i < numDofs; ++i )
    {
      dofColIndices[i] = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
      localColIndices[i] = sei( iconn, i );
    }

    for( localIndex i = 0; i < numFluxElems; ++i )
    {
      if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
      {
        globalIndex const globalRow = pressureDofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
        localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                          dofColIndices.data(),
                                                                          localFluxJacobian[i].dataIfContiguous(),
                                                                          stencilSize );

        dR_dAper.addToRowBinarySearch< parallelDeviceAtomic >( sei( iconn, i ),
                                                               localColIndices.data(),
                                                               dFlux_dAper[i].dataIfContiguous(),
                                                               stencilSize );
      }
    }
  } );

}



GEOSX_HOST_DEVICE
void
FaceElementFluxKernel::compute( localIndex const numFluxElems,
                                arraySlice1d< localIndex const > const & seri,
                                arraySlice1d< localIndex const > const & sesri,
                                arraySlice1d< localIndex const > const & sei,
                                real64 const (&transmissibility)[2],
                                real64 const (&dTrans_dPres)[2],
                                real64 const (&dTrans_dAper)[2],
                                ElementViewConst< arrayView1d< real64 const > > const & pres,
                                ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                                ElementViewConst< arrayView2d< real64 const > > const & dens,
                                ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                ElementViewConst< arrayView1d< real64 const > > const & mob,
                                ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                                real64 const dt,
                                arraySlice1d< real64 > const & flux,
                                arraySlice2d< real64 > const & fluxJacobian,
                                arraySlice2d< real64 > const & dFlux_dAperture )
{

  GEOSX_UNUSED_VAR( numFluxElems );
  real64 fluxVal = 0.0;
  real64 dFlux_dTrans[2] = {0.0, 0.0};
  real64 dFlux_dP[2] = {0.0, 0.0};

  computeSinglePhaseFlux( seri, sesri, sei,
                          transmissibility,
                          dTrans_dPres,
                          pres,
                          dPres,
                          gravCoef,
                          dens,
                          dDens_dPres,
                          mob,
                          dMob_dPres,
                          fluxVal,
                          dFlux_dP,
                          dFlux_dTrans );

  real64 dFlux_dAper[2] = {0.0, 0.0};
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    dFlux_dAper[ke] = dFlux_dTrans[ke] * dTrans_dAper[ke];
  }

  // populate local flux vector and derivatives
  flux[0] =  dt * fluxVal;
  flux[1] = -dt * fluxVal;

  fluxJacobian[0][0] += dFlux_dP[0];
  fluxJacobian[0][1] += dFlux_dP[1];
  fluxJacobian[1][0] -= dFlux_dP[0];
  fluxJacobian[1][1] -= dFlux_dP[1];

  dFlux_dAperture[0][0] += dFlux_dAper[0];
  dFlux_dAperture[0][1] += dFlux_dAper[1];
  dFlux_dAperture[1][0] -= dFlux_dAper[0];
  dFlux_dAperture[1][1] -= dFlux_dAper[1];
}


}// namespace SinglePhaseFVMKernels

} // namespace geosx
