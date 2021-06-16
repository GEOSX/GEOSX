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
 * @file ProppantTransportKernels.cpp
 */

#include "ProppantTransportKernels.hpp"

#include "constitutive/fluid/ParticleFluidBase.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

namespace ProppantTransportKernels
{

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
AccumulationKernel::
  compute( localIndex const NC,
           real64 const proppantConcOld,
           real64 const proppantConcNew,
           arraySlice1d< real64 const > const & componentDensOld,
           arraySlice1d< real64 const > const & componentDensNew,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dCompDens_dPres ),
           arraySlice2d< real64 const > const & dCompDens_dCompConc,
           real64 const volume,
           real64 const packPoreVolume,
           real64 const proppantLiftVolume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian )
{

  // proppant mass conservation
  localAccum[0] = (proppantConcNew - proppantConcOld) * volume - proppantLiftVolume;

  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {
    for( localIndex c2 = 0; c2 < NC; ++c2 )
    {
      localAccumJacobian[c1][c2] = 0.0;
    }
  }

  localAccumJacobian[0][0] = volume;

  // component mass conservation
  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {

    localAccum[c1+1] = ( componentDensNew[c1] * (1.0 - proppantConcNew) - componentDensOld[c1] * (1.0 - proppantConcOld) ) * volume +
                       (componentDensNew[c1] - componentDensOld[c1]) * packPoreVolume;

    for( localIndex c2 = 0; c2 < NC; ++c2 )
    {
      localAccumJacobian[c1 + 1][c2 + 1] = dCompDens_dCompConc[c1][c2] * ( 1.0 - proppantConcNew ) * volume
                                           + dCompDens_dCompConc[c1][c2] * packPoreVolume;
    }

    localAccumJacobian[c1+1][0] = -componentDensNew[c1] * volume;
  }
}

void
AccumulationKernel::
  launch( localIndex const size,
          localIndex const NC,
          localIndex const NDOF,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & proppantConc,
          arrayView1d< real64 const > const & dProppantConc,
          arrayView2d< real64 const > const & componentDensOld,
          arrayView3d< real64 const > const & componentDens,
          arrayView3d< real64 const > const & dCompDens_dPres,
          arrayView4d< real64 const > const & dCompDens_dCompConc,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & proppantPackVf,
          arrayView1d< real64 const > const & proppantLiftFlux,
          real64 const dt,
          real64 const maxProppantConcentration,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] < 0 )
    {
      localIndex constexpr MAX_NC = ProppantTransport::MAX_NUM_COMPONENTS;
      stackArray1d< globalIndex, MAX_NC > localAccumDOF( NDOF );
      stackArray1d< real64, MAX_NC > localAccum( NDOF );
      stackArray2d< real64, MAX_NC * MAX_NC > localAccumJacobian( NDOF, NDOF );

      real64 effectiveVolume = volume[ei];
      real64 packPoreVolume = 0.0;

      if( proppantPackVf[ei] < 1.0 )
      {
        effectiveVolume = volume[ei] * ( 1.0 - proppantPackVf[ei] );
        packPoreVolume = volume[ei] * proppantPackVf[ei] * ( 1.0 - maxProppantConcentration );
      }

      real64 const proppantLiftVolume = proppantLiftFlux[ei] * dt;

      compute( NC,
               proppantConc[ei],
               proppantConc[ei] + dProppantConc[ei],
               componentDensOld[ei],
               componentDens[ei][0],
               dCompDens_dPres[ei][0],
               dCompDens_dCompConc[ei][0],
               effectiveVolume,
               packPoreVolume,
               proppantLiftVolume,
               localAccum,
               localAccumJacobian );

      globalIndex const elemDOF = dofNumber[ei];

      for( localIndex idof = 0; idof < NDOF; ++idof )
      {
        localAccumDOF[idof] = elemDOF + idof;
      }
      // add contribution to global residual and dRdP

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
GEOSX_FORCE_INLINE
void
FluxKernel::
  computeJunction( localIndex const numElems,
                   localIndex const numDofPerCell,
                   arraySlice1d< localIndex const > const & stencilElementIndices,
                   arraySlice1d< real64 const > const & stencilWeights,
                   arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                   arrayView1d< real64 const > const & pres,
                   arrayView1d< real64 const > const & dPres,
                   arrayView1d< real64 const > const & proppantConc,
                   arrayView1d< real64 const > const & dProppantConc,
                   arrayView3d< real64 const > const & componentDens,
                   arrayView3d< real64 const > const & dComponentDens_dPres,
                   arrayView4d< real64 const > const & dComponentDens_dComponentConc,
                   arrayView1d< real64 const > const & gravDepth,
                   arrayView2d< real64 const > const & dens,
                   arrayView2d< real64 const > const & dDens_dPres,
                   arrayView2d< real64 const > const & dDens_dProppantConc,
                   arrayView3d< real64 const > const & dDens_dComponentConc,
                   arrayView2d< real64 const > const & visc,
                   arrayView2d< real64 const > const & dVisc_dPres,
                   arrayView2d< real64 const > const & dVisc_dProppantConc,
                   arrayView3d< real64 const > const & dVisc_dComponentConc,
                   arrayView2d< real64 const > const & fluidDensity,
                   arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDens_dPres ),
                   arrayView3d< real64 const > const & GEOSX_UNUSED_PARAM( dFluidDens_dComponentConc ),
                   arrayView1d< real64 const > const & settlingFactor,
                   arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dPres ),
                   arrayView1d< real64 const > const & dSettlingFactor_dProppantConc,
                   arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                   arrayView1d< real64 const > const & collisionFactor,
                   arrayView1d< real64 const > const & dCollisionFactor_dProppantConc,
                   arrayView1d< integer const > const & isProppantMobile,
                   arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                   arrayView1d< real64 const > const & aperture,
                   R1Tensor const & unitGravityVector,
                   arrayView2d< real64 const > const & transTMultiplier,
                   real64 const dt,
                   arraySlice1d< real64 > const & localFlux,
                   arraySlice2d< real64 > const & localFluxJacobian )
{

  // We assume numElems == stencilSize;

  real64 const TINY = 1e-10;

  // working array
  constexpr localIndex maxNumFluxElems = SurfaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxNumComponents = ProppantTransport::MAX_NUM_COMPONENTS;

  localIndex const NC = numDofPerCell - 1;

  stackArray1d< real64, maxNumFluxElems > weight( numElems );

  // mixture density and fluid density in each face
  stackArray1d< real64, maxNumFluxElems > mixDens( numElems );
  stackArray1d< real64, maxNumFluxElems > fluidDens( numElems );

  // realted to slip velocity calculation
  stackArray1d< real64, maxNumFluxElems > transT( numElems );
  stackArray1d< real64, maxNumFluxElems > coefs( numElems );

  real64 edgeDensity = 0.0;
  stackArray1d< real64, maxNumFluxElems > dEdgeDens_dP( numElems );
  stackArray1d< real64, maxNumFluxElems > dEdgeDens_dProppantC( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dEdgeDens_dComponentC( numElems, NC );

  real64 edgeViscosity = 0.0;
  stackArray1d< real64, maxNumFluxElems > dEdgeVisc_dP( numElems );
  stackArray1d< real64, maxNumFluxElems > dEdgeVisc_dProppantC( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dEdgeVisc_dComponentC( numElems, NC );

  stackArray1d< real64, maxNumFluxElems > dPe_dP( numElems );
  stackArray1d< real64, maxNumFluxElems > dPe_dProppantC( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dPe_dComponentC( numElems, NC );

  stackArray1d< real64, maxNumFluxElems > edgeToFaceFlux( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFlux_dP( numElems, numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceFlux_dComponentC( numElems, numElems, NC );

  stackArray1d< real64, maxNumFluxElems > edgeToFaceProppantFlux( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceProppantFlux_dP( numElems, numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceProppantFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceProppantFlux_dComponentC( numElems, numElems, NC );

  stackArray1d< real64, maxNumFluxElems > edgeToFaceFluidFlux( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFluidFlux_dP( numElems, numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumFluxElems > dEdgeToFaceFluidFlux_dProppantC( numElems, numElems );
  stackArray3d< real64, maxNumFluxElems * maxNumFluxElems * maxNumComponents > dEdgeToFaceFluidFlux_dComponentC( numElems, numElems, NC );

  stackArray1d< real64, maxNumFluxElems > proppantC( numElems );

  real64 sumOfWeights = 0.0;

  for( localIndex i=0; i<numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const aperTerm = aperture[ei] * aperture[ei] * aperture[ei];

    sumOfWeights += stencilWeights[i];
    weight[i] = stencilWeights[i];

    transT[i] = aperTerm * stencilWeights[i];

    real64 const stencilCellToEdgeDistance = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[i] );
    real64 const edgeLength = 12.0 * stencilWeights[i] * stencilCellToEdgeDistance;

    real64 const stencilEdgeToFaceDownDistance =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[i], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance;

    coefs[i] = stencilEdgeToFaceDownDistance * aperture[ei];

    if( fabs( stencilEdgeToFaceDownDistance ) > TINY )
    {
      // vertical flow component
      transT[i] *= transTMultiplier[ei][1];
    }
    else
    {
      // horizontal flow component
      transT[i] *= transTMultiplier[ei][0];
    }
  }

  for( localIndex i=0; i<numElems; ++i )
  {
    weight[i] /= sumOfWeights;
  }

  localIndex numberOfMobileProppantElems = 0;

  //get averaged edgeDensity and edgeViscosity
  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    edgeDensity += weight[i] * dens[ei][0];
    dEdgeDens_dP[i] = weight[i] * dDens_dPres[ei][0];
    dEdgeDens_dProppantC[i] = weight[i] * dDens_dProppantConc[ei][0];

    edgeViscosity += weight[i] * visc[ei][0];
    dEdgeVisc_dP[i] = weight[i] * dVisc_dPres[ei][0];
    dEdgeVisc_dProppantC[i] = weight[i] * dVisc_dProppantConc[ei][0];

    proppantC[i] = proppantConc[ei] + dProppantConc[ei];

    mixDens[i] = dens[ei][0];
    fluidDens[i] = fluidDensity[ei][0];

    if( isProppantMobile[ei] == 1 )
    {
      numberOfMobileProppantElems++;
    }

    for( localIndex c = 0; c < NC; ++c )
    {
      dEdgeDens_dComponentC[i][c] = weight[i] * dDens_dComponentConc[ei][0][c];
      dEdgeVisc_dComponentC[i][c] = weight[i] * dVisc_dComponentConc[ei][0][c];
    }
  }

  real64 const proppantFluxCoef = ( numberOfMobileProppantElems > 1 ) ? 1.0 : 0.0;

  real64 transTSum = 0.0;
  real64 Pe = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    Pe += transT[i] * (pres[ei] + dPres[ei] - gravTerm);
    transTSum += transT[i];
    dPe_dP[i] += transT[i];

    for( localIndex j = 0; j < numElems; ++j )
    {
      dPe_dP[j] += -transT[i] * gravD * dEdgeDens_dP[j];
      dPe_dProppantC[j] += -transT[i] * gravD * dEdgeDens_dProppantC[j];

      for( localIndex c = 0; c < NC; ++c )
      {
        dPe_dComponentC[j][c] += -transT[i] * gravD * dEdgeDens_dComponentC[j][c];
      }
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    dPe_dP[i] /= transTSum;
    dPe_dProppantC[i] /= transTSum;

    for( localIndex c = 0; c < NC; ++c )
    {
      dPe_dComponentC[i][c] /= transTSum;
    }
  }

  Pe /= transTSum;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    real64 const fluxTerm = Pe - (pres[ei] + dPres[ei] - gravTerm);

    edgeToFaceFlux[i] = transT[i] * fluxTerm / edgeViscosity;
    dEdgeToFaceFlux_dP[i][i] += -transT[i] / edgeViscosity;

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFlux_dP[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dP[j] / (edgeViscosity * edgeViscosity) + transT[i] *
                                  (dPe_dP[j] + dEdgeDens_dP[j] * gravD) / edgeViscosity;

      dEdgeToFaceFlux_dProppantC[i][j] += -transT[i] * fluxTerm * dEdgeVisc_dProppantC[j] / (edgeViscosity * edgeViscosity) + transT[i] *
                                          (dPe_dProppantC[j] + dEdgeDens_dProppantC[j] * gravD) / edgeViscosity;

      for( localIndex c = 0; c < NC; ++c )
      {
        dEdgeToFaceFlux_dComponentC[i][j][c] += -transT[i] * fluxTerm * dEdgeVisc_dComponentC[j][c] / (edgeViscosity * edgeViscosity) + transT[i] *
                                                (dPe_dComponentC[j][c] + dEdgeDens_dComponentC[j][c] * gravD) / edgeViscosity;
      }
    }

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceProppantFlux_dProppantC[i][j] = 0.0;
      for( localIndex c = 0; c < NC; ++c )
      {
        dEdgeToFaceProppantFlux_dComponentC[i][j][c] = 0.0;
      }
    }

    if( fabs( coefs[i] ) > TINY )
    {
      // vertical
      edgeToFaceProppantFlux[i] = (1.0 - proppantC[i]) * settlingFactor[ei] * coefs[i] * fluidDens[i] / mixDens[i];

      dEdgeToFaceProppantFlux_dProppantC[i][i] = (-settlingFactor[ei] + (1 - proppantC[i]) * dSettlingFactor_dProppantConc[ei]) *
                                                 coefs[i] * fluidDens[i] / mixDens[i];

      edgeToFaceProppantFlux[i] += proppantFluxCoef * edgeToFaceFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dEdgeToFaceProppantFlux_dProppantC[i][j] += proppantFluxCoef * dEdgeToFaceFlux_dProppantC[i][j];

        for( localIndex c = 0; c < NC; ++c )
        {
          dEdgeToFaceProppantFlux_dComponentC[i][j][c] += proppantFluxCoef * dEdgeToFaceFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // horizontal
      edgeToFaceProppantFlux[i] = (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) * proppantFluxCoef * edgeToFaceFlux[i];

      dEdgeToFaceProppantFlux_dProppantC[i][i] = -fluidDens[i] / mixDens[i] * (collisionFactor[ei] - (1.0 - proppantC[i]) * dCollisionFactor_dProppantConc[ei]) *
                                                 proppantFluxCoef * edgeToFaceFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dEdgeToFaceProppantFlux_dProppantC[i][j] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) *
                                                    proppantFluxCoef * dEdgeToFaceFlux_dProppantC[i][j];

        for( localIndex c = 0; c < NC; ++c )
        {
          dEdgeToFaceProppantFlux_dComponentC[i][j][c] += (1.0 + fluidDens[i] / mixDens[i] * (1.0 - proppantC[i]) * collisionFactor[ei]) *
                                                          proppantFluxCoef * dEdgeToFaceFlux_dComponentC[i][j][c];
        }
      }
    }

    // fluid flux for component transport

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dP[i][j] = 0.0;
      dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;

      for( localIndex c = 0; c < NC; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;
      }
    }

    // note that all the fluid properties are from previous time step

    real64 const fluidFluxCoef = ( isProppantMobile[ei] == 0 || numElems == 1 ) ? 0.0 : 1.0;

    edgeToFaceFluidFlux[i] = mixDens[i] / fluidDens[i] * edgeToFaceFlux[i] - fluidFluxCoef * (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) /
                             fluidDens[i] * edgeToFaceProppantFlux[i];

    dEdgeToFaceFluidFlux_dProppantC[i][i] = -fluidFluxCoef * edgeToFaceProppantFlux[i];

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dProppantC[i][j] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dProppantC[i][j] - fluidFluxCoef *
                                               (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] * dEdgeToFaceProppantFlux_dProppantC[i][j];
      for( localIndex c = 0; c < NC; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] += mixDens[i] / fluidDens[i] * dEdgeToFaceFlux_dComponentC[i][j][c] - fluidFluxCoef *
                                                     (mixDens[i] -  fluidDens[i] * (1.0 - proppantC[i])) / fluidDens[i] *
                                                     dEdgeToFaceProppantFlux_dComponentC[i][j][c];
      }
    }

    for( localIndex j = 0; j < numElems; ++j )
    {
      dEdgeToFaceFluidFlux_dProppantC[i][j] = 0.0;
      for( localIndex c = 0; c < NC; ++c )
      {
        dEdgeToFaceFluidFlux_dComponentC[i][j][c] = 0.0;
      }
    }
  }

  // get proppantCe

  real64 proppantCe = 0.0;
  stackArray1d< real64, maxNumFluxElems > dProppantCe_dP( numElems );
  stackArray1d< real64, maxNumFluxElems > dProppantCe_dProppantC( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dProppantCe_dComponentC( numElems, NC );

  real64 downStreamFlux = 0.0;
  stackArray1d< real64, maxNumFluxElems > dDownStreamFlux_dP( numElems );
  stackArray1d< real64, maxNumFluxElems > dDownStreamFlux_dProppantC( numElems );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dDownStreamFlux_dComponentC( numElems, NC );

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    if( isProppantMobile[ei] == 0 )
      continue;

    if( edgeToFaceProppantFlux[i] >= 0.0 )
    {
      // downstream
      downStreamFlux += edgeToFaceProppantFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dDownStreamFlux_dP[j] += dEdgeToFaceProppantFlux_dP[i][j];
        dDownStreamFlux_dProppantC[j] += dEdgeToFaceProppantFlux_dProppantC[i][j];

        for( localIndex c = 0; c < NC; ++c )
        {
          dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceProppantFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // upstream
      proppantCe += -edgeToFaceProppantFlux[i] * proppantC[i];
      dProppantCe_dProppantC[i] += -edgeToFaceProppantFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dProppantCe_dP[j] += -dEdgeToFaceProppantFlux_dP[i][j] * proppantC[i];
        dProppantCe_dProppantC[j] += -dEdgeToFaceProppantFlux_dProppantC[i][j] * proppantC[i];

        for( localIndex c = 0; c < NC; ++c )
        {
          dProppantCe_dComponentC[j][c] += -dEdgeToFaceProppantFlux_dComponentC[i][j][c] * proppantC[i];
        }
      }
    }
  }

  if( downStreamFlux > 0.0 )
  {
    for( localIndex i = 0; i < numElems; ++i )
    {
      localIndex const ei = stencilElementIndices[i];

      if( isProppantMobile[ei] == 0 )
        continue;

      dProppantCe_dP[i] =  dProppantCe_dP[i] / downStreamFlux - proppantCe * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
      dProppantCe_dProppantC[i] =  dProppantCe_dProppantC[i] / downStreamFlux - proppantCe * dDownStreamFlux_dProppantC[i] /
                                  (downStreamFlux * downStreamFlux);;

      for( localIndex c = 0; c < NC; ++c )
      {
        dProppantCe_dComponentC[i][c] =  dProppantCe_dComponentC[i][c] / downStreamFlux - proppantCe * dDownStreamFlux_dComponentC[i][c] /
                                        (downStreamFlux * downStreamFlux);
      }
    }

    proppantCe = proppantCe / downStreamFlux;
  }
  else
  {
    proppantCe = 0.0;
    for( localIndex i = 0; i < numElems; ++i )
    {
      localIndex const ei = stencilElementIndices[i];

      if( isProppantMobile[ei] == 0 )
        continue;

      dProppantCe_dP[i] = 0.0;
      dProppantCe_dProppantC[i] =  weight[i];
      proppantCe += proppantC[i] * weight[i];
    }
  }

  // get componentCe

  stackArray1d< real64, maxNumComponents > componentCe( NC );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentCe_dP( numElems, NC );
  stackArray2d< real64, maxNumFluxElems * maxNumComponents > dComponentCe_dProppantC( numElems, NC );
  stackArray3d< real64, maxNumFluxElems * maxNumComponents * maxNumComponents > dComponentCe_dComponentC( numElems, NC, NC );

  downStreamFlux = 0.0;
  for( localIndex i = 0; i < numElems; ++i )
  {
    dDownStreamFlux_dP[i] = 0.0;
    dDownStreamFlux_dProppantC[i] = 0.0;
    for( localIndex c = 0; c < NC; ++c )
    {
      dDownStreamFlux_dComponentC[i][c] = 0.0;
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    if( edgeToFaceFluidFlux[i] >= 0.0 )
    {
      // downstream
      downStreamFlux += edgeToFaceFluidFlux[i];

      for( localIndex j = 0; j < numElems; ++j )
      {
        dDownStreamFlux_dP[j] += dEdgeToFaceFluidFlux_dP[i][j];
        dDownStreamFlux_dProppantC[j] += dEdgeToFaceFluidFlux_dProppantC[i][j];

        for( localIndex c = 0; c < NC; ++c )
        {
          dDownStreamFlux_dComponentC[j][c] += dEdgeToFaceFluidFlux_dComponentC[i][j][c];
        }
      }
    }
    else
    {
      // upstream
      for( localIndex c1 = 0; c1 < NC; ++c1 )
      {
        componentCe[c1] += -edgeToFaceFluidFlux[i] * componentDens[ei][0][c1];
        dComponentCe_dP[i][c1] += -edgeToFaceFluidFlux[i] * dComponentDens_dPres[ei][0][c1];

        for( localIndex c2 = 0; c2 < NC; ++c2 )
        {
          dComponentCe_dComponentC[i][c1][c2] += -edgeToFaceFluidFlux[i] * dComponentDens_dComponentConc[ei][0][c1][c2];
        }

        for( localIndex j = 0; j < numElems; ++j )
        {
          dComponentCe_dP[j][c1] += -dEdgeToFaceFluidFlux_dP[i][j] * componentDens[ei][0][c1];
          dComponentCe_dProppantC[j][c1] += -dEdgeToFaceFluidFlux_dProppantC[i][j] * componentDens[ei][0][c1];

          for( localIndex c2 = 0; c2 < NC; ++c2 )
          {
            dComponentCe_dComponentC[j][c1][c2] += -dEdgeToFaceFluidFlux_dComponentC[i][j][c2] * componentDens[ei][0][c1];
          }
        }
      }
    }
  }

  if( downStreamFlux > 0.0 )
  {
    for( localIndex c1 = 0; c1 < NC; ++c1 )
    {
      for( localIndex i = 0; i < numElems; ++i )
      {
        dComponentCe_dP[i][c1] =  dComponentCe_dP[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dP[i] / (downStreamFlux * downStreamFlux);
        dComponentCe_dProppantC[i][c1] =
          dComponentCe_dProppantC[i][c1] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dProppantC[i] / (downStreamFlux * downStreamFlux);

        for( localIndex c2 = 0; c2 < NC; ++c2 )
        {
          dComponentCe_dComponentC[i][c1][c2] =
            dComponentCe_dComponentC[i][c1][c2] / downStreamFlux - componentCe[c1] * dDownStreamFlux_dComponentC[i][c2] / (downStreamFlux * downStreamFlux);

        }
      }

      componentCe[c1] = componentCe[c1] / downStreamFlux;
    }
  }
  else
  {
    for( localIndex c = 0; c < NC; ++c )
    {
      componentCe[c] = 0.0;

      for( localIndex i = 0; i < numElems; ++i )
      {
        localIndex const ei = stencilElementIndices[i];

        componentCe[c] += componentDens[ei][0][c] * weight[i];
        dComponentCe_dP[i][c] = dComponentDens_dPres[ei][0][c] * weight[i];

        for( localIndex c2 = 0; c2 < NC; ++c2 )
        {
          dComponentCe_dComponentC[i][c][c2] = dComponentDens_dComponentConc[ei][0][c][c2] * weight[i];
        }
      }
    }
  }

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei = stencilElementIndices[i];

    localIndex idx1 = i * numDofPerCell; // proppant

    if( isProppantMobile[ei] == 1 && !(numElems == 1 && coefs[i] > TINY) )
    {
      if( edgeToFaceProppantFlux[i] >= 0.0 )
      {
        localFlux[idx1] = -proppantCe * edgeToFaceProppantFlux[i] * dt;
      }
      else
      {
        localFlux[idx1] = -proppantC[i] * edgeToFaceProppantFlux[i] * dt;
      }

      for( localIndex j = 0; j < numElems; ++j )
      {
        localIndex idx2 = j * numDofPerCell;

        if( edgeToFaceProppantFlux[i] >= 0.0 )
        {
          localFluxJacobian[idx1][idx2] =
            -(dProppantCe_dProppantC[j] * edgeToFaceProppantFlux[i] + proppantCe * dEdgeToFaceProppantFlux_dProppantC[i][j]) * dt;
        }
        else
        {
          localFluxJacobian[idx1][idx2] = -proppantC[i] * dEdgeToFaceProppantFlux_dProppantC[i][j] * dt;

          if( i == j )
          {
            localFluxJacobian[idx1][idx2] += -edgeToFaceProppantFlux[i] * dt;
          }
        }
      }
    }
    else
    {
      localFlux[idx1] = 0.0;

      for( localIndex j = 0; j < numElems; ++j )
      {
        for( localIndex c = 0; c < numDofPerCell; ++c )
        {
          localIndex idx2 = j * numDofPerCell + c;
          localFluxJacobian[idx1][idx2] = 0.0;
        }
      }
    }

    // component

    if( numElems > 1 )
    {
      for( localIndex c1 = 0; c1 < NC; ++c1 )
      {
        idx1 = i * numDofPerCell + 1 + c1;

        if( edgeToFaceFluidFlux[i] >= 0.0 )
        {
          localFlux[idx1] = -componentCe[c1] * edgeToFaceFluidFlux[i] * dt;
        }
        else
        {
          localFlux[idx1] = -componentDens[ei][0][c1] * edgeToFaceFluidFlux[i] * dt;
        }

        for( localIndex j = 0; j < numElems; ++j )
        {
          localIndex idx2 = j * numDofPerCell;

          if( edgeToFaceFluidFlux[i] >= 0.0 )
          {
            for( localIndex c2 = 0; c2 < NC; ++c2 )
            {
              localFluxJacobian[idx1][idx2 + 1 +
                                      c2] =
                -(dComponentCe_dComponentC[j][c1][c2] * edgeToFaceFluidFlux[i] + 0*componentCe[c1] * dEdgeToFaceFluidFlux_dComponentC[i][j][c2]) * dt;
            }
          }
          else
          {
            if( i == j )
            {
              for( localIndex c2 = 0; c2 < NC; ++c2 )
              {
                localFluxJacobian[idx1][idx2 + 1 + c2] += -dComponentDens_dComponentConc[ei][0][c1][c2] * edgeToFaceFluidFlux[i] * dt;
              }
            }
          }
        }
      }
    }
  }
}

template<>
void FluxKernel::
  launch< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                    localIndex const GEOSX_UNUSED_PARAM( numDofPerCell ),
                                    real64 const GEOSX_UNUSED_PARAM( dt ),
                                    globalIndex const GEOSX_UNUSED_PARAM( rankOffset ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                    integer const GEOSX_UNUSED_PARAM( updateProppantPacking ),
                                    R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                    ElementViewConst< arrayView1d< globalIndex const > > const & GEOSX_UNUSED_PARAM( dofNumber ),
                                    ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( ghostRank ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantConc ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dProppantConc ),
                                    ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( componentDens ),
                                    ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dPres ),
                                    ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dComponentConc ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dPres ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dProppantConc ),
                                    ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dComponentConc ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dPres ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dProppantConc ),
                                    ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dComponentConc ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dPres ),
                                    ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dComponentConc ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dPres ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dProppantConc ),
                                    ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( collisionFactor ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dCollisionFactor_dProppantConc ),
                                    ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                    ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                    CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                    arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR( "Not implemented" );
}

template<>
void FluxKernel::
  launch< EmbeddedSurfaceToCellStencil >( EmbeddedSurfaceToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                          localIndex const GEOSX_UNUSED_PARAM( numDofPerCell ),
                                          real64 const GEOSX_UNUSED_PARAM( dt ),
                                          globalIndex const GEOSX_UNUSED_PARAM( rankOffset ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                          integer const GEOSX_UNUSED_PARAM( updateProppantPacking ),
                                          R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                          ElementViewConst< arrayView1d< globalIndex const > > const & GEOSX_UNUSED_PARAM( dofNumber ),
                                          ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( ghostRank ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantConc ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dProppantConc ),
                                          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( componentDens ),
                                          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dPres ),
                                          ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dComponentConc ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dPres ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dProppantConc ),
                                          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dComponentConc ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dPres ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dProppantConc ),
                                          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dComponentConc ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dPres ),
                                          ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dComponentConc ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dPres ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dProppantConc ),
                                          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( collisionFactor ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dCollisionFactor_dProppantConc ),
                                          ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                          ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                          CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                          arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR( "Not implemented" );
}

template<>
void FluxKernel::
  launch< FaceElementToCellStencil >( FaceElementToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                      localIndex const GEOSX_UNUSED_PARAM( numDofPerCell ),
                                      real64 const GEOSX_UNUSED_PARAM( dt ),
                                      globalIndex const GEOSX_UNUSED_PARAM( rankOffset ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                      integer const GEOSX_UNUSED_PARAM( updateProppantPacking ),
                                      R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                      ElementViewConst< arrayView1d< globalIndex const > > const & GEOSX_UNUSED_PARAM( dofNumber ),
                                      ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( ghostRank ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantConc ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dProppantConc ),
                                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( componentDens ),
                                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dPres ),
                                      ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dComponentDens_dComponentConc ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dPres ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dProppantConc ),
                                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dDens_dComponentConc ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dPres ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dProppantConc ),
                                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dVisc_dComponentConc ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dPres ),
                                      ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dFluidDens_dComponentConc ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dPres ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dProppantConc ),
                                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dSettlingFactor_dComponentConc ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( collisionFactor ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dCollisionFactor_dProppantConc ),
                                      ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                      CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
                                      arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR( "Not implemented" );
}


template<>
void FluxKernel::
  launch< SurfaceElementStencil >( SurfaceElementStencil const & stencil,
                                   localIndex const numDofPerCell,
                                   real64 const dt,
                                   globalIndex const rankOffset,
                                   ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                   integer const updateProppantPacking,
                                   R1Tensor const & unitGravityVector,
                                   ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
                                   ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                                   ElementViewConst< arrayView1d< real64 const > > const & pres,
                                   ElementViewConst< arrayView1d< real64 const > > const & dPres,
                                   ElementViewConst< arrayView1d< real64 const > > const & proppantConc,
                                   ElementViewConst< arrayView1d< real64 const > > const & dProppantConc,
                                   ElementViewConst< arrayView3d< real64 const > > const & componentDens,
                                   ElementViewConst< arrayView3d< real64 const > > const & dComponentDens_dPres,
                                   ElementViewConst< arrayView4d< real64 const > > const & dComponentDens_dComponentConc,
                                   ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                   ElementViewConst< arrayView2d< real64 const > > const & dens,
                                   ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                                   ElementViewConst< arrayView2d< real64 const > > const & dDens_dProppantConc,
                                   ElementViewConst< arrayView3d< real64 const > > const & dDens_dComponentConc,
                                   ElementViewConst< arrayView2d< real64 const > > const & visc,
                                   ElementViewConst< arrayView2d< real64 const > > const & dVisc_dPres,
                                   ElementViewConst< arrayView2d< real64 const > > const & dVisc_dProppantConc,
                                   ElementViewConst< arrayView3d< real64 const > > const & dVisc_dComponentConc,
                                   ElementViewConst< arrayView2d< real64 const > > const & fluidDensity,
                                   ElementViewConst< arrayView2d< real64 const > > const & dFluidDens_dPres,
                                   ElementViewConst< arrayView3d< real64 const > > const & dFluidDens_dComponentConc,
                                   ElementViewConst< arrayView1d< real64 const > > const & settlingFactor,
                                   ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dPres,
                                   ElementViewConst< arrayView1d< real64 const > > const & dSettlingFactor_dProppantConc,
                                   ElementViewConst< arrayView2d< real64 const > > const & dSettlingFactor_dComponentConc,
                                   ElementViewConst< arrayView1d< real64 const > > const & collisionFactor,
                                   ElementViewConst< arrayView1d< real64 const > > const & dCollisionFactor_dProppantConc,
                                   ElementViewConst< arrayView1d< integer const > > const & isProppantMobile,
                                   ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
                                   ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                   arrayView1d< real64 > const & localRhs )
{
  constexpr localIndex maxNumFluxElems = SurfaceElementStencil::NUM_POINT_IN_FLUX;
  constexpr localIndex maxStencilSize = SurfaceElementStencil::MAX_STENCIL_SIZE;

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  constexpr localIndex DOF1 = maxNumFluxElems * constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;
  constexpr localIndex DOF2 = maxStencilSize * constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    if( ( numFluxElems > 1 || updateProppantPacking != 0 ) ) //isGhostConnectors[iconn][0] < 0 )
    {
      localIndex const stencilSize  = numFluxElems;
      localIndex const DOF = numFluxElems * numDofPerCell;

      // working arrays
      stackArray1d< globalIndex, DOF2 > dofColIndices( DOF );

      stackArray1d< real64, DOF1 > localFlux( DOF );
      stackArray2d< real64, DOF1 * DOF2 > localFluxJacobian( DOF, DOF );

      localIndex const er = seri[iconn][0];
      localIndex const esr = sesri[iconn][0];

      computeJunction( numFluxElems,
                       numDofPerCell,
                       sei[iconn],
                       weights[iconn],
                       cellCenterToEdgeCenters[iconn],
                       pres[er][esr],
                       dPres[er][esr],
                       proppantConc[er][esr],
                       dProppantConc[er][esr],
                       componentDens[er][esr],
                       dComponentDens_dPres[er][esr],
                       dComponentDens_dComponentConc[er][esr],
                       gravDepth[er][esr],
                       dens[er][esr],
                       dDens_dPres[er][esr],
                       dDens_dProppantConc[er][esr],
                       dDens_dComponentConc[er][esr],
                       visc[er][esr],
                       dVisc_dPres[er][esr],
                       dVisc_dProppantConc[er][esr],
                       dVisc_dComponentConc[er][esr],
                       fluidDensity[er][esr],
                       dFluidDens_dPres[er][esr],
                       dFluidDens_dComponentConc[er][esr],
                       settlingFactor[er][esr],
                       dSettlingFactor_dPres[er][esr],
                       dSettlingFactor_dProppantConc[er][esr],
                       dSettlingFactor_dComponentConc[er][esr],
                       collisionFactor[er][esr],
                       dCollisionFactor_dProppantConc[er][esr],
                       isProppantMobile[er][esr],
                       proppantPackVf[er][esr],
                       aperture[er][esr],
                       unitGravityVector,
                       transTMultiplier[er][esr],
                       dt,
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
    }
  } );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
FluxKernel::
  computeCellBasedFlux( localIndex const numElems,
                        arraySlice1d< localIndex const > const & stencilElementIndices,
                        arraySlice1d< real64 const > const & stencilWeights,
                        arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                        arrayView2d< real64 const > const & transMultiplier,
                        R1Tensor const & unitGravityVector,
                        arrayView1d< real64 const > const & pres,
                        arrayView1d< real64 const > const & gravDepth,
                        arrayView2d< real64 const > const & dens,
                        arrayView2d< real64 const > const & visc,
                        arrayView1d< real64 const > const & aperture,
                        arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                        arrayView2d< real64 > const & cellBasedFlux )
{

  real64 constexpr TINY = 1e-10;

  localIndex constexpr maxNumFluxElems = SurfaceElementStencil::NUM_POINT_IN_FLUX;

  stackArray1d< real64, maxNumFluxElems > weight( numElems );
  stackArray1d< real64, maxNumFluxElems > transT( numElems );

  // clear working arrays

  real64 aperTerm;
  real64 sumOfWeights = 0;

  for( localIndex i=0; i<numElems; ++i )
  {

    localIndex const ei  = stencilElementIndices[i];

    aperTerm = aperture[ei] * aperture[ei] * aperture[ei];

    sumOfWeights += stencilWeights[i];
    weight[i] = stencilWeights[i];

    transT[i] = aperTerm * stencilWeights[i];

    real64 const stencilCellToEdgeDistance = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[i] );
    real64 const edgeLength = 12.0 * stencilWeights[i] * stencilCellToEdgeDistance;
    real64 const stencilEdgeToFaceDownDistance =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[i], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance;

    if( fabs( stencilEdgeToFaceDownDistance ) > TINY )
    {
      // vertical flow component
      transT[i] *= transMultiplier[ei][1];
    }
    else
    {
      // horizontal flow component
      transT[i] *= transMultiplier[ei][0];
    }
  }

  for( localIndex i=0; i<numElems; ++i )
  {
    weight[i] /= sumOfWeights;
  }

  //get averaged edgeDensity and edgeViscosity
  real64 edgeDensity = 0.0;
  real64 edgeViscosity = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];
    edgeDensity += weight[i] * dens[ei][0];
    edgeViscosity += weight[i] * visc[ei][0];
  }

  real64 transTSum = 0.0;
  real64 Pe = 0.0;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    Pe += transT[i] * (pres[ei] - gravTerm);
    transTSum += transT[i];
  }

  Pe /= transTSum;

  for( localIndex i = 0; i < numElems; ++i )
  {
    localIndex const ei  = stencilElementIndices[i];

    real64 const gravD    = gravDepth[ei];
    real64 const gravTerm = edgeDensity * gravD;

    real64 const fluxTerm = Pe - (pres[ei] - gravTerm);
    real64 const edgeToFaceFlux = transT[i] * fluxTerm / edgeViscosity;

    for( localIndex k = 0; k < 3; ++k )
    {
      RAJA::atomicAdd( parallelDeviceAtomic{}, &cellBasedFlux[ei][k], -edgeToFaceFlux * stencilCellCenterToEdgeCenters[i][k] );
    }
  }
}

template<>
void FluxKernel::
  launchCellBasedFluxCalculation< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                            ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                                            R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                            ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                                            ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                                            ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                                            ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                                            ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                            ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                            ElementView< arrayView2d< real64 > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ) )
{}

template<>
void FluxKernel::
  launchCellBasedFluxCalculation< EmbeddedSurfaceToCellStencil >( EmbeddedSurfaceToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                                  ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                                                  R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                  ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                                                  ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                                                  ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                                                  ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                                                  ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                                  ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                                  ElementView< arrayView2d< real64 > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ) )
{}

template<>
void FluxKernel::
  launchCellBasedFluxCalculation< FaceElementToCellStencil >( FaceElementToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                              ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
                                                              R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                              ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                                                              ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravDepth ),
                                                              ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dens ),
                                                              ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( visc ),
                                                              ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                              ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                              ElementView< arrayView2d< real64 > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ) )
{}


template<>
void FluxKernel::
  launchCellBasedFluxCalculation< SurfaceElementStencil >( SurfaceElementStencil const & stencil,
                                                           FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & transTMultiplier,
                                                           R1Tensor const & unitGravityVector,
                                                           FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & pres,
                                                           FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & gravDepth,
                                                           FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & dens,
                                                           FluxKernel::ElementViewConst< arrayView2d< real64 const > > const & visc,
                                                           FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & aperture,
                                                           FluxKernel::ElementViewConst< arrayView1d< real64 const > > const & proppantPackVf,
                                                           FluxKernel::ElementView< arrayView2d< real64 > > const & cellBasedFlux )
{

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    computeCellBasedFlux( numFluxElems,
                          sei[iconn],
                          weights[iconn],
                          cellCenterToEdgeCenters[iconn],
                          transTMultiplier[er][esr],
                          unitGravityVector,
                          pres[er][esr],
                          gravDepth[er][esr],
                          dens[er][esr],
                          visc[er][esr],
                          aperture[er][esr],
                          proppantPackVf[er][esr],
                          cellBasedFlux[er][esr] );
  } );
}


template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeCalculation< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                                 real64 const GEOSX_UNUSED_PARAM( dt ),
                                                                 real64 const GEOSX_UNUSED_PARAM( proppantDensity ),
                                                                 real64 const GEOSX_UNUSED_PARAM( proppantDiameter ),
                                                                 real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                                 R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                 real64 const GEOSX_UNUSED_PARAM( criticalShieldsNumber ),
                                                                 real64 const GEOSX_UNUSED_PARAM( frictionCoefficient ),
                                                                 ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                                                 ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( density ),
                                                                 ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                                                 ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidViscosity ),
                                                                 ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                                 ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantBoundaryElement ),
                                                                 ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                                 ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( volume ),
                                                                 ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
                                                                 ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ),
                                                                 ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                                 ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                                 ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                                 ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantLiftFlux ) )
{}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeCalculation< EmbeddedSurfaceToCellStencil >( EmbeddedSurfaceToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                                       real64 const GEOSX_UNUSED_PARAM( dt ),
                                                                       real64 const GEOSX_UNUSED_PARAM( proppantDensity ),
                                                                       real64 const GEOSX_UNUSED_PARAM( proppantDiameter ),
                                                                       real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                                       R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                       real64 const GEOSX_UNUSED_PARAM( criticalShieldsNumber ),
                                                                       real64 const GEOSX_UNUSED_PARAM( frictionCoefficient ),
                                                                       ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                                                       ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( density ),
                                                                       ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                                                       ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidViscosity ),
                                                                       ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                                       ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantBoundaryElement ),
                                                                       ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                                       ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( volume ),
                                                                       ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
                                                                       ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ),
                                                                       ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                                       ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                                       ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                                       ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantLiftFlux ) )
{}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeCalculation< FaceElementToCellStencil >( FaceElementToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                                   real64 const GEOSX_UNUSED_PARAM( dt ),
                                                                   real64 const GEOSX_UNUSED_PARAM( proppantDensity ),
                                                                   real64 const GEOSX_UNUSED_PARAM( proppantDiameter ),
                                                                   real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                                   R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                   real64 const GEOSX_UNUSED_PARAM( criticalShieldsNumber ),
                                                                   real64 const GEOSX_UNUSED_PARAM( frictionCoefficient ),
                                                                   ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( settlingFactor ),
                                                                   ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( density ),
                                                                   ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidDensity ),
                                                                   ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( fluidViscosity ),
                                                                   ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                                   ElementViewConst< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantBoundaryElement ),
                                                                   ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( aperture ),
                                                                   ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( volume ),
                                                                   ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( elemGhostRank ),
                                                                   ElementView< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( cellBasedFlux ),
                                                                   ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                                   ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ),
                                                                   ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                                   ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantLiftFlux ) )
{}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ProppantPackVolumeKernel::
  computeProppantPackVolume( localIndex const numElems,
                             real64 const dt,
                             real64 const proppantDensity,
                             real64 const proppantDiameter,
                             real64 const maxProppantConcentration,
                             R1Tensor const & unitGravityVector,
                             real64 const criticalShieldsNumber,
                             real64 const frictionCoefficient,
                             arraySlice1d< localIndex const > const & stencilElementIndices,
                             arraySlice1d< real64 const > const & stencilWeights,
                             arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                             arrayView1d< real64 const > const & settlingFactor,
                             arrayView2d< real64 const > const & density,
                             arrayView2d< real64 const > const & fluidDensity,
                             arrayView2d< real64 const > const &,
                             arrayView1d< real64 const > const & volume,
                             arrayView1d< real64 const > const & aperture,
                             arrayView1d< integer const > const & elemGhostRank,
                             arrayView1d< integer const > const & isProppantBoundaryElement,
                             arrayView1d< integer const > const & isProppantMobile,
                             arrayView2d< real64 const > const & cellBasedFlux,
                             arrayView1d< real64 > const & conc,
                             arrayView1d< real64 > const & proppantPackVf,
                             arrayView1d< real64 > const & proppantExcessPackV,
                             arrayView1d< real64 > const & proppantLiftFlux )
{

  integer faceIndex = -1;

  real64 constexpr TINY = 1e-10;

  real64 const stencilCellToEdgeDistance0 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[0] );
  real64 const edgeLength = 12.0 * stencilWeights[0] * stencilCellToEdgeDistance0;

  if( numElems == 1 )
  {
    localIndex const ei  = stencilElementIndices[0];

    real64 const stencilEdgeToFaceDownDistance =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    if( stencilEdgeToFaceDownDistance < -TINY && isProppantMobile[ei] == 1 )
    {
      faceIndex = 0;
    }
  }
  else if( numElems == 2 )
  {
    localIndex const ei0  = stencilElementIndices[0];
    localIndex const ei1  = stencilElementIndices[1];

    real64 const stencilCellToEdgeDistance1 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[1] );

    real64 const stencilEdgeToFaceDownDistance0 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    real64 const stencilEdgeToFaceDownDistance1 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance1;

    //0: top  1: bottom

    if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 )
    {
      faceIndex = 0;
    }

    //0: bottom  1: top

    if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 )
    {
      faceIndex = 1;
    }
  }

  if( faceIndex >= 0 )
  {
    localIndex const ei = stencilElementIndices[faceIndex];

    if( elemGhostRank[ei] < 0 && isProppantBoundaryElement[ei] == 0 )
    {
      real64 const L = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[faceIndex] ) * 2.0;

      real64 velocity[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( cellBasedFlux[ei] );

      real64 const downVelocity = LvArray::tensorOps::AiBi< 3 >( cellBasedFlux[ei], unitGravityVector );

      LvArray::tensorOps::scaledAdd< 3 >( velocity, unitGravityVector, -downVelocity );

      real64 const velocityMag = LvArray::tensorOps::l2Norm< 3 >( velocity ) / volume[ei];

      real64 dH = fluidDensity[ei][0] / density[ei][0] * (1.0 - conc[ei]) * settlingFactor[ei] * conc[ei] / maxProppantConcentration * dt;

      real64 const tau = 1.0/8.0 * frictionCoefficient * fluidDensity[ei][0] * velocityMag * velocityMag;

      real64 const ShieldsNumber = tau / (proppantDensity - fluidDensity[ei][0]) / 9.81 / proppantDiameter;

      proppantLiftFlux[ei] = 0.0;

      if( ShieldsNumber > criticalShieldsNumber )
      {
        proppantLiftFlux[ei] = aperture[ei] *
                               (proppantDiameter * sqrt( 9.81 * proppantDiameter * (proppantDensity - fluidDensity[ei][0]) / fluidDensity[ei][0] )) *
                               (9.64 * pow( ShieldsNumber, 0.166 )) * pow( ShieldsNumber - criticalShieldsNumber, 1.5 );
      }

      if( proppantLiftFlux[ei] < 0.0 )
      {
        proppantLiftFlux[ei] = 0.0;
      }

      real64 liftH =  proppantLiftFlux[ei] / edgeLength / aperture[ei] / maxProppantConcentration * dt;

      real64 Vf = proppantPackVf[ei] + (dH - liftH) / L;

      if( Vf < 0.0 )
      {
        Vf = 0.0;
        liftH = dH + proppantPackVf[ei] * L;
        proppantLiftFlux[ei] = liftH * edgeLength * aperture[ei] * maxProppantConcentration / dt;
      }

      if( Vf >= 1.0 )
      {
        proppantExcessPackV[ei] = (Vf - 1.0) * L;
        proppantPackVf[ei] = 1.0;
        conc[ei] = maxProppantConcentration;
      }
      else
      {
        proppantPackVf[ei] = Vf;
      }
    }
  }
}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeCalculation< SurfaceElementStencil >( SurfaceElementStencil const & stencil,
                                                                real64 const dt,
                                                                real64 const proppantDensity,
                                                                real64 const proppantDiameter,
                                                                real64 const maxProppantConcentration,
                                                                R1Tensor const & unitGravityVector,
                                                                real64 const criticalShieldsNumber,
                                                                real64 const frictionCoefficient,
                                                                ElementView< arrayView1d< real64 const > > const & settlingFactor,
                                                                ElementView< arrayView2d< real64 const > > const & density,
                                                                ElementView< arrayView2d< real64 const > > const & fluidDensity,
                                                                ElementView< arrayView2d< real64 const > > const & fluidViscosity,
                                                                ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                                                ElementView< arrayView1d< integer const > > const & isProppantBoundaryElement,
                                                                ElementView< arrayView1d< real64 const > > const & aperture,
                                                                ElementView< arrayView1d< real64 const > > const & volume,
                                                                ElementView< arrayView1d< integer const > > const & elemGhostRank,
                                                                ElementView< arrayView2d< real64 const > > const & cellBasedFlux,
                                                                ElementView< arrayView1d< real64 > > const & conc,
                                                                ElementView< arrayView1d< real64 > > const & proppantPackVf,
                                                                ElementView< arrayView1d< real64 > > const & proppantExcessPackV,
                                                                ElementView< arrayView1d< real64 > > const & proppantLiftFlux )
{

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    computeProppantPackVolume( numFluxElems,
                               dt,
                               proppantDensity,
                               proppantDiameter,
                               maxProppantConcentration,
                               unitGravityVector,
                               criticalShieldsNumber,
                               frictionCoefficient,
                               sei[iconn],
                               weights[iconn],
                               cellCenterToEdgeCenters[iconn],
                               settlingFactor[er][esr],
                               density[er][esr],
                               fluidDensity[er][esr],
                               fluidViscosity[er][esr],
                               volume[er][esr],
                               aperture[er][esr],
                               elemGhostRank[er][esr],
                               isProppantBoundaryElement[er][esr],
                               isProppantMobile[er][esr],
                               cellBasedFlux[er][esr],
                               conc[er][esr],
                               proppantPackVf[er][esr],
                               proppantExcessPackV[er][esr],
                               proppantLiftFlux[er][esr] );
  } );
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ProppantPackVolumeKernel::
  updateProppantPackVolume( localIndex const numElems,
                            arraySlice1d< localIndex const > const & stencilElementIndices,
                            arraySlice1d< real64 const > const & stencilWeights,
                            arraySlice1d< R1Tensor const > const & stencilCellCenterToEdgeCenters,
                            R1Tensor const & unitGravityVector,
                            real64 const maxProppantConcentration,
                            arrayView1d< integer const > const & isProppantMobile,
                            arrayView1d< real64 const > const & proppantExcessPackV,
                            arrayView1d< real64 > const & conc,
                            arrayView1d< real64 > const & proppantPackVf )
{
  integer faceIndex = -1;
  real64 excessV = 0.0;

  real64 constexpr TINY = 1e-10;

  real64 const stencilCellToEdgeDistance0 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[0] );
  real64 const edgeLength = 12.0 * stencilWeights[0] * stencilCellToEdgeDistance0;

  if( numElems == 2 )
  {
    localIndex const ei0  = stencilElementIndices[0];
    localIndex const ei1  = stencilElementIndices[1];

    real64 const stencilCellToEdgeDistance1 = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[1] );

    real64 const stencilEdgeToFaceDownDistance0 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[0], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance0;

    real64 const stencilEdgeToFaceDownDistance1 =
      -LvArray::tensorOps::AiBi< 3 >( stencilCellCenterToEdgeCenters[1], unitGravityVector ) * edgeLength / stencilCellToEdgeDistance1;

    //0: top  1: bottom

    if( stencilEdgeToFaceDownDistance0 < -TINY && stencilEdgeToFaceDownDistance1 > TINY && isProppantMobile[ei0] == 1 && isProppantMobile[ei1] == 0 &&
        proppantExcessPackV[ei1] > 0.0 )
    {
      faceIndex = 0;
      excessV = proppantExcessPackV[ei1];
    }

    //0: bottom  1: top

    if( stencilEdgeToFaceDownDistance0 > TINY && stencilEdgeToFaceDownDistance1 < -TINY && isProppantMobile[ei1] == 1 && isProppantMobile[ei0] == 0 &&
        proppantExcessPackV[ei0] > 0.0 )
    {
      faceIndex = 1;
      excessV = proppantExcessPackV[ei0];
    }
  }

  if( faceIndex >= 0 )
  {
    localIndex const ei = stencilElementIndices[faceIndex];
    real64 const L = LvArray::tensorOps::l2Norm< 3 >( stencilCellCenterToEdgeCenters[faceIndex] ) * 2.0;
    real64 const Vf = proppantPackVf[ei] + excessV / L;

    if( Vf >= 1.0 )
    {
      proppantPackVf[ei] = 1.0;
      conc[ei] = maxProppantConcentration;
    }
    else
    {
      proppantPackVf[ei] = Vf;
    }
  }
}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeUpdate< CellElementStencilTPFA >( CellElementStencilTPFA const & GEOSX_UNUSED_PARAM( stencil ),
                                                            R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                            real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                            ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                            ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                            ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                            ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ) )
{}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeUpdate< EmbeddedSurfaceToCellStencil >( EmbeddedSurfaceToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                                  R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                                  real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                                  ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                                  ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                                  ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                                  ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ) )
{}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeUpdate< FaceElementToCellStencil >( FaceElementToCellStencil const & GEOSX_UNUSED_PARAM( stencil ),
                                                              R1Tensor const & GEOSX_UNUSED_PARAM( unitGravityVector ),
                                                              real64 const GEOSX_UNUSED_PARAM( maxProppantConcentration ),
                                                              ElementView< arrayView1d< integer const > > const & GEOSX_UNUSED_PARAM( isProppantMobile ),
                                                              ElementView< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( proppantExcessPackV ),
                                                              ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( conc ),
                                                              ElementView< arrayView1d< real64 > > const & GEOSX_UNUSED_PARAM( proppantPackVf ) )
{}

template<>
void ProppantPackVolumeKernel::
  launchProppantPackVolumeUpdate< SurfaceElementStencil >( SurfaceElementStencil const & stencil,
                                                           R1Tensor const & unitGravityVector,
                                                           real64 const maxProppantConcentration,
                                                           ElementView< arrayView1d< integer const > > const & isProppantMobile,
                                                           ElementView< arrayView1d< real64 const > > const & proppantExcessPackV,
                                                           ElementView< arrayView1d< real64 > > const & conc,
                                                           ElementView< arrayView1d< real64 > > const & proppantPackVf )
{

  typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename SurfaceElementStencil::WeightContainerViewConstType const & weights = stencil.getWeights();

  ArrayOfArraysView< R1Tensor const > const & cellCenterToEdgeCenters = stencil.getCellCenterToEdgeCenters();

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
  {
    localIndex const numFluxElems = seri.sizeOfArray( iconn );

    localIndex const er = seri[iconn][0];
    localIndex const esr = sesri[iconn][0];

    updateProppantPackVolume( numFluxElems,
                              sei[iconn],
                              weights[iconn],
                              cellCenterToEdgeCenters[iconn],
                              unitGravityVector,
                              maxProppantConcentration,
                              isProppantMobile[er][esr],
                              proppantExcessPackV[er][esr],
                              conc[er][esr],
                              proppantPackVf[er][esr] );
  } );
}

} // namespace ProppantTransportKernels

} // namespace geosx
