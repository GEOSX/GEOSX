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
 * @file CompositionalMultiphaseFVMKernels.cpp
 */

#include "CompositionalMultiphaseFVMKernels.hpp"
#include "CompositionalMultiphaseFlowUpwindHelperKernels.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/FaceElementStencil.hpp"

namespace geosx
{

namespace CompositionalMultiphaseFVMKernels
{

/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
PhaseMobilityKernel::
  compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( phaseDens ),
           arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM( dPhaseDens_dPres ),
           arraySlice2d< real64 const > const & GEOSX_UNUSED_PARAM( dPhaseDens_dComp ),
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp )
{
  real64 dRelPerm_dC[NC];
  real64 dVisc_dC[NC];

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const viscosity = phaseVisc[ip];
    real64 const dVisc_dP = dPhaseVisc_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseVisc_dComp[ip], dVisc_dC );

    real64 const relPerm = phaseRelPerm[ip];
    real64 dRelPerm_dP = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dRelPerm_dC[ic] = 0.0;
    }

    for( localIndex jp = 0; jp < NP; ++jp )
    {
      real64 const dRelPerm_dS = dPhaseRelPerm_dPhaseVolFrac[ip][jp];
      dRelPerm_dP += dRelPerm_dS * dPhaseVolFrac_dPres[jp];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dRelPerm_dC[jc] += dRelPerm_dS * dPhaseVolFrac_dComp[jp][jc];
      }
    }

    real64 const mobility = relPerm  / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP / viscosity
                          - mobility * dVisc_dP / viscosity;

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] / viscosity
                                - mobility * dVisc_dC[jc] / viscosity;
    }
  }
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        compute< NC, NP >( dCompFrac_dCompDens[a],
                           phaseDens[a][0],
                           dPhaseDens_dPres[a][0],
                           dPhaseDens_dComp[a][0],
                           phaseVisc[a][0],
                           dPhaseVisc_dPres[a][0],
                           dPhaseVisc_dComp[a][0],
                           phaseRelPerm[a][0],
                           dPhaseRelPerm_dPhaseVolFrac[a][0],
                           dPhaseVolFrac_dPres[a],
                           dPhaseVolFrac_dComp[a],
                           phaseMob[a],
                           dPhaseMob_dPres[a],
                           dPhaseMob_dComp[a] );
      } );
}

template< localIndex NC, localIndex NP >
void PhaseMobilityKernel::
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseDens,
          arrayView3d< real64 const > const & dPhaseDens_dPres,
          arrayView4d< real64 const > const & dPhaseDens_dComp,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
      {
        localIndex const a = targetSet[ i ];
        compute< NC, NP >( dCompFrac_dCompDens[a],
                           phaseDens[a][0],
                           dPhaseDens_dPres[a][0],
                           dPhaseDens_dComp[a][0],
                           phaseVisc[a][0],
                           dPhaseVisc_dPres[a][0],
                           dPhaseVisc_dComp[a][0],
                           phaseRelPerm[a][0],
                           dPhaseRelPerm_dPhaseVolFrac[a][0],
                           dPhaseVolFrac_dPres[a],
                           dPhaseVolFrac_dComp[a],
                           phaseMob[a],
                           dPhaseMob_dPres[a],
                           dPhaseMob_dComp[a] );
      } );
}

#define INST_PhaseMobilityKernel( NC, NP ) \
  template \
  void \
  PhaseMobilityKernel:: \
    launch< NC, NP >( localIndex const size, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp ); \
  template \
  void \
  PhaseMobilityKernel:: \
    launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseVisc, \
                      arrayView3d< real64 const > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const > const & phaseRelPerm, \
                      arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64 > const & phaseMob, \
                      arrayView2d< real64 > const & dPhaseMob_dPres, \
                      arrayView3d< real64 > const & dPhaseMob_dComp )

INST_PhaseMobilityKernel( 1, 1 );
INST_PhaseMobilityKernel( 2, 1 );
INST_PhaseMobilityKernel( 3, 1 );
INST_PhaseMobilityKernel( 4, 1 );
INST_PhaseMobilityKernel( 5, 1 );

INST_PhaseMobilityKernel( 1, 2 );
INST_PhaseMobilityKernel( 2, 2 );
INST_PhaseMobilityKernel( 3, 2 );
INST_PhaseMobilityKernel( 4, 2 );
INST_PhaseMobilityKernel( 5, 2 );

INST_PhaseMobilityKernel( 1, 3 );
INST_PhaseMobilityKernel( 2, 3 );
INST_PhaseMobilityKernel( 3, 3 );
INST_PhaseMobilityKernel( 4, 3 );
INST_PhaseMobilityKernel( 5, 3 );

#undef INST_PhaseMobilityKernel

/******************************** FluxKernel ********************************/

template< CompositionalMultiphaseFlowUpwindHelperKernels::term T >
using UpwindSchemeType = CompositionalMultiphaseFlowUpwindHelperKernels::HybridUpwind< T >;

template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, bool IS_UT_FORM >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
FluxKernel::
  compute( localIndex const numPhases,
           localIndex const stencilSize,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           arraySlice1d< real64 const > const stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const localFlux,
           arraySlice2d< real64 > const localFluxJacobian )
{
  localIndex constexpr NDOF = NC + 1;
  localIndex const NP = numPhases;

  real64 compFlux[NC]{};
  real64 dCompFlux_dP[MAX_STENCIL][NC]{};
  real64 dCompFlux_dC[MAX_STENCIL][NC][NC]{};

  real64 totFlux_unw{};
  real64 dTotFlux_dP[MAX_STENCIL]{};
  real64 dTotFlux_dC[MAX_STENCIL][NC]{};


  using UpwindHelpers = geosx::CompositionalMultiphaseFlowUpwindHelperKernels::UpwindHelpers;

  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( localIndex ip = 0; ip < NP; ++ip )
  {

    real64 phaseFlux{};
    real64 dPhaseFlux_dP[MAX_STENCIL]{};
    real64 dPhaseFlux_dC[MAX_STENCIL][NC]{};

    localIndex k_up = -1;

    UpwindHelpers::formPPUVelocity< NC, NUM_ELEMS, MAX_STENCIL >(
      NP,
      ip,
      stencilSize,
      seri,
      sesri,
      sei,
      stencilWeights,
      pres,
      dPres,
      gravCoef,
      phaseMob,
      dPhaseMob_dPres,
      dPhaseMob_dComp,
      dPhaseVolFrac_dPres,
      dPhaseVolFrac_dComp,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens_dPres,
      dPhaseMassDens_dComp,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      capPressureFlag,
      k_up,
      phaseFlux,
      dPhaseFlux_dP,
      dPhaseFlux_dC
      );

    // updateing phase Flux
    totFlux_unw += phaseFlux;

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dTotFlux_dP[ke] += dPhaseFlux_dP[ke];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dTotFlux_dC[ke][jc] += dPhaseFlux_dC[ke][jc];
      }
    }

//    mdensmultiply
    UpwindHelpers::mdensMultiply(
      ip,
      k_up,
      stencilSize,
      seri,
      sesri,
      sei,
      dCompFrac_dCompDens,
      phaseDens,
      dPhaseDens_dPres,
      dPhaseDens_dComp,
      phaseFlux,
      dPhaseFlux_dP,
      dPhaseFlux_dC );

    if( !IS_UT_FORM ) // skip  if you intend to use fixed total velocity formulation
    {
      UpwindHelpers::formPhaseComp( ip,
                                    k_up,
                                    stencilSize,
                                    seri,
                                    sesri,
                                    sei,
                                    phaseCompFrac,
                                    dPhaseCompFrac_dPres,
                                    dPhaseCompFrac_dComp,
                                    dCompFrac_dCompDens,
                                    phaseFlux,
                                    dPhaseFlux_dP,
                                    dPhaseFlux_dC,
                                    compFlux,
                                    dCompFlux_dP,
                                    dCompFlux_dC );

    }
  }
  // *** end of upwinding

  //if total flux formulation
  if( IS_UT_FORM )
  {
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      // choose upstream cell
      // create local work arrays
      real64 phaseFlux{};
      real64 dPhaseFlux_dP[MAX_STENCIL]{};
      real64 dPhaseFlux_dC[MAX_STENCIL][NC]{};

      real64 phaseFluxV{};
      real64 dPhaseFluxV_dP[MAX_STENCIL]{};
      real64 dPhaseFluxV_dC[MAX_STENCIL][NC]{};

      real64 fflow{};
      real64 dFflow_dP[MAX_STENCIL]{};
      real64 dFflow_dC[MAX_STENCIL][NC]{};


      // and the fractional flow for viscous part as \lambda_i^{up}/\sum_{NP}(\lambda_j^{up}) with up decided upon
      // the Upwind strategy
      localIndex k_up = -1;
      UpwindHelpers::formFracFlow< NC, NUM_ELEMS, MAX_STENCIL,
                                   CompositionalMultiphaseFlowUpwindHelperKernels::term::Viscous,
                                   UpwindSchemeType >( NP,
                                                       ip,
                                                       stencilSize,
                                                       seri,
                                                       sesri,
                                                       sei,
                                                       stencilWeights,
                                                       totFlux_unw,
                                                       pres,
                                                       dPres,
                                                       gravCoef,
                                                       dCompFrac_dCompDens,
                                                       phaseMassDens,
                                                       dPhaseMassDens_dPres,
                                                       dPhaseMassDens_dComp,
                                                       phaseMob,
                                                       dPhaseMob_dPres,
                                                       dPhaseMob_dComp,
                                                       dPhaseVolFrac_dPres,
                                                       dPhaseVolFrac_dComp,
                                                       phaseCapPressure,
                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                       capPressureFlag,
                                                       k_up,
                                                       fflow,
                                                       dFflow_dP,
                                                       dFflow_dC );

      //mdensmultiply
      UpwindHelpers::mdensMultiply(
        ip,
        k_up,
        stencilSize,
        seri,
        sesri,
        sei,
        dCompFrac_dCompDens,
        phaseDens,
        dPhaseDens_dPres,
        dPhaseDens_dComp,
        fflow,
        dFflow_dP,
        dFflow_dC );

//      std::cout << k_up << " | ";
      // Assembling the viscous flux (and derivatives) from fractional flow and total velocity as \phi_{\mu} = f_i^{up,\mu} uT
      phaseFluxV = fflow * totFlux_unw;
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dPhaseFluxV_dP[ke] += dFflow_dP[ke] * totFlux_unw;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPhaseFluxV_dC[ke][jc] += dFflow_dC[ke][jc] * totFlux_unw;
        }
      }

      //NON-FIXED UT -- to be canceled out if considered fixed
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dPhaseFluxV_dP[ke] += fflow * dTotFlux_dP[ke];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPhaseFluxV_dC[ke][jc] += fflow * dTotFlux_dC[ke][jc];
        }
      }

      // accumulate in the flux and its derivatives
      phaseFlux += phaseFluxV;
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] += dPhaseFluxV_dP[ke];
        for( localIndex ic = 0; ic < NC; ++ic )
          dPhaseFlux_dC[ke][ic] += dPhaseFluxV_dC[ke][ic];
      }

      // Distributing the viscous flux of phase i onto component
      UpwindHelpers::formPhaseComp( ip,
                                    k_up,
                                    stencilSize,
                                    seri,
                                    sesri,
                                    sei,
                                    phaseCompFrac,
                                    dPhaseCompFrac_dPres,
                                    dPhaseCompFrac_dComp,
                                    dCompFrac_dCompDens,
                                    phaseFluxV,
                                    dPhaseFluxV_dP,
                                    dPhaseFluxV_dC,
                                    compFlux,
                                    dCompFlux_dP,
                                    dCompFlux_dC );

      /***           GRAVITY TERM                ***/
      localIndex k_up_g = -1;
      localIndex k_up_og = -1;

      real64 phaseFluxG{};
      real64 dPhaseFluxG_dP[MAX_STENCIL]{};
      real64 dPhaseFluxG_dC[MAX_STENCIL][NC]{};

      UpwindHelpers::formPotFluxes< NC,
                                    CompositionalMultiphaseFlowUpwindHelperKernels::term::Gravity,
                                    NUM_ELEMS, MAX_STENCIL, NUM_ELEMS, UpwindSchemeType >(
        NP,
        ip,
        stencilSize,
        seri,
        sesri,
        sei,
        stencilWeights,
        totFlux_unw,
        pres,
        dPres,
        gravCoef,
        phaseMob,
        dPhaseMob_dPres,
        dPhaseMob_dComp,
        dPhaseVolFrac_dPres,
        dPhaseVolFrac_dComp,
        dCompFrac_dCompDens,
        phaseDens,
        dPhaseDens_dPres,
        dPhaseDens_dComp,
        phaseMassDens,
        dPhaseMassDens_dPres,
        dPhaseMassDens_dComp,
        phaseCapPressure,
        dPhaseCapPressure_dPhaseVolFrac,
        capPressureFlag,
        k_up_g,
        k_up_og,
        phaseFlux,
        dPhaseFlux_dP,
        dPhaseFlux_dC
        );

      // Distributing the gravitational flux of phase i onto component
      UpwindHelpers::formPhaseComp( ip,
                                    k_up_g,
                                    stencilSize,
                                    seri,
                                    sesri,
                                    sei,
                                    phaseCompFrac,
                                    dPhaseCompFrac_dPres,
                                    dPhaseCompFrac_dComp,
                                    dCompFrac_dCompDens,
                                    phaseFluxG,
                                    dPhaseFluxG_dP,
                                    dPhaseFluxG_dC,
                                    compFlux,
                                    dCompFlux_dP,
                                    dCompFlux_dC );

      //update phaseFlux from gravitational
      phaseFlux += phaseFluxG;
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dPhaseFlux_dP[ke] += dPhaseFluxG_dP[ke];
        for( localIndex ic = 0; ic < NC; ++ic )
          dPhaseFlux_dC[ke][ic] += dPhaseFluxG_dC[ke][ic];
      }


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////Adding capillary
      /// contribution
      /***           CAPILLARY TERM                ***/
      if( capPressureFlag )
      {
        localIndex k_up_pc = -1;
        localIndex k_up_opc = -1;

        real64 phaseFluxCap{};
        real64 dPhaseFluxCap_dP[MAX_STENCIL]{};
        real64 dPhaseFluxCap_dC[MAX_STENCIL][NC]{};

        UpwindHelpers::formPotFluxes< NC,
                                      CompositionalMultiphaseFlowUpwindHelperKernels::term::Capillary,
                                      NUM_ELEMS, MAX_STENCIL, MAX_STENCIL, UpwindSchemeType >(
          NP,
          ip,
          stencilSize,
          seri,
          sesri,
          sei,
          stencilWeights,
          totFlux_unw,
          pres,
          dPres,
          gravCoef,
          phaseMob,
          dPhaseMob_dPres,
          dPhaseMob_dComp,
          dPhaseVolFrac_dPres,
          dPhaseVolFrac_dComp,
          dCompFrac_dCompDens,
          phaseDens,
          dPhaseDens_dPres,
          dPhaseDens_dComp,
          phaseMassDens,
          dPhaseMassDens_dPres,
          dPhaseMassDens_dComp,
          phaseCapPressure,
          dPhaseCapPressure_dPhaseVolFrac,
          capPressureFlag,
          k_up_pc,
          k_up_opc,
          phaseFluxCap,
          dPhaseFluxCap_dP,
          dPhaseFluxCap_dC
          );

        // Distributing the gravitational flux of phase i onto component
        UpwindHelpers::formPhaseComp( ip,
                                      k_up_pc,
                                      stencilSize,
                                      seri,
                                      sesri,
                                      sei,
                                      phaseCompFrac,
                                      dPhaseCompFrac_dPres,
                                      dPhaseCompFrac_dComp,
                                      dCompFrac_dCompDens,
                                      phaseFluxCap,
                                      dPhaseFluxCap_dP,
                                      dPhaseFluxCap_dC,
                                      compFlux,
                                      dCompFlux_dP,
                                      dCompFlux_dC );

        //update phaseFlux from gravitational and capillary part
        phaseFlux += phaseFluxCap;
        for( localIndex ke = 0; ke < stencilSize; ++ke )
        {
          dPhaseFlux_dP[ke] += dPhaseFluxCap_dP[ke];
          for( localIndex ic = 0; ic < NC; ++ic )
            dPhaseFlux_dC[ke][ic] += dPhaseFluxCap_dC[ke][ic];
        }

      }//end if capPressureFlag

//      std::cout << " \n ------ \n";

    }
  }//end if UT_FORM

  UpwindHelpers::fillLocalJacobi< NC, MAX_STENCIL, NDOF >( compFlux,
                                                           dCompFlux_dP,
                                                           dCompFlux_dC,
                                                           stencilSize,
                                                           dt,
                                                           localFlux,
                                                           localFluxJacobian );
}

template< localIndex NC, typename STENCIL_TYPE, bool IS_UT_FORM >
void
FluxKernel::
  launch( localIndex const numPhases,
          STENCIL_TYPE const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  typename STENCIL_TYPE::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  typename STENCIL_TYPE::WeightContainerViewConstType const & weights = stencil.getWeights();

  localIndex constexpr NUM_ELEMS   = STENCIL_TYPE::NUM_POINT_IN_FLUX;
  localIndex constexpr MAX_STENCIL = STENCIL_TYPE::MAX_STENCIL_SIZE;

  forAll< parallelDevicePolicy<> >( stencil.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
      {
        // TODO: hack! for MPFA, etc. must obtain proper size from e.g. seri
        localIndex const stencilSize = MAX_STENCIL;
        localIndex constexpr NDOF = NC + 1;

        stackArray1d< real64, NUM_ELEMS * NC >                      localFlux( NUM_ELEMS * NC );
        stackArray2d< real64, NUM_ELEMS * NC * MAX_STENCIL * NDOF > localFluxJacobian( NUM_ELEMS * NC, stencilSize * NDOF );

        FluxKernel::compute< NC, NUM_ELEMS, MAX_STENCIL, IS_UT_FORM >( numPhases,
                                                                       stencilSize,
                                                                       seri[iconn],
                                                                       sesri[iconn],
                                                                       sei[iconn],
                                                                       weights[iconn],
                                                                       pres,
                                                                       dPres,
                                                                       gravCoef,
                                                                       phaseMob,
                                                                       dPhaseMob_dPres,
                                                                       dPhaseMob_dComp,
                                                                       dPhaseVolFrac_dPres,
                                                                       dPhaseVolFrac_dComp,
                                                                       dCompFrac_dCompDens,
                                                                       phaseDens,
                                                                       dPhaseDens_dPres,
                                                                       dPhaseDens_dComp,
                                                                       phaseMassDens,
                                                                       dPhaseMassDens_dPres,
                                                                       dPhaseMassDens_dComp,
                                                                       phaseCompFrac,
                                                                       dPhaseCompFrac_dPres,
                                                                       dPhaseCompFrac_dComp,
                                                                       phaseCapPressure,
                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                       capPressureFlag,
                                                                       dt,
                                                                       localFlux,
                                                                       localFluxJacobian );

        // populate dof indices
        globalIndex dofColIndices[ MAX_STENCIL * NDOF ];
        for( localIndex i = 0; i < stencilSize; ++i )
        {
          globalIndex const offset = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

          for( localIndex jdof = 0; jdof < NDOF; ++jdof )
          {
            dofColIndices[i * NDOF + jdof] = offset + jdof;
          }
        }

        // TODO: apply equation/variable change transformation(s)

        // Add to residual/jacobian
        for( localIndex i = 0; i < NUM_ELEMS; ++i )
        {
          if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
          {
            globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
            localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
            GEOSX_ASSERT_GE( localRow, 0 );
            GEOSX_ASSERT_GT( localMatrix.numRows(), localRow + NC );

            for( localIndex ic = 0; ic < NC; ++ic )
            {
              RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow + ic], localFlux[i * NC + ic] );
              localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow + ic,
                                                                                dofColIndices,
                                                                                localFluxJacobian[i * NC + ic].dataIfContiguous(),
                                                                                stencilSize * NDOF );
            }
          }
        }
      } );
}

#define INST_FluxKernel( NC, STENCIL_TYPE, IS_UT_FORM ) \
  template \
  void FluxKernel:: \
    launch< NC, STENCIL_TYPE, IS_UT_FORM >( localIndex const numPhases, \
                                            STENCIL_TYPE const & stencil, \
                                            globalIndex const rankOffset, \
                                            ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber, \
                                            ElementViewConst< arrayView1d< integer const > > const & ghostRank, \
                                            ElementViewConst< arrayView1d< real64 const > > const & pres, \
                                            ElementViewConst< arrayView1d< real64 const > > const & dPres, \
                                            ElementViewConst< arrayView1d< real64 const > > const & gravCoef, \
                                            ElementViewConst< arrayView2d< real64 const > > const & phaseMob, \
                                            ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres, \
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp, \
                                            ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres, \
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp, \
                                            ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens, \
                                            ElementViewConst< arrayView3d< real64 const > > const & phaseDens, \
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres, \
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp, \
                                            ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens, \
                                            ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres, \
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp, \
                                            ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac, \
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres, \
                                            ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp, \
                                            ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure, \
                                            ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac, \
                                            integer const capPressureFlag, \
                                            real64 const dt, \
                                            CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                            arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1, CellElementStencilTPFA, false );
INST_FluxKernel( 2, CellElementStencilTPFA, false );
INST_FluxKernel( 3, CellElementStencilTPFA, false );
INST_FluxKernel( 4, CellElementStencilTPFA, false );
INST_FluxKernel( 5, CellElementStencilTPFA, false );

INST_FluxKernel( 1, CellElementStencilTPFA, true );
INST_FluxKernel( 2, CellElementStencilTPFA, true );
INST_FluxKernel( 3, CellElementStencilTPFA, true );
INST_FluxKernel( 4, CellElementStencilTPFA, true );
INST_FluxKernel( 5, CellElementStencilTPFA, true );

INST_FluxKernel( 1, FaceElementStencil, false );
INST_FluxKernel( 2, FaceElementStencil, false );
INST_FluxKernel( 3, FaceElementStencil, false );
INST_FluxKernel( 4, FaceElementStencil, false );
INST_FluxKernel( 5, FaceElementStencil, false );

INST_FluxKernel( 1, FaceElementStencil, true );
INST_FluxKernel( 2, FaceElementStencil, true );
INST_FluxKernel( 3, FaceElementStencil, true );
INST_FluxKernel( 4, FaceElementStencil, true );
INST_FluxKernel( 5, FaceElementStencil, true );

#undef INST_FluxKernel

} // namespace CompositionalMultiphaseFVMKernels

} // namespace geosx
