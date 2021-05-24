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
 * @file CompositionalMultiphaseFlowUpwindHelperKernels.hpp
 */
#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HPP

namespace geosx
{

namespace CompositionalMultiphaseFlowUpwindHelperKernels
{
/**
 * @brief The type for element-based non-constitutive data parameters.
 * Consists entirely of ArrayView's.
 *
 * Can be converted from ElementRegionManager::ElementViewAccessor
 * by calling .toView() or .toViewConst() on an accessor instance
 */
template< typename VIEWTYPE >
using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;
/// This enum to select the proper physics in Upwind class specialization
enum class term
{
  Viscous, Gravity, Capillary
};

/*** struct holding actual static's ***/
struct UpwindHelpers
{
/******************************** generic templated  ********************************/

/**
 * @brief Form the PhasePotentialUpwind from pressure gradient and gravitational heads (legacy)
 * @tparam NC number of components
 * @tparam NUM_ELEMS numberof elements involve in the stencil connexion
 * @tparam MAX_STENCIL maximum number of points in the stencil
 * @param numPhase number of phases
 * @param ip index of the treated phase
 * @param stencilSize number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 * @param stencilWeights weights associated with elements in the stencil
 */
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void formPPUVelocity( localIndex const numPhase,
                               localIndex const ip,
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
                               ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                               ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                               ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                               ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                               ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                               integer const capPressureFlag,
                               localIndex( &k_up ),
                               real64( &phaseFlux ),
                               real64 ( & dPhaseFlux_dP )[MAX_STENCIL],
                               real64 ( & dPhaseFlux_dC )[MAX_STENCIL][NC]
                               )
  {

    real64 densMean{};
    real64 dDensMean_dP[NUM_ELEMS]{};
    real64 dDensMean_dC[NUM_ELEMS][NC]{};

    real64 presGrad{};
    real64 dPresGrad_dP[MAX_STENCIL]{};
    real64 dPresGrad_dC[MAX_STENCIL][NC]{};

    real64 gravHead{};
    real64 dGravHead_dP[NUM_ELEMS]{};
    real64 dGravHead_dC[NUM_ELEMS][NC]{};

    real64 dCapPressure_dC[NC]{};
// Working array
    real64 dProp_dC[NC]{};

// calculate quantities on primary connected cells
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // density
      real64 const density = phaseMassDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseMassDens_dPres[er][esr][ei][0][ip];

      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens_dComp[er][esr][ei][0][ip],
                      dProp_dC
                      );

// average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dP;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
      }
    }
// compute potential difference MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];
      real64 const weight = stencilWeights[i];

      // capillary pressure
      real64 capPressure = 0.0;
      real64 dCapPressure_dP = 0.0;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dCapPressure_dC[ic] = 0.0;
      }

      if( capPressureFlag )
      {
        capPressure = phaseCapPressure[er][esr][ei][0][ip];

        for( localIndex jp = 0; jp < numPhase; ++jp )
        {
          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP += dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCapPressure_dC[jc] += dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }
        }
      }

      presGrad += weight * ( pres[er][esr][ei] + dPres[er][esr][ei] - capPressure );
      dPresGrad_dP[i] += weight * ( 1 - dCapPressure_dP );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPresGrad_dC[i][jc] += -weight * dCapPressure_dC[jc];
      }

      real64 const gravD = weight * gravCoef[er][esr][ei];

      // the density used in the potential difference is always a mass density
      // unlike the density used in the phase mobility, which is a mass density
      // if useMass == 1 and a molar density otherwise
      gravHead += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < NUM_ELEMS; ++j )
      {
        dGravHead_dP[j] += dDensMean_dP[j] * gravD;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dGravHead_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
        }
      }
    }

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    k_up = ( potGrad >= 0 ) ? 0 : 1;

    localIndex er_up = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];
// pressure gradient depends on all points in the stencil
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dPhaseFlux_dP[ke] += dPresGrad_dP[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] += dPresGrad_dC[ke][jc];
      }

    }

    // gravitational head depends only on the two cells connected (same as mean density)
    for( localIndex ke = 0; ke < NUM_ELEMS; ++ke )
    {
      dPhaseFlux_dP[ke] -= dGravHead_dP[ke];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] -= dGravHead_dC[ke][jc];
      }
    }

    // compute the phase flux and derivatives using upstream cell mobility
    phaseFlux = mobility * potGrad;

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dPhaseFlux_dP[ke] *= mobility;

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPhaseFlux_dC[ke][jc] *= mobility;
      }
    }

    real64 const dMob_dP = dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
    arraySlice1d< real64 const > const dMob_dC = dPhaseMob_dComp[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dMob_dC[jc] * potGrad;
    }

  }

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void formWAVelocity( localIndex const numPhase,
                              localIndex const ip,
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
                              ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                              ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                              ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                              integer const capPressureFlag,
                              localIndex( &k_up ),
                              real64( &phaseFlux ),
                              real64 ( & dPhaseFlux_dP )[MAX_STENCIL],
                              real64 ( & dPhaseFlux_dC )[MAX_STENCIL][NC]
                              )
  {
    // to be implemented (Bosma et al.)



  }

/**
 * @brief Helper function to fill the Jacobi local entries from the compositional fluxes (and derivatives)
 * @tparam NC number of components
 * @tparam MAX_STENCIL maximum number of points in the stencil
 * @tparam NDOF number of degrees of freedom for this element
 * @param compFlux the component fluxes table
 * @param dCompFlux_dP the component fluxes derivatives wrt pressure
 * @param dCompFlux_dC the component fluxes derivatives wrt component
 * @param stencilSize number of points in the stencil
 * @param dt the time step increment
 * @param localFlux arraySlice of localFlux to be filled
 * @param localFluxJacobian arraySlice of local Jacobian to be filled
 */
  template< localIndex NC, localIndex MAX_STENCIL, localIndex NDOF >
  GEOSX_HOST_DEVICE
  static void fillLocalJacobi( real64 const (&compFlux)[NC],
                               real64 const (&dCompFlux_dP)[MAX_STENCIL][NC],
                               real64 const (&dCompFlux_dC)[MAX_STENCIL][NC][NC],
                               localIndex const stencilSize,
                               real64 const dt,
                               arraySlice1d< real64 > const localFlux,
                               arraySlice2d< real64 > const localFluxJacobian )
  {
    // populate jacobian from compnent fluxes (and derivatives)
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      localFlux[ic] = dt * compFlux[ic];
      localFlux[NC + ic] = -dt * compFlux[ic];

      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        localIndex const localDofIndexPres = ke * NDOF;
        localFluxJacobian[ic][localDofIndexPres] = dt * dCompFlux_dP[ke][ic];
        localFluxJacobian[NC + ic][localDofIndexPres] = -dt * dCompFlux_dP[ke][ic];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
          localFluxJacobian[ic][localDofIndexComp] = dt * dCompFlux_dC[ke][ic][jc];
          localFluxJacobian[NC + ic][localDofIndexComp] = -dt * dCompFlux_dC[ke][ic][jc];
        }
      }
    }
  }

  //Form potential-related parts of fluxes
  template< localIndex NC, term T, localIndex NUM_ELEMS, localIndex MAX_STENCIL, localIndex SIZE, template< term > class UPWIND >
  GEOSX_HOST_DEVICE
  static void formPotFluxes(
    localIndex const numPhase,
    localIndex const ip,
    localIndex const stencilSize,
    arraySlice1d< localIndex const > const seri,
    arraySlice1d< localIndex const > const sesri,
    arraySlice1d< localIndex const > const sei,
    arraySlice1d< real64 const > const stencilWeights,
    real64 const totFlux_unw,
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
    ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
    ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
    localIndex const capPressureFlag,
    localIndex (&k_up),
    localIndex (&k_up_o),
    real64 & phaseFlux,
    real64 (& dPhaseFlux_dP) [MAX_STENCIL],
    real64 (& dPhaseFlux_dC) [MAX_STENCIL][NC]
    )
  {

    real64 fflow{};
    real64 dFflow_dP[MAX_STENCIL]{};
    real64 dFflow_dC[MAX_STENCIL][NC]{};

    real64 pot{};
    real64 dPot_dP[SIZE]{};
    real64 dPot_dC[SIZE][NC]{};
    real64 dProp_dC[NC]{};

    //SIZE depends on T (if gravity then MAX_STENCIL, if Cap then NUM_ELEMS)
    UpwindHelpers::formPotential< NC, T, NUM_ELEMS, MAX_STENCIL >::compute(
      numPhase,
      ip,
      stencilSize,
      seri,
      sesri,
      sei,
      stencilWeights,
      totFlux_unw,
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens_dPres,
      dPhaseMassDens_dComp,
      dPhaseVolFrac_dPres,
      dPhaseVolFrac_dComp,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      pot,
      dPot_dP,
      dPot_dC,
      dProp_dC
      );

    // and the fractional flow for gravitational part as \lambda_i^{up}/\sum_{numPhase}(\lambda_k^{up}) with up decided upon
    // the Upwind strategy
    UpwindHelpers::formFracFlow< NC, NUM_ELEMS, MAX_STENCIL, T, UPWIND >( numPhase,
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
      dFflow_dC
      );

//      std::cout << k_up_g << " ";

    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( ip != jp )
      {

        real64 potOther{};
        real64 dPotOther_dP[MAX_STENCIL]{};
        real64 dPotOther_dC[MAX_STENCIL][NC]{};
        real64 dPropOther_dC[NC]{};

        //Fetch pot for phase j!=i defined as \rho_j g dz/dx
        UpwindHelpers::formPotential< NC, T, NUM_ELEMS, MAX_STENCIL >::compute(
          numPhase,
          jp,
          stencilSize,
          seri,
          sesri,
          sei,
          stencilWeights,
          totFlux_unw,
          gravCoef,
          dCompFrac_dCompDens,
          phaseMassDens,
          dPhaseMassDens_dPres,
          dPhaseMassDens_dComp,
          dPhaseVolFrac_dPres,
          dPhaseVolFrac_dComp,
          phaseCapPressure,
          dPhaseCapPressure_dPhaseVolFrac,
          potOther,
          dPotOther_dP,
          dPotOther_dC,
          dPropOther_dC
          );

        //Eventually get the mobility of the second phase
        real64 mobOther{};
        real64 dMobOther_dP{};
        real64 dMobOther_dC[NC]{};

        // and the other mobility for gravitational part as \lambda_j^{up} with up decided upon
        // the Upwind strategy - Note that it should be the same as the gravitational fractional flow

        UpwindHelpers::upwindMob< NC, NUM_ELEMS, MAX_STENCIL, T, UPWIND >( numPhase,
                                                                           jp,
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
                                                                           k_up_o,
                                                                           mobOther,
                                                                           dMobOther_dP,
                                                                           dMobOther_dC );


//          std::cout << k_up_og << " ";

        // Assembling gravitational flux phase-wise as \phi_{i,g} = \sum_{k\nei} \lambda_k^{up,g} f_k^{up,g} (G_i - G_k)
        phaseFlux -= fflow * mobOther * ( pot - potOther );
        dPhaseFlux_dP[k_up_o] -= fflow * dMobOther_dP * ( pot - potOther );
        for( localIndex jc = 0; jc < NC; ++jc )
          dPhaseFlux_dC[k_up_o][jc] -= fflow * dMobOther_dC[jc] * ( pot - potOther );

        //mob related part of dFflow_dP is only upstream defined but totMob related is defined everywhere
        for( localIndex ke = 0; ke < stencilSize; ++ke )
        {
          dPhaseFlux_dP[ke] -= dFflow_dP[ke] * mobOther * ( pot - potOther );

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= dFflow_dC[ke][jc] * mobOther * ( pot - potOther );
          }
        }

        for( localIndex ke = 0; ke < SIZE; ++ke )
        {
          dPhaseFlux_dP[ke] -= fflow * mobOther * ( dPot_dP[ke] - dPotOther_dP[ke] );
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dPhaseFlux_dC[ke][jc] -= fflow * mobOther * ( dPot_dC[ke][jc] - dPotOther_dC[ke][jc] );
          }
        }
      }
    }

  }

  /**
   * @brief  Struct defining formation of potential from different Physics (flagged by enum type T) to be used
   *            in Upwind discretization schemes
   * @tparam NC
   * @tparam T the concerned physics (Viscou,Graviaty or Capillary)
   * @tparam NUM_ELEMS
   * @tparam MAX_STENCIL
   */
  template< localIndex NC, term T, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  struct  formPotential
  {

    GEOSX_HOST_DEVICE
    static void compute( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const stencilSize,
                         arraySlice1d< localIndex const > const seri,
                         arraySlice1d< localIndex const > const sesri,
                         arraySlice1d< localIndex const > const sei,
                         arraySlice1d< real64 const > const stencilWeights,
                         real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                         ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                         ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                         ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                         ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                         ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                         real64 & potHead,
                         real64 ( & dPotHead_dPres)[MAX_STENCIL],
                         real64 (& dPotHead_dComp)[MAX_STENCIL][NC],
                         real64 ( & dProp_dComp)[NC] ) {};
  };
/*****/
// TODO specializing for viscous ?
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  struct formPotential< NC, term::Viscous, NUM_ELEMS, MAX_STENCIL >
  {

    GEOSX_HOST_DEVICE
    static void compute( localIndex const GEOSX_UNUSED_PARAM( numPhase ),
                         localIndex const GEOSX_UNUSED_PARAM( ip ),
                         localIndex const GEOSX_UNUSED_PARAM( stencilSize ),
                         arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM( seri ),
                         arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM( sesri ),
                         arraySlice1d< localIndex const > const GEOSX_UNUSED_PARAM( sei ),
                         arraySlice1d< real64 const > const GEOSX_UNUSED_PARAM( stencilWeights ),
                         real64 const totFlux,
                         ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravCoef ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dCompFrac_dCompDens ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseMassDens ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseMassDens_dPres ),
                         ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseMassDens_dComp ),
                         ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseVolFrac_dPres ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseVolFrac_dComp ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseCapPressure ),
                         ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                         real64 & pot,
                         real64( &GEOSX_UNUSED_PARAM( dPot_dPres ))[MAX_STENCIL],
                         real64( &GEOSX_UNUSED_PARAM( dPot_dComp ) )[MAX_STENCIL][NC],
                         real64( &GEOSX_UNUSED_PARAM( dProp_dComp ))[NC] )
    {
      pot = totFlux;
      //could be relevant for symmetry to include derivative

    }
  };

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  struct formPotential< NC, term::Gravity, NUM_ELEMS, MAX_STENCIL >
  {
/**
 * @brief Form gravitational head for phase from gravity and massDensities
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 * @param ip phase concerned
 * @param stencilSize number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 * @param stencilWeights weights associated with elements in the stencil
 */
    GEOSX_HOST_DEVICE
    static void compute( localIndex const GEOSX_UNUSED_PARAM( numPhase ),
                         localIndex const ip,
                         localIndex const stencilSize,
                         arraySlice1d< localIndex const >
                         const seri,
                         arraySlice1d< localIndex const > const sesri,
                         arraySlice1d< localIndex const > const sei,
                         arraySlice1d< real64 const > const stencilWeights,
                         real64 const GEOSX_UNUSED_PARAM( totFlux ), //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                         ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                         ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                         ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                         ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                         ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseVolFrac_dPres ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseVolFrac_dComp ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseCapPressure ),
                         ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ),
                         real64 & pot,
                         real64 ( & dPot_dPres )[NUM_ELEMS],
                         real64 (& dPot_dComp )[NUM_ELEMS][NC],
                         real64 ( & dProp_dComp )[NC] )
    {
      //working arrays
      real64 densMean{};
      real64 dDensMean_dPres[NUM_ELEMS]{};
      real64 dDensMean_dComp[NUM_ELEMS][NC]{};

      //init
      pot = 0.0;
      for( localIndex i = 0; i < MAX_STENCIL; ++i )
      {
        dPot_dPres[i] = 0.0;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPot_dComp[i][jc] = 0.0;
          dProp_dComp[jc] = 0.0;
        }
      }

      //inner loop to get average density
      for( localIndex i = 0; i < NUM_ELEMS; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];

        // density
        real64 const density = phaseMassDens[er][esr][ei][0][ip];
        real64 const dDens_dPres = dPhaseMassDens_dPres[er][esr][ei][0][ip];

        applyChainRule( NC,
                        dCompFrac_dCompDens[er][esr][ei],
                        dPhaseMassDens_dComp[er][esr][ei][0][ip],
                        dProp_dComp );

        // average density and derivatives
        densMean += 0.5 * density;
        dDensMean_dPres[i] = 0.5 * dDens_dPres;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dDensMean_dComp[i][jc] = 0.5 * dProp_dComp[jc];
        }
      }

      // compute potential difference MPFA-style
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];
        real64 const weight = stencilWeights[i];

        real64 const gravD = weight * gravCoef[er][esr][ei];
        pot += densMean * gravD;

        // need to add contributions from both cells the mean density depends on
        for( localIndex j = 0; j < NUM_ELEMS; ++j )
        {
          dPot_dPres[j] += dDensMean_dPres[j] * gravD;
          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dPot_dComp[j][jc] += dDensMean_dComp[j][jc] * gravD;
          }
        }
      }

    }
  };

  /**
   * @brief Form capillary head
   * @tparam NC
   * @tparam NUM_ELEMS
   * @param ip
   * @param stencilSize
   * @param seri
   * @param sesri
   * @param sei
   * @param stencilWeights
   */

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  struct formPotential< NC, term::Capillary, NUM_ELEMS, MAX_STENCIL >
  {

    GEOSX_HOST_DEVICE
    static void compute( localIndex const numPhase,
                         localIndex const ip,
                         localIndex const stencilSize,
                         arraySlice1d< localIndex const > const seri,
                         arraySlice1d< localIndex const > const sesri,
                         arraySlice1d< localIndex const > const sei,
                         arraySlice1d< real64 const > const stencilWeights,
                         real64 const GEOSX_UNUSED_PARAM( totFlux ),  //in fine should be a ElemnetViewConst once seq form are in place
                         ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( gravCoef ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dCompFrac_dCompDens ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseMassDens ),
                         ElementViewConst< arrayView3d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseMassDens_dPres ),
                         ElementViewConst< arrayView4d< real64 const > > const & GEOSX_UNUSED_PARAM( dPhaseMassDens_dComp ),
                         ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                         ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                         ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                         ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                         real64 & pot,
                         real64 ( & dPot_dPres)[MAX_STENCIL],
                         real64 (& dPot_dComp)[MAX_STENCIL][NC],
                         real64( &GEOSX_UNUSED_PARAM( dProp_dComp ))[NC] )
    {


      for( localIndex i = 0; i < stencilSize; ++i )
      {
        localIndex const er = seri[i];
        localIndex const esr = sesri[i];
        localIndex const ei = sei[i];
        real64 const weight = stencilWeights[i];

        pot += weight * phaseCapPressure[er][esr][ei][0][ip];
        // need to add contributions from both cells
        for( localIndex jp = 0; jp < numPhase; ++jp )
        {

          real64 const dCapPressure_dS = dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dPot_dPres[i] += weight * dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dPot_dComp[i][jc] += weight * dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
          }

        }

      }

    }
  };

/**
 * @brief Remultiply field by molar density (and derivatives)
 * @tparam NC number of components
 * @tparam MAX_STENCIL maximum number of points in the stencil
 * @param ip concerned phase index
 * @param k_up upwind direction of the phase
 * @param stencilSize  number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 * @param dCompFrac_dCompDens table of derivatives of composition fraction wrt their molar densities
 * @param phaseDens table of molar density by elements and phase
 * @param dPhaseDens_dPres table of molar density derivatives wrt pressure by elements and phase
 * @param dPhaseDens_dComp table of molar density derivatives wrt pressure by elements and phase
 * @param field a scalar field to be rescaled by molar density
 * @param dField_dPres the scalar field's pressure derivative table
 * @param dField_dComp the scalar field's component derivative table
 */

  template< localIndex NC, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void mdensMultiply( localIndex const ip,
                             localIndex const k_up,
                             localIndex const stencilSize,
                             arraySlice1d< localIndex const > const seri,
                             arraySlice1d< localIndex const > const sesri,
                             arraySlice1d< localIndex const > const sei,
                             ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                             ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                             ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                             ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
                             real64 & field,
                             real64 ( & dField_dPres)[MAX_STENCIL],
                             real64 (& dField_dComp)[MAX_STENCIL][NC]
                             )
  {
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dField_dPres[ke] = dField_dPres[ke] * phaseDens[er_up][esr_up][ei_up][0][ip];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dField_dComp[ke][jc] = dField_dComp[ke][jc] * phaseDens[er_up][esr_up][ei_up][0][ip];
      }
    }

    dField_dPres[k_up] += dPhaseDens_dPres[er_up][esr_up][ei_up][0][ip] * field;
    real64 dPhaseDens_dCompDens[NC] = { 0.0 };
    applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseDens_dComp[er_up][esr_up][ei_up][0][ip],
                    dPhaseDens_dCompDens );

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dField_dComp[k_up][jc] += dPhaseDens_dCompDens[jc] * field;
    }

    //last as multiplicative use in the second part of derivatives
    field = field * phaseDens[er_up][esr_up][ei_up][0][ip];
  }

/**
 * @brief Distribute phaseFlux onto component fluxes
 * @tparam NC number of components
 * @tparam MAX_STENCIL maximum number of points in the stencil
 * @param ip concerned phase index
 * @param k_up upwind direction of the phase
 * @param stencilSize number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 */
  template< localIndex NC, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static void formPhaseComp( localIndex const ip,
                             localIndex const k_up,
                             localIndex const stencilSize,
                             arraySlice1d< localIndex const > const seri,
                             arraySlice1d< localIndex const > const sesri,
                             arraySlice1d< localIndex const > const sei,
                             ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                             ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                             ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
                             ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                             real64 const & phaseFlux,
                             real64 const (&dPhaseFlux_dPres)[MAX_STENCIL],
                             real64 const (&dPhaseFlux_dComp)[MAX_STENCIL][NC],
                             real64 ( & compFlux)[NC],
                             real64 (& dCompFlux_dPres)[MAX_STENCIL][NC],
                             real64 ( & dCompFlux_dComp)[MAX_STENCIL][NC][NC] )
  {
    /*update phaseComp from grav part*/
    localIndex const er_up = seri[k_up];
    localIndex const esr_up = sesri[k_up];
    localIndex const ei_up = sei[k_up];

    arraySlice1d< real64 const > phaseCompFracSub = phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice1d< real64 const > dPhaseCompFrac_dPresSub = dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const > dPhaseCompFrac_dCompSub = dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][ip];

    real64 dProp_dC[NC]{};

    // compute component fluxes and derivatives using upstream cell composition
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dCompFlux_dPres[ke][ic] += dPhaseFlux_dPres[ke] * ycp;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dComp[ke][ic][jc] += dPhaseFlux_dComp[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dPres[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

      // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
      // densities
      applyChainRule( NC, dCompFrac_dCompDens[er_up][esr_up][ei_up], dPhaseCompFrac_dCompSub[ic], dProp_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dComp[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }



/*********** upwinded mobilities and fractional flow (upwind templated)  ************/

/**
 * @brief Function returning upwinded mobility (and derivatives)  of specified phase as well as uppwind direction
 *        according to specified upwind scheme
 * @tparam NC number of components
 * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
 * @param numPhase total number of phases
 * @param ip concerned phase index
 * @param stencilSizei number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 * @param stencilWeights weights associated with elements in the stencil
 * @param totFlux total flux signed value
 */
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, term T, template< term > class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  upwindMob( localIndex const numPhase,
             localIndex const ip,
             localIndex const stencilSize,
             arraySlice1d< localIndex const > const seri,
             arraySlice1d< localIndex const > const sesri,
             arraySlice1d< localIndex const > const sei,
             arraySlice1d< real64 const > const stencilWeights,
             real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
             ElementViewConst< arrayView1d< real64 const > > const & pres,
             ElementViewConst< arrayView1d< real64 const > > const & dPres,
             ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
             ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
             ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
             ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
             ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
             ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
             ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
             ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
             ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
             ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
             ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
             ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
             integer const capPressureFlag,
             localIndex & upwindDir,
             real64 & mob,
             real64 (&dMob_dP),
             real64 ( & dMob_dC)[NC] )
  {

    //reinit
    mob = 0.0;
    dMob_dP = 0.0;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dMob_dC[ic] = 0.0;
    }

    UPWIND< T > scheme;
    scheme.template getUpwindDir< NC, NUM_ELEMS, MAX_STENCIL, UPWIND >( numPhase,
                                                                        ip,
                                                                        stencilSize,
                                                                        seri,
                                                                        sesri,
                                                                        sei,
                                                                        stencilWeights,
                                                                        totFlux, //in fine should be a ElementViewConst once seq form are in
                                                                                 // place
                                                                        pres,
                                                                        dPres,
                                                                        gravCoef,
                                                                        phaseMob,
                                                                        dCompFrac_dCompDens,
                                                                        phaseMassDens,
                                                                        dPhaseMassDens_dPres,
                                                                        dPhaseMassDens_dComp,
                                                                        dPhaseVolFrac_dPres,
                                                                        dPhaseVolFrac_dComp,
                                                                        phaseCapPressure,
                                                                        dPhaseCapPressure_dPhaseVolFrac,
                                                                        capPressureFlag,
                                                                        upwindDir );

    localIndex const er_up = seri[upwindDir];
    localIndex const esr_up = sesri[upwindDir];
    localIndex const ei_up = sei[upwindDir];

    if( std::fabs( phaseMob[er_up][esr_up][ei_up][ip] ) > 1e-20 )
    {
      mob = phaseMob[er_up][esr_up][ei_up][ip];
      dMob_dP = dPhaseMob_dPres[er_up][esr_up][ei_up][ip];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dMob_dC[ic] = dPhaseMob_dComp[er_up][esr_up][ei_up][ip][ic];
      }
    }
  }

/**
 * @brief Function returning upwinded fractional flow (and derivatives) as well as upwind direction of specified phase
 *        according to specified upwind scheme
 * @tparam NC number of components
 * @tparam UpwindScheme Desciption of how to construct potential used to decide upwind direction
 * @param numPhase total number of phases
 * @param ip concerned phase index
 * @param stencilSizei number of points in the stencil
 * @param seri arraySlice of the stencil implied element region index
 * @param sesri arraySlice of the stencil implied element subregion index
 * @param sei arraySlice of the stencil implied element index
 * @param stencilWeights weights associated with elements in the stencil
 * @param totFlux total flux signed value
 */
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, term T, template< term > class UPWIND >
  GEOSX_HOST_DEVICE
  static void
  formFracFlow( localIndex const numPhase,
                localIndex const ip,
                localIndex const stencilSize,
                arraySlice1d< localIndex const > const seri,
                arraySlice1d< localIndex const > const sesri,
                arraySlice1d< localIndex const > const sei,
                arraySlice1d< real64 const > const stencilWeights,
                real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                ElementViewConst< arrayView1d< real64 const > > const & pres,
                ElementViewConst< arrayView1d< real64 const > > const & dPres,
                ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
                ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                integer const capPressureFlag,
                localIndex & k_up_main,
                real64 & fflow,
                real64 (& dFflow_dP)[MAX_STENCIL],
                real64 ( & dFflow_dC)[MAX_STENCIL][NC] )
  {
    // get var to memorized the numerator mobility properly upwinded
    real64 mainMob{};
    real64 dMMob_dP{};
    real64 dMMob_dC[NC]{};

    real64 totMob{};
    real64 dTotMob_dP[MAX_STENCIL]{};
    real64 dTotMob_dC[MAX_STENCIL][NC]{};

    //reinit
    //fractional flow too low to let the upstream phase flow
    k_up_main = -1; //to throw error if unmodified
    fflow = 0;
    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      dFflow_dP[ke] = 0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dFflow_dC[ke][jc] = 0;
      }
    }

    //Form totMob
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {

      localIndex k_up;
      real64 mob{};
      real64 dMob_dP{};
      real64 dMob_dC[NC]{};

      upwindMob< NC, NUM_ELEMS, MAX_STENCIL, T, UPWIND >( numPhase,
                                                          jp,
                                                          stencilSize,
                                                          seri,
                                                          sesri,
                                                          sei,
                                                          stencilWeights,
                                                          totFlux, //in fine should be a ElemnetViewConst once seq form are in place
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
                                                          mob,
                                                          dMob_dP,
                                                          dMob_dC );


      totMob += mob;
      dTotMob_dP[k_up] += dMob_dP;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dTotMob_dC[k_up][ic] += dMob_dC[ic];
      }

      if( jp == ip )
      {
        k_up_main = k_up;
        mainMob = mob;
        dMMob_dP = dMob_dP;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dMMob_dC[ic] = dMob_dC[ic];
        }
      }
    }

    //guard against no flow region
    if( std::fabs( mainMob ) > 1e-20 )
    {
      fflow = mainMob / totMob;
      dFflow_dP[k_up_main] = dMMob_dP / totMob;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dFflow_dC[k_up_main][jc] = dMMob_dC[jc] / totMob;

      }

      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dFflow_dP[ke] -= fflow * dTotMob_dP[ke] / totMob;

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dFflow_dC[ke][jc] -= fflow * dTotMob_dC[ke][jc] / totMob;
        }
      }
    }
  }

};

/************************* UPWIND ******************/

/**
 * @brief Template base class for different upwind Scheme
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 * @tparam T physics concerned by the scheme if specialized
 */
template< term T >
class UpwindScheme
{

public:

  //default ctor
  UpwindScheme() = default;

  //usual copy ctor
  UpwindScheme( UpwindScheme const & scheme ) = default;

  //default move ctor
  UpwindScheme( UpwindScheme && ) = default;

  //deleted copy and move assignement
  UpwindScheme & operator=( UpwindScheme const & ) = delete;

  UpwindScheme & operator=( UpwindScheme && ) = delete;

  virtual ~UpwindScheme() = default;

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, template< term > class UPWIND >
  GEOSX_HOST_DEVICE
  void getUpwindDir( localIndex const numPhase,
                     localIndex const ip,
                     localIndex const stencilSize,
                     arraySlice1d< localIndex const > const seri,
                     arraySlice1d< localIndex const > const sesri,
                     arraySlice1d< localIndex const > const sei,
                     arraySlice1d< real64 const > const stencilWeights,
                     real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                     ElementViewConst< arrayView1d< real64 const > > const & pres,
                     ElementViewConst< arrayView1d< real64 const > > const & dPres,
                     ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                     ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                     ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                     ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                     ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                     ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                     ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                     ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                     ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                     ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                     integer const capPressureFlag,
                     localIndex & upwindDir
                     )
  {
    real64 pot{};
    localIndex source{};

    UPWIND< T >::template calcPotential< NC, NUM_ELEMS, MAX_STENCIL >( numPhase,
                                                                       ip,
                                                                       stencilSize,
                                                                       seri,
                                                                       sesri,
                                                                       sei,
                                                                       stencilWeights,
                                                                       totFlux, //in fine should be a ElemnetViewConst once seq form are in
                                                                                // place
                                                                       pres,
                                                                       dPres,
                                                                       gravCoef,
                                                                       phaseMob,
                                                                       dCompFrac_dCompDens,
                                                                       phaseMassDens,
                                                                       dPhaseMassDens_dPres,
                                                                       dPhaseMassDens_dComp,
                                                                       dPhaseVolFrac_dPres,
                                                                       dPhaseVolFrac_dComp,
                                                                       phaseCapPressure,
                                                                       dPhaseCapPressure_dPhaseVolFrac,
                                                                       capPressureFlag,
                                                                       source,
                                                                       pot );

    //all definition has been changed to fit pot>0 => first cell is upstream
    upwindDir = ( pot > 0 ) ? source : ( ( source == 0 ) ? 1 : 0 );
  }

  //refactor getPotential - 3 overload// by phase

  //TODO Settle passing lambda MAX_STENCIL or NUM_ELEMS
  template< localIndex NC, localIndex SIZE, typename LAMBDA >
  GEOSX_HOST_DEVICE
  static void potential( localIndex numPhase,
                         localIndex ip,
                         arraySlice1d< localIndex const > const seri,
                         arraySlice1d< localIndex const > const sesri,
                         arraySlice1d< localIndex const > const sei,
                         ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                         real64 & wpot,
                         LAMBDA && fn )
  {
    //getPhase Pot
    real64 pot {};
    real64 pot_dP[SIZE] {};
    real64 pot_dC[SIZE][NC] {};
    real64 dProp_dC[NC] {};

    fn( ip, pot, pot_dP, pot_dC, dProp_dC );

    localIndex const k_up = 0;
    localIndex const k_dw = 1;

    //loop other other phases to form
    for( localIndex jp = 0; jp < numPhase; ++jp )
    {
      if( jp != ip )
      {
        localIndex const er_up = seri[k_up];
        localIndex const esr_up = sesri[k_up];
        localIndex const ei_up = sei[k_up];

        localIndex const er_dw = seri[k_dw];
        localIndex const esr_dw = sesri[k_dw];
        localIndex const ei_dw = sei[k_dw];

        real64 potOther {};
        real64 potOther_dP[SIZE] {};
        real64 potOther_dC[SIZE][NC] {};
        real64 dPropOther_dC[NC] {};

        fn( jp, potOther, potOther_dP, potOther_dC, dPropOther_dC );

        real64 const mob_up = phaseMob[er_up][esr_up][ei_up][jp];
        real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][jp];

        wpot +=  ( pot - potOther >= 0 ) ? mob_dw * ( potOther - pot ) : mob_up * ( potOther - pot );

      }
      ;
    }
  }

};

/**
 * @brief Class describing the classical Phase Potential Upwind Scheme as studied in Sammon. "An analysis of upstream
 *        differencing." SPE reservoir engineering (1988)
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 * @tparam T physics concerned by the scheme if specialized
 */
template< term T >
class PhasePotentialUpwind : public UpwindScheme< T >
{
public:

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static
  void calcPotential( localIndex const numPhase,
                      localIndex const ip,
                      localIndex const stencilSize,
                      arraySlice1d< localIndex const > const seri,
                      arraySlice1d< localIndex const > const sesri,
                      arraySlice1d< localIndex const > const sei,
                      arraySlice1d< real64 const > const stencilWeights,
                      real64 const totFlux,
                      ElementViewConst< arrayView1d< real64 const > > const & pres,
                      ElementViewConst< arrayView1d< real64 const > > const & dPres,
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseMob ),
                      ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                      ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                      integer const capPressureFlag,
                      localIndex & GEOSX_UNUSED_PARAM( source ),
                      real64 & pot
                      )
  {
    //compute presGrad
    real64 presGrad{};

    for( localIndex i = 0; i < stencilSize; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];
      real64 const weight = stencilWeights[i];

      //TODO add capillary
      presGrad += weight * ( pres[er][esr][ei] + dPres[er][esr][ei] );

      if( capPressureFlag ) // could also rely on capHead
        presGrad -= weight * phaseCapPressure[er][esr][ei][0][ip];

    }

    // then form GravHead
    real64 gravHead{};
    real64 dGravHead_dP[MAX_STENCIL]{};
    real64 dGravHead_dC[MAX_STENCIL][NC]{};
    real64 dProp_dC[NC]{};

    UpwindHelpers::formPotential< NC, term::Gravity, NUM_ELEMS, MAX_STENCIL >::compute(
      numPhase,
      ip,
      stencilSize,
      seri,
      sesri,
      sei,
      stencilWeights,
      totFlux,
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens_dPres,
      dPhaseMassDens_dComp,
      dPhaseVolFrac_dPres,
      dPhaseVolFrac_dComp,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      gravHead,
      dGravHead_dP,
      dGravHead_dC,
      dProp_dC );



    // presGrad and gravHead are including (-K) by weithing , will then fit pot>0 => first-cell upwind
    pot = presGrad - gravHead;
  }

};

/**
 * @brief Class describing Phase Potential Upwind scheme as defined in "Upstream differencing for multiphase flow
 *           in reservoir simulation " (Brenier and Jaffre,1986)
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
template< term T >
class PhaseUpwind : public UpwindScheme< T >
{

public:

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static
  void calcPotential( localIndex const numPhase,
                      localIndex const ip,
                      localIndex const stencilSize,
                      arraySlice1d< localIndex const > const seri,
                      arraySlice1d< localIndex const > const sesri,
                      arraySlice1d< localIndex const > const sei,
                      arraySlice1d< real64 const > const stencilWeights,
                      real64 const totFlux, //in fine should be a ElemnetViewConst once seq form are in place
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                      ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                      ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                      integer const capPressureFlag,
                      localIndex & GEOSX_UNUSED_PARAM( source ),
                      real64 & pot
                      )
  {
    // Allow different treatment
    real64 dPot_dP[MAX_STENCIL] {};
    real64 dPot_dC[MAX_STENCIL][NC] {};
    real64 dProp_dC[NC] {};

    //Form total velocity
    UpwindHelpers::formPotential< NC, term::Viscous, NUM_ELEMS, MAX_STENCIL >::compute(
      numPhase,
      ip,
      stencilSize,
      seri,
      sesri,
      sei,
      stencilWeights,
      totFlux, //in fine should be a ElemnetViewConst once seq form are in place
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens_dPres,
      dPhaseMassDens_dComp,
      dPhaseVolFrac_dPres,
      dPhaseVolFrac_dComp,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      pot,
      dPot_dP,
      dPot_dC,
      dProp_dC );


    UpwindScheme< T >::template potential< NC, MAX_STENCIL >( numPhase, ip, seri, sesri, sei, phaseMob, pot,
                                                              [&]( localIndex ipp, real64 & potk, real64 (& dPotk_dP)[MAX_STENCIL], real64 (& dPotk_dC)[MAX_STENCIL][NC],
                                                                   real64 (& dPropk_dC)[NC] ){

          real64 potu = 0;

          UpwindHelpers::formPotential< NC, term::Gravity, NUM_ELEMS, MAX_STENCIL >::compute(
            numPhase,
            ipp,
            stencilSize,
            seri,
            sesri,
            sei,
            stencilWeights,
            totFlux,                                                   //in fine should be a ElemnetViewConst once seq form are in place
            gravCoef,
            dCompFrac_dCompDens,
            phaseMassDens,
            dPhaseMassDens_dPres,
            dPhaseMassDens_dComp,
            dPhaseVolFrac_dPres,
            dPhaseVolFrac_dComp,
            phaseCapPressure,
            dPhaseCapPressure_dPhaseVolFrac,
            potk,
            dPotk_dP,
            dPotk_dC,
            dPropk_dC );



          potk += potu;
          if( capPressureFlag )
          {
            potu = 0;
            UpwindHelpers::formPotential< NC, term::Capillary, NUM_ELEMS, MAX_STENCIL >::compute(
              numPhase,
              ipp,
              stencilSize,
              seri,
              sesri,
              sei,
              stencilWeights,
              totFlux,                                                   //in fine should be a ElemnetViewConst once seq form are in place
              gravCoef,
              dCompFrac_dCompDens,
              phaseMassDens,
              dPhaseMassDens_dPres,
              dPhaseMassDens_dComp,
              dPhaseVolFrac_dPres,
              dPhaseVolFrac_dComp,
              phaseCapPressure,
              dPhaseCapPressure_dPhaseVolFrac,
              potk,
              dPotk_dP,
              dPotk_dC,
              dPropk_dC );

            potk += potu;
          }



        } );
  }

};

/**
 * @brief  Class describing the Hybrid Upwind scheme as defined in "Consistent upwinding for sequential fully implicit
 *         multiscale compositional simulation" (Moncorge,2020)
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
template< term T >
class HybridUpwind : public UpwindScheme< T >
{

public:
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static
  void calcPotential( localIndex const numPhase,
                      localIndex const ip,
                      localIndex const stencilSize,
                      arraySlice1d< localIndex const > const seri,
                      arraySlice1d< localIndex const > const sesri,
                      arraySlice1d< localIndex const > const sei,
                      arraySlice1d< real64 const > const stencilWeights,
                      real64 const totFlux,
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                      ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                      ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                      integer const GEOSX_UNUSED_PARAM( capPressureFlag ),
                      localIndex & GEOSX_UNUSED_PARAM( source ),
                      real64 & pot
                      )
  {
    //Form total velocity
    pot = 0;


    //TODO check if MAX_STENCIL or NUM_ELEMS
    UpwindScheme< T >::template potential< NC, MAX_STENCIL >( numPhase, ip, seri, sesri, sei, phaseMob, pot,
                                                              [&]( localIndex ipp,
                                                                   real64 & potk,
                                                                   real64 (& dPotk_dP)[MAX_STENCIL],
                                                                   real64 (& dPotk_dC)[MAX_STENCIL][NC],
                                                                   real64 (& dProp_dC)[NC] )
        {

          UpwindHelpers::formPotential< NC, T, NUM_ELEMS, MAX_STENCIL >::compute(
            numPhase,
            ipp,
            stencilSize,
            seri,
            sesri,
            sei,
            stencilWeights,
            totFlux,                                              //in fine should be a ElemnetViewConst once seq form are in place
            gravCoef,
            dCompFrac_dCompDens,
            phaseMassDens,
            dPhaseMassDens_dPres,
            dPhaseMassDens_dComp,
            dPhaseVolFrac_dPres,
            dPhaseVolFrac_dComp,
            phaseCapPressure,
            dPhaseCapPressure_dPhaseVolFrac,
            potk,
            dPotk_dP,
            dPotk_dC,
            dProp_dC );

        } );
  }
//    // reinit source and target depending on gravity orientation
//
//    //form gravHead for the currentPhase
//    // then form GravHead
//    real64 gravHead{};
//    real64 dGravHead_dP[NUM_ELEMS]{};
//    real64 dGravHead_dC[NUM_ELEMS][NC]{};
//    real64 dProp_dC[NC]{};
//
//    UpwindHelpers::formGravHead( ip,
//                                 stencilSize,
//                                 seri,
//                                 sesri,
//                                 sei,
//                                 stencilWeights,
//                                 gravCoef,
//                                 dCompFrac_dCompDens,
//                                 phaseMassDens,
//                                 dPhaseMassDens_dPres,
//                                 dPhaseMassDens_dComp,
//                                 gravHead,
//                                 dGravHead_dP,
//                                 dGravHead_dC,
//                                 dProp_dC );
//
//    localIndex const k_up = 0;
//    localIndex const k_dw = 1;
//
//    //loop other other phases to form
//    for( localIndex jp = 0; jp < numPhase; ++jp )
//    {
//
//      if( jp != ip )
//      {
//        // then form GravHead
//        real64 gravHeadOther{};
//        real64 dGravHeadOther_dP[NUM_ELEMS]{};
//        real64 dGravHeadOther_dC[NUM_ELEMS][NC]{};
//        real64 dPropOther_dC[NC]{};
//
//        localIndex const er_up = seri[k_up];
//        localIndex const esr_up = sesri[k_up];
//        localIndex const ei_up = sei[k_up];
//
//        localIndex const er_dw = seri[k_dw];
//        localIndex const esr_dw = sesri[k_dw];
//        localIndex const ei_dw = sei[k_dw];
//
//        UpwindHelpers::formGravHead( jp,
//                                     stencilSize,
//                                     seri,
//                                     sesri,
//                                     sei,
//                                     stencilWeights,
//                                     gravCoef,
//                                     dCompFrac_dCompDens,
//                                     phaseMassDens,
//                                     dPhaseMassDens_dPres,
//                                     dPhaseMassDens_dComp,
//                                     gravHeadOther,
//                                     dGravHeadOther_dP,
//                                     dGravHeadOther_dC,
//                                     dPropOther_dC );
//
//        real64 const mob_up = phaseMob[er_up][esr_up][ei_up][jp];
//        real64 const mob_dw = phaseMob[er_dw][esr_dw][ei_dw][jp];
//
//        pot += ( gravHead - gravHeadOther >= 0 ) ? mob_dw * ( gravHeadOther - gravHead ) :
//               mob_up * ( gravHeadOther - gravHead );
//
//
//      }
//    }

};
//Special case for Viscous as the assembly of potential is different

/**
 * @brief Specialization of the Viscous term as it does not fit the generic framework summing over other phases
 * @tparam NC number of components
 * @tparam NUM_ELEMS number of elements neighbors of considered face
 */
template<>
class HybridUpwind< term::Viscous > : public UpwindScheme< term::Viscous >
{

public:
  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  static
  void calcPotential( localIndex const numPhase,
                      localIndex const ip,
                      localIndex const stencilSize,
                      arraySlice1d< localIndex const > const seri,
                      arraySlice1d< localIndex const > const sesri,
                      arraySlice1d< localIndex const > const sei,
                      arraySlice1d< real64 const > const stencilWeights,
                      real64 const totFlux,
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( pres ),
                      ElementViewConst< arrayView1d< real64 const > > const & GEOSX_UNUSED_PARAM( dPres ),
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( phaseMob ),
                      ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dComp,
                      ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
                      ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
                      ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
                      ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
                      integer const GEOSX_UNUSED_PARAM( capPressureFlag ),
                      localIndex & GEOSX_UNUSED_PARAM( source ),
                      real64 & pot
                      )
  {
    real64 dPot_dP[MAX_STENCIL] {};
    real64 dPot_dC[MAX_STENCIL][NC]{};
    real64 dProp_dC[NC] {};


    UpwindHelpers::formPotential< NC, term::Viscous, NUM_ELEMS, MAX_STENCIL >::compute(
      numPhase,
      ip,
      stencilSize,
      seri,
      sesri,
      sei,
      stencilWeights,
      totFlux,                                                    //in fine should be a ElemnetViewConst once seq form are in place
      gravCoef,
      dCompFrac_dCompDens,
      phaseMassDens,
      dPhaseMassDens_dPres,
      dPhaseMassDens_dComp,
      dPhaseVolFrac_dPres,
      dPhaseVolFrac_dComp,
      phaseCapPressure,
      dPhaseCapPressure_dPhaseVolFrac,
      pot,
      dPot_dP,
      dPot_dC,
      dProp_dC );
  }

};

}//end namespace CompositionalMultiphaseFlowUpwindHelperKernels

}//end namespace geosx
#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEFLOWUPWINDHELPERKERNELS_HP
