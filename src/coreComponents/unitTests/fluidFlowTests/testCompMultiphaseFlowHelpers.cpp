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

// Source includes
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mainInterface/initialization.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "physicsSolvers/PhysicsSolverManager.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVM.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFVMKernels.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowUpwindHelperKernels.hpp"
#include "unitTests/fluidFlowTests/testCompFlowUtils.hpp"

#include "unitTests/fluidFlowTests/testFlowKernelHelpers.hpp"
// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::CompositionalMultiphaseFlowUpwindHelperKernels;
using UHelpers = geosx::CompositionalMultiphaseFlowUpwindHelperKernels::UpwindHelpers;

CommandLineOptions g_commandLineOptions;

// Sphinx start after input XML
template< localIndex NC, localIndex NP, localIndex MAX_STENCIL >
void mdensMultiplyTest( localIndex const ip,
                        localIndex const k_up,
                        localIndex const stencilSize,
                        real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                        real64 const (*phaseDens)[MAX_STENCIL][1][NP],
                        real64 const (*dPhaseDens_dP)[MAX_STENCIL][1][NP],
                        real64 const (*dPhaseDens_dC)[MAX_STENCIL][1][NP][NC],
                        real64 & field,
                        real64 (& dField_dP)[MAX_STENCIL],
                        real64 (& dField_dC)[MAX_STENCIL][NC] )
{
  for( localIndex ke = 0; ke < stencilSize; ++ke )
  {
    dField_dP[ke] = dField_dP[ke] * ( *phaseDens )[k_up][0][ip];
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dField_dC[ke][jc] = dField_dC[ke][jc] * ( *phaseDens )[k_up][0][ip];
    }
  }

  dField_dP[k_up] += ( *dPhaseDens_dP )[k_up][0][ip] * field;
  real64 dPhaseDens_dCompDens[NC] = { 0.0 };
  applyChainRule( NC, ( *dCompFrac_dCompDens )[k_up], ( *dPhaseDens_dC )[k_up][0][ip],
                  dPhaseDens_dCompDens );

  for( localIndex jc = 0; jc < NC; ++jc )
  {
    dField_dC[k_up][jc] += dPhaseDens_dCompDens[jc] * field;
  }

//last as multiplicative use in the second part of derivatives
  field = field * ( *phaseDens )[k_up][0][ip];
}

template< localIndex NC, localIndex NP, localIndex MAX_STENCIL >
void formPhaseCompTest( localIndex const ip,
                        localIndex const k_up,
                        localIndex const stencilSize,
                        real64 const (*phaseCompFrac)[MAX_STENCIL][1][NP][NC],
                        real64 const (*dPhaseCompFrac_dPres)[MAX_STENCIL][1][NP][NC],
                        real64 const (*dPhaseCompFrac_dComp)[MAX_STENCIL][1][NP][NC][NC],
                        real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                        real64 const & phaseFlux,
                        real64 const (&dPhaseFlux_dPres)[MAX_STENCIL],
                        real64 const (&dPhaseFlux_dComp)[MAX_STENCIL][NC],
                        real64 (& compFlux)[NC],
                        real64 (& dCompFlux_dPres)[MAX_STENCIL][NC],
                        real64 (& dCompFlux_dComp)[MAX_STENCIL][NC][NC] )
{

  real64 dProp_dC[NC] = { 0.0 };

  // compute component fluxes and derivatives using upstream cell composition
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    real64 const ycp = ( *phaseCompFrac )[k_up][0][ip][ic];
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
    dCompFlux_dPres[k_up][ic] += phaseFlux * ( *dPhaseCompFrac_dPres )[k_up][0][ip][ic];

    // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
    // densities
    applyChainRule( NC, ( *dCompFrac_dCompDens )[k_up], ( *dPhaseCompFrac_dComp )[k_up][0][ip][ic], dProp_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dCompFlux_dComp[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
    }
  }

}

// --
template< localIndex NC, localIndex NP, term T, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
struct formPotTest
{

  static void compute( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const stencilSize,
                       real64 const (*weights),
                       real64 const (*gravCoef)[MAX_STENCIL],
                       real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                       real64 const (*phaseMassDens)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseMassDens_dP)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseMassDens_dC)[MAX_STENCIL][1][NP][NC],
                       real64 const (*dPhaseVolFrac_dP)[MAX_STENCIL][NP],
                       real64 const (*dPhaseVolFrac_dC)[MAX_STENCIL][NP][NC],
                       real64 const (*phaseCapPressure)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseCapPressure_dPhaseVolFrac)[MAX_STENCIL][1][NP][NP],
                       real64 & pot,
                       real64 (& dPot_dP)[NUM_ELEMS],
                       real64 (& dPot_dC)[NUM_ELEMS][NC],
                       real64 (& dProp_dC)[NC] )
  { }
};

// --
template< localIndex NC, localIndex NP, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
struct formPotTest< NC, NP, term::Gravity, NUM_ELEMS, MAX_STENCIL >
{
  static void compute( localIndex const GEOSX_UNUSED_PARAM( numPhase ),
                       localIndex const ip,
                       localIndex const stencilSize,
                       real64 const (*weights),
                       real64 const (*gravCoef)[MAX_STENCIL],
                       real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                       real64 const (*phaseMassDens)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseMassDens_dP)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseMassDens_dC)[MAX_STENCIL][1][NP][NC],
                       real64 const (*GEOSX_UNUSED_PARAM( dPhaseVolFrac_dP ))[MAX_STENCIL][NP],
                       real64 const (*GEOSX_UNUSED_PARAM( dPhaseVolFrac_dC ))[MAX_STENCIL][NP][NC],
                       real64 const (*GEOSX_UNUSED_PARAM( phaseCapPressure ))[MAX_STENCIL][1][NP],
                       real64 const (*GEOSX_UNUSED_PARAM( dPhaseCapPressure_dPhaseVolFrac ))[MAX_STENCIL][1][NP][NP],
                       real64 & pot,
                       real64 (& dPot_dP)[NUM_ELEMS],
                       real64 (& dPot_dC)[NUM_ELEMS][NC],
                       real64 (& dProp_dC)[NC] )
  {
    //working arrays
    real64 densMean{};
    real64 dDensMean_dP[NUM_ELEMS]{};
    real64 dDensMean_dC[NUM_ELEMS][NC]{};

    //init
    pot = 0.0;
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      dPot_dP[i] = 0.0;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dPot_dC[i][jc] = 0.0;
        dProp_dC[jc] = 0.0;
      }
    }

    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      // density
      real64 const density = ( *phaseMassDens )[i][0][ip];
      real64 const dDens_dP = ( *dPhaseMassDens_dP )[i][0][ip];

      applyChainRule( NC,
                      ( *dCompFrac_dCompDens )[i],
                      ( *dPhaseMassDens_dC )[i][0][ip],
                      dProp_dC );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dP;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
      }
    }

    // compute potential diffeirence MPFA-style
    for( localIndex i = 0; i < stencilSize; ++i )
    {
      real64 const weight = weights[i];

      real64 const gravD = weight * ( *gravCoef )[i];
      pot += densMean * gravD;

      // need to add contributions from both cells the mean density depends on
      for( localIndex j = 0; j < NUM_ELEMS; ++j )
      {
        dPot_dP[j] += dDensMean_dP[j] * gravD;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPot_dC[j][jc] += dDensMean_dC[j][jc] * gravD;
        }
      }
    }
  }

};

// --
template< localIndex NC, localIndex NP, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
struct formPotTest< NC, NP, term::Capillary, NUM_ELEMS, MAX_STENCIL >
{
  static void compute( localIndex const numPhase,
                       localIndex const ip,
                       localIndex const stencilSize,
                       real64 const (*weights),
                       real64 const (*GEOSX_UNUSED_PARAM( gravCoef ))[MAX_STENCIL],
                       real64 const (*GEOSX_UNUSED_PARAM( dCompFrac_dCompDens ))[MAX_STENCIL][NC][NC],
                       real64 const (*GEOSX_UNUSED_PARAM( phaseMassDens ))[MAX_STENCIL][1][NP],
                       real64 const (*GEOSX_UNUSED_PARAM( dPhaseMassDens_dP ))[MAX_STENCIL][1][NP],
                       real64 const (*GEOSX_UNUSED_PARAM( dPhaseMassDens_dC ))[MAX_STENCIL][1][NP][NC],
                       real64 const (*dPhaseVolFrac_dP)[MAX_STENCIL][NP],
                       real64 const (*dPhaseVolFrac_dC)[MAX_STENCIL][NP][NC],
                       real64 const (*phaseCapPressure)[MAX_STENCIL][1][NP],
                       real64 const (*dPhaseCapPressure_dPhaseVolFrac)[MAX_STENCIL][1][NP][NP],
                       real64 & pot,
                       real64 (& dPot_dP)[NUM_ELEMS],
                       real64 (& dPot_dC)[NUM_ELEMS][NC],
                       real64 (&GEOSX_UNUSED_PARAM( dProp_dC ))[NC] )
  {
    for( localIndex i = 0; i < stencilSize; ++i )
    {

      pot += weights[i] * (*phaseCapPressure)[i][0][ip];
      // need to add contributions from both cells
      for( localIndex jp = 0; jp < numPhase; ++jp )
      {

        real64 const dCapPressure_dS = (*dPhaseCapPressure_dPhaseVolFrac)[i][0][ip][jp];
        dPot_dP[i] += weights[i] * dCapPressure_dS * (*dPhaseVolFrac_dP)[i][jp];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dPot_dC[i][jc] += weights[i] * dCapPressure_dS * (*dPhaseVolFrac_dC)[i][jp][jc];
        }

      }
    }
  }
};

template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, localIndex NDOF >
void fillLocalJacobiTest( real64 const (&compFlux)[NC],
                          real64 const (&dCompFlux_dP)[MAX_STENCIL][NC],
                          real64 const (&dCompFlux_dC)[MAX_STENCIL][NC][NC],
                          localIndex const stencilSize,
                          real64 const dt,
                          real64 (& localFlux)[NUM_ELEMS * NC],
                          real64 (& localFluxJacobian)[NUM_ELEMS * NC][MAX_STENCIL * NDOF] )
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
// comparators

template< localIndex NC, localIndex NP, localIndex MAX_STENCIL, bool FULL >
void testCompositionalUpwindDensMult( CellElementStencilTPFA const & stencil,
                                      real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                                      real64 const (*phaseDens)[MAX_STENCIL][1][NP],
                                      real64 const (*dPhaseDens_dP)[MAX_STENCIL][1][NP],
                                      real64 const (*dPhaseDens_dC)[MAX_STENCIL][1][NP][NC] )
{
  localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();

  //some type rearganging
  auto phaseDensView = AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( phaseDens[0][0][0][0] ),
                                                                                  MAX_STENCIL,
                                                                                  seri[0],
                                                                                  sesri[0],
                                                                                  sei[0], 1, NP );

  auto dPhaseDens_dPView = AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dPhaseDens_dP[0][0][0][0] ),
                                                                                      MAX_STENCIL,
                                                                                      seri[0],
                                                                                      sesri[0],
                                                                                      sei[0], 1, NP );

  auto dPhaseDens_dCView = AccessorHelper< FULL >::template makeElementAccessor< 4 >( &( dPhaseDens_dC[0][0][0][0][0] ),
                                                                                      MAX_STENCIL,
                                                                                      seri[0],
                                                                                      sesri[0],
                                                                                      sei[0], 1, NP, NC );

  auto dCompFrac_dCompDensView =
    AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dCompFrac_dCompDens[0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], NC, NC );
  for( localIndex ke = 0; ke < MAX_STENCIL; ++ke )
  {
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 field = 5.650e+03;
      real64 dField_dP[NUM_ELEMS] = { 6.952e-03, 4.991e+00 };
      real64 dField_dC[NUM_ELEMS][NC] = { { 5.523e-01, 6.299e-03, 0.320e+00, 6.147e-02 },
        { 1.231e-03, 2.055e+01, 1.465e-02, 1.891e+00 } };

      UHelpers::mdensMultiply( ip,
                               ke,
                               MAX_STENCIL,
                               seri[0],
                               sesri[0],
                               sei[0],
                               dCompFrac_dCompDensView.toNestedViewConst(),
                               phaseDensView.toNestedViewConst(),
                               dPhaseDens_dPView.toNestedViewConst(),
                               dPhaseDens_dCView.toNestedViewConst(),
                               field,
                               dField_dP,
                               dField_dC );

      real64 fieldRef = 5.650e+03;
      real64 dFieldRef_dP[NUM_ELEMS] = { 6.952e-03, 4.991e+00 };
      real64 dFieldRef_dC[NUM_ELEMS][NC] = { { 5.523e-01, 6.299e-03, 0.320e+00, 6.147e-02 },
        { 1.231e-03, 2.055e+01, 1.465e-02, 1.891e+00 } };


      mdensMultiplyTest( ip,
                         ke,
                         MAX_STENCIL,
                         dCompFrac_dCompDens,
                         phaseDens,
                         dPhaseDens_dP,
                         dPhaseDens_dC,
                         fieldRef,
                         dFieldRef_dP,
                         dFieldRef_dC );

      EXPECT_DOUBLE_EQ( fieldRef, field );
      for( localIndex kl = 0; kl < MAX_STENCIL; ++kl )
      {

        EXPECT_DOUBLE_EQ( dFieldRef_dP[kl], dField_dP[kl] );
        for( localIndex ic = 0; ic < NC; ++ic )
          EXPECT_DOUBLE_EQ( dFieldRef_dC[kl][ic], dField_dC[kl][ic] );
      }

    }
  }


}

template< localIndex NC, localIndex NP, localIndex MAX_STENCIL, bool FULL >
void testCompositionalUpwindFormPhaseComp( CellElementStencilTPFA const & stencil,
                                           real64 const (*dCompFrac_dCompDens )[MAX_STENCIL][NC][NC],
                                           real64 const (*phaseCompFrac )[MAX_STENCIL][1][NP][NC],
                                           real64 const (*dPhaseCompFrac_dPres )[MAX_STENCIL][1][NP][NC],
                                           real64 const (*dPhaseCompFrac_dComp )[MAX_STENCIL][1][NP][NC][NC] )
{
  localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();

  //some type rearganging
  auto phaseCompFracView = AccessorHelper< FULL >::template makeElementAccessor< 4 >( &( phaseCompFrac[0][0][0][0][0] ),
                                                                                      MAX_STENCIL,
                                                                                      seri[0],
                                                                                      sesri[0],
                                                                                      sei[0], 1, NP, NC );

  auto dPhaseCompFrac_dPView =
    AccessorHelper< FULL >::template makeElementAccessor< 4 >( &( dPhaseCompFrac_dPres[0][0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], 1, NP, NC );

  auto dPhaseCompFrac_dCView =
    AccessorHelper< FULL >::template makeElementAccessor< 5 >( &( dPhaseCompFrac_dComp[0][0][0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], 1, NP, NC, NC );

  auto dCompFrac_dCompDensView =
    AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dCompFrac_dCompDens[0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], NC, NC );

  for( localIndex ke = 0; ke < MAX_STENCIL; ++ke )
  {
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      real64 const phaseFlux = 9.283e+03;
      real64 const dPhaseFlux_dP[NUM_ELEMS] = { 0.835e+03, 6.260e+03 };
      real64 const dPhaseFlux_dC[NUM_ELEMS][NC] = { { 0.005e+00, 8.654e+00, 6.126e+02, 9.900e-02 },
        { 4.981e+02, 9.009e+01, 5.747e-02, 8.452e+01 } };

      real64 compFlux[NC]{};
      real64 dCompFlux_dP[MAX_STENCIL][NC]{};
      real64 dCompFlux_dC[MAX_STENCIL][NC][NC]{};

      UHelpers::formPhaseComp( ip,
                               ke,
                               MAX_STENCIL,
                               seri[0],
                               sesri[0],
                               sei[0],
                               phaseCompFracView.toNestedViewConst(),
                               dPhaseCompFrac_dPView.toNestedViewConst(),
                               dPhaseCompFrac_dCView.toNestedViewConst(),
                               dCompFrac_dCompDensView.toNestedViewConst(),
                               phaseFlux,
                               dPhaseFlux_dP,
                               dPhaseFlux_dC,
                               compFlux,
                               dCompFlux_dP,
                               dCompFlux_dC );


      real64 compFluxRef[NC]{};
      real64 dCompFluxRef_dP[MAX_STENCIL][NC]{};
      real64 dCompFluxRef_dC[MAX_STENCIL][NC][NC]{};

      formPhaseCompTest( ip,
                         ke,
                         MAX_STENCIL,
                         phaseCompFrac,
                         dPhaseCompFrac_dPres,
                         dPhaseCompFrac_dComp,
                         dCompFrac_dCompDens,
                         phaseFlux,
                         dPhaseFlux_dP,
                         dPhaseFlux_dC,
                         compFluxRef,
                         dCompFluxRef_dP,
                         dCompFluxRef_dC );



      /* Actual testing (and derivatives) */
      for( localIndex kl = 0; kl < MAX_STENCIL; ++kl )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          EXPECT_DOUBLE_EQ( compFluxRef[ic], compFlux[ic] );
          EXPECT_DOUBLE_EQ( dCompFluxRef_dP[kl][ic], dCompFlux_dP[kl][ic] );
          for( localIndex jc = 0; jc < NC; ++jc )
            EXPECT_DOUBLE_EQ( dCompFluxRef_dC[kl][ic][jc], dCompFlux_dC[kl][ic][jc] );
        }
      }


    }
  }


}

template< localIndex NC, localIndex NP, term T, localIndex MAX_STENCIL, bool FULL >
void testCompositionalUpwindFormPotential( CellElementStencilTPFA const & stencil,
                                           real64 const (*gravCoef)[MAX_STENCIL],
                                           real64 const (*dCompFrac_dCompDens)[MAX_STENCIL][NC][NC],
                                           real64 const (*phaseMassDens)[MAX_STENCIL][1][NP],
                                           real64 const (*dPhaseMassDens_dP)[MAX_STENCIL][1][NP],
                                           real64 const (*dPhaseMassDens_dC)[MAX_STENCIL][1][NP][NC],
                                           real64 const (*dPhaseVolFrac_dP)[MAX_STENCIL][NP],
                                           real64 const (*dPhaseVolFrac_dC)[MAX_STENCIL][NP][NC],
                                           real64 const (*phaseCapPressure)[MAX_STENCIL][1][NP],
                                           real64 const (*dPhaseCapPressure_dPhaseVolFrac)[MAX_STENCIL][1][NP][NP]
                                           )
{
  localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;

  typename CellElementStencilTPFA::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
  typename CellElementStencilTPFA::IndexContainerViewConstType const & sei = stencil.getElementIndices();
  auto const & weights = stencil.getWeights();

  auto gravCoefView = AccessorHelper< FULL >::template makeElementAccessor< 1 >( &( gravCoef[0][0] ),
                                                                                 MAX_STENCIL,
                                                                                 seri[0],
                                                                                 sesri[0],
                                                                                 sei[0] );

  auto phaseMassDensView = AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( phaseMassDens[0][0][0][0] ),
                                                                                      MAX_STENCIL,
                                                                                      seri[0],
                                                                                      sesri[0],
                                                                                      sei[0], 1, NP );

  auto dPhaseMassDens_dPView =
    AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dPhaseMassDens_dP[0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], 1, NP );

  auto dPhaseMassDens_dCView =
    AccessorHelper< FULL >::template makeElementAccessor< 4 >( &( dPhaseMassDens_dC[0][0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], 1, NP, NC );
  auto dCompFrac_dCompDensView =
    AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dCompFrac_dCompDens[0][0][0][0] ),
                                                               MAX_STENCIL,
                                                               seri[0],
                                                               sesri[0],
                                                               sei[0], NC, NC );

  //capillary part
  auto phaseCapPressureView = AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( phaseCapPressure[0][0][0][0] ),
                                                                                         MAX_STENCIL,
                                                                                         seri[0],
                                                                                         sesri[0],
                                                                                         sei[0], 1, NP );

  auto dPhaseCapPressure_dPhaseVolFracView = AccessorHelper< FULL >::template makeElementAccessor< 4 >( &( dPhaseCapPressure_dPhaseVolFrac[0][0][0][0][0] ),
                                                                                                        MAX_STENCIL,
                                                                                                        seri[0],
                                                                                                        sesri[0],
                                                                                                        sei[0], 1, NP, NP );

  auto dPhaseVolFrac_dPView = AccessorHelper< FULL >::template makeElementAccessor< 2 >( &( dPhaseVolFrac_dP[0][0][0] ),
                                                                                         MAX_STENCIL,
                                                                                         seri[0],
                                                                                         sesri[0],
                                                                                         sei[0], NP );

  //NP NC
  auto dPhaseVolFrac_dCView = AccessorHelper< FULL >::template makeElementAccessor< 3 >( &( dPhaseVolFrac_dC[0][0][0][0] ),
                                                                                         MAX_STENCIL,
                                                                                         seri[0],
                                                                                         sesri[0],
                                                                                         sei[0], NP, NC );

  for( localIndex ip = 0; ip < NP; ++ip )
  {

    real64 pot{};
    real64 dPot_dP[NUM_ELEMS]{};
    real64 dPot_dC[NUM_ELEMS][NC]{};
    real64 dProp_dC[NC]{};

    real64 totFlux = 0.0;

    UHelpers::formPotential< NC, T, NUM_ELEMS, MAX_STENCIL >::compute(
      NP,
      ip,
      MAX_STENCIL,
      seri[0],
      sesri[0],
      sei[0],
      weights[0],
      totFlux,
      gravCoefView.toNestedViewConst(),
      dCompFrac_dCompDensView.toNestedViewConst(),
      phaseMassDensView.toNestedViewConst(),
      dPhaseMassDens_dPView.toNestedViewConst(),
      dPhaseMassDens_dCView.toNestedViewConst(),
      dPhaseVolFrac_dPView.toNestedViewConst(),
      dPhaseVolFrac_dCView.toNestedViewConst(),
      phaseCapPressureView.toNestedViewConst(),
      dPhaseCapPressure_dPhaseVolFracView.toNestedViewConst(),
      pot,
      dPot_dP,
      dPot_dC,
      dProp_dC
      );

    real64 potRef{};
    real64 dPotRef_dP[NUM_ELEMS]{};
    real64 dPotRef_dC[NUM_ELEMS][NC]{};
    real64 dPropRef_dC[NC]{};


    formPotTest< NC, NP, T, NUM_ELEMS, MAX_STENCIL >::compute( NP,
                                                               ip,
                                                               MAX_STENCIL,
                                                               weights.data(),
                                                               gravCoef,
                                                               dCompFrac_dCompDens,
                                                               phaseMassDens,
                                                               dPhaseMassDens_dP,
                                                               dPhaseMassDens_dC,
                                                               dPhaseVolFrac_dP,
                                                               dPhaseVolFrac_dC,
                                                               phaseCapPressure,
                                                               dPhaseCapPressure_dPhaseVolFrac,
                                                               potRef,
                                                               dPotRef_dP,
                                                               dPotRef_dC,
                                                               dPropRef_dC );


    EXPECT_DOUBLE_EQ( pot, potRef );
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      EXPECT_DOUBLE_EQ( dPot_dP[i], dPotRef_dP[i] );

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        EXPECT_DOUBLE_EQ( dPot_dC[i][ic], dPotRef_dC[i][ic] );
        EXPECT_DOUBLE_EQ( dProp_dC[ic], dPropRef_dC[ic] );
      }
    }
  }
}


template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL, localIndex NDOF >
void testCompositionalUpwindFillLocalJacobi( real64 const (*compFlux)[NC],
                                             real64 const (*dCompFlux_dP)[MAX_STENCIL][NC],
                                             real64 const (*dCompFlux_dC)[MAX_STENCIL][NC][NC],
                                             localIndex const stencilSize,
                                             real64 const * dt )
{

  //todo type manip to  get from localFlux as arraySlice
  stackArray1d< real64, NUM_ELEMS * NC > localFlux( NUM_ELEMS * NC );
  stackArray2d< real64, NUM_ELEMS * NC * MAX_STENCIL * NDOF > localFluxJacobian( NUM_ELEMS * NC, MAX_STENCIL * NDOF );

  UHelpers::fillLocalJacobi< NC, MAX_STENCIL, NDOF >( ( *compFlux ),
                                                      ( *dCompFlux_dP ),
                                                      ( *dCompFlux_dC ),
                                                      MAX_STENCIL,
                                                      ( *dt ),
                                                      localFlux,
                                                      localFluxJacobian );


  real64 localFluxRef[NUM_ELEMS * NC]{};
  real64 localFluxJacobianRef[NUM_ELEMS * NC][MAX_STENCIL * NDOF]{};

  fillLocalJacobiTest< NC, MAX_STENCIL, NUM_ELEMS, NDOF >( ( *compFlux ),
                                                           ( *dCompFlux_dP ),
                                                           ( *dCompFlux_dC ),
                                                           MAX_STENCIL,
                                                           ( *dt ),
                                                           localFluxRef,
                                                           localFluxJacobianRef );


  //expecting
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    EXPECT_DOUBLE_EQ( localFlux[ic], localFluxRef[ic] );

    for( localIndex ke = 0; ke < stencilSize; ++ke )
    {
      localIndex const localDofIndexPres = ke * NDOF;
      EXPECT_DOUBLE_EQ( localFluxJacobian[ic][localDofIndexPres],
                        localFluxJacobianRef[ic][localDofIndexPres] );
      EXPECT_DOUBLE_EQ( localFluxJacobian[ic + NC][localDofIndexPres],
                        localFluxJacobianRef[ic + NC][localDofIndexPres] );

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        localIndex const localDofIndexComp = localDofIndexPres + jc + 1;
        EXPECT_DOUBLE_EQ( localFluxJacobian[ic][localDofIndexComp],
                          localFluxJacobianRef[ic][localDofIndexComp] );
        EXPECT_DOUBLE_EQ( localFluxJacobian[ic + NC][localDofIndexComp],
                          localFluxJacobianRef[ic + NC][localDofIndexComp] );

      }
    }
  }

}

// ------------------------- gtest
TEST( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_mDens )
{
  localIndex constexpr MAX_STENCIL = 2;
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;

  CellElementStencilTPFA stencil;

  localIndex elemReg[2] = { 0, 1 };
  localIndex elemSubReg[2] = { 0, 0 };
  localIndex elemIndex[2] = { 1, 0 };
  real64 weight[] = { 1e-12, -1e-12 };
  stencil.add( MAX_STENCIL,
               elemReg,
               elemSubReg,
               elemIndex,
               weight,
               0 );


  int constexpr NTEST = 1;

  // we keep these around for easy aggregate initialization
  real64 const densData[NTEST][MAX_STENCIL][1][NP] = {
    { { { 1.562e+3, 2e+3 } },
      { { 2e+0, 4e+0 } } }
  };
  real64 const dDensData_dP[NTEST][MAX_STENCIL][1][NP] = {
    { { { 1e+5, 1e+5 } },
      { { 1e+2, 2e+2 } } }
  };
  real64 const dDensData_dC[NTEST][MAX_STENCIL][1][NP][NC] = {
    {
      { { { 8.244e+01, 9.827e-03, 7.302e+03, 3.439e+03 },
        { 8.178e-01, 2.607e-01, 5.944e-02, 0.225e-02 }
      } },
      { { { 5.313e+02, 3.251e-01, 1.056e-03, 6.110e-02 },
        { 1.537e+00, 2.810e+03, 4.401e+00, 5.271e+03 },
      } }
    }
  };

  real64 const dCompFracData_dCompDens[NTEST][MAX_STENCIL][NC][NC] = {
    {
      {
        { 8.055e+03, 5.767e-03, 1.829e+00, 2.399e-02 },
        { 9.787e-03, 7.127e+01, 5.005e-03, 4.711e-03 },
        { 5.216e+02, 0.967e-02, 8.181e+01, 8.175e+00 },
        { 9.730e+00, 6.490e+02, 8.003e-03, 4.538e-03 }
      },
      {
        { 1.734e-03, 3.909e-01, 8.314e+00, 8.034e-01 },
        { 6.569e-03, 6.280e+03, 2.920e-02, 4.317e-03 },
        { 3.724e+03, 1.981e+03, 4.897e-03, 3.395e+02 },
        { 2.691e-01, 4.228e+03, 5.479e-01, 9.427e+01 }
      },
    }
  };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );
    testCompositionalUpwindDensMult< NC, NP, MAX_STENCIL, true >( stencil,
                                                                  &dCompFracData_dCompDens[i],
                                                                  &densData[i],
                                                                  &dDensData_dP[i],
                                                                  &dDensData_dC[i]
                                                                  );

  }

}

TEST( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_formPhaseComp )
{
  localIndex constexpr MAX_STENCIL = 2;
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;

  CellElementStencilTPFA stencil;

  localIndex elemReg[2] = { 0, 1 };
  localIndex elemSubReg[2] = { 0, 0 };
  localIndex elemIndex[2] = { 1, 0 };
  real64 weight[] = { 1e-12, -1e-12 };
  stencil.add( MAX_STENCIL,
               elemReg,
               elemSubReg,
               elemIndex,
               weight,
               0 );

  int constexpr NTEST = 1;

  // we keep these around for easy aggregate initialization
  real64 const phaseCompFrac[NTEST][MAX_STENCIL][1][NP][NC] = {
    {
      { {
        { 6.663e-02, 5.391e-03, 6.981e+03, 6.665e-02 },
        { 0.326e-02, 5.612e-01, 8.819e+00, 6.692e+03 }
      } },
      { {
        { 1.564e-02, 8.555e-01, 6.448e+00, 3.763e-03 },
        { 5.895e-02, 2.262e-01, 3.846e+01, 5.830e-02 }
      } }
    }
  };
  real64 const dPhaseCompFrac_dPres[NTEST][MAX_STENCIL][1][NP][NC] = {
    {
      { { { 6.022e+00, 3.868e-01, 9.160e+00, 0.012e+02 },
        { 3.225e-02, 7.847e+02, 4.714e+00, 0.358e-02 } } },
      { { { 3.411e-02, 6.074e+03, 1.917e-02, 7.384e+02 },
        { 1.887e+01, 2.875e+00, 0.911e-01, 5.762e+01 } } }
    }
  };

  real64 const dPhaseCompFrac_dComp[NTEST][MAX_STENCIL][1][NP][NC][NC] = {
    {//MAX_STENCIL
      {//1
        {//NP
          {
            { 6.476e-02, 6.790e+01, 6.358e-02, 9.452e-03 },
            { 6.073e+02, 4.501e-01, 4.587e+01, 6.619e-01 },
            { 8.419e+01, 8.329e+00, 2.564e+03, 6.135e-02 },
            { 3.181e+00, 1.192e+01, 9.398e+00, 6.456e+01 }
          },
          {
            { 5.439e-02, 7.210e-03, 5.225e-03, 9.937e-03 },
            { 4.046e+01, 4.484e+02, 3.658e+03, 7.635e+03 },
            { 1.920e+00, 1.389e+00, 6.963e+03, 0.938e+00 },
            { 3.935e-01, 6.714e-02, 7.413e+01, 5.201e-02 }
          },
        }
      },
      {
        {//NP
          {
            { 0.445e+01, 7.549e-01, 2.428e+02, 4.424e-01 },
            { 6.834e-01, 7.040e-01, 4.423e-02, 0.196e-02 },
            { 8.217e+02, 4.299e-01, 8.878e+02, 3.912e+02 },
            { 3.774e-01, 2.160e+01, 7.904e+00, 9.493e+02 }
          },
          {
            { 7.689e+00, 1.673e+03, 8.620e+01, 9.899e-02 },
            { 1.999e+02, 4.070e-01, 7.487e+00, 8.256e-03 },
            { 1.117e-02, 1.363e+00, 6.787e-02, 4.952e-03 },
            { 8.507e+01, 5.606e+02, 9.296e+03, 6.967e+03 }
          },
        }
      }
    }
  };

  real64 const dCompFracData_dCompDens[NTEST][MAX_STENCIL][NC][NC] = {
    {
      {
        { 8.055e+03, 5.767e-03, 1.829e+00, 2.399e-02 },
        { 9.787e-03, 7.127e+01, 5.005e-03, 4.711e-03 },
        { 5.216e+02, 0.967e-02, 8.181e+01, 8.175e+00 },
        { 9.730e+00, 6.490e+02, 8.003e-03, 4.538e-03 }
      },
      {
        { 1.734e-03, 3.909e-01, 8.314e+00, 8.034e-01 },
        { 6.569e-03, 6.280e+03, 2.920e-02, 4.317e-03 },
        { 3.724e+03, 1.981e+03, 4.897e-03, 3.395e+02 },
        { 2.691e-01, 4.228e+03, 5.479e-01, 9.427e+01 }
      },
    }
  };


  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );
    testCompositionalUpwindFormPhaseComp< NC, NP, MAX_STENCIL, true >( stencil,
                                                                       &dCompFracData_dCompDens[i],
                                                                       &phaseCompFrac[i],
                                                                       &dPhaseCompFrac_dPres[i],
                                                                       &dPhaseCompFrac_dComp[i] );

  }

}

template< term T >
void testPotential()
{
  localIndex constexpr MAX_STENCIL = 2;
  localIndex constexpr NP = 2;
  localIndex constexpr NC = 4;

  CellElementStencilTPFA stencil;

  localIndex elemReg[2] = { 0, 1 };
  localIndex elemSubReg[2] = { 0, 0 };
  localIndex elemIndex[2] = { 1, 0 };
  real64 weight[] = { 1e-12, -1e-12 };
  stencil.add( MAX_STENCIL,
               elemReg,
               elemSubReg,
               elemIndex,
               weight,
               0 );

  int constexpr NTEST = 1;
  //here goes synthetic data
  real64 const gravCoef[NTEST][MAX_STENCIL] = { -9.81, -0.981 };

  real64 const phaseMassDensData[NTEST][MAX_STENCIL][1][NP] = {
    { { { 1.562e+3, 2e+3 } }, { { 2e+0, 4e+0 } } }
  };

  real64 const dPhaseMassDensData_dP[NTEST][MAX_STENCIL][1][NP] = {
    { { { 1e+5, 1e+5 } },
      { { 1e+2, 2e+2 } } }
  };
  real64 const dPhaseMassDensData_dC[NTEST][MAX_STENCIL][1][NP][NC] = {
    {
      { { { 1e+2, 2e+2, 3e+3, 4e+1 }, { 2e+2, 3e+3, 4e+1, 5e+2 } } },
      { { { 5e+1, 6e+1, 7e+1, 8e+1 }, { 7e+1, 8e+1, 9e+1, 1e+0 } } }
    }
  };

  real64 const dCompFracData_dCompDens[NTEST][MAX_STENCIL][NC][NC] = {
    {
      {
        { 8.055e+03, 5.767e-03, 1.829e+00, 2.399e-02 },
        { 9.787e-03, 7.127e+01, 5.005e-03, 4.711e-03 },
        { 5.216e+02, 0.967e-02, 8.181e+01, 8.175e+00 },
        { 9.730e+00, 6.490e+02, 8.003e-03, 4.538e-03 }
      },
      {
        { 1.734e-03, 3.909e-01, 8.314e+00, 8.034e-01 },
        { 6.569e-03, 6.280e+03, 2.920e-02, 4.317e-03 },
        { 3.724e+03, 1.981e+03, 4.897e-03, 3.395e+02 },
        { 2.691e-01, 4.228e+03, 5.479e-01, 9.427e+01 }
      },
    }
  };

  real64 const dPhaseVolFracData_dP[NTEST][MAX_STENCIL][NP] = {
    {
      { 8.055e+03, 5.767e-03 },
      { 1.981e+03, 4.897e-03 }
    }
  };

  real64 const dPhaseVolFracData_dC[NTEST][MAX_STENCIL][NP][NC] = {
    {
      { { 8.055e+03, 5.767e-03, 2.920e-02, 4.317e-03 },
        { 3.724e+03, 1.981e+03, 4.897e-03, 3.395e+02 } },
      { { 2.691e-01, 4.228e+03, 5.479e-01, 9.427e+01 },
        { 5.216e+02, 0.967e-02, 8.181e+01, 8.175e+00 } }
    }
  };

  real64 const phaseCapPressure[NTEST][MAX_STENCIL][1][NP] = {
    {
      {{ 9.787e-03, 7.127e+01 }},
      {{  1.981e+03, 4.897e-03 }}
    }
  };

  real64 const dPhaseCapPressureData_dPhaseVolFrac[NTEST][MAX_STENCIL][1][NP][NP] = {
    {
      {
        {
          { 8.055e+03, 5.767e-03 },
          { 9.787e-03, 7.127e+01 }
        }
      },
      {
        {
          { 1.734e-03, 3.909e-01 },
          { 6.569e-03, 6.280e+03 }
        }
      }
    }
  };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );
    testCompositionalUpwindFormPotential< NC, NP, T, MAX_STENCIL, true >( stencil,
                                                                          &gravCoef[i],
                                                                          &dCompFracData_dCompDens[i],
                                                                          &phaseMassDensData[i],
                                                                          &dPhaseMassDensData_dP[i],
                                                                          &dPhaseMassDensData_dC[i],
                                                                          &dPhaseVolFracData_dP[i],
                                                                          &dPhaseVolFracData_dC[i],
                                                                          &phaseCapPressure[i],
                                                                          &dPhaseCapPressureData_dPhaseVolFrac[i] );

  }
}

TEST( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_formPotGrav )
{
  testPotential< term::Gravity >();
}
TEST( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_formPotCap )
{
  testPotential< term::Capillary >();
}


TEST( CompositionalMultiphaseFlowUpwindHelperKernelsTest, test_fillLocalJacobi )
{
  localIndex constexpr NC = 4;
  localIndex constexpr NDOF = NC + 1;
  localIndex constexpr NUM_ELEMS = CellElementStencilTPFA::NUM_POINT_IN_FLUX;
  localIndex constexpr MAX_STENCIL = CellElementStencilTPFA::MAX_STENCIL_SIZE;

  localIndex constexpr NTEST = 1;

  //synthetic data compFluex (and derivatives)  and dt
  real64 const compFlux[NTEST][NC] = { 5.085e+01, 5.108e-01, 8.176e+02, 7.948e+00 };

  real64 const dCompFlux_dP[NTEST][MAX_STENCIL][NC] = {
    {
      { 3.507e+01, 9.390e+01, 8.759e-02, 5.502e-01 },
      { 4.709e-02, 2.305e-02, 8.443e-02, 1.948e+00 }
    }
  };

  real64 const dCompFlux_dC[NTEST][MAX_STENCIL][NC][NC] = {
    {
      {
        { 3.111e+03, 9.234e+03, 4.302e+00, 1.848e-03 },
        { 2.581e+01, 4.087e+01, 5.949e-02, 2.622e-03 },
        { 2.967e-03, 3.188e-02, 4.242e+02, 5.079e-03 },
        { 9.289e-02, 7.303e+00, 4.886e+03, 5.785e+00 },
      },
      {
        { 5.211e+01, 2.316e-01, 4.889e-01, 6.241e+03 },
        { 0.377e-03, 8.852e-02, 9.133e-01, 7.962e+01 },
        { 1.366e+00, 7.212e+02, 1.068e+02, 6.538e+03 },
        { 8.909e-03, 3.342e+02, 6.987e+00, 1.978e+00 }
      }
    }
  };

  real64 const dt[NTEST] = { 1e4 };

  for( int i = 0; i < NTEST; ++i )
  {
    SCOPED_TRACE( "Input # " + std::to_string( i ) );
    testCompositionalUpwindFillLocalJacobi< NC, NUM_ELEMS, MAX_STENCIL, NDOF >( &compFlux[i],
                                                                                &dCompFlux_dP[i],
                                                                                &dCompFlux_dC[i],
                                                                                MAX_STENCIL,
                                                                                &dt[i] );
  }

}


int main( int argc,
          char * * argv )
{
  ::testing::InitGoogleTest( &argc, argv );
  g_commandLineOptions = *geosx::basicSetup( argc, argv );
  int const result = RUN_ALL_TESTS();
  geosx::basicCleanup();
  return result;
}
