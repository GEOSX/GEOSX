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
 * @file CompositionalMultiphaseFlowKernels.cpp
 */

#include "CompositionalMultiphaseFlowKernels.hpp"

#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "finiteVolume/FaceElementStencil.hpp"

namespace geosx
{
namespace CompositionalMultiphaseFlowKernels
{
/******************************** ComponentFractionKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
ComponentFractionKernel::Compute(
  arraySlice1d< real64 const > const compDens,
  arraySlice1d< real64 const > const dCompDens,
  arraySlice1d< real64 > const compFrac,
  arraySlice2d< real64 > const dCompFrac_dCompDens )
{
  real64 totalDensity = 0.0;

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  real64 const totalDensityInv = 1.0 / totalDensity;

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    compFrac[ic] = ( compDens[ic] + dCompDens[ic] ) * totalDensityInv;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
    }
    dCompFrac_dCompDens[ic][ic] += totalDensityInv;
  }
}

template< localIndex NC >
void
ComponentFractionKernel::Launch( localIndex const size,
                                 arrayView2d< real64 const > const & compDens,
                                 arrayView2d< real64 const > const & dCompDens,
                                 arrayView2d< real64 > const & compFrac,
                                 arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex const a ) {
    Compute< NC >( compDens[a], dCompDens[a], compFrac[a], dCompFrac_dCompDens[a] );
  } );
}

template< localIndex NC >
void
ComponentFractionKernel::Launch(
  SortedArrayView< localIndex const > const & targetSet,
  arrayView2d< real64 const > const & compDens,
  arrayView2d< real64 const > const & dCompDens,
  arrayView2d< real64 > const & compFrac,
  arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >(
    targetSet.size(),
    [=] GEOSX_HOST_DEVICE( localIndex const i ) {
      localIndex const a = targetSet[i];
      Compute< NC >( compDens[a], dCompDens[a], compFrac[a], dCompFrac_dCompDens[a] );
    } );
}

#define INST_ComponentFractionKernel( NC )                 \
  template void ComponentFractionKernel::Launch< NC >(     \
    localIndex const size,                                 \
    arrayView2d< real64 const > const & compDens,          \
    arrayView2d< real64 const > const & dCompDens,         \
    arrayView2d< real64 > const & compFrac,                \
    arrayView3d< real64 > const & dCompFrac_dCompDens );   \
  template void ComponentFractionKernel::Launch< NC >(     \
    SortedArrayView< localIndex const > const & targetSet, \
    arrayView2d< real64 const > const & compDens,          \
    arrayView2d< real64 const > const & dCompDens,         \
    arrayView2d< real64 > const & compFrac,                \
    arrayView3d< real64 > const & dCompFrac_dCompDens )

INST_ComponentFractionKernel( 1 );
INST_ComponentFractionKernel( 2 );
INST_ComponentFractionKernel( 3 );
INST_ComponentFractionKernel( 4 );
INST_ComponentFractionKernel( 5 );

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
PhaseVolumeFractionKernel::Compute(
  arraySlice1d< real64 const > const & compDens,
  arraySlice1d< real64 const > const & dCompDens,
  arraySlice2d< real64 const > const & dCompFrac_dCompDens,
  arraySlice1d< real64 const > const & phaseDens,
  arraySlice1d< real64 const > const & dPhaseDens_dPres,
  arraySlice2d< real64 const > const & dPhaseDens_dComp,
  arraySlice1d< real64 const > const & phaseFrac,
  arraySlice1d< real64 const > const & dPhaseFrac_dPres,
  arraySlice2d< real64 const > const & dPhaseFrac_dComp,
  arraySlice1d< real64 > const & phaseVolFrac,
  arraySlice1d< real64 > const & dPhaseVolFrac_dPres,
  arraySlice2d< real64 > const & dPhaseVolFrac_dComp )
{
  real64 work[NC];

  // compute total density from component partial densities
  real64 totalDensity = 0.0;
  real64 const dTotalDens_dCompDens = 1.0;
  for( localIndex ic = 0; ic < NC; ++ic )
  {
    totalDensity += compDens[ic] + dCompDens[ic];
  }

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // Expression for volume fractions: S_p = (nu_p / rho_p) * rho_t
    real64 const phaseDensInv = 1.0 / phaseDens[ip];

    // compute saturation and derivatives except multiplying by the total density
    phaseVolFrac[ip] = phaseFrac[ip] * phaseDensInv;

    dPhaseVolFrac_dPres[ip] =
      ( dPhaseFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseDens_dPres[ip] ) *
      phaseDensInv;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseVolFrac_dComp[ip][jc] =
        ( dPhaseFrac_dComp[ip][jc] - phaseVolFrac[ip] * dPhaseDens_dComp[ip][jc] ) *
        phaseDensInv;
    }

    // apply chain rule to convert derivatives from global component fractions to densities
    applyChainRuleInPlace( NC, dCompFrac_dCompDens, dPhaseVolFrac_dComp[ip], work );

    // now finalize the computation by multiplying by total density
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseVolFrac_dComp[ip][jc] *= totalDensity;
      dPhaseVolFrac_dComp[ip][jc] += phaseVolFrac[ip] * dTotalDens_dCompDens;
    }

    phaseVolFrac[ip] *= totalDensity;
    dPhaseVolFrac_dPres[ip] *= totalDensity;
  }
}

template< localIndex NC, localIndex NP >
void
PhaseVolumeFractionKernel::Launch(
  localIndex const size,
  arrayView2d< real64 const > const & compDens,
  arrayView2d< real64 const > const & dCompDens,
  arrayView3d< real64 const > const & dCompFrac_dCompDens,
  arrayView3d< real64 const > const & phaseDens,
  arrayView3d< real64 const > const & dPhaseDens_dPres,
  arrayView4d< real64 const > const & dPhaseDens_dComp,
  arrayView3d< real64 const > const & phaseFrac,
  arrayView3d< real64 const > const & dPhaseFrac_dPres,
  arrayView4d< real64 const > const & dPhaseFrac_dComp,
  arrayView2d< real64 > const & phaseVolFrac,
  arrayView2d< real64 > const & dPhaseVolFrac_dPres,
  arrayView3d< real64 > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex const a ) {
    Compute< NC, NP >( compDens[a],
                       dCompDens[a],
                       dCompFrac_dCompDens[a],
                       phaseDens[a][0],
                       dPhaseDens_dPres[a][0],
                       dPhaseDens_dComp[a][0],
                       phaseFrac[a][0],
                       dPhaseFrac_dPres[a][0],
                       dPhaseFrac_dComp[a][0],
                       phaseVolFrac[a],
                       dPhaseVolFrac_dPres[a],
                       dPhaseVolFrac_dComp[a] );
  } );
}

template< localIndex NC, localIndex NP >
void
PhaseVolumeFractionKernel::Launch(
  SortedArrayView< localIndex const > const & targetSet,
  arrayView2d< real64 const > const & compDens,
  arrayView2d< real64 const > const & dCompDens,
  arrayView3d< real64 const > const & dCompFrac_dCompDens,
  arrayView3d< real64 const > const & phaseDens,
  arrayView3d< real64 const > const & dPhaseDens_dPres,
  arrayView4d< real64 const > const & dPhaseDens_dComp,
  arrayView3d< real64 const > const & phaseFrac,
  arrayView3d< real64 const > const & dPhaseFrac_dPres,
  arrayView4d< real64 const > const & dPhaseFrac_dComp,
  arrayView2d< real64 > const & phaseVolFrac,
  arrayView2d< real64 > const & dPhaseVolFrac_dPres,
  arrayView3d< real64 > const & dPhaseVolFrac_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(),
                                    [=] GEOSX_HOST_DEVICE( localIndex const i ) {
                                      localIndex const a = targetSet[i];
                                      Compute< NC, NP >( compDens[a],
                                                         dCompDens[a],
                                                         dCompFrac_dCompDens[a],
                                                         phaseDens[a][0],
                                                         dPhaseDens_dPres[a][0],
                                                         dPhaseDens_dComp[a][0],
                                                         phaseFrac[a][0],
                                                         dPhaseFrac_dPres[a][0],
                                                         dPhaseFrac_dComp[a][0],
                                                         phaseVolFrac[a],
                                                         dPhaseVolFrac_dPres[a],
                                                         dPhaseVolFrac_dComp[a] );
                                    } );
}

#define INST_PhaseVolumeFractionKernel( NC, NP )             \
  template void PhaseVolumeFractionKernel::Launch< NC, NP >( \
    localIndex const size,                                   \
    arrayView2d< real64 const > const & compDens,            \
    arrayView2d< real64 const > const & dCompDens,           \
    arrayView3d< real64 const > const & dCompFrac_dCompDens, \
    arrayView3d< real64 const > const & phaseDens,           \
    arrayView3d< real64 const > const & dPhaseDens_dPres,    \
    arrayView4d< real64 const > const & dPhaseDens_dComp,    \
    arrayView3d< real64 const > const & phaseFrac,           \
    arrayView3d< real64 const > const & dPhaseFrac_dPres,    \
    arrayView4d< real64 const > const & dPhaseFrac_dComp,    \
    arrayView2d< real64 > const & phaseVolFrac,              \
    arrayView2d< real64 > const & dPhaseVolFrac_dPres,       \
    arrayView3d< real64 > const & dPhaseVolFrac_dComp );     \
  template void PhaseVolumeFractionKernel::Launch< NC, NP >( \
    SortedArrayView< localIndex const > const & targetSet,   \
    arrayView2d< real64 const > const & compDens,            \
    arrayView2d< real64 const > const & dCompDens,           \
    arrayView3d< real64 const > const & dCompFrac_dCompDens, \
    arrayView3d< real64 const > const & phaseDens,           \
    arrayView3d< real64 const > const & dPhaseDens_dPres,    \
    arrayView4d< real64 const > const & dPhaseDens_dComp,    \
    arrayView3d< real64 const > const & phaseFrac,           \
    arrayView3d< real64 const > const & dPhaseFrac_dPres,    \
    arrayView4d< real64 const > const & dPhaseFrac_dComp,    \
    arrayView2d< real64 > const & phaseVolFrac,              \
    arrayView2d< real64 > const & dPhaseVolFrac_dPres,       \
    arrayView3d< real64 > const & dPhaseVolFrac_dComp )

INST_PhaseVolumeFractionKernel( 1, 1 );
INST_PhaseVolumeFractionKernel( 2, 1 );
INST_PhaseVolumeFractionKernel( 3, 1 );
INST_PhaseVolumeFractionKernel( 4, 1 );
INST_PhaseVolumeFractionKernel( 5, 1 );

INST_PhaseVolumeFractionKernel( 1, 2 );
INST_PhaseVolumeFractionKernel( 2, 2 );
INST_PhaseVolumeFractionKernel( 3, 2 );
INST_PhaseVolumeFractionKernel( 4, 2 );
INST_PhaseVolumeFractionKernel( 5, 2 );

INST_PhaseVolumeFractionKernel( 1, 3 );
INST_PhaseVolumeFractionKernel( 2, 3 );
INST_PhaseVolumeFractionKernel( 3, 3 );
INST_PhaseVolumeFractionKernel( 4, 3 );
INST_PhaseVolumeFractionKernel( 5, 3 );

#undef INST_PhaseVolumeFractionKernel

/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
PhaseMobilityKernel::Compute(
  arraySlice2d< real64 const > const & dCompFrac_dCompDens,
  arraySlice1d< real64 const > const & phaseDens,
  arraySlice1d< real64 const > const & dPhaseDens_dPres,
  arraySlice2d< real64 const > const & dPhaseDens_dComp,
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
  real64 dDens_dC[NC];
  real64 dVisc_dC[NC];

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const density = phaseDens[ip];
    real64 const dDens_dP = dPhaseDens_dPres[ip];
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dDens_dC );

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

    real64 const mobility = relPerm * density / viscosity;

    phaseMob[ip] = mobility;
    dPhaseMob_dPres[ip] = dRelPerm_dP * density / viscosity +
      mobility * ( dDens_dP / density - dVisc_dP / viscosity );

    // compositional derivatives
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseMob_dComp[ip][jc] = dRelPerm_dC[jc] * density / viscosity +
        mobility * ( dDens_dC[jc] / density - dVisc_dC[jc] / viscosity );
    }
  }
}

template< localIndex NC, localIndex NP >
void
PhaseMobilityKernel::Launch(
  localIndex const size,
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
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex const a ) {
    Compute< NC, NP >( dCompFrac_dCompDens[a],
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
void
PhaseMobilityKernel::Launch(
  SortedArrayView< localIndex const > const & targetSet,
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
  forAll< parallelDevicePolicy<> >( targetSet.size(),
                                    [=] GEOSX_HOST_DEVICE( localIndex const i ) {
                                      localIndex const a = targetSet[i];
                                      Compute< NC, NP >(
                                        dCompFrac_dCompDens[a],
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

#define INST_PhaseMobilityKernel( NC, NP )                           \
  template void PhaseMobilityKernel::Launch< NC, NP >(               \
    localIndex const size,                                           \
    arrayView3d< real64 const > const & dCompFrac_dCompDens,         \
    arrayView3d< real64 const > const & phaseDens,                   \
    arrayView3d< real64 const > const & dPhaseDens_dPres,            \
    arrayView4d< real64 const > const & dPhaseDens_dComp,            \
    arrayView3d< real64 const > const & phaseVisc,                   \
    arrayView3d< real64 const > const & dPhaseVisc_dPres,            \
    arrayView4d< real64 const > const & dPhaseVisc_dComp,            \
    arrayView3d< real64 const > const & phaseRelPerm,                \
    arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres,         \
    arrayView3d< real64 const > const & dPhaseVolFrac_dComp,         \
    arrayView2d< real64 > const & phaseMob,                          \
    arrayView2d< real64 > const & dPhaseMob_dPres,                   \
    arrayView3d< real64 > const & dPhaseMob_dComp );                 \
  template void PhaseMobilityKernel::Launch< NC, NP >(               \
    SortedArrayView< localIndex const > const & targetSet,           \
    arrayView3d< real64 const > const & dCompFrac_dCompDens,         \
    arrayView3d< real64 const > const & phaseDens,                   \
    arrayView3d< real64 const > const & dPhaseDens_dPres,            \
    arrayView4d< real64 const > const & dPhaseDens_dComp,            \
    arrayView3d< real64 const > const & phaseVisc,                   \
    arrayView3d< real64 const > const & dPhaseVisc_dPres,            \
    arrayView4d< real64 const > const & dPhaseVisc_dComp,            \
    arrayView3d< real64 const > const & phaseRelPerm,                \
    arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac, \
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres,         \
    arrayView3d< real64 const > const & dPhaseVolFrac_dComp,         \
    arrayView2d< real64 > const & phaseMob,                          \
    arrayView2d< real64 > const & dPhaseMob_dPres,                   \
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

/******************************** AccumulationKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
AccumulationKernel::Compute(
  localIndex const numPhases,
  real64 const & volume,
  real64 const & porosityOld,
  real64 const & porosityRef,
  real64 const & pvMult,
  real64 const & dPvMult_dPres,
  arraySlice2d< real64 const > const & dCompFrac_dCompDens,
  arraySlice1d< real64 const > const & phaseVolFracOld,
  arraySlice1d< real64 const > const & phaseVolFrac,
  arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
  arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
  arraySlice1d< real64 const > const & phaseDensOld,
  arraySlice1d< real64 const > const & phaseDens,
  arraySlice1d< real64 const > const & dPhaseDens_dPres,
  arraySlice2d< real64 const > const & dPhaseDens_dComp,
  arraySlice2d< real64 const > const & phaseCompFracOld,
  arraySlice2d< real64 const > const & phaseCompFrac,
  arraySlice2d< real64 const > const & dPhaseCompFrac_dPres,
  arraySlice3d< real64 const > const & dPhaseCompFrac_dComp,
  real64 ( &localAccum )[NC],
  real64 ( &localAccumJacobian )[NC][NC + 1] )
{
  localIndex constexpr NDOF = NC + 1;
  localIndex const NP = numPhases;

  // temporary work arrays
  real64 dPhaseAmount_dC[NC];
  real64 dPhaseCompFrac_dC[NC];

  // reset the local values
  for( localIndex i = 0; i < NC; ++i )
  {
    localAccum[i] = 0.0;
    for( localIndex j = 0; j < NDOF; ++j )
    {
      localAccumJacobian[i][j] = 0.0;
    }
  }

  // compute fluid-independent (pore volume) part
  real64 const volNew = volume;
  real64 const volOld = volume;
  real64 const dVol_dP = 0.0;  // used in poroelastic solver

  real64 const poroNew = porosityRef * pvMult;
  real64 const poroOld = porosityOld;
  real64 const dPoro_dP = porosityRef * dPvMult_dPres;

  real64 const poreVolNew = volNew * poroNew;
  real64 const poreVolOld = volOld * poroOld;
  real64 const dPoreVol_dP = dVol_dP * poroNew + volNew * dPoro_dP;

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
    real64 const phaseAmountOld =
      poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

    real64 const dPhaseAmount_dP =
      dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip] +
      poreVolNew *
        ( dPhaseVolFrac_dPres[ip] * phaseDens[ip] +
          phaseVolFrac[ip] * dPhaseDens_dPres[ip] );

    // assemble density dependence
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip] +
        phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
      dPhaseAmount_dC[jc] *= poreVolNew;
    }

    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
      real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

      real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic] +
        phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

      localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
      localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

      // jc - index of component w.r.t. whose compositional var the derivative is being taken
      // (i.e. col number in local matrix)

      // assemble phase composition dependence
      applyChainRule( NC,
                      dCompFrac_dCompDens,
                      dPhaseCompFrac_dComp[ip][ic],
                      dPhaseCompFrac_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseCompAmount_dC =
          dPhaseCompFrac_dC[jc] * phaseAmountNew +
          phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
        localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
      }
    }
  }
}

template< localIndex NC >
void
AccumulationKernel::Launch(
  localIndex const numPhases,
  localIndex const size,
  globalIndex const rankOffset,
  arrayView1d< globalIndex const > const & dofNumber,
  arrayView1d< integer const > const & elemGhostRank,
  arrayView1d< real64 const > const & volume,
  arrayView1d< real64 const > const & porosityOld,
  arrayView1d< real64 const > const & porosityRef,
  arrayView2d< real64 const > const & pvMult,
  arrayView2d< real64 const > const & dPvMult_dPres,
  arrayView3d< real64 const > const & dCompFrac_dCompDens,
  arrayView2d< real64 const > const & phaseVolFracOld,
  arrayView2d< real64 const > const & phaseVolFrac,
  arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
  arrayView2d< real64 const > const & phaseDensOld,
  arrayView3d< real64 const > const & phaseDens,
  arrayView3d< real64 const > const & dPhaseDens_dPres,
  arrayView4d< real64 const > const & dPhaseDens_dComp,
  arrayView3d< real64 const > const & phaseCompFracOld,
  arrayView4d< real64 const > const & phaseCompFrac,
  arrayView4d< real64 const > const & dPhaseCompFrac_dPres,
  arrayView5d< real64 const > const & dPhaseCompFrac_dComp,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
    if( elemGhostRank[ei] >= 0 ) return;

    localIndex constexpr NDOF = NC + 1;

    real64 localAccum[NC];
    real64 localAccumJacobian[NC][NDOF];

    Compute< NC >( numPhases,
                   volume[ei],
                   porosityOld[ei],
                   porosityRef[ei],
                   pvMult[ei][0],
                   dPvMult_dPres[ei][0],
                   dCompFrac_dCompDens[ei],
                   phaseVolFracOld[ei],
                   phaseVolFrac[ei],
                   dPhaseVolFrac_dPres[ei],
                   dPhaseVolFrac_dCompDens[ei],
                   phaseDensOld[ei],
                   phaseDens[ei][0],
                   dPhaseDens_dPres[ei][0],
                   dPhaseDens_dComp[ei][0],
                   phaseCompFracOld[ei],
                   phaseCompFrac[ei][0],
                   dPhaseCompFrac_dPres[ei][0],
                   dPhaseCompFrac_dComp[ei][0],
                   localAccum,
                   localAccumJacobian );

    // set DOF indices for this block
    localIndex const localRow = dofNumber[ei] - rankOffset;
    globalIndex dofIndices[NDOF];
    for( localIndex idof = 0; idof < NDOF; ++idof )
    {
      dofIndices[idof] = dofNumber[ei] + idof;
    }

    // TODO: apply equation/variable change transformation(s)

    // add contribution to residual and jacobian
    for( localIndex i = 0; i < NC; ++i )
    {
      localRhs[localRow + i] += localAccum[i];
      localMatrix.addToRow< serialAtomic >( localRow + i,
                                            dofIndices,
                                            localAccumJacobian[i],
                                            NDOF );
    }
  } );
}

#define INST_AccumulationKernel( NC )                               \
  template void AccumulationKernel::Launch< NC >(                   \
    localIndex const numPhases,                                     \
    localIndex const size,                                          \
    globalIndex const rankOffset,                                   \
    arrayView1d< globalIndex const > const & dofNumber,             \
    arrayView1d< integer const > const & elemGhostRank,             \
    arrayView1d< real64 const > const & volume,                     \
    arrayView1d< real64 const > const & porosityOld,                \
    arrayView1d< real64 const > const & porosityRef,                \
    arrayView2d< real64 const > const & pvMult,                     \
    arrayView2d< real64 const > const & dPvMult_dPres,              \
    arrayView3d< real64 const > const & dCompFrac_dCompDens,        \
    arrayView2d< real64 const > const & phaseVolFracOld,            \
    arrayView2d< real64 const > const & phaseVolFrac,               \
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres,        \
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,    \
    arrayView2d< real64 const > const & phaseDensOld,               \
    arrayView3d< real64 const > const & phaseDens,                  \
    arrayView3d< real64 const > const & dPhaseDens_dPres,           \
    arrayView4d< real64 const > const & dPhaseDens_dComp,           \
    arrayView3d< real64 const > const & phaseCompFracOld,           \
    arrayView4d< real64 const > const & phaseCompFrac,              \
    arrayView4d< real64 const > const & dPhaseCompFrac_dPres,       \
    arrayView5d< real64 const > const & dPhaseCompFrac_dComp,       \
    CRSMatrixView< real64, globalIndex const > const & localMatrix, \
    arrayView1d< real64 > const & localRhs )

INST_AccumulationKernel( 1 );
INST_AccumulationKernel( 2 );
INST_AccumulationKernel( 3 );
INST_AccumulationKernel( 4 );
INST_AccumulationKernel( 5 );

#undef INST_AccumulationKernel

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
FluxKernel::Compute(
  localIndex const numPhases,
  localIndex const stencilSize,
  arraySlice1d< localIndex const > const & seri,
  arraySlice1d< localIndex const > const & sesri,
  arraySlice1d< localIndex const > const & sei,
  arraySlice1d< real64 const > const & stencilWeights,
  ElementView< arrayView1d< real64 const > > const & pres,
  ElementView< arrayView1d< real64 const > > const & dPres,
  ElementView< arrayView1d< real64 const > > const & gravCoef,
  ElementView< arrayView2d< real64 const > > const & phaseMob,
  ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
  ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp,
  ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
  ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
  ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
  ElementView< arrayView3d< real64 const > > const & phaseDens,
  ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
  ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp,
  ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
  ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
  ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
  ElementView< arrayView3d< real64 const > > const & phaseCapPressure,
  ElementView< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
  integer const capPressureFlag,
  real64 const dt,
  arraySlice1d< real64 > const & localFlux,
  arraySlice2d< real64 > const & localFluxJacobian )
{
  localIndex constexpr NDOF = NC + 1;
  localIndex const NP = numPhases;

  real64 compFlux[NC] {};
  real64 dCompFlux_dP[MAX_STENCIL][NC] {};
  real64 dCompFlux_dC[MAX_STENCIL][NC][NC] {};

  // loop over phases, compute and upwind phase flux and sum contributions to each component's flux
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // clear working arrays
    real64 densMean {};
    real64 dDensMean_dP[NUM_ELEMS] {};
    real64 dDensMean_dC[NUM_ELEMS][NC] {};

    // create local work arrays
    real64 phaseFlux {};
    real64 dPhaseFlux_dP[MAX_STENCIL] {};
    real64 dPhaseFlux_dC[MAX_STENCIL][NC] {};

    real64 presGrad {};
    real64 dPresGrad_dP[MAX_STENCIL] {};
    real64 dPresGrad_dC[MAX_STENCIL][NC] {};

    real64 gravHead {};
    real64 dGravHead_dP[NUM_ELEMS] {};
    real64 dGravHead_dC[NUM_ELEMS][NC] {};

    real64 dCapPressure_dC[NC] {};

    // Working array
    real64 dProp_dC[NC] {};

    // calculate quantities on primary connected cells
    for( localIndex i = 0; i < NUM_ELEMS; ++i )
    {
      localIndex const er = seri[i];
      localIndex const esr = sesri[i];
      localIndex const ei = sei[i];

      // density
      real64 const density = phaseDens[er][esr][ei][0][ip];
      real64 const dDens_dP = dPhaseDens_dPres[er][esr][ei][0][ip];

      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseDens_dComp[er][esr][ei][0][ip],
                      dProp_dC );

      // average density and derivatives
      densMean += 0.5 * density;
      dDensMean_dP[i] = 0.5 * dDens_dP;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDensMean_dC[i][jc] = 0.5 * dProp_dC[jc];
      }
    }

    //***** calculation of flux *****

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

        for( localIndex jp = 0; jp < NP; ++jp )
        {
          real64 const dCapPressure_dS =
            dPhaseCapPressure_dPhaseVolFrac[er][esr][ei][0][ip][jp];
          dCapPressure_dP +=
            dCapPressure_dS * dPhaseVolFrac_dPres[er][esr][ei][jp];

          for( localIndex jc = 0; jc < NC; ++jc )
          {
            dCapPressure_dC[jc] +=
              dCapPressure_dS * dPhaseVolFrac_dComp[er][esr][ei][jp][jc];
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

    // *** upwinding ***

    // use PPU currently; advanced stuff like IHU would go here
    // TODO isolate into a kernel?

    // compute phase potential gradient
    real64 const potGrad = presGrad - gravHead;

    // choose upstream cell
    localIndex const k_up = ( potGrad >= 0 ) ? 0 : 1;

    localIndex er_up = seri[k_up];
    localIndex esr_up = sesri[k_up];
    localIndex ei_up = sei[k_up];

    real64 const mobility = phaseMob[er_up][esr_up][ei_up][ip];

    // skip the phase flux if phase not present or immobile upstream
    if( std::fabs( mobility ) < 1e-20 )  // TODO better constant
    {
      continue;
    }

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
    arraySlice1d< real64 const > dPhaseMob_dCompSub =
      dPhaseMob_dComp[er_up][esr_up][ei_up][ip];

    // add contribution from upstream cell mobility derivatives
    dPhaseFlux_dP[k_up] += dMob_dP * potGrad;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseFlux_dC[k_up][jc] += dPhaseMob_dCompSub[jc] * potGrad;
    }

    // slice some constitutive arrays to avoid too much indexing in component loop
    arraySlice1d< real64 const > phaseCompFracSub =
      phaseCompFrac[er_up][esr_up][ei_up][0][ip];
    arraySlice1d< real64 const > dPhaseCompFrac_dPresSub =
      dPhaseCompFrac_dPres[er_up][esr_up][ei_up][0][ip];
    arraySlice2d< real64 const > dPhaseCompFrac_dCompSub =
      dPhaseCompFrac_dComp[er_up][esr_up][ei_up][0][ip];

    // compute component fluxes and derivatives using upstream cell composition
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const ycp = phaseCompFracSub[ic];
      compFlux[ic] += phaseFlux * ycp;

      // derivatives stemming from phase flux
      for( localIndex ke = 0; ke < stencilSize; ++ke )
      {
        dCompFlux_dP[ke][ic] += dPhaseFlux_dP[ke] * ycp;
        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dCompFlux_dC[ke][ic][jc] += dPhaseFlux_dC[ke][jc] * ycp;
        }
      }

      // additional derivatives stemming from upstream cell phase composition
      dCompFlux_dP[k_up][ic] += phaseFlux * dPhaseCompFrac_dPresSub[ic];

      // convert derivatives of component fraction w.r.t. component fractions to derivatives w.r.t. component
      // densities
      applyChainRule( NC,
                      dCompFrac_dCompDens[er_up][esr_up][ei_up],
                      dPhaseCompFrac_dCompSub[ic],
                      dProp_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dCompFlux_dC[k_up][ic][jc] += phaseFlux * dProp_dC[jc];
      }
    }
  }

  // *** end of upwinding

  // populate local flux vector and derivatives
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
        localFluxJacobian[NC + ic][localDofIndexComp] =
          -dt * dCompFlux_dC[ke][ic][jc];
      }
    }
  }
}

template< localIndex NC, typename STENCIL_TYPE >
void
FluxKernel::Launch(
  localIndex const numPhases,
  STENCIL_TYPE const & stencil,
  globalIndex const rankOffset,
  ElementView< arrayView1d< globalIndex const > > const & dofNumber,
  ElementView< arrayView1d< integer const > > const & ghostRank,
  ElementView< arrayView1d< real64 const > > const & pres,
  ElementView< arrayView1d< real64 const > > const & dPres,
  ElementView< arrayView1d< real64 const > > const & gravCoef,
  ElementView< arrayView2d< real64 const > > const & phaseMob,
  ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,
  ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp,
  ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
  ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
  ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
  ElementView< arrayView3d< real64 const > > const & phaseDens,
  ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,
  ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp,
  ElementView< arrayView4d< real64 const > > const & phaseCompFrac,
  ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
  ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
  ElementView< arrayView3d< real64 const > > const & phaseCapPressure,
  ElementView< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
  integer const capPressureFlag,
  real64 const dt,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  typename STENCIL_TYPE::IndexContainerViewConstType const & seri =
    stencil.getElementRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sesri =
    stencil.getElementSubRegionIndices();
  typename STENCIL_TYPE::IndexContainerViewConstType const & sei =
    stencil.getElementIndices();
  typename STENCIL_TYPE::WeightContainerViewConstType const & weights =
    stencil.getWeights();

  localIndex constexpr NUM_ELEMS = STENCIL_TYPE::NUM_POINT_IN_FLUX;
  localIndex constexpr MAX_STENCIL = STENCIL_TYPE::MAX_STENCIL_SIZE;
  localIndex constexpr NDOF = NC + 1;

  forAll< parallelDevicePolicy<> >(
    stencil.size(),
    [=] GEOSX_HOST_DEVICE( localIndex const iconn ) {
      // TODO: hack! for MPFA, etc. must obtain proper size from e.g. seri
      localIndex const stencilSize = MAX_STENCIL;

      stackArray1d< real64, NUM_ELEMS * NC > localFlux( NUM_ELEMS * NC );
      stackArray2d< real64, NUM_ELEMS * NC * MAX_STENCIL * NDOF > localFluxJacobian(
        NUM_ELEMS * NC,
        stencilSize * NDOF );

      FluxKernel::Compute< NC, NUM_ELEMS, MAX_STENCIL >(
        numPhases,
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
      globalIndex dofColIndices[MAX_STENCIL * NDOF];
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        globalIndex const offset =
          dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

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
          globalIndex const globalRow =
            dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow =
            LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow + NC );

          for( localIndex ic = 0; ic < NC; ++ic )
          {
            RAJA::atomicAdd( parallelDeviceAtomic {},
                             &localRhs[localRow + ic],
                             localFlux[i * NC + ic] );
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
              localRow + ic,
              dofColIndices,
              localFluxJacobian[i * NC + ic].dataIfContiguous(),
              stencilSize * NDOF );
          }
        }
      }
    } );
}

#define INST_FluxKernel( NC, STENCIL_TYPE )                                             \
  template void FluxKernel::Launch< NC, STENCIL_TYPE >(                                 \
    localIndex const numPhases,                                                         \
    STENCIL_TYPE const & stencil,                                                       \
    globalIndex const rankOffset,                                                       \
    ElementView< arrayView1d< globalIndex const > > const & dofNumber,                  \
    ElementView< arrayView1d< integer const > > const & ghostRank,                      \
    ElementView< arrayView1d< real64 const > > const & pres,                            \
    ElementView< arrayView1d< real64 const > > const & dPres,                           \
    ElementView< arrayView1d< real64 const > > const & gravCoef,                        \
    ElementView< arrayView2d< real64 const > > const & phaseMob,                        \
    ElementView< arrayView2d< real64 const > > const & dPhaseMob_dPres,                 \
    ElementView< arrayView3d< real64 const > > const & dPhaseMob_dComp,                 \
    ElementView< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,             \
    ElementView< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,             \
    ElementView< arrayView3d< real64 const > > const & dCompFrac_dCompDens,             \
    ElementView< arrayView3d< real64 const > > const & phaseDens,                       \
    ElementView< arrayView3d< real64 const > > const & dPhaseDens_dPres,                \
    ElementView< arrayView4d< real64 const > > const & dPhaseDens_dComp,                \
    ElementView< arrayView4d< real64 const > > const & phaseCompFrac,                   \
    ElementView< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,            \
    ElementView< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,            \
    ElementView< arrayView3d< real64 const > > const & phaseCapPressure,                \
    ElementView< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac, \
    integer const capPressureFlag,                                                      \
    real64 const dt,                                                                    \
    CRSMatrixView< real64, globalIndex const > const & localMatrix,                     \
    arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 1, CellElementStencilTPFA );
INST_FluxKernel( 2, CellElementStencilTPFA );
INST_FluxKernel( 3, CellElementStencilTPFA );
INST_FluxKernel( 4, CellElementStencilTPFA );
INST_FluxKernel( 5, CellElementStencilTPFA );

INST_FluxKernel( 1, FaceElementStencil );
INST_FluxKernel( 2, FaceElementStencil );
INST_FluxKernel( 3, FaceElementStencil );
INST_FluxKernel( 4, FaceElementStencil );
INST_FluxKernel( 5, FaceElementStencil );

#undef INST_FluxKernel

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
VolumeBalanceKernel::Compute(
  real64 const & volume,
  real64 const & porosityRef,
  real64 const & pvMult,
  real64 const & dPvMult_dPres,
  arraySlice1d< real64 const > const & phaseVolFrac,
  arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
  arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
  real64 & localVolBalance,
  real64 * const localVolBalanceJacobian )
{
  localIndex constexpr NDOF = NC + 1;

  real64 const poro = porosityRef * pvMult;
  real64 const dPoro_dP = porosityRef * dPvMult_dPres;

  real64 const poreVol = volume * poro;
  real64 const dPoreVol_dP = volume * dPoro_dP;

  localVolBalance = 1.0;
  for( localIndex i = 0; i < NDOF; ++i )
  {
    localVolBalanceJacobian[i] = 0.0;
  }

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localVolBalance -= phaseVolFrac[ip];
    localVolBalanceJacobian[0] -= dPhaseVolFrac_dPres[ip];

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      localVolBalanceJacobian[jc + 1] -= dPhaseVolFrac_dCompDens[ip][jc];
    }
  }

  // scale saturation-based volume balance by pore volume (for better scaling w.r.t. other equations)
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    localVolBalanceJacobian[idof] *= poreVol;
  }
  localVolBalanceJacobian[0] += dPoreVol_dP * localVolBalance;
  localVolBalance *= poreVol;
}

template< localIndex NC, localIndex NP >
void
VolumeBalanceKernel::Launch(
  localIndex const size,
  globalIndex const rankOffset,
  arrayView1d< globalIndex const > const & dofNumber,
  arrayView1d< integer const > const & elemGhostRank,
  arrayView1d< real64 const > const & volume,
  arrayView1d< real64 const > const & porosityRef,
  arrayView2d< real64 const > const & pvMult,
  arrayView2d< real64 const > const & dPvMult_dPres,
  arrayView2d< real64 const > const & phaseVolFrac,
  arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
    if( elemGhostRank[ei] >= 0 ) return;

    localIndex constexpr NDOF = NC + 1;

    real64 localVolBalance;
    real64 localVolBalanceJacobian[NDOF];

    Compute< NC, NP >( volume[ei],
                       porosityRef[ei],
                       pvMult[ei][0],
                       dPvMult_dPres[ei][0],
                       phaseVolFrac[ei],
                       dPhaseVolFrac_dPres[ei],
                       dPhaseVolFrac_dCompDens[ei],
                       localVolBalance,
                       localVolBalanceJacobian );

    // get equation/dof indices
    localIndex const localRow = dofNumber[ei] + NC - rankOffset;
    globalIndex dofIndices[NDOF];
    for( localIndex jdof = 0; jdof < NDOF; ++jdof )
    {
      dofIndices[jdof] = dofNumber[ei] + jdof;
    }

    // TODO: apply equation/variable change transformation(s)

    // add contribution to residual and jacobian
    localRhs[localRow] += localVolBalance;
    localMatrix.addToRow< serialAtomic >( localRow,
                                          dofIndices,
                                          localVolBalanceJacobian,
                                          NDOF );
  } );
}

#define INST_VolumeBalanceKernel( NC, NP )                          \
  template void VolumeBalanceKernel::Launch< NC, NP >(              \
    localIndex const size,                                          \
    globalIndex const rankOffset,                                   \
    arrayView1d< globalIndex const > const & dofNumber,             \
    arrayView1d< integer const > const & elemGhostRank,             \
    arrayView1d< real64 const > const & volume,                     \
    arrayView1d< real64 const > const & porosityRef,                \
    arrayView2d< real64 const > const & pvMult,                     \
    arrayView2d< real64 const > const & dPvMult_dPres,              \
    arrayView2d< real64 const > const & phaseVolFrac,               \
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres,        \
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,    \
    CRSMatrixView< real64, globalIndex const > const & localMatrix, \
    arrayView1d< real64 > const & localRhs )

INST_VolumeBalanceKernel( 1, 1 );
INST_VolumeBalanceKernel( 2, 1 );
INST_VolumeBalanceKernel( 3, 1 );
INST_VolumeBalanceKernel( 4, 1 );
INST_VolumeBalanceKernel( 5, 1 );

INST_VolumeBalanceKernel( 1, 2 );
INST_VolumeBalanceKernel( 2, 2 );
INST_VolumeBalanceKernel( 3, 2 );
INST_VolumeBalanceKernel( 4, 2 );
INST_VolumeBalanceKernel( 5, 2 );

INST_VolumeBalanceKernel( 1, 3 );
INST_VolumeBalanceKernel( 2, 3 );
INST_VolumeBalanceKernel( 3, 3 );
INST_VolumeBalanceKernel( 4, 3 );
INST_VolumeBalanceKernel( 5, 3 );

#undef INST_VolumeBalanceKernel

}  // namespace CompositionalMultiphaseFlowKernels

}  // namespace geosx
