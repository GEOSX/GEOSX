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
 * @file CompositionalMultiphaseBaseKernels.cpp
 */

#include "CompositionalMultiphaseBaseKernels.hpp"

namespace geosx
{

namespace CompositionalMultiphaseBaseKernels
{

/******************************** ComponentFractionKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ComponentFractionKernel::
  compute( arraySlice1d< real64 const > const compDens,
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
    compFrac[ic] = (compDens[ic] + dCompDens[ic]) * totalDensityInv;
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dCompFrac_dCompDens[ic][jc] = -compFrac[ic] * totalDensityInv;
    }
    dCompFrac_dCompDens[ic][ic] += totalDensityInv;
  }
}

template< localIndex NC >
void
ComponentFractionKernel::
  launch( localIndex const size,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

template< localIndex NC >
void
ComponentFractionKernel::
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC >( compDens[a],
                   dCompDens[a],
                   compFrac[a],
                   dCompFrac_dCompDens[a] );
  } );
}

#define INST_ComponentFractionKernel( NC ) \
  template \
  void ComponentFractionKernel:: \
    launch< NC >( localIndex const size, \
                  arrayView2d< real64 const > const & compDens, \
                  arrayView2d< real64 const > const & dCompDens, \
                  arrayView2d< real64 > const & compFrac, \
                  arrayView3d< real64 > const & dCompFrac_dCompDens ); \
  template \
  void ComponentFractionKernel:: \
    launch< NC >( SortedArrayView< localIndex const > const & targetSet, \
                  arrayView2d< real64 const > const & compDens, \
                  arrayView2d< real64 const > const & dCompDens, \
                  arrayView2d< real64 > const & compFrac, \
                  arrayView3d< real64 > const & dCompFrac_dCompDens )

INST_ComponentFractionKernel( 1 );
INST_ComponentFractionKernel( 2 );
INST_ComponentFractionKernel( 3 );
INST_ComponentFractionKernel( 4 );
INST_ComponentFractionKernel( 5 );

#undef INST_ComponentFractionKernel

/******************************** PhaseVolumeFractionKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
PhaseVolumeFractionKernel::
  compute( arraySlice1d< real64 const > const & compDens,
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
      (dPhaseFrac_dPres[ip] - phaseVolFrac[ip] * dPhaseDens_dPres[ip]) * phaseDensInv;

    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseVolFrac_dComp[ip][jc] =
        (dPhaseFrac_dComp[ip][jc] - phaseVolFrac[ip] * dPhaseDens_dComp[ip][jc]) * phaseDensInv;
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
void PhaseVolumeFractionKernel::
  launch( localIndex const size,
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
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC, NP >( compDens[a],
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
void PhaseVolumeFractionKernel::
  launch( SortedArrayView< localIndex const > const & targetSet,
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
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC, NP >( compDens[a],
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

#define INST_PhaseVolumeFractionKernel( NC, NP ) \
  template \
  void \
  PhaseVolumeFractionKernel:: \
    launch< NC, NP >( localIndex const size, \
                      arrayView2d< real64 const > const & compDens, \
                      arrayView2d< real64 const > const & dCompDens, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseFrac, \
                      arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                      arrayView2d< real64 > const & phaseVolFrac, \
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 > const & dPhaseVolFrac_dComp ); \
  template \
  void \
  PhaseVolumeFractionKernel:: \
    launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView2d< real64 const > const & compDens, \
                      arrayView2d< real64 const > const & dCompDens, \
                      arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const > const & phaseDens, \
                      arrayView3d< real64 const > const & dPhaseDens_dPres, \
                      arrayView4d< real64 const > const & dPhaseDens_dComp, \
                      arrayView3d< real64 const > const & phaseFrac, \
                      arrayView3d< real64 const > const & dPhaseFrac_dPres, \
                      arrayView4d< real64 const > const & dPhaseFrac_dComp, \
                      arrayView2d< real64 > const & phaseVolFrac, \
                      arrayView2d< real64 > const & dPhaseVolFrac_dPres, \
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


/******************************** AccumulationKernel ********************************/

template< localIndex NC >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
AccumulationKernel::
  compute( localIndex const numPhases,
           real64 const & poreVolOld,
           real64 const & poreVolNew,
           real64 const & dPoreVol_dP,
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
           real64 ( & localAccum )[NC],
           real64 ( & localAccumJacobian )[NC][NC + 1] )
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

  // sum contributions to component accumulation from each phase
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const phaseAmountNew = poreVolNew * phaseVolFrac[ip] * phaseDens[ip];
    real64 const phaseAmountOld = poreVolOld * phaseVolFracOld[ip] * phaseDensOld[ip];

    real64 const dPhaseAmount_dP = dPoreVol_dP * phaseVolFrac[ip] * phaseDens[ip]
                                   + poreVolNew * (dPhaseVolFrac_dPres[ip] * phaseDens[ip]
                                                   + phaseVolFrac[ip] * dPhaseDens_dPres[ip]);

    // assemble density dependence
    applyChainRule( NC, dCompFrac_dCompDens, dPhaseDens_dComp[ip], dPhaseAmount_dC );
    for( localIndex jc = 0; jc < NC; ++jc )
    {
      dPhaseAmount_dC[jc] = dPhaseAmount_dC[jc] * phaseVolFrac[ip]
                            + phaseDens[ip] * dPhaseVolFrac_dCompDens[ip][jc];
      dPhaseAmount_dC[jc] *= poreVolNew;
    }

    // ic - index of component whose conservation equation is assembled
    // (i.e. row number in local matrix)
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      real64 const phaseCompAmountNew = phaseAmountNew * phaseCompFrac[ip][ic];
      real64 const phaseCompAmountOld = phaseAmountOld * phaseCompFracOld[ip][ic];

      real64 const dPhaseCompAmount_dP = dPhaseAmount_dP * phaseCompFrac[ip][ic]
                                         + phaseAmountNew * dPhaseCompFrac_dPres[ip][ic];

      localAccum[ic] += phaseCompAmountNew - phaseCompAmountOld;
      localAccumJacobian[ic][0] += dPhaseCompAmount_dP;

      // jc - index of component w.r.t. whose compositional var the derivative is being taken
      // (i.e. col number in local matrix)

      // assemble phase composition dependence
      applyChainRule( NC, dCompFrac_dCompDens, dPhaseCompFrac_dComp[ip][ic], dPhaseCompFrac_dC );
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        real64 const dPhaseCompAmount_dC = dPhaseCompFrac_dC[jc] * phaseAmountNew
                                           + phaseCompFrac[ip][ic] * dPhaseAmount_dC[jc];
        localAccumJacobian[ic][jc + 1] += dPhaseCompAmount_dC;
      }
    }
  }
}

template< localIndex NC >
void
AccumulationKernel::
  launch( localIndex const numPhases,
          localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosityOld,
          arrayView2d< real64 const > const & porosityNew,
          arrayView2d< real64 const > const & dPoro_dPres,
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
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] >= 0 )
      return;

    localIndex constexpr NDOF = NC + 1;

    real64 localAccum[NC];
    real64 localAccumJacobian[NC][NDOF];

    real64 const poreVolumeNew = volume[ei] * porosityNew[ei][0];
    real64 const poreVolumeOld = volume[ei] * porosityOld[ei][0];
    real64 const dPoreVolume_dPres = volume[ei] * dPoro_dPres[ei][0];

    compute< NC >( numPhases,
                   poreVolumeOld,
                   poreVolumeNew,
                   dPoreVolume_dPres,
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

#define INST_AccumulationKernel( NC ) \
  template \
  void \
  AccumulationKernel:: \
    launch< NC >( localIndex const numPhases, \
                  localIndex const size, \
                  globalIndex const rankOffset, \
                  arrayView1d< globalIndex const > const & dofNumber, \
                  arrayView1d< integer const > const & elemGhostRank, \
                  arrayView1d< real64 const > const & volume, \
                  arrayView2d< real64 const > const & porosityOld, \
                  arrayView2d< real64 const > const & porosityNew, \
                  arrayView2d< real64 const > const & dPoro_dPres, \
                  arrayView3d< real64 const > const & dCompFrac_dCompDens, \
                  arrayView2d< real64 const > const & phaseVolFracOld, \
                  arrayView2d< real64 const > const & phaseVolFrac, \
                  arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens, \
                  arrayView2d< real64 const > const & phaseDensOld, \
                  arrayView3d< real64 const > const & phaseDens, \
                  arrayView3d< real64 const > const & dPhaseDens_dPres, \
                  arrayView4d< real64 const > const & dPhaseDens_dComp, \
                  arrayView3d< real64 const > const & phaseCompFracOld, \
                  arrayView4d< real64 const > const & phaseCompFrac, \
                  arrayView4d< real64 const > const & dPhaseCompFrac_dPres, \
                  arrayView5d< real64 const > const & dPhaseCompFrac_dComp, \
                  CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                  arrayView1d< real64 > const & localRhs )

INST_AccumulationKernel( 1 );
INST_AccumulationKernel( 2 );
INST_AccumulationKernel( 3 );
INST_AccumulationKernel( 4 );
INST_AccumulationKernel( 5 );

#undef INST_AccumulationKernel

/******************************** VolumeBalanceKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
VolumeBalanceKernel::
  compute( real64 const & volume,
           real64 const & porosity,
           real64 const & dPoro_dPres,
           arraySlice1d< real64 const > const & phaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           real64 * const localVolBalanceJacobian )
{
  localIndex constexpr NDOF = NC + 1;

  real64 const poreVol     = volume * porosity;
  real64 const dPoreVol_dP = volume * dPoro_dPres;

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
      localVolBalanceJacobian[jc+1] -= dPhaseVolFrac_dCompDens[ip][jc];
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
VolumeBalanceKernel::
  launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & volume,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const > const & dPoro_dPres,
          arrayView2d< real64 const > const & phaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const ei )
  {
    if( elemGhostRank[ei] >= 0 )
      return;

    localIndex constexpr NDOF = NC + 1;

    real64 localVolBalance;
    real64 localVolBalanceJacobian[NDOF];

    compute< NC, NP >( volume[ei],
                       porosity[ei][0],
                       dPoro_dPres[ei][0],
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

#define INST_VolumeBalanceKernel( NC, NP ) \
  template \
  void VolumeBalanceKernel:: \
    launch< NC, NP >( localIndex const size, \
                      globalIndex const rankOffset, \
                      arrayView1d< globalIndex const > const & dofNumber, \
                      arrayView1d< integer const > const & elemGhostRank, \
                      arrayView1d< real64 const > const & volume, \
                      arrayView2d< real64 const > const & porosity, \
                      arrayView2d< real64 const > const & dPoro_dPres, \
                      arrayView2d< real64 const > const & phaseVolFrac, \
                      arrayView2d< real64 const > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens, \
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

} // namespace CompositionalMultiphaseFlowKernels

} // namespace geosx
