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
 * @file CompositionalMultiphaseHybridFVMKernels.cpp
 */

#include "CompositionalMultiphaseHybridFVMKernels.hpp"

#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/BdVLMInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"

#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{

/******************************** UpwindingHelper ********************************/

template< localIndex NC, localIndex NP >
void
UpwindingHelper::
  upwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                            localIndex const (&neighborIds)[ 3 ],
                            ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                            ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
                            ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                            ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
                            ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
                            ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                            real64 const & oneSidedVolFlux,
                            real64 ( & upwPhaseViscCoef )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                            real64 ( & dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                            globalIndex & upwViscDofNumber )
{
  real64 dUpwMobRatio_dCompDens[ NC ]{};
  real64 dUpwDensMobRatio_dCompDens[ NC ]{};
  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  // 1) Upwind
  localIndex const er  = ( oneSidedVolFlux > 0 ) ? localIds[0] : neighborIds[0];
  localIndex const esr = ( oneSidedVolFlux > 0 ) ? localIds[1] : neighborIds[1];
  localIndex const ei  = ( oneSidedVolFlux > 0 ) ? localIds[2] : neighborIds[2];

  // 2) Compute total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0;
  real64 dTotalMob_dPres = 0;
  real64 dTotalMob_dCompDens[ NC ]{};
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres = dTotalMob_dPres + dPhaseMob_dPres[er][esr][ei][ip];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[ic] = dTotalMob_dCompDens[ic] + dPhaseMob_dCompDens[er][esr][ei][ip][ic];
    }
  }
  real64 const totalMobInv = 1.0 / totalMob;

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    // 3) Compute viscous mobility ratio: \frac{\lambda_{\ell}}{\lambda_T}
    real64 const upwMobRatio = phaseMob[er][esr][ei][ip] * totalMobInv;
    real64 const dUpwMobRatio_dPres = ( dPhaseMob_dPres[er][esr][ei][ip] - upwMobRatio * dTotalMob_dPres )
                                      * totalMobInv;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dUpwMobRatio_dCompDens[ic] = ( dPhaseMob_dCompDens[er][esr][ei][ip][ic] - upwMobRatio * dTotalMob_dCompDens[ic] )
                                   * totalMobInv;
    }

    // 4) Multiply mobility ratio by phase density: \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseDens_dCompFrac[er][esr][ei][0][ip],
                    dPhaseDens_dC );
    real64 const upwDensMobRatio = phaseDens[er][esr][ei][0][ip] * upwMobRatio;
    real64 const dUpwDensMobRatio_dPres = dPhaseDens_dPres[er][esr][ei][0][ip] * upwMobRatio
                                          + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dPres;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dUpwDensMobRatio_dCompDens[ic] = dPhaseDens_dC[ic] * upwMobRatio
                                       + phaseDens[er][esr][ei][0][ip] * dUpwMobRatio_dCompDens[ic];
    }

    // 5) Multiply density mobility ratio by phase comp fraction: x_{c,\ell} \rho^{up}_{\ell} \frac{\lambda_{\ell}}{\lambda_T}
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseCompFrac_dCompFrac[er][esr][ei][0][ip][ic],
                      dPhaseCompFrac_dC );
      upwPhaseViscCoef[ip][ic] = phaseCompFrac[er][esr][ei][0][ip][ic] * upwDensMobRatio;
      dUpwPhaseViscCoef_dPres[ip][ic] = dPhaseCompFrac_dPres[er][esr][ei][0][ip][ic] * upwDensMobRatio
                                        + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dPres;
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dUpwPhaseViscCoef_dCompDens[ip][ic][jc] = dPhaseCompFrac_dC[jc] * upwDensMobRatio
                                                  + phaseCompFrac[er][esr][ei][0][ip][ic] * dUpwDensMobRatio_dCompDens[jc];
      }
    }
  }
  // 6) Save the dof number of the upwind cell
  upwViscDofNumber = elemDofNumber[er][esr][ei];

}

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::
  upwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                             localIndex const (&neighborIds)[ 3 ],
                             real64 const & transGravCoef,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                             ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac,
                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                             ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
                             ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                             ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
                             ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
                             real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                             real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                             real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                             real64 ( & upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                             real64 ( & dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                             real64 ( & dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] )
{
  // 1) Compute the driving force: T ( \rho^{avg}_{\ell} - \rho^{avg}_m ) g \Delta z
  computePhaseGravTerm( localIds,
                        neighborIds,
                        transGravCoef,
                        phaseMassDens,
                        dPhaseMassDens_dPres,
                        dPhaseMassDens_dCompFrac,
                        dCompFrac_dCompDens,
                        phaseGravTerm,
                        dPhaseGravTerm_dPres,
                        dPhaseGravTerm_dCompDens );

  // 2) Compute the total mobility: \lambda_T = \sum_{\ell} \lambda_{\ell}
  real64 totalMob = 0.0;
  real64 dTotalMob_dPres[ 2 ]{};
  real64 dTotalMob_dCompDens[ 2 ][ NC ]{};

  computeUpwindedTotalMobility( localIds,
                                neighborIds,
                                phaseMob,
                                dPhaseMob_dPres,
                                dPhaseMob_dCompDens,
                                phaseGravTerm,
                                totalMob,
                                dTotalMob_dPres,
                                dTotalMob_dCompDens );
  real64 const totalMobInv = 1.0 / totalMob;

  // 3) Compute the quantities \x_{up}_{c,p} \rho_p \frac{\lambda_p \lambda_m}{\lambda_T}
  real64 dMobRatio_dPres[ 2 ]{};
  real64 dMobRatio_dCompDens[ 2 ][ NC ]{};
  real64 dDensMobRatio_dPres[ 2 ]{};
  real64 dDensMobRatio_dCompDens[ 2 ][ NC ]{};

  real64 dPhaseDens_dC[ NC ]{};
  real64 dPhaseCompFrac_dC[ NC ]{};

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localIndex k = 0;
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      // 3.a) Upwinding using the gravity term
      localIndex eru, esru, eiu, posu; // upwind
      localIndex erd, esrd, eid, posd; // downwind
      setIndicesForMobilityRatioUpwinding( localIds, neighborIds,
                                           phaseGravTerm[ip][k],
                                           eru, esru, eiu, posu,
                                           erd, esrd, eid, posd );

      // 3.b) Compute mobility ratio \frac{\lambda_l \lambda_m}{\lambda_T}
      real64 const mobRatio = phaseMob[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                              * totalMobInv;
      dMobRatio_dPres[posu] = ( dPhaseMob_dPres[eru][esru][eiu][ip] * phaseMob[erd][esrd][eid][jp]
                                - mobRatio * dTotalMob_dPres[posu] ) * totalMobInv;
      dMobRatio_dPres[posd] = ( dPhaseMob_dPres[erd][esrd][eid][jp] * phaseMob[eru][esru][eiu][ip]
                                - mobRatio * dTotalMob_dPres[posd] ) * totalMobInv;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dMobRatio_dCompDens[posu][ic] = ( dPhaseMob_dCompDens[eru][esru][eiu][ip][ic] * phaseMob[erd][esrd][eid][jp]
                                          - mobRatio * dTotalMob_dCompDens[posu][ic] ) * totalMobInv;
        dMobRatio_dCompDens[posd][ic] = ( dPhaseMob_dCompDens[erd][esrd][eid][jp][ic] * phaseMob[eru][esru][eiu][ip]
                                          - mobRatio * dTotalMob_dCompDens[posd][ic] ) * totalMobInv;
      }

      // 3.c) Compute mobility ratio multiplied by upwinded phase density \rho_l \frac{\lambda_l \lambda_m}{\lambda_T}
      applyChainRule( NC,
                      dCompFrac_dCompDens[eru][esru][eiu],
                      dPhaseDens_dCompFrac[eru][esru][eiu][0][ip],
                      dPhaseDens_dC );
      real64 const densMobRatio = phaseDens[eru][esru][eiu][0][ip] * mobRatio;
      dDensMobRatio_dPres[posu] = dPhaseDens_dPres[eru][esru][eiu][0][ip] * mobRatio
                                  + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posu];
      dDensMobRatio_dPres[posd] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dPres[posd];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dDensMobRatio_dCompDens[posu][ic] = dPhaseDens_dC[ic] * mobRatio
                                            + phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posu][ic];
        dDensMobRatio_dCompDens[posd][ic] = phaseDens[eru][esru][eiu][0][ip] * dMobRatio_dCompDens[posd][ic];
      }

      // 3.d) Compute the final gravity coefficient \x_{up}_{c,p} \rho_p \frac{\lambda_l \lambda_m}{\lambda_T}
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        applyChainRule( NC,
                        dCompFrac_dCompDens[eru][esru][eiu],
                        dPhaseCompFrac_dCompFrac[eru][esru][eiu][0][ip][ic],
                        dPhaseCompFrac_dC );
        upwPhaseGravCoef[ip][k][ic] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * densMobRatio;
        dUpwPhaseGravCoef_dPres[ip][k][ic][posu] = dPhaseCompFrac_dPres[eru][esru][eiu][0][ip][ic] * densMobRatio
                                                   + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posu];
        dUpwPhaseGravCoef_dPres[ip][k][ic][posd] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dPres[posd];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posu][jc] = dPhaseCompFrac_dC[jc] * densMobRatio
                                                             + phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posu][jc];
          dUpwPhaseGravCoef_dCompDens[ip][k][ic][posd][jc] = phaseCompFrac[eru][esru][eiu][0][ip][ic] * dDensMobRatio_dCompDens[posd][jc];
        }
      }
      ++k;
    }
  }
}

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::
  computePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                        localIndex const (&neighborIds)[ 3 ],
                        real64 const & transGravCoef,
                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac,
                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                        real64 ( & phaseGravTerm )[ NP ][ NP-1 ],
                        real64 ( & dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                        real64 ( & dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] )
{
  localIndex const er   = localIds[0];
  localIndex const esr  = localIds[1];
  localIndex const ei   = localIds[2];
  localIndex const ern  = neighborIds[0];
  localIndex const esrn = neighborIds[1];
  localIndex const ein  = neighborIds[2];

  real64 dPhaseMassDens_dCLoc[ NC ]{};
  real64 dPhaseMassDens_dCNeighbor[ NC ]{};
  real64 dPhaseMassDens_dC[ NC ]{};

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dCompFrac_dCompDens[er][esr][ei],
                    dPhaseMassDens_dCompFrac[er][esr][ei][0][ip],
                    dPhaseMassDens_dCLoc );
    applyChainRule( NC,
                    dCompFrac_dCompDens[ern][esrn][ein],
                    dPhaseMassDens_dCompFrac[ern][esrn][ein][0][ip],
                    dPhaseMassDens_dCNeighbor );

    localIndex k = 0;
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      if( ip == jp )
      {
        continue;
      }

      phaseGravTerm[ip][k] = -( phaseMassDens[er][esr][ei][0][ip] + phaseMassDens[ern][esrn][ein][0][ip] );
      phaseGravTerm[ip][k] += ( phaseMassDens[er][esr][ei][0][jp] + phaseMassDens[ern][esrn][ein][0][jp] );
      phaseGravTerm[ip][k] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] = ( -dPhaseMassDens_dPres[er][esr][ei][0][ip] + dPhaseMassDens_dPres[er][esr][ei][0][jp] );
      dPhaseGravTerm_dPres[ip][k][Pos::LOCAL] *= 0.5 * transGravCoef;

      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] = ( -dPhaseMassDens_dPres[ern][esrn][ein][0][ip] + dPhaseMassDens_dPres[ern][esrn][ein][0][jp] );
      dPhaseGravTerm_dPres[ip][k][Pos::NEIGHBOR] *= 0.5 * transGravCoef;

      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCLoc[ic];
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] = -0.5 * transGravCoef * dPhaseMassDens_dCNeighbor[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[er][esr][ei],
                      dPhaseMassDens_dCompFrac[er][esr][ei][0][jp],
                      dPhaseMassDens_dC );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::LOCAL][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      applyChainRule( NC,
                      dCompFrac_dCompDens[ern][esrn][ein],
                      dPhaseMassDens_dCompFrac[ern][esrn][ein][0][jp],
                      dPhaseMassDens_dC );
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dPhaseGravTerm_dCompDens[ip][k][Pos::NEIGHBOR][ic] += 0.5 * transGravCoef * dPhaseMassDens_dC[ic];
      }
      ++k;
    }
  }
}

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::
  computeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                localIndex const (&neighborIds)[ 3 ],
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                                ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
                                ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
                                real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                real64 & totalMob,
                                real64 ( & dTotalMob_dPres )[ 2 ],
                                real64 ( & dTotalMob_dCompDens )[ 2 ][ NC ] )
{
  localIndex totalMobIds[ NP ][ 3 ]{};
  localIndex totalMobPos[ NP ]{};
  setIndicesForTotalMobilityUpwinding< NP >( localIds,
                                             neighborIds,
                                             phaseGravTerm,
                                             totalMobIds,
                                             totalMobPos );
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    localIndex const er  = totalMobIds[ip][0];
    localIndex const esr = totalMobIds[ip][1];
    localIndex const ei  = totalMobIds[ip][2];
    localIndex const pos = totalMobPos[ip];
    totalMob = totalMob + phaseMob[er][esr][ei][ip];
    dTotalMob_dPres[pos] = dTotalMob_dPres[pos] + dPhaseMob_dPres[er][esr][ei][pos];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dTotalMob_dCompDens[pos][ic] = dTotalMob_dCompDens[pos][ic] + dPhaseMob_dCompDens[er][esr][ei][ip][ic];
    }
  }
  if( totalMob < 1e-12 )
  {
    totalMob = 1e-12;
  }
}

GEOSX_HOST_DEVICE
void
UpwindingHelper::
  setIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd )
{
  if( gravTerm > 0 )
  {
    eru  = localIds[0];
    esru = localIds[1];
    eiu  = localIds[2];
    posu = Pos::LOCAL;
    erd  = neighborIds[0];
    esrd = neighborIds[1];
    eid  = neighborIds[2];
    posd = Pos::NEIGHBOR;
  }
  else
  {
    eru  = neighborIds[0];
    esru = neighborIds[1];
    eiu  = neighborIds[2];
    posu = Pos::NEIGHBOR;
    erd  = localIds[0];
    esrd = localIds[1];
    eid  = localIds[2];
    posd = Pos::LOCAL;
  }
}

template< localIndex NP >
GEOSX_HOST_DEVICE
void
UpwindingHelper::
  setIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const (&gravTerm)[ NP ][ NP-1 ],
                                       localIndex ( & totalMobIds )[ NP ][ 3 ],
                                       localIndex ( & totalMobPos )[ NP ] )
{
  if( NP == 2 )
  {
    if( gravTerm[0][0] > 0 )
    {
      totalMobIds[0][0] = localIds[0];
      totalMobIds[0][1] = localIds[1];
      totalMobIds[0][2] = localIds[2];
      totalMobPos[0] = Pos::LOCAL;
      totalMobIds[1][0] = neighborIds[0];
      totalMobIds[1][1] = neighborIds[1];
      totalMobIds[1][2] = neighborIds[2];
      totalMobPos[1] = Pos::NEIGHBOR;
    }
    else
    {
      totalMobIds[0][0] = neighborIds[0];
      totalMobIds[0][1] = neighborIds[1];
      totalMobIds[0][2] = neighborIds[2];
      totalMobPos[0] = Pos::NEIGHBOR;
      totalMobIds[1][0] = localIds[0];
      totalMobIds[1][1] = localIds[1];
      totalMobIds[1][2] = localIds[2];
      totalMobPos[1] = Pos::LOCAL;
    }
  }
  else if( NP == 3 )
  {
    // TODO Francois: this should be improved
    // currently this implements the algorithm proposed by SH Lee
    for( localIndex ip = 0; ip < NP; ++ip )
    {
      if( ( gravTerm[ip][0] >= 0 && gravTerm[ip][1] >= 0 ) || // includes the no-buoyancy case
          ( ( fabs( gravTerm[ip][0] ) >= fabs( gravTerm[ip][1] ) ) && gravTerm[ip][1] >= 0 ) ||
          ( ( fabs( gravTerm[ip][1] ) >= fabs( gravTerm[ip][0] ) ) && gravTerm[ip][0] >= 0 ) )
      {
        totalMobIds[ip][0] = localIds[0];
        totalMobIds[ip][1] = localIds[1];
        totalMobIds[ip][2] = localIds[2];
        totalMobPos[ip] = Pos::LOCAL;
      }
      else
      {
        totalMobIds[ip][0] = neighborIds[0];
        totalMobIds[ip][1] = neighborIds[1];
        totalMobIds[ip][2] = neighborIds[2];
        totalMobPos[ip] = Pos::NEIGHBOR;
      }
    }
  }
}


#define INST_UpwindingHelperNCNP( NC, NP ) \
  template \
  void \
  UpwindingHelper:: \
    upwindViscousCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                        localIndex const (&neighborIds)[ 3 ], \
                                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                        ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac, \
                                        ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                        ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                                        ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                        ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                                        ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac, \
                                        ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                        real64 const & oneSidedVolFlux, \
                                        real64 ( &upwPhaseViscCoef )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ], \
                                        real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ], \
                                        globalIndex & upwViscDofNumber ); \
  template \
  void \
  UpwindingHelper:: \
    upwindBuoyancyCoefficient< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                         localIndex const (&neighborIds)[ 3 ], \
                                         real64 const & transGravCoef,  \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac, \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                         ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac, \
                                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                         ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                                         ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                         ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                                         ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac, \
                                         real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                         real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                         real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ], \
                                         real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ], \
                                         real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ], \
                                         real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] ); \
  template \
  void \
  UpwindingHelper:: \
    computePhaseGravTerm< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                    localIndex const (&neighborIds)[ 3 ], \
                                    real64 const & transGravCoef, \
                                    ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                    ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                                    ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac, \
                                    ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                    real64 ( &phaseGravTerm )[ NP ][ NP-1 ], \
                                    real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ], \
                                    real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] ); \
  template \
  void \
  UpwindingHelper:: \
    computeUpwindedTotalMobility< NC, NP >( localIndex const (&localIds)[ 3 ], \
                                            localIndex const (&neighborIds)[ 3 ], \
                                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                            ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                            ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                                            real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                            real64 & totalMob,    \
                                            real64 ( &dTotalMob_dPres )[ 2 ], \
                                            real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] )

INST_UpwindingHelperNCNP( 1, 2 );
INST_UpwindingHelperNCNP( 2, 2 );
INST_UpwindingHelperNCNP( 3, 2 );
INST_UpwindingHelperNCNP( 4, 2 );
INST_UpwindingHelperNCNP( 5, 2 );

INST_UpwindingHelperNCNP( 1, 3 );
INST_UpwindingHelperNCNP( 2, 3 );
INST_UpwindingHelperNCNP( 3, 3 );
INST_UpwindingHelperNCNP( 4, 3 );
INST_UpwindingHelperNCNP( 5, 3 );

#undef INST_UpwindingHelperNCNP

#define INST_UpwindingHelperNP( NP ) \
  template \
  void \
  UpwindingHelper:: \
    setIndicesForTotalMobilityUpwinding< NP >( localIndex const (&localIds)[ 3 ], \
                                               localIndex const (&neighborIds)[ 3 ], \
                                               real64 const (&gravTerm)[ NP ][ NP-1 ], \
                                               localIndex ( &totalMobIds )[ NP ][ 3 ], \
                                               localIndex ( &totalMobPos )[ NP ] )

INST_UpwindingHelperNP( 2 );
INST_UpwindingHelperNP( 3 );

#undef INST_UpwindingHelperNP

/******************************** AssemblerKernelHelper ********************************/

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::
  applyGradient( arrayView1d< real64 const > const & facePres,
                 arrayView1d< real64 const > const & dFacePres,
                 arrayView1d< real64 const > const & faceGravCoef,
                 arraySlice1d< localIndex const > const & elemToFaces,
                 real64 const & elemPres,
                 real64 const & dElemPres,
                 real64 const & elemGravCoef,
                 arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & elemPhaseMassDens,
                 arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dElemPhaseMassDens_dPres,
                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dElemPhaseMassDens_dCompFrac,
                 arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & elemPhaseMob,
                 arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dElemPhaseMob_dPres,
                 arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dElemPhaseMob_dCompDens,
                 arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dElemCompFrac_dCompDens,
                 arraySlice2d< real64 const > const & transMatrix,
                 real64 ( & oneSidedVolFlux )[ NF ],
                 real64 ( & dOneSidedVolFlux_dPres )[ NF ],
                 real64 ( & dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                 real64 ( & dOneSidedVolFlux_dCompDens )[ NF ][ NC ] )
{
  real64 dPhaseMassDens_dC[ NP ][ NC ]{};
  real64 dPresDif_dCompDens[ NC ]{};
  real64 dPhaseGravDif_dCompDens[ NC ]{};
  real64 dPhaseMobPotDif_dCompDens[ NC ]{};

  // 0) precompute dPhaseDens_dC since it is always computed at the element center
  for( localIndex ip = 0; ip < NP; ++ip )
  {
    applyChainRule( NC,
                    dElemCompFrac_dCompDens,
                    dElemPhaseMassDens_dCompFrac[ip],
                    dPhaseMassDens_dC[ip] );
  }

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {

      // depth difference between element center and face center
      real64 const ccGravCoef = elemGravCoef;
      real64 const fGravCoef = faceGravCoef[elemToFaces[jfaceLoc]];
      real64 const gravCoefDif = ccGravCoef - fGravCoef;

      for( localIndex ip = 0; ip < NP; ++ip )
      {

        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres + dElemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

        // pressure difference
        real64 const presDif = ccPres - fPres;
        real64 const dPresDif_dPres = 1;
        real64 const dPresDif_dFacePres = -1;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPresDif_dCompDens[ic] = 0.0; // no capillary pressure
        }

        // gravity term
        real64 const phaseGravDif = elemPhaseMassDens[ip] * gravCoefDif;
        real64 const dPhaseGravDif_dPres = dElemPhaseMassDens_dPres[ip] * gravCoefDif;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseGravDif_dCompDens[ic] = dPhaseMassDens_dC[ip][ic] * gravCoefDif;
        }
        // no density evaluated at the face center

        // potential difference
        real64 const phasePotDif = presDif - phaseGravDif;
        real64 const phaseMobPotDif = elemPhaseMob[ip] * phasePotDif;
        real64 const dPhaseMobPotDif_dPres = dElemPhaseMob_dPres[ip] * phasePotDif
                                             + elemPhaseMob[ip] * (dPresDif_dPres - dPhaseGravDif_dPres);
        real64 const dPhaseMobPotDif_dFacePres = elemPhaseMob[ip] * dPresDif_dFacePres;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dPhaseMobPotDif_dCompDens[ic] = dElemPhaseMob_dCompDens[ip][ic] * phasePotDif
                                          + elemPhaseMob[ip] * (dPresDif_dCompDens[ic] - dPhaseGravDif_dCompDens[ic]);
        }

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        oneSidedVolFlux[ifaceLoc] = oneSidedVolFlux[ifaceLoc]
                                    + transMatrix[ifaceLoc][jfaceLoc] * phaseMobPotDif;
        dOneSidedVolFlux_dPres[ifaceLoc] = dOneSidedVolFlux_dPres[ifaceLoc]
                                           + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dPres;
        dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc]
                                                         + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dFacePres;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dOneSidedVolFlux_dCompDens[ifaceLoc][ic] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic]
                                                     + transMatrix[ifaceLoc][jfaceLoc] * dPhaseMobPotDif_dCompDens[ic];
        }
      }
    }
  }
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arrayView1d< real64 const > const & mimFaceGravCoef,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          real64 const & elemGravCoef,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
{
  localIndex constexpr NDOF = NC+1;

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF*(NF+1) ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + idof;
  }

  // divergence of fluxes
  real64 divMassFluxes[ NC ]{};
  real64 dDivMassFluxes_dElemVars[ NC ][ NDOF*(NF+1) ]{};
  real64 dDivMassFluxes_dFaceVars[ NC ][ NF ]{};

  // auxiliary variables for upwinding

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseViscCoef[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dPres[ NP ][ NC ]{};
  real64 dUpwPhaseViscCoef_dCompDens[ NP ][ NC ][ NC ]{};
  globalIndex upwViscDofNumber = 0;

  // gravity term: ( \rho_l - \rho_m ) g \Delta z
  real64 phaseGravTerm[ NP ][ NP-1 ]{};
  real64 dPhaseGravTerm_dPres[ NP ][ NP-1 ][ 2 ]{};
  real64 dPhaseGravTerm_dCompDens[ NP ][ NP-1 ][ 2 ][ NC ]{};

  // upwinding phase buoyancy transport coefficients
  real64 upwPhaseGravCoef[ NP ][ NP-1 ][ NC ]{};
  real64 dUpwPhaseGravCoef_dPres[ NP ][ NP-1 ][ NC ][ 2 ]{};
  real64 dUpwPhaseGravCoef_dCompDens[ NP ][ NP-1 ][ NC ][ 2 ][ NC ]{};

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {

    // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element

    localIndex neighborIds[ 3 ] = { localIds[0], localIds[1], localIds[2] };
    hybridFVMKernels::CellConnectivity::findNeighbor( localIds,
                                                      ifaceLoc,
                                                      elemRegionList,
                                                      elemSubRegionList,
                                                      elemList,
                                                      regionFilter,
                                                      elemToFaces,
                                                      neighborIds );
    localIndex const neighborDofNumber = elemDofNumber[neighborIds[0]][neighborIds[1]][neighborIds[2]];

    // 2) *************** Assemble viscous terms ******************

    // 2.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} for each phase at this face
    UpwindingHelper::upwindViscousCoefficient< NC, NP >( localIds,
                                                         neighborIds,
                                                         phaseDens,
                                                         dPhaseDens_dPres,
                                                         dPhaseDens_dCompFrac,
                                                         phaseMob,
                                                         dPhaseMob_dPres,
                                                         dPhaseMob_dCompDens,
                                                         dCompFrac_dCompDens,
                                                         phaseCompFrac,
                                                         dPhaseCompFrac_dPres,
                                                         dPhaseCompFrac_dCompFrac,
                                                         elemDofNumber,
                                                         oneSidedVolFlux[ifaceLoc],
                                                         upwPhaseViscCoef,
                                                         dUpwPhaseViscCoef_dPres,
                                                         dUpwPhaseViscCoef_dCompDens,
                                                         upwViscDofNumber );

    // 2.b) Add the \x_{c,\ell} \rho_{\ell} \frac{\lambda_{\ell}}{\lambda_T} q_T of this face to the divergence of the flux in this cell
    assembleViscousFlux< NF, NC, NP >( ifaceLoc,
                                       oneSidedVolFlux,
                                       dOneSidedVolFlux_dPres,
                                       dOneSidedVolFlux_dFacePres,
                                       dOneSidedVolFlux_dCompDens,
                                       upwPhaseViscCoef,
                                       dUpwPhaseViscCoef_dPres,
                                       dUpwPhaseViscCoef_dCompDens,
                                       elemDofNumber[localIds[0]][localIds[1]][localIds[2]],
                                       neighborDofNumber,
                                       upwViscDofNumber,
                                       faceDofNumber[elemToFaces[ifaceLoc]],
                                       dt,
                                       divMassFluxes,
                                       dDivMassFluxes_dElemVars,
                                       dDivMassFluxes_dFaceVars,
                                       dofColIndicesElemVars,
                                       dofColIndicesFaceVars );

    // 3) *************** Assemble buoyancy terms ******************

    real64 const transGravCoef = (localIds[0] != neighborIds[0] || localIds[1] != neighborIds[1] || localIds[2] != neighborIds[2])
                                 * transMatrixGrav[ifaceLoc][ifaceLoc] * (elemGravCoef - mimFaceGravCoef[elemToFaces[ifaceLoc]]);

    // 3.a) Compute the upwinded x_{c, \ell} \rho_{\ell} \frac{\lambda_{\ell}\lambda_m}{\lambda_T}
    //      and (\rho_{\ell} - \rho_m) g \Delta z for each phase at this face
    UpwindingHelper::upwindBuoyancyCoefficient< NC, NP >( localIds,
                                                          neighborIds,
                                                          transGravCoef,
                                                          phaseDens,
                                                          dPhaseDens_dPres,
                                                          dPhaseDens_dCompFrac,
                                                          phaseMassDens,
                                                          dPhaseMassDens_dPres,
                                                          dPhaseMassDens_dCompFrac,
                                                          phaseMob,
                                                          dPhaseMob_dPres,
                                                          dPhaseMob_dCompDens,
                                                          dCompFrac_dCompDens,
                                                          phaseCompFrac,
                                                          dPhaseCompFrac_dPres,
                                                          dPhaseCompFrac_dCompFrac,
                                                          phaseGravTerm,
                                                          dPhaseGravTerm_dPres,
                                                          dPhaseGravTerm_dCompDens,
                                                          upwPhaseGravCoef,
                                                          dUpwPhaseGravCoef_dPres,
                                                          dUpwPhaseGravCoef_dCompDens );

    // 3.b) Add the buoyancy term of this face to the divergence of the flux in this cell
    assembleBuoyancyFlux< NF, NC, NP >( ifaceLoc,
                                        phaseGravTerm,
                                        dPhaseGravTerm_dPres,
                                        dPhaseGravTerm_dCompDens,
                                        upwPhaseGravCoef,
                                        dUpwPhaseGravCoef_dPres,
                                        dUpwPhaseGravCoef_dCompDens,
                                        dt,
                                        divMassFluxes,
                                        dDivMassFluxes_dElemVars );

  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  for( localIndex ic = 0; ic < NC; ++ic )
  {
    localIndex const eqnRowLocalIndex =
      LvArray::integerConversion< localIndex >( elemDofNumber[localIds[0]][localIds[1]][localIds[2]] + ic - rankOffset );

    GEOSX_ASSERT_GE( eqnRowLocalIndex, 0 );
    GEOSX_ASSERT_GT( localMatrix.numRows(), eqnRowLocalIndex );

    // residual
    localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + divMassFluxes[ic];

    // jacobian -- derivative wrt elem centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesElemVars[0],
                                                              &dDivMassFluxes_dElemVars[0][0] + ic * NDOF * (NF+1),
                                                              NDOF * (NF+1) );

    // jacobian -- derivatives wrt face centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &dofColIndicesFaceVars[0],
                                                              &dDivMassFluxes_dFaceVars[0][0] + ic * NF,
                                                              NF );
  }
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::
  assembleViscousFlux( localIndex const ifaceLoc,
                       real64 const (&oneSidedVolFlux)[ NF ],
                       real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ],
                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ],
                       globalIndex const elemDofNumber,
                       globalIndex const neighborDofNumber,
                       globalIndex const upwViscDofNumber,
                       globalIndex const faceDofNumber,
                       real64 const & dt,
                       real64 ( & divMassFluxes )[ NC ],
                       real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                       real64 ( & dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                       globalIndex ( & dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                       globalIndex ( & dofColIndicesFaceVars )[ NF ] )
{
  localIndex constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      // compute the mass flux at the one-sided face plus its derivatives
      // add the newly computed flux to the sum

      real64 const dt_upwPhaseViscCoef = dt * upwPhaseViscCoef[ip][ic];

      // residual
      divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseViscCoef * oneSidedVolFlux[ifaceLoc];

      // local derivatives
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + dt_upwPhaseViscCoef * dOneSidedVolFlux_dPres[ifaceLoc];
      dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                        + ( elemDofNumber == upwViscDofNumber )
                                        * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];
      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + dt_upwPhaseViscCoef * dOneSidedVolFlux_dCompDens[ifaceLoc][jc];
        dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                             + ( elemDofNumber == upwViscDofNumber )
                                             * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      // neighbor derivatives
      dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                     + ( elemDofNumber != upwViscDofNumber )
                                                     * dt * dUpwPhaseViscCoef_dPres[ip][ic] * oneSidedVolFlux[ifaceLoc];

      for( localIndex jc = 0; jc < NC; ++jc )
      {
        dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                            + ( elemDofNumber != upwViscDofNumber )
                                                            * dt * dUpwPhaseViscCoef_dCompDens[ip][ic][jc] * oneSidedVolFlux[ifaceLoc];
      }

      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dDivMassFluxes_dFaceVars[ic][jfaceLoc] = dDivMassFluxes_dFaceVars[ic][jfaceLoc]
                                                 + dt_upwPhaseViscCoef * dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      }
    }
  }

  // collect the relevant dof numbers
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[elemVarsOffset+idof] = neighborDofNumber + idof;
  }
  dofColIndicesFaceVars[ifaceLoc] = faceDofNumber;
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::
  assembleBuoyancyFlux( localIndex const ifaceLoc,
                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ],
                        real64 const & dt,
                        real64 ( & divMassFluxes )[ NC ],
                        real64 ( & dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] )
{
  localIndex constexpr NDOF = NC+1;
  localIndex const elemVarsOffset = NDOF*(ifaceLoc+1);

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex jp = 0; jp < NP - 1; ++jp )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        real64 const dt_upwPhaseGravCoef = dt * upwPhaseGravCoef[ip][jp][ic];

        // residual
        divMassFluxes[ic] = divMassFluxes[ic] + dt_upwPhaseGravCoef * phaseGravTerm[ip][jp];

        // local derivatives
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::LOCAL] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][0] = dDivMassFluxes_dElemVars[ic][0]
                                          + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::LOCAL];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::LOCAL][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][jc+1] = dDivMassFluxes_dElemVars[ic][jc+1]
                                               + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::LOCAL][jc];
        }

        // neighbor derivatives
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt * dUpwPhaseGravCoef_dPres[ip][jp][ic][Pos::NEIGHBOR] * phaseGravTerm[ip][jp];
        dDivMassFluxes_dElemVars[ic][elemVarsOffset] = dDivMassFluxes_dElemVars[ic][elemVarsOffset]
                                                       + dt_upwPhaseGravCoef * dPhaseGravTerm_dPres[ip][jp][Pos::NEIGHBOR];

        for( localIndex jc = 0; jc < NC; ++jc )
        {
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt * dUpwPhaseGravCoef_dCompDens[ip][jp][ic][Pos::NEIGHBOR][jc] * phaseGravTerm[ip][jp];
          dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1] = dDivMassFluxes_dElemVars[ic][elemVarsOffset+jc+1]
                                                              + dt_upwPhaseGravCoef * dPhaseGravTerm_dCompDens[ip][jp][Pos::NEIGHBOR][jc];
        }
      }
    }
  }
}

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernelHelper::
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[ NF ],
                           real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
{
  localIndex constexpr NDOF = NC+1;

  // fluxes
  real64 dFlux_dElemVars[ NDOF ]{};
  real64 dFlux_dFaceVars[ NF ]{};

  // dof numbers
  globalIndex dofColIndicesElemVars[ NDOF ]{};
  globalIndex dofColIndicesFaceVars[ NF ]{};
  for( localIndex idof = 0; idof < NDOF; ++idof )
  {
    dofColIndicesElemVars[idof] = elemDofNumber + idof;
  }

  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
    {
      continue;
    }

    // flux at this face
    real64 const flux = oneSidedVolFlux[ifaceLoc];
    dFlux_dElemVars[0] = dOneSidedVolFlux_dPres[ifaceLoc];
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dFlux_dElemVars[ic+1] = dOneSidedVolFlux_dCompDens[ifaceLoc][ic];
    }

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dFaceVars[jfaceLoc] = dOneSidedVolFlux_dFacePres[ifaceLoc][jfaceLoc];
      dofColIndicesFaceVars[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

    GEOSX_ASSERT_GE( eqnLocalRowIndex, 0 );
    GEOSX_ASSERT_GT( localMatrix.numRows(), eqnLocalRowIndex );

    // residual
    atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivatives wrt elem-centered terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesElemVars[0],
                                                                      &dFlux_dElemVars[0],
                                                                      NDOF );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                      &dofColIndicesFaceVars[0],
                                                                      &dFlux_dFaceVars[0],
                                                                      NF );
  }
}


#define INST_AssemblerKernelHelper( NF, NC, NP ) \
  template \
  void \
  AssemblerKernelHelper:: \
    applyGradient< NF, NC, NP >( arrayView1d< real64 const > const & facePres, \
                                 arrayView1d< real64 const > const & dFacePres, \
                                 arrayView1d< real64 const > const & faceGravCoef, \
                                 arraySlice1d< localIndex const > const & elemToFaces, \
                                 real64 const & elemPres, \
                                 real64 const & dElemPres, \
                                 real64 const & elemGravCoef, \
                                 arraySlice1d< real64 const, multifluid::USD_PHASE-2 > const & elemPhaseMassDens, \
                                 arraySlice1d< real64 const, multifluid::USD_PHASE-2 > const & dElemPhaseMassDens_dPres, \
                                 arraySlice2d< real64 const, multifluid::USD_PHASE_DC-2 > const & dElemPhaseMassDens_dCompFrac, \
                                 arraySlice1d< real64 const, compflow::USD_PHASE-1 > const & elemPhaseMob, \
                                 arraySlice1d< real64 const, compflow::USD_PHASE-1 > const & dElemPhaseMob_dPres, \
                                 arraySlice2d< real64 const, compflow::USD_PHASE_DC-1 > const & dElemPhaseMob_dCompDens, \
                                 arraySlice2d< real64 const, compflow::USD_COMP_DC-1 > const & dElemCompFrac_dCompDens, \
                                 arraySlice2d< real64 const > const & transMatrix, \
                                 real64 ( &oneSidedVolFlux )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dPres )[ NF ], \
                                 real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ], \
                                 real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] ); \
  template \
  void \
  AssemblerKernelHelper:: \
    assembleFluxDivergence< NF, NC, NP >( localIndex const (&localIds)[ 3 ], \
                                          globalIndex const rankOffset, \
                                          arrayView2d< localIndex const > const & elemRegionList, \
                                          arrayView2d< localIndex const > const & elemSubRegionList, \
                                          arrayView2d< localIndex const > const & elemList, \
                                          SortedArrayView< localIndex const > const & regionFilter, \
                                          arrayView1d< globalIndex const > const & faceDofNumber, \
                                          arrayView1d< real64 const > const & mimFaceGravCoef, \
                                          arraySlice1d< localIndex const > const & elemToFaces, \
                                          real64 const & elemGravCoef, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac, \
                                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                                          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                                          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac, \
                                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                          arraySlice2d< real64 const > const & transMatrixGrav, \
                                          real64 const (&oneSidedVolFlux)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                          real64 const & dt, \
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                          arrayView1d< real64 > const & localRhs ); \
  template \
  void \
  AssemblerKernelHelper:: \
    assembleViscousFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                       real64 const (&oneSidedVolFlux)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                       real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                       real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                       real64 const (&upwPhaseViscCoef)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dPres)[ NP ][ NC ], \
                                       real64 const (&dUpwPhaseViscCoef_dCompDens)[ NP ][ NC ][ NC ], \
                                       globalIndex const elemDofNumber, \
                                       globalIndex const neighborDofNumber, \
                                       globalIndex const upwViscDofNumber, \
                                       globalIndex const faceDofNumber, \
                                       real64 const & dt, \
                                       real64 ( &divMassFluxes )[ NC ], \
                                       real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ], \
                                       real64 ( &dDivMassFluxes_dFaceVars )[ NC ][ NF ], \
                                       globalIndex ( &dofColIndicesElemVars )[ (NC+1)*(NF+1) ], \
                                       globalIndex ( &dofColIndicesFaceVars )[ NF ] ); \
  template \
  void \
  AssemblerKernelHelper:: \
    assembleBuoyancyFlux< NF, NC, NP >( localIndex const ifaceLoc, \
                                        real64 const (&phaseGravTerm)[ NP ][ NP-1 ], \
                                        real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ], \
                                        real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ], \
                                        real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ], \
                                        real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ], \
                                        real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ], \
                                        real64 const & dt, \
                                        real64 ( &divMassFluxes )[ NC ], \
                                        real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] ); \
  template \
  void \
  AssemblerKernelHelper:: \
    assembleFaceConstraints< NF, NC, NP >( arrayView1d< globalIndex const > const & faceDofNumber, \
                                           arrayView1d< integer const > const & faceGhostRank, \
                                           arraySlice1d< localIndex const > const & elemToFaces, \
                                           globalIndex const elemDofNumber, \
                                           globalIndex const rankOffset,  \
                                           real64 const (&oneSidedVolFlux)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dPres)[ NF ], \
                                           real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ], \
                                           real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ], \
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernelHelper( 4, 1, 2 );
INST_AssemblerKernelHelper( 4, 2, 2 );
INST_AssemblerKernelHelper( 4, 3, 2 );
INST_AssemblerKernelHelper( 4, 4, 2 );
INST_AssemblerKernelHelper( 4, 5, 2 );

INST_AssemblerKernelHelper( 4, 1, 3 );
INST_AssemblerKernelHelper( 4, 2, 3 );
INST_AssemblerKernelHelper( 4, 3, 3 );
INST_AssemblerKernelHelper( 4, 4, 3 );
INST_AssemblerKernelHelper( 4, 5, 3 );

INST_AssemblerKernelHelper( 5, 1, 2 );
INST_AssemblerKernelHelper( 5, 2, 2 );
INST_AssemblerKernelHelper( 5, 3, 2 );
INST_AssemblerKernelHelper( 5, 4, 2 );
INST_AssemblerKernelHelper( 5, 5, 2 );

INST_AssemblerKernelHelper( 5, 1, 3 );
INST_AssemblerKernelHelper( 5, 2, 3 );
INST_AssemblerKernelHelper( 5, 3, 3 );
INST_AssemblerKernelHelper( 5, 4, 3 );
INST_AssemblerKernelHelper( 5, 5, 3 );

INST_AssemblerKernelHelper( 6, 1, 2 );
INST_AssemblerKernelHelper( 6, 2, 2 );
INST_AssemblerKernelHelper( 6, 3, 2 );
INST_AssemblerKernelHelper( 6, 4, 2 );
INST_AssemblerKernelHelper( 6, 5, 2 );

INST_AssemblerKernelHelper( 6, 1, 3 );
INST_AssemblerKernelHelper( 6, 2, 3 );
INST_AssemblerKernelHelper( 6, 3, 3 );
INST_AssemblerKernelHelper( 6, 4, 3 );
INST_AssemblerKernelHelper( 6, 5, 3 );

#undef INST_AssemblerKernelHelper

/******************************** AssemblerKernel ********************************/

template< localIndex NF, localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
AssemblerKernel::
  compute( localIndex const er, localIndex const esr, localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arrayView1d< real64 const > const & mimFaceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  real64 oneSidedVolFlux[ NF ]{};
  real64 dOneSidedVolFlux_dPres[ NF ]{};
  real64 dOneSidedVolFlux_dFacePres[ NF ][ NF ]{};
  real64 dOneSidedVolFlux_dCompDens[ NF ][ NC ]{};

  localIndex const localIds[3] = { er, esr, ei };

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::applyGradient< NF, NC, NP >( facePres,
                                                      dFacePres,
                                                      faceGravCoef,
                                                      elemToFaces,
                                                      elemPres,
                                                      dElemPres,
                                                      elemGravCoef,
                                                      phaseMassDens[er][esr][ei][0],
                                                      dPhaseMassDens_dPres[er][esr][ei][0],
                                                      dPhaseMassDens_dCompFrac[er][esr][ei][0],
                                                      phaseMob[er][esr][ei],
                                                      dPhaseMob_dPres[er][esr][ei],
                                                      dPhaseMob_dCompDens[er][esr][ei],
                                                      dCompFrac_dCompDens[er][esr][ei],
                                                      transMatrix,
                                                      oneSidedVolFlux,
                                                      dOneSidedVolFlux_dPres,
                                                      dOneSidedVolFlux_dFacePres,
                                                      dOneSidedVolFlux_dCompDens );

  // at this point, we know the local flow direction in the element
  // so we can upwind the transport coefficients (mobilities) at the one sided faces
  // ** this function needs non-local information **
  if( elemGhostRank < 0 )
  {
    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::assembleFluxDivergence< NF, NC, NP >( localIds,
                                                                 rankOffset,
                                                                 elemRegionList,
                                                                 elemSubRegionList,
                                                                 elemList,
                                                                 regionFilter,
                                                                 faceDofNumber,
                                                                 mimFaceGravCoef,
                                                                 elemToFaces,
                                                                 elemGravCoef,
                                                                 phaseDens,
                                                                 dPhaseDens_dPres,
                                                                 dPhaseDens_dCompFrac,
                                                                 phaseMassDens,
                                                                 dPhaseMassDens_dPres,
                                                                 dPhaseMassDens_dCompFrac,
                                                                 phaseMob,
                                                                 dPhaseMob_dPres,
                                                                 dPhaseMob_dCompDens,
                                                                 dCompFrac_dCompDens,
                                                                 phaseCompFrac,
                                                                 dPhaseCompFrac_dPres,
                                                                 dPhaseCompFrac_dCompFrac,
                                                                 elemDofNumber,
                                                                 transMatrixGrav,
                                                                 oneSidedVolFlux,
                                                                 dOneSidedVolFlux_dPres,
                                                                 dOneSidedVolFlux_dFacePres,
                                                                 dOneSidedVolFlux_dCompDens,
                                                                 dt,
                                                                 localMatrix,
                                                                 localRhs );
  }

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::assembleFaceConstraints< NF, NC, NP >( faceDofNumber,
                                                                faceGhostRank,
                                                                elemToFaces,
                                                                elemDofNumber[er][esr][ei],
                                                                rankOffset,
                                                                oneSidedVolFlux,
                                                                dOneSidedVolFlux_dPres,
                                                                dOneSidedVolFlux_dFacePres,
                                                                dOneSidedVolFlux_dCompDens,
                                                                localMatrix,
                                                                localRhs );

}

#define INST_AssemblerKernel( NF, NC, NP ) \
  template \
  void \
  AssemblerKernel:: \
    compute< NF, NC, NP >( localIndex const er, localIndex const esr, localIndex const ei, \
                           SortedArrayView< localIndex const > const & regionFilter, \
                           arrayView2d< localIndex const > const & elemRegionList, \
                           arrayView2d< localIndex const > const & elemSubRegionList, \
                           arrayView2d< localIndex const > const & elemList, \
                           arrayView1d< globalIndex const > const & faceDofNumber, \
                           arrayView1d< integer const > const & faceGhostRank, \
                           arrayView1d< real64 const > const & facePres, \
                           arrayView1d< real64 const > const & dFacePres, \
                           arrayView1d< real64 const > const & faceGravCoef, \
                           arrayView1d< real64 const > const & mimFaceGravCoef, \
                           arraySlice1d< localIndex const > const & elemToFaces, \
                           real64 const & elemPres, \
                           real64 const & dElemPres, \
                           real64 const & elemGravCoef, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                           ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac, \
                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                           ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                           ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                           ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                           ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac, \
                           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                           integer const elemGhostRank, \
                           globalIndex const rankOffset, \
                           real64 const & dt, \
                           arraySlice2d< real64 const > const & transMatrix, \
                           arraySlice2d< real64 const > const & transMatrixGrav, \
                           CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                           arrayView1d< real64 > const & localRhs )

INST_AssemblerKernel( 4, 1, 2 );
INST_AssemblerKernel( 4, 2, 2 );
INST_AssemblerKernel( 4, 3, 2 );
INST_AssemblerKernel( 4, 4, 2 );
INST_AssemblerKernel( 4, 5, 2 );

INST_AssemblerKernel( 4, 1, 3 );
INST_AssemblerKernel( 4, 2, 3 );
INST_AssemblerKernel( 4, 3, 3 );
INST_AssemblerKernel( 4, 4, 3 );
INST_AssemblerKernel( 4, 5, 3 );

INST_AssemblerKernel( 5, 1, 2 );
INST_AssemblerKernel( 5, 2, 2 );
INST_AssemblerKernel( 5, 3, 2 );
INST_AssemblerKernel( 5, 4, 2 );
INST_AssemblerKernel( 5, 5, 2 );

INST_AssemblerKernel( 5, 1, 3 );
INST_AssemblerKernel( 5, 2, 3 );
INST_AssemblerKernel( 5, 3, 3 );
INST_AssemblerKernel( 5, 4, 3 );
INST_AssemblerKernel( 5, 5, 3 );

INST_AssemblerKernel( 6, 1, 2 );
INST_AssemblerKernel( 6, 2, 2 );
INST_AssemblerKernel( 6, 3, 2 );
INST_AssemblerKernel( 6, 4, 2 );
INST_AssemblerKernel( 6, 5, 2 );

INST_AssemblerKernel( 6, 1, 3 );
INST_AssemblerKernel( 6, 2, 3 );
INST_AssemblerKernel( 6, 3, 3 );
INST_AssemblerKernel( 6, 4, 3 );
INST_AssemblerKernel( 6, 5, 3 );

#undef INST_AssemblerKernel



/******************************** FluxKernel ********************************/

template< localIndex NF, localIndex NC, localIndex NP, typename IP_TYPE >
void
FluxKernel::
  launch( localIndex er, localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::PermeabilityBase const & permeabilityModel,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView1d< globalIndex const > const & faceDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          arrayView1d< real64 const > const & faceGravCoef,
          arrayView1d< real64 const > const & mimFaceGravCoef,
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens,
          ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
{
  // get the cell-centered DOF numbers and ghost rank for the assembly
  arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres  =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & dElemPres =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::deltaPressureString() );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >( CellBlock::viewKeyStruct::elementCenterString() );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >( CellBlock::viewKeyStruct::elementVolumeString() );

  arrayView2d< real64 const > const & elemPerm = permeabilityModel.permeability();

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >( CompositionalMultiphaseBase::viewKeyStruct::gravityCoefString() );

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
  {

    // transmissibility matrix
    stackArray2d< real64, NF *NF > transMatrix( NF, NF );
    stackArray2d< real64, NF *NF > transMatrixGrav( NF, NF );

    real64 const perm[ 3 ] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    IP_TYPE::template compute< NF >( nodePosition,
                                     transMultiplier,
                                     faceToNodes,
                                     elemToFaces[ei],
                                     elemCenter[ei],
                                     elemVolume[ei],
                                     perm,
                                     lengthTolerance,
                                     transMatrix );

    // currently the gravity term in the transport scheme is treated as in MRST, that is, always with TPFA
    // this is why below we have to recompute the TPFA transmissibility in addition to the transmissibility matrix above
    // TODO: treat the gravity term with a consistent inner product
    mimeticInnerProduct::TPFAInnerProduct::compute< NF >( nodePosition,
                                                          transMultiplier,
                                                          faceToNodes,
                                                          elemToFaces[ei],
                                                          elemCenter[ei],
                                                          elemVolume[ei],
                                                          perm,
                                                          lengthTolerance,
                                                          transMatrixGrav );

    // perform flux assembly in this element
    CompositionalMultiphaseHybridFVMKernels::AssemblerKernel::compute< NF, NC, NP >( er, esr, ei,
                                                                                     regionFilter,
                                                                                     elemRegionList,
                                                                                     elemSubRegionList,
                                                                                     elemList,
                                                                                     faceDofNumber,
                                                                                     faceGhostRank,
                                                                                     facePres,
                                                                                     dFacePres,
                                                                                     faceGravCoef,
                                                                                     mimFaceGravCoef,
                                                                                     elemToFaces[ei],
                                                                                     elemPres[ei],
                                                                                     dElemPres[ei],
                                                                                     elemGravCoef[ei],
                                                                                     phaseDens,
                                                                                     dPhaseDens_dPres,
                                                                                     dPhaseDens_dCompFrac,
                                                                                     phaseMassDens,
                                                                                     dPhaseMassDens_dPres,
                                                                                     dPhaseMassDens_dCompFrac,
                                                                                     phaseMob,
                                                                                     dPhaseMob_dPres,
                                                                                     dPhaseMob_dCompDens,
                                                                                     dCompFrac_dCompDens,
                                                                                     phaseCompFrac,
                                                                                     dPhaseCompFrac_dPres,
                                                                                     dPhaseCompFrac_dCompFrac,
                                                                                     elemDofNumber,
                                                                                     elemGhostRank[ei],
                                                                                     rankOffset,
                                                                                     dt,
                                                                                     transMatrix,
                                                                                     transMatrixGrav,
                                                                                     localMatrix,
                                                                                     localRhs );
  } );
}

#define INST_FluxKernel( NF, NC, NP, IP_TYPE ) \
  template \
  void \
  FluxKernel:: \
    launch< NF, NC, NP, IP_TYPE >( localIndex er, localIndex esr, \
                                   CellElementSubRegion const & subRegion, \
                                   constitutive::PermeabilityBase const & permeabilityModel, \
                                   SortedArrayView< localIndex const > const & regionFilter, \
                                   arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
                                   arrayView2d< localIndex const > const & elemRegionList, \
                                   arrayView2d< localIndex const > const & elemSubRegionList, \
                                   arrayView2d< localIndex const > const & elemList, \
                                   ArrayOfArraysView< localIndex const > const & faceToNodes, \
                                   arrayView1d< globalIndex const > const & faceDofNumber, \
                                   arrayView1d< integer const > const & faceGhostRank, \
                                   arrayView1d< real64 const > const & facePres, \
                                   arrayView1d< real64 const > const & dFacePres, \
                                   arrayView1d< real64 const > const & faceGravCoef, \
                                   arrayView1d< real64 const > const & mimFaceGravCoef, \
                                   arrayView1d< real64 const > const & transMultiplier, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseDens_dPres, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseDens_dCompFrac, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & phaseMassDens, \
                                   ElementViewConst< arrayView3d< real64 const, multifluid::USD_PHASE > > const & dPhaseMassDens_dPres, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_DC > > const & dPhaseMassDens_dCompFrac, \
                                   ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & phaseMob, \
                                   ElementViewConst< arrayView2d< real64 const, compflow::USD_PHASE > > const & dPhaseMob_dPres, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_PHASE_DC > > const & dPhaseMob_dCompDens, \
                                   ElementViewConst< arrayView3d< real64 const, compflow::USD_COMP_DC > > const & dCompFrac_dCompDens, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & phaseCompFrac, \
                                   ElementViewConst< arrayView4d< real64 const, multifluid::USD_PHASE_COMP > > const & dPhaseCompFrac_dPres, \
                                   ElementViewConst< arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > > const & dPhaseCompFrac_dCompFrac, \
                                   ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber, \
                                   globalIndex const rankOffset, \
                                   real64 const lengthTolerance, \
                                   real64 const dt, \
                                   CRSMatrixView< real64, globalIndex const > const & localMatrix, \
                                   arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::TPFAInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::TPFAInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::TPFAInnerProduct const );


INST_FluxKernel( 4, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 4, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 4, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 5, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 5, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 2, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 2, mimeticInnerProduct::BdVLMInnerProduct const );

INST_FluxKernel( 6, 1, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 2, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 3, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 4, 3, mimeticInnerProduct::BdVLMInnerProduct const );
INST_FluxKernel( 6, 5, 3, mimeticInnerProduct::BdVLMInnerProduct const );


#undef INST_FluxKernel


/******************************** PhaseMobilityKernel ********************************/

template< localIndex NC, localIndex NP >
GEOSX_HOST_DEVICE
void
PhaseMobilityKernel::
  compute( arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseVisc,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const, relperm::USD_RELPERM - 2 > const & phaseRelPerm,
           arraySlice2d< real64 const, relperm::USD_RELPERM_DS - 2 > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseMob,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseMob_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseMob_dComp )
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

    real64 const mobility = relPerm / viscosity;

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
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( size, [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    compute< NC, NP >( dCompFrac_dCompDens[a],
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
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp,
          arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm,
          arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseMob,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )
{
  forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const i )
  {
    localIndex const a = targetSet[ i ];
    compute< NC, NP >( dCompFrac_dCompDens[a],
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
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm, \
                      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseMob, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp ); \
  template \
  void \
  PhaseMobilityKernel:: \
    launch< NC, NP >( SortedArrayView< localIndex const > const & targetSet, \
                      arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseVisc, \
                      arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseVisc_dPres, \
                      arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseVisc_dComp, \
                      arrayView3d< real64 const, relperm::USD_RELPERM > const & phaseRelPerm, \
                      arrayView4d< real64 const, relperm::USD_RELPERM_DS > const & dPhaseRelPerm_dPhaseVolFrac, \
                      arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres, \
                      arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp, \
                      arrayView2d< real64, compflow::USD_PHASE > const & phaseMob, \
                      arrayView2d< real64, compflow::USD_PHASE > const & dPhaseMob_dPres, \
                      arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseMob_dComp )

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


} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx
