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
 * @file CompositionalMultiphaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"

namespace geosx
{

namespace CompositionalMultiphaseHybridFVMKernels
{



// struct to specify local and neighbor derivatives
struct Pos
{
  static constexpr integer LOCAL = 0;
  static constexpr integer NEIGHBOR = 1;
};

/******************************** UpwindingHelper ********************************/

struct UpwindingHelper
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

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] oneSidedVolFlux the total volumetric flux at this face
   * @param[out] upwPhaseViscCoef the upwind viscous transport coef at this face
   * @param[out] dUpwPhaseViscCoef_dPres the derivative of the upwind viscous transport coef wrt pressure at this face
   * @param[out] dUpwPhaseViscCoef_dCompDens the derivatives of the upwind viscous transport coef wrt component density at this face
   * @param[out] upwViscDofNumber the dof number of the upwind cell at this face
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    upwindViscousCoefficient( localIndex const (&localIds)[ 3 ],
                              localIndex const (&neighborIds)[ 3 ],
                              ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                              ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                              ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                              ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                              ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                              ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                              ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                              ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                              ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                              ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                              real64 const & oneSidedVolFlux,
                              real64 ( &upwPhaseViscCoef )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dPres )[ NP ][ NC ],
                              real64 ( &dUpwPhaseViscCoef_dCompDens )[ NP ][ NC ][ NC ],
                              globalIndex & upwViscDofNumber );

  /**
   * @brief At a given one-sided face, compute the upwind viscous transport coefficient
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[inout] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[inout] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[inout] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    upwindBuoyancyCoefficient( localIndex const (&localIds)[ 3 ],
                               localIndex const (&neighborIds)[ 3 ],
                               real64 const & transGravCoef,
                               ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                               ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                               ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                               ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                               ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                               ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
                               ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                               ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                               ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                               ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                               ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                               ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                               ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                               real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                               real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                               real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ],
                               real64 ( &upwPhaseGravCoef )[ NP ][ NP-1 ][ NC ],
                               real64 ( &dUpwPhaseGravCoef_dPres )[ NP ][ NP-1 ][ NC ][ 2 ],
                               real64 ( &dUpwPhaseGravCoef_dCompDens )[ NP ][ NP-1 ][ NC ][ 2 ][ NC ] );


  /**
   * @brief At a given one-sided face, compute the gravCoef multiplied by the difference in phase densities
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] transGravCoef
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[inout] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[inout] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    computePhaseGravTerm( localIndex const (&localIds)[ 3 ],
                          localIndex const (&neighborIds)[ 3 ],
                          real64 const & transGravCoef,
                          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
                          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                          real64 ( &phaseGravTerm )[ NP ][ NP-1 ],
                          real64 ( &dPhaseGravTerm_dPres )[ NP ][ NP-1 ][ 2 ],
                          real64 ( &dPhaseGravTerm_dCompDens )[ NP ][ NP-1 ][ 2 ][ NC ] );

  /**
   * @brief At a given one-sided face, compute the upwinded total mobility
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[inout] totalMob the upwinded total mobility
   * @param[inout] dTotalMob_dPres the derivative of the upwinded total mobility wrt pressure
   * @param[inout] dTotalMob_dCompDens the derivative of the upwinded total mobility wrt component density
   */
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    computeUpwindedTotalMobility( localIndex const (&localIds)[ 3 ],
                                  localIndex const (&neighborIds)[ 3 ],
                                  ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                                  ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                                  ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                                  real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                                  real64 & totalMob,
                                  real64 ( &dTotalMob_dPres )[ 2 ],
                                  real64 ( &dTotalMob_dCompDens )[ 2 ][ NC ] );

  /**
   * @brief Set the element indices used to evaluate the mobility ratios of the buoyancy term in hybrid upwinding
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] neighborIds region, subRegion, and element indices of the neigbhbor element
   * @param[in] gravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] totalMob the upwinded total mobility
   * @param[in] eru region index of the upwind element
   * @param[in] esru subRegion index of the upwind element
   * @param[in] eiu element index of the upwind element
   * @param[in] posu position (local or neighbor) of the upwind element
   * @param[in] erd region index of the downwind element
   * @param[in] esrd subRegion index of the downwind element
   * @param[in] eid element index of the downwind element
   * @param[in] posd position (local or neighbor) of the downwind element
   */
  GEOSX_HOST_DEVICE
  static void
  setIndicesForMobilityRatioUpwinding( localIndex const (&localIds)[ 3 ],
                                       localIndex const (&neighborIds)[ 3 ],
                                       real64 const & gravTerm,
                                       localIndex & eru, localIndex & esru, localIndex & eiu, localIndex & posu,
                                       localIndex & erd, localIndex & esrd, localIndex & eid, localIndex & posd );

  /**
   * @brief Set the element indices used to evaluate the total mobility of the buoyancy term in hybrid upwinding
   * @param[in] localIds triplet of indices for the local element
   * @param[in] neighborIds triplet of indices for the neighbor element
   * @param[in] gravTerm gravity term used to upwind
   * @param[out] totalMobIds for each phase, triplet of indices of the upwind element
   * @param[out] totalMobPos for each phase, flag specifying with the upwind element is local or neighbor
   */
  template< localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    setIndicesForTotalMobilityUpwinding( localIndex const (&localIds)[ 3 ],
                                         localIndex const (&neighborIds)[ 3 ],
                                         real64 const (&gravTerm)[ NP ][ NP-1 ],
                                         localIndex ( &totalMobIds )[ NP ][ 3 ],
                                         localIndex ( &totalMobPos )[ NP ] );

};


/******************************** AssemblerKernelHelper ********************************/

struct AssemblerKernelHelper
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

  /**
   * @brief In a given element, compute the transmissibility-weighted pressure gradients in the cell
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the element
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the element wrt pressure
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the element wrt component fraction
   * @param[in] phaseMob the phase mobilities in the element
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the element wrt pressure
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the element wrt component fraction
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the element wrt component density
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[out] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    applyGradient( arrayView1d< real64 const > const & facePres,
                   arrayView1d< real64 const > const & dFacePres,
                   arrayView1d< real64 const > const & faceGravCoef,
                   arraySlice1d< localIndex const > const & elemToFaces,
                   real64 const & elemPres,
                   real64 const & dElemPres,
                   real64 const & elemGravCoef,
                   arraySlice1d< real64 const > const & elemPhaseMassDens,
                   arraySlice1d< real64 const > const & dElemPhaseMassDens_dPres,
                   arraySlice2d< real64 const > const & dElemPhaseMassDens_dCompFrac,
                   arraySlice1d< real64 const > const & elemPhaseMob,
                   arraySlice1d< real64 const > const & dElemPhaseMob_dPres,
                   arraySlice2d< real64 const > const & dElemPhaseMob_dCompDens,
                   arraySlice2d< real64 const > const & dElemCompFrac_dCompDens,
                   arraySlice2d< real64 const > const & transMatrix,
                   real64 ( &oneSidedVolFlux )[ NF ],
                   real64 ( &dOneSidedVolFlux_dPres )[ NF ],
                   real64 ( &dOneSidedVolFlux_dFacePres )[ NF ][ NF ],
                   real64 ( &dOneSidedVolFlux_dCompDens )[ NF ][ NC ] );

  /**
   * @brief In a given element, compute the flux divergence, i.e, sum the fluxes at this element's faces
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] rankOffset offset of this rank
   * @param[in] faceDofNumber the dof numbers of the face pressures in the domain
   * @param[in] elemRegionList map from face to element region index
   * @param[in] elemSubRegionList map from face to element subRegion index
   * @param[in] elemList map from face to element index
   * @param[in] regionFilter set of target regions of the solver
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] phaseDens phaseDens the phase densities in the domain
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities wrt pressure
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities wrt component fraction
   * @param[in] phaseMob the phase mobilities in the domain
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the element wrt pressure
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the element wrt component density
   * @param[in] dCompFrac_dCompDens the derivatives of the component fractions wrt component density
   * @param[in] phaseCompFrac the phase component fractions in the domain
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions wrt pressure
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions wrt component fractions
   * @param[in] elemDofNumber the dof numbers of the element in the domain
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] dt the time step
   * @param[inout] localMatrix the Jacobian matrix
   * @param[inout] localRhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
                          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
                          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
                          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
                          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
                          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
                          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
                          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
                          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
                          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
                          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
                          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
                          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
                          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          arraySlice2d< real64 const > const & transMatrixGrav,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dPres)[ NF ],
                          real64 const (&dOneSidedVolFlux_dFacePres)[ NF ][ NF ],
                          real64 const (&dOneSidedVolFlux_dCompDens)[ NF ][ NC ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**
   * @brief In a given element, compute the viscous flux divergence, i.e, sum the viscous fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face centered pressure
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[in] upwPhaseViscCoef the upwinded viscous transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseViscCoef_dPres the derivative of the upwinded viscous transport coefficient wrt pressure
   * @param[in] dUpwPhaseViscCoef_dCompDens the derivative of the upwinded viscous transport coefficient wrt component density
   * @param[in] upwViscDofNumber degree of freedom number of the upwind element
   * @param[in] faceDofNumber degree of freedom number of the face
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dDivMassFluxes_dFaceVars the derivatives of the flux divergence wrt the face centered vars
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   * @param[inout] dofColIndicesFaceVars degrees of freedom of the faces involved in the flux divergence
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
                         real64 ( &divMassFluxes )[ NC ],
                         real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ],
                         real64 ( &dDivMassFluxes_dFaceVars )[ NC ][ NF ],
                         globalIndex ( &dofColIndicesElemVars )[ (NC+1)*(NF+1) ],
                         globalIndex ( &dofColIndicesFaceVars )[ NF ] );

  /**
   * @brief In a given element, compute the buoyancy flux divergence, i.e, sum the buoyancy fluxes at this element's faces
   * @param[in] ifaceLoc the local index of the face
   * @param[in] phaseGravTerm the gravCoef multiplied by the difference in phase densities
   * @param[in] dPhaseGravTerm_dPres the derivatives of the gravCoef multiplied by the difference in phase densities wrt pressure
   * @param[in] dPhaseGravTerm_dCompDens the derivatives of the gravCoef multiplied by the difference in phase densities wrt component
   * density
   * @param[in] upwPhaseGravCoef the upwinded buoyancy transport coefficient at this face (ifaceLoc)
   * @param[in] dUpwPhaseGravCoef_dPres the derivative of the upwinded buoyancy transport coefficient wrt pressure
   * @param[in] dUpwPhaseGravCoef_dCompDens the derivative of the upwinded buoyancy transport coefficient wrt component density
   * @param[in] neighborDofNumber the degree of freedom number of the neighbor element
   * @param[in] dt the time step
   * @param[inout] divMassFluxes the divergence of the fluxes in the element
   * @param[inout] dDivMassFluxes_dElemVars the derivatives of the flux divergence wrt the element centered vars (pres and comp dens)
   * @param[inout] dofColIndicesElemVars degrees of freedom of the cells involved in the flux divergence
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
    assembleBuoyancyFlux( localIndex const ifaceLoc,
                          real64 const (&phaseGravTerm)[ NP ][ NP-1 ],
                          real64 const (&dPhaseGravTerm_dPres)[ NP ][ NP-1 ][ 2 ],
                          real64 const (&dPhaseGravTerm_dCompDens)[ NP ][ NP-1 ][ 2 ][ NC ],
                          real64 const (&upwPhaseGravCoef)[ NP ][ NP-1 ][ NC ],
                          real64 const (&dUpwPhaseGravCoef_dPres)[ NP ][ NP-1 ][ NC ][ 2 ],
                          real64 const (&dUpwPhaseGravCoef_dCompDens)[ NP ][ NP-1 ][ NC ][ 2 ][ NC ],
                          real64 const & dt,
                          real64 ( &divMassFluxes )[ NC ],
                          real64 ( &dDivMassFluxes_dElemVars )[ NC ][ (NC+1)*(NF+1) ] );

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures in the subRegion
   * @param[in] faceGhostRank the ghost ranks of the face pressures in the subRegion
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] rankOffset the offset of this rank
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dPres the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dFacePres the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dOneSidedVolFlux_dCompDens the derivatives of the vol fluxes wrt to this element's component density
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
                           arrayView1d< real64 > const & localRhs );

};


/******************************** AssemblerKernel ********************************/

struct AssemblerKernel
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

  /**
   * @brief In a given element, assemble the mass conservation equations and the contribution of this element to the face
   * constraints
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] ei index of this element
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank ghost rank of each face
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] elemGhostRank the ghost rank of the element in which we assemble the fluxes
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF, localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
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
           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
           ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           arraySlice2d< real64 const > const & transMatrixGrav,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs );
};


/******************************** FluxKernel ********************************/

struct FluxKernel
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

  /**
   * @brief In a given subRegion, assemble the mass conservation equations and the contribution of the elements of this subRegion  to the
   * face constraints.
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion the subRegion in which we are going to assemble the fluxes
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] faceGhostRank  the ghost ranks of the face pressures
   * @param[in] facePres the pressures at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravDepth the depth at the mesh faces
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravDepth the depth at this element's center
   * @param[in] phaseDens the phase densities in the domain (non-local)
   * @param[in] dPhaseDens_dPres the derivatives of the phase densities in the domain wrt pressure (non-local)
   * @param[in] dPhaseDens_dCompFrac the derivatives of the phase densities in the domain wrt component fraction (non-local)
   * @param[in] phaseMob the phase mobilities in the domain (non-local)
   * @param[in] dPhaseMob_dPres the derivatives of the phase mobilities in the domain wrt pressure (non-local)
   * @param[in] dPhaseMob_dCompDens the derivatives of the phase mobilities in the domain wrt component fraction (non-local)
   * @param[in] dCompFrac_dCompDens the derivatives of the component fraction in the domain wrt component density (non-local)
   * @param[in] phaseCompFrac the phase component fractions in the domain (non-local)
   * @param[in] dPhaseCompFrac_dPres the derivatives of the phase component fractions in the domain wrt pressure (non-local)
   * @param[in] dPhaseCompFrac_dCompFrac the derivatives of the phase component fractions in the domain wrt component fraction (non-local)
   * @param[in] elemDofNumber the dof number of the cell centered pressures (non-local)
   * @param[in] rankOffset the offset of this rank
   * @param[in] lengthTolerance tolerance used in the transmissibility matrix computation
   * @param[in] dt time step size
   * @param[inout] matrix the system matrix
   * @param[inout] rhs the system right-hand side vector
   */
  template< localIndex NF, localIndex NC, localIndex NP, typename IP_TYPE >
  static void
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
          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dCompFrac,
          ElementViewConst< arrayView3d< real64 const > > const & phaseMassDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMassDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseMassDens_dCompFrac,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dCompFrac,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          globalIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );

};


/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const > const & phaseVisc,
           arraySlice1d< real64 const > const & dPhaseVisc_dPres,
           arraySlice2d< real64 const > const & dPhaseVisc_dComp,
           arraySlice1d< real64 const > const & phaseRelPerm,
           arraySlice2d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dComp,
           arraySlice1d< real64 > const & phaseMob,
           arraySlice1d< real64 > const & dPhaseMob_dPres,
           arraySlice2d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView3d< real64 const > const & dCompFrac_dCompDens,
          arrayView3d< real64 const > const & phaseVisc,
          arrayView3d< real64 const > const & dPhaseVisc_dPres,
          arrayView4d< real64 const > const & dPhaseVisc_dComp,
          arrayView3d< real64 const > const & phaseRelPerm,
          arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac,
          arrayView2d< real64 const > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const > const & dPhaseVolFrac_dComp,
          arrayView2d< real64 > const & phaseMob,
          arrayView2d< real64 > const & dPhaseMob_dPres,
          arrayView3d< real64 > const & dPhaseMob_dComp );
};


/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY, typename REDUCE_POLICY >
  static void
  launch( arrayView1d< real64 const > const & localResidual,
          globalIndex const rankOffset,
          localIndex const numPhases,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          SortedArrayView< localIndex const > const & regionFilter,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMobOld,
          real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( facePresDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      // if not ghost face and if adjacent to target region, increment the residual norm
      if( faceGhostRank[iface] < 0 && facePresDofNumber[iface] >= 0 )
      {
        real64 normalizer = 0;
        localIndex elemCounter = 0;
        for( localIndex k=0; k<elemRegionList.size( 1 ); ++k )
        {
          localIndex const er  = elemRegionList[iface][k];
          localIndex const esr = elemSubRegionList[iface][k];
          localIndex const ei  = elemList[iface][k];

          bool const onBoundary = (er == -1 || esr == -1 || ei == -1);
          bool const isInTarget = regionFilter.contains( er );

          if( !onBoundary && isInTarget )
          {
            // compute a normalizer to obtain a dimensionless norm
            real64 sumMobOld = 0.0;
            for( localIndex ip = 0; ip < numPhases; ++ip )
            {
              sumMobOld += phaseMobOld[er][esr][ei][ip];
            }
            real64 const totalMobOld = ( sumMobOld < 1e-3 ) ? 1e-3 : sumMobOld;
            normalizer += elemVolume[er][esr][ei] / totalMobOld;
            elemCounter++;
          }
        }
        normalizer /= elemCounter;

        localIndex const lid = LvArray::integerConversion< localIndex >( facePresDofNumber[iface] - rankOffset );
        // note: unit of localResidual[lid] * totalMobOld: m^3, so this is dimensionless
        real64 const val = localResidual[lid] / normalizer;
        sumScaled += val * val;
      }
    } );

    localResidualNorm = localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & facePres,
          arrayView1d< real64 const > const & dFacePres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      if( ghostRank[iface] < 0 && dofNumber[iface] >= 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofNumber[iface] - rankOffset );
        {
          real64 const newFacePres = facePres[iface] + dFacePres[iface] + scalingFactor * localSolution[localRow];
          check.min( newFacePres >= 0.0 );
        }
      }
    } );
    return check.get();
  }

};

/******************************** PrecomputeKernel ********************************/

struct PrecomputeKernel
{

  template< typename IP_TYPE, localIndex NF >
  static void
  launch( localIndex const subRegionSize,
          localIndex const faceManagerSize,
          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
          ArrayOfArraysView< localIndex const > const & faceToNodes,
          arrayView2d< real64 const > const & elemCenter,
          arrayView1d< real64 const > const & elemVolume,
          arrayView3d< real64 const > const & elemPerm,
          arrayView1d< real64 const > const & elemGravCoef,
          arrayView2d< localIndex const > const & elemToFaces,
          arrayView1d< real64 const > const & transMultiplier,
          real64 const & lengthTolerance,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefNumerator,
          arrayView1d< RAJA::ReduceSum< serialReduce, real64 > > const & mimFaceGravCoefDenominator,
          arrayView1d< real64 > const & mimFaceGravCoef )
  {
    forAll< serialPolicy >( subRegionSize, [=] ( localIndex const ei )
    {
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

      real64 const perm[ 3 ] = { elemPerm[ei][0][0], elemPerm[ei][0][1], elemPerm[ei][0][2] };

      IP_TYPE::template compute< NF >( nodePosition,
                                       transMultiplier,
                                       faceToNodes,
                                       elemToFaces[ei],
                                       elemCenter[ei],
                                       elemVolume[ei],
                                       perm,
                                       lengthTolerance,
                                       transMatrix );

      for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
      {
        mimFaceGravCoefNumerator[elemToFaces[ei][ifaceLoc]] += elemGravCoef[ei] * transMatrix[ifaceLoc][ifaceLoc];
        mimFaceGravCoefDenominator[elemToFaces[ei][ifaceLoc]] += transMatrix[ifaceLoc][ifaceLoc];
      }
    } );

    forAll< serialPolicy >( faceManagerSize, [=] ( localIndex const iface )
    {
      if( !isZero( mimFaceGravCoefDenominator[iface].get() ) )
      {
        mimFaceGravCoef[iface] = mimFaceGravCoefNumerator[iface].get() / mimFaceGravCoefDenominator[iface].get();
      }
    } );
  }
};

/******************************** Kernel switches ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorFaceSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorFaceSwitch: type should be integral" );

  switch( value )
  {
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return;}
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return;}
    case 6:
    { lambda( std::integral_constant< T, 6 >() ); return;}
    default: GEOSX_ERROR( "Unknown numFacesInElem value: " << value );
  }
}

} // namespace internal

template< typename KERNELWRAPPER, typename IP_TYPE, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, localIndex numComps, localIndex numPhases, ARGS && ... args )
{
  internal::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    if( numPhases == 2 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template launch< NF(), 2, 2, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template launch< NF(), 3, 2, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template launch< NF(), 4, 2, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template launch< NF(), 5, 2, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else if( numPhases == 3 )
    {
      if( numComps == 2 )
      {
        KERNELWRAPPER::template launch< NF(), 2, 3, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 3 )
      {
        KERNELWRAPPER::template launch< NF(), 3, 3, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 4 )
      {
        KERNELWRAPPER::template launch< NF(), 4, 3, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else if( numComps == 5 )
      {
        KERNELWRAPPER::template launch< NF(), 5, 3, IP_TYPE >( std::forward< ARGS >( args )... );
      }
      else
      {
        GEOSX_ERROR( "Unsupported number of components: " << numComps );
      }
    }
    else
    {
      GEOSX_ERROR( "Unsupported number of phases: " << numPhases );
    }
  } );
}

} // namespace CompositionalMultiphaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVMKERNELS_HPP
