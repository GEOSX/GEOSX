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
 * @file SinglePhaseHybridFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/mimeticInnerProducts/MimeticInnerProductBase.hpp"
#include "finiteVolume/mimeticInnerProducts/TPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiTPFAInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/QuasiRTInnerProduct.hpp"
#include "finiteVolume/mimeticInnerProducts/SimpleInnerProduct.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/MeshLevel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/HybridFVMHelperKernels.hpp"

namespace geosx
{
namespace SinglePhaseHybridFVMKernels
{

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
   * @brief In a given element, compute the one-sided volumetric fluxes at this element's faces
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh facesb
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] elemPres the pressure at this element's center
   * @param[in] dElemPres the accumulated pressure updates at this element's center
   * @param[in] elemGravCoef the depth at this element's center
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[out] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[out] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[out] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  applyGradient( arrayView1d< real64 const > const & facePres,
                 arrayView1d< real64 const > const & dFacePres,
                 arrayView1d< real64 const > const & faceGravCoef,
                 arraySlice1d< localIndex const > const & elemToFaces,
                 real64 const & elemPres,
                 real64 const & dElemPres,
                 real64 const & elemGravCoef,
                 real64 const & elemDens,
                 real64 const & dElemDens_dp,
                 arraySlice2d< real64 const > const & transMatrix,
                 real64 ( & oneSidedVolFlux )[ NF ],
                 real64 ( & dOneSidedVolFlux_dp )[ NF ],
                 real64 ( & dOneSidedVolFlux_dfp )[ NF ][ NF ] )
  {
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      // now in the following nested loop,
      // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        // 1) compute the potential diff between the cell center and the face center
        real64 const ccPres = elemPres + dElemPres;
        real64 const fPres  = facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

        real64 const ccGravCoef = elemGravCoef;
        real64 const fGravCoef  = faceGravCoef[elemToFaces[jfaceLoc]];

        real64 const ccDens     = elemDens;
        real64 const dCcDens_dp = dElemDens_dp;
        // no density evaluated at the face center

        // pressure difference
        real64 const presDif      = ccPres - fPres;
        real64 const dPresDif_dp  =  1;
        real64 const dPresDif_dfp = -1;

        // gravity term
        real64 const gravCoefDif  = ccGravCoef - fGravCoef;
        real64 const gravTerm     = ccDens     * gravCoefDif;
        real64 const dGravTerm_dp = dCcDens_dp * gravCoefDif;

        // potential difference
        real64 const potDif      = presDif     - gravTerm;
        real64 const dPotDif_dp  = dPresDif_dp - dGravTerm_dp;
        real64 const dPotDif_dfp = dPresDif_dfp;

        // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
        oneSidedVolFlux[ifaceLoc]                = oneSidedVolFlux[ifaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * potDif;
        dOneSidedVolFlux_dp[ifaceLoc]            = dOneSidedVolFlux_dp[ifaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp;
        dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc]
                                                   + transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp;
      }
    }
  }

  /**
   * @brief In a given element, assemble the mass conservation equation
   * @param[in] localIds region, subRegion, and element indices of the local element
   * @param[in] rankOffset offset of this rank
   * @param[in] faceDofNumber the dof numbers of the face pressures in the domain
   * @param[in] elemRegionList map from face to element region index
   * @param[in] elemSubRegionList map from face to element subRegion index
   * @param[in] elemList map from face to element index
   * @param[in] regionFilter set of target regions of the solver
   * @param[in] elemToFaces the map from one-sided face to face
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof numbers of all the cell centered pressures (non-local)
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[in] dt the time step size
   * @param[inout] matrix the jacobian matrix
   * @param[inout] rhs the residual
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  assembleFluxDivergence( localIndex const (&localIds)[ 3 ],
                          globalIndex const rankOffset,
                          arrayView2d< localIndex const > const & elemRegionList,
                          arrayView2d< localIndex const > const & elemSubRegionList,
                          arrayView2d< localIndex const > const & elemList,
                          SortedArrayView< localIndex const > const & regionFilter,
                          arrayView1d< globalIndex const > const & faceDofNumber,
                          arraySlice1d< localIndex const > const & elemToFaces,
                          ElementViewConst< arrayView1d< real64 const > > const & mob,
                          ElementViewConst< arrayView1d< real64 const > > const & dMob_dp,
                          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
                          real64 const (&oneSidedVolFlux)[ NF ],
                          real64 const (&dOneSidedVolFlux_dp)[ NF ],
                          real64 const (&dOneSidedVolFlux_dfp)[ NF ][ NF ],
                          real64 const & dt,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs )
  {
    // fluxes
    real64 divMassFluxes = 0;
    real64 dDivMassFluxes_dElemVars[ NF+1 ]{};
    real64 dDivMassFluxes_dFaceVars[ NF ]{};

    // dof numbers
    globalIndex const eqnRowLocalIndex = elemDofNumber[localIds[0]][localIds[1]][localIds[2]] - rankOffset;
    globalIndex elemDofColIndices[ NF+1 ]{};
    globalIndex faceDofColIndices[ NF ]{};
    elemDofColIndices[0] = elemDofNumber[localIds[0]][localIds[1]][localIds[2]];

    // upwinded mobility
    real64 upwMob = 0.0;
    real64 dUpwMob_dp = 0.0;
    globalIndex upwDofNumber = 0;

    // for each element, loop over the one-sided faces
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {

      // 1) Find if there is a neighbor, and if there is, grab the indices of the neighbor element
      localIndex neighborIds[ 3 ] = { localIds[0], localIds[1], localIds[2] };
      bool const hasNeighbor = hybridFVMKernels::CellConnectivity::findNeighbor( localIds,
                                                                                 ifaceLoc,
                                                                                 elemRegionList,
                                                                                 elemSubRegionList,
                                                                                 elemList,
                                                                                 regionFilter,
                                                                                 elemToFaces,
                                                                                 neighborIds );

      // 2) Upwind the mobility at this face
      if( oneSidedVolFlux[ifaceLoc] >= 0 || !hasNeighbor )
      {
        upwMob       = mob[localIds[0]][localIds[1]][localIds[2]];
        dUpwMob_dp   = dMob_dp[localIds[0]][localIds[1]][localIds[2]];
        upwDofNumber = elemDofNumber[localIds[0]][localIds[1]][localIds[2]];
      }
      else
      {
        upwMob       = mob[neighborIds[0]][neighborIds[1]][neighborIds[2]];
        dUpwMob_dp   = dMob_dp[neighborIds[0]][neighborIds[1]][neighborIds[2]];
        upwDofNumber = elemDofNumber[neighborIds[0]][neighborIds[1]][neighborIds[2]];
      }

      // 3) Add to the flux divergence
      real64 const dt_upwMob = dt * upwMob;

      // compute the mass flux at the one-sided face plus its derivatives and add the newly computed flux to the sum
      divMassFluxes                        = divMassFluxes + dt_upwMob * oneSidedVolFlux[ifaceLoc];
      dDivMassFluxes_dElemVars[0]          = dDivMassFluxes_dElemVars[0] + dt * upwMob * dOneSidedVolFlux_dp[ifaceLoc];
      dDivMassFluxes_dElemVars[ifaceLoc+1] = dt * dUpwMob_dp * oneSidedVolFlux[ifaceLoc];
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dDivMassFluxes_dFaceVars[jfaceLoc] = dDivMassFluxes_dFaceVars[jfaceLoc]
                                             + dt_upwMob * dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
      }

      // collect the relevant dof numbers
      elemDofColIndices[ifaceLoc+1] = upwDofNumber;
      faceDofColIndices[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
    }

    // we are ready to assemble the local flux and its derivatives
    // no need for atomic adds - each row is assembled by a single thread

    // residual
    localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + divMassFluxes;

    // jacobian -- derivative wrt elem centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &elemDofColIndices[0],
                                                              &dDivMassFluxes_dElemVars[0],
                                                              NF+1 );

    // jacobian -- derivatives wrt face centered vars
    localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowLocalIndex,
                                                              &faceDofColIndices[0],
                                                              &dDivMassFluxes_dFaceVars[0],
                                                              NF );
  }

  /**
   * @brief In a given element, assemble the constraints at this element's faces
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] elemToFaces the map from one-sided face to face to access face Dof numbers
   * @param[in] elemDofNumber the dof number of this element's cell centered pressure
   * @param[in] rankOffset the offset of this rank
   * @param[in] oneSidedVolFlux the volumetric fluxes at this element's faces
   * @param[in] dOneSidedVolFlux_dp the derivatives of the vol fluxes wrt to this element's cell centered pressure
   * @param[in] dOneSidedVolFlux_dfp the derivatives of the vol fluxes wrt to this element's face pressures
   * @param[inout] matrix the local Jacobian matrix
   * @param[inout] rhs the local residual
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  assembleFaceConstraints( arrayView1d< globalIndex const > const & faceDofNumber,
                           arrayView1d< integer const > const & faceGhostRank,
                           arraySlice1d< localIndex const > const & elemToFaces,
                           globalIndex const elemDofNumber,
                           globalIndex const rankOffset,
                           real64 const (&oneSidedVolFlux)[ NF ],
                           real64 const (&dOneSidedVolFlux_dp)[ NF ],
                           real64 const (&dOneSidedVolFlux_dfp)[ NF ][ NF ],
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs )
  {
    // fluxes
    real64 dFlux_dfp[ NF ]{};

    // dof numbers
    globalIndex dofColIndicesFacePres[ NF ]{};
    globalIndex const dofColIndexElemPres = elemDofNumber;    // fluxes

    // for each element, loop over the local (one-sided) faces
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
      {
        continue;
      }

      // flux at this face
      real64 const flux      = oneSidedVolFlux[ifaceLoc];
      real64 const dFlux_dp  = dOneSidedVolFlux_dp[ifaceLoc];

      // dof number of this face constraint
      localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >( faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        dFlux_dfp[jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
        dofColIndicesFacePres[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
      }

      // residual
      atomicAdd( parallelDeviceAtomic{}, &localRhs[eqnLocalRowIndex], flux );

      // jacobian -- derivative wrt local cell centered pressure term
      localMatrix.addToRow< parallelDeviceAtomic >( eqnLocalRowIndex, &dofColIndexElemPres, &dFlux_dp, 1 );

      // jacobian -- derivatives wrt face pressure terms
      localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                                        &dofColIndicesFacePres[0],
                                                                        &dFlux_dfp[0],
                                                                        NF );
    }
  }

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
   * @brief In a given element, assemble the mass conservation equation and the contribution of this element to the face
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
   * @param[in] elemDens the density at this elenent's center
   * @param[in] dElemDens_dp the derivative of the density at this element's center
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof numbers of the cells in the domain (non-local)
   * @param[in] elemGhostRank the ghost rank of the cell
   * @param[in] rankOffset the offset of this rank
   * @param[in] dt time step size
   * @param[in] transMatrix the transmissibility matrix in this element
   * @param[inout] localMatrix the local Jacobian matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const er,
           localIndex const esr,
           localIndex const ei,
           SortedArrayView< localIndex const > const & regionFilter,
           arrayView2d< localIndex const > const & elemRegionList,
           arrayView2d< localIndex const > const & elemSubRegionList,
           arrayView2d< localIndex const > const & elemList,
           arrayView1d< globalIndex const > const & faceDofNumber,
           arrayView1d< integer const > const & faceGhostRank,
           arrayView1d< real64 const > const & facePres,
           arrayView1d< real64 const > const & dFacePres,
           arrayView1d< real64 const > const & faceGravCoef,
           arraySlice1d< localIndex const > const & elemToFaces,
           real64 const & elemPres,
           real64 const & dElemPres,
           real64 const & elemGravCoef,
           real64 const & elemDens,
           real64 const & dElemDens_dp,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dp,
           ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
           integer const elemGhostRank,
           globalIndex const rankOffset,
           real64 const & dt,
           arraySlice2d< real64 const > const & transMatrix,
           CRSMatrixView< real64, globalIndex const > const & localMatrix,
           arrayView1d< real64 > const & localRhs )
  {
    // one sided flux
    real64 oneSidedVolFlux[ NF ]{};
    real64 dOneSidedVolFlux_dp[ NF ]{};
    real64 dOneSidedVolFlux_dfp[ NF ][ NF ]{};

    localIndex const localIds[3] = { er, esr, ei };

    /*
     * compute auxiliary quantities at the one sided faces of this element:
     * 1) One-sided volumetric fluxes
     * 2) Upwinded mobilities
     */

    // for each one-sided face of the elem,
    // compute the volumetric flux using transMatrix
    AssemblerKernelHelper::applyGradient< NF >( facePres,
                                                dFacePres,
                                                faceGravCoef,
                                                elemToFaces,
                                                elemPres,
                                                dElemPres,
                                                elemGravCoef,
                                                elemDens,
                                                dElemDens_dp,
                                                transMatrix,
                                                oneSidedVolFlux,
                                                dOneSidedVolFlux_dp,
                                                dOneSidedVolFlux_dfp );

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
      AssemblerKernelHelper::assembleFluxDivergence< NF >( localIds,
                                                           rankOffset,
                                                           elemRegionList,
                                                           elemSubRegionList,
                                                           elemList,
                                                           regionFilter,
                                                           faceDofNumber,
                                                           elemToFaces,
                                                           mob,
                                                           dMob_dp,
                                                           elemDofNumber,
                                                           oneSidedVolFlux,
                                                           dOneSidedVolFlux_dp,
                                                           dOneSidedVolFlux_dfp,
                                                           dt,
                                                           localMatrix,
                                                           localRhs );
    }

    // use the computed one sided vol fluxes to assemble the constraints
    // enforcing flux continuity at this element's faces
    AssemblerKernelHelper::assembleFaceConstraints< NF >( faceDofNumber,
                                                          faceGhostRank,
                                                          elemToFaces,
                                                          elemDofNumber[er][esr][ei],
                                                          rankOffset,
                                                          oneSidedVolFlux,
                                                          dOneSidedVolFlux_dp,
                                                          dOneSidedVolFlux_dfp,
                                                          localMatrix,
                                                          localRhs );
  }


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
   * @brief Assemble the mass conservation equations and face constraints in the cell subregion
   * @param[in] er index of this element's region
   * @param[in] esr index of this element's subregion
   * @param[in] subRegion pointer to the cell element subregion
   * @param[in] fluid the (single-phase) fluid model associated with this subRegion
   * @param[in] regionFilter set containing the indices of the target regions
   * @param[in] nodePosition position of the nodes
   * @param[in] elemRegionList face-to-elemRegions map
   * @param[in] elemSubRegionList face-to-elemSubRegions map
   * @param[in] elemList face-to-elemIds map
   * @param[in] faceToNodes map from face to nodes
   * @param[in] faceDofNumber the dof numbers of the face pressures
   * @param[in] facePres the pressure at the mesh faces at the beginning of the time step
   * @param[in] dFacePres the accumulated pressure updates at the mesh face
   * @param[in] faceGravCoef the depth at the mesh faces
   * @param[in] transMultiplier the transmissibility multiplier at the mesh faces
   * @param[in] mob the mobilities in the domain (non-local)
   * @param[in] dMob_dp the derivatives of the mobilities in the domain wrt cell-centered pressure (non-local)
   * @param[in] elemDofNumber the dof numbers of the cells in the domain (non-local)
   * @param[in] rankOffset the offset of this rank
   * @param[in] dt time step size
   * @param[inout] localMatrix the local Jacobian matrix
   * @param[inout] localRhs the local right-hand side vector
   */
  template< typename IP_TYPE, localIndex NF >
  static void
  launch( localIndex er,
          localIndex esr,
          CellElementSubRegion const & subRegion,
          constitutive::SingleFluidBase const & fluid,
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
          arrayView1d< real64 const > const & transMultiplier,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dp,
          ElementViewConst< arrayView1d< globalIndex const > > const & elemDofNumber,
          localIndex const rankOffset,
          real64 const lengthTolerance,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    // get the cell-centered DOF numbers and ghost rank for the assembly
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    // get the map from elem to faces
    arrayView2d< localIndex const > const elemToFaces = subRegion.faceList().toViewConst();

    // get the cell-centered pressures
    arrayView1d< real64 const > const elemPres  =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const dElemPres =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() );

    // get the element data needed for transmissibility computation
    arrayView2d< real64 const > const elemCenter =
      subRegion.getReference< array2d< real64 > >( CellElementSubRegion::viewKeyStruct::elementCenterString() );
    arrayView1d< real64 const > const elemVolume =
      subRegion.getReference< array1d< real64 > >( CellElementSubRegion::viewKeyStruct::elementVolumeString() );
    arrayView2d< real64 const > const elemPerm =
      subRegion.getReference< array2d< real64 > >( SinglePhaseBase::viewKeyStruct::permeabilityString() );

    // get the cell-centered depth
    arrayView1d< real64 const > const elemGravCoef =
      subRegion.getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::gravityCoefString() );

    // get the fluid data
    arrayView2d< real64 const > const elemDens = fluid.density();
    arrayView2d< real64 const > const dElemDens_dp = fluid.dDensity_dPressure();

    // assemble the residual and Jacobian element by element
    // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
    using KERNEL_POLICY = parallelDevicePolicy< 32 >;
    forAll< KERNEL_POLICY >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex const ei )
    {

      // transmissibility matrix
      stackArray2d< real64, NF *NF > transMatrix( NF, NF );

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

      // perform flux assembly in this element
      SinglePhaseHybridFVMKernels::AssemblerKernel::compute< NF >( er, esr, ei,
                                                                   regionFilter,
                                                                   elemRegionList,
                                                                   elemSubRegionList,
                                                                   elemList,
                                                                   faceDofNumber,
                                                                   faceGhostRank,
                                                                   facePres,
                                                                   dFacePres,
                                                                   faceGravCoef,
                                                                   elemToFaces[ei],
                                                                   elemPres[ei],
                                                                   dElemPres[ei],
                                                                   elemGravCoef[ei],
                                                                   elemDens[ei][0],
                                                                   dElemDens_dp[ei][0],
                                                                   mob,
                                                                   dMob_dp,
                                                                   elemDofNumber,
                                                                   elemGhostRank[ei],
                                                                   rankOffset,
                                                                   dt,
                                                                   transMatrix,
                                                                   localMatrix,
                                                                   localRhs );

    } );
  }


};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void
  launch( LOCAL_VECTOR const localResidual,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & facePresDofNumber,
          arrayView1d< integer const > const & faceGhostRank,
          arrayView2d< localIndex const > const & elemRegionList,
          arrayView2d< localIndex const > const & elemSubRegionList,
          arrayView2d< localIndex const > const & elemList,
          ElementViewConst< arrayView1d< real64 const > > const & elemVolume,
          real64 const & defaultViscosity,
          real64 * localResidualNorm )
  {

    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( facePresDofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
    {
      // if not ghost face and if adjacent to target region
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

          // if not on boundary, save the mobility and the upwDofNumber
          if( !onBoundary )
          {
            normalizer += elemVolume[er][esr][ei];
            elemCounter++;
          }
        }
        normalizer /= elemCounter;
        normalizer /= defaultViscosity;

        localIndex const lid = LvArray::integerConversion< localIndex >( facePresDofNumber[iface] - rankOffset );
        real64 const val = localResidual[lid] / normalizer; // to get something dimensionless
        sumScaled += val * val;
      }
    } );

    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

/******************************** Kernel switches ********************************/

namespace helpers
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

} // namespace helpers

template< typename IP_TYPE, typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector( localIndex numFacesInElem, ARGS && ... args )
{
  helpers::KernelLaunchSelectorFaceSwitch( numFacesInElem, [&] ( auto NF )
  {
    KERNELWRAPPER::template launch< IP_TYPE, NF() >( std::forward< ARGS >( args )... );
  } );
}


} // namespace SinglePhaseHybridFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVMKERNELS_HPP
