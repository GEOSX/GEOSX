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
 * @file SinglePhaseHybridFVMKernels.cpp
 */

#include "SinglePhaseHybridFVMKernels.hpp"

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{
namespace SinglePhaseHybridFVMKernels
{
/******************************** AssemblerKernelHelper ********************************/

template< localIndex NF >
GEOSX_HOST_DEVICE void
AssemblerKernelHelper::ComputeOneSidedVolFluxes(
  arrayView1d< real64 const > const & facePres,
  arrayView1d< real64 const > const & dFacePres,
  arrayView1d< real64 const > const & faceGravCoef,
  arraySlice1d< localIndex const > const & elemToFaces,
  real64 const & elemPres,
  real64 const & dElemPres,
  real64 const & elemGravCoef,
  real64 const & elemDens,
  real64 const & dElemDens_dp,
  arraySlice2d< real64 const > const & transMatrix,
  arraySlice1d< real64 > const & oneSidedVolFlux,
  arraySlice1d< real64 > const & dOneSidedVolFlux_dp,
  arraySlice2d< real64 > const & dOneSidedVolFlux_dfp )
{
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // now in the following nested loop,
    // we compute the contribution of face jfaceLoc to the one sided total volumetric flux at face iface
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      // 1) compute the potential diff between the cell center and the face center
      real64 const ccPres = elemPres + dElemPres;
      real64 const fPres =
        facePres[elemToFaces[jfaceLoc]] + dFacePres[elemToFaces[jfaceLoc]];

      real64 const ccGravCoef = elemGravCoef;
      real64 const fGravCoef = faceGravCoef[elemToFaces[jfaceLoc]];

      real64 const ccDens = elemDens;
      real64 const dCcDens_dp = dElemDens_dp;
      // no density evaluated at the face center

      // pressure difference
      real64 const presDif = ccPres - fPres;
      real64 const dPresDif_dp = 1;
      real64 const dPresDif_dfp = -1;

      // gravity term
      real64 const gravCoefDif = ccGravCoef - fGravCoef;
      real64 const gravTerm = ccDens * gravCoefDif;
      real64 const dGravTerm_dp = dCcDens_dp * gravCoefDif;

      // potential difference
      real64 const potDif = presDif - gravTerm;
      real64 const dPotDif_dp = dPresDif_dp - dGravTerm_dp;
      real64 const dPotDif_dfp = dPresDif_dfp;

      // this is going to store T \sum_p \lambda_p (\nabla p - \rho_p g \nabla d)
      oneSidedVolFlux[ifaceLoc] =
        oneSidedVolFlux[ifaceLoc] + transMatrix[ifaceLoc][jfaceLoc] * potDif;
      dOneSidedVolFlux_dp[ifaceLoc] = dOneSidedVolFlux_dp[ifaceLoc] +
        transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dp;
      dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] =
        dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc] +
        transMatrix[ifaceLoc][jfaceLoc] * dPotDif_dfp;
    }
  }
}

template< localIndex NF >
GEOSX_HOST_DEVICE void
AssemblerKernelHelper::UpdateUpwindedCoefficients(
  localIndex const er,
  localIndex const esr,
  localIndex const ei,
  arrayView2d< localIndex const > const & elemRegionList,
  arrayView2d< localIndex const > const & elemSubRegionList,
  arrayView2d< localIndex const > const & elemList,
  SortedArrayView< localIndex const > const & regionFilter,
  arraySlice1d< localIndex const > const & elemToFaces,
  ElementView< arrayView1d< real64 const > > const & mob,
  ElementView< arrayView1d< real64 const > > const & dMob_dp,
  ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
  arraySlice1d< real64 const > const & oneSidedVolFlux,
  arraySlice1d< real64 > const & upwMobility,
  arraySlice1d< real64 > const & dUpwMobility_dp,
  arraySlice1d< globalIndex > const & upwDofNumber )
{
  // for this element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    // we initialize these upw quantities with the values of the local elem
    upwMobility[ifaceLoc] = mob[er][esr][ei];
    dUpwMobility_dp[ifaceLoc] = dMob_dp[er][esr][ei];
    upwDofNumber[ifaceLoc] = elemDofNumber[er][esr][ei];

    // if the local elem if upstream, we are done, we can proceed to the next one-sided face
    // otherwise, we have to access the properties of the neighbor element
    // this is done on the fly below
    if( oneSidedVolFlux[ifaceLoc] < 0 )
    {
      // the face has at most two adjacent elements
      // one of these two elements is the current element indexed by er, esr, ei
      // but here we are interested in the indices of the other element
      // this other element is "the neighbor" for this one-sided face
      for( localIndex k = 0; k < elemRegionList.size( 1 ); ++k )
      {
        localIndex const erNeighbor = elemRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const esrNeighbor =
          elemSubRegionList[elemToFaces[ifaceLoc]][k];
        localIndex const eiNeighbor = elemList[elemToFaces[ifaceLoc]][k];

        // this element is not the current element
        // we have found the neighbor or we are at the boundary
        if( erNeighbor != er || esrNeighbor != esr || eiNeighbor != ei )
        {
          bool const onBoundary =
            ( erNeighbor == -1 || esrNeighbor == -1 || eiNeighbor == -1 );
          bool const neighborInTarget = regionFilter.contains( erNeighbor );

          // if not on boundary, save the mobility and the upwDofNumber
          if( !onBoundary && neighborInTarget )
          {
            upwMobility[ifaceLoc] = mob[erNeighbor][esrNeighbor][eiNeighbor];
            dUpwMobility_dp[ifaceLoc] =
              dMob_dp[erNeighbor][esrNeighbor][eiNeighbor];
            upwDofNumber[ifaceLoc] =
              elemDofNumber[erNeighbor][esrNeighbor][eiNeighbor];
          }
          // if the face is on the boundary, use the properties of the local elem
        }
      }
    }
  }
}

template< localIndex NF >
GEOSX_HOST_DEVICE void
AssemblerKernelHelper::AssembleOneSidedMassFluxes(
  arrayView1d< globalIndex const > const & faceDofNumber,
  arraySlice1d< localIndex const > const & elemToFaces,
  globalIndex const elemDofNumber,
  globalIndex const rankOffset,
  arraySlice1d< real64 const > const & oneSidedVolFlux,
  arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
  arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
  arraySlice1d< real64 const > const & upwMobility,
  arraySlice1d< real64 const > const & dUpwMobility_dp,
  arraySlice1d< globalIndex const > const & upwDofNumber,
  real64 const & dt,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  // fluxes
  real64 sumLocalMassFluxes = 0;
  stackArray1d< real64, 1 + NF > dSumLocalMassFluxes_dElemVars( 1 + NF );
  stackArray1d< real64, NF > dSumLocalMassFluxes_dFaceVars( NF );
  for( localIndex i = 0; i < NF + 1; ++i )
  {
    dSumLocalMassFluxes_dElemVars( i ) = 0.;
  }
  for( localIndex i = 0; i < NF; ++i )
  {
    dSumLocalMassFluxes_dFaceVars( i ) = 0.;
  }

  // dof numbers
  globalIndex const eqnRowLocalIndex = elemDofNumber - rankOffset;
  stackArray1d< globalIndex, 1 + NF > elemDofColIndices( 1 + NF );
  stackArray1d< globalIndex, NF > faceDofColIndices( NF );
  elemDofColIndices[0] = elemDofNumber;

  // for each element, loop over the one-sided faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    real64 const dt_upwMob = dt * upwMobility[ifaceLoc];
    // compute the mass flux at the one-sided face plus its derivatives
    // add the newly computed flux to the sum
    sumLocalMassFluxes =
      sumLocalMassFluxes + dt_upwMob * oneSidedVolFlux[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[0] = dSumLocalMassFluxes_dElemVars[0] +
      dt * upwMobility[ifaceLoc] * dOneSidedVolFlux_dp[ifaceLoc];
    dSumLocalMassFluxes_dElemVars[ifaceLoc + 1] =
      dt * dUpwMobility_dp[ifaceLoc] * oneSidedVolFlux[ifaceLoc];
    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dSumLocalMassFluxes_dFaceVars[jfaceLoc] =
        dSumLocalMassFluxes_dFaceVars[jfaceLoc] +
        dt_upwMob * dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
    }

    // collect the relevant dof numbers
    elemDofColIndices[ifaceLoc + 1] =
      upwDofNumber[ifaceLoc];  // if upwDofNumber == elemDofNumber, the derivative is zero
    faceDofColIndices[ifaceLoc] = faceDofNumber[elemToFaces[ifaceLoc]];
  }

  // we are ready to assemble the local flux and its derivatives
  // no need for atomic adds - each row is assembled by a single thread

  // residual
  localRhs[eqnRowLocalIndex] = localRhs[eqnRowLocalIndex] + sumLocalMassFluxes;

  // jacobian -- derivative wrt elem centered vars
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >(
    eqnRowLocalIndex,
    elemDofColIndices.data(),
    dSumLocalMassFluxes_dElemVars.data(),
    elemDofColIndices.size() );

  // jacobian -- derivatives wrt face centered vars
  localMatrix.addToRowBinarySearchUnsorted< serialAtomic >(
    eqnRowLocalIndex,
    faceDofColIndices.data(),
    dSumLocalMassFluxes_dFaceVars.data(),
    faceDofColIndices.size() );
}

template< localIndex NF >
GEOSX_HOST_DEVICE void
AssemblerKernelHelper::AssembleConstraints(
  arrayView1d< globalIndex const > const & faceDofNumber,
  arrayView1d< integer const > const & faceGhostRank,
  arraySlice1d< localIndex const > const & elemToFaces,
  globalIndex const elemDofNumber,
  globalIndex const rankOffset,
  arraySlice1d< real64 const > const & oneSidedVolFlux,
  arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,
  arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  // fluxes
  stackArray1d< real64, NF > dFlux_dfp( NF );

  // dof numbers
  stackArray1d< globalIndex, NF > dofColIndicesFacePres( NF );
  globalIndex const dofColIndexElemPres = elemDofNumber;

  // for each element, loop over the local (one-sided) faces
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    if( faceGhostRank[elemToFaces[ifaceLoc]] >= 0 )
    {
      continue;
    }

    // flux at this face
    real64 const flux = oneSidedVolFlux[ifaceLoc];
    real64 const dFlux_dp = dOneSidedVolFlux_dp[ifaceLoc];

    // dof number of this face constraint
    localIndex const eqnLocalRowIndex = LvArray::integerConversion< localIndex >(
      faceDofNumber[elemToFaces[ifaceLoc]] - rankOffset );

    for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
    {
      dFlux_dfp[jfaceLoc] = dOneSidedVolFlux_dfp[ifaceLoc][jfaceLoc];
      dofColIndicesFacePres[jfaceLoc] = faceDofNumber[elemToFaces[jfaceLoc]];
    }

    // residual
    atomicAdd( parallelDeviceAtomic {}, &localRhs[eqnLocalRowIndex], flux );

    // jacobian -- derivative wrt local cell centered pressure term
    localMatrix.addToRow< parallelDeviceAtomic >( eqnLocalRowIndex,
                                                  &dofColIndexElemPres,
                                                  &dFlux_dp,
                                                  1 );

    // jacobian -- derivatives wrt face pressure terms
    localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
      eqnLocalRowIndex,
      dofColIndicesFacePres.data(),
      dFlux_dfp.data(),
      NF );
  }
}

/******************************** AssemblerKernel ********************************/

template< localIndex NF >
GEOSX_HOST_DEVICE void
AssemblerKernel::Compute(
  localIndex const er,
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
  ElementView< arrayView1d< real64 const > > const & mobility,
  ElementView< arrayView1d< real64 const > > const & dMobility_dp,
  ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
  integer const elemGhostRank,
  globalIndex const rankOffset,
  real64 const & dt,
  arraySlice2d< real64 const > const & transMatrix,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  // one sided flux
  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF * NF > dOneSidedVolFlux_dfp( NF, NF );
  for( localIndex i = 0; i < NF; ++i )
  {
    oneSidedVolFlux( i ) = 0.;
    dOneSidedVolFlux_dp( i ) = 0.;
    for( localIndex j = 0; j < NF; ++j )
    {
      dOneSidedVolFlux_dfp( i, j ) = 0.;  // assume row major
    }
  }

  // upwinded mobility
  stackArray1d< real64, NF > upwMobility( NF );
  stackArray1d< real64, NF > dUpwMobility_dp( NF );
  stackArray1d< globalIndex, NF > upwDofNumber( NF );

  /*
   * compute auxiliary quantities at the one sided faces of this element:
   * 1) One-sided volumetric fluxes
   * 2) Upwinded mobilities
   */

  // for each one-sided face of the elem,
  // compute the volumetric flux using transMatrix
  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
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
    AssemblerKernelHelper::UpdateUpwindedCoefficients< NF >( er,
                                                             esr,
                                                             ei,
                                                             elemRegionList,
                                                             elemSubRegionList,
                                                             elemList,
                                                             regionFilter,
                                                             elemToFaces,
                                                             mobility,
                                                             dMobility_dp,
                                                             elemDofNumber,
                                                             oneSidedVolFlux,
                                                             upwMobility,
                                                             dUpwMobility_dp,
                                                             upwDofNumber );

    /*
     * perform assembly in this element in two steps:
     * 1) mass conservation equations
     * 2) face constraints
     */

    // use the computed one sided vol fluxes and the upwinded mobilities
    // to assemble the upwinded mass fluxes in the mass conservation eqn of the elem
    AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >(
      faceDofNumber,
      elemToFaces,
      elemDofNumber[er][esr][ei],
      rankOffset,
      oneSidedVolFlux,
      dOneSidedVolFlux_dp,
      dOneSidedVolFlux_dfp,
      upwMobility,
      dUpwMobility_dp,
      upwDofNumber,
      dt,
      localMatrix,
      localRhs );
  }

  // use the computed one sided vol fluxes to assemble the constraints
  // enforcing flux continuity at this element's faces
  AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
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

/******************************** FluxKernel ********************************/

template< localIndex NF >
void
FluxKernel::Launch(
  localIndex er,
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
  ElementView< arrayView1d< real64 const > > const & mobility,
  ElementView< arrayView1d< real64 const > > const & dMobility_dp,
  ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,
  localIndex const rankOffset,
  real64 const lengthTolerance,
  real64 const dt,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  // get the cell-centered DOF numbers and ghost rank for the assembly
  arrayView1d< integer const > const & elemGhostRank =
    subRegion.getReference< array1d< integer > >(
      ObjectManagerBase::viewKeyStruct::ghostRankString );

  // get the map from elem to faces
  arrayView2d< localIndex const > const & elemToFaces = subRegion.faceList();

  // get the cell-centered pressures
  arrayView1d< real64 const > const & elemPres =
    subRegion.getReference< array1d< real64 > >(
      SinglePhaseBase::viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dElemPres =
    subRegion.getReference< array1d< real64 > >(
      SinglePhaseBase::viewKeyStruct::deltaPressureString );

  // get the element data needed for transmissibility computation
  arrayView2d< real64 const > const & elemCenter =
    subRegion.getReference< array2d< real64 > >(
      CellBlock::viewKeyStruct::elementCenterString );
  arrayView1d< real64 const > const & elemVolume =
    subRegion.getReference< array1d< real64 > >(
      CellBlock::viewKeyStruct::elementVolumeString );
  arrayView1d< R1Tensor const > const & elemPerm =
    subRegion.getReference< array1d< R1Tensor > >(
      SinglePhaseBase::viewKeyStruct::permeabilityString );

  // get the cell-centered depth
  arrayView1d< real64 const > const & elemGravCoef =
    subRegion.getReference< array1d< real64 > >(
      SinglePhaseBase::viewKeyStruct::gravityCoefString );

  // get the fluid data
  arrayView2d< real64 const > const & elemDens = fluid.density();
  arrayView2d< real64 const > const & dElemDens_dp = fluid.dDensity_dPressure();

  // assemble the residual and Jacobian element by element
  // in this loop we assemble both equation types: mass conservation in the elements and constraints at the faces
  using KERNEL_POLICY = parallelDevicePolicy< 32 >;
  forAll< KERNEL_POLICY >( subRegion.size(), [=] GEOSX_DEVICE( localIndex const ei ) {
    // transmissibility matrix
    stackArray2d< real64, NF * NF > transMatrix( NF, NF );

    real64 const perm[3] = { elemPerm[ei][0], elemPerm[ei][1], elemPerm[ei][2] };

    // recompute the local transmissibility matrix at each iteration
    // we can decide later to precompute transMatrix if needed
    HybridFVMInnerProduct::QTPFACellInnerProductKernel::Compute< NF >(
      nodePosition,
      faceToNodes,
      elemToFaces[ei],
      elemCenter[ei],
      elemVolume[ei],
      perm,
      2,
      lengthTolerance,
      transMatrix );

    // perform flux assembly in this element
    SinglePhaseHybridFVMKernels::AssemblerKernel::Compute< NF >( er,
                                                                 esr,
                                                                 ei,
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
                                                                 mobility,
                                                                 dMobility_dp,
                                                                 elemDofNumber,
                                                                 elemGhostRank[ei],
                                                                 rankOffset,
                                                                 dt,
                                                                 transMatrix,
                                                                 localMatrix,
                                                                 localRhs );
  } );
}

#define INST_AssembleKernelHelper( NF )                                  \
  template void AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >(   \
    arrayView1d< real64 const > const & facePres,                        \
    arrayView1d< real64 const > const & dFacePres,                       \
    arrayView1d< real64 const > const & faceGravCoef,                    \
    arraySlice1d< localIndex const > const & elemToFaces,                \
    real64 const & elemPres,                                             \
    real64 const & dElemPres,                                            \
    real64 const & elemGravDepth,                                        \
    real64 const & elemDens,                                             \
    real64 const & dElemDens_dp,                                         \
    arraySlice2d< real64 const > const & transMatrix,                    \
    arraySlice1d< real64 > const & oneSidedVolFlux,                      \
    arraySlice1d< real64 > const & dOneSidedVolFlux_dp,                  \
    arraySlice2d< real64 > const & dOneSidedVolFlux_dfp );               \
  template void AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >( \
    arrayView1d< globalIndex const > const & faceDofNumber,              \
    arraySlice1d< localIndex const > const & elemToFaces,                \
    globalIndex const elemDofNumber,                                     \
    globalIndex const rankOffset,                                        \
    arraySlice1d< real64 const > const & oneSidedVolFlux,                \
    arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,            \
    arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,           \
    arraySlice1d< real64 const > const & upwMobility,                    \
    arraySlice1d< real64 const > const & dUpwMobility_dp,                \
    arraySlice1d< globalIndex const > const & upwDofNumber,              \
    real64 const & dt,                                                   \
    CRSMatrixView< real64, globalIndex const > const & localMatrix,      \
    arrayView1d< real64 > const & localRhs );                            \
  template void AssemblerKernelHelper::AssembleConstraints< NF >(        \
    arrayView1d< globalIndex const > const & faceDofNumber,              \
    arrayView1d< integer const > const & faceGhostRank,                  \
    arraySlice1d< localIndex const > const & elemToFaces,                \
    globalIndex const elemDofNumber,                                     \
    globalIndex const rankOffset,                                        \
    arraySlice1d< real64 const > const & oneSidedVolFlux,                \
    arraySlice1d< real64 const > const & dOneSidedVolFlux_dp,            \
    arraySlice2d< real64 const > const & dOneSidedVolFlux_dfp,           \
    CRSMatrixView< real64, globalIndex const > const & localMatrix,      \
    arrayView1d< real64 > const & localRhs )

INST_AssembleKernelHelper( 4 );
INST_AssembleKernelHelper( 5 );
INST_AssembleKernelHelper( 6 );

#undef INST_AssembleKernelHelper

#define INST_FluxKernel( NF )                                                        \
  template void FluxKernel::Launch< NF >(                                            \
    localIndex er,                                                                   \
    localIndex esr,                                                                  \
    CellElementSubRegion const & subRegion,                                          \
    constitutive::SingleFluidBase const & fluid,                                     \
    SortedArrayView< localIndex const > const & regionFilter,                        \
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition, \
    arrayView2d< localIndex const > const & elemRegionList,                          \
    arrayView2d< localIndex const > const & elemSubRegionList,                       \
    arrayView2d< localIndex const > const & elemList,                                \
    ArrayOfArraysView< localIndex const > const & faceToNodes,                       \
    arrayView1d< globalIndex const > const & faceDofNumber,                          \
    arrayView1d< integer const > const & faceGhostRank,                              \
    arrayView1d< real64 const > const & facePres,                                    \
    arrayView1d< real64 const > const & dFacePres,                                   \
    arrayView1d< real64 const > const & faceGravCoef,                                \
    ElementView< arrayView1d< real64 const > > const & mobility,                     \
    ElementView< arrayView1d< real64 const > > const & dMobility_dp,                 \
    ElementView< arrayView1d< globalIndex const > > const & elemDofNumber,           \
    localIndex const rankOffset,                                                     \
    real64 const lengthTolerance,                                                    \
    real64 const dt,                                                                 \
    CRSMatrixView< real64, globalIndex const > const & localMatrix,                  \
    arrayView1d< real64 > const & localRhs )

INST_FluxKernel( 4 );
INST_FluxKernel( 5 );
INST_FluxKernel( 6 );

#undef INST_FluxKernel

}  // namespace SinglePhaseHybridFVMKernels

}  // namespace geosx
