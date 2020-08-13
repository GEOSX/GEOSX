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

// Source includes
#include "managers/initialization.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"
#include "codingUtilities/UnitTestUtilities.hpp"

// TPL includes
#include <gtest/gtest.h>

using namespace geosx;
using namespace geosx::SinglePhaseHybridFVMKernels;
using namespace geosx::HybridFVMInnerProduct;
using namespace geosx::testing;

static real64 constexpr relTol = 1e-5;
static localIndex constexpr NF = 4;  // we consider a tetrahedron in this file

// helpers

void
updateDensity( real64 const & refPres,
               real64 const & elemPres,
               real64 const & dElemPres,
               real64 & elemDens,
               real64 & dElemDens_dp )
{
  // we assume (very) compressible flow to catch wrong derivatives
  real64 const compressibility = 1e-3;
  real64 const refDens = 1000;

  // hard-coded relationship between pressure and density
  elemDens = refDens * exp( compressibility * ( elemPres + dElemPres - refPres ) );
  dElemDens_dp = compressibility * elemDens;
}

void
updateUpwindedMobilities( globalIndex const elemDofNumber,
                          real64 const & elemDens,
                          real64 const & dElemDens_dp,
                          arraySlice1d< real64 > const & upwMobility,
                          arraySlice1d< real64 > const & dUpwMobility_dp,
                          arraySlice1d< globalIndex > const & upwDofNumber )
{
  // we assume that viscosity is independent of pressure
  real64 const elemVisc = 0.001;
  real64 elemMobility = elemDens / elemVisc;
  real64 dElemMobility_dp = dElemDens_dp / elemVisc;

  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    upwMobility( ifaceLoc ) = elemMobility;
    dUpwMobility_dp( ifaceLoc ) = dElemMobility_dp;
    upwDofNumber( ifaceLoc ) = elemDofNumber;
  }
}

// in this function, we set up a 1-element problem in a tetrahedron
// the QTPFA transmissibility matrix comes from testHybridFVMInnerProducts
void
setupProblemForTetra( array1d< localIndex > & elemToFaces,
                      array1d< real64 > & facePres,
                      array1d< real64 > & dFacePres,
                      array1d< real64 > & faceGravCoef,
                      real64 & refPres,
                      real64 & elemPres,
                      real64 & dElemPres,
                      real64 & elemGravCoef,
                      real64 & elemDens,
                      real64 & dElemDens_dp,
                      arraySlice2d< real64 > const & transMatrix )
{
  facePres.resize( NF );
  dFacePres.resize( NF );
  facePres( 0 ) = 1e5;
  dFacePres( 0 ) = 0;
  facePres( 1 ) = 2e5;
  dFacePres( 1 ) = 0;
  facePres( 2 ) = 1e4;
  dFacePres( 2 ) = 0;
  facePres( 3 ) = 5e4;
  dFacePres( 3 ) = 0;

  elemGravCoef = 5e1;
  faceGravCoef.resize( NF );
  faceGravCoef( 0 ) = 1e1;
  faceGravCoef( 1 ) = 1e2;
  faceGravCoef( 2 ) = 3e1;
  faceGravCoef( 3 ) = 7e1;

  elemToFaces.resize( NF );
  elemToFaces( 0 ) = 0;
  elemToFaces( 1 ) = 1;
  elemToFaces( 2 ) = 2;
  elemToFaces( 3 ) = 3;

  elemPres = 1.2e5;
  refPres = elemPres;
  dElemPres = 0;
  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );

  // the transmissibility matrix comes from the HybridFVMInnerProduct unit tests
  transMatrix( 0, 0 ) = 5.25e-12;
  transMatrix( 0, 1 ) = 3.75e-12;
  transMatrix( 0, 2 ) = 2.25e-12;
  transMatrix( 0, 3 ) = 3.75e-12;

  transMatrix( 1, 0 ) = 3.75e-12;
  transMatrix( 1, 1 ) = 12.75e-12;
  transMatrix( 1, 2 ) = -5.25e-12;
  transMatrix( 1, 3 ) = 3.75e-12;

  transMatrix( 2, 0 ) = 2.25e-12;
  transMatrix( 2, 1 ) = -5.25e-12;
  transMatrix( 2, 2 ) = 18.75e-12;
  transMatrix( 2, 3 ) = -0.75e-12;

  transMatrix( 3, 0 ) = 3.75e-12;
  transMatrix( 3, 1 ) = 3.75e-12;
  transMatrix( 3, 2 ) = -0.75e-12;
  transMatrix( 3, 3 ) = 8.25e-12;
}

void
setupMatrixAndRhsForTetra( globalIndex & elemDofNumber,
                           array1d< globalIndex > & faceDofNumber,
                           array1d< integer > & faceGhostRank,
                           CRSMatrix< real64, globalIndex > & matrix,
                           CRSMatrix< real64, globalIndex > & matrixPerturb,
                           CRSMatrix< real64, globalIndex > & matrixFD,
                           array1d< real64 > & rhs,
                           array1d< real64 > & rhsPerturb )
{
  elemDofNumber = 0;
  faceDofNumber.resize( NF );
  faceDofNumber( 0 ) = 1;
  faceDofNumber( 1 ) = 2;
  faceDofNumber( 2 ) = 3;
  faceDofNumber( 3 ) = 4;

  faceGhostRank.resize( NF );
  faceGhostRank.setValues< serialPolicy >( -1 );

  matrix.resize( NF + 1, NF + 1, NF + 1 );
  matrixPerturb.resize( NF + 1, NF + 1, NF + 1 );
  matrixFD.resize( NF + 1, NF + 1, NF + 1 );

  rhs.resize( NF + 1 );
  rhsPerturb.resize( NF + 1 );

  // the matrices are full for the one-cell problem
  for( globalIndex i = 0; i < NF + 1; ++i )
  {
    for( globalIndex j = 0; j < NF + 1; ++j )
    {
      matrix.insertNonZero( i, j, 1 );
      matrixPerturb.insertNonZero( i, j, 1 );
      matrixFD.insertNonZero( i, j, 1 );
    }
  }
}

TEST( SinglePhaseHybridFVMKernels, assembleConstraints )
{
  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );

  //////////////////////
  // 1) problem setup //
  //////////////////////

  array1d< localIndex > elemToFaces;
  array1d< real64 > facePres;
  array1d< real64 > dFacePres;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  real64 elemDens;
  real64 dElemDens_dp;
  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  setupProblemForTetra( elemToFaces,
                        facePres,
                        dFacePres,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemDens,
                        dElemDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  array1d< integer > faceGhostRank;
  CRSMatrix< real64, globalIndex > jacobian;
  CRSMatrix< real64, globalIndex > jacobianPerturb;
  CRSMatrix< real64, globalIndex > jacobianFD;
  array1d< real64 > rhs;
  array1d< real64 > rhsPerturb;

  setupMatrixAndRhsForTetra( elemDofNumber,
                             faceDofNumber,
                             faceGhostRank,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF * NF > dOneSidedVolFlux_dfp( NF, NF );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

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

  jacobianFD.setValues< parallelHostPolicy >( 0.0 );
  jacobian.setValues< parallelHostPolicy >( 0.0 );
  rhs.setValues< parallelHostPolicy >( 0 );

  AssemblerKernelHelper::AssembleConstraints< NF >( faceDofNumber,
                                                    faceGhostRank,
                                                    elemToFaces,
                                                    elemDofNumber,
                                                    0,
                                                    oneSidedVolFlux,
                                                    dOneSidedVolFlux_dp,
                                                    dOneSidedVolFlux_dfp,
                                                    jacobian.toViewConstSizes(),
                                                    rhs.toView() );

  ///////////////////////////////////////////////////////////////////////////////////
  // 3) Compute finite-difference derivatives with respect to the element pressure //
  ///////////////////////////////////////////////////////////////////////////////////

  real64 const dElemPresPerturb =
    perturbParameter * ( elemPres + perturbParameter );
  // we need to update density to account for the perturbation
  updateDensity( refPres, elemPres, dElemPresPerturb, elemDens, dElemDens_dp );

  LvArray::tensorOps::fill< NF >( oneSidedVolFlux, 0 );
  LvArray::tensorOps::fill< NF >( dOneSidedVolFlux_dp, 0 );
  LvArray::tensorOps::fill< NF, NF >( dOneSidedVolFlux_dfp, 0 );

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPresPerturb,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianPerturb.setValues< parallelHostPolicy >( 0.0 );
  rhsPerturb.setValues< parallelHostPolicy >( 0.0 );

  AssemblerKernelHelper::AssembleConstraints< NF >(
    faceDofNumber,
    faceGhostRank,
    elemToFaces,
    elemDofNumber,
    0,
    oneSidedVolFlux,
    dOneSidedVolFlux_dp,
    dOneSidedVolFlux_dfp,
    jacobianPerturb.toViewConstSizes(),
    rhsPerturb.toView() );

  for( localIndex row = 0; row < rhs.size(); ++row )
  {
    real64 const dR_dp = ( rhsPerturb[row] - rhs[row] ) / dElemPresPerturb;
    if( std::fabs( dR_dp ) > 0.0 )
    {
      jacobianFD.addToRow< serialAtomic >( row, &elemDofNumber, &dR_dp, 1 );
    }
  }

  /////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the face pressures //
  /////////////////////////////////////////////////////////////////////////////////

  // we need to revert the density to its initial value (it is independent of face pressure)
  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    dFacePres.setValues< serialPolicy >( 0 );
    dFacePres[ifaceLoc] =
      perturbParameter * ( facePres[ifaceLoc] + perturbParameter );

    LvArray::tensorOps::fill< NF >( oneSidedVolFlux, 0 );
    LvArray::tensorOps::fill< NF >( dOneSidedVolFlux_dp, 0 );
    LvArray::tensorOps::fill< NF, NF >( dOneSidedVolFlux_dfp, 0 );

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

    jacobianPerturb.setValues< parallelHostPolicy >( 0.0 );
    rhsPerturb.setValues< parallelHostPolicy >( 0.0 );

    AssemblerKernelHelper::AssembleConstraints< NF >(
      faceDofNumber,
      faceGhostRank,
      elemToFaces,
      elemDofNumber,
      0,
      oneSidedVolFlux,
      dOneSidedVolFlux_dp,
      dOneSidedVolFlux_dfp,
      jacobianPerturb.toViewConstSizes(),
      rhsPerturb.toView() );

    for( localIndex row = 0; row < rhs.size(); ++row )
    {
      real64 const dR_dfp = ( rhsPerturb[row] - rhs[row] ) / dFacePres[ifaceLoc];
      if( std::fabs( dR_dfp ) > 0.0 )
      {
        jacobianFD.addToRow< serialAtomic >( row, &faceDofNumber[ifaceLoc], &dR_dfp, 1 );
      }
    }
  }

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

TEST( SinglePhaseHybridFVMKernels, assembleOneSidedMassFluxes )
{
  real64 const perturbParameter = sqrt( std::numeric_limits< real64 >::epsilon() );
  real64 const dt = 1e3;

  //////////////////////
  // 1) problem setup //
  //////////////////////

  array1d< localIndex > elemToFaces;
  array1d< real64 > facePres;
  array1d< real64 > dFacePres;
  array1d< real64 > faceGravCoef;
  real64 refPres;
  real64 elemPres;
  real64 dElemPres;
  real64 elemGravCoef;
  real64 elemDens;
  real64 dElemDens_dp;
  stackArray2d< real64, NF * NF > transMatrix( NF, NF );

  setupProblemForTetra( elemToFaces,
                        facePres,
                        dFacePres,
                        faceGravCoef,
                        refPres,
                        elemPres,
                        dElemPres,
                        elemGravCoef,
                        elemDens,
                        dElemDens_dp,
                        transMatrix );

  globalIndex elemDofNumber;
  array1d< globalIndex > faceDofNumber;
  array1d< integer > faceGhostRank;
  CRSMatrix< real64, globalIndex > jacobian;
  CRSMatrix< real64, globalIndex > jacobianPerturb;
  CRSMatrix< real64, globalIndex > jacobianFD;
  array1d< real64 > rhs;
  array1d< real64 > rhsPerturb;

  setupMatrixAndRhsForTetra( elemDofNumber,
                             faceDofNumber,
                             faceGhostRank,
                             jacobian,
                             jacobianPerturb,
                             jacobianFD,
                             rhs,
                             rhsPerturb );

  stackArray1d< real64, NF > upwMobility( NF );
  stackArray1d< real64, NF > dUpwMobility_dp( NF );
  stackArray1d< globalIndex, NF > upwDofNumber( NF );

  // no upwinding to do since we are in a one-cell problem
  updateUpwindedMobilities( elemDofNumber,
                            elemDens,
                            dElemDens_dp,
                            upwMobility,
                            dUpwMobility_dp,
                            upwDofNumber );

  stackArray1d< real64, NF > oneSidedVolFlux( NF );
  stackArray1d< real64, NF > dOneSidedVolFlux_dp( NF );
  stackArray2d< real64, NF * NF > dOneSidedVolFlux_dfp( NF, NF );

  ///////////////////////////////////////
  // 2) Compute analytical derivatives //
  ///////////////////////////////////////

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

  jacobianFD.setValues< parallelHostPolicy >( 0.0 );
  jacobian.setValues< parallelHostPolicy >( 0.0 );
  rhs.setValues< parallelHostPolicy >( 0 );

  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >(
    faceDofNumber,
    elemToFaces,
    elemDofNumber,
    0,
    oneSidedVolFlux,
    dOneSidedVolFlux_dp,
    dOneSidedVolFlux_dfp,
    upwMobility,
    dUpwMobility_dp,
    upwDofNumber,
    dt,
    jacobian.toViewConstSizes(),
    rhs.toView() );

  ///////////////////////////////////////////////////////////////////////////////////
  // 3) Compute finite-difference derivatives with respect to the element pressure //
  ///////////////////////////////////////////////////////////////////////////////////

  real64 const dElemPresPerturb =
    perturbParameter * ( elemPres + perturbParameter );
  // we need to update density and mobility to account for the perturbation
  updateDensity( refPres, elemPres, dElemPresPerturb, elemDens, dElemDens_dp );
  updateUpwindedMobilities( elemDofNumber,
                            elemDens,
                            dElemDens_dp,
                            upwMobility,
                            dUpwMobility_dp,
                            upwDofNumber );

  LvArray::tensorOps::fill< NF >( oneSidedVolFlux, 0 );
  LvArray::tensorOps::fill< NF >( dOneSidedVolFlux_dp, 0 );
  LvArray::tensorOps::fill< NF, NF >( dOneSidedVolFlux_dfp, 0 );

  AssemblerKernelHelper::ComputeOneSidedVolFluxes< NF >( facePres,
                                                         dFacePres,
                                                         faceGravCoef,
                                                         elemToFaces,
                                                         elemPres,
                                                         dElemPresPerturb,
                                                         elemGravCoef,
                                                         elemDens,
                                                         dElemDens_dp,
                                                         transMatrix,
                                                         oneSidedVolFlux,
                                                         dOneSidedVolFlux_dp,
                                                         dOneSidedVolFlux_dfp );

  jacobianPerturb.setValues< parallelHostPolicy >( 0.0 );
  rhsPerturb.setValues< serialPolicy >( 0.0 );

  AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >(
    faceDofNumber,
    elemToFaces,
    elemDofNumber,
    0,
    oneSidedVolFlux,
    dOneSidedVolFlux_dp,
    dOneSidedVolFlux_dfp,
    upwMobility,
    dUpwMobility_dp,
    upwDofNumber,
    dt,
    jacobianPerturb.toViewConstSizes(),
    rhsPerturb.toView() );

  for( localIndex row = 0; row < rhs.size(); ++row )
  {
    real64 const dR_dp = ( rhsPerturb[row] - rhs[row] ) / dElemPresPerturb;
    if( std::fabs( dR_dp ) > 0.0 )
    {
      jacobianFD.addToRow< serialAtomic >( row, &elemDofNumber, &dR_dp, 1 );
    }
  }

  /////////////////////////////////////////////////////////////////////////////////
  // 4) Compute finite-difference derivatives with respect to the face pressures //
  /////////////////////////////////////////////////////////////////////////////////

  updateDensity( refPres, elemPres, dElemPres, elemDens, dElemDens_dp );
  updateUpwindedMobilities( elemDofNumber,
                            elemDens,
                            dElemDens_dp,
                            upwMobility,
                            dUpwMobility_dp,
                            upwDofNumber );
  for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
  {
    dFacePres.setValues< serialPolicy >( 0 );
    dFacePres[ifaceLoc] =
      perturbParameter * ( facePres[ifaceLoc] + perturbParameter );

    LvArray::tensorOps::fill< NF >( oneSidedVolFlux, 0 );
    LvArray::tensorOps::fill< NF >( dOneSidedVolFlux_dp, 0 );
    LvArray::tensorOps::fill< NF, NF >( dOneSidedVolFlux_dfp, 0 );

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

    jacobianPerturb.setValues< parallelHostPolicy >( 0.0 );
    rhsPerturb.setValues< parallelHostPolicy >( 0.0 );

    AssemblerKernelHelper::AssembleOneSidedMassFluxes< NF >(
      faceDofNumber,
      elemToFaces,
      elemDofNumber,
      0,
      oneSidedVolFlux,
      dOneSidedVolFlux_dp,
      dOneSidedVolFlux_dfp,
      upwMobility,
      dUpwMobility_dp,
      upwDofNumber,
      dt,
      jacobianPerturb.toViewConstSizes(),
      rhsPerturb.toView() );

    for( localIndex row = 0; row < rhs.size(); ++row )
    {
      real64 const dR_dfp = ( rhsPerturb[row] - rhs[row] ) / dFacePres[ifaceLoc];
      if( std::fabs( dR_dfp ) > 0.0 )
      {
        jacobianFD.addToRow< serialAtomic >( row, &faceDofNumber[ifaceLoc], &dR_dfp, 1 );
      }
    }
  }

  compareLocalMatrices( jacobian.toViewConst(), jacobianFD.toViewConst(), relTol );
}

int
main( int argc, char ** argv )
{
  ::testing::InitGoogleTest( &argc, argv );

  geosx::basicSetup( argc, argv );

  int const result = RUN_ALL_TESTS();

  geosx::basicCleanup();

  return result;
}
