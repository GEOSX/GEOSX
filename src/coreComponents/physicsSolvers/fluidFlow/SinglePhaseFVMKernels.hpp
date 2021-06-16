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
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/FluxKernelsHelper.hpp"

namespace geosx
{

namespace SinglePhaseFVMKernels
{

using namespace FluxKernelsHelper;

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
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   * @param[in] stencil The stencil object.
   * @param[in] dt The timestep for the integration step.
   * @param[in] dofNumber The dofNumbers for each element
   * @param[in] pres The pressures in each element
   * @param[in] dPres The change in pressure for each element
   * @param[in] gravCoef The factor for gravity calculations (g*H)
   * @param[in] dens The material density in each element
   * @param[in] dDens_dPres The change in material density for each element
   * @param[in] mob The fluid mobility in each element
   * @param[in] dMob_dPres The derivative of mobility wrt pressure in each element
   * @param[in] permeability
   * @param[in] dPerm_dPres The derivative of permeability wrt pressure in each element
   * @param[in] transTMultiplier
   * @param[in] gravityVector
   * @param[out] localMatrix The linear system matrix
   * @param[out] localRhs The linear system residual
   */
  template< typename STENCILWRAPPER_TYPE >
  static void
  launch( STENCILWRAPPER_TYPE const & stencilWrapper,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & dens,
          ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
          ElementViewConst< arrayView1d< real64 const > > const & mob,
          ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & permeability,
          ElementViewConst< arrayView3d< real64 const > > const & dPerm_dPres,
          ElementViewConst< arrayView2d< real64 const > > const & GEOSX_UNUSED_PARAM( transTMultiplier ),
          R1Tensor const & GEOSX_UNUSED_PARAM ( gravityVector ),
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
//    constexpr localIndex maxNumFluxElems = STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX;
//    constexpr localIndex maxStencilSize  = STENCILWRAPPER_TYPE::MAX_STENCIL_SIZE;

    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & seri = stencilWrapper.getElementRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sesri = stencilWrapper.getElementSubRegionIndices();
    typename STENCILWRAPPER_TYPE::IndexContainerViewConstType const & sei = stencilWrapper.getElementIndices();


    forAll< parallelDevicePolicy<> >( stencilWrapper.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      localIndex const stencilSize = stencilWrapper.stencilSize( iconn );
      localIndex const numFluxElems = stencilWrapper.numPointsInFlux( iconn );

      // working arrays
      stackArray1d< globalIndex, STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX > dofColIndices( stencilSize );
      stackArray1d< real64, STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX > localFlux( numFluxElems );
      stackArray2d< real64, STENCILWRAPPER_TYPE::NUM_POINT_IN_FLUX * STENCILWRAPPER_TYPE::MAX_STENCIL_SIZE > localFluxJacobian( numFluxElems, stencilSize );

      // compute transmissibility
      real64 transmissiblity[2], dTrans_dPres[2];
      stencilWrapper.computeTransmissibility( iconn,
                                              permeability,
                                              dPerm_dPres,
                                              transmissiblity,
                                              dTrans_dPres );

      compute( stencilSize,
               seri[iconn],
               sesri[iconn],
               sei[iconn],
               transmissiblity,
               dTrans_dPres,
               pres,
               dPres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               dt,
               localFlux,
               localFluxJacobian );


      // extract DOF numbers
      for( localIndex i = 0; i < stencilSize; ++i )
      {
        dofColIndices[i] = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];

      }

      for( localIndex i = 0; i < numFluxElems; ++i )
      {

        if( ghostRank[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )] < 0 )
        {
          globalIndex const globalRow = dofNumber[seri( iconn, i )][sesri( iconn, i )][sei( iconn, i )];
          localIndex const localRow = LvArray::integerConversion< localIndex >( globalRow - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GT( localMatrix.numRows(), localRow );

          RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], localFlux[i] );
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( localRow,
                                                                            dofColIndices.data(),
                                                                            localFluxJacobian[i].dataIfContiguous(),
                                                                            stencilSize );

        }
      }

    } );
  }

  /**
   * @brief Compute flux and its derivatives for a given tpfa connector.
   *
   *
   */
  GEOSX_HOST_DEVICE
  static void
  compute( localIndex const numFluxElems,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           real64 const (&transmissibility)[2],
           real64 const (&dTrans_dPres)[2],
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & dens,
           ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
           ElementViewConst< arrayView1d< real64 const > > const & mob,
           ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian )
  {
    GEOSX_UNUSED_VAR( numFluxElems );

    real64 fluxVal = 0.0;
    real64 dFlux_dTrans[2] = {0.0, 0.0};
    real64 dFlux_dP[2] = {0.0, 0.0};

    computeSinglePhaseFlux( seri, sesri, sei,
                            transmissibility,
                            dTrans_dPres,
                            pres,
                            dPres,
                            gravCoef,
                            dens,
                            dDens_dPres,
                            mob,
                            dMob_dPres,
                            fluxVal,
                            dFlux_dP,
                            dFlux_dTrans );

    // populate local flux vector and derivatives
    flux[0] =  dt * fluxVal;
    flux[1] = -dt * fluxVal;


    for( localIndex ke = 0; ke < 2; ++ke )
    {
      fluxJacobian[0][ke] =  dt * dFlux_dP[ke];
      fluxJacobian[1][ke] = -dt * dFlux_dP[ke];
    }
  }
};

struct FaceDirichletBCKernel
{
  template< typename VIEWTYPE >
  using ElementViewConst = FluxKernel::ElementViewConst< VIEWTYPE >;

  template< typename FLUID_WRAPPER >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void compute( arraySlice1d< localIndex const > const & seri,
                       arraySlice1d< localIndex const > const & sesri,
                       arraySlice1d< localIndex const > const & sefi,
                       arraySlice1d< real64 const > const & trans,
                       ElementViewConst< arrayView1d< real64 const > > const & pres,
                       ElementViewConst< arrayView1d< real64 const > > const & dPres,
                       ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                       ElementViewConst< arrayView2d< real64 const > > const & dens,
                       ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                       ElementViewConst< arrayView1d< real64 const > > const & mob,
                       ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                       arrayView1d< real64 const > const & presFace,
                       arrayView1d< real64 const > const & gravCoefFace,
                       FLUID_WRAPPER const & fluidWrapper,
                       real64 const dt,
                       real64 & flux,
                       real64 & dFlux_dP )
  {
    using Order = BoundaryStencil::Order;
    localIndex constexpr numElems = BoundaryStencil::NUM_POINT_IN_FLUX;

    stackArray1d< real64, numElems > mobility( numElems );
    stackArray1d< real64, numElems > dMobility_dP( numElems );

    localIndex const er  = seri[ Order::ELEM ];
    localIndex const esr = sesri[ Order::ELEM ];
    localIndex const ei  = sefi[ Order::ELEM ];
    localIndex const kf  = sefi[ Order::FACE ];

    // Get flow quantities on the elem/face
    real64 faceDens, faceVisc;
    fluidWrapper.compute( presFace[kf], faceDens, faceVisc );

    mobility[Order::ELEM] = mob[er][esr][ei];
    SinglePhaseBaseKernels::MobilityKernel::compute( faceDens, faceVisc, mobility[Order::FACE] );

    dMobility_dP[Order::ELEM] = dMob_dPres[er][esr][ei];
    dMobility_dP[Order::FACE] = 0.0;

    // Compute average density
    real64 const densMean = 0.5 * ( dens[er][esr][ei][0] + faceDens );
    real64 const dDens_dP = 0.5 * dDens_dPres[er][esr][ei][0];

    // Evaluate potential difference
    real64 const potDif = trans[ Order::ELEM ] * ( pres[er][esr][ei] + dPres[er][esr][ei] - densMean * gravCoef[er][esr][ei] )
                          + trans[ Order::FACE ] * ( presFace[kf] - densMean * gravCoefFace[kf] );

    real64 const dPotDif_dP = trans[ Order::ELEM ] * ( 1.0 - dDens_dP * gravCoef[er][esr][ei] );

    // Upwind mobility
    localIndex const k_up = ( potDif >= 0 ) ? Order::ELEM : Order::FACE;

    flux = dt * mobility[k_up] * potDif;
    dFlux_dP = dt * ( mobility[k_up] * dPotDif_dP + dMobility_dP[k_up] * potDif );
  }

  template< typename FLUID_WRAPPER >
  static void launch( BoundaryStencil::IndexContainerViewConstType const & seri,
                      BoundaryStencil::IndexContainerViewConstType const & sesri,
                      BoundaryStencil::IndexContainerViewConstType const & sefi,
                      BoundaryStencil::WeightContainerViewConstType const & trans,
                      ElementViewConst< arrayView1d< integer const > > const & ghostRank,
                      ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
                      globalIndex const rankOffset,
                      ElementViewConst< arrayView1d< real64 const > > const & pres,
                      ElementViewConst< arrayView1d< real64 const > > const & dPres,
                      ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
                      ElementViewConst< arrayView2d< real64 const > > const & dens,
                      ElementViewConst< arrayView2d< real64 const > > const & dDens_dPres,
                      ElementViewConst< arrayView1d< real64 const > > const & mob,
                      ElementViewConst< arrayView1d< real64 const > > const & dMob_dPres,
                      arrayView1d< real64 const > const & presFace,
                      arrayView1d< real64 const > const & gravCoefFace,
                      FLUID_WRAPPER const & fluidWrapper,
                      real64 const dt,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs )
  {
    forAll< parallelDevicePolicy<> >( seri.size( 0 ), [=] GEOSX_HOST_DEVICE ( localIndex const iconn )
    {
      real64 flux, fluxJacobian;

      compute( seri[iconn],
               sesri[iconn],
               sefi[iconn],
               trans[iconn],
               pres,
               dPres,
               gravCoef,
               dens,
               dDens_dPres,
               mob,
               dMob_dPres,
               presFace,
               gravCoefFace,
               fluidWrapper,
               dt,
               flux,
               fluxJacobian );

      localIndex const er  = seri( iconn, BoundaryStencil::Order::ELEM );
      localIndex const esr = sesri( iconn, BoundaryStencil::Order::ELEM );
      localIndex const ei  = sefi( iconn, BoundaryStencil::Order::ELEM );

      if( ghostRank[er][esr][ei] < 0 )
      {
        // Add to global residual/jacobian
        globalIndex const dofIndex = dofNumber[er][esr][ei];
        localIndex const localRow = LvArray::integerConversion< localIndex >( dofIndex - rankOffset );

        RAJA::atomicAdd( parallelDeviceAtomic{}, &localRhs[localRow], flux );
        localMatrix.addToRow< parallelDeviceAtomic >( localRow, &dofIndex, &fluxJacobian, 1 );
      }
    } );
  }
};

} // namespace SinglePhaseFVMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
