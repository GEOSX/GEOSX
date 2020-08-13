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
 * @file SinglePhaseFVMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

namespace geosx
{
namespace SinglePhaseFVMKernels
{
/******************************** FluxKernel ********************************/

/**
 * @struct Structure to contain helper functions for the FluxKernel struct.
 */
struct FluxKernelHelper
{
  /**
   * @tparam INTEGRATION_OPTION This specifies the choice of integration rule for the aperture term
   *         in a lubrication permeability.
   * @param[in] aper0 The beginning of step aperture
   * @param[in] aper The current approximation to the end of step aperture
   * @param[out] aperTerm The resulting
   * @param[out] dAperTerm_dAper
   *
   * Typically in lubrication theory, the permeabilty involves a \f$ aperture^3 \f$ term. The
   * flow residual equation assumes a constant value for all parameters over \f$ dt \f$, which
   * may introduce significant errors given the highly nonlinear nature of the cubic aperture term.
   * The template parameter provides options:
   *  - (0) Forward Euler. This results in no non-linearity since the beginning of step aperture
   *    does not change.
   *  - (1) Exact/Simpson's Rule. This is the result of taking
   *    \f$ \int_{0}^{1} (aperture0 + (aperture-aperture0)x)x^3 dx \f$. This results
   *    in a cubic non-linearity in the resulting set of equations.
   *  .
   *  @note The use of option (1) does not imply that the time integration of the residual
   *        equation is exact, or applying Simpson's Rule. It only means that the integral of
   *        the cubic aperture term in the permeablity is exact. All other components of the
   *        residual equation are assumed constant over the integral, or use a backward
   *        Euler as the case may be. Also, we omit the use of a backward Euler option as
   *        it offers no benefit over the exact integration.
   */
  template< int INTEGRATION_OPTION >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  apertureForPermeablityCalculation( real64 const aper0,
                                     real64 const aper,
                                     real64 & aperTerm,
                                     real64 & dAperTerm_dAper );
};

template<>
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
FluxKernelHelper::apertureForPermeablityCalculation< 0 >( real64 const aper0,
                                                          real64 const,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = aper0 * aper0 * aper0;
  dAperTerm_dAper = 0.0;
}

template<>
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
FluxKernelHelper::apertureForPermeablityCalculation< 1 >( real64 const aper0,
                                                          real64 const aper,
                                                          real64 & aperTerm,
                                                          real64 & dAperTerm_dAper )
{
  aperTerm = 0.25 *
    ( aper0 * aper0 * aper0 + aper0 * aper0 * aper + aper0 * aper * aper +
      aper * aper * aper );

  dAperTerm_dAper = 0.25 * ( aper0 * aper0 + 2 * aper0 * aper + 3 * aper * aper );

  //printf( "aper0, aper, Kf = %4.2e, %4.2e, %4.2e \n", aper0, aper, aperTerm );
}

template<>
inline void
FluxKernelHelper::apertureForPermeablityCalculation< 2 >(
  real64 const,
  real64 const aper,
  real64 & aperTerm,
  real64 & dAperTerm_dAper )
{
  aperTerm = aper * aper * aper;
  dAperTerm_dAper = 3.0 * aper * aper;
}

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
  using ElementView =
    typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

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
   * @param[out] jacobian The linear system matrix
   * @param[out] residual The linear system residual
   */
  template< typename STENCIL_TYPE >
  static void
  Launch( STENCIL_TYPE const & stencil,
          real64 const dt,
          globalIndex const rankOffset,
          ElementView< arrayView1d< globalIndex const > > const & dofNumber,
          ElementView< arrayView1d< integer const > > const & ghostRank,
          ElementView< arrayView1d< real64 const > > const & pres,
          ElementView< arrayView1d< real64 const > > const & dPres,
          ElementView< arrayView1d< real64 const > > const & gravCoef,
          ElementView< arrayView2d< real64 const > > const & dens,
          ElementView< arrayView2d< real64 const > > const & dDens_dPres,
          ElementView< arrayView1d< real64 const > > const & mob,
          ElementView< arrayView1d< real64 const > > const & dMob_dPres,
          ElementView< arrayView1d< real64 const > > const & aperture0,
          ElementView< arrayView1d< real64 const > > const & aperture,
          ElementView< arrayView1d< R1Tensor const > > const & transTMultiplier,
          R1Tensor const gravityVector,
          real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
          ElementView< arrayView1d< real64 const > > const & s,
          ElementView< arrayView1d< real64 const > > const & dSdAper,
#endif
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs,
          CRSMatrixView< real64, localIndex const > const & dR_dAper );

  /**
   * @brief Compute flux and its derivatives for a given connection
   *
   * This is a general version that assumes different element regions.
   * See below for a specialized version for fluxes within a region.
   */
  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementView< arrayView1d< real64 const > > const & pres,
           ElementView< arrayView1d< real64 const > > const & dPres,
           ElementView< arrayView1d< real64 const > > const & gravCoef,
           ElementView< arrayView2d< real64 const > > const & dens,
           ElementView< arrayView2d< real64 const > > const & dDens_dPres,
           ElementView< arrayView1d< real64 const > > const & mob,
           ElementView< arrayView1d< real64 const > > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian );

  /**
   * @brief Compute flux and its derivatives for a given connection
   *.
   * This is a specialized version for fluxes within the same region.
   * See above for a general version.
   */
  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const stencilSize,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const &,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           arrayView1d< real64 const > const & pres,
           arrayView1d< real64 const > const & dPres,
           arrayView1d< real64 const > const & gravCoef,
           arrayView2d< real64 const > const & dens,
           arrayView2d< real64 const > const & dDens_dPres,
           arrayView1d< real64 const > const & mob,
           arrayView1d< real64 const > const & dMob_dPres,
           real64 const dt,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian );

  /**
   * @brief Compute flux and its derivatives for a given multi-element connector.
   *
   * This is a specialized version that flux in a single region, and uses
   * element pairing instead of a proper junction.
   */
  GEOSX_HOST_DEVICE
  static void
  ComputeJunction(
    localIndex const numFluxElems,
    arraySlice1d< localIndex const > const & stencilElementIndices,
    arraySlice1d< real64 const > const & stencilWeights,
    arrayView1d< real64 const > const & pres,
    arrayView1d< real64 const > const & dPres,
    arrayView1d< real64 const > const & gravCoef,
    arrayView2d< real64 const > const & dens,
    arrayView2d< real64 const > const & dDens_dPres,
    arrayView1d< real64 const > const & mob,
    arrayView1d< real64 const > const & dMob_dPres,
    arrayView1d< real64 const > const & aperture0,
    arrayView1d< real64 const > const & aperture,
    real64 const meanPermCoeff,
#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    arrayView1d< real64 const > const & GEOSX_GEOSX_UNUSED_PARAM( s ),
    arrayView1d< real64 const > const & GEOSX_GEOSX_UNUSED_PARAM( dSdAper ),
#endif
    real64 const dt,
    arraySlice1d< real64 > const & flux,
    arraySlice2d< real64 > const & fluxJacobian,
    arraySlice2d< real64 > const & dFlux_dAperture );
};

struct FaceDirichletBCKernel
{
  template< typename VIEWTYPE >
  using ElementView = FluxKernel::ElementView< VIEWTYPE >;

  template< typename FLUID_WRAPPER >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
    arraySlice1d< localIndex const > const & seri,
    arraySlice1d< localIndex const > const & sesri,
    arraySlice1d< localIndex const > const & sefi,
    arraySlice1d< real64 const > const & trans,
    ElementView< arrayView1d< real64 const > > const & pres,
    ElementView< arrayView1d< real64 const > > const & dPres,
    ElementView< arrayView1d< real64 const > > const & gravCoef,
    ElementView< arrayView2d< real64 const > > const & dens,
    ElementView< arrayView2d< real64 const > > const & dDens_dPres,
    ElementView< arrayView1d< real64 const > > const & mob,
    ElementView< arrayView1d< real64 const > > const & dMob_dPres,
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

    localIndex const er = seri[Order::ELEM];
    localIndex const esr = sesri[Order::ELEM];
    localIndex const ei = sefi[Order::ELEM];
    localIndex const kf = sefi[Order::FACE];

    // Get flow quantities on the elem/face
    real64 faceDens, faceVisc;
    fluidWrapper.Compute( presFace[kf], faceDens, faceVisc );

    mobility[Order::ELEM] = mob[er][esr][ei];
    SinglePhaseBaseKernels::MobilityKernel::Compute( faceDens,
                                                     faceVisc,
                                                     mobility[Order::FACE] );

    dMobility_dP[Order::ELEM] = dMob_dPres[er][esr][ei];
    dMobility_dP[Order::FACE] = 0.0;

    // Compute average density
    real64 const densMean = 0.5 * ( dens[er][esr][ei][0] + faceDens );
    real64 const dDens_dP = 0.5 * dDens_dPres[er][esr][ei][0];

    // Evaluate potential difference
    real64 const potDif = trans[Order::ELEM] *
        ( pres[er][esr][ei] + dPres[er][esr][ei] +
          densMean * gravCoef[er][esr][ei] ) +
      trans[Order::FACE] * ( presFace[kf] + densMean * gravCoefFace[kf] );

    real64 const dPotDif_dP =
      trans[Order::ELEM] * ( 1.0 + dDens_dP * gravCoef[er][esr][ei] );

    // Upwind mobility
    localIndex const k_up = ( potDif >= 0 ) ? Order::ELEM : Order::FACE;

    flux = dt * mobility[k_up] * potDif;
    dFlux_dP = dt * ( mobility[k_up] * dPotDif_dP + dMobility_dP[k_up] * potDif );
  }

  template< typename FLUID_WRAPPER >
  static void
  Launch( BoundaryStencil::IndexContainerViewConstType const & seri,
          BoundaryStencil::IndexContainerViewConstType const & sesri,
          BoundaryStencil::IndexContainerViewConstType const & sefi,
          BoundaryStencil::WeightContainerViewConstType const & trans,
          ElementView< arrayView1d< integer const > > const & ghostRank,
          ElementView< arrayView1d< globalIndex const > > const & dofNumber,
          globalIndex const rankOffset,
          ElementView< arrayView1d< real64 const > > const & pres,
          ElementView< arrayView1d< real64 const > > const & dPres,
          ElementView< arrayView1d< real64 const > > const & gravCoef,
          ElementView< arrayView2d< real64 const > > const & dens,
          ElementView< arrayView2d< real64 const > > const & dDens_dPres,
          ElementView< arrayView1d< real64 const > > const & mob,
          ElementView< arrayView1d< real64 const > > const & dMob_dPres,
          arrayView1d< real64 const > const & presFace,
          arrayView1d< real64 const > const & gravCoefFace,
          FLUID_WRAPPER const & fluidWrapper,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    forAll< parallelDevicePolicy<> >(
      seri.size( 0 ),
      [=] GEOSX_HOST_DEVICE( localIndex const iconn ) {
        real64 flux, fluxJacobian;

        Compute( seri[iconn],
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

        localIndex const er = seri( iconn, BoundaryStencil::Order::ELEM );
        localIndex const esr = sesri( iconn, BoundaryStencil::Order::ELEM );
        localIndex const ei = sefi( iconn, BoundaryStencil::Order::ELEM );

        if( ghostRank[er][esr][ei] < 0 )
        {
          // Add to global residual/jacobian
          globalIndex const dofIndex = dofNumber[er][esr][ei];
          localIndex const localRow =
            LvArray::integerConversion< localIndex >( dofIndex - rankOffset );

          RAJA::atomicAdd( parallelDeviceAtomic {}, &localRhs[localRow], flux );
          localMatrix.addToRow< parallelDeviceAtomic >( localRow,
                                                        &dofIndex,
                                                        &fluxJacobian,
                                                        1 );
        }
      } );
  }
};

}  // namespace SinglePhaseFVMKernels

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVMKERNELS_HPP
