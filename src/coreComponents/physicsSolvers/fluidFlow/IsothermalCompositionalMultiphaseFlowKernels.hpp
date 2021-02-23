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
 * @file IsothermalCompositionalMultiphaseFlowKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_ISOTHERMALCOMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_ISOTHERMALCOMPOSITIONALMULTIPHASEFLOWKERNELS_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace IsothermalCompositionalMultiphaseFlowKernels
{

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
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
           arraySlice2d< real64 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
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
          arrayView3d< real64 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
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
          arrayView3d< real64 > const & dPhaseVolFrac_dComp );
};

/******************************** PhaseMobilityKernel ********************************/

/**
 * @brief Functions to compute phase mobilities and derivatives from density, viscosity and relperm
 */
struct PhaseMobilityKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( arraySlice2d< real64 const > const & dCompFrac_dCompDens,
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
           arraySlice2d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
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
          arrayView3d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
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
          arrayView3d< real64 > const & dPhaseMob_dComp );
};

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( size, [=] ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( size, [=] ( localIndex const k )
    {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }
};

/******************************** AccumulationKernel ********************************/

/**
 * @brief Functions to assemble accumulation term contributions to residual and Jacobian
 */
struct AccumulationKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
    compute( localIndex const numPhases,
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
             real64 ( &localAccumJacobian )[NC][NC+1] );

  template< localIndex NC >
  static void
  launch( localIndex const numPhases,
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
          arrayView1d< real64 > const & localRhs );
};

/******************************** FluxKernel ********************************/

/**
 * @brief Functions to assemble flux term contributions to residual and Jacobian
 */
struct FluxKernel
{

  /**
   * @brief The type for element-based data. Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewConstAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( localIndex const stencilSize,
           localIndex const numPhases,
           arraySlice1d< localIndex const > const seri,
           arraySlice1d< localIndex const > const sesri,
           arraySlice1d< localIndex const > const sei,
           arraySlice1d< real64 const > const stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
           ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
           ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
           ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
           ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
           ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
           ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
           ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
           ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
           integer const capPressureFlag,
           real64 const dt,
           arraySlice1d< real64 > const localFlux,
           arraySlice2d< real64 > const localFluxJacobian );

  template< localIndex NC, typename STENCIL_TYPE >
  static void
  launch( localIndex const numPhases,
          STENCIL_TYPE const & stencil,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView1d< real64 const > > const & gravCoef,
          ElementViewConst< arrayView2d< real64 const > > const & phaseMob,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseMob_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseMob_dComp,
          ElementViewConst< arrayView2d< real64 const > > const & dPhaseVolFrac_dPres,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseVolFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & dCompFrac_dCompDens,
          ElementViewConst< arrayView3d< real64 const > > const & phaseDens,
          ElementViewConst< arrayView3d< real64 const > > const & dPhaseDens_dPres,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseDens_dComp,
          ElementViewConst< arrayView4d< real64 const > > const & phaseCompFrac,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCompFrac_dPres,
          ElementViewConst< arrayView5d< real64 const > > const & dPhaseCompFrac_dComp,
          ElementViewConst< arrayView3d< real64 const > > const & phaseCapPressure,
          ElementViewConst< arrayView4d< real64 const > > const & dPhaseCapPressure_dPhaseVolFrac,
          integer const capPressureFlag,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** VolumeBalanceKernel ********************************/

/**
 * @brief Functions to assemble volume balance contributions to residual and Jacobian
 */
struct VolumeBalanceKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  compute( real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const > const & phaseVolFrac,
           arraySlice1d< real64 const > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const > const & dPhaseVolFrac_dCompDens,
           real64 & localVolBalance,
           real64 * const localVolBalanceJacobian );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
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
          arrayView1d< real64 > const & localRhs );
};


/******************************** Kernel launch machinery ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernellaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernellaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    { lambda( std::integral_constant< T, 1 >() ); return; }
    case 2:
    { lambda( std::integral_constant< T, 2 >() ); return; }
    case 3:
    { lambda( std::integral_constant< T, 3 >() ); return; }
    case 4:
    { lambda( std::integral_constant< T, 4 >() ); return; }
    case 5:
    { lambda( std::integral_constant< T, 5 >() ); return; }
    default:
    { GEOSX_ERROR( "Unsupported number of components: " << value ); }
  }
}

} // namespace helpers

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector1( localIndex numComp, ARGS && ... args )
{
  internal::KernellaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS && ... args )
{
  internal::KernellaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    switch( numPhase )
    {
      case 2:
        { KERNELWRAPPER::template launch< NC(), 2 >( std::forward< ARGS >( args )... ); return; }
      case 3:
        { KERNELWRAPPER::template launch< NC(), 3 >( std::forward< ARGS >( args )... ); return; }
      default:
        { GEOSX_ERROR( "Unsupported number of phases: " << numPhase ); }
    }
  } );
}


} // namespace IsotherCompositionalMultiphaseFlowKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_ISOTHERMALCOMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
