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
 * @file CompositionalMultiphaseFlowKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP

#include "common/DataTypes.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
namespace CompositionalMultiphaseFlowKernels
{
/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
    arraySlice1d< real64 const > compDens,
    arraySlice1d< real64 const > dCompDens,
    arraySlice1d< real64 > compFrac,
    arraySlice2d< real64 > dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  Launch( localIndex const size,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const > const & compDens,
          arrayView2d< real64 const > const & dCompDens,
          arrayView2d< real64 > const & compFrac,
          arrayView3d< real64 > const & dCompFrac_dCompDens );
};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
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
    arraySlice2d< real64 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  Launch( localIndex const size,
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
  Launch( SortedArrayView< localIndex const > const & targetSet,
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
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
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
    arraySlice2d< real64 > const & dPhaseMob_dComp );

  template< localIndex NC, localIndex NP >
  static void
  Launch( localIndex const size,
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
  Launch( SortedArrayView< localIndex const > const & targetSet,
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
  Launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( size, [=]( localIndex const k ) {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.Update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  Launch( localIndex const size,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( size, [=]( localIndex const k ) {
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.Update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=]( localIndex const a ) {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.Update( k, q, pres[k] + dPres[k], temp, compFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename FLUID_WRAPPER >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          real64 const temp,
          arrayView2d< real64 const > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=]( localIndex const a ) {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.Update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }
};

/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  Launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const k ) {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.Update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE( localIndex const a ) {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.Update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  Launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const k ) {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.Update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  Launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE( localIndex const a ) {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.Update( k, q, phaseVolFrac[k] );
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
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
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
    real64 ( &localAccumJacobian )[NC][NC + 1] );

  template< localIndex NC >
  static void
  Launch( localIndex const numPhases,
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
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementView =
    typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< localIndex NC, localIndex NUM_ELEMS, localIndex MAX_STENCIL >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
    localIndex const stencilSize,
    localIndex const numPhases,
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
    arraySlice2d< real64 > const & localFluxJacobian );

  template< localIndex NC, typename STENCIL_TYPE >
  static void
  Launch(
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
    arrayView1d< real64 > const & localRhs );
};

/******************************** VolumeBalanceKernel ********************************/

/**
 * @brief Functions to assemble volume balance contributions to residual and Jacobian
 */
struct VolumeBalanceKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
    real64 const & volume,
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
  Launch( localIndex const size,
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
void
KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value,
                 "KernelLaunchSelectorCompSwitch: type should be integral" );

  switch( value )
  {
    case 1:
    {
      lambda( std::integral_constant< T, 1 >() );
      return;
    }
    case 2:
    {
      lambda( std::integral_constant< T, 2 >() );
      return;
    }
    case 3:
    {
      lambda( std::integral_constant< T, 3 >() );
      return;
    }
    case 4:
    {
      lambda( std::integral_constant< T, 4 >() );
      return;
    }
    case 5:
    {
      lambda( std::integral_constant< T, 5 >() );
      return;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported number of components: " << value );
    }
  }
}

}  // namespace internal

template< typename KERNELWRAPPER, typename... ARGS >
void
KernelLaunchSelector1( localIndex numComp, ARGS &&... args )
{
  internal::KernelLaunchSelectorCompSwitch( numComp, [&]( auto NC ) {
    KERNELWRAPPER::template Launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename... ARGS >
void
KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS &&... args )
{
  internal::KernelLaunchSelectorCompSwitch( numComp, [&]( auto NC ) {
    switch( numPhase )
    {
      case 2:
      {
        KERNELWRAPPER::template Launch< NC(), 2 >( std::forward< ARGS >( args )... );
        return;
      }
      case 3:
      {
        KERNELWRAPPER::template Launch< NC(), 3 >( std::forward< ARGS >( args )... );
        return;
      }
      default:
      {
        GEOSX_ERROR( "Unsupported number of phases: " << numPhase );
      }
    }
  } );
}

}  // namespace CompositionalMultiphaseFlowKernels

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOWKERNELS_HPP
