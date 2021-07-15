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
 * @file CompositionalMultiphaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP

#include "common/DataLayouts.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/fluid/layouts.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace CompositionalMultiphaseBaseKernels
{

using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

/******************************** ComponentFractionKernel ********************************/

/**
 * @brief Functions to compute component fractions from global component densities (mass or molar)
 */
struct ComponentFractionKernel
{
  template< localIndex NC >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > dCompDens,
           arraySlice1d< real64, compflow::USD_COMP - 1 > compFrac,
           arraySlice2d< real64, compflow::USD_COMP_DC - 1 > dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );

  template< localIndex NC >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView2d< real64, compflow::USD_COMP > const & compFrac,
          arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens );
};

/******************************** PhaseVolumeFractionKernel ********************************/

/**
 * @brief Functions to compute phase volume fractions (saturations) and derivatives
 */
struct PhaseVolumeFractionKernel
{
  template< localIndex NC, localIndex NP >
  GEOSX_HOST_DEVICE
  static void
  compute( arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & compDens,
           arraySlice1d< real64 const, compflow::USD_COMP - 1 > const & dCompDens,
           arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseFrac,
           arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseFrac_dPres,
           arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseFrac_dComp,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( localIndex const size,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );

  template< localIndex NC, localIndex NP >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp,
          arrayView2d< real64, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dComp );
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
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
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
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
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
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
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
          arrayView2d< real64 const, compflow::USD_COMP > const & compFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
      {
        fluidWrapper.update( k, q, pres[k], temp, compFrac[k] );
      }
    } );
  }
};

/******************************** RelativePermeabilityUpdateKernel ********************************/

struct RelativePermeabilityUpdateKernel
{
  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( localIndex const size,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename RELPERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          RELPERM_WRAPPER const & relPermWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < relPermWrapper.numGauss(); ++q )
      {
        relPermWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }
};

/******************************** CapillaryPressureUpdateKernel ********************************/

struct CapillaryPressureUpdateKernel
{
  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( localIndex const size,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
      }
    } );
  }

  template< typename POLICY, typename CAPPRES_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          CAPPRES_WRAPPER const & capPresWrapper,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < capPresWrapper.numGauss(); ++q )
      {
        capPresWrapper.update( k, q, phaseVolFrac[k] );
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
  static void
    compute( localIndex const numPhases,
             real64 const & volume,
             real64 const & porosityOld,
             real64 const & porosityRef,
             real64 const & pvMult,
             real64 const & dPvMult_dPres,
             arraySlice2d< real64 const, compflow::USD_COMP_DC - 1 > const & dCompFrac_dCompDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFracOld,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
             arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
             arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseDensOld,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & phaseDens,
             arraySlice1d< real64 const, multifluid::USD_PHASE - 2 > const & dPhaseDens_dPres,
             arraySlice2d< real64 const, multifluid::USD_PHASE_DC - 2 > const & dPhaseDens_dComp,
             arraySlice2d< real64 const, compflow::USD_PHASE_COMP-1 > const & phaseCompFracOld,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > const & phaseCompFrac,
             arraySlice2d< real64 const, multifluid::USD_PHASE_COMP-2 > const & dPhaseCompFrac_dPres,
             arraySlice3d< real64 const, multifluid::USD_PHASE_COMP_DC-2 > const & dPhaseCompFrac_dComp,
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
          arrayView1d< real64 >       const & porosity,
          arrayView1d< real64 const > const & porosityOld,
          arrayView1d< real64 const > const & porosityRef,
          arrayView2d< real64 const > const & pvMult,
          arrayView2d< real64 const > const & dPvMult_dPres,
          arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFracOld,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseDensOld,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens,
          arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres,
          arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp,
          arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & phaseCompFracOld,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFrac,
          arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dPres,
          arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFrac_dComp,
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
  static void
  compute( real64 const & volume,
           real64 const & porosityRef,
           real64 const & pvMult,
           real64 const & dPvMult_dPres,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & phaseVolFrac,
           arraySlice1d< real64 const, compflow::USD_PHASE - 1 > const & dPhaseVolFrac_dPres,
           arraySlice2d< real64 const, compflow::USD_PHASE_DC - 1 > const & dPhaseVolFrac_dCompDens,
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
          arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac,
          arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres,
          arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY >
  static void launch( arrayView1d< real64 const > const & localResidual,
                      globalIndex const rankOffset,
                      localIndex const numComponents,
                      arrayView1d< globalIndex const > const & dofNumber,
                      arrayView1d< integer const > const & ghostRank,
                      arrayView1d< real64 const > const & refPoro,
                      arrayView1d< real64 const > const & volume,
                      arrayView1d< real64 const > const & totalDensOld,
                      real64 & localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        real64 const normalizer = totalDensOld[ei] * refPoro[ei] * volume[ei];

        for( localIndex idof = 0; idof < numComponents + 1; ++idof )
        {
          real64 const val = localResidual[localRow + idof] / normalizer;
          localSum += val * val;
        }
      }
    } );
    localResidualNorm += localSum.get();
  }

};



/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY >
  static localIndex
  launch( arrayView1d< real64 const > const & localSolution,
          globalIndex const rankOffset,
          localIndex const numComponents,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView2d< real64 const, compflow::USD_COMP > const & compDens,
          arrayView2d< real64 const, compflow::USD_COMP > const & dCompDens,
          integer const allowCompDensChopping,
          real64 const scalingFactor )
  {
    real64 constexpr eps = minDensForDivision;

    RAJA::ReduceMin< REDUCE_POLICY, integer > check( 1 );

    forAll< POLICY >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        {
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
          check.min( newPres >= 0.0 );
        }

        // if component density chopping is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in ApplySystemSolution)
        if( !allowCompDensChopping )
        {
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newDens >= 0.0 );
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( localIndex ic = 0; ic < numComponents; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          check.min( totalDens >= eps );
        }
      }
    } );
    return check.get();
  }

};


/******************************** Kernel launch machinery ********************************/

namespace internal
{

template< typename T, typename LAMBDA >
void KernelLaunchSelectorCompSwitch( T value, LAMBDA && lambda )
{
  static_assert( std::is_integral< T >::value, "KernelLaunchSelectorCompSwitch: type should be integral" );

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
  internal::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
  {
    KERNELWRAPPER::template launch< NC() >( std::forward< ARGS >( args )... );
  } );
}

template< typename KERNELWRAPPER, typename ... ARGS >
void KernelLaunchSelector2( localIndex numComp, localIndex numPhase, ARGS && ... args )
{
  internal::KernelLaunchSelectorCompSwitch( numComp, [&] ( auto NC )
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

} // namespace CompositionalMultiphaseBaseKernels

} // namespace geosx


#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEBASEKERNELS_HPP
