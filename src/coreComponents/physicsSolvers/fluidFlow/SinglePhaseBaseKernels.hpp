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
 * @file SinglePhaseBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{
namespace SinglePhaseBaseKernels
{
/******************************** MobilityKernel ********************************/

struct MobilityKernel
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( real64 const & dens,
           real64 const & dDens_dPres,
           real64 const & visc,
           real64 const & dVisc_dPres,
           real64 & mob,
           real64 & dMob_dPres )
  {
    mob = dens / visc;
    dMob_dPres = dDens_dPres / visc - mob / visc * dVisc_dPres;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( real64 const & dens, real64 const & visc, real64 & mob )
  {
    mob = dens / visc;
  }

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & dDens_dPres,
          arrayView2d< real64 const > const & visc,
          arrayView2d< real64 const > const & dVisc_dPres,
          arrayView1d< real64 > const & mob,
          arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const a ) {
      Compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  template< typename POLICY >
  static void
  Launch( SortedArrayView< localIndex const > targetSet,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & dDens_dPres,
          arrayView2d< real64 const > const & visc,
          arrayView2d< real64 const > const & dVisc_dPres,
          arrayView1d< real64 > const & mob,
          arrayView1d< real64 > const & dMob_dPres )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE( localIndex const i ) {
      localIndex const a = targetSet[i];
      Compute( dens[a][0],
               dDens_dPres[a][0],
               visc[a][0],
               dVisc_dPres[a][0],
               mob[a],
               dMob_dPres[a] );
    } );
  }

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & visc,
          arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const a ) {
      Compute( dens[a][0], visc[a][0], mob[a] );
    } );
  }

  template< typename POLICY >
  static void
  Launch( SortedArrayView< localIndex const > targetSet,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & visc,
          arrayView1d< real64 > const & mob )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE( localIndex const i ) {
      localIndex const a = targetSet[i];
      Compute( dens[a][0], visc[a][0], mob[a] );
    } );
  }
};

/******************************** AccumulationKernel ********************************/

template< bool ISPORO >
struct AssembleAccumulationTermsHelper;

template<>
struct AssembleAccumulationTermsHelper< true >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const biotCoefficient,
                  real64 const poroOld,
                  real64 const bulkModulus,
                  real64 const totalMeanStress,
                  real64 const oldTotalMeanStress,
                  real64 const dPres,
                  real64 const GEOSX_UNUSED_PARAM( poroRef ),
                  real64 const GEOSX_UNUSED_PARAM( pvmult ),
                  real64 const GEOSX_UNUSED_PARAM( dPVMult_dPres ) )
  {
    dPoro_dPres = ( biotCoefficient - poroOld ) / bulkModulus;
    poro = poroOld + dPoro_dPres * ( totalMeanStress - oldTotalMeanStress + dPres );
  }
};

template<>
struct AssembleAccumulationTermsHelper< false >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  porosityUpdate( real64 & poro,
                  real64 & dPoro_dPres,
                  real64 const GEOSX_UNUSED_PARAM( biotCoefficient ),
                  real64 const GEOSX_UNUSED_PARAM( poroOld ),
                  real64 const GEOSX_UNUSED_PARAM( bulkModulus ),
                  real64 const GEOSX_UNUSED_PARAM( totalMeanStress ),
                  real64 const GEOSX_UNUSED_PARAM( oldTotalMeanStress ),
                  real64 const GEOSX_UNUSED_PARAM( dPres ),
                  real64 const poroRef,
                  real64 const pvmult,
                  real64 const dPVMult_dPres )
  {
    poro = poroRef * pvmult;
    dPoro_dPres = dPVMult_dPres * poroRef;
  }
};

template< typename REGIONTYPE >
struct AccumulationKernel
{};

template<>
struct AccumulationKernel< CellElementSubRegion >
{
  template< bool COUPLED >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Compute(
    real64 const & dPres,
    real64 const & densNew,
    real64 const & densOld,
    real64 const & dDens_dPres,
    real64 const & volume,
    real64 const & dVol,
    real64 const & poroRef,
    real64 const & poroOld,
    real64 const & pvMult,
    real64 const & dPVMult_dPres,
    real64 const & biotCoefficient,
    real64 const & bulkModulus,
    real64 const & totalMeanStress,
    real64 const & oldTotalMeanStress,
    real64 & poroNew,
    real64 & localAccum,
    real64 & localAccumJacobian )
  {
    real64 const volNew = volume + dVol;

    // TODO porosity update needs to be elsewhere...
    real64 dPoro_dPres;
    AssembleAccumulationTermsHelper< COUPLED >::porosityUpdate( poroNew,
                                                                dPoro_dPres,
                                                                biotCoefficient,
                                                                poroOld,
                                                                bulkModulus,
                                                                totalMeanStress,
                                                                oldTotalMeanStress,
                                                                dPres,
                                                                poroRef,
                                                                pvMult,
                                                                dPVMult_dPres );

    // Residual contribution is mass conservation in the cell
    localAccum = poroNew * densNew * volNew - poroOld * densOld * volume;

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian = ( dPoro_dPres * densNew + dDens_dPres * poroNew ) * volNew;
  }

  template< bool COUPLED, typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & dPres,
          arrayView1d< real64 const > const & densOld,
          arrayView1d< real64 > const & poro,
          arrayView1d< real64 const > const & poroOld,
          arrayView1d< real64 const > const & poroRef,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & dVol,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & dDens_dPres,
          arrayView2d< real64 const > const & pvmult,
          arrayView2d< real64 const > const & dPVMult_dPres,
          arrayView1d< real64 const > const & oldTotalMeanStress,
          arrayView1d< real64 const > const & totalMeanStress,
          arrayView1d< real64 const > const & bulkModulus,
          real64 const biotCoefficient,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
      if( elemGhostRank[ei] < 0 )
      {
        real64 localAccum, localAccumJacobian;

        Compute< COUPLED >( dPres[ei],
                            dens[ei][0],
                            densOld[ei],
                            dDens_dPres[ei][0],
                            volume[ei],
                            dVol[ei],
                            poroRef[ei],
                            poroOld[ei],
                            pvmult[ei][0],
                            dPVMult_dPres[ei][0],
                            biotCoefficient,
                            bulkModulus[ei],
                            totalMeanStress[ei],
                            oldTotalMeanStress[ei],
                            poro[ei],
                            localAccum,
                            localAccumJacobian );

        globalIndex const elemDOF = dofNumber[ei];
        localIndex const localElemDof = elemDOF - rankOffset;

        // add contribution to global residual and jacobian (no need for atomics here)
        localMatrix.addToRow< serialAtomic >( localElemDof,
                                              &elemDOF,
                                              &localAccumJacobian,
                                              1 );
        localRhs[localElemDof] += localAccum;
      }
    } );
  }
};

template<>
struct AccumulationKernel< FaceElementSubRegion >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( real64 const & densNew,
           real64 const & densOld,
           real64 const & dDens_dPres,
           real64 const & volume,
           real64 const & dVol,
           real64 & localAccum,
           real64 & localAccumJacobian )
  {
    real64 const volNew = volume + dVol;

    // Residual contribution is mass conservation in the cell
    localAccum = densNew * volNew - densOld * volume;

    // Derivative of residual wrt to pressure in the cell
    localAccumJacobian = dDens_dPres * volNew;
  }

  template< bool COUPLED, typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView1d< real64 const > const & densOld,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & dVol,
          arrayView2d< real64 const > const & dens,
          arrayView2d< real64 const > const & dDens_dPres,
          arrayView1d< real64 const > const & poroMultiplier,
#if ALLOW_CREATION_MASS
          arrayView1d< real64 const > const & creationMass,
#endif
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
      if( elemGhostRank[ei] < 0 )
      {
        real64 localAccum, localAccumJacobian;

        real64 const effectiveVolume = volume[ei] * poroMultiplier[ei];

        Compute( dens[ei][0],
                 densOld[ei],
                 dDens_dPres[ei][0],
                 effectiveVolume,
                 dVol[ei],
                 localAccum,
                 localAccumJacobian );

#if ALLOW_CREATION_MASS
        if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
        {
          localAccum += creationMass[ei] * 0.25;
        }
#endif

        globalIndex const elemDOF = dofNumber[ei];
        localIndex const localElemDof = elemDOF - rankOffset;

        // add contribution to global residual and jacobian (no need for atomics here)
        localMatrix.addToRow< serialAtomic >( localElemDof,
                                              &elemDOF,
                                              &localAccumJacobian,
                                              1 );
        localRhs[localElemDof] += localAccum;
      }
    } );
  }
};

/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void
  Launch( FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres )
  {
    forAll< parallelDevicePolicy<> >(
      fluidWrapper.numElems(),
      [=] GEOSX_HOST_DEVICE( localIndex const k ) {
        for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
        {
          fluidWrapper.Update( k, q, pres[k] + dPres[k] );
        }
      } );
  }
};

/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void
  Launch( LOCAL_VECTOR const localResidual,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & presDofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & refPoro,
          arrayView1d< real64 const > const & volume,
          arrayView1d< real64 const > const & densOld,
          real64 * localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > localSum( 0.0 );
    RAJA::ReduceSum< REDUCE_POLICY, real64 > normSum( 0.0 );
    RAJA::ReduceSum< REDUCE_POLICY, localIndex > count( 0 );

    forAll< POLICY >( presDofNumber.size(),
                      [=] GEOSX_HOST_DEVICE( localIndex const a ) {
                        if( ghostRank[a] < 0 )
                        {
                          localIndex const lid = presDofNumber[a] - rankOffset;
                          real64 const val = localResidual[lid];
                          localSum += val * val;
                          normSum += refPoro[a] * densOld[a] * volume[a];
                          count += 1;
                        }
                      } );

    localResidualNorm[0] += localSum.get();
    localResidualNorm[1] += normSum.get();
    localResidualNorm[2] += count.get();
  }
};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static localIndex
  Launch( LOCAL_VECTOR const localSolution,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & presDofNumber,
          arrayView1d< integer const > const & ghostRank,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( presDofNumber.size(),
                      [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                        if( ghostRank[ei] < 0 && presDofNumber[ei] >= 0 )
                        {
                          localIndex const lid = presDofNumber[ei] - rankOffset;
                          real64 const newPres = pres[ei] + dPres[ei] +
                            scalingFactor * localSolution[lid];

                          if( newPres < 0.0 )
                          {
                            minVal.min( 0 );
                          }
                        }
                      } );
    return minVal.get();
  }
};

}  // namespace SinglePhaseBaseKernels

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASEKERNELS_HPP
