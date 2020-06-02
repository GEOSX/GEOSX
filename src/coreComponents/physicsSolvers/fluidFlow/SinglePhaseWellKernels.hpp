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
 * @file SinglePhaseWellKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEWELLKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEWELLKERNELS_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseWell.hpp"
#include "wells/WellControls.hpp"

namespace geosx
{

namespace SinglePhaseWellKernels
{


/******************************** ControlEquationHelper ********************************/

struct ControlEquationHelper
{

  static void
  ComputeJacobianEntry( globalIndex const rankOffset,
                        WellControls const & wellControls,
                        globalIndex const wellElemDofNumber,
                        real64 const & wellElemPressure,
                        real64 const & dWellElemPressure,
                        real64 const & connRate,
                        real64 const & dConnRate,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs )
  {
    globalIndex eqnRowIndex = 0;
    globalIndex dofColIndex = 0;
    real64 controlEqn = 0;
    real64 dControlEqn_dX = 0;

    // get well control and type
    WellControls::Control const currentControl = wellControls.GetControl();

    // BHP control
    if( currentControl == WellControls::Control::BHP )
    {
      // get pressures and compute normalizer
      real64 const currentBHP = wellElemPressure + dWellElemPressure;
      real64 const targetBHP  = wellControls.GetTargetBHP();
      real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                                ? 1.0 / targetBHP
                                : 1.0;

      // control equation is a normalized difference between current pressure and target pressure
      controlEqn = ( currentBHP - targetBHP ) * normalizer;
      dControlEqn_dX = normalizer;
      dofColIndex = wellElemDofNumber + SinglePhaseWell::ColOffset::DPRES;
      eqnRowIndex = wellElemDofNumber + SinglePhaseWell::RowOffset::CONTROL - rankOffset;
    }
    // rate control
    else
    {

      // get rates and compute normalizer
      real64 const currentConnRate = connRate + dConnRate;
      real64 const & targetConnRate = wellControls.GetTargetRate();
      real64 const normalizer = fabs( targetConnRate ) > std::numeric_limits< real64 >::min()
                                ? 1.0 / ( 1e-2 * fabs( targetConnRate ) ) // hard-coded value comes from AD-GPRS
                                : 1.0;

      // for a producer, the actual (target) rate is negative

      // control equation is a normalized difference between current rate and target rate
      controlEqn = ( currentConnRate - targetConnRate ) * normalizer;
      dControlEqn_dX = normalizer;
      dofColIndex = wellElemDofNumber + SinglePhaseWell::ColOffset::DRATE;
      eqnRowIndex = wellElemDofNumber + SinglePhaseWell::RowOffset::CONTROL - rankOffset;
    }

    localMatrix.addToRow< serialAtomic >( eqnRowIndex,
                                          &dofColIndex,
                                          &dControlEqn_dX,
                                          1 );
    localRhs[eqnRowIndex] += controlEqn;
  }
};


/******************************** FluxKernel ********************************/

struct FluxKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & connRate,
          arrayView1d< real64 const > const & dConnRate,
          real64 const & dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    // loop over the well elements to compute the fluxes between elements
    forAll< POLICY >( size, [=]( localIndex const iwelem )
    {

      // 1) Compute the flux and its derivatives

      /*  currentConnRate < 0 flow from iwelem to iwelemNext
       *  currentConnRate > 0 flow from iwelemNext to iwelem
       *  With this convention, currentConnRate < 0 at the last connection for a producer
       *                        currentConnRate > 0 at the last connection for a injector
       */

      // get next well element index
      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      // there is nothing to upwind for single-phase flow
      real64 const currentConnRate = connRate[iwelem] + dConnRate[iwelem];
      real64 const flux = dt * currentConnRate;
      real64 const dFlux_dRate = dt;

      // 2) Assemble the flux into residual and Jacobian
      if( iwelemNext < 0 )
      {
        // flux terms
        real64 const oneSidedLocalFlux = -flux;
        real64 const oneSidedLocalFluxJacobian_dRate = -dFlux_dRate;

        // jacobian indices
        globalIndex const offset = wellElemDofNumber[iwelem];
        globalIndex const oneSidedEqnRowIndex = offset + SinglePhaseWell::RowOffset::MASSBAL - rankOffset;
        globalIndex const oneSidedDofColIndex_dRate = offset + SinglePhaseWell::ColOffset::DRATE;

        if( oneSidedEqnRowIndex >= 0 && oneSidedEqnRowIndex < localMatrix.numRows() )
        {
          localMatrix.addToRow< serialAtomic >( oneSidedEqnRowIndex,
                                                &oneSidedDofColIndex_dRate,
                                                &oneSidedLocalFluxJacobian_dRate,
                                                1 );
          atomicAdd( serialAtomic{}, &localRhs[oneSidedEqnRowIndex], oneSidedLocalFlux );
        }
      }
      else
      {
        // local working variables and arrays
        stackArray1d< globalIndex, 2 > eqnRowIndices( 2 );

        stackArray1d< real64, 2 > localFlux( 2 );
        stackArray1d< real64, 2 > localFluxJacobian_dRate( 2 );

        // flux terms
        localFlux[SinglePhaseWell::ElemTag::NEXT] = flux;
        localFlux[SinglePhaseWell::ElemTag::CURRENT] = -flux;

        localFluxJacobian_dRate[SinglePhaseWell::ElemTag::NEXT] = dFlux_dRate;
        localFluxJacobian_dRate[SinglePhaseWell::ElemTag::CURRENT] = -dFlux_dRate;

        // indices
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];
        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];
        eqnRowIndices[SinglePhaseWell::ElemTag::CURRENT] = offsetCurrent
                                                           + SinglePhaseWell::RowOffset::MASSBAL - rankOffset;
        eqnRowIndices[SinglePhaseWell::ElemTag::NEXT] = offsetNext
                                                        + SinglePhaseWell::RowOffset::MASSBAL - rankOffset;
        globalIndex const dofColIndex_dRate = offsetCurrent + SinglePhaseWell::ColOffset::DRATE;

        for( localIndex i = 0; i < localFlux.size(); ++i )
        {
          if( eqnRowIndices[i] >= 0 && eqnRowIndices[i] < localMatrix.numRows() )
          {
            localMatrix.addToRow< serialAtomic >( eqnRowIndices[i],
                                                  &dofColIndex_dRate,
                                                  &localFluxJacobian_dRate[i],
                                                  1 );
            atomicAdd( serialAtomic{}, &localRhs[eqnRowIndices[i]], localFlux[i] );
          }
        }
      }
    } );
  }
};


/******************************** PressureRelationKernel ********************************/

struct PressureRelationKernel
{

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & wellElemDofNumber,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< localIndex const > const & nextWellElemIndex,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          real64 const & targetBHP,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs )
  {
    // loop over the well elements to compute the pressure relations between well elements
    forAll< POLICY >( size, [=]( localIndex const iwelem )
    {

      localIndex const iwelemNext = nextWellElemIndex[iwelem];

      if( iwelemNext >= 0 ) // if iwelemNext < 0, form control equation, not momentum
      {

        // local working variables and arrays
        stackArray1d< globalIndex, 2 > dofColIndices( 2 );
        stackArray1d< real64, 2 > localPresRelJacobian( 2 );

        dofColIndices = -1;
        localPresRelJacobian = 0;

        // compute avg density
        real64 const avgDensity = 0.5 * ( wellElemDensity[iwelem][0] + wellElemDensity[iwelemNext][0] );
        real64 const dAvgDensity_dPresNext = 0.5 * dWellElemDensity_dPres[iwelemNext][0];
        real64 const dAvgDensity_dPresCurrent = 0.5 * dWellElemDensity_dPres[iwelem][0];

        // compute depth diff times acceleration
        real64 const gravD = wellElemGravCoef[iwelemNext] - wellElemGravCoef[iwelem];

        // compute the current pressure in the two well elements
        real64 const pressureCurrent = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
        real64 const pressureNext = wellElemPressure[iwelemNext] + dWellElemPressure[iwelemNext];

        // compute a coefficient to normalize the momentum equation
        real64 const normalizer = targetBHP > std::numeric_limits< real64 >::epsilon()
                                  ? 1.0 / targetBHP
                                  : 1.0;

        // compute momentum flux and derivatives
        real64 const localPresRel = ( pressureNext - pressureCurrent - avgDensity * gravD ) * normalizer;
        localPresRelJacobian[SinglePhaseWell::ElemTag::NEXT] = ( 1 - dAvgDensity_dPresNext * gravD ) * normalizer;
        localPresRelJacobian[SinglePhaseWell::ElemTag::CURRENT] = ( -1 - dAvgDensity_dPresCurrent * gravD ) * normalizer;

        // TODO: add friction and acceleration terms

        // jacobian indices
        globalIndex const offsetNext = wellElemDofNumber[iwelemNext];
        globalIndex const offsetCurrent = wellElemDofNumber[iwelem];
        globalIndex const eqnRowIndex = offsetCurrent + SinglePhaseWell::RowOffset::CONTROL - rankOffset;
        dofColIndices[SinglePhaseWell::ElemTag::NEXT] = offsetNext + SinglePhaseWell::ColOffset::DPRES;
        dofColIndices[SinglePhaseWell::ElemTag::CURRENT] = offsetCurrent + SinglePhaseWell::ColOffset::DPRES;

        if( eqnRowIndex >= 0 && eqnRowIndex < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< serialAtomic >( eqnRowIndex,
                                                                    dofColIndices.data(),
                                                                    localPresRelJacobian.data(),
                                                                    2 );
          atomicAdd( serialAtomic{}, &localRhs[eqnRowIndex], localPresRel );
        }
      }
    } );
  }

};


/******************************** PerforationKernel ********************************/

struct PerforationKernel
{

  /**
   * @brief The type for element-based non-constitutive data parameters.
   * Consists entirely of ArrayView's.
   *
   * Can be converted from ElementRegionManager::ElementViewAccessor
   * by calling .toView() or .toViewConst() on an accessor instance
   */
  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementViewAccessor< VIEWTYPE >::ViewTypeConst;

  template< typename POLICY >
  static void
  Launch( localIndex const size,
          ElementView< arrayView1d< real64 const > > const & resPressure,
          ElementView< arrayView1d< real64 const > > const & dResPressure,
          ElementView< arrayView2d< real64 const > > const & resDensity,
          ElementView< arrayView2d< real64 const > > const & dResDensity_dPres,
          ElementView< arrayView2d< real64 const > > const & resViscosity,
          ElementView< arrayView2d< real64 const > > const & dResViscosity_dPres,
          arrayView1d< real64 const > const & wellElemGravCoef,
          arrayView1d< real64 const > const & wellElemPressure,
          arrayView1d< real64 const > const & dWellElemPressure,
          arrayView2d< real64 const > const & wellElemDensity,
          arrayView2d< real64 const > const & dWellElemDensity_dPres,
          arrayView2d< real64 const > const & wellElemViscosity,
          arrayView2d< real64 const > const & dWellElemViscosity_dPres,
          arrayView1d< real64 const > const & perfGravCoef,
          arrayView1d< localIndex const > const & perfWellElemIndex,
          arrayView1d< real64 const > const & perfTransmissibility,
          arrayView1d< localIndex const > const & resElementRegion,
          arrayView1d< localIndex const > const & resElementSubRegion,
          arrayView1d< localIndex const > const & resElementIndex,
          arrayView1d< real64 > const & perfRate,
          arrayView2d< real64 > const & dPerfRate_dPres )
  {
    forAll< POLICY >( size, [=]( localIndex const iperf )
    {
      // local working variables and arrays
      stackArray1d< real64, 2 > pressure( 2 );
      stackArray1d< real64, 2 > dPressure_dP( 2 );
      stackArray1d< localIndex, 2 > multiplier( 2 );
      pressure = 0;
      dPressure_dP = 0;
      multiplier = 0;

      // 1) Reservoir side

      // get the reservoir (sub)region and element indices
      localIndex const er  = resElementRegion[iperf];
      localIndex const esr = resElementSubRegion[iperf];
      localIndex const ei  = resElementIndex[iperf];

      // get reservoir variables
      pressure[SinglePhaseWell::SubRegionTag::RES] = resPressure[er][esr][ei] + dResPressure[er][esr][ei];
      dPressure_dP[SinglePhaseWell::SubRegionTag::RES] = 1;

      // TODO: add a buoyancy term for the reservoir side here

      // multiplier for reservoir side in the flux
      multiplier[SinglePhaseWell::SubRegionTag::RES] = 1;

      // 2) Well side

      // get the local index of the well element
      localIndex const iwelem = perfWellElemIndex[iperf];

      // get well variables
      pressure[SinglePhaseWell::SubRegionTag::WELL] = wellElemPressure[iwelem] + dWellElemPressure[iwelem];
      dPressure_dP[SinglePhaseWell::SubRegionTag::WELL] = 1.0;

      real64 const gravD = ( perfGravCoef[iperf] - wellElemGravCoef[iwelem] );
      pressure[SinglePhaseWell::SubRegionTag::WELL] += wellElemDensity[iwelem][0] * gravD;
      dPressure_dP[SinglePhaseWell::SubRegionTag::WELL] += dWellElemDensity_dPres[iwelem][0] * gravD;

      // multiplier for well side in the flux
      multiplier[SinglePhaseWell::SubRegionTag::WELL] = -1;

      // get transmissibility at the interface
      real64 const trans = perfTransmissibility[iperf];

      // compute potential difference
      real64 potDif = 0.0;
      for( localIndex i = 0; i < 2; ++i )
      {
        potDif += multiplier[i] * trans * pressure[i];
        dPerfRate_dPres[iperf][i] = multiplier[i] * trans * dPressure_dP[i];
      }

      // choose upstream cell based on potential difference
      localIndex const k_up = (potDif >= 0)
                            ? SinglePhaseWell::SubRegionTag::RES
                            : SinglePhaseWell::SubRegionTag::WELL;

      // compute upstream density, viscosity, and mobility
      real64 densityUp       = 0.0;
      real64 dDensityUp_dP   = 0.0;
      real64 viscosityUp     = 0.0;
      real64 dViscosityUp_dP = 0.0;

      // upwinding the variables
      if( k_up == SinglePhaseWell::SubRegionTag::RES ) // use reservoir vars
      {
        densityUp     = resDensity[er][esr][ei][0];
        dDensityUp_dP = dResDensity_dPres[er][esr][ei][0];

        viscosityUp     = resViscosity[er][esr][ei][0];
        dViscosityUp_dP = dResViscosity_dPres[er][esr][ei][0];
      }
      else // use well vars
      {
        densityUp = wellElemDensity[iwelem][0];
        dDensityUp_dP = dWellElemDensity_dPres[iwelem][0];

        viscosityUp = wellElemViscosity[iwelem][0];
        dViscosityUp_dP = dWellElemViscosity_dPres[iwelem][0];
      }

      // compute mobility
      real64 const mobilityUp     = densityUp / viscosityUp;
      real64 const dMobilityUp_dP = dDensityUp_dP / viscosityUp
                                    - mobilityUp / viscosityUp * dViscosityUp_dP;

      perfRate[iperf] = mobilityUp * potDif;
      for( localIndex ke = 0; ke < 2; ++ke )
      {
        dPerfRate_dPres[iperf][ke] *= mobilityUp;
      }
      dPerfRate_dPres[iperf][k_up] += dMobilityUp_dP * potDif;
    } );
  }
};


/******************************** ResidualNormKernel ********************************/

struct ResidualNormKernel
{

  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static void Launch( LOCAL_VECTOR const localResidual,
                      globalIndex const rankOffset,
                      arrayView1d< globalIndex const > const & wellElemDofNumber,
                      arrayView1d< integer const > const & wellElemGhostRank,
                      arrayView1d< real64 const > const & wellElemVolume,
                      arrayView2d< real64 const > const & wellElemDensity,
                      real64 * localResidualNorm )
  {
    RAJA::ReduceSum< REDUCE_POLICY, real64 > sumScaled( 0.0 );

    forAll< POLICY >( wellElemDofNumber.size(), [=] ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex idof = 0; idof < 2; ++idof )
        {
          real64 const normalizer = ( idof == SinglePhaseWell::RowOffset::MASSBAL )
                                    ? wellElemDensity[iwelem][0] * wellElemVolume[iwelem]
                                    : 1;
          localIndex const lid = wellElemDofNumber[iwelem] + idof - rankOffset;
          real64 const val = localResidual[lid] / normalizer;
          sumScaled += val * val;
        }
      }
    } );
    *localResidualNorm = *localResidualNorm + sumScaled.get();
  }

};

/******************************** SolutionCheckKernel ********************************/

struct SolutionCheckKernel
{
  template< typename POLICY, typename REDUCE_POLICY, typename LOCAL_VECTOR >
  static localIndex Launch( LOCAL_VECTOR const localSolution,
                            globalIndex const rankOffset,
                            arrayView1d< globalIndex const > const & presDofNumber,
                            arrayView1d< integer const > const & ghostRank,
                            arrayView1d< real64 const > const & pres,
                            arrayView1d< real64 const > const & dPres,
                            real64 const scalingFactor )
  {
    RAJA::ReduceMin< REDUCE_POLICY, localIndex > minVal( 1 );

    forAll< POLICY >( presDofNumber.size(), [=] ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 && presDofNumber[ei] >= 0 )
      {
        localIndex const lid = presDofNumber[ei] - rankOffset;
        real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[lid];

        if( newPres < 0.0 )
        {
          minVal.min( 0 );
        }
      }
    } );
    return minVal.get();
  }
};

} // end namespace SinglePhaseWellKernels

} // end namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEWELLKERNELS_HPP
