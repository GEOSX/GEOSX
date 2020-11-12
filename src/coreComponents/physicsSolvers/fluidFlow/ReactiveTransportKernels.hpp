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
 * @file ReactiveTransportKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVETRANSPORTKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVETRANSPORTKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "physicsSolvers/fluidFlow/ReactiveTransport.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace ReactiveTransportKernels
{

/******************************** AccumulationKernel ********************************/

struct AccumulationKernel
{
  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const NC,
           arraySlice1d< real64 const > const & componentConc,
           arraySlice1d< real64 const > const & dComponentConc,
           real64 const volume,
           arraySlice1d< real64 > const & localAccum,
           arraySlice2d< real64 > const & localAccumJacobian );

  static void
  Launch( localIndex const size,
          localIndex const NC,
          localIndex const NDOF,
          globalIndex const rankOffset,
          arrayView1d< globalIndex const > const & dofNumber,
          arrayView1d< integer const > const & elemGhostRank,
          arrayView2d< real64 const > const & componentConc,
          arrayView2d< real64 const > const & dComponentConc,
          arrayView1d< real64 const > const & porosity,
          arrayView1d< real64 const > const & volume,
          real64 const dt,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );
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
  using ElementViewConst = typename ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using ElementView = typename ElementRegionManager::ElementView< VIEWTYPE >;

  /**
   * @brief launches the kernel to assemble the flux contributions to the linear system.
   * @tparam STENCIL_TYPE The type of the stencil that is being used.
   */
  template< typename STENCIL_TYPE >
  static void
  Launch( STENCIL_TYPE const & stencil,
          localIndex const NC,
          real64 const dt,
          globalIndex const rankOffset,
          ElementViewConst< arrayView1d< globalIndex const > > const & dofNumber,
          ElementViewConst< arrayView1d< integer const > > const & ghostRank,
          ElementViewConst< arrayView1d< real64 const > > const & pres,
          ElementViewConst< arrayView1d< real64 const > > const & dPres,
          ElementViewConst< arrayView2d< real64 const > > const & componentConc,
          ElementViewConst< arrayView2d< real64 const > > const & dComponentConc,
          real64 const viscosity,
          CRSMatrixView< real64, globalIndex const > const & localMatrix,
          arrayView1d< real64 > const & localRhs );


  GEOSX_HOST_DEVICE
  static void
  Compute( localIndex const stencilSize,
           localIndex const NC,
           arraySlice1d< localIndex const > const & seri,
           arraySlice1d< localIndex const > const & sesri,
           arraySlice1d< localIndex const > const & sei,
           arraySlice1d< real64 const > const & stencilWeights,
           ElementViewConst< arrayView1d< real64 const > > const & pres,
           ElementViewConst< arrayView1d< real64 const > > const & dPres,
           ElementViewConst< arrayView2d< real64 const > > const & componentConc,
           ElementViewConst< arrayView2d< real64 const > > const & dComponentConc,
           real64 const viscosity,
           arraySlice1d< real64 > const & flux,
           arraySlice2d< real64 > const & fluxJacobian );
};

} // namespace ReactiveTransportKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_REACTIVETRANSPORTKERNELS_HPP_
