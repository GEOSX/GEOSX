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
 * @file SinglePhaseProppantBaseKernels.hpp
 */
#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASEKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASEKERNELS_HPP_

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
namespace SinglePhaseProppantBaseKernels
{
/******************************** FluidUpdateKernel ********************************/

struct FluidUpdateKernel
{
  template< typename FLUID_WRAPPER >
  static void
  Launch( FLUID_WRAPPER const & fluidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres,
          arrayView1d< real64 const > const & proppantConcentration,
          arrayView1d< real64 const > const & dProppantConcentration,
          arrayView2d< real64 const > const & componentConcentration,
          arrayView1d< R1Tensor const > const & cellBasedFlux,
          arrayView1d< integer const > const & isProppantBoundaryElement )
  {
    forAll< parallelDevicePolicy<> >(
      fluidWrapper.numElems(),
      [=] GEOSX_HOST_DEVICE( localIndex const a ) {
        for( localIndex q = 0; q < fluidWrapper.numGauss(); ++q )
        {
          fluidWrapper.Update(
            a,
            q,
            pres[a] + dPres[a],
            proppantConcentration[a] + dProppantConcentration[a],
            componentConcentration[a],
            cellBasedFlux[a].L2_Norm(),
            isProppantBoundaryElement[a] );
        }
      } );
  }
};

}  //namespace SinglePhaseProppantBaseKernels

}  //namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEPROPPANTBASEKERNELS_HPP_
