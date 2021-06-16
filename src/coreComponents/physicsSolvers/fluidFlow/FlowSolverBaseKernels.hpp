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
 * @file FlowSolverBaseKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP

#include "common/DataTypes.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"

namespace geosx
{

namespace FlowSolverBaseKernels
{

/******************************** PorosityKernel ********************************/
struct PorosityKernel
{
  template< typename POLICY, typename SOLID_WRAPPER >
  static void
  launch( localIndex const size,
          SOLID_WRAPPER const & solidWrapper,
          arrayView1d< real64 const > const & pres,
          arrayView1d< real64 const > const & dPres )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < solidWrapper.numGauss(); ++q )
      {
        solidWrapper.updatePorosity( k, q, pres[k] + dPres[k] );
      }
    } );
  }

};

/******************************** PermeabilityKernel ********************************/

template< typename REGIONTYPE >
struct PermeabilityKernel
{};

template<>
struct PermeabilityKernel< CellElementSubRegion >
{
  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( localIndex const size,
          PERM_WRAPPER & permWrapper,
          arrayView2d< real64 const > const & porosity )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.updatePorosity( k, q, porosity[k][q] );
      }
    } );
  }

  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          PERM_WRAPPER & permWrapper,
          arrayView2d< real64 const > const & porosity )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.updatePorosity( k, q, porosity[k][q] );
      }
    } );
  }
};

template<>
struct PermeabilityKernel< SurfaceElementSubRegion >
{
  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( localIndex const size,
          PERM_WRAPPER & permWrapper,
          arrayView1d< real64 const > const & effectiveAperture )
  {
    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.updateAperture( k, q, effectiveAperture[k] );
      }
    } );
  }

  template< typename POLICY, typename PERM_WRAPPER >
  static void
  launch( SortedArrayView< localIndex const > const & targetSet,
          PERM_WRAPPER & permWrapper,
          arrayView1d< real64 const > const & effectiveAperture )
  {
    forAll< POLICY >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const k = targetSet[a];
      for( localIndex q = 0; q < permWrapper.numGauss(); ++q )
      {
        permWrapper.updateAperture( k, q, effectiveAperture[k] );
      }
    } );
  }
};

} // FlowSolverBaseKernels

} // geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_FLOWSOLVERBASEKERNELS_HPP
