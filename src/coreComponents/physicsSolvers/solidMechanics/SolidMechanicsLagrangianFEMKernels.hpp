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
 * @file SolidMechanicsLagrangianFEMKernels.hpp
 */

#pragma once

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveBase.hpp"
#include "finiteElement/Kinematics.h"
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "TimeIntegrationOption.hpp"

namespace geosx
{
namespace SolidMechanicsLagrangianFEMKernels
{
inline void
velocityUpdate(
  arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
  arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
  real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = acceleration.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE( localIndex const i ) {
    LvArray::tensorOps::scaledAdd< 3 >( velocity[i], acceleration[i], dt );
    LvArray::tensorOps::fill< 3 >( acceleration[i], 0 );
  } );
}

inline void
velocityUpdate(
  arrayView2d< real64, nodes::ACCELERATION_USD > const & acceleration,
  arrayView1d< real64 const > const & mass,
  arrayView2d< real64, nodes::VELOCITY_USD > const & velocity,
  real64 const dt,
  SortedArrayView< localIndex const > const & indices )
{
  GEOSX_MARK_FUNCTION;

  forAll< parallelDevicePolicy<> >(
    indices.size(),
    [=] GEOSX_DEVICE( localIndex const i ) {
      localIndex const a = indices[i];
      LvArray::tensorOps::scale< 3 >( acceleration[a], 1.0 / mass[a] );
      LvArray::tensorOps::scaledAdd< 3 >( velocity[a], acceleration[a], dt );
    } );
}

inline void
displacementUpdate(
  arrayView2d< real64 const, nodes::VELOCITY_USD > const & velocity,
  arrayView2d< real64, nodes::INCR_DISPLACEMENT_USD > const & uhat,
  arrayView2d< real64, nodes::TOTAL_DISPLACEMENT_USD > const & u,
  real64 const dt )
{
  GEOSX_MARK_FUNCTION;

  localIndex const N = velocity.size( 0 );
  forAll< parallelDevicePolicy<> >( N, [=] GEOSX_DEVICE( localIndex const i ) {
    LvArray::tensorOps::scaledCopy< 3 >( uhat[i], velocity[i], dt );
    LvArray::tensorOps::add< 3 >( u[i], uhat[i] );
  } );
}

/**
 * @brief Function to select which templated kernel function to call.
 * @tparam KERNELWRAPPER A struct or class that contains the following method
 *  "Launch<NUM_NODES_PER_ELEM, NUM_QUADRATURE_POINTS, CONSTITUTIVE_TYPE>( CONSTITUTIVE_TYPE *, PARAMS... )"
 * @tparam PARAMS Variadic parameter pack to pass arguments to Launch function.
 * @param NUM_NODES_PER_ELEM The number of nodes in an element.
 * @param NUM_QUADRATURE_POINTS The number of quadrature points in an element.
 * @param params Variadic parameter list to hold all parameters that are forwarded to the kernel function.
 * @return Depends on the kernel.
 */
//template< typename KERNELWRAPPER, typename ... PARAMS >
//inline real64
//ElementKernelLaunchSelector( localIndex NUM_NODES_PER_ELEM,
//                             localIndex NUM_QUADRATURE_POINTS,
//                             constitutive::ConstitutiveBase * const constitutiveRelation,
//                             PARAMS && ... params )
//{
//  real64 rval = 0;
//
//  using namespace constitutive;
//
//  ConstitutivePassThru< SolidBase >::Execute( constitutiveRelation,
//                                              [&]( auto * const constitutive )
//  {
//    rval = finiteElementLaunchDispatch< KERNELWRAPPER >( NUM_NODES_PER_ELEM, NUM_QUADRATURE_POINTS, &constitutive, std::forward< PARAMS >(
// params )... );
//  } );
//  return rval;
//}

/**
 * @struct Structure to wrap templated function that implements the explicit time integration kernel.
 */
struct ExplicitKernel
{
#if defined( GEOSX_USE_CUDA )
  #define CALCFEMSHAPE
#endif

  template< int N, int USD >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE static void
  Integrate(
    arraySlice1d< real64 const, USD > const & fieldVar,
#if defined( CALCFEMSHAPE )
    real64 const ( &dNdX )[N][3],
#else
    arraySlice2d< real64 const > const & dNdX,
#endif
    real64 const detJ,
    real64 const detF,
    real64 const ( &fInv )[3][3],
    real64 ( &result )[N][3] )
  {
    GEOSX_ASSERT_EQ( fieldVar.size(), 6 );

    real64 const integrationFactor = -detJ * detF;

    real64 P[3][3];
    LvArray::tensorOps::symAikBjk< 3 >( P, fieldVar, fInv );
    LvArray::tensorOps::scale< 3, 3 >( P, integrationFactor );

    for( int a = 0; a < N; ++a )  // loop through all shape functions in element
    {
      LvArray::tensorOps::plusAijBj< 3, 3 >( result[a], P, dNdX[a] );
    }
  }

  static inline real64
  CalculateSingleNodalForce(
    localIndex const k,
    localIndex const targetNode,
    localIndex const numQuadraturePoints,
    arrayView4d< real64 const > const & dNdX,
    arrayView2d< real64 const > const & detJ,
    arrayView3d< real64 const, solid::STRESS_USD > const & stress,
    R1Tensor & force )
  {
    GEOSX_MARK_FUNCTION;
    localIndex const & a = targetNode;

    //Compute Quadrature
    for( localIndex q = 0; q < numQuadraturePoints; ++q )
    {
      force[0] -= ( stress( k, q, 0 ) * dNdX( k, q, a, 0 ) +
                    stress( k, q, 5 ) * dNdX( k, q, a, 1 ) +
                    stress( k, q, 4 ) * dNdX( k, q, a, 2 ) ) *
        detJ( k, q );
      force[1] -= ( stress( k, q, 5 ) * dNdX( k, q, a, 0 ) +
                    stress( k, q, 1 ) * dNdX( k, q, a, 1 ) +
                    stress( k, q, 3 ) * dNdX( k, q, a, 2 ) ) *
        detJ( k, q );
      force[2] -= ( stress( k, q, 4 ) * dNdX( k, q, a, 0 ) +
                    stress( k, q, 3 ) * dNdX( k, q, a, 1 ) +
                    stress( k, q, 2 ) * dNdX( k, q, a, 2 ) ) *
        detJ( k, q );

    }  //quadrature loop

    return 0;
  }
};

}  // namespace SolidMechanicsLagrangianFEMKernels

}  // namespace geosx
