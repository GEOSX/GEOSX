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
 * @file SolidMechanicsSmallStrainExplicitNewmarkKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"


namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

#if defined(GEOSX_USE_CUDA)
  #define CALCFEMSHAPE
#endif
// If UPDATE_STRESS is undef, uses total displacement and stress is not updated at all.
//  #define UPDATE_STRESS 1 // uses total displacement to and adds material
// stress state to integral for nodalforces.
#define UPDATE_STRESS 2 // uses velocity*dt and updates material stress state.


/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### Explicit Small Strain Description
 * Implements the KernelBase interface functions required for explicit time
 * integration of the equations of motion using the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the interface for KernelBase is used, but
 * ExplicitSmallStrain only conforms to the interface set by KernelBase, and
 * does not inherit from KernelBase.
 * The number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3`.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          int NUM_NODES_PER_ELEM,
          int UNUSED >
class ExplicitSmallStrain : public finiteElement::KernelBase< SUBREGION_TYPE,
                                                              CONSTITUTIVE_TYPE,
                                                              NUM_NODES_PER_ELEM,
                                                              NUM_NODES_PER_ELEM,
                                                              3,
                                                              3 >
{
public:

  /// Alias for the base class;
  using Base = finiteElement::KernelBase< SUBREGION_TYPE,
                                          CONSTITUTIVE_TYPE,
                                          NUM_NODES_PER_ELEM,
                                          NUM_NODES_PER_ELEM,
                                          3,
                                          3 >;

  static constexpr int numTestDofPerSP = 3;
  static constexpr int numTrialDofPerSP = 3;
  static constexpr int numNodesPerElem = NUM_NODES_PER_ELEM;

  using Base::m_elemsToNodes;
  using Base::m_elemGhostRank;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;




//*****************************************************************************

  ExplicitSmallStrain( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       SUBREGION_TYPE const & elementSubRegion,
                       FiniteElementBase const * const finiteElementSpace,
                       CONSTITUTIVE_TYPE * const inputConstitutiveType,
                       real64 const dt,
                       string const & elementListName ):
    Base( elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType ),
#if !defined(CALCFEMSHAPE)
    m_dNdX( elementSubRegion.template getReference< array3d< R1Tensor > >( dataRepository::keys::dNdX )),
    m_detJ( elementSubRegion.template getReference< array2d< real64 > >( dataRepository::keys::detJ ) ),
#endif
    m_X( nodeManager.referencePosition()),
    m_u( nodeManager.totalDisplacement()),
    m_vel( nodeManager.velocity()),
    m_acc( nodeManager.acceleration() ),
    m_dt( dt ),
    m_elementList( elementSubRegion.template getReference< SortedArray< localIndex > >( elementListName ).toViewConst() )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
  }

  //*****************************************************************************
  struct StackVariables
  {
public:
    GEOSX_HOST_DEVICE
    StackVariables():
      fLocal{ { 0.0} },
      varLocal{ {0.0} }
  #if defined(CALCFEMSHAPE)
      ,
      xLocal(),
      dNdX(),
      detJ()
  #endif
    {}

    real64 fLocal[ numNodesPerElem ][ numTrialDofPerSP ];
    real64 varLocal[ numNodesPerElem ][ numTestDofPerSP ];
  #if defined(CALCFEMSHAPE)
// This needs to be returned to service when the FEM kernels are expanded properly
//    real64 xLocal[ numNodesPerElem ][ numTestDofPerSP ];
//    real64 dNdX[ numNodesPerElem ][ numTestDofPerSP ];
    real64 xLocal[ 8 ][ numTestDofPerSP ];
    real64 dNdX[ 8 ][ numTestDofPerSP ];
    real64 detJ;
  #endif
  };



  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i=0; i<numTrialDofPerSP; ++i )
      {
#if defined(CALCFEMSHAPE)
        stack.xLocal[ a ][ i ] = m_X[ nodeIndex ][ i ];
#endif

#if UPDATE_STRESS==2
        stack.varLocal[ a ][ i ] = m_vel[ nodeIndex ][ i ] * m_dt;
#else
        stack.varLocal[ a ][ i ] = m_u[ nodeIndex ][ i ];
#endif
      }
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointStateUpdate( localIndex const k,
                                   localIndex const q,
                                   StackVariables & stack ) const
  {

#if defined(CALCFEMSHAPE)
    real64 dNdX[ 8 ][ 3 ];
    real64 const detJ = FiniteElementShapeKernel::shapeFunctionDerivatives( q, stack.xLocal, dNdX );
#define DNDX dNdX
  #define DETJ detJ
#else //defined(CALCFEMSHAPE)
  #define DNDX m_dNdX[k][q]
  #define DETJ m_detJ( k, q )
#endif //defined(CALCFEMSHAPE)

    real64 stressLocal[ 6 ] = {0};
    real64 strain[6] = {0};
    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      strain[0] = strain[0] + DNDX[ a ][0] * stack.varLocal[ a ][0];
      strain[1] = strain[1] + DNDX[ a ][1] * stack.varLocal[ a ][1];
      strain[2] = strain[2] + DNDX[ a ][2] * stack.varLocal[ a ][2];
      strain[3] = strain[3] + DNDX[ a ][2] * stack.varLocal[ a ][1] + DNDX[ a ][1] * stack.varLocal[ a ][2];
      strain[4] = strain[4] + DNDX[ a ][2] * stack.varLocal[ a ][0] + DNDX[ a ][0] * stack.varLocal[ a ][2];
      strain[5] = strain[5] + DNDX[ a ][1] * stack.varLocal[ a ][0] + DNDX[ a ][0] * stack.varLocal[ a ][1];
    }

#if UPDATE_STRESS == 2
    m_constitutiveUpdate.SmallStrain( k, q, strain );
#else
    m_constitutiveUpdate.SmallStrainNoState( k, strain, stressLocal );
#endif

    for( localIndex c = 0; c < 6; ++c )
    {
#if UPDATE_STRESS == 2
      stressLocal[ c ] =  m_constitutiveUpdate.m_stress( k, q, c ) * (-DETJ);
#elif UPDATE_STRESS == 1
      stressLocal[ c ] = ( stressLocal[ c ] + m_constitutiveUpdate.m_stress( k, q, c ) ) *(-DETJ);
#else
      stressLocal[ c ] *= -DETJ;
#endif
    }

    for( localIndex a=0; a< NUM_NODES_PER_ELEM; ++a )
    {
      stack.fLocal[ a ][ 0 ] = stack.fLocal[ a ][ 0 ] + stressLocal[ 0 ] * DNDX[ a ][ 0 ] + stressLocal[ 5 ] * DNDX[ a ][ 1 ] + stressLocal[ 4 ] * DNDX[ a ][ 2 ];
      stack.fLocal[ a ][ 1 ] = stack.fLocal[ a ][ 1 ] + stressLocal[ 5 ] * DNDX[ a ][ 0 ] + stressLocal[ 1 ] * DNDX[ a ][ 1 ] + stressLocal[ 3 ] * DNDX[ a ][ 2 ];
      stack.fLocal[ a ][ 2 ] = stack.fLocal[ a ][ 2 ] + stressLocal[ 4 ] * DNDX[ a ][ 0 ] + stressLocal[ 3 ] * DNDX[ a ][ 1 ] + stressLocal[ 2 ] * DNDX[ a ][ 2 ];
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointJacobianContribution( localIndex const,
                                            localIndex const,
                                            StackVariables & ) const
  {}

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointResidualContribution( localIndex const GEOSX_UNUSED_PARAM( k ),
                                            localIndex const GEOSX_UNUSED_PARAM( q ),
                                            StackVariables & GEOSX_UNUSED_PARAM( stack ) ) const
  {}

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables const & stack ) const
  {
    for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int b = 0; b < numTestDofPerSP; ++b )
      {
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_acc( nodeIndex, b ), stack.fLocal[ a ][ b ] );
      }
    }
    return 0;
  }


  template< typename POLICY,
            int NUM_QUADRATURE_POINTS,
            typename KERNEL_TYPE >
  static real64
  Launch( localIndex const, //numElems,
          KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    localIndex const numElems = kernelComponent.m_elementList.size();
    forAll< POLICY >( numElems,
                      [=] GEOSX_DEVICE ( localIndex const index )
    {
      localIndex const k = kernelComponent.m_elementList[ index ];

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<NUM_QUADRATURE_POINTS; ++q )
      {
        kernelComponent.quadraturePointStateUpdate( k, q, stack );

        kernelComponent.quadraturePointJacobianContribution( k, q, stack );

        kernelComponent.quadraturePointResidualContribution( k, q, stack );
      }
      kernelComponent.complete( k, stack );
    } );
    return 0;
  }


protected:
  #if !defined(CALCFEMSHAPE)
  arrayView3d< R1Tensor const > const m_dNdX;
  arrayView2d< real64 const > const m_detJ;
  #endif
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_u;
  arrayView2d< real64 const, nodes::VELOCITY_USD > const m_vel;
  arrayView2d< real64, nodes::ACCELERATION_USD > const m_acc;
  real64 const m_dt;
  SortedArrayView< localIndex const > const m_elementList;


};
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_