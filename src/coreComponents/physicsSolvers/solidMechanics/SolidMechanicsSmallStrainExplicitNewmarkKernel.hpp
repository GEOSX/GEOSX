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
 * @file SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_

#include "finiteElement/kernelInterface/KernelBase.hpp"

namespace geosx
{
/// Namespace to contain the solid mechanics kernels.
namespace SolidMechanicsLagrangianFEMKernels
{
#if defined( GEOSX_USE_CUDA )
  /// Macro variable to indicate whether or not to calculate the shape function
  /// derivatives in the kernel instead of using a pre-calculated value.
  #define CALCFEMSHAPE
#endif
/// If UPDATE_STRESS is undef, uses total displacement and stress is not
/// updated at all.
/// If UPDATE_STRESS 1, uses total displacement to and adds material stress
/// state to integral for nodalforces.
/// If UPDATE_STRESS 2 then velocity*dt is used to update material stress state
#define UPDATE_STRESS 2

/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the small strain assumption.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam SUBREGION_TYPE The type of subregion that the kernel will act on.
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
template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
class ExplicitSmallStrain
  : public finiteElement::KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 3, 3 >
{
public:
  /// Alias for the base class;
  using Base =
    finiteElement::KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 3, 3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  using Base::m_constitutiveUpdate;
  using Base::m_elemGhostRank;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;

  //*****************************************************************************
  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param dt The time interval for the step.
   * @param elementListName The name of the entry that holds the list of
   *   elements to be processed during this kernel launch.
   */
  ExplicitSmallStrain( NodeManager & nodeManager,
                       EdgeManager const & edgeManager,
                       FaceManager const & faceManager,
                       SUBREGION_TYPE const & elementSubRegion,
                       FE_TYPE const & finiteElementSpace,
                       CONSTITUTIVE_TYPE * const inputConstitutiveType,
                       real64 const dt,
                       string const & elementListName ) :
    Base( elementSubRegion, finiteElementSpace, inputConstitutiveType ),
#if !defined( CALCFEMSHAPE )
    m_dNdX( elementSubRegion.dNdX() ),
    m_detJ( elementSubRegion.detJ() ),
#endif
    m_X( nodeManager.referencePosition() ),
    m_u( nodeManager.totalDisplacement() ),
    m_vel( nodeManager.velocity() ),
    m_acc( nodeManager.acceleration() ),
    m_dt( dt ),
    m_elementList(
      elementSubRegion
        .template getReference< SortedArray< localIndex > >( elementListName )
        .toViewConst() )
  {
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
  }

  //*****************************************************************************
  /**
   * @copydoc geosx::finiteElement::KernelBase::StackVariables
   *
   * ### ExplicitSmallStrain Description
   * Adds a stack arrays for the nodal force, primary displacement variable, etc.
   */
  struct StackVariables : Base::StackVariables
  {
  public:
    GEOSX_HOST_DEVICE
    StackVariables() :
      // clang-format off
      fLocal{ { 0.0} },
      varLocal{ {0.0} }
    // clang-format on
#if defined( CALCFEMSHAPE )
    , xLocal(), dNdX(), detJ()
#endif
    {}

    /// C-array stack storage for the element local force
    real64 fLocal[numNodesPerElem][numDofPerTrialSupportPoint];

    /// C-array stack storage for element local primary variable values.
    real64 varLocal[numNodesPerElem][numDofPerTestSupportPoint];
#if defined( CALCFEMSHAPE )
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][numDofPerTestSupportPoint];

    /// C-array stack storage for shape function derivatives at a point.
    real64 dNdX[numNodesPerElem][numDofPerTestSupportPoint];

    /// C-array stack storage for the jacobian of the parent space mapping.
    real64 detJ;
#endif
  };
  //***************************************************************************

  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * Copies the primary variable, and position into the local stack array.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  setup( localIndex const k, StackVariables & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int i = 0; i < numDofPerTrialSupportPoint; ++i )
      {
#if defined( CALCFEMSHAPE )
        stack.xLocal[a][i] = m_X[nodeIndex][i];
#endif

#if UPDATE_STRESS == 2
        stack.varLocal[a][i] = m_vel[nodeIndex][i] * m_dt;
#else
        stack.varLocal[a][i] = m_u[nodeIndex][i];
#endif
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointStateUpdate
   *
   * ### ExplicitSmallStrain Description
   * Calculates the shape function derivatives, and the strain tensor. Then
   * calls the constitutive update, and also performs the integration of
   * the stress divergence, rather than using the dedicated component function
   * to allow for some variable reuse.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  quadraturePointStateUpdate( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
#if defined( CALCFEMSHAPE )
    real64 dNdX[numNodesPerElem][3];
    real64 const detJ = FE_TYPE::shapeFunctionDerivatives( q, stack.xLocal, dNdX );

  /// Macro to substitute in the shape function derivatives.
  #define DNDX dNdX

  /// Macro to substitute the determinant of the jacobian transformation to the parent space.
  #define DETJ detJ
#else  //defined(CALCFEMSHAPE)
  /// @cond DOXYGEN_SKIP
  #define DNDX m_dNdX[k][q]
  #define DETJ m_detJ( k, q )
      /// @endcond DOXYGEN_SKIP
#endif  //defined(CALCFEMSHAPE)

    real64 stressLocal[6] = { 0 };
    real64 strain[6] = { 0 };
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      strain[0] = strain[0] + DNDX[a][0] * stack.varLocal[a][0];
      strain[1] = strain[1] + DNDX[a][1] * stack.varLocal[a][1];
      strain[2] = strain[2] + DNDX[a][2] * stack.varLocal[a][2];
      strain[3] = strain[3] + DNDX[a][2] * stack.varLocal[a][1] +
        DNDX[a][1] * stack.varLocal[a][2];
      strain[4] = strain[4] + DNDX[a][2] * stack.varLocal[a][0] +
        DNDX[a][0] * stack.varLocal[a][2];
      strain[5] = strain[5] + DNDX[a][1] * stack.varLocal[a][0] +
        DNDX[a][0] * stack.varLocal[a][1];
    }

#if UPDATE_STRESS == 2
    m_constitutiveUpdate.SmallStrain( k, q, strain );
#else
    m_constitutiveUpdate.SmallStrainNoState( k, strain, stressLocal );
#endif

    for( localIndex c = 0; c < 6; ++c )
    {
#if UPDATE_STRESS == 2
      stressLocal[c] = m_constitutiveUpdate.m_stress( k, q, c ) * ( -DETJ );
#elif UPDATE_STRESS == 1
      stressLocal[c] =
        ( stressLocal[c] + m_constitutiveUpdate.m_stress( k, q, c ) ) * ( -DETJ );
#else
      stressLocal[c] *= -DETJ;
#endif
    }

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      stack.fLocal[a][0] = stack.fLocal[a][0] + stressLocal[0] * DNDX[a][0] +
        stressLocal[5] * DNDX[a][1] + stressLocal[4] * DNDX[a][2];
      stack.fLocal[a][1] = stack.fLocal[a][1] + stressLocal[5] * DNDX[a][0] +
        stressLocal[1] * DNDX[a][1] + stressLocal[3] * DNDX[a][2];
      stack.fLocal[a][2] = stack.fLocal[a][2] + stressLocal[4] * DNDX[a][0] +
        stressLocal[3] * DNDX[a][1] + stressLocal[2] * DNDX[a][2];
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * ### ExplicitSmallStrain Description
   * Performs the distribution of the nodal force out to the rank local arrays.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64
  complete( localIndex const k, StackVariables const & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = m_elemsToNodes( k, a );
      for( int b = 0; b < numDofPerTestSupportPoint; ++b )
      {
        RAJA::atomicAdd< parallelDeviceAtomic >( &m_acc( nodeIndex, b ),
                                                 stack.fLocal[a][b] );
      }
    }
    return 0;
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::kernelLaunch
   *
   * ### ExplicitSmallStrain Description
   * Copy of the KernelBase::kernelLaunch function without the exclusion of ghost
   * elements.
   */
  template< typename POLICY, typename KERNEL_TYPE >
  static real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    GEOSX_UNUSED_VAR( numElems );

    localIndex const numProcElems = kernelComponent.m_elementList.size();
    forAll< POLICY >( numProcElems, [=] GEOSX_DEVICE( localIndex const index ) {
      localIndex const k = kernelComponent.m_elementList[index];

      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q = 0; q < KERNEL_TYPE::numQuadraturePointsPerElem; ++q )
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
#if !defined( CALCFEMSHAPE )
  /// The shape function derivative for each quadrature point.
  arrayView4d< real64 const > const m_dNdX;
  /// The parent->physical jacobian determinant for each quadrature point.
  arrayView2d< real64 const > const m_detJ;
#endif
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The array containing the nodal displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_u;

  /// The array containing the nodal velocity array.
  arrayView2d< real64 const, nodes::VELOCITY_USD > const m_vel;

  /// The array containing the nodal acceleration array, which is used to store
  /// the force.
  arrayView2d< real64, nodes::ACCELERATION_USD > const m_acc;

  /// The time increment for this time integration step.
  real64 const m_dt;

  /// The list of elements to process for the kernel launch.
  SortedArrayView< localIndex const > const m_elementList;
};
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

}  // namespace SolidMechanicsLagrangianFEMKernels

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINEXPLICITNEWMARK_HPP_
