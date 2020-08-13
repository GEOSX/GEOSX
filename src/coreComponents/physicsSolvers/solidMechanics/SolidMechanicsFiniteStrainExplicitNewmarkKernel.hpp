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
 * @file SolidMechanicsFiniteStrainExplicitNewmarkKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_

#include "SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp"
#include "finiteElement/Kinematics.h"

namespace geosx
{
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
 * @brief Integrate the divergence of a rank-2 tensor in Voigt notation
 * @tparam N The number of support points.
 * @tparam USD The stride one dimensions.
 * @param fieldVar The rank-2 tensor in Voigt notation.
 * @param dNdX The shape function derivatives at this quadrature point.
 * @param detJ The jacobian for the parent space transformation.
 * @param detF The determinant of the deformation gradient.
 * @param fInv The inverse of the deformation gradient.
 * @param result The resulting vector.
 */
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

/**
 * @brief Implements kernels for solving the equations of motion using the
 *   explicit Newmark method under the finite strain assumption.
 * @copydoc ExplicitSmallStrain
 *
 * ### Explicit Small Strain Description
 * Finite strain implementation.
 */
template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
class ExplicitFiniteStrain
  : public ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = ExplicitSmallStrain< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;

  using Base::m_constitutiveUpdate;
  using Base::m_elemGhostRank;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::numNodesPerElem;

  using Base::m_acc;
  using Base::m_dt;
  using Base::m_u;
  using Base::m_vel;
#if !defined( CALCFEMSHAPE )
  using Base::m_detJ;
  using Base::m_dNdX;
#else
  using Base::m_X;
#endif

  /**
   * @copydoc ExplicitSmallStrain
   */
  ExplicitFiniteStrain( NodeManager & nodeManager,
                        EdgeManager const & edgeManager,
                        FaceManager const & faceManager,
                        SUBREGION_TYPE const & elementSubRegion,
                        FE_TYPE const & finiteElementSpace,
                        CONSTITUTIVE_TYPE * const inputConstitutiveType,
                        real64 const dt,
                        string const & elementListName ) :
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          dt,
          elementListName )
  {}

  //*****************************************************************************
  /**
   * @copydoc ExplicitSmallStrain::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    using Base::StackVariables::fLocal;
    using Base::StackVariables::varLocal;

#if defined( CALCFEMSHAPE )
    using Base::StackVariables::detJ;
    using Base::StackVariables::dNdX;
    using Base::StackVariables::xLocal;
#endif

    /**
     * @brief constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables() :
      Base::StackVariables(),
      uLocal { { 0.0 } }
    {}

    /// Local stack storage for nodal displacements.
    real64 uLocal[numNodesPerElem][numDofPerTrialSupportPoint];
  };
  //*****************************************************************************

  /**
   * @copydoc ExplicitSmallStrain::setup
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
        stack.uLocal[a][i] = m_u[nodeIndex][i];
        stack.varLocal[a][i] = m_vel[nodeIndex][i];
      }
    }
  }

  /**
   * @copydoc ExplicitSmallStrain::quadraturePointStateUpdate
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
#else
  /// @cond DOXYGEN_SKIP
  #define DNDX m_dNdX[k][q]
  #define DETJ m_detJ( k, q )
      /// @endcond DOXYGEN_SKIP
#endif
    real64 dUhatdX[3][3], dUdX[3][3];
    CalculateGradients< numNodesPerElem >( dUhatdX,
                                           dUdX,
                                           stack.varLocal,
                                           stack.uLocal,
                                           DNDX );
    LvArray::tensorOps::scale< 3, 3 >( dUhatdX, m_dt );

    real64 F[3][3], Ldt[3][3], fInv[3][3];

    // calculate du/dX
    LvArray::tensorOps::scaledCopy< 3, 3 >( F, dUhatdX, 0.5 );
    LvArray::tensorOps::add< 3, 3 >( F, dUdX );
    LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
    LvArray::tensorOps::invert< 3 >( fInv, F );

    // chain rule: calculate dv/du = dv/dX * dX/du
    LvArray::tensorOps::AikBkj< 3, 3, 3 >( Ldt, dUhatdX, fInv );

    // calculate gradient (end of step)
    LvArray::tensorOps::copy< 3, 3 >( F, dUhatdX );
    LvArray::tensorOps::add< 3, 3 >( F, dUdX );
    LvArray::tensorOps::addIdentity< 3 >( F, 1.0 );
    real64 const detF = LvArray::tensorOps::invert< 3 >( fInv, F );

    real64 Rot[3][3];
    real64 Dadt[6];
    HughesWinget( Rot, Dadt, Ldt );

    m_constitutiveUpdate.HypoElastic( k, q, Dadt, Rot );

    Integrate< numNodesPerElem >( m_constitutiveUpdate.m_stress[k][q].toSliceConst(),
                                  DNDX,
                                  DETJ,
                                  detF,
                                  fInv,
                                  stack.fLocal );
  }
};
#undef CALCFEMSHAPE
#undef DNDX
#undef DETJ
#undef UPDATE_STRESS

}  // namespace SolidMechanicsLagrangianFEMKernels

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSFINITESTRAINEXPLICITNEWMARK_HPP_
