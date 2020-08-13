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
 * @file SolidMechanicsSmallStrainQuasiStaticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{
namespace SolidMechanicsLagrangianFEMKernels
{
/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### QuasiStatic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static equilibrium equations using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `3` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
class QuasiStatic
  : public finiteElement::
      ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 3, 3 >
{
public:
  /// Alias for the base class;
  using Base =
    finiteElement::ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 3, 3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  using Base::m_constitutiveUpdate;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_elemsToNodes;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  QuasiStatic( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE * const inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDofNumber,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const & inputMatrix,
               arrayView1d< real64 > const & inputRhs,
               real64 const ( &inputGravityVector )[3] ) :
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_disp( nodeManager.totalDisplacement() ),
    m_uhat( nodeManager.incrementalDisplacement() ),
    m_dNdX( elementSubRegion.dNdX() ),
    m_detJ( elementSubRegion.detJ() ),
    m_gravityVector { inputGravityVector[0],
                      inputGravityVector[1],
                      inputGravityVector[2] },
    m_density( inputConstitutiveType->getDensity() )
  {}

  //*****************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the displacement, incremental displacement, and the
   * constitutive stiffness.
   */
  struct StackVariables : public Base::StackVariables
  {
  public:
    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables() :
      Base::StackVariables(),
      u_local(),
      uhat_local(),
      constitutiveStiffness { { 0.0 } }
    {}

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[6][6];
  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the QuasiStatic implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  setup( localIndex const k, StackVariables & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i = 0; i < 3; ++i )
      {
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a * 3 + i] = m_dofNumber[localNodeIndex] + i;
        stack.localColDofIndex[a * 3 + i] = m_dofNumber[localNodeIndex] + i;
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointStateUpdate
   *
   * For solid mechanics kernels, the strain increment is calculated, and the
   * constitutive update is called. In addition, the constitutive stiffness
   * stack variable is filled by the constitutive model.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  quadraturePointStateUpdate( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 strainInc[6] = { 0 };
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      strainInc[0] = strainInc[0] + m_dNdX( k, q, a, 0 ) * stack.uhat_local[a][0];
      strainInc[1] = strainInc[1] + m_dNdX( k, q, a, 1 ) * stack.uhat_local[a][1];
      strainInc[2] = strainInc[2] + m_dNdX( k, q, a, 2 ) * stack.uhat_local[a][2];
      strainInc[3] = strainInc[3] + m_dNdX( k, q, a, 2 ) * stack.uhat_local[a][1] +
        m_dNdX( k, q, a, 1 ) * stack.uhat_local[a][2];

      strainInc[4] = strainInc[4] + m_dNdX( k, q, a, 2 ) * stack.uhat_local[a][0] +
        m_dNdX( k, q, a, 0 ) * stack.uhat_local[a][2];

      strainInc[5] = strainInc[5] + m_dNdX( k, q, a, 1 ) * stack.uhat_local[a][0] +
        m_dNdX( k, q, a, 0 ) * stack.uhat_local[a][1];
    }

    m_constitutiveUpdate.SmallStrain( k, q, strainInc );

    GEOSX_UNUSED_VAR( q )
    m_constitutiveUpdate.GetStiffness( k, stack.constitutiveStiffness );
  }

  /**
   * @brief Internal struct to provide no-op defaults used in the inclusion
   *   of lambda functions into kernel component functions.
   * @struct NoOpFunctors
   */
  struct NoOpFunctors
  {
    /**
     * @brief operator() no-op used for adding an additional dynamics term
     *   inside the jacobian assembly loop.
     * @param a Node index for the row.
     * @param b Node index for the col.
     */
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr void
    operator()(
      localIndex const a,
      localIndex const b )
    {
      GEOSX_UNUSED_VAR( a );
      GEOSX_UNUSED_VAR( b );
    }

    /**
     * @brief operator() no-op used for modifying the stress tensor prior to
     *   integrating the divergence to produce nodal forces.
     * @param stress The stress array.
     */
    GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE constexpr void operator()(
      real64 ( &stress )[6] )
    {
      GEOSX_UNUSED_VAR( stress );
    }
  };

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointJacobianContribution
   * @tparam DYNAMICS_LAMBDA The type of the lambda that will generate dynamics
   *   terms inside the Jacobian loop.
   * @param dynamicsTerms The lambda that generates dynamics terms.
   * For solid mechanics kernels, the derivative of the force residual wrt
   * the incremental displacement is filled into the local element jacobian.
   */
  template< typename DYNAMICS_LAMBDA = NoOpFunctors >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
  quadraturePointJacobianContribution(
    localIndex const k,
    localIndex const q,
    StackVariables & stack,
    DYNAMICS_LAMBDA && dynamicsTerms = NoOpFunctors {} ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        real64 const( &c )[6][6] = stack.constitutiveStiffness;
        stack.localJacobian[a * 3 + 0][b * 3 + 0] -=
          ( c[0][0] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 0 ) +
            c[5][5] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 1 ) +
            c[4][4] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 0][b * 3 + 1] -=
          ( c[5][5] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 0 ) +
            c[0][1] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 1 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 0][b * 3 + 2] -=
          ( c[4][4] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 0 ) +
            c[0][2] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 1][b * 3 + 1] -=
          ( c[5][5] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 0 ) +
            c[1][1] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 1 ) +
            c[3][3] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 1][b * 3 + 0] -=
          ( c[0][1] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 0 ) +
            c[5][5] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 1 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 1][b * 3 + 2] -=
          ( c[3][3] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 1 ) +
            c[1][2] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 2][b * 3 + 0] -=
          ( c[0][2] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 0 ) +
            c[4][4] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 2][b * 3 + 1] -=
          ( c[1][2] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 1 ) +
            c[3][3] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        stack.localJacobian[a * 3 + 2][b * 3 + 2] -=
          ( c[4][4] * m_dNdX( k, q, a, 0 ) * m_dNdX( k, q, b, 0 ) +
            c[3][3] * m_dNdX( k, q, a, 1 ) * m_dNdX( k, q, b, 1 ) +
            c[2][2] * m_dNdX( k, q, a, 2 ) * m_dNdX( k, q, b, 2 ) ) *
          m_detJ( k, q );

        dynamicsTerms( a, b );
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointResidualContribution
   * @tparam STRESS_MODIFIER The type of the function to modify the stress
   *   prior to integration.
   * @param stressModifier The stress modifier functor/lambda.
   *
   * The divergence of the stress is integrated over the volume of the element,
   * yielding the nodal force (residual) contributions.
   */
  template< typename STRESS_MODIFIER = NoOpFunctors >
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void
  quadraturePointResidualContribution(
    localIndex const k,
    localIndex const q,
    StackVariables & stack,
    STRESS_MODIFIER && stressModifier = NoOpFunctors {} ) const
  {
    real64 stress[6] = { m_constitutiveUpdate.m_stress( k, q, 0 ),
                         m_constitutiveUpdate.m_stress( k, q, 1 ),
                         m_constitutiveUpdate.m_stress( k, q, 2 ),
                         m_constitutiveUpdate.m_stress( k, q, 3 ),
                         m_constitutiveUpdate.m_stress( k, q, 4 ),
                         m_constitutiveUpdate.m_stress( k, q, 5 ) };

    stressModifier( stress );

    real64 const gravityForce[3] = { m_gravityVector[0] * m_density( k, q ),
                                     m_gravityVector[1] * m_density( k, q ),
                                     m_gravityVector[2] * m_density( k, q ) };

    real64 N[numNodesPerElem];
    FE_TYPE::shapeFunctionValues( q, N );
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      stack.localResidual[a * 3 + 0] -=
        ( stress[0] * m_dNdX( k, q, a, 0 ) + stress[5] * m_dNdX( k, q, a, 1 ) +
          stress[4] * m_dNdX( k, q, a, 2 ) - gravityForce[0] * N[a] ) *
        m_detJ( k, q );
      stack.localResidual[a * 3 + 1] -=
        ( stress[5] * m_dNdX( k, q, a, 0 ) + stress[1] * m_dNdX( k, q, a, 1 ) +
          stress[3] * m_dNdX( k, q, a, 2 ) - gravityForce[1] * N[a] ) *
        m_detJ( k, q );
      stack.localResidual[a * 3 + 2] -=
        ( stress[4] * m_dNdX( k, q, a, 0 ) + stress[3] * m_dNdX( k, q, a, 1 ) +
          stress[2] * m_dNdX( k, q, a, 2 ) - gravityForce[2] * N[a] ) *
        m_detJ( k, q );
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64
  complete( localIndex const k, StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >(
          stack.localRowDofIndex[numDofPerTestSupportPoint * localNode + dim] -
          m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
          dof,
          stack.localRowDofIndex,
          stack.localJacobian[numDofPerTestSupportPoint * localNode + dim],
          numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >(
          &m_rhs[dof],
          stack.localResidual[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax(
          maxForce,
          fabs( stack.localResidual[numDofPerTestSupportPoint * localNode + dim] ) );
      }
    }

    return maxForce;
  }

protected:
  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The shape function derivative for each quadrature point.
  arrayView4d< real64 const > const m_dNdX;

  /// The parent->physical jacobian determinant for each quadrature point.
  arrayView2d< real64 const > const m_detJ;

  /// The gravity vector.
  real64 const m_gravityVector[3];

  /// The rank global density
  arrayView2d< real64 const > const m_density;
};

}  // namespace SolidMechanicsLagrangianFEMKernels

}  // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif  // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINQUASISTATIC_HPP_
