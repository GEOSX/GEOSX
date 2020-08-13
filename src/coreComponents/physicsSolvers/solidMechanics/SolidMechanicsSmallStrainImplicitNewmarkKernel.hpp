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
 * @file SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_

#include "SolidMechanicsSmallStrainQuasiStaticKernel.hpp"

namespace geosx
{
namespace SolidMechanicsLagrangianFEMKernels
{
/**
 * @brief Implements kernels for solving the equations of motion using an
 *   implicit Newmark's method..
 * @copydoc QuasiStatic
 *
 * ### Implicit Newmark Description
 * Implements the KernelBase interface functions required for solving the
 * equations of motion using with an Implicit Newmark's Method with one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
class ImplicitNewmark
  : public QuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >
{
public:
  /// Alias for the base class;
  using Base = QuasiStatic< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::numNodesPerElem;
  using Base::numTestSupportPointsPerElem;
  using Base::numTrialSupportPointsPerElem;

  using Base::m_density;
  using Base::m_detJ;
  using Base::m_disp;
  using Base::m_dofNumber;
  using Base::m_elemsToNodes;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_uhat;

  /**
   * @brief Constructor
   * @copydoc QuasiStatic
   * @param inputNewmarkGamma The Gamma parameter of the Newmark method.
   * @param inputNewmarkBeta The Beta parameter for the Newmark method.
   * @param inputMassDamping The mass damping coefficient.
   * @param inputStiffnessDamping The stiffness damping coefficient.
   * @param inputDt The timestep for the physics update.
   */
  ImplicitNewmark( NodeManager const & nodeManager,
                   EdgeManager const & edgeManager,
                   FaceManager const & faceManager,
                   SUBREGION_TYPE const & elementSubRegion,
                   FE_TYPE const & finiteElementSpace,
                   CONSTITUTIVE_TYPE * const inputConstitutiveType,
                   arrayView1d< globalIndex const > const & inputDofNumber,
                   globalIndex const rankOffset,
                   CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                   arrayView1d< real64 > const & inputRhs,
                   real64 const ( &inputGravityVector )[3],
                   real64 const inputNewmarkGamma,
                   real64 const inputNewmarkBeta,
                   real64 const inputMassDamping,
                   real64 const inputStiffnessDamping,
                   real64 const inputDt ) :
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector ),
    m_vtilde( nodeManager.totalDisplacement() ),
    m_uhattilde( nodeManager.totalDisplacement() ),
    m_newmarkGamma( inputNewmarkGamma ),
    m_newmarkBeta( inputNewmarkBeta ),
    m_massDamping( inputMassDamping ),
    m_stiffnessDamping( inputStiffnessDamping ),
    m_dt( inputDt )
  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc QuasiStatic::StackVariables
   *
   * Adds a stack array for the vtilde, uhattilde, and the
   * Inertial mass damping.
   */
  struct StackVariables : public Base::StackVariables
  {
  public:
    using Base::StackVariables::numCols;
    using Base::StackVariables::numRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables() :
      Base::StackVariables(),
      dRdU_InertiaMassDamping { { 0.0 } },
      vtilde_local(),
      uhattilde_local()
    {}

    /// Stack storage for the Inertial damping contributions to the Jacobian
    real64 dRdU_InertiaMassDamping[numRows][numCols];

    /// Stack storage for the velocity predictor.
    real64 vtilde_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the incremental displacement predictor.
    real64 uhattilde_local[numNodesPerElem][numDofPerTrialSupportPoint];
  };
  //***************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc QuasiStatic::setup
   *
   * For the ImplicitNewmark implementation, global values from the velocity
   * predictor, and the incremental displacement predictor are placed into
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
      for( localIndex i = 0; i < numDofPerTrialSupportPoint; ++i )
      {
        stack.vtilde_local[a][i] = m_vtilde[localNodeIndex][i];
        stack.uhattilde_local[a][i] = m_uhattilde[localNodeIndex][i];
      }
    }
    Base::setup( k, stack );
  }

  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointJacobianContribution
   *
   * The ImplcitNewmark kernel adds the calculation of the inertia damping,
   * jacobian and residual contributions.
   */
  GEOSX_DEVICE
  GEOSX_FORCE_INLINE
  void
  quadraturePointJacobianContribution( localIndex const k,
                                       localIndex const q,
                                       StackVariables & stack ) const
  {
    real64 N[numNodesPerElem];
    FE_TYPE::shapeFunctionValues( q, N );

    Base::quadraturePointJacobianContribution(
      k,
      q,
      stack,
      [&] GEOSX_DEVICE( localIndex const a, localIndex const b ) mutable {
        real64 const integrationFactor =
          m_density( k, q ) * N[a] * N[b] * m_detJ( k, q );
        real64 const temp1 =
          ( m_massDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) +
            1.0 / ( m_newmarkBeta * m_dt * m_dt ) ) *
          integrationFactor;

        constexpr int nsdof = numDofPerTestSupportPoint;
        for( int i = 0; i < nsdof; ++i )
        {
          realT const acc = 1.0 / ( m_newmarkBeta * m_dt * m_dt ) *
            ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i] );
          realT const vel = stack.vtilde_local[b][i] +
            m_newmarkGamma / ( m_newmarkBeta * m_dt ) *
              ( stack.uhat_local[b][i] - stack.uhattilde_local[b][i] );

          stack.dRdU_InertiaMassDamping[a * nsdof + i][b * nsdof + i] -= temp1;
          stack.localResidual[a * nsdof + i] -=
            ( m_massDamping * vel + acc ) * integrationFactor;
        }
      } );
  }

  /**
   * @copydoc QuasiStatic::complete
   *
   * The ImplicitNewmark implementation adds residual and jacobian
   * contributions from  stiffness based damping.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64
  complete( localIndex const k, StackVariables & stack ) const
  {
    for( int a = 0; a < numNodesPerElem; ++a )
    {
      for( int b = 0; b < numNodesPerElem; ++b )
      {
        for( int i = 0; i < numDofPerTestSupportPoint; ++i )
        {
          for( int j = 0; j < numDofPerTrialSupportPoint; ++j )
          {
            stack.localResidual[a * numDofPerTestSupportPoint + i] =
              stack.localResidual[a * numDofPerTestSupportPoint + i] +
              m_stiffnessDamping *
                stack.localJacobian[a * numDofPerTestSupportPoint + i]
                                   [b * numDofPerTrialSupportPoint + j] *
                ( stack.vtilde_local[b][j] +
                  m_newmarkGamma / ( m_newmarkBeta * m_dt ) *
                    ( stack.uhat_local[b][j] - stack.uhattilde_local[b][j] ) );

            stack.localJacobian[a * numDofPerTestSupportPoint + i]
                               [b * numDofPerTrialSupportPoint + j] =
              stack.localJacobian[a * numDofPerTestSupportPoint + i]
                                 [b * numDofPerTrialSupportPoint + j] +
              stack.localJacobian[a][b] *
                ( 1.0 +
                  m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) ) +
              stack.dRdU_InertiaMassDamping[a * numDofPerTestSupportPoint + i]
                                           [b * numDofPerTrialSupportPoint + j];
          }
        }
      }
    }

    for( int a = 0; a < stack.numRows; ++a )
    {
      for( int b = 0; b < stack.numCols; ++b )
      {
        stack.localJacobian[a][b] += stack.localJacobian[a][b] *
            ( 1.0 + m_stiffnessDamping * m_newmarkGamma / ( m_newmarkBeta * m_dt ) ) +
          stack.dRdU_InertiaMassDamping[a][b];
      }
    }

    return Base::complete( k, stack );
  }

protected:
  /// The rank-global velocity predictor
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_vtilde;

  /// The rank-global incremental displacement predictor
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhattilde;

  /// The Gamma parameter for Newmark's method.
  real64 const m_newmarkGamma;

  /// The Beta parameter for Newmark's method.
  real64 const m_newmarkBeta;

  /// The mass damping coefficient.
  real64 const m_massDamping;

  /// The stiffness damping coefficient.
  real64 const m_stiffnessDamping;

  /// The timestep for the update.
  real64 const m_dt;
};

}  // namespace SolidMechanicsLagrangianFEMKernels

}  // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif  //GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSSMALLSTRAINIMPLICITNEWMARK_HPP_
