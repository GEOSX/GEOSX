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
 * @file ThermalDiffusionKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSIONKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSIONKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

namespace geosx
{

template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class ThermalDiffusionKernel :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            1,
                                            1 >
{
public:

  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  1,
                                                  1 >;

  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_finiteElementSpace;

  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  ThermalDiffusionKernel( NodeManager const & nodeManager,
                          EdgeManager const & edgeManager,
                          FaceManager const & faceManager,
                          localIndex const targetRegionIndex,
                          SUBREGION_TYPE const & elementSubRegion,
                          FE_TYPE const & finiteElementSpace,
                          CONSTITUTIVE_TYPE & inputConstitutiveType,
                          arrayView1d< globalIndex const > const & inputDofNumber,
                          globalIndex const rankOffset,
                          CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                          arrayView1d< real64 > const & inputRhs,
                          real64 const & thermalDiffusion,
                          real64 const dt ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition() ),
    m_temperature( nodeManager.template getReference< array1d< real64 > >( dataRepository::keys::Temperature ) ),
    m_deltaTemperature( nodeManager.template getReference< array1d< real64 > >( dataRepository::keys::IncrementalTemperature ) ),
    m_thermalDiffusion( thermalDiffusion ),
    m_dt( dt )
  {}

  struct StackVariables : Base::StackVariables
  {
public:

    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            temperature_local{ 0.0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// C-array storage for the element local temperature field variable.
    real64 temperature_local[numNodesPerElem];

    /// C-array storage for the element local incremental temperature field variable.
    real64 deltaTemperature_local[numNodesPerElem];
  };

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

#if defined(CALC_FEM_SHAPE_IN_KERNEL)
      for( int i=0; i<3; ++i )
      {
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
      }
#endif

      stack.temperature_local[a] = m_temperature[localNodeIndex];
      stack.deltaTemperature_local[a] = m_deltaTemperature[localNodeIndex];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    // Shape function
    real64 N[numNodesPerElem];
    FE_TYPE::calcN( q, N );

    // Gradient of shape function
    real64 dNdX[numNodesPerElem][3];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      for( localIndex b=0; b<numNodesPerElem; ++b )
      {
        real64 tmpValue1 = N[a] * N[b] * detJ / m_dt;
        real64 tmpValue2 = LvArray::tensorOps::AiBi< 3 >( dNdX[a], dNdX[b] ) * detJ * m_thermalDiffusion;

        // Update local Jacobian
        stack.localJacobian[a][b] += tmpValue1 + tmpValue2;

        // Update local Rhs, the sign of Rhs is inverted here
        stack.localResidual[a] += tmpValue1 * stack.deltaTemperature_local[b];
        stack.localResidual[a] += tmpValue2 * stack.temperature_local[b];
      }
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[a] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.localColDofIndex,
                                                                              stack.localJacobian[a],
                                                                              numNodesPerElem );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidual[a] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[a] ) );
    }

    return maxForce;
  }

protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The global temperature field array.
  arrayView1d< real64 const > const m_temperature;

  /// The global incremental temperature field array.
  arrayView1d< real64 const > const m_deltaTemperature;

  /// Thermal diffusion coefficient
  real64 const m_thermalDiffusion;

  /// Time step
  real64 const m_dt;

};

/// The factory used to construct a ThermalDiffusionKernel.
using ThermalDiffusionKernelFactory =
  finiteElement::KernelFactory< ThermalDiffusionKernel,
                                arrayView1d< globalIndex const > const &,
                                globalIndex, CRSMatrixView< real64, globalIndex const > const &,
                                arrayView1d< real64 > const &,
                                real64 const &,
                                real64 const >;

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSIONKERNELS_HPP_*/
