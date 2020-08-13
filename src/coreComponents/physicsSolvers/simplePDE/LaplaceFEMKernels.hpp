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

///*
// * @file SolidMechanicsLagrangianFEMKernels.hpp
// */
//
//#pragma once
//
//#include "common/DataTypes.hpp"
//#include "common/TimingMacros.hpp"
//#include "finiteElement/FiniteElementShapeFunctionKernel.hpp"
//#include "rajaInterface/GEOS_RAJA_Interface.hpp"
//
//namespace geosx
//{
//namespace LaplaceFEMKernels
//{
//
//struct ImplicitKernel
//{
//
//  template< localIndex NUM_NODES_PER_ELEM, localIndex NUM_QUADRATURE_POINTS >
//  static real64 Launch( arrayView4d< real64 const > const & dNdX,
//                        arrayView2d< real64 const > const & detJ,
//                        arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemNodes,
//                        arrayView1d< globalIndex const > const & dofIndex,
//                        globalIndex const dofRankOffset,
//                        CRSMatrixView< real64, globalIndex const > const & matrix )
//  {
//    localIndex const numElems = dNdX.size( 0 );
//    GEOSX_ERROR_IF_NE( dNdX.size( 0 ), numElems );
//    GEOSX_ERROR_IF_NE( dNdX.size( 1 ), NUM_QUADRATURE_POINTS );
//    GEOSX_ERROR_IF_NE( dNdX.size( 2 ), NUM_NODES_PER_ELEM );
//
//    GEOSX_ERROR_IF_NE( detJ.size( 0 ), numElems );
//    GEOSX_ERROR_IF_NE( detJ.size( 1 ), NUM_QUADRATURE_POINTS );
//
//    GEOSX_ERROR_IF_NE( elemNodes.size( 0 ), numElems );
//    GEOSX_ERROR_IF_NE( elemNodes.size( 1 ), NUM_NODES_PER_ELEM );
//
//    // begin element loop, skipping ghost elements
//    forAll< parallelDevicePolicy< 32 > >( numElems, [=] GEOSX_HOST_DEVICE ( localIndex const k )
//    {
//      globalIndex dofIndices[ NUM_NODES_PER_ELEM ];
//      real64 elementMatrix[ NUM_NODES_PER_ELEM ][ NUM_NODES_PER_ELEM ] = { { 0 } };
//
//      for( localIndex q = 0; q < NUM_QUADRATURE_POINTS; ++q )
//      {
//        for( localIndex a = 0; a < NUM_NODES_PER_ELEM; ++a )
//        {
//          dofIndices[ a ] = dofIndex[ elemNodes( k, a ) ];
//
//          real64 diffusion = 1.0;
//          for( localIndex b = 0; b < NUM_NODES_PER_ELEM; ++b )
//          {
//            elementMatrix[ a ][ b ] += detJ( k, q ) *
//                                       diffusion *
//                                       + LvArray::tensorOps::AiBi< 3 >( dNdX[ k ][ q ][ a ], dNdX[ k ][ q ][ b ] );
//          }
//
//        }
//      }
//
//      for( localIndex i = 0; i < NUM_NODES_PER_ELEM; ++i )
//      {
//        globalIndex const dof = dofIndices[ i ] - dofRankOffset;
//        if( dof < 0 || dof >= matrix.numRows() ) continue;
//        matrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
//                                                                     dofIndices,
//                                                                     elementMatrix[ i ],
//                                                                     NUM_NODES_PER_ELEM );
//      }
//    } );
//
//    return 0;
//  }
//};
//
//} // namespace LaplaceFEMKernels
//} // namespace geosx
//=======
/*
 * @file LaplaceFEMKernels.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_

#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"

namespace geosx
{
//*****************************************************************************
/**
 * @brief Implements kernels for solving Laplace's equation.
 * @copydoc geosx::finiteElement::KernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### LaplaceFEMKernel Description
 * Implements the KernelBase interface functions required for solving Laplace's
 * equation using on of the finite element kernel application functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 * In this implementation, the template parameter @p NUM_NODES_PER_ELEM is used
 * in place of both @p NUM_TEST_SUPPORT_POINTS_PER_ELEM and
 * @p NUM_TRIAL_SUPPORT_POINTS_PER_ELEM, which are assumed to be equal. This
 * results in the @p UNUSED template parameter as only the NUM_NODES_PER_ELEM
 * is passed to the ImplicitKernelBase template to form the base class.
 *
 * Additionally, the number of degrees of freedom per support point for both
 * the test and trial spaces are specified as `1` when specifying the base
 * class.
 */
template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
class LaplaceFEMKernel
  : public finiteElement::
      ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 1, 1 >
{
public:
  /// An alias for the base class.
  using Base =
    finiteElement::ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, 1, 1 >;

  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_elemsToNodes;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;

  /// The number of nodes per element.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param fieldName The name of the primary field
   *                  (i.e. Temperature, Pressure, etc.)
   */
  LaplaceFEMKernel( NodeManager const & nodeManager,
                    EdgeManager const & edgeManager,
                    FaceManager const & faceManager,
                    SUBREGION_TYPE const & elementSubRegion,
                    FE_TYPE const & finiteElementSpace,
                    CONSTITUTIVE_TYPE * const inputConstitutiveType,
                    arrayView1d< globalIndex const > const & inputDofNumber,
                    globalIndex const rankOffset,
                    CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                    arrayView1d< real64 > const & inputRhs,
                    string const & fieldName ) :
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
    m_primaryField( nodeManager.template getReference< array1d< real64 > >( fieldName ) ),
    m_dNdX( elementSubRegion.dNdX() ),
    m_detJ( elementSubRegion.detJ() )
  {}

  //***************************************************************************
  /**
   * @class StackVariables
   * @copydoc geosx::finiteElement::ImplicitKernelBase::StackVariables
   *
   * Adds a stack array for the primary field.
   */
  struct StackVariables : Base::StackVariables
  {
  public:
    /**
     * @brief Constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables() :
      Base::StackVariables(),
      primaryField_local { 0.0 }
    {}

    /// C-array storage for the element local primary field variable.
    real64 primaryField_local[numNodesPerElem];
  };

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the LaplaceFEMKernel implementation, global values from the
   * primaryField, and degree of freedom numbers are placed into element local
   * stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  setup( localIndex const k, StackVariables & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      stack.primaryField_local[a] = m_primaryField[localNodeIndex];
      stack.localRowDofIndex[a] = m_dofNumber[localNodeIndex];
      stack.localColDofIndex[a] = m_dofNumber[localNodeIndex];
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::quadraturePointJacobianContribution
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  quadraturePointJacobianContribution( localIndex const k,
                                       localIndex const q,
                                       StackVariables & stack ) const
  {
    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        stack.localJacobian[a][b] +=
          LvArray::tensorOps::AiBi< 3 >( m_dNdX[k][q][a], m_dNdX[k][q][b] ) *
          m_detJ( k, q );
      }
    }
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   *
   * Form element residual from the fully formed element Jacobian dotted with
   * the primary field and map the element local Jacobian/Residual to the
   * global matrix/vector.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64
  complete( localIndex const k, StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );
    real64 maxForce = 0;

    for( localIndex a = 0; a < numNodesPerElem; ++a )
    {
      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        stack.localResidual[a] +=
          stack.localJacobian[a][b] * stack.primaryField_local[b];
      }
    }

    for( int a = 0; a < numNodesPerElem; ++a )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >(
        stack.localRowDofIndex[a] - m_dofRankOffset );
      if( dof < 0 || dof >= m_matrix.numRows() ) continue;
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
        dof,
        stack.localColDofIndex,
        stack.localJacobian[a],
        numNodesPerElem );

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidual[a] );
      maxForce = fmax( maxForce, fabs( stack.localResidual[a] ) );
    }

    return maxForce;
  }

protected:
  /// The global primary field array.
  arrayView1d< real64 const > const m_primaryField;

  /// The global shape function derivatives array.
  arrayView4d< real64 const > const m_dNdX;

  /// The global determinant of the parent/physical Jacobian.
  arrayView2d< real64 const > const m_detJ;
};

}  // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif  // GEOSX_PHYSICSSOLVERS_SIMPLEPDE_LAPLACEFEMKERNELS_HPP_
