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

#include "KernelBase.hpp"

/**
 * @file ImplicitKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_
  #define GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_

namespace geosx
{
namespace finiteElement
{
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class ImplicitKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geosx::finiteElement::KernelBase
 *
 * ### ImplicitKernelBase Description
 * Provides a common base for kernels that require the assembly of a system of
 * equations. The types required to assemble the system, such as DOF
 * information, the Matrix and Vector object, etc., are declared and set here.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class ImplicitKernelBase
  : public KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, NUM_DOF_PER_TEST_SP, NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geosx::finiteElement::KernelBase)
  using Base =
    KernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, NUM_DOF_PER_TEST_SP, NUM_DOF_PER_TRIAL_SP >;

  using Base::m_elemsToNodes;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::numTestSupportPointsPerElem;
  using Base::numTrialSupportPointsPerElem;

  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param inputDofNumber The dof number for the primary field.
   * @param rankOffset dof index offset of current rank
   * @param inputMatrix Reference to the Jacobian matrix.
   * @param inputRhs Reference to the RHS vector.
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   */
  ImplicitKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      SUBREGION_TYPE const & elementSubRegion,
                      FE_TYPE const & finiteElementSpace,
                      CONSTITUTIVE_TYPE * const inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      globalIndex const rankOffset,
                      CRSMatrixView< real64, globalIndex const > const & inputMatrix,
                      arrayView1d< real64 > const & inputRhs ) :
    Base( elementSubRegion, finiteElementSpace, inputConstitutiveType ),
    m_dofNumber( inputDofNumber ),
    m_dofRankOffset( rankOffset ),
    m_matrix( inputMatrix ),
    m_rhs( inputRhs )
  {
    GEOSX_UNUSED_VAR( nodeManager );
    GEOSX_UNUSED_VAR( edgeManager );
    GEOSX_UNUSED_VAR( faceManager );
  }

  //***************************************************************************
  /**
   * @copydoc finiteElement::KernelBase::StackVariables
   */
  struct StackVariables : public Base::StackVariables
  {
    /// The number of rows in the element local jacobian matrix.
    static constexpr int numRows =
      numTestSupportPointsPerElem * numDofPerTestSupportPoint;

    /// The number of columns in the element local jacobian matrix.
    static constexpr int numCols =
      numTrialSupportPointsPerElem * numDofPerTrialSupportPoint;

    /**
     * Default constructor
     */
    GEOSX_HOST_DEVICE
    StackVariables() :
      localRowDofIndex { 0 },
      localColDofIndex { 0 },
      localResidual { 0.0 },
      localJacobian { { 0.0 } }
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localRowDofIndex[numRows];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex localColDofIndex[numCols];

    /// C-array storage for the element local residual vector.
    real64 localResidual[numRows];

    /// C-array storage for the element local Jacobian matrix.
    real64 localJacobian[numRows][numCols];
  };
  //***************************************************************************

  /**
   * @copydoc geosx::finiteElement::KernelBase::setup
   *
   * ### ImplicitKernelBase::setup() Description
   *
   * In this implementation, the element local Row and Column DOF stack arrays
   * are filled for when we fill the global matrix and rhs.
   *
   * @note This seems like a waste of register space. We should do this in
   *       complete() unless we actually need these dof somewhere else in the kernel.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  setup( localIndex const k, StackVariables & stack ) const
  {
    for( localIndex a = 0; a < numTestSupportPointsPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i = 0; i < numDofPerTestSupportPoint; ++i )
      {
        stack.localRowDofIndex[a * numDofPerTestSupportPoint + i] =
          m_dofNumber[localNodeIndex] + i;
      }
    }

    for( localIndex a = 0; a < numTrialSupportPointsPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes[k][a];
      for( int i = 0; i < numDofPerTrialSupportPoint; ++i )
      {
        stack.localColDofIndex[a * numDofPerTrialSupportPoint + i] =
          m_dofNumber[localNodeIndex] + i;
      }
    }
  }

protected:
  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_dofNumber;

  /// The global rank offset
  globalIndex const m_dofRankOffset;

  /// The global Jacobian matrix.
  CRSMatrixView< real64, globalIndex const > const m_matrix;

  /// The global residaul vector.
  arrayView1d< real64 > const m_rhs;
};

}  // namespace finiteElement
}  // namespace geosx

#endif /* GEOSX_FINITEELEMENT_IMPLICITKERNELBASE_HPP_ */
