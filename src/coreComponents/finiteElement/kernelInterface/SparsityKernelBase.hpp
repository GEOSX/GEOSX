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

#include "ImplicitKernelBase.hpp"

/**
 * @file SparsityKernelBase.hpp
 */

#ifndef GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_
  #define GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_

namespace geosx
{
namespace finiteElement
{
//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @class SparsityKernelBase
 * @brief Define the base interface for implicit finite element kernels.
 * @copydoc geosx::finiteElement::KernelBase
 *
 * ### SparsityKernelBase Description
 * Provide common kernels for generation of CRS Sparsity patterns.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          int NUM_DOF_PER_TEST_SP,
          int NUM_DOF_PER_TRIAL_SP >
class SparsityKernelBase
  : public ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, NUM_DOF_PER_TEST_SP, NUM_DOF_PER_TRIAL_SP >
{
public:
  /// Alias for the base class. (i.e. #geosx::finiteElement::ImplicitKernelBase)
  using Base =
    ImplicitKernelBase< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE, NUM_DOF_PER_TEST_SP, NUM_DOF_PER_TRIAL_SP >;

  using Base::m_dofRankOffset;
  using typename Base::StackVariables;

  using Base::setup;

  /**
   * @brief Constructor
   * @param nodeManager Reference to the NodeManager object.
   * @param edgeManager Reference to the EdgeManager object.
   * @param faceManager Reference to the FaceManager object.
   * @param inputDofNumber The dof number for the primary field.
   * @param rankOffset dof index offset of current rank
   * @param inputSparsity The sparsity pattern to fill.
   * @param rowSizes The array that will be filled with row sizes.
   * @copydoc geosx::finiteElement::KernelBase::KernelBase
   */
  SparsityKernelBase( NodeManager const & nodeManager,
                      EdgeManager const & edgeManager,
                      FaceManager const & faceManager,
                      SUBREGION_TYPE const & elementSubRegion,
                      FE_TYPE const & finiteElementSpace,
                      CONSTITUTIVE_TYPE * const inputConstitutiveType,
                      arrayView1d< globalIndex const > const & inputDofNumber,
                      globalIndex const rankOffset,
                      SparsityPattern< globalIndex > & inputSparsity,
                      arrayView1d< localIndex > const & rowSizes ) :
    Base( nodeManager,
          edgeManager,
          faceManager,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDofNumber,
          rankOffset,
          CRSMatrix< real64, globalIndex >().toViewConstSizes(),
          array1d< real64 >().toView() ),
    m_sparsity( inputSparsity ),
    m_rowSizes( rowSizes )
  {}

  /**
   * @copydoc geosx::finiteElement::KernelBase::complete
   *
   * In this implementation, only the matrix values are inserted, making this
   * implementation appropriate for generating the sparsity pattern.
   */
  GEOSX_FORCE_INLINE
  real64
  complete( localIndex const k, StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );

    for( localIndex r = 0; r < stack.numRows; ++r )
    {
      localIndex const row = stack.localRowDofIndex[r] - m_dofRankOffset;
      if( row < 0 || row >= m_sparsity.numRows() ) continue;
      for( localIndex c = 0; c < stack.numCols; ++c )
      {
        m_sparsity.insertNonZero( row, stack.localColDofIndex[c] );
      }
    }
    return 0;
  }

private:
  SparsityPattern< globalIndex > & m_sparsity;

  arrayView1d< localIndex > const & m_rowSizes;
};

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Helper struct to define a specialization of
 *   #::geosx::finiteElement::SparsityKernelBase that may be used to generate the sparsity pattern.
 * @tparam KERNEL_TEMPLATE Templated class that defines the physics kernel.
 *                         Most likely derives from SparsityKernelBase.
 */
template< template< typename, typename, typename > class KERNEL_TEMPLATE >
struct SparsityHelper
{
  /**
   * Defines an alias for the specialization of
   * #geosx::finiteElement::SparsityKernelBase from the compile time constants
   * defined in @p KERNEL_TEMPLATE. Specifically, the
   * NUM_TEST_SUPPORT_POINTS_PER_ELEM and NUM_TRIAL_SUPPORT_POINTS_PER_ELEM
   * parameters are specified by the physics solver.
   */
  template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
  using Kernel = SparsityKernelBase<
    SUBREGION_TYPE,
    CONSTITUTIVE_TYPE,
    FE_TYPE,
    KERNEL_TEMPLATE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::numDofPerTestSupportPoint,
    KERNEL_TEMPLATE< SUBREGION_TYPE, CONSTITUTIVE_TYPE, FE_TYPE >::numDofPerTrialSupportPoint >;
};

//*****************************************************************************
//*****************************************************************************
//*****************************************************************************
/**
 * @brief Fills matrix sparsity.
 * @tparam POLICY The RAJA launch policy to pass to the kernel launch.
 * @tparam REGION_TYPE The type of region to loop over. TODO make this a
 *                     parameter pack?
 * @tparam KERNEL_TEMPLATE The type of template for the physics kernel, which
 *                         conforms to the interface specified by KernelBase.
 * @param mesh The MeshLevel object.
 * @param targetRegions The names of the target regions(of type @p REGION_TYPE)
 *                      to apply the @p KERNEL_TEMPLATE.
 * @param discretizationName The name of the finite element discretization.
 * @param inputDofNumber The global degree of freedom numbers.
 * @param rankOffset Offset of dof indices on curren rank.
 * @param inputSparsityPattern The local sparsity pattern to fill.
 * @param rowSizes The array of local row sizes to be populated
 * @return 0
 *
 * Fills matrix sparsity using information from physics specific implementation
 * of #geosx::finiteElement::KernelBase interface using the
 * #geosx::finiteElement::ImplicitKernelBase alias to specialize
 * #geosx::finiteElement::SparsityKernelBase to conform with the template
 * pattern specified in the physics kernels.
 */
template< typename POLICY,
          typename REGION_TYPE,
          template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
          class KERNEL_TEMPLATE >
static real64
fillSparsity( MeshLevel & mesh,
              arrayView1d< string const > const & targetRegions,
              string const & discretizationName,
              arrayView1d< globalIndex const > const & inputDofNumber,
              globalIndex const rankOffset,
              SparsityPattern< globalIndex > & inputSparsityPattern,
              arrayView1d< localIndex > const & rowSizes )
{
  regionBasedKernelApplication< POLICY,
                                constitutive::NullModel,
                                REGION_TYPE,
                                SparsityHelper< KERNEL_TEMPLATE >::template Kernel >(
    mesh,
    targetRegions,
    discretizationName,
    array1d< string >(),
    inputDofNumber,
    rankOffset,
    inputSparsityPattern,
    rowSizes );

  return 0;
}

}  // namespace finiteElement
}  // namespace geosx

#endif /* GEOSX_FINITEELEMENT_SPARSITYKERNELBASE_HPP_ */
