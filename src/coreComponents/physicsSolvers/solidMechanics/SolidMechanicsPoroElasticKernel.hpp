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
 * @file SolidMechanicsPoroElasticKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOROELASTICKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOROELASTICKERNEL_HPP_

#include "physicsSolvers/solidMechanics/SolidMechanicsSmallStrainQuasiStaticKernel.hpp"

namespace geosx
{

namespace SolidMechanicsLagrangianFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static equilibrium.
 * @copydoc QuasiStatic
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 * @tparam BASE The base class kernel for this kernel. Likely one of the
 *   existing SolidMechanicsKernel implementations.
 *
 * ### PoroElastic Description
 * Implements the KernelBase interface functions required for using the
 * effective stress for the integration of the stress divergence. This is
 * templated on one of the existing solid mechanics kernel implementations
 * such as QuasiStatic.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE,
          template< typename,
                    typename,
                    typename > class BASE >
class PoroElastic : public BASE< SUBREGION_TYPE,
                                 CONSTITUTIVE_TYPE,
                                 FE_TYPE >
{
public:
  /// Alias for the base class.
  using Base = BASE< SUBREGION_TYPE,
                     CONSTITUTIVE_TYPE,
                     FE_TYPE >;

  using Base::m_constitutiveUpdate;
  using typename Base::StackVariables;

  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
  PoroElastic( NodeManager const & nodeManager,
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
               real64 const (&inputGravityVector)[3] ):
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
          inputRhs,
          inputGravityVector ),
    m_fluidPressure( elementSubRegion.template getReference< array1d< real64 > >( "pressure" ) ),
    m_deltaFluidPressure( elementSubRegion.template getReference< array1d< real64 > >( "deltaPressure" ) )
  {}



  /**
   * @copydoc geosx::finiteElement::KernelBase::quadraturePointKernel
   *
   * The divergence of the effective stress is integrated over the volume of
   * the element, yielding the nodal force (residual) contributions.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {
    real64 const biotCoefficient = 1.0;//m_constitutiveUpdate.getBiotCoefficient();
    Base::quadraturePointKernel( k, q, stack, [=] GEOSX_HOST_DEVICE ( real64 (& stress)[6] )
    {
      real64 const pressure = biotCoefficient * ( m_fluidPressure[k] + m_deltaFluidPressure[k] );
      stress[0] -= pressure;
      stress[1] -= pressure;
      stress[2] -= pressure;
    } );
  }


protected:
  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;
};

/**
 * Alias for the PoroElastic kernel with the QuasiStatic base class.
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
using QuasiStaticPoroElastic = PoroElastic< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            QuasiStatic >;

/// The factory used to construct a QuasiStaticPoroElastic kernel.
using QuasiStaticPoroElasticFactory = finiteElement::KernelFactory< QuasiStaticPoroElastic,
                                                                    arrayView1d< globalIndex const > const &,
                                                                    globalIndex,
                                                                    CRSMatrixView< real64, globalIndex const > const &,
                                                                    arrayView1d< real64 > const &,
                                                                    real64 const (&)[3] >;

} // namespace SolidMechanicsLagrangianFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSPOROELASTICKERNEL_HPP_
