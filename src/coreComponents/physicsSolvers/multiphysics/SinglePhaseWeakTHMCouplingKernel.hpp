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
 * @file SinglePhaseWeakTHMCouplingKernel.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEWEAKTHMCOUPLINGKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEWEAKTHMCOUPLINGKERNEL_HPP_

#include "SinglePhasePoromechanicsKernel.hpp"
#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

namespace geosx
{

namespace THMKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase THM problems with weak thermal coupling.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhase :
  public PoromechanicsKernels::SinglePhase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE >
{
public:

  using Base = PoromechanicsKernels::SinglePhase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE >;

  using Base::m_constitutiveUpdate;
  using typename Base::StackVariables;
  using Base::m_elemsToNodes;


  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               localIndex const targetRegionIndex,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE & inputConstitutiveType,
               arrayView1d< globalIndex const > const & inputDispDofNumber,
               string const & inputFlowDofKey,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const & inputMatrix,
               arrayView1d< real64 > const & inputRhs,
               real64 const (&inputGravityVector)[3],
               arrayView1d< string const > const fluidModelNames ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          inputDispDofNumber,
          inputFlowDofKey,
          rankOffset,
          inputMatrix,
          inputRhs,
          inputGravityVector,
          fluidModelNames ),
    m_temperature( nodeManager.template getReference< array1d< real64 > >( dataRepository::keys::Temperature ) ),
    m_deltaTemperature( nodeManager.template getReference< array1d< real64 > >( "newDeltaTemperature" ) )
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
    real64 const thermalStressCoefficient = m_constitutiveUpdate.getThermalStressCoefficient();
    real64 const thermalPorosityCoefficient = m_constitutiveUpdate.getThermalPorosityCoefficient();
    Base::quadraturePointKernel( k, q, stack, [=] GEOSX_HOST_DEVICE ( real64 (& stress)[6], real64 & porosity, real64 & porosityOld )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, q );
      real64 const thermalStress = thermalStressCoefficient * m_temperature[localNodeIndex];

      stress[0] -= thermalStress;
      stress[1] -= thermalStress;
      stress[2] -= thermalStress;



      porosityOld -= thermalPorosityCoefficient * ( m_temperature[localNodeIndex] - m_deltaTemperature[localNodeIndex] );


      porosity -= thermalPorosityCoefficient * m_temperature[localNodeIndex];

    } );

  }

protected:
  /// The rank-global temperature array.
  arrayView1d< real64 const > const m_temperature;

  /// The rank-global incremental temperature array.
  arrayView1d< real64 const > const m_deltaTemperature;

};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               arrayView1d< globalIndex const > const &,
                                                               string const &,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const &,
                                                               arrayView1d< real64 > const &,
                                                               real64 const (&)[3],
                                                               arrayView1d< string const > const >;

} /* namespace THMKernels */

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEWEAKTHMCOUPLINGKERNEL_HPP_ */
