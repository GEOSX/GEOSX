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
 * @file SinglePhasePoroelasticKernel.hpp
 */


#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{

namespace PoromechanicsKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase THM problems with weak thermal coupling.
 * @copydoc geosx::finiteElement::ImplicitKernelBase
 * @tparam NUM_NODES_PER_ELEM The number of nodes per element for the
 *                            @p SUBREGION_TYPE.
 * @tparam UNUSED An unused parameter since we are assuming that the test and
 *                trial space have the same number of support points.
 *
 * ### SinglePhasePoroelastic Description
 * Implements the KernelBase interface functions required for solving the
 * quasi-static single-phase poromechanics problem using one of the
 * "finite element kernel application" functions such as
 * geosx::finiteElement::RegionBasedKernelApplication.
 *
 */
template< typename SUBREGION_TYPE,
          typename CONSTITUTIVE_TYPE,
          typename FE_TYPE >
class SinglePhase :
  public finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                            CONSTITUTIVE_TYPE,
                                            FE_TYPE,
                                            3,
                                            3 >
{
public:
  /// Alias for the base class;
  using Base = finiteElement::ImplicitKernelBase< SUBREGION_TYPE,
                                                  CONSTITUTIVE_TYPE,
                                                  FE_TYPE,
                                                  3,
                                                  3 >;

  /// Number of nodes per element...which is equal to the
  /// numTestSupportPointPerElem and numTrialSupportPointPerElem by definition.
  static constexpr int numNodesPerElem = Base::numTestSupportPointsPerElem;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  /**
   * @brief Constructor
   * @copydoc geosx::finiteElement::ImplicitKernelBase::ImplicitKernelBase
   * @param inputGravityVector The gravity vector.
   */
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
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_uhat( nodeManager.incrementalDisplacement()),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_fluidDensity( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).density() ),
    m_fluidDensityOld( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::densityOldString() ) ),
    m_dFluidDensity_dPressure( elementSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).dDensity_dPressure() ),
    m_flowDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey )),
    m_fluidPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() ) ),
    m_deltaFluidPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() ) ),
    m_poroRef( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::referencePorosityString() ) ),
    m_biotCoefficient( m_constitutiveUpdate.getBiotCoefficient() )
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

    static constexpr int numDispDofPerElem =  Base::StackVariables::numRows;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            xLocal(),
            u_local(),
            uhat_local(),
            localFlowResidual{ 0.0 },
      localDispFlowJacobian{ {0.0} },
      localFlowDispJacobian{ {0.0} },
      localFlowFlowJacobian{ {0.0} },
      localFlowDofIndex{ 0 }
    {}

#if !defined(CALC_FEM_SHAPE_IN_KERNEL)
    /// Dummy
    int xLocal;
#else
    /// C-array stack storage for element local the nodal positions.
    real64 xLocal[numNodesPerElem][3];
#endif

    /// Stack storage for the element local nodal displacement
    real64 u_local[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for the element local nodal incremental displacement
    real64 uhat_local[numNodesPerElem][numDofPerTrialSupportPoint];

    real64 localFlowResidual[1];
    real64 localDispFlowJacobian[numDispDofPerElem][1];
    real64 localFlowDispJacobian[1][numDispDofPerElem];
    real64 localFlowFlowJacobian[1][1];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex localFlowDofIndex[1];

  };
  //*****************************************************************************

  /**
   * @brief Copy global values from primary field to a local stack array.
   * @copydoc ::geosx::finiteElement::ImplicitKernelBase::setup
   *
   * For the SinglePhase implementation, global values from the displacement,
   * incremental displacement, and degree of freedom numbers are placed into
   * element local stack storage.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void setup( localIndex const k,
              StackVariables & stack ) const
  {
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
#if defined(CALC_FEM_SHAPE_IN_KERNEL)
        stack.xLocal[a][i] = m_X[localNodeIndex][i];
#endif
        stack.u_local[a][i] = m_disp[localNodeIndex][i];
        stack.uhat_local[a][i] = m_uhat[localNodeIndex][i];
        stack.localRowDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
        stack.localColDofIndex[a*3+i] = m_dofNumber[localNodeIndex]+i;
      }
    }

    for( int flowDofIndex=0; flowDofIndex<1; ++flowDofIndex )
    {
      stack.localFlowDofIndex[flowDofIndex] = m_flowDofNumber[k] + flowDofIndex;
    }

  }

  /**
   * @brief Internal struct to provide no-op defaults used in the inclusion
   *   of lambda functions into kernel component functions.
   * @struct NoOpFunctors
   */
  struct NoOpFunctors
  {
    /**
     * @brief operator() no-op used for modifying the stress tensor ( prior to
     *   integrating the divergence to produce nodal forces) and the porosity.
     * @param stress The stress array.
     * @param porosity The porosity.
     */
    GEOSX_HOST_DEVICE
    GEOSX_FORCE_INLINE constexpr
    void operator() ( real64 (& stress)[6], real64 & porosity, real64 & porosityOld )
    {
      GEOSX_UNUSED_VAR( stress );
      GEOSX_UNUSED_VAR( porosity );
      GEOSX_UNUSED_VAR( porosityOld );
    }
  };



  template< typename MODIFIER = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack,
                              MODIFIER && modifier = NoOpFunctors{} ) const
  {
    // Get displacement: (i) basis functions (N), (ii) basis function
    // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
    // matrix times the quadrature weight (detJxW)
    real64 N[numNodesPerElem];
    real64 dNdX[numNodesPerElem][3];
    FE_TYPE::calcN( q, N );
    real64 const detJxW = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // Evaluate total stress tensor
    real64 strainIncrement[6] = {0};
    real64 totalStress[6] = {0};

    // --- Update effective stress tensor (stored in totalStress)
    typename CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps stiffness;
    FE_TYPE::symmetricGradient( dNdX, stack.uhat_local, strainIncrement );
    m_constitutiveUpdate.smallStrainUpdate( k, q, strainIncrement, totalStress, stiffness );

    // --- Subtract pressure term
    real64 const biotTimesPressure = m_biotCoefficient * ( m_fluidPressure[k] + m_deltaFluidPressure[k] );
    totalStress[0] -= biotTimesPressure;
    totalStress[1] -= biotTimesPressure;
    totalStress[2] -= biotTimesPressure;

    // Evaluate fluid phase mass content
    // --- TODO: temporary solution -----------------------------------------------------------------------------------//
    //           update PoroElastic Constitutive Model (PECM) such that:                                               //
    //           (  i) the kernel has a view m_porosityOld to the field stored in PECM                                 //
    //           ( ii) the newPororosity, stored in PECM, is updated as                                                //
    //                 m_constitutiveUpdate.porosityUpdate( k, q, m_deltaFluidPressure, volStrainIncrement );          //
    //           (iii) get methods for porosityNew, dPorosity_dFluidPressure, dPorosity_dVolStrainIncrement            //
    real64 const biotSkeletonModulusInverse = 0.0; //TODO: 1/N = 0 correct only for biotCoefficient = 1                //
    real64 const volumetricStrainNew = FE_TYPE::symmetricGradientTrace( dNdX, stack.u_local );                         //
    real64 const volumetricStrainOld = volumetricStrainNew - FE_TYPE::symmetricGradientTrace( dNdX, stack.uhat_local );//
    real64 porosityOld = m_poroRef( k ) + m_biotCoefficient * volumetricStrainOld;// +  DeltaPoro                //
    real64 const dPorosity_dPressure = biotSkeletonModulusInverse;                                                     //
    real64 const dPorosity_dVolStrainIncrement =  m_biotCoefficient;                                                   //
                                                                                                                       //
    GEOSX_ERROR_IF_GT_MSG( fabs( m_biotCoefficient - 1.0 ),                                                            //
                           1e-10,                                                                                      //
                           "Correct only for Biot's coefficient equal to 1" );                                         //
    // --------------------------------------------------------------------------------------------------------------- //
    real64 porosityNew = porosityOld
                         + m_biotCoefficient * (strainIncrement[0] + strainIncrement[1] + strainIncrement[2] )
                         + biotSkeletonModulusInverse * m_deltaFluidPressure[k];

    // Modify stress and porosity in case of thermal coupling
    modifier( totalStress, porosityNew, porosityOld );

    // Evaluate body force vector
    real64 bodyForce[3] = { m_gravityVector[0],
                            m_gravityVector[1],
                            m_gravityVector[2]};
    if( m_gravityAcceleration > 0.0 )
    {
      real64 mixtureDensity = ( 1.0 - porosityNew ) * m_solidDensity( k, q ) + porosityNew * m_fluidDensity( k, q );
      mixtureDensity *= detJxW;
      bodyForce[0] *= mixtureDensity;
      bodyForce[1] *= mixtureDensity;
      bodyForce[2] *= mixtureDensity;
    }

    // Assemble local jacobian and residual

    // --- Momentum balance
    for( localIndex i=0; i<6; ++i )
    {
      totalStress[i] *= -detJxW;
    }

    FE_TYPE::plusGradNajAijPlusNaFi( dNdX,
                                     totalStress,
                                     N,
                                     bodyForce,
                                     reinterpret_cast< real64 (&)[numNodesPerElem][3] >(stack.localResidual) );

    stiffness.template upperBTDB< numNodesPerElem >( dNdX, -detJxW, stack.localJacobian );

    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localDispFlowJacobian[a*3+0][0] += dNdX[a][0] * m_biotCoefficient * detJxW;
      stack.localDispFlowJacobian[a*3+1][0] += dNdX[a][1] * m_biotCoefficient * detJxW;
      stack.localDispFlowJacobian[a*3+2][0] += dNdX[a][2] * m_biotCoefficient * detJxW;
    }

    if( m_gravityAcceleration > 0.0 )
    {
      // Assumptions: ( i) dMixtureDens_dVolStrain contribution is neglected
      //              (ii) grains are assumed incompressible

      real64 const dMixtureDens_dPressure = dPorosity_dPressure * ( -m_solidDensity( k, q ) + m_fluidDensity( k, q ) )
                                            + porosityNew * m_dFluidDensity_dPressure( k, q );
      for( integer a = 0; a < numNodesPerElem; ++a )
      {
        stack.localDispFlowJacobian[a*3+0][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[0] * detJxW;
        stack.localDispFlowJacobian[a*3+1][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[1] * detJxW;
        stack.localDispFlowJacobian[a*3+2][0] += N[a] * dMixtureDens_dPressure * m_gravityVector[2] * detJxW;
      }
    }

    // --- Mass balance accumulation
    for( integer a = 0; a < numNodesPerElem; ++a )
    {
      stack.localFlowDispJacobian[0][a*3+0] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][0] * detJxW;
      stack.localFlowDispJacobian[0][a*3+1] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][1] * detJxW;
      stack.localFlowDispJacobian[0][a*3+2] += dPorosity_dVolStrainIncrement * m_fluidDensity( k, q ) * dNdX[a][2] * detJxW;
    }

    stack.localFlowResidual[0] += ( porosityNew * m_fluidDensity( k, q ) - porosityOld * m_fluidDensityOld( k ) ) * detJxW;
    stack.localFlowFlowJacobian[0][0] += ( dPorosity_dPressure * m_fluidDensity( k, q ) + porosityNew * m_dFluidDensity_dPressure( k, q ) ) * detJxW;

  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  template< typename MODIFIER = NoOpFunctors >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    GEOSX_UNUSED_VAR( k );

    real64 maxForce = 0;

    CONSTITUTIVE_TYPE::KernelWrapper::DiscretizationOps::template fillLowerBTDB< numNodesPerElem >( stack.localJacobian );

    constexpr int nUDof = numNodesPerElem * numDofPerTestSupportPoint;

    for( int localNode = 0; localNode < numNodesPerElem; ++localNode )
    {
      for( int dim = 0; dim < numDofPerTestSupportPoint; ++dim )
      {
        localIndex const dof = LvArray::integerConversion< localIndex >( stack.localRowDofIndex[numDofPerTestSupportPoint*localNode + dim] - m_dofRankOffset );
        if( dof < 0 || dof >= m_matrix.numRows() ) continue;
        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localRowDofIndex,
                                                                                stack.localJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                numNodesPerElem * numDofPerTrialSupportPoint );

        RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localResidual[numDofPerTestSupportPoint * localNode + dim] );
        maxForce = fmax( maxForce, fabs( stack.localResidual[numDofPerTestSupportPoint * localNode + dim] ) );

        m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                                stack.localFlowDofIndex,
                                                                                stack.localDispFlowJacobian[numDofPerTestSupportPoint * localNode + dim],
                                                                                1 );

      }
    }


    localIndex const dof = LvArray::integerConversion< localIndex >( stack.localFlowDofIndex[0] - m_dofRankOffset );
    if( 0 <= dof && dof < m_matrix.numRows() )
    {
      m_matrix.template addToRowBinarySearchUnsorted< serialAtomic >( dof,
                                                                      stack.localRowDofIndex,
                                                                      stack.localFlowDispJacobian[0],
                                                                      nUDof );
      m_matrix.template addToRow< serialAtomic >( dof,
                                                  stack.localFlowDofIndex,
                                                  stack.localFlowFlowJacobian[0],
                                                  1 );
      RAJA::atomicAdd< serialAtomic >( &m_rhs[dof], stack.localFlowResidual[0] );
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_uhat;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView1d< real64 const > const m_fluidDensityOld;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_flowDofNumber;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_fluidPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaFluidPressure;

  /// The rank-global reference porosity array
  arrayView1d< real64 const > const m_poroRef;

  /// Biot's coefficient
  real64 const m_biotCoefficient;


};

using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               arrayView1d< globalIndex const > const &,
                                                               string const &,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const &,
                                                               arrayView1d< real64 > const &,
                                                               real64 const (&)[3],
                                                               arrayView1d< string const > const >;

} // namespace PoroelasticKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROMECHANICSKERNEL_HPP_
