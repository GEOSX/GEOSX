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

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICEFEMKERNEL_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICEFEMKERNEL_HPP_
#include "finiteElement/kernelInterface/ImplicitKernelBase.hpp"
#include "SinglePhasePoroelasticKernel.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"


namespace geosx
{

namespace PoroelasticEFEMKernels
{

/**
 * @brief Implements kernels for solving quasi-static single-phase poromechanics.
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
  /// Compile time value for the number of quadrature points per element.
  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  using Base::numDofPerTestSupportPoint;
  using Base::numDofPerTrialSupportPoint;
  using Base::m_dofNumber;
  using Base::m_dofRankOffset;
  using Base::m_matrix;
  using Base::m_rhs;
  using Base::m_elemsToNodes;
  using Base::m_constitutiveUpdate;
  using Base::m_finiteElementSpace;


  SinglePhase( NodeManager const & nodeManager,
               EdgeManager const & edgeManager,
               FaceManager const & faceManager,
               localIndex const targetRegionIndex,
               SUBREGION_TYPE const & elementSubRegion,
               FE_TYPE const & finiteElementSpace,
               CONSTITUTIVE_TYPE & inputConstitutiveType,
               EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion,
               arrayView1d< globalIndex const > const & dispDofNumber,
               arrayView1d< globalIndex const > const & jumpDofNumber,
               string const & inputFlowDofKey,
               globalIndex const rankOffset,
               CRSMatrixView< real64, globalIndex const > const & inputMatrix,
               arrayView1d< real64 > const & inputRhs,
               real64 const (&inputGravityVector)[3],
               arrayView1d< string const > const & fluidModelNames ):
    Base( nodeManager,
          edgeManager,
          faceManager,
          targetRegionIndex,
          elementSubRegion,
          finiteElementSpace,
          inputConstitutiveType,
          dispDofNumber,
          rankOffset,
          inputMatrix,
          inputRhs ),
    m_X( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement()),
    m_deltaDisp( nodeManager.incrementalDisplacement()),
    m_w( embeddedSurfSubRegion.displacementJump() ),
    m_matrixPresDofNumber( elementSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
    m_fracturePresDofNumber( embeddedSurfSubRegion.template getReference< array1d< globalIndex > >( inputFlowDofKey ) ),
    m_wDofNumber( jumpDofNumber ),
    m_solidDensity( inputConstitutiveType.getDensity() ),
    m_fluidDensity( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).density() ),
    m_fluidDensityOld( embeddedSurfSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::densityOldString() ) ),
    m_dFluidDensity_dPressure( embeddedSurfSubRegion.template getConstitutiveModel< constitutive::SingleFluidBase >( fluidModelNames[targetRegionIndex] ).dDensity_dPressure() ),
    m_matrixPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::pressureString() ) ),
    m_deltaMatrixPressure( elementSubRegion.template getReference< array1d< real64 > >( SinglePhaseBase::viewKeyStruct::deltaPressureString() ) ),
    m_oldPorosity( inputConstitutiveType.getOldPorosity() ),
    m_tractionVec( embeddedSurfSubRegion.tractionVector() ),
    m_dTraction_dJump( embeddedSurfSubRegion.dTraction_dJump() ),
    m_dTraction_dPressure( embeddedSurfSubRegion.dTraction_dPressure() ),
    m_nVec( embeddedSurfSubRegion.getNormalVector() ),
    m_tVec1( embeddedSurfSubRegion.getTangentVector1() ),
    m_tVec2( embeddedSurfSubRegion.getTangentVector2() ),
    m_surfaceCenter( embeddedSurfSubRegion.getElementCenter() ),
    m_surfaceArea( embeddedSurfSubRegion.getElementArea() ),
    m_elementVolume( elementSubRegion.getElementVolume() ),
    m_deltaVolume( elementSubRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString() ) ),
    m_fracturedElems( elementSubRegion.fracturedElementsList() ),
    m_cellsToEmbeddedSurfaces( elementSubRegion.embeddedSurfacesList().toViewConst() ),
    m_gravityVector{ inputGravityVector[0], inputGravityVector[1], inputGravityVector[2] },
    m_gravityAcceleration( LvArray::tensorOps::l2Norm< 3 >( inputGravityVector ) )
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

    /// The number of displacement dofs per element.
    static constexpr int numUdofs = numNodesPerElem * 3;


    /// The number of jump dofs per element.
    static constexpr int numWdofs = 3;

    /// Constructor.
    GEOSX_HOST_DEVICE
    StackVariables():
      Base::StackVariables(),
            dispEqnRowIndices{ 0 },
      dispColIndices{ 0 },
      jumpEqnRowIndices{ 0 },
      jumpColIndices{ 0 },
      localDispResidual{ 0.0 },
      localJumpResidual{ 0.0 },
      localKww{ { 0.0 } },
      localKwu{ { 0.0 } },
      localKuw{ { 0.0 } },
      localKwpm{ 0.0 },
      localKwpf( 0.0 ),
      wLocal(),
      dispLocal(),
      deltaDispLocal(),
      hInv(),
      xLocal(),
      tractionVec(),
      dTractiondw{ { 0.0 } },
      constitutiveStiffness()
    {}

    /// C-array storage for the element local row degrees of freedom.
    globalIndex dispEqnRowIndices[numUdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex dispColIndices[numUdofs];

    /// C-array storage for the element local row degrees of freedom.
    globalIndex jumpEqnRowIndices[numWdofs];

    /// C-array storage for the element local column degrees of freedom.
    globalIndex jumpColIndices[numWdofs];

    /// C-array storage for the element local Ru residual vector.
    real64 localDispResidual[numUdofs];

    /// C-array storage for the element local Rw residual vector.
    real64 localJumpResidual[numWdofs];

    /// C-array storage for the element local Kww matrix.
    real64 localKww[numWdofs][numWdofs];

    /// C-array storage for the element local Kwu matrix.
    real64 localKwu[numWdofs][numUdofs];

    /// C-array storage for the element local Kuw matrix.
    real64 localKuw[numUdofs][numWdofs];

    /// C-array storage for the element local Kwpm matrix.
    real64 localKwpm[numWdofs];

    /// C-array storage for the element local Kwpf matrix.
    real64 localKwpf;

    /// Stack storage for the element local jump vector
    real64 wLocal[3];

    /// Stack storage for the element displacement vector.
    real64 dispLocal[numUdofs];

    // Stack storage for incremental displacement
    real64 deltaDispLocal[numNodesPerElem][numDofPerTrialSupportPoint];

    /// Stack storage for Area/Volume
    real64 hInv;

    /// local nodal coordinates
    real64 xLocal[ numNodesPerElem ][ 3 ];

    /// Stack storage for the traction
    real64 tractionVec[3];

    /// Stack storage for the derivative of the traction
    real64 dTractiondw[3][3];

    /// Stack storage for the constitutive stiffness at a quadrature point.
    real64 constitutiveStiffness[ 6 ][ 6 ];
  };
  //*****************************************************************************

  //START_kernelLauncher
  template< typename POLICY,
            typename KERNEL_TYPE >
  static
  real64
  kernelLaunch( localIndex const numElems,
                KERNEL_TYPE const & kernelComponent )
  {
    GEOSX_MARK_FUNCTION;

    GEOSX_UNUSED_VAR( numElems );

    // Define a RAJA reduction variable to get the maximum residual contribution.
    RAJA::ReduceMax< ReducePolicy< POLICY >, real64 > maxResidual( 0 );

    forAll< POLICY >( kernelComponent.m_fracturedElems.size(),
                      [=] GEOSX_HOST_DEVICE ( localIndex const i )
    {
      localIndex k = kernelComponent.m_fracturedElems[i];
      typename KERNEL_TYPE::StackVariables stack;

      kernelComponent.setup( k, stack );
      for( integer q=0; q<numQuadraturePointsPerElem; ++q )
      {
        kernelComponent.quadraturePointKernel( k, q, stack );
      }
      maxResidual.max( kernelComponent.complete( k, stack ) );
    } );

    return maxResidual.get();
  }
  //END_kernelLauncher


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
    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    stack.hInv = m_surfaceArea[embSurfIndex] / m_elementVolume[k];
    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        stack.dispEqnRowIndices[a*3+i] = m_dofNumber[localNodeIndex]+i-m_dofRankOffset;
        stack.dispColIndices[a*3+i]    = m_dofNumber[localNodeIndex]+i;
        stack.xLocal[ a ][ i ] = m_X[ localNodeIndex ][ i ];
        stack.dispLocal[ a*3 + i ] = m_disp[ localNodeIndex ][ i ];
        stack.deltaDispLocal[ a ][ i ] = m_deltaDisp[ localNodeIndex ][ i ];
      }
    }

    for( int i=0; i<3; ++i )
    {
      // need to grab the index.
      stack.jumpEqnRowIndices[i] = m_wDofNumber[embSurfIndex] + i - m_dofRankOffset;
      stack.jumpColIndices[i]    = m_wDofNumber[embSurfIndex] + i;
      stack.wLocal[ i ] = m_w[ embSurfIndex ][i];
      stack.tractionVec[ i ] = m_tractionVec[ embSurfIndex ][i] * m_surfaceArea[embSurfIndex];
      for( int ii=0; ii < 3; ++ii )
      {
        stack.dTractiondw[ i ][ ii ] = m_dTraction_dJump[embSurfIndex][i][ii] * m_surfaceArea[embSurfIndex];
      }
    }
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void quadraturePointKernel( localIndex const k,
                              localIndex const q,
                              StackVariables & stack ) const
  {

    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    // Get displacement: (i) basis functions (N), (ii) basis function
    // derivatives (dNdX), and (iii) determinant of the Jacobian transformation
    // matrix times the quadrature weight (detJxW)
    real64 dNdX[numNodesPerElem][3];
    real64 const detJ = m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, stack.xLocal, dNdX );

    // EFEM part starts here
    constexpr int nUdof = numNodesPerElem*3;

    // Gauss contribution to Kww, Kwu and Kuw blocks
    real64 Kww_gauss[3][3], Kwu_gauss[3][nUdof], Kuw_gauss[nUdof][3], Kwpm_gauss[3];

    //  Compatibility, equilibrium and strain operators. The compatibility operator is constructed as
    //  a 3 x 6 because it is more convenient for construction purposes (reduces number of local var).
    real64 compMatrix[3][6], strainMatrix[6][nUdof], eqMatrix[3][6];
    real64 matBD[nUdof][6], matED[3][6];
    real64 const biotCoefficient = 1.0;

    int Heaviside[ numNodesPerElem ];

    // TODO: asking for the stiffness here will only work for elastic models.  most other models
    //       need to know the strain increment to compute the current stiffness value.
    m_constitutiveUpdate.getElasticStiffness( k, stack.constitutiveStiffness );

    SolidMechanicsEFEMKernelsHelper::computeHeavisideFunction< numNodesPerElem >( Heaviside,
                                                                                  stack.xLocal,
                                                                                  m_nVec[embSurfIndex],
                                                                                  m_surfaceCenter[embSurfIndex] );


    SolidMechanicsEFEMKernelsHelper::assembleEquilibriumOperator( eqMatrix,
                                                                  m_nVec[embSurfIndex],
                                                                  m_tVec1[embSurfIndex],
                                                                  m_tVec2[embSurfIndex],
                                                                  stack.hInv );

    SolidMechanicsEFEMKernelsHelper::assembleCompatibilityOperator< numNodesPerElem >( compMatrix,
                                                                                       m_nVec[embSurfIndex],
                                                                                       m_tVec1[embSurfIndex],
                                                                                       m_tVec2[embSurfIndex],
                                                                                       Heaviside,
                                                                                       dNdX );

    SolidMechanicsEFEMKernelsHelper::assembleStrainOperator< 6, nUdof, numNodesPerElem >( strainMatrix, dNdX );

    // transp(B)D
    LvArray::tensorOps::Rij_eq_AkiBkj< nUdof, 6, 6 >( matBD, strainMatrix, stack.constitutiveStiffness );
    // ED
    LvArray::tensorOps::Rij_eq_AikBkj< 3, 6, 6 >( matED, eqMatrix, stack.constitutiveStiffness );
    // EDC
    LvArray::tensorOps::Rij_eq_AikBjk< 3, 3, 6 >( Kww_gauss, matED, compMatrix );
    // EDB
    LvArray::tensorOps::Rij_eq_AikBkj< 3, nUdof, 6 >( Kwu_gauss, matED, strainMatrix );
    // transp(B)DB
    LvArray::tensorOps::Rij_eq_AikBjk< nUdof, 3, 6 >( Kuw_gauss, matBD, compMatrix );

    LvArray::tensorOps::fill< 3 >( Kwpm_gauss, 0 );
    for( int i=0; i < 3; ++i )
    {
      Kwpm_gauss[0] += eqMatrix[0][i];
      Kwpm_gauss[1] += eqMatrix[1][i];
      Kwpm_gauss[2] += eqMatrix[2][i];
    }

    // multiply by determinant and add to element matrix
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, Kww_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< 3, nUdof >( stack.localKwu, Kwu_gauss, -detJ );
    LvArray::tensorOps::scaledAdd< nUdof, 3 >( stack.localKuw, Kuw_gauss, -detJ );
    // No neg coz the effective stress is total stress - porePressure
    // and all signs are flipped here.
    LvArray::tensorOps::scaledAdd< 3 >( stack.localKwpm, Kwpm_gauss, detJ*biotCoefficient );
  }

  /**
   * @copydoc geosx::finiteElement::ImplicitKernelBase::complete
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  real64 complete( localIndex const k,
                   StackVariables & stack ) const
  {
    real64 maxForce = 0;
    constexpr int nUdof = numNodesPerElem*3;

    globalIndex matrixPressureColIndex = m_matrixPresDofNumber[k];

    // Compute the local residuals
    LvArray::tensorOps::Ri_add_AijBj< 3, 3 >( stack.localJumpResidual, stack.localKww, stack.wLocal );
    LvArray::tensorOps::Ri_add_AijBj< 3, nUdof >( stack.localJumpResidual, stack.localKwu, stack.dispLocal );
    LvArray::tensorOps::Ri_add_AijBj< nUdof, 3 >( stack.localDispResidual, stack.localKuw, stack.wLocal );

    // add pore pressure contribution
    LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.localKwpm, m_matrixPressure[ k ] + m_deltaMatrixPressure[ k ] );

    localIndex const embSurfIndex = m_cellsToEmbeddedSurfaces[k][0];

    // Add traction contribution tranction
    LvArray::tensorOps::scaledAdd< 3 >( stack.localJumpResidual, stack.tractionVec, -1 );
    LvArray::tensorOps::scaledAdd< 3, 3 >( stack.localKww, stack.dTractiondw, -1 );

    // JumpFractureFlowJacobian
    real64 const localJumpFracPressureJacobian = -m_dTraction_dPressure[embSurfIndex] * m_surfaceArea[embSurfIndex];

    // Mass balance accumulation
    real64 const newVolume = m_elementVolume( embSurfIndex ) + m_deltaVolume( embSurfIndex );
    real64 const newMass =  m_fluidDensity( embSurfIndex, 0 ) * newVolume;
    real64 const oldMass =  m_fluidDensityOld( embSurfIndex ) * m_elementVolume( embSurfIndex );
    real64 const localFlowResidual = ( newMass - oldMass );
    real64 const localFlowJumpJacobian = m_fluidDensity( embSurfIndex, 0 ) * m_surfaceArea[ embSurfIndex ];
    real64 const localFlowFlowJacobian = m_dFluidDensity_dPressure( embSurfIndex, 0 ) * newVolume;

    for( localIndex i = 0; i < nUdof; ++i )
    {
      localIndex const uDof = LvArray::integerConversion< localIndex >( stack.dispEqnRowIndices[ i ] );
      if( uDof < 0 || uDof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[uDof], stack.localDispResidual[i] );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( uDof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKuw[i],
                                                                              3 );

    }

    for( localIndex i=0; i < 3; ++i )
    {
      localIndex const dof = LvArray::integerConversion< localIndex >( stack.jumpEqnRowIndices[ i ] );

      if( dof < 0 || dof >= m_matrix.numRows() ) continue;

      RAJA::atomicAdd< parallelDeviceAtomic >( &m_rhs[dof], stack.localJumpResidual[i] );

      // fill in matrix
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.jumpColIndices,
                                                                              stack.localKww[i],
                                                                              3 );
      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              stack.dispColIndices,
                                                                              stack.localKwu[i],
                                                                              numNodesPerElem*3 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( dof,
                                                                              &matrixPressureColIndex,
                                                                              &stack.localKwpm[i],
                                                                              1 );
    }

//    // it only affects the normal jump

    if( stack.jumpEqnRowIndices[0] >= 0 && stack.jumpEqnRowIndices[0] < m_matrix.numRows() )
    {

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( stack.jumpEqnRowIndices[0],
                                                                              &m_fracturePresDofNumber[ embSurfIndex ],
                                                                              &localJumpFracPressureJacobian,
                                                                              1 );
    }

    localIndex const fracturePressureDof = m_fracturePresDofNumber[ embSurfIndex ] - m_dofRankOffset;
    if( fracturePressureDof >= 0 && fracturePressureDof < m_matrix.numRows() )
    {

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( fracturePressureDof,
                                                                              &stack.jumpColIndices[0],
                                                                              &localFlowJumpJacobian,
                                                                              1 );

      m_matrix.template addToRowBinarySearchUnsorted< parallelDeviceAtomic >( fracturePressureDof,
                                                                              &m_fracturePresDofNumber[ embSurfIndex ],
                                                                              &localFlowFlowJacobian,
                                                                              1 );

      RAJA::atomicAdd< serialAtomic >( &m_rhs[ fracturePressureDof ], localFlowResidual );
    }

    return maxForce;
  }



protected:
  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_X;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;

  /// The rank-global incremental displacement array.
  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const m_deltaDisp;

  arrayView2d< real64 const > const m_w;

  /// The global degree of freedom number
  arrayView1d< globalIndex const > const m_matrixPresDofNumber;

  arrayView1d< globalIndex const > const m_fracturePresDofNumber;

  arrayView1d< globalIndex const > const m_wDofNumber;

  /// The rank global densities
  arrayView2d< real64 const > const m_solidDensity;
  arrayView2d< real64 const > const m_fluidDensity;
  arrayView1d< real64 const > const m_fluidDensityOld;
  arrayView2d< real64 const > const m_dFluidDensity_dPressure;

  /// The rank-global fluid pressure array.
  arrayView1d< real64 const > const m_matrixPressure;

  /// The rank-global delta-fluid pressure array.
  arrayView1d< real64 const > const m_deltaMatrixPressure;

  arrayView2d< real64 const > const m_oldPorosity;

  arrayView2d< real64 const > const m_tractionVec;

  arrayView3d< real64 const > const m_dTraction_dJump;

  arrayView1d< real64 const > const m_dTraction_dPressure;

  arrayView2d< real64 const > const m_nVec;

  arrayView2d< real64 const > const m_tVec1;

  arrayView2d< real64 const > const m_tVec2;

  arrayView2d< real64 const > const m_surfaceCenter;

  arrayView1d< real64 const > const m_surfaceArea;

  arrayView1d< real64 const > const m_elementVolume;

  arrayView1d< real64 const > const m_deltaVolume;

  SortedArrayView< localIndex const > const m_fracturedElems;

  ArrayOfArraysView< localIndex const > const m_cellsToEmbeddedSurfaces;

  /// The gravity vector.
  real64 const m_gravityVector[3];
  real64 const m_gravityAcceleration;

};


using SinglePhaseKernelFactory = finiteElement::KernelFactory< SinglePhase,
                                                               EmbeddedSurfaceSubRegion const &,
                                                               arrayView1d< globalIndex const > const &,
                                                               arrayView1d< globalIndex const > const &,
                                                               string const &,
                                                               globalIndex const,
                                                               CRSMatrixView< real64, globalIndex const > const &,
                                                               arrayView1d< real64 > const &,
                                                               real64 const (&)[3],
                                                               arrayView1d< string const > const >;


} // namespace PoroelasticEFEMKernels

} // namespace geosx

#include "finiteElement/kernelInterface/SparsityKernelBase.hpp"

#endif // GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_SINGLEPHASEPOROELASTICKERNEL_HPP_
