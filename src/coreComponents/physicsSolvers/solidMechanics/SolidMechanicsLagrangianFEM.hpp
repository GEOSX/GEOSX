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
 * @file SolidMechanicsLagrangianFEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_
#define GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_

#include "common/TimingMacros.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/SolverBase.hpp"

#include "SolidMechanicsLagrangianFEMKernels.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;
}
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

/**
 * @class SolidMechanicsLagrangianFEM
 *
 * This class implements a finite element solution to the equations of motion.
 */
class SolidMechanicsLagrangianFEM : public SolverBase
{
public:
  /**
   * Constructor
   * @param name The name of the solver instance
   * @param parent the parent group of the solver
   */
  SolidMechanicsLagrangianFEM( const std::string & name, Group * const parent );

  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM( SolidMechanicsLagrangianFEM && ) = default;

  SolidMechanicsLagrangianFEM &
  operator=( SolidMechanicsLagrangianFEM const & ) = delete;
  SolidMechanicsLagrangianFEM &
  operator=( SolidMechanicsLagrangianFEM && ) = delete;

  /**
   * destructor
   */
  virtual ~SolidMechanicsLagrangianFEM() override;

  /**
   * @return The string that may be used to generate a new instance from the SolverBase::CatalogInterface::CatalogType
   */
  static string
  CatalogName()
  {
    return "SolidMechanics_LagrangianFEM";
  }

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  RegisterDataOnMesh( Group * const MeshBody ) override final;

  void
  updateIntrinsicNodalData( DomainPartition * const domain );

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  virtual real64
  ExplicitStep( real64 const & time_n,
                real64 const & dt,
                integer const cycleNumber,
                DomainPartition & domain ) override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = false ) override;

  virtual void
  AssembleSystem(
    real64 const time,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  ApplyBoundaryConditions(
    real64 const time,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm(
    DomainPartition const & domain,
    DofManager const & dofManager,
    arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  void
  ResetStressToBeginningOfStep( DomainPartition & domain );

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**@}*/

  template< typename CONSTITUTIVE_BASE,
            template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
            class KERNEL_TEMPLATE,
            typename... PARAMS >
  void
  AssemblyLaunch( DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs,
                  PARAMS &&... params );

  template< typename... PARAMS >
  real64
  explicitKernelDispatch( PARAMS &&... params );

  /**
   * Applies displacement boundary conditions to the system for implicit time integration
   * @param time The time to use for any lookups associated with this BC
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain The DomainPartition.
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   * @param solution the solution vector
   */
  void
  ApplyDisplacementBC_implicit(
    real64 const time,
    DofManager const & dofManager,
    DomainPartition & domain,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs );

  void
  CRSApplyTractionBC( real64 const time,
                      DofManager const & dofManager,
                      DomainPartition & domain,
                      arrayView1d< real64 > const & localRhs );

  void
  ApplyChomboPressure( DofManager const & dofManager,
                       DomainPartition & domain,
                       arrayView1d< real64 > const & localRhs );

  void
  ApplyContactConstraint(
    DofManager const & dofManager,
    DomainPartition & domain,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs );

  virtual real64
  ScalingForSystemSolution(
    DomainPartition const & domain,
    DofManager const & dofManager,
    arrayView1d< real64 const > const & localSolution ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr auto vTildeString = "velocityTilde";
    static constexpr auto uhatTildeString = "uhatTilde";
    static constexpr auto cflFactorString = "cflFactor";
    static constexpr auto newmarkGammaString = "newmarkGamma";
    static constexpr auto newmarkBetaString = "newmarkBeta";
    static constexpr auto massDampingString = "massDamping";
    static constexpr auto stiffnessDampingString = "stiffnessDamping";
    static constexpr auto useVelocityEstimateForQSString = "useVelocityForQS";
    static constexpr auto timeIntegrationOptionString = "timeIntegrationOption";
    static constexpr auto maxNumResolvesString = "maxNumResolves";
    static constexpr auto strainTheoryString = "strainTheory";
    static constexpr auto solidMaterialNamesString = "solidMaterialNames";
    static constexpr auto stress_n = "beginningOfStepStress";
    static constexpr auto forceExternal = "externalForce";
    static constexpr auto contactRelationNameString = "contactRelationName";
    static constexpr auto noContactRelationNameString = "NOCONTACT";
    static constexpr auto contactForceString = "contactForce";
    static constexpr auto maxForce = "maxForce";
    static constexpr auto elemsAttachedToSendOrReceiveNodes =
      "elemsAttachedToSendOrReceiveNodes";
    static constexpr auto elemsNotAttachedToSendOrReceiveNodes =
      "elemsNotAttachedToSendOrReceiveNodes";
    static constexpr auto effectiveStress = "effectiveStress";

    dataRepository::ViewKey vTilde = { vTildeString };
    dataRepository::ViewKey uhatTilde = { uhatTildeString };
    dataRepository::ViewKey newmarkGamma = { newmarkGammaString };
    dataRepository::ViewKey newmarkBeta = { newmarkBetaString };
    dataRepository::ViewKey massDamping = { massDampingString };
    dataRepository::ViewKey stiffnessDamping = { stiffnessDampingString };
    dataRepository::ViewKey useVelocityEstimateForQS = {
      useVelocityEstimateForQSString };
    dataRepository::ViewKey timeIntegrationOption = { timeIntegrationOptionString };
  } solidMechanicsViewKeys;

  struct groupKeyStruct
  {
  } solidMechanicsGroupKeys;

  arrayView1d< string const > const &
  solidMaterialNames() const
  {
    return m_solidMaterialNames;
  }

  SortedArray< localIndex > &
  getElemsAttachedToSendOrReceiveNodes(
    ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >(
      viewKeyStruct::elemsAttachedToSendOrReceiveNodes );
  }

  SortedArray< localIndex > &
  getElemsNotAttachedToSendOrReceiveNodes(
    ElementSubRegionBase & subRegion )
  {
    return subRegion.getReference< SortedArray< localIndex > >(
      viewKeyStruct::elemsNotAttachedToSendOrReceiveNodes );
  }

  void
  setEffectiveStress( integer const input )
  {
    m_effectiveStress = input;
  }

  real64 &
  getMaxForce()
  {
    return m_maxForce;
  }

protected:
  virtual void
  PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    dataRepository::Group * const problemManager ) override final;

  real64 m_newmarkGamma;
  real64 m_newmarkBeta;
  real64 m_massDamping;
  real64 m_stiffnessDamping;
  TimeIntegrationOption m_timeIntegrationOption;
  integer m_useVelocityEstimateForQS;
  real64 m_maxForce = 0.0;
  integer m_maxNumResolves;
  integer m_strainTheory;
  array1d< string > m_solidMaterialNames;
  string m_contactRelationName;
  SortedArray< localIndex > m_sendOrReceiveNodes;
  SortedArray< localIndex > m_nonSendOrReceiveNodes;
  MPI_iCommData m_iComm;

  /// Indicates whether or not to use effective stress when integrating the
  /// stress divergence in the kernels. This means calling the poroelastic
  /// variant of the solid mechanics kernels.
  integer m_effectiveStress;

  SolidMechanicsLagrangianFEM();
};

//**********************************************************************************************************************
//**********************************************************************************************************************
//**********************************************************************************************************************

template< typename CONSTITUTIVE_BASE,
          template< typename SUBREGION_TYPE, typename CONSTITUTIVE_TYPE, typename FE_TYPE >
          class KERNEL_TEMPLATE,
          typename... PARAMS >
void
SolidMechanicsLagrangianFEM::AssemblyLaunch(
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs,
  PARAMS &&... params )
{
  GEOSX_MARK_FUNCTION;
  MeshLevel & mesh =
    *( domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 ) );

  NodeManager const & nodeManager = *( mesh.getNodeManager() );

  string const dofKey =
    dofManager.getKey( dataRepository::keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & dofNumber =
    nodeManager.getReference< globalIndex_array >( dofKey );

  ResetStressToBeginningOfStep( domain );

  real64 const gravityVectorData[3] = { gravityVector().Data()[0],
                                        gravityVector().Data()[1],
                                        gravityVector().Data()[2] };

  m_maxForce =
    finiteElement::regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                                 CONSTITUTIVE_BASE,
                                                 CellElementSubRegion,
                                                 KERNEL_TEMPLATE >(
      mesh,
      targetRegionNames(),
      this->getDiscretizationName(),
      m_solidMaterialNames,
      dofNumber,
      dofManager.rankOffset(),
      localMatrix,
      localRhs,
      gravityVectorData,
      std::forward< PARAMS >( params )... );

  ApplyContactConstraint( dofManager, domain, localMatrix, localRhs );
}

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSLAGRANGIANFEM_HPP_ */
