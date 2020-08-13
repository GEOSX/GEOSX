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
 * @file SinglePhaseBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

namespace geosx
{
namespace constitutive
{
class ConstitutiveBase;

}

/**
 * @class SinglePhaseBase
 *
 * base class to perform a single phase finite volume solve.
 */
class SinglePhaseBase : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseBase( const std::string & name, Group * const parent );

  /// deleted default constructor
  SinglePhaseBase() = delete;

  /// deleted copy constructor
  SinglePhaseBase( SinglePhaseBase const & ) = delete;

  /// default move constructor
  SinglePhaseBase( SinglePhaseBase && ) = default;

  /// deleted assignment operator
  SinglePhaseBase &
  operator=( SinglePhaseBase const & ) = delete;

  /// deleted move operator
  SinglePhaseBase &
  operator=( SinglePhaseBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseBase() override = default;

  virtual void
  RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual void
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  ///@{

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  AssembleSystem(
    real64 const time_n,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  virtual void
  ApplyBoundaryConditions(
    real64 const time_n,
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
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  template< bool ISPORO, typename POLICY >
  void
  AccumulationLaunch( localIndex const targetIndex,
                      CellElementSubRegion & subRegion,
                      DofManager const & dofManager,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs );

  template< bool ISPORO, typename POLICY >
  void
  AccumulationLaunch( localIndex const targetIndex,
                      FaceElementSubRegion const & subRegion,
                      DofManager const & dofManager,
                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                      arrayView1d< real64 > const & localRhs );

  ///@}

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  template< bool ISPORO, typename POLICY >
  void
  AssembleAccumulationTerms(
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  AssembleFluxTerms(
    real64 const time_n,
    real64 const dt,
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) = 0;

  void
  ApplyDiricletBC( real64 const time_n,
                   real64 const dt,
                   DomainPartition & domain,
                   DofManager const & dofManager,
                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                   arrayView1d< real64 > const & localRhs ) const;

  void
  ApplySourceFluxBC( real64 const time_n,
                     real64 const dt,
                     DomainPartition & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Function to update all constitutive state and dependent variables
   * @param dataGroup group that contains the fields
   */
  virtual void
  UpdateState( Group & dataGroup, localIndex const targetIndex ) const;

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    // used for face-based BC
    static constexpr auto facePressureString = "facePressure";

    // intermediate fields
    static constexpr auto mobilityString = "mobility";
    static constexpr auto dMobility_dPressureString = "dMobility_dPressure";
    static constexpr auto porosityString = "porosity";

    // backup fields
    static constexpr auto porosityOldString = "porosityOld";
    static constexpr auto densityOldString = "densityOld";
    static constexpr auto transTMultString = "transTMult";
    static constexpr auto poroMultString = "poroMult";

  } viewKeysSinglePhaseBase;

  viewKeyStruct &
  viewKeys()
  {
    return viewKeysSinglePhaseBase;
  }

  viewKeyStruct const &
  viewKeys() const
  {
    return viewKeysSinglePhaseBase;
  }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseBase;

  groupKeyStruct &
  groupKeys()
  {
    return groupKeysSinglePhaseBase;
  }

  groupKeyStruct const &
  groupKeys() const
  {
    return groupKeysSinglePhaseBase;
  }

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void
  ResetViews( MeshLevel & mesh ) override;

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    dataRepository::Group * const rootGroup ) override;

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param mesh the mesh to operate on
   */
  void
  BackupFields( MeshLevel & mesh ) const;

protected:
  virtual void
  ValidateFluidModels( DomainPartition const & domain ) const;

  /**
   * @brief Structure holding views into fluid properties used by the base solver.
   */
  struct FluidPropViews
  {
    arrayView2d< real64 const > const & dens;         ///< density
    arrayView2d< real64 const > const & dDens_dPres;  ///< derivative of density w.r.t. pressure
    arrayView2d< real64 const > const & visc;         ///< viscosity
    arrayView2d< real64 const > const & dVisc_dPres;  ///< derivative of viscosity w.r.t. pressure
    real64 const defaultDensity;                      ///< default density to use for new elements
    real64 const defaulViscosity;                     ///< default vi to use for new elements
  };

  /**
   * @brief Extract properties from a fluid.
   * @param fluid base reference to the fluid object
   * @return structure with property views
   *
   * This function enables derived solvers to substitute SingleFluidBase for a different,
   * unrelated fluid class, and customize property extraction. For example, it is used by
   * SinglePhaseProppantBase to allow using SlurryFluidBase, which does not inherit from
   * SingleFluidBase currently (but this design may need to be revised).
   */
  virtual FluidPropViews
  getFluidProperties(
    constitutive::ConstitutiveBase const & fluid ) const;

  /**
   * @brief Extract pore volume multiplier array from a subregion.
   * @param subRegion the subregion reference
   * @return array view for pore volume multiplier
   *
   * This function allows derived solvers to specialize access to pore volume multiplier that
   * is used in accumulation kernel. For example, it is used by SinglePhaseProppantBase to use
   * multiplier produced by ProppantTransport solver, which we otherwise don't know about.
   * This design should DEFINITELY be revisited.
   */
  virtual arrayView1d< real64 const > const &
  getPoreVolumeMult(
    ElementSubRegionBase const & subRegion ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  virtual void
  UpdateFluidModel( Group & dataGroup,
                    localIndex const targetIndex ) const;

  /**
   * @brief Function to update all constitutive models
   * @param dataGroup group that contains the fields
   */
  void
  UpdateSolidModel( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Function to update fluid mobility
   * @param dataGroup group that contains the fields
   */
  void
  UpdateMobility( Group & dataGroup, localIndex const targetIndex ) const;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaVolume;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_mobility;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_dMobility_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_density;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dDens_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_viscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dVisc_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > m_transTMultiplier;

private:
  virtual void
  ResetViewsPrivate( ElementRegionManager const & elemManager );
};

} /* namespace geosx */

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEBASE_HPP_
