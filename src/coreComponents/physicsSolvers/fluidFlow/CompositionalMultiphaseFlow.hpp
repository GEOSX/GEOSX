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
 * @file CompositionalMultiphaseFlow.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;

namespace keys
{
string const compositionalMultiphaseFlow = "CompositionalMultiphaseFlow";
}
}  // namespace dataRepository
class FieldSpecificationBase;
class FiniteElementBase;
class DomainPartition;

namespace constitutive
{
class MultiFluidBase;
}

//START_SPHINX_INCLUDE_00
/**
 * @class CompositionalMultiphaseFlow
 *
 * A compositional multiphase solver
 */
class CompositionalMultiphaseFlow : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseFlow( const string & name, Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseFlow() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseFlow( CompositionalMultiphaseFlow && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseFlow &
  operator=( CompositionalMultiphaseFlow const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseFlow &
  operator=( CompositionalMultiphaseFlow && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseFlow() override = default;

  //START_SPHINX_INCLUDE_01

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string
  CatalogName()
  {
    return dataRepository::keys::compositionalMultiphaseFlow;
  }

  virtual void
  RegisterDataOnMesh( Group * const MeshBodies ) override;

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

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

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

  virtual real64
  CalculateResidualNorm(
    DomainPartition const & domain,
    DofManager const & dofManager,
    arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  SolveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual real64
  ScalingForSystemSolution(
    DomainPartition const & domain,
    DofManager const & dofManager,
    arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  CheckSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param dataGroup the group storing the required fields
   */
  void
  UpdateComponentFraction( Group & dataGroup ) const;

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param dataGroup the group storing the required fields
   */
  void
  UpdatePhaseVolumeFraction( Group & dataGroup,
                             localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param dataGroup the group storing the required fields
   */
  void
  UpdateFluidModel( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant solid models using current values of pressure
   * @param dataGroup the group storing the required fields
   */
  void
  UpdateSolidModel( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param castedRelPerm the group storing the required fields
   */
  void
  UpdateRelPermModel( Group & castedRelPerm,
                      localIndex const targetIndex ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param castedCapPres the group storing the required fields
   */
  void
  UpdateCapPressureModel( Group & castedCapPres,
                          localIndex const targetIndex ) const;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  void
  UpdatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  void
  UpdateState( Group & dataGroup, localIndex const targetIndex ) const;

  /**
   * @brief Get the number of fluid components (species)
   * @return the number of components
   */
  localIndex
  numFluidComponents() const
  {
    return m_numComponents;
  }

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  localIndex
  numFluidPhases() const
  {
    return m_numPhases;
  }

  /**
   * @brief assembles the accumulation terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void
  AssembleAccumulationTerms(
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void
  AssembleFluxTerms( real64 const dt,
                     DomainPartition const & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief assembles the volume balance terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void
  AssembleVolumeBalanceTerms(
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) const;

  /**@}*/

  arrayView1d< string const > const &
  relPermModelNames() const
  {
    return m_relPermModelNames;
  }

  arrayView1d< string const > const &
  capPresModelNames() const
  {
    return m_capPressureModelNames;
  }

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "compositionalVariables";

    // inputs
    static constexpr auto temperatureString = "temperature";
    static constexpr auto useMassFlagString = "useMass";

    static constexpr auto relPermNamesString = "relPermNames";
    static constexpr auto capPressureNamesString = "capPressureNames";

    static constexpr auto maxCompFracChangeString = "maxCompFractionChange";
    static constexpr auto allowLocalCompDensChoppingString =
      "allowLocalCompDensityChopping";

    static constexpr auto facePressureString = "facePressure";
    static constexpr auto bcPressureString = "bcPressure";

    static constexpr auto globalCompDensityString = "globalCompDensity";
    static constexpr auto deltaGlobalCompDensityString =
      "deltaGlobalCompDensity";

    // intermediate values for constitutive model input
    static constexpr auto globalCompFractionString = "globalCompFraction";
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString =
      "dGlobalCompFraction_dGlobalCompDensity";

    // intermediate values for saturations
    static constexpr auto phaseVolumeFractionString = "phaseVolumeFraction";
    static constexpr auto dPhaseVolumeFraction_dPressureString =
      "dPhaseVolumeFraction_dPressure";
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString =
      "dPhaseVolumeFraction_dGlobalCompDensity";

    // intermediate values for mobilities
    static constexpr auto phaseMobilityString = "phaseMobility";
    static constexpr auto dPhaseMobility_dPressureString =
      "dPhaseMobility_dPressure";
    static constexpr auto dPhaseMobility_dGlobalCompDensityString =
      "dPhaseMobility_dGlobalCompDensity";

    // these are used to store last converged time step values
    static constexpr auto phaseVolumeFractionOldString =
      "phaseVolumeFractionOld";
    static constexpr auto phaseDensityOldString = "phaseDensityOld";
    static constexpr auto phaseComponentFractionOldString =
      "phaseComponentFractionOld";
    static constexpr auto porosityOldString = "porosityOld";

    // these are allocated on faces for BC application until we can get constitutive models on faces
    static constexpr auto phaseViscosityString = "phaseViscosity";
    static constexpr auto phaseRelativePermeabilityString =
      "phaseRelativePermeability";
    static constexpr auto phaseCapillaryPressureString =
      "phaseCapillaryPressure";

  } viewKeysCompMultiphaseFlow;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseFlow;

  /**
   * @brief Initialize all variables from initial conditions
   * @param domain the domain containing the mesh and fields
   *
   * Initialize all variables from initial conditions. This calculating primary variable values
   * from prescribed intermediate values (i.e. global densities from global fractions)
   * and any applicable hydrostatic equilibration of the domain
   */
  void
  InitializeFluidState( MeshLevel & mesh ) const;

  /**
   * @brief Backup current values of all constitutive fields that participate in the accumulation term
   * @param domain the domain containing the mesh and fields
   */
  void
  BackupFields( MeshLevel & mesh ) const;

  /**
   * @brief Function to perform the Application of Dirichlet type BC's
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  ApplyDirichletBC( real64 const time,
                    real64 const dt,
                    DofManager const & dofManager,
                    DomainPartition & domain,
                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                    arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Apply source flux boundary conditions to the system
   * @param time current time
   * @param dt time step
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param localMatrix local system matrix
   * @param localRhs local system right-hand side vector
   */
  void
  ApplySourceFluxBC( real64 const time,
                     real64 const dt,
                     DofManager const & dofManager,
                     DomainPartition & domain,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const;

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void
  ChopNegativeDensities( DomainPartition & domain );

protected:
  virtual void
  PostProcessInput() override;

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    dataRepository::Group * const rootGroup ) override;

  /**
   * @brief Checks constitutive models for consistency
   * @param cm        reference to the global constitutive model manager
   */
  void
  ValidateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const;

private:
  /**
   * @brief Resize the allocated multidimensional fields
   * @param domain the domain containing the mesh and fields
   *
   * Resize fields along dimensions 1 and 2 (0 is the size of containing object, i.e. element subregion)
   * once the number of phases/components is known (e.g. component fractions)
   */
  void
  ResizeFields( MeshLevel & meshLevel ) const;

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void
  ResetViews( MeshLevel & mesh ) override;

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// name of the rel perm constitutive model
  array1d< string > m_relPermModelNames;

  /// flag to determine whether or not to apply capillary pressure
  integer m_capPressureFlag;

  /// name of the cap pressure constitutive model
  array1d< string > m_capPressureModelNames;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaPressure;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
    m_dPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_phaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dPhaseMob_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_phaseDens;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dPhaseDens_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > m_dPhaseDens_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > m_phaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > >
    m_dPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 const > >
    m_dPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_phaseCapPressure;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > >
    m_dPhaseCapPressure_dPhaseVolFrac;
};

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_COMPOSITIONALMULTIPHASEFLOW_HPP_
