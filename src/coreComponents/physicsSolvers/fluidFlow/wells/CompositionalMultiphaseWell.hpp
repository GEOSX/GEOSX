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
 * @file CompositionalMultiphaseWell.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_

#include "WellSolverBase.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlow.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;
}

namespace constitutive
{
class ConstitutiveManager;
class MultiFluidBase;
}  // namespace constitutive
class WellElementSubRegion;

/**
 * @class CompositionalMultiphaseWell
 *
 * A compositional multiphase well solver
 */
class CompositionalMultiphaseWell : public WellSolverBase
{
public:
  // define the column offset of the derivatives
  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DCOMP = 1;
  };

  // define the row offset of the residual equations
  struct RowOffset
  {
    static constexpr integer CONTROL = 0;
    static constexpr integer MASSBAL = 1;
  };

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseWell( const string & name, Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseWell() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseWell( CompositionalMultiphaseWell && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseWell &
  operator=( CompositionalMultiphaseWell const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseWell &
  operator=( CompositionalMultiphaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "CompositionalMultiphaseWell";
  }

  virtual void
  RegisterDataOnMesh( Group * const meshBodies ) override;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual real64
  CalculateResidualNorm(
    DomainPartition const & domain,
    DofManager const & dofManager,
    arrayView1d< real64 const > const & localRhs ) override;

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

  /**@}*/

  /**
   * @brief Recompute component fractions from primary variables (component densities)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void
  UpdateComponentFraction( WellElementSubRegion & subRegion ) const;

  /**
   * @brief Update all relevant fluid models using current values of pressure and composition
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void
  UpdateFluidModel( WellElementSubRegion & subRegion,
                    localIndex const targetIndex );

  /**
   * @brief Recompute phase volume fractions (saturations) from constitutive and primary variables
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  void
  UpdatePhaseVolumeFraction( WellElementSubRegion & subRegion,
                             localIndex const targetIndex ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param subRegion the well subregion containing all the primary and dependent fields
   */
  virtual void
  UpdateState( WellElementSubRegion & subRegion,
               localIndex const targetIndex ) override;

  virtual string
  WellElementDofName() const override
  {
    return viewKeyStruct::dofFieldString;
  }

  virtual string
  ResElementDofName() const override
  {
    return CompositionalMultiphaseFlow::viewKeyStruct::dofFieldString;
  }

  virtual localIndex
  NumFluidComponents() const override
  {
    return m_numComponents;
  }

  virtual localIndex
  NumFluidPhases() const override
  {
    return m_numPhases;
  }

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleFluxTerms(
    real64 const time_n,
    real64 const dt,
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleVolumeBalanceTerms(
    real64 const time_n,
    real64 const dt,
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  FormPressureRelations(
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief Sets all the negative component densities (if any) to zero.
   * @param domain the physical domain object
   */
  void
  ChopNegativeDensities( DomainPartition & domain );

  arrayView1d< string const > const &
  relPermModelNames() const
  {
    return m_relPermModelNames;
  }

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "compositionalWellVars";

    // inputs
    static constexpr auto temperatureString = "wellTemperature";
    static constexpr auto useMassFlagString =
      CompositionalMultiphaseFlow::viewKeyStruct::useMassFlagString;

    static constexpr auto relPermNamesString =
      CompositionalMultiphaseFlow::viewKeyStruct::relPermNamesString;

    static constexpr auto maxCompFracChangeString =
      CompositionalMultiphaseFlow::viewKeyStruct::maxCompFracChangeString;
    static constexpr auto allowLocalCompDensChoppingString =
      CompositionalMultiphaseFlow::viewKeyStruct::allowLocalCompDensChoppingString;

    // primary solution field
    static constexpr auto pressureString =
      CompositionalMultiphaseFlow::viewKeyStruct::pressureString;
    static constexpr auto deltaPressureString =
      CompositionalMultiphaseFlow::viewKeyStruct::deltaPressureString;
    static constexpr auto globalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::globalCompDensityString;
    static constexpr auto deltaGlobalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::deltaGlobalCompDensityString;
    static constexpr auto mixtureConnRateString =
      "wellElementMixtureConnectionRate";
    static constexpr auto deltaMixtureConnRateString =
      "deltaWellElementMixtureConnectionRate";

    // saturations
    static constexpr auto phaseVolumeFractionString =
      CompositionalMultiphaseFlow::viewKeyStruct::phaseVolumeFractionString;
    static constexpr auto dPhaseVolumeFraction_dPressureString =
      CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dPressureString;
    static constexpr auto dPhaseVolumeFraction_dGlobalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString;

    // global component fractions
    static constexpr auto globalCompFractionString =
      CompositionalMultiphaseFlow::viewKeyStruct::globalCompFractionString;
    static constexpr auto dGlobalCompFraction_dGlobalCompDensityString =
      CompositionalMultiphaseFlow::viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString;

    // perforation rates and derivatives
    static constexpr auto compPerforationRateString = "compPerforationRate";
    static constexpr auto dCompPerforationRate_dPresString =
      "dCompPerforationRate_dPres";
    static constexpr auto dCompPerforationRate_dCompString =
      "dCompPerforationRate_dComp";

  } viewKeysCompMultiphaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysCompMultiphaseWell;

protected:
  virtual void
  PostProcessInput() override;

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    Group * const rootGroup ) override;

  /**
   * @brief Checks constitutive models for consistency
   * @param meshLevel reference to the mesh
   * @param cm        reference to the global constitutive model manager
   *
   * This function will produce an error if one of the well constitutive models
   * (fluid, relperm) is incompatible with the corresponding models in reservoir
   * regions connected to that particular well.
   */
  void
  ValidateConstitutiveModels( MeshLevel const & meshLevel,
                              constitutive::ConstitutiveManager const & cm ) const;

  /**
   * @brief Checks injection streams for validity (compositions sum to one)
   * @param meshLevel reference to the mesh
   */
  void
  ValidateInjectionStreams( MeshLevel const & meshLevel ) const;

private:
  /**
   * @brief Compute all the perforation rates for this well
   * @param well the well with its perforations
   */
  void
  ComputePerforationRates( WellElementSubRegion & subRegion,
                           localIndex const targetIndex );

  /**
   * @brief Setup stored reservoir views into domain data for the current step
   */
  void
  ResetViews( DomainPartition & domain ) override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void
  InitializeWells( DomainPartition & domain ) override;

  /**
   * @brief Resize the allocated multidimensional fields
   * @param well the well for which the fields are resized
   */
  void
  ResizeFields( WellElementSubRegion & subRegion );

  /// the max number of fluid phases
  localIndex m_numPhases;

  /// the number of fluid components
  localIndex m_numComponents;

  /// the (uniform) temperature
  real64 m_temperature;

  /// flag indicating whether mass or molar formulation should be used
  integer m_useMass;

  /// list of relative permeability model names per target region
  array1d< string > m_relPermModelNames;

  /// maximum (absolute) change in a component fraction between two Newton iterations
  real64 m_maxCompFracChange;

  /// minimum value of the scaling factor obtained by enforcing maxCompFracChange
  real64 m_minScalingFactor;

  /// flag indicating whether local (cell-wise) chopping of negative compositions is allowed
  integer m_allowCompDensChopping;

  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaResPressure;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
    m_resGlobalCompDensity;

  /// views into other reservoir variable fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_resPhaseVolFrac;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > >
    m_dResPhaseVolFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dResPhaseVolFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dResCompFrac_dCompDens;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_resPhaseMob;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dResPhaseMob_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dResPhaseMob_dCompDens;

  /// views into reservoir material fields

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_resPhaseDens;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_resPhaseVisc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > >
    m_dResPhaseVisc_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > >
    m_dResPhaseVisc_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > m_resPhaseCompFrac;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > >
    m_dResPhaseCompFrac_dPres;
  ElementRegionManager::ElementViewAccessor< arrayView5d< real64 const > >
    m_dResPhaseCompFrac_dComp;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_resPhaseRelPerm;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > >
    m_dResPhaseRelPerm_dPhaseVolFrac;
};

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_COMPOSITIONALMULTIPHASEWELL_HPP_
