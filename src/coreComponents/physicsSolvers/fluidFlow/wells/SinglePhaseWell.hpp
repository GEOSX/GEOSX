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
 * @file SinglePhaseWell.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_

#include "WellSolverBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;
}

namespace constitutive
{
class SingleFluidBase;
}
class WellElementSubRegion;

/**
 * @class SinglePhaseWell
 *
 * A single-phase well solver
 */
class SinglePhaseWell : public WellSolverBase
{
public:
  // define the column offset of the derivatives
  struct ColOffset
  {
    static constexpr integer DPRES = 0;
    static constexpr integer DRATE = 1;
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
  SinglePhaseWell( const string & name, Group * const parent );

  /// deleted default constructor
  SinglePhaseWell() = delete;

  /// deleted copy constructor
  SinglePhaseWell( SinglePhaseWell const & ) = delete;

  /// default move constructor
  SinglePhaseWell( SinglePhaseWell && ) = default;

  /// deleted assignment operator
  SinglePhaseWell &
  operator=( SinglePhaseWell const & ) = delete;

  /// deleted move operator
  SinglePhaseWell &
  operator=( SinglePhaseWell && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseWell() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "SinglePhaseWell";
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

  virtual string
  WellElementDofName() const override
  {
    return viewKeyStruct::dofFieldString;
  }

  virtual string
  ResElementDofName() const override
  {
    return SinglePhaseBase::viewKeyStruct::pressureString;
  }

  virtual localIndex
  NumFluidComponents() const override
  {
    return 1;
  }

  virtual localIndex
  NumFluidPhases() const override
  {
    return 1;
  }

  /**
   * @brief Update fluid constitutive model state
   * @param dataGroup group that contains the fields
   */
  virtual void
  UpdateFluidModel( WellElementSubRegion & subRegion,
                    localIndex const targetIndex ) const;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models) on the well
   * @param subRegion the well subRegion containing the well elements and their associated fields
   */
  virtual void
  UpdateState( WellElementSubRegion & subRegion,
               localIndex const targetIndex ) override;

  /**
   * @brief assembles the flux terms for all connections between well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void
  AssembleFluxTerms( real64 const time_n,
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

  struct viewKeyStruct : WellSolverBase::viewKeyStruct
  {
    static constexpr auto dofFieldString = "singlePhaseWellVars";

    // primary solution field
    static constexpr auto pressureString =
      SinglePhaseBase::viewKeyStruct::pressureString;
    static constexpr auto deltaPressureString =
      SinglePhaseBase::viewKeyStruct::deltaPressureString;
    static constexpr auto connRateString = "connectionRate";
    static constexpr auto deltaConnRateString = "deltaConnectionRate";

    // perforation rates
    static constexpr auto perforationRateString = "perforationRate";
    static constexpr auto dPerforationRate_dPresString =
      "dPerforationRate_dPres";
  } viewKeysSinglePhaseWell;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseWell;

protected:
  virtual void
  PostProcessInput() override;

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

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
   * @param domain the domain containing the well manager to access individual wells
   */
  void
  ResetViews( DomainPartition & domain ) override;

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  void
  InitializeWells( DomainPartition & domain ) override;

private:
  /// views into reservoir primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_resPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaResPressure;

  /// views into reservoir material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_resDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dResDens_dPres;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_resViscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dResVisc_dPres;
};

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_SINGLEPHASEWELL_HPP_
