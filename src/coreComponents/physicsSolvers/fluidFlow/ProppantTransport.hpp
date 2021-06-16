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
 * @file ProppantTransport.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_HPP_

#include <constitutive/fluid/ParticleFluidBase.hpp>
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "constitutive/fluid/SlurryFluidBase.hpp"

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
 * @class ProppantTransport
 *
 * class to perform a proppant finite volume solve.
 */
class ProppantTransport : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ProppantTransport( const string & name,
                     Group * const parent );


  /// deleted default constructor
  ProppantTransport() = delete;

  /// deleted copy constructor
  ProppantTransport( ProppantTransport const & ) = delete;

  /// default move constructor
  ProppantTransport( ProppantTransport && ) = default;

  /// deleted assignment operator
  ProppantTransport & operator=( ProppantTransport const & ) = delete;

  /// deleted move operator
  ProppantTransport & operator=( ProppantTransport && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ProppantTransport() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "ProppantTransport"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  void preStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      DomainPartition & domain );

  void postStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       DomainPartition & domain );

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  applyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void
  solveSystem( DofManager const & dofManager,
               ParallelMatrix & matrix,
               ParallelVector & rhs,
               ParallelVector & solution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  implicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void assembleAccumulationTerms( real64 const dt,
                                  DomainPartition const & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs );

  /**
   * @brief assembles the flux terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */
  void assembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**@}*/

  void resizeFractureFields( MeshLevel & mesh );

  arrayView1d< string const > const proppantModelNames() const { return m_proppantModelNames; }

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {
    static constexpr char const * proppantNamesString() { return "proppantNames"; }

    // primary solution field
    static constexpr char const * proppantConcentrationString() { return "proppantConcentration"; }
    static constexpr char const * deltaProppantConcentrationString() { return "deltaProppantConcentration"; }
    static constexpr char const * componentConcentrationString() { return "componentConcentration"; }
    static constexpr char const * deltaComponentConcentrationString() { return "deltaComponentConcentration"; }
    static constexpr char const * bcComponentConcentrationString() { return "bcComponentConcentration"; }

    // these are used to store last converged time step values

    static constexpr char const * oldComponentDensityString() { return "oldComponentDensity"; }

    static constexpr char const * updateProppantPackingString() { return "updateProppantPacking"; }
    static constexpr char const * cellBasedFluxString() { return "cellBasedFlux"; }

    static constexpr char const * isProppantBoundaryString() { return "isProppantBoundary"; }
    static constexpr char const * isProppantMobileString() { return "isProppantMobile"; }

    static constexpr char const * proppantPackVolumeFractionString() { return "proppantPackVolumeFraction"; }
    static constexpr char const * proppantExcessPackVolumeString() { return "proppantExcessPackVolume"; }
    static constexpr char const * proppantLiftFluxString() { return "proppantLiftFlux"; }

    static constexpr char const * poroMultiplierString() { return "poroMultiplier"; }
    static constexpr char const * transTMultiplierString() { return "transTMultiplier"; }
    static constexpr char const * bridgingFactorString() { return "bridgingFactor"; }

    static constexpr char const * maxProppantConcentrationString() { return "maxProppantConcentration"; }

    static constexpr char const * proppantDiameterString() { return "proppantDiameter"; }
    static constexpr char const * proppantDensityString() { return "proppantDensity"; }

    static constexpr char const * criticalShieldsNumberString() { return "criticalShieldsNumber"; }
    static constexpr char const * frictionCoefficientString() { return "frictionCoefficient"; }
  };

  static constexpr localIndex MAX_NUM_COMPONENTS = constitutive::ParticleFluidBase::MAX_NUM_COMPONENTS;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  void updateProppantMobility( Group & dataGroup );

  /**
   * @brief Function to update proppant pack volume fraction
   */
  void updateProppantPackVolume( real64 const time_n, real64 const dt, DomainPartition & domain );

  virtual void updateState ( DomainPartition & domain ) override final { GEOSX_UNUSED_VAR( domain ) };

  /**
   * @brief Function to update fluid and proppant properties
   * @param domain the domain
   */
  void updateState( Group & dataGroup, localIndex const targetIndex );

protected:

  virtual void postProcessInput() override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void resetViews( MeshLevel & mesh ) override;

  /**
   * @brief Function to update fluid properties
   * @param domain the domain
   */
  void updateFluidModel( Group & dataGroup, localIndex const targetIndex );

  void updateComponentDensity( Group & dataGroup, localIndex const targetIndex );

  void updateProppantModel( Group & dataGroup, localIndex const targetIndex );

  /**
   * @brief Function to update cell-based fluid flux
   */
  void updateCellBasedFlux( real64 const time_n,
                            DomainPartition & domain );

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_proppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_cellBasedFlux;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_proppantLiftFlux;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_proppantPackVolumeFraction;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_proppantExcessPackVolume;
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > m_isProppantBoundaryElement;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_transTMultiplier;
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > m_isProppantMobile;

  /// views into material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_density;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dDensity_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_componentDensity;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dComponentDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView4d< real64 const > > m_dComponentDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_fluidDensity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dFluidDensity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dFluidDensity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_fluidViscosity;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_viscosity;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dViscosity_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dViscosity_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dViscosity_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_settlingFactor;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_dSettlingFactor_dPressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_dSettlingFactor_dProppantConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dSettlingFactor_dComponentConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_collisionFactor;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_dCollisionFactor_dProppantConcentration;

  array1d< string > m_proppantModelNames;

  integer m_numComponents;

  integer m_updateProppantPacking;
  real64 m_proppantPackPermeability;
  real64 m_bridgingFactor;
  real64 m_minAperture;
  real64 m_maxProppantConcentration;
  real64 m_proppantDiameter;
  real64 m_proppantDensity;
  real64 m_criticalShieldsNumber;
  real64 m_frictionCoefficient;
};


} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_PROPPANTTRANSPORT_HPP_
