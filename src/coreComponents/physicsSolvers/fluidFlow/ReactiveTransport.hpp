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
 * @file ReactiveTransport.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_REACTIVETRANSPORT_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_REACTIVETRANSPORT_HPP_

#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"
#include "physicsSolvers/fluidFlow/GeochemicalModel.hpp"

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
 * @class ReactiveTransport
 *
 * class to perform a reactive transport finite volume solve.
 */
class ReactiveTransport : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  ReactiveTransport( const std::string & name,
                     Group * const parent );


  /// deleted default constructor
  ReactiveTransport() = delete;

  /// deleted copy constructor
  ReactiveTransport( ReactiveTransport const & ) = delete;

  /// default move constructor
  ReactiveTransport( ReactiveTransport && ) = default;

  /// deleted assignment operator
  ReactiveTransport & operator=( ReactiveTransport const & ) = delete;

  /// deleted move operator
  ReactiveTransport & operator=( ReactiveTransport && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~ReactiveTransport() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "ReactiveTransport"; }

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              integer const cycleNumber,
              DomainPartition & domain ) override;

  void PreStepUpdate( real64 const & time_n,
                      real64 const & dt,
                      DomainPartition & domain );

  void PostStepUpdate( real64 const & time_n,
                       real64 const & dt,
                       DomainPartition & domain );

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  AssembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

  virtual void
  ApplyBoundaryConditions( real64 const time,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

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
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleAccumulationTerms( real64 const dt,
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
  void AssembleFluxTerms( real64 const time_n,
                          real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  /**@}*/

  void ResizeFields( MeshLevel & mesh );

  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {

    static constexpr auto reactiveTransportModelString      = "reactiveTransportModel";
    static constexpr auto viscosityString      = "viscosity";

    static constexpr auto bcComponentConcentrationString      = "bcComponentConcentration";

    static constexpr auto deltaComponentConcentrationString      = "deltaComponentConcentration";

  } viewKeysReactiveTransport;

  viewKeyStruct & viewKeys() { return viewKeysReactiveTransport; }
  viewKeyStruct const & viewKeys() const { return viewKeysReactiveTransport; }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysReactiveTransport;

  groupKeyStruct & groupKeys() { return groupKeysReactiveTransport; }
  groupKeyStruct const & groupKeys() const { return groupKeysReactiveTransport; }

  static constexpr localIndex MAX_NUM_COMPONENTS = constitutive::ReactiveFluidBase::MAX_NUM_SPECIES;

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

protected:

  virtual void PostProcessInput() override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( MeshLevel & mesh ) override;

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_deltaPressure;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_componentConcentration;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_deltaComponentConcentration;

  array1d< string > m_reactiveFluidNames;

  localIndex m_numComponents;

  real64 m_viscosity;
};


} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_REACTIVETRANSPORT_HPP_
