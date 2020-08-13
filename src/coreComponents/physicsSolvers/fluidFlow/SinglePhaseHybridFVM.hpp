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
 * @file SinglePhaseHybridFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseHybridFVMKernels.hpp"

namespace geosx
{
/**
 * @class SinglePhaseHybridFVM
 *
 * class to assemble the single-phase flow equations based
 * on a mixed hybrid formulation known as the mimetic method
 */
class SinglePhaseHybridFVM : public SinglePhaseBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseHybridFVM( const std::string & name, Group * const parent );

  /// deleted default constructor
  SinglePhaseHybridFVM() = delete;

  /// deleted copy constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM const & ) = delete;

  /// default move constructor
  SinglePhaseHybridFVM( SinglePhaseHybridFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseHybridFVM &
  operator=( SinglePhaseHybridFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseHybridFVM &
  operator=( SinglePhaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseHybridFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "SinglePhaseHybridFVM";
  }

  virtual void
  RegisterDataOnMesh( Group * const MeshBodies ) override;

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
   * @brief assembles the flux terms for all cells
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

  /**@}*/

  struct viewKeyStruct : SinglePhaseBase::viewKeyStruct
  {
    // primary face-based field
    static constexpr auto deltaFacePressureString = "deltaFacePressure";

  } viewKeysSinglePhaseHybridFVM;

  viewKeyStruct &
  viewKeys()
  {
    return viewKeysSinglePhaseHybridFVM;
  }

  viewKeyStruct const &
  viewKeys() const
  {
    return viewKeysSinglePhaseHybridFVM;
  }

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysSinglePhaseHybridFVM;

  groupKeyStruct &
  groupKeys()
  {
    return groupKeysSinglePhaseHybridFVM;
  }

  groupKeyStruct const &
  groupKeys() const
  {
    return groupKeysSinglePhaseHybridFVM;
  }

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    dataRepository::Group * const rootGroup ) override;

private:
  /// Dof key for the member functions that do not have access to the coupled Dof manager
  string m_faceDofKey;

  /// relative tolerance (redundant with FluxApproximationBase)
  real64 m_areaRelTol;

  /// region filter used in flux assembly
  SortedArray< localIndex > m_regionFilter;
};

} /* namespace geosx */

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEHYBRIDFVM_HPP_
