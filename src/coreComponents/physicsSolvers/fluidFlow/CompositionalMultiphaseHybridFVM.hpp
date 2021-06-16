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
 * @file CompositionalMultiphaseHybridFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_

#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBase.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp"

namespace geosx
{

namespace dataRepository
{
class Group;
}

namespace constitutive
{
class MultiFluidBase;
}

/**
 * @class CompositionalMultiphaseHybridFVM
 *
 * A compositional multiphase solver
 * using only cell-centered and face-centered variables
 */
class CompositionalMultiphaseHybridFVM : public CompositionalMultiphaseBase
{
public:

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  CompositionalMultiphaseHybridFVM( const string & name,
                                    Group * const parent );

  /// deleted default constructor
  CompositionalMultiphaseHybridFVM() = delete;

  /// deleted copy constructor
  CompositionalMultiphaseHybridFVM( CompositionalMultiphaseHybridFVM const & ) = delete;

  /// default move constructor
  CompositionalMultiphaseHybridFVM( CompositionalMultiphaseHybridFVM && ) = default;

  /// deleted assignment operator
  CompositionalMultiphaseHybridFVM & operator=( CompositionalMultiphaseHybridFVM const & ) = delete;

  /// deleted move operator
  CompositionalMultiphaseHybridFVM & operator=( CompositionalMultiphaseHybridFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~CompositionalMultiphaseHybridFVM() override = default;

  /**
   * @brief name of the solver in the object catalog
   * @return string that contains the catalog name to generate a new object through the object catalog.
   */
  static string catalogName() { return "CompositionalMultiphaseHybridFVM"; }

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

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
  applyBoundaryConditions( real64 const time_n,
                           real64 const dt,
                           DomainPartition & domain,
                           DofManager const & dofManager,
                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                           arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  calculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual bool
  checkSystemSolution( DomainPartition const & domain,
                       DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor ) override;

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
   * @brief assembles the flux terms for all cells
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assembleFluxTerms( real64 const dt,
                     DomainPartition & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) const override;

  /**
   * @brief Recompute phase mobility from constitutive and primary variables
   * @param domain the domain containing the mesh and fields
   */
  virtual void
  updatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const override;

  /**@}*/

  struct viewKeyStruct : CompositionalMultiphaseBase::viewKeyStruct
  {
    static constexpr char const * faceDofFieldString() { return "faceCenteredVariables"; }

    // inputs
    static constexpr char const * maxRelativePresChangeString() { return "maxRelativePressureChange"; }

    // primary face-based field
    static constexpr char const * deltaFacePressureString() { return "deltaFacePressure"; }

    // auxiliary data for the buoyancy term
    // TODO: change the name
    static constexpr char const * mimGravityCoefString() { return "mimGravityCoefficient"; }

  };

  virtual void initializePostInitialConditionsPreSubGroups() override;

  virtual void initializePreSubGroups() override;

protected:

  /// precompute the minGravityCoefficient for the buoyancy term
  void precomputeData( MeshLevel & mesh ) override;

private:

  /// maximum relative face pressure change between two Newton iterations
  real64 m_maxRelativePresChange;

  /// tolerance used in the  computation of the transmissibility matrix
  real64 m_lengthTolerance;

  /// name of the transmissibility multiplier field
  string m_transMultName;

  /// region filter used in flux assembly
  SortedArray< localIndex > m_regionFilter;

};

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_COMPOSITIONALMULTIPHASEHYBRIDFVM_HPP_
