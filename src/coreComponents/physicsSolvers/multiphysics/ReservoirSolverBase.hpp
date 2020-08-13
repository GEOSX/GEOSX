/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file ReservoirSolverBase.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
class FlowSolverBase;
class WellSolverBase;

class ReservoirSolverBase : public SolverBase
{
public:
  /**
   * @brief main constructor for ManagedGroup Objects
   * @param name the name of this instantiation of ManagedGroup in the repository
   * @param parent the parent group of this instantiation of ManagedGroup
   */
  ReservoirSolverBase( const std::string & name, Group * const parent );

  /**
   * @brief default destructor
   */
  virtual ~ReservoirSolverBase() override;

  /// deleted copy constructor
  ReservoirSolverBase( ReservoirSolverBase const & ) = delete;

  /// default move constructor
  ReservoirSolverBase( ReservoirSolverBase && ) = default;

  /// deleted assignment operator
  ReservoirSolverBase &
  operator=( ReservoirSolverBase const & ) = delete;

  /// deleted move operator
  ReservoirSolverBase &
  operator=( ReservoirSolverBase && ) = delete;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "Reservoir";
  }

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual void
  AssembleSystem(
    real64 const time,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

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

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  /**@}*/

  /**
   * @Brief assembles the perforation rate terms
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleCouplingTerms(
    real64 const time_n,
    real64 const dt,
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) = 0;

  FlowSolverBase *
  GetFlowSolver() const
  {
    return m_flowSolver;
  }

  WellSolverBase *
  GetWellSolver() const
  {
    return m_wellSolver;
  }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // solver that assembles the reservoir equations
    constexpr static auto flowSolverNameString = "flowSolverName";

    // solver that assembles the well
    constexpr static auto wellSolverNameString = "wellSolverName";

  } reservoirWellsSolverViewKeys;

protected:
  virtual void
  InitializePostInitialConditions_PreSubGroups(
    Group * const rootGroup ) override;

  virtual void
  PostProcessInput() override;

  void
  AddCouplingNumNonzeros( DomainPartition & domain,
                          DofManager & dofManager,
                          arrayView1d< localIndex > const & rowLengths ) const;

  /**
   * @Brief add the sparsity pattern induced by the perforations
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  virtual void
  AddCouplingSparsityPattern(
    DomainPartition const & domain,
    DofManager const & dofManager,
    SparsityPatternView< globalIndex > const & pattern ) const = 0;

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void
  ResetViews( DomainPartition * const domain );

  /// solver that assembles the reservoir equations
  string m_flowSolverName;

  /// solver that assembles the well equations and compute perforation rates
  string m_wellSolverName;

  /// pointer to the flow sub-solver
  FlowSolverBase * m_flowSolver;

  /// pointer to the well sub-solver
  WellSolverBase * m_wellSolver;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_RESERVOIRSOLVERBASE_HPP_ */
