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
 * @file PoroelasticSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
class SolidMechanicsLagrangianFEM;
class FlowSolverBase;

class PoroelasticSolver : public SolverBase
{
public:
  PoroelasticSolver( const std::string & name, Group * const parent );
  ~PoroelasticSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "Poroelastic";
  }

  virtual void
  RegisterDataOnMesh(
    dataRepository::Group * const MeshBodies ) override final;

  virtual void
  SetupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  AssembleSystem(
    real64 const time,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  void
  AssembleCouplingTerms(
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs );

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

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  ImplicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override final;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  void
  UpdateDeformationForCoupling( DomainPartition & domain );

  real64
  SplitOperatorStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain );

  enum class couplingTypeOption : int
  {
    FIM,
    SIM_FixedStress
  };

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto couplingTypeOptionString = "couplingTypeOptionEnum";
    constexpr static auto couplingTypeOptionStringString = "couplingTypeOption";

    constexpr static auto totalMeanStressString = "totalMeanStress";
    constexpr static auto oldTotalMeanStressString = "oldTotalMeanStress";

    constexpr static auto solidSolverNameString = "solidSolverName";
    constexpr static auto fluidSolverNameString = "fluidSolverName";
  } poroElasticSolverViewKeys;

  SolidMechanicsLagrangianFEM *
  getSolidSolver()
  {
    return this->getParent()
      ->GetGroup( m_solidSolverName )
      ->group_cast< SolidMechanicsLagrangianFEM * >();
  }
  SolidMechanicsLagrangianFEM const *
  getSolidSolver() const
  {
    return this->getParent()
      ->GetGroup( m_solidSolverName )
      ->group_cast< SolidMechanicsLagrangianFEM const * >();
  }

  FlowSolverBase *
  getFlowSolver()
  {
    return this->getParent()
      ->GetGroup( m_flowSolverName )
      ->group_cast< FlowSolverBase * >();
  }
  FlowSolverBase const *
  getFlowSolver() const
  {
    return this->getParent()
      ->GetGroup( m_flowSolverName )
      ->group_cast< FlowSolverBase const * >();
  }

protected:
  virtual void
  PostProcessInput() override final;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    dataRepository::Group * const problemManager ) override final;

private:
  void
  CreatePreconditioner();

  string m_solidSolverName;
  string m_flowSolverName;
  string m_couplingTypeOptionString;
  couplingTypeOption m_couplingTypeOption;

  // pointer to the flow sub-solver
  FlowSolverBase * m_flowSolver;

  // pointer to the solid mechanics sub-solver
  SolidMechanicsLagrangianFEM * m_solidSolver;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_COUPLEDSOLVERS_POROELASTICSOLVER_HPP_ */
