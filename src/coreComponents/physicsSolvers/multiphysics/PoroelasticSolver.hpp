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
 * @file PoroelasticSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROELASTICSOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROELASTICSOLVER_HPP_

#include "codingUtilities/EnumStrings.hpp"
#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{


class SolidMechanicsLagrangianFEM;
class SinglePhaseBase;

class PoroelasticSolver : public SolverBase
{
public:
  PoroelasticSolver( const string & name,
                     Group * const parent );
  ~PoroelasticSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName() { return "Poroelastic"; }

  virtual void registerDataOnMesh( dataRepository::Group & MeshBodies ) override;


  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = true ) override;

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  virtual void
  assembleSystem( real64 const time,
                  real64 const dt,
                  DomainPartition & domain,
                  DofManager const & dofManager,
                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                  arrayView1d< real64 > const & localRhs ) override;

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

  virtual void updateState( DomainPartition & domain ) override;

  void updatePermeability( NodeManager const & nodeManger,
                           CellElementSubRegion & subRegion,
                           localIndex const targetIndex ) const;

  virtual void
  implicitStepComplete( real64 const & time_n,
                        real64 const & dt,
                        DomainPartition & domain ) override;

  virtual void
  resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64
  solverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  void updateDeformationForCoupling( DomainPartition & domain );

  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );


  enum class CouplingTypeOption : integer
  {
    FIM,
    SIM_FixedStress
  };



  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * couplingTypeOptionString() { return "couplingTypeOptionEnum"; }
    constexpr static char const * couplingTypeOptionStringString() { return "couplingTypeOption"; }

    constexpr static char const * totalMeanStressString() { return "totalMeanStress"; }
    constexpr static char const * oldTotalMeanStressString() { return "oldTotalMeanStress"; }
    constexpr static char const * dPerm_dDisplacementString() { return "dPerm_dDisplacement"; }

    constexpr static char const * solidSolverNameString() { return "solidSolverName"; }
    constexpr static char const * fluidSolverNameString() { return "fluidSolverName"; }
  };

protected:

  virtual void postProcessInput() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  string m_solidSolverName;
  string m_flowSolverName;

  CouplingTypeOption m_couplingTypeOption;

  // pointer to the flow sub-solver
  SinglePhaseBase * m_flowSolver;

  // pointer to the solid mechanics sub-solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

private:

  void createPreconditioner();

};

ENUM_STRINGS( PoroelasticSolver::CouplingTypeOption,
              "FIM",
              "SIM_FixedStress" );

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_POROELASTICSOLVER_HPP_ */
