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

#ifndef GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSION_HPP_
#define GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSION_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mainInterface/GeosxState.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

class ThermalDiffusion : public SolverBase
{
public:

  ThermalDiffusion( string const & name,
                    Group * const parent );

  virtual ~ThermalDiffusion() override = default;

  static string catalogName() { return "ThermalDiffusion"; }

private:

  struct viewKeyStruct
  {
    static constexpr char const * thermalDiffusionString() { return "thermalDiffusion"; }
    static constexpr char const * newDeltaTemperatureString() { return "newDeltaTemperature"; }
  };

  virtual void registerDataOnMesh( Group & meshBodies ) override;

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = false ) override;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void solveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void applySystemSolution( DofManager const & dofManager,
                                    arrayView1d< real64 const > const & localSolution,
                                    real64 const scalingFactor,
                                    DomainPartition & domain ) override;

  virtual void implicitStepSetup( real64 const & time_n,
                                  real64 const & dt,
                                  DomainPartition & domain ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  virtual real64 calculateResidualNorm( DomainPartition const & domain,
                                        DofManager const & dofManager,
                                        arrayView1d< real64 const > const & localRhs ) override;

  virtual void implicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void applyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  // Diffusion coefficient
  real64 m_thermalDiffusion;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SIMPLEPDE_THERMALDIFFUSION_HPP_ */
