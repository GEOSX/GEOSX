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
 * @file HydrofractureSolver.hpp
 *
 */

#ifndef GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_
#define GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_

#include "physicsSolvers/multiphysics/PoroelasticSolver.hpp"

namespace geosx
{

class SurfaceGenerator;


class HydrofractureSolver : public PoroelasticSolver
{
public:
  HydrofractureSolver( const string & name,
                       Group * const parent );

  ~HydrofractureSolver() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string catalogName()
  {
    return "Hydrofracture";
  }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  virtual void RegisterDataOnMesh( Group & MeshBodies ) override final;
#endif

  virtual void setupDofs( DomainPartition const & domain,
                          DofManager & dofManager ) const override;

  virtual void setupSystem( DomainPartition & domain,
                            DofManager & dofManager,
                            CRSMatrix< real64, globalIndex > & localMatrix,
                            array1d< real64 > & localRhs,
                            array1d< real64 > & localSolution,
                            bool const setSparsity = true ) override;

  virtual void
  implicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void assembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void applyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  scalingForSystemSolution( DomainPartition const & domain,
                            DofManager const & dofManager,
                            arrayView1d< real64 const > const & localSolution ) override;

  virtual void
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void resetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             int const cycleNumber,
                             DomainPartition & domain ) override;

  virtual void setNextDt( real64 const & currentDt,
                          real64 & nextDt ) override;


  virtual real64 explicitStep( real64 const & time_n,
                               real64 const & dt,
                               integer const cycleNumber,
                               DomainPartition & domain ) override;

  void updateDeformationForCoupling( DomainPartition & domain );

//  void ApplyFractureFluidCoupling( DomainPartition * const domain,
//                                   systemSolverInterface::EpetraBlockSystem & blockSystem );

  void assembleForceResidualDerivativeWrtPressure( DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs );

  void assembleFluidMassResidualDerivativeWrtDisplacement( DomainPartition const & domain,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix );


  real64 splitOperatorStep( real64 const & time_n,
                            real64 const & dt,
                            integer const cycleNumber,
                            DomainPartition & domain );

  void initializeNewFaceElements( DomainPartition const & domain );


  std::unique_ptr< CRSMatrix< real64, localIndex > > & getRefDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture;
  }

  CRSMatrixView< real64, localIndex const > getDerivativeFluxResidual_dAperture()
  {
    return m_derivativeFluxResidual_dAperture->toViewConstSizes();
  }

  CRSMatrixView< real64 const, localIndex const > getDerivativeFluxResidual_dAperture() const
  {
    return m_derivativeFluxResidual_dAperture->toViewConst();
  }


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * contactRelationNameString() { return "contactRelationName"; }

    constexpr static char const * surfaceGeneratorNameString() { return "surfaceGeneratorName"; }

    constexpr static char const * maxNumResolvesString() { return "maxNumResolves"; }

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
    constexpr static char const * separationCoeff0String() { return "separationCoeff0"; }
    constexpr static char const * apertureAtFailureString() { return "apertureAtFailure"; }
#endif
  };

protected:
  virtual void postProcessInput() override final;

  virtual void
  initializePostInitialConditionsPreSubGroups() override final;

  /**
   * @Brief add the nnz induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param rowLenghts the nnz in each row
   */
  void addFluxApertureCouplingNNZ( DomainPartition & domain,
                                   DofManager & dofManager,
                                   arrayView1d< localIndex > const & rowLengths ) const;


  /**
   * @Brief add the sparsity pattern induced by the flux-aperture coupling
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param pattern the sparsity pattern
   */
  void addFluxApertureCouplingSparsityPattern( DomainPartition & domain,
                                               DofManager & dofManager,
                                               SparsityPatternView< globalIndex > const & pattern ) const;


  void setUpDflux_dApertureMatrix( DomainPartition & domain,
                                   DofManager const & dofManager,
                                   CRSMatrix< real64, globalIndex > & localMatrix );



private:

  // name of the contact relation
  string m_contactRelationName;

  /// name of the surface generator
  string m_surfaceGeneratorName;

  /// pointer to the surface generator
  SurfaceGenerator * m_surfaceGenerator;

  std::unique_ptr< ParallelMatrix > m_blockDiagUU;

  // it is only important for this case.
  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

  integer m_maxNumResolves;
  integer m_numResolves[2];
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_MULTIPHYSICS_HYDROFRACTURESOLVER_HPP_ */
