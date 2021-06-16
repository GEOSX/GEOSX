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
 * @file SinglePhaseFVM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_

#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseProppantBase.hpp"

namespace geosx
{


/**
 * @class SinglePhaseFVM
 *
 * class to perform a single phase finite volume solve
 * using only cell-centered variables
 * works with both TPFA and MPFA
 */
template< typename BASE = SinglePhaseBase >
class SinglePhaseFVM : public BASE
{
public:


  // Aliasing public/protected members/methods of Group so we don't
  // have to use this->member etc.
  using BASE::getLogLevel;

  // Aliasing public/protected members/methods of SolverBase so we don't
  // have to use this->member etc.
  using BASE::m_cflFactor;
  using BASE::m_maxStableDt;
  using BASE::m_nextDt;
  using BASE::m_discretizationName;
  using BASE::targetRegionNames;
  using BASE::forTargetRegions;
  using BASE::forTargetRegionsComplete;
  using BASE::forTargetSubRegions;
  using BASE::forTargetSubRegionsComplete;
  using BASE::m_dofManager;
  using BASE::m_matrix;
  using BASE::m_rhs;
  using BASE::m_solution;
  using BASE::m_localMatrix;
  using BASE::m_localRhs;
  using BASE::m_localSolution;
  using BASE::m_linearSolverParameters;
  using BASE::m_nonlinearSolverParameters;

  // Aliasing public/protected members/methods of FlowSolverBase so we don't
  // have to use this->member etc.
  using BASE::m_fluidModelNames;
  using BASE::m_solidModelNames;
  using BASE::m_permeabilityModelNames;
  using BASE::m_poroElasticFlag;
  using BASE::m_coupledWellsFlag;
  using BASE::m_numDofPerCell;
  using BASE::m_fluxEstimate;
  using BASE::m_elemGhostRank;
  using BASE::m_volume;
  using BASE::m_gravCoef;


  // Aliasing public/protected members/methods of SinglePhaseBase so we don't
  // have to use this->member etc.
  using BASE::m_pressure;
  using BASE::m_deltaPressure;
  using BASE::m_deltaVolume;
  using BASE::m_permeability;
  using BASE::m_dPerm_dPressure;
  using BASE::m_mobility;
  using BASE::m_dMobility_dPres;
  using BASE::m_density;
  using BASE::m_dDens_dPres;
  using BASE::m_viscosity;
  using BASE::m_dVisc_dPres;
  using BASE::m_transTMultiplier;

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  SinglePhaseFVM( const string & name,
                  dataRepository::Group * const parent );


  /// deleted default constructor
  SinglePhaseFVM() = delete;

  /// deleted copy constructor
  SinglePhaseFVM( SinglePhaseFVM const & ) = delete;

  /// default move constructor
  SinglePhaseFVM( SinglePhaseFVM && ) = default;

  /// deleted assignment operator
  SinglePhaseFVM & operator=( SinglePhaseFVM const & ) = delete;

  /// deleted move operator
  SinglePhaseFVM & operator=( SinglePhaseFVM && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~SinglePhaseFVM() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  template< typename _BASE=BASE >
  static
  typename std::enable_if< std::is_same< _BASE, SinglePhaseBase >::value, string >::type
  catalogName()
  {
    return "SinglePhaseFVM";
  }

  template< typename _BASE=BASE >
  static
  typename std::enable_if< std::is_same< _BASE, SinglePhaseProppantBase >::value, string >::type
  catalogName()
  {
    return "SinglePhaseProppantFVM";
  }

  /**
   * @defgroup Solver Interface Function
   *
   * This function provides the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  setupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  setupSystem( DomainPartition & domain,
               DofManager & dofManager,
               CRSMatrix< real64, globalIndex > & localMatrix,
               array1d< real64 > & localRhs,
               array1d< real64 > & localSolution,
               bool const setSparsity = true ) override;

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
  applySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
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
  assembleFluxTerms( real64 const time_n,
                     real64 const dt,
                     DomainPartition & domain,
                     DofManager const & dofManager,
                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                     arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the flux terms for all cells for the poroelastic case
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param localMatrix the system matrix
   * @param localRhs the system right-hand side vector
   */
  virtual void
  assemblePoroelasticFluxTerms( real64 const time_n,
                                real64 const dt,
                                DomainPartition & domain,
                                DofManager const & dofManager,
                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                arrayView1d< real64 > const & localRhs,
                                string const & jumpDofKey ) override;

  void assemblePoroelasticFluxTerms( real64 const time_n,
                                     real64 const dt,
                                     DomainPartition & domain,
                                     DofManager const & dofManager,
                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                     arrayView1d< real64 > const & localRhs,
                                     CRSMatrixView< real64, localIndex const > const & dR_dAper );

  /**@}*/

  virtual void initializePreSubGroups() override;

private:

  /**
   * @brief Function to perform the application of Dirichlet BCs on faces
   * @param time_n current time
   * @param dt time step
   * @param faceSet degree-of-freedom manager associated with the linear system
   * @param domain the domain
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  void applyFaceDirichletBC( real64 const time_n,
                             real64 const dt,
                             DofManager const & faceSet,
                             DomainPartition & domain,
                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                             arrayView1d< real64 > const & localRhs );

  // no data needed here, see SinglePhaseBase

};

} /* namespace geosx */

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_SINGLEPHASEFVM_HPP_
