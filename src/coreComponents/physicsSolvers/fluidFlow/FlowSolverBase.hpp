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
 * @file FlowSolverBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{

/**
 * @class FlowSolverBase
 *
 * Base class for finite volume fluid flow solvers.
 * Provides some common features
 */
class FlowSolverBase : public SolverBase
{
public:
/**
 * @brief main constructor for Group Objects
 * @param name the name of this instantiation of Group in the repository
 * @param parent the parent group of this instantiation of Group
 */
  FlowSolverBase( const string & name,
                  Group * const parent );


  /// deleted default constructor
  FlowSolverBase() = delete;

  /// deleted copy constructor
  FlowSolverBase( FlowSolverBase const & ) = delete;

  /// default move constructor
  FlowSolverBase( FlowSolverBase && ) = default;

  /// deleted assignment operator
  FlowSolverBase & operator=( FlowSolverBase const & ) = delete;

  /// deleted move operator
  FlowSolverBase & operator=( FlowSolverBase && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~FlowSolverBase() override;

  virtual void registerDataOnMesh( Group & MeshBodies ) override;

  void setPoroElasticCoupling() { m_poroElasticFlag = 1; }

  void setReservoirWellsCoupling() { m_coupledWellsFlag = 1; }

  arrayView1d< string const > fluidModelNames() const { return m_fluidModelNames; }

  arrayView1d< string const > solidModelNames() const { return m_solidModelNames; }

  virtual std::vector< string > getConstitutiveRelations( string const & regionName ) const override;


  localIndex numDofPerCell() const { return m_numDofPerCell; }

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // input data
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * permeabilityString() { return "permeability"; }

    // gravity term precomputed values
    static constexpr char const * gravityCoefString() { return "gravityCoefficient"; }

    // misc inputs
    static constexpr char const * fluidNamesString() { return "fluidNames"; }
    static constexpr char const * solidNamesString() { return "solidNames"; }
    static constexpr char const * permeabilityNamesString() { return "permeabilityNames"; }
    static constexpr char const * pressureString() { return "pressure"; }
    static constexpr char const * deltaPressureString() { return "deltaPressure"; }
    static constexpr char const * deltaVolumeString() { return "deltaVolume"; }
    static constexpr char const * aperture0String() { return "aperture_n"; }
    static constexpr char const * effectiveApertureString() { return "effectiveAperture"; }
    static constexpr char const * inputFluxEstimateString() { return "inputFluxEstimate"; }
    static constexpr char const * meanPermCoeffString() { return "meanPermCoeff"; }
  };

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void resetViews( MeshLevel & mesh );

  virtual void setUpDflux_dApertureMatrix( DomainPartition & domain,
                                           DofManager const & dofManager,
                                           CRSMatrix< real64, globalIndex > & localMatrix );


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

private:

  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain partition
   */


protected:

  virtual void precomputeData( MeshLevel & mesh );

  virtual void postProcessInput() override;

  virtual void initializePreSubGroups() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

  /// name of the fluid constitutive model
  array1d< string > m_fluidModelNames;

  /// name of the solid constitutive model
  array1d< string > m_solidModelNames;

  /// name of the permeability constituive model
  array1d< string > m_permeabilityModelNames;

  /// flag to determine whether or not coupled with solid solver
  integer m_poroElasticFlag;

  /// flag to determine whether or not coupled with wells
  integer m_coupledWellsFlag;

  /// the number of Degrees of Freedom per cell
  integer m_numDofPerCell;

  std::unique_ptr< CRSMatrix< real64, localIndex > > m_derivativeFluxResidual_dAperture;

  real64 m_fluxEstimate;

  real64 m_meanPermCoeff;

  /// views into constant data fields
  ElementRegionManager::ElementViewAccessor< arrayView1d< integer const > > m_elemGhostRank;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_volume;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_gravCoef;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_porosityRef;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_elementArea;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_elementAperture0;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_elementAperture;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > >  m_effectiveAperture;

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  m_elementSeparationCoefficient;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > >  m_element_dSeparationCoefficient_dAperture;
#endif

};

}

#endif //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_FLOWSOLVERBASE_HPP_
