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
 * @file GeochemicalModel.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_
#define SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_

#include "FlowSolverBase.hpp"
#include "constitutive/fluid/ReactiveFluidBase.hpp"

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
 * @class GeochemicalModel
 *
 * class to perform a proppant finite volume solve.
 */
class GeochemicalModel : public FlowSolverBase
{
public:
  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  GeochemicalModel( const std::string & name,
                    Group * const parent );


  /// deleted default constructor
  GeochemicalModel() = delete;

  /// deleted copy constructor
  GeochemicalModel( GeochemicalModel const & ) = delete;

  /// default move constructor
  GeochemicalModel( GeochemicalModel && ) = default;

  /// deleted assignment operator
  GeochemicalModel & operator=( GeochemicalModel const & ) = delete;

  /// deleted move operator
  GeochemicalModel & operator=( GeochemicalModel && ) = delete;

  /**
   * @brief default destructor
   */
  virtual ~GeochemicalModel() override = default;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string CatalogName() { return "GeochemicalModel"; }

  virtual void InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void RegisterDataOnMesh( Group * const MeshBodies ) override;

  virtual real64 SolverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

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

  virtual void AssembleSystem( real64 const time,
                               real64 const dt,
                               DomainPartition & domain,
                               DofManager const & dofManager,
                               CRSMatrixView< real64, globalIndex const > const & localMatrix,
                               arrayView1d< real64 > const & localRhs ) override;

  virtual void ApplyBoundaryConditions( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs ) override;

  virtual real64
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs ) override;

  virtual void SolveSystem( DofManager const & dofManager,
                            ParallelMatrix & matrix,
                            ParallelVector & rhs,
                            ParallelVector & solution ) override;

  virtual void
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;


  virtual void ResetStateToBeginningOfStep( DomainPartition & domain ) override;

  virtual void ImplicitStepComplete( real64 const & time,
                                     real64 const & dt,
                                     DomainPartition & domain ) override;

  /**
   * @brief assembles the accumulation terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleAccumulationTerms( DomainPartition & domain,
                                  DofManager const & dofManager,
                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                  arrayView1d< real64 > const & localRhs,
                                  bool isInitialization );

  /**
   * @brief assembles the flux terms for all cells
   * @param domain the physical domain object
   * @param blockSystem the entire block system
   * @param time_n previous time value
   * @param dt time step
   */

  void AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          real64 const dt,
                          DomainPartition const & domain,
                          DofManager const & dofManager,
                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                          arrayView1d< real64 > const & localRhs );

  void ResizeFields( MeshLevel * const meshLevel );

  void WriteSpeciesToFile( DomainPartition * const domain, real64 const & time );

  void CalculateTotalConcentration( DomainPartition & domain );

  /**@}*/


  struct viewKeyStruct : FlowSolverBase::viewKeyStruct
  {

    static constexpr auto reactiveFluidNamesString      = "reactiveFluidNames";
    // primary solution field
    static constexpr auto temperatureString      = "temperature";
    static constexpr auto deltaTemperatureString = "deltaTemperature";

    static constexpr auto concentrationString      = "concentration";
    static constexpr auto concentrationOutString      = "concentrationOut";
    static constexpr auto deltaConcentrationString      = "deltaConcentration";

    static constexpr auto totalConcentrationString      = "totalConcentration";
    static constexpr auto concentrationNewString      = "concentrationNew";
    static constexpr auto outputSpeciesFileNameString      = "outputSpeciesFileName";
    static constexpr auto geochemicalModelString      = "geochemicalModel";

  } viewKeysGeochemicalModel;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {} groupKeysGeochemicalModel;

protected:

  virtual void InitializePostInitialConditions_PreSubGroups( dataRepository::Group * const rootGroup ) override;

private:

  /**
   * @brief Setup stored views into domain data for the current step
   */
  void ResetViews( MeshLevel & mesh ) override;

  /**
   * @brief Function to update all constitutive models
   * @param domain the domain
   */
  void UpdateReactiveFluidModel( Group * const dataGroup, localIndex const targetIndex ) const;

  void UpdateState( Group * dataGroup, localIndex const targetIndex ) const;

  /// views into primary variable fields

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_pressure;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaPressure;

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_temperature;
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_deltaTemperature;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_concentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_concentrationOut;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_deltaConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_totalConcentration;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_totalConcentrationOut;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_concentrationNew;

  /// views into material fields

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 const > > m_dependentConc;
  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > m_dDependentConc_dConc;

  array1d< string > m_reactiveFluidNames;
  array1d< localIndex > m_numBasisSpecies;
  array1d< localIndex > m_numDependentSpecies;

  string m_outputSpeciesFileName;

  //Below is not used in GeochemicalModel model

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > m_porosity;

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_pvMult;
  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > m_dPvMult_dPres;

};


} /* namespace geosx */

#endif //SRC_COMPONENTS_CORE_SRC_PHYSICSSOLVERS_GEOCHEMICALMODEL_HPP_
