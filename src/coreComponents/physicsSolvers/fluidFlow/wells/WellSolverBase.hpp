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
 * @file WellSolverBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_

#include "physicsSolvers/SolverBase.hpp"
#include "physicsSolvers/fluidFlow/FlowSolverBase.hpp"

namespace geosx
{
namespace dataRepository
{
class Group;
}
class DomainPartition;
class WellControls;
class WellElementSubRegion;

/**
 * @class WellSolverBase
 *
 * Base class for well solvers.
 * Provides some common features
 */
class WellSolverBase : public SolverBase
{
public:
  // tag to access well and reservoir elements in perforation rates computation
  struct SubRegionTag
  {
    static constexpr integer RES = 0;
    static constexpr integer WELL = 1;
  };

  // tag to access the next and current well elements of a connection
  struct ElemTag
  {
    static constexpr integer CURRENT = 0;
    static constexpr integer NEXT = 1;
  };

  /**
   * @brief main constructor for Group Objects
   * @param name the name of this instantiation of Group in the repository
   * @param parent the parent group of this instantiation of Group
   */
  WellSolverBase( const std::string & name, Group * const parent );

  /// default destructor
  virtual ~WellSolverBase() override;

  /// deleted default constructor
  WellSolverBase() = delete;

  /// deleted copy constructor
  WellSolverBase( WellSolverBase const & ) = delete;

  /// default move constructor
  WellSolverBase( WellSolverBase && ) = default;

  /// deleted assignment operator
  WellSolverBase &
  operator=( WellSolverBase const & ) = delete;

  /// deleted move operator
  WellSolverBase &
  operator=( WellSolverBase && ) = delete;

  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /// Expand catalog for schema generation
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief setter for the name of the flow solver (needed to use the flow kernels like UpdateFluid)
   * @param name the name of the flow solver
   */
  void
  SetFlowSolverName( string const & name )
  {
    m_flowSolverName = name;
  }

  /**
   * @brief getter for the name of the flow solver (used in UpdateState)
   * @return a string containing the name of the flow solver
   */
  string const &
  GetFlowSolverName() const
  {
    return m_flowSolverName;
  }

  /**
   * @brief getter for the number of degrees of freedom per well element
   * @return the number of dofs
   */
  localIndex
  NumDofPerWellElement() const
  {
    return m_numDofPerWellElement;
  }

  /**
   * @brief getter for the number of degrees of freedom per mesh element
   * @return the number of dofs
   */
  localIndex
  NumDofPerResElement() const
  {
    return m_numDofPerResElement;
  }

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string
  WellElementDofName() const = 0;

  /**
   * @brief get the name of DOF defined on well elements
   * @return name of the DOF field used by derived solver type
   */
  virtual string
  ResElementDofName() const = 0;

  /**
   * @brief const getter for the number of fluid components
   * @return the number of fluid components
   */
  virtual localIndex
  NumFluidComponents() const = 0;

  /**
   * @brief Get the number of fluid phases
   * @return the number of phases
   */
  virtual localIndex
  NumFluidPhases() const = 0;

  /**
   * @brief getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a pointer to the controls
   */
  WellControls &
  GetWellControls( WellElementSubRegion const & subRegion );

  /**
   * @brief const getter for the well controls associated to this well subRegion
   * @param subRegion the well subRegion whose controls are requested
   * @return a pointer to the const controls
   */
  WellControls const &
  GetWellControls( WellElementSubRegion const & subRegion ) const;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual void
  RegisterDataOnMesh( Group * const meshBodies ) override;

  virtual void
  SetupDofs( DomainPartition const & domain,
             DofManager & dofManager ) const override;

  virtual void
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override;

  /**@}*/

  /**
   * @brief function to assemble the linear system matrix and rhs
   * @param time the time at the beginning of the step
   * @param dt the desired timestep
   * @param domain the domain partition
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleSystem(
    real64 const time,
    real64 const dt,
    DomainPartition & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) override;

  /**
   * @brief assembles the flux terms for all connections between well elements
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
    arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the volume balance terms for all well elements
   * @param time_n previous time value
   * @param dt time step
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  AssembleVolumeBalanceTerms(
    real64 const time_n,
    real64 const dt,
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief assembles the pressure relations at all connections between well elements except at the well head
   * @param domain the physical domain object
   * @param dofManager degree-of-freedom manager associated with the linear system
   * @param matrix the system matrix
   * @param rhs the system right-hand side vector
   */
  virtual void
  FormPressureRelations(
    DomainPartition const & domain,
    DofManager const & dofManager,
    CRSMatrixView< real64, globalIndex const > const & localMatrix,
    arrayView1d< real64 > const & localRhs ) = 0;

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param domain the domain containing the mesh and fields
   */
  virtual void
  UpdateStateAll( DomainPartition & domain );

  /**
   * @brief Recompute all dependent quantities from primary variables (including constitutive models)
   * @param well the well containing all the primary and dependent fields
   */
  virtual void
  UpdateState( WellElementSubRegion & subRegion,
               localIndex const targetIndex ) = 0;

  arrayView1d< string const > const &
  fluidModelNames() const
  {
    return m_fluidModelNames;
  }

  virtual std::vector< string >
  getConstitutiveRelations(
    string const & regionName ) const override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    // gravity term precomputed values
    static constexpr auto gravityCoefString =
      FlowSolverBase::viewKeyStruct::gravityCoefString;

    // misc inputs
    static constexpr auto fluidNamesString = "fluidNames";

  } viewKeysWellSolverBase;

  struct groupKeyStruct : SolverBase::groupKeyStruct
  {
  } groupKeysWellSolverBase;

private:
  /**
   * @brief This function generates various discretization information for later use.
   * @param domain the domain parition
   */
  void
  PrecomputeData( DomainPartition & domain );

protected:
  virtual void
  PostProcessInput() override;

  virtual void
  InitializePreSubGroups( Group * const rootGroup ) override;

  virtual void
  InitializePostInitialConditions_PreSubGroups(
    Group * const rootGroup ) override;

  /**
   * @brief Setup stored views into domain data for the current step
   */
  virtual void
  ResetViews( DomainPartition & domain );

  /**
   * @brief Initialize all the primary and secondary variables in all the wells
   * @param domain the domain containing the well manager to access individual wells
   */
  virtual void
  InitializeWells( DomainPartition & domain ) = 0;

  /**
   * @brief Check if the controls are viable; if not, switch the controls
   * @param domain the domain containing the well manager to access individual wells
   */
  //virtual void CheckWellControlSwitch( DomainPartition & domain ) = 0;

  /// name of the flow solver
  string m_flowSolverName;

  /// names of the fluid constitutive models
  array1d< string > m_fluidModelNames;

  /// the number of Degrees of Freedom per well element
  localIndex m_numDofPerWellElement;

  /// the number of Degrees of Freedom per reservoir element
  localIndex m_numDofPerResElement;

  /// views into reservoir constant data fields
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > m_resGravCoef;
};

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLSOLVERBASE_HPP_
