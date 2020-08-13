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

/*
 * SolidMechanicsEmbeddedFractures.hpp
 *
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_

#include "physicsSolvers/SolverBase.hpp"

namespace geosx
{
using namespace constitutive;

class SolidMechanicsLagrangianFEM;

class SolidMechanicsEmbeddedFractures : public SolverBase
{
public:
  SolidMechanicsEmbeddedFractures( const std::string & name, Group * const parent );

  ~SolidMechanicsEmbeddedFractures() override;

  /**
   * @brief name of the node manager in the object catalog
   * @return string that contains the catalog name to generate a new NodeManager object through the object catalog.
   */
  static string
  CatalogName()
  {
    return "SolidMechanicsEmbeddedFractures";
  }

  virtual void
  RegisterDataOnMesh(
    dataRepository::Group * const MeshBodies ) override final;

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
  ImplicitStepSetup( real64 const & time_n,
                     real64 const & dt,
                     DomainPartition & domain ) override final;

  virtual void
  ImplicitStepComplete( real64 const & time_n,
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
  ApplySystemSolution( DofManager const & dofManager,
                       arrayView1d< real64 const > const & localSolution,
                       real64 const scalingFactor,
                       DomainPartition & domain ) override;

  virtual void
  ResetStateToBeginningOfStep( DomainPartition & domain ) override final;

  virtual real64
  SolverStep( real64 const & time_n,
              real64 const & dt,
              int const cycleNumber,
              DomainPartition & domain ) override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static auto solidSolverNameString = "solidSolverName";

    constexpr static auto contactRelationNameString = "contactRelationName";

    constexpr static auto dispJumpString = "displacementJump";

    constexpr static auto deltaDispJumpString = "deltaDisplacementJump";

  } SolidMechanicsEmbeddedFracturesViewKeys;

protected:
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
  void
  AddCouplingSparsityPattern(
    DomainPartition const & domain,
    DofManager const & dofManager,
    SparsityPatternView< globalIndex > const & pattern ) const;

  /*
   * @brief Assemble Equilibrium operator
   * @param eqMatrix Equilibrium operator
   * @param embeddedSurfaceSubRegion subRegion
   * @param k cell index
   * @param hInv scaling coefficient
   */
  void
  AssembleEquilibriumOperator(
    array2d< real64 > & eqMatrix,
    EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
    const localIndex k,
    const real64 hInv );
  /*
   * @brief Assemble Compatibility operator
   * @param compMatrix
   * @param embeddedSurfaceSubRegion
   * @param k cell index
   * @param q quadrature point index
   * @param elemsToNodes element to node map
   * @param nodesCoord nodes coordinates
   * @param embeddedSurfaceToCell embedded surface to cell maps
   * @param numNodesPerElement number of nodes per element
   * @param dNdX shape functions derivatives
   */
  void
  AssembleCompatibilityOperator(
    array2d< real64 > & compMatrix,
    EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
    localIndex const k,
    localIndex const q,
    CellBlock::NodeMapType const & elemsToNodes,
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord,
    arrayView1d< localIndex const > const & embeddedSurfaceToCell,
    localIndex const numNodesPerElement,
    arrayView4d< real64 const > const & dNdX );

  /*
   * @brief Assemble Compatibility operator
   * @param strainMatrix strain matrix (B)
   * @param elIndex element index
   * @param q quadrature point index
   * @param numNodesPerElement number of nodes per element
   * @param dNdX shape functions derivatives
   */
  void
  AssembleStrainOperator( array2d< real64 > & strainMatrix,
                          localIndex const elIndex,
                          localIndex const q,
                          localIndex const numNodesPerElement,
                          arrayView4d< real64 const > const & dNdX );
  /*
   * @brief Computes traction and derivative on each fracture segment.
   * @param constitutiveManager constant pointer to the constitutive mamanger
   * @param dispJump displacement jump
   * @param tractionVector traction vector
   * @param dTdw Derivative of the traction w.r.t. the jump.
   */
  void
  ComputeTraction( ConstitutiveManager const * const constitutiveManager,
                   array1d< real64 > const & dispJump,
                   array1d< real64 > & tractionVector,
                   array2d< real64 > & dTdw );

private:
  /// Solid mechanics solver name
  string m_solidSolverName;

  /// pointer to the solid mechanics solver
  SolidMechanicsLagrangianFEM * m_solidSolver;

  /// contact relation name string
  string m_contactRelationName;
};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_SOLIDMECHANICS_SOLIDMECHANICSEMBEDDEDFRACTURES_HPP_ */
