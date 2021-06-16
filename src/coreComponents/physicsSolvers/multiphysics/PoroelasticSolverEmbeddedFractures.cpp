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
 * @file PoroelasticSolverEmbeddedFractures.cpp
 *
 */


#include "PoroelasticSolverEmbeddedFractures.hpp"

#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/interfaces/dense/BlasLapackLA.hpp"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsEFEMKernelsHelper.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsEmbeddedFractures.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "SinglePhasePoroelasticKernel.hpp"
#include "SinglePhasePoroelasticEFEMKernel.hpp"
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

PoroelasticSolverEmbeddedFractures::PoroelasticSolverEmbeddedFractures( const std::string & name,
                                                                        Group * const parent ):
  PoroelasticSolver( name, parent ),
  m_fracturesSolverName()
{
  registerWrapper( viewKeyStruct::fracturesSolverNameString(), &m_fracturesSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fractures solver to use in the fractured poroelastic solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

PoroelasticSolverEmbeddedFractures::~PoroelasticSolverEmbeddedFractures()
{}

void PoroelasticSolverEmbeddedFractures::postProcessInput()
{
  PoroelasticSolver::postProcessInput();

  m_fracturesSolver  = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_fracturesSolverName );
}

void PoroelasticSolverEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  PoroelasticSolver::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
    {
      region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
      {
        subRegion.registerWrapper< array1d< real64 > >( EmbeddedSurfaceSubRegion::viewKeyStruct::dTraction_dPressureString() );
      } );
    } );
  } );
}

void PoroelasticSolverEmbeddedFractures::initializePostInitialConditionsPreSubGroups()
{
  updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
}

void PoroelasticSolverEmbeddedFractures::setupDofs( DomainPartition const & domain,
                                                    DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_fracturesSolver->setupDofs( domain, dofManager );
  m_flowSolver->setupDofs( domain, dofManager );

  // Add coupling between displacement and cell pressures
  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );

  // Add coupling between fracture pressure and displacement jump
  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = meshLevel.getElemManager();

  array1d< string > regions;
  elemManager.forElementRegions< SurfaceElementRegion >( [&]( SurfaceElementRegion const & region ) {
    regions.emplace_back( region.getName() );
  } );

  dofManager.addCoupling( FlowSolverBase::viewKeyStruct::pressureString(),
                          SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString(),
                          DofManager::Connector::Elem,
                          regions );
}

void PoroelasticSolverEmbeddedFractures::setupSystem( DomainPartition & domain,
                                                      DofManager & dofManager,
                                                      CRSMatrix< real64, globalIndex > & localMatrix,
                                                      array1d< real64 > & localRhs,
                                                      array1d< real64 > & localSolution,
                                                      bool const setSparsity )
{
  // Add missing couplings ( matrix pressure with displacement jump and jump - displacement )

  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );

  dofManager.setMesh( domain.getMeshBody( 0 ).getMeshLevel( 0 ) );
  setupDofs( domain, dofManager );
  dofManager.reorderByRank();

  // Set the sparsity pattern without the Kwu and Kuw blocks.
  SparsityPattern< globalIndex > patternDiag;
  dofManager.setSparsityPattern( patternDiag );

  // Get the original row lengths (diagonal blocks only)
  array1d< localIndex > rowLengths( patternDiag.numRows() );
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    rowLengths[localRow] = patternDiag.numNonZeros( localRow );
  }

  // Add the number of nonzeros induced by coupling jump-pm
  addCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(), patternDiag.numColumns(), rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow, cols, cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  addCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localRhs.resize( localMatrix.numRows() );
  localSolution.resize( localMatrix.numRows() );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );
}

void PoroelasticSolverEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                                 DofManager & dofManager,
                                                                 arrayView1d< localIndex > const & rowLengths ) const
{
  // 1. Add the number of nonzeros induced by coupling jump-displacement
  m_fracturesSolver->addCouplingNumNonzeros( domain, dofManager, rowLengths );

  // 2. Add the number of nonzeros induced by coupling jump - matrix pressure
  MeshLevel const & mesh                   = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
  {
    localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

    FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

    arrayView1d< globalIndex const > const &
    embeddedElementDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
    arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

    for( localIndex k=0; k<numEmbeddedElems; ++k )
    {
      // Get rock matrix element subregion
      CellElementSubRegion const & subRegion =
        elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
          getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

      arrayView1d< globalIndex const > const &
      pressureDofNumber = subRegion.getReference< globalIndex_array >( pressureDofKey );

      localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

      if( ghostRank[k] < 0 )
      {
        localIndex const localRow = LvArray::integerConversion< localIndex >( embeddedElementDofNumber[k] - rankOffset );
        GEOSX_ASSERT_GE( localRow, 0 );
        GEOSX_ASSERT_GE( rowLengths.size(), localRow + embeddedSurfaceSubRegion.numOfJumpEnrichments()  );

        for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
        {
          rowLengths[localRow + i] += 1;
        }

        localIndex const localPressureRow = LvArray::integerConversion< localIndex >( pressureDofNumber[cellElementIndex] - rankOffset );
        GEOSX_ASSERT_GE( localPressureRow, 0 );
        GEOSX_ASSERT_GE( rowLengths.size(), localPressureRow + embeddedSurfaceSubRegion.numOfJumpEnrichments() );

        rowLengths[ localPressureRow ] += embeddedSurfaceSubRegion.numOfJumpEnrichments();
      }
    }
  } );

  // 3. Add the number of nonzeros induced by coupling jump (aperture) - fracture pressure due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

      arrayView1d< globalIndex const > const &
      pressureDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( pressureDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = pressureDofNumber[sei[iconn][k0]];
        globalIndex const rowNumber = activeFlowDOF - rankOffset;

        if( rowNumber >= 0 && rowNumber < rowLengths.size() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the jump of the cell itself has already been added by the dofManager
            // so we only add the coupling with the jumps of the neighbours.
            if( k1 != k0 )
            {
              rowLengths[ rowNumber ] += 1;
            }
          }
        }
      }
    }
  } );

}

void PoroelasticSolverEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                     DofManager const & dofManager,
                                                                     SparsityPatternView< globalIndex > const & pattern ) const
{
  // 1. Add sparsity pattern induced by coupling jump-displacement
  m_fracturesSolver->addCouplingSparsityPattern( domain, dofManager, pattern );

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  // 2. Add the sparsity pattern induced by coupling jump - matrix pressure
  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion )
  {
    localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();

    FixedToManyElementRelation const & embeddedSurfacesToCells = embeddedSurfaceSubRegion.getToCellRelation();

    arrayView1d< globalIndex const > const &
    jumpDofNumber = embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
    arrayView1d< integer const > const & ghostRank = embeddedSurfaceSubRegion.ghostRank();

    for( localIndex k=0; k<numEmbeddedElems; ++k )
    {
      // Get rock matrix element subregion
      CellElementSubRegion const & subRegion =
        elemManager.getRegion( embeddedSurfacesToCells.m_toElementRegion[k][0] ).
          getSubRegion< CellElementSubRegion >( embeddedSurfacesToCells.m_toElementSubRegion[k][0] );

      arrayView1d< globalIndex const > const &
      pressureDofNumber = subRegion.getReference< globalIndex_array >( pressureDofKey );

      localIndex cellElementIndex = embeddedSurfacesToCells.m_toElementIndex[k][0];

      if( ghostRank[k] < 0 ) /// TODO is this really necessary?
      {
        localIndex const localJumpRow = LvArray::integerConversion< localIndex >( jumpDofNumber[k] - rankOffset );
        localIndex const localPressureRow = LvArray::integerConversion< localIndex >( pressureDofNumber[cellElementIndex] - rankOffset );

        for( localIndex i=0; i<embeddedSurfaceSubRegion.numOfJumpEnrichments(); ++i )
        {
          if( localJumpRow + i >= 0 && localJumpRow + i < pattern.numRows() )
            pattern.insertNonZero( localJumpRow + i, pressureDofNumber[cellElementIndex] );
          if( localPressureRow >= 0 && localPressureRow < pattern.numRows() )
            pattern.insertNonZero( localPressureRow, jumpDofNumber[k] + i );
        }
      }
    }
  } );

  // 3. Add the sparsity pattern induced by coupling jump (aperture) - fracture pressure due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( SurfaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename SurfaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename SurfaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

      EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion =
        elemManager.getRegion( seri[iconn][0] ).getSubRegion< EmbeddedSurfaceSubRegion >( sesri[iconn][0] );

      arrayView1d< globalIndex const > const &
      pressureDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( pressureDofKey );
      arrayView1d< globalIndex const > const &
      jumpDofNumber =  embeddedSurfaceSubRegion.getReference< globalIndex_array >( jumpDofKey );

      for( localIndex k0=0; k0<numFluxElems; ++k0 )
      {
        globalIndex const activeFlowDOF = pressureDofNumber[sei[iconn][k0]];
        globalIndex const rowIndex = activeFlowDOF - rankOffset;

        if( rowIndex >= 0 && rowIndex < pattern.numRows() )
        {
          for( localIndex k1=0; k1<numFluxElems; ++k1 )
          {
            // The coupling with the jump of the cell itself has already been added by the dofManager
            // so we only add the coupling with the jumps of the neighbours.
            if( k1 != k0 )
            {
              globalIndex const colIndex = jumpDofNumber[sei[iconn][k1]];
              pattern.insertNonZero( rowIndex, colIndex );
            }
          }
        }
      }
    }
  } );

}

void PoroelasticSolverEmbeddedFractures::assembleSystem( real64 const time_n,
                                                         real64 const dt,
                                                         DomainPartition & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  //updateState( domain );

  MeshLevel & mesh = domain.getMeshBodies().getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  NodeManager const & nodeManager = mesh.getNodeManager();
  ElementRegionManager const & elemManager = mesh.getElemManager();
  SurfaceElementRegion const & region = elemManager.getRegion< SurfaceElementRegion >( m_fracturesSolver->getFractureRegionName() );
  EmbeddedSurfaceSubRegion const & subRegion = region.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  string const dofKey = dofManager.getKey( dataRepository::keys::TotalDisplacement );
  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );

  arrayView1d< globalIndex const > const & dispDofNumber = nodeManager.getReference< globalIndex_array >( dofKey );
  arrayView1d< globalIndex const > const & jumpDofNumber = subRegion.getReference< globalIndex_array >( jumpDofKey );

  string const pDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  real64 const gravityVectorData[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );


  // 1. Cell-based contributions of standard poroelasticity
  PoroelasticKernels::SinglePhaseKernelFactory kernelFactory( dispDofNumber,
                                                              pDofKey,
                                                              dofManager.rankOffset(),
                                                              localMatrix,
                                                              localRhs,
                                                              gravityVectorData,
                                                              m_flowSolver->fluidModelNames() );

  m_solidSolver->getMaxForce() =
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::SolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            m_solidSolver->getDiscretizationName(),
                                                            m_solidSolver->solidMaterialNames(),
                                                            kernelFactory );

  // 2.  Add EFEM poroelastic contribution
  PoroelasticEFEMKernels::SinglePhaseKernelFactory EFEMkernelFactory( subRegion,
                                                                      dispDofNumber,
                                                                      jumpDofNumber,
                                                                      pDofKey,
                                                                      dofManager.rankOffset(),
                                                                      localMatrix,
                                                                      localRhs,
                                                                      gravityVectorData,
                                                                      m_flowSolver->fluidModelNames() );

  real64 maxTraction =
    finiteElement::
      regionBasedKernelApplication< parallelDevicePolicy< 32 >,
                                    constitutive::SolidBase,
                                    CellElementSubRegion >( mesh,
                                                            targetRegionNames(),
                                                            m_solidSolver->getDiscretizationName(),
                                                            m_solidSolver->solidMaterialNames(),
                                                            EFEMkernelFactory );

  GEOSX_UNUSED_VAR( maxTraction );

  // 3. Assemble poroelastic fluxes and all derivatives
  m_flowSolver->assemblePoroelasticFluxTerms( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs,
                                              jumpDofKey );

}



void PoroelasticSolverEmbeddedFractures::applyBoundaryConditions( real64 const time_n,
                                                                  real64 const dt,
                                                                  DomainPartition & domain,
                                                                  DofManager const & dofManager,
                                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                  arrayView1d< real64 > const & localRhs )
{
  m_fracturesSolver->applyBoundaryConditions( time_n, dt,
                                              domain,
                                              dofManager,
                                              localMatrix,
                                              localRhs );

  m_flowSolver->applyBoundaryConditions( time_n, dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

void PoroelasticSolverEmbeddedFractures::implicitStepSetup( real64 const & time_n,
                                                            real64 const & dt,
                                                            DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_fracturesSolver->implicitStepSetup( time_n, dt, domain );
}

void PoroelasticSolverEmbeddedFractures::implicitStepComplete( real64 const & time_n,
                                                               real64 const & dt,
                                                               DomainPartition & domain )
{
  m_fracturesSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}

void PoroelasticSolverEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_fracturesSolver->resetStateToBeginningOfStep( domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&] ( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const & oldTotalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString() );
    arrayView1d< real64 > const & totalMeanStress =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMeanStressString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      totalMeanStress[ei] = oldTotalMeanStress[ei];
    } );
  } );
}

real64 PoroelasticSolverEmbeddedFractures::solverStep( real64 const & time_n,
                                                       real64 const & dt,
                                                       int const cycleNumber,
                                                       DomainPartition & domain )
{
  real64 dtReturn = dt;

  /// TODO
  // for (integer outerIter = 0; outerIter < m_maxOuterIter; outerIter++)
  {
    setupSystem( domain,
                 m_dofManager,
                 m_localMatrix,
                 m_localRhs,
                 m_localSolution );

    implicitStepSetup( time_n, dt, domain );

    // Given a fracture state we solve the system
    dtReturn = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    implicitStepComplete( time_n, dtReturn, domain );

    // check the fracture state
//    bool fractureStateUnchaged = true; // TODO
//    if ( fractureStateUnChanged )
//    {
//      break;
//    }
  }

  return dtReturn;
}

real64 PoroelasticSolverEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
                                                                  DofManager const & dofManager,
                                                                  arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm = m_fracturesSolver->calculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm = m_flowSolver->calculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank==0 )
  {
    char output[200] = {0};
    sprintf( output, "    ( Rsolid, Rfluid ) = ( %4.2e, %4.2e )", momementumResidualNorm, massResidualNorm );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm + massResidualNorm * massResidualNorm );
}

void PoroelasticSolverEmbeddedFractures::applySystemSolution( DofManager const & dofManager,
                                                              arrayView1d< real64 const > const & localSolution,
                                                              real64 const scalingFactor,
                                                              DomainPartition & domain )
{
  // update displacement and jump
  m_fracturesSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );
}

void PoroelasticSolverEmbeddedFractures::updateState( DomainPartition & domain )
{
  /// 1. update the reservoir
  PoroelasticSolver::updateState( domain );

  /// 2. update the fractures

  m_fracturesSolver->updateState( domain );

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  //ElementRegionManager & elemManager = meshLevel.getElemManager();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();

  ContactRelationBase const &
  contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_fracturesSolver->getContactRelationName() );

  this->template forTargetSubRegions< EmbeddedSurfaceSubRegion >( meshLevel, [&] ( localIndex const targetIndex,
                                                                                   auto & subRegion )
  {
    arrayView2d< real64 const > const dispJump =
      subRegion.template getReference< array2d< real64 > >( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );

    arrayView1d< real64 > const aperture = subRegion.getElementAperture();

    arrayView1d< real64 > const effectiveAperture =
      subRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString() );

    arrayView1d< real64 const > const volume = subRegion.getElementVolume();

    arrayView1d< real64 > const deltaVolume =
      subRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString() );
    arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();

    arrayView2d< real64 > const & fractureTraction =
      subRegion.template getReference< array2d< real64 > >( SolidMechanicsEmbeddedFractures::viewKeyStruct::fractureTractionString() );

    arrayView1d< real64 >  const & dTdpf =
      subRegion.template getReference< array1d< real64 > >( viewKeyStruct::dTraction_dPressureString() );

    arrayView1d< real64 const > const & pressure =
      subRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );

    arrayView1d< real64 const > const & deltaPressure =
      subRegion.template getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );

    forAll< serialPolicy >( subRegion.size(), [=, &contactRelation] ( localIndex const k )
    {
      // update aperture to be equal to the normal displacement jump
      aperture[k] = dispJump[k][0];   // the first component of the jump is the normal one.

      effectiveAperture[k] = contactRelation.effectiveAperture( aperture[k] );

      deltaVolume[k] = effectiveAperture[k] * area[k] - volume[k];

      // traction on the fracture to include the pressure contribution
      contactRelation.addPressureToTraction( pressure[k] + deltaPressure[k], fractureTraction[k] );

      bool open = aperture[k] >= 0 ? true : false;
      contactRelation.dTraction_dPressure( dTdpf[k], open );
    } );

    // update fracture's permeability and porosity
    m_flowSolver->updateSolidFlowProperties( subRegion, targetIndex );
    // update fluid model
    m_flowSolver->updateFluidState( subRegion, targetIndex );

  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, PoroelasticSolverEmbeddedFractures, std::string const &, Group * const )

} /* namespace geosx */
