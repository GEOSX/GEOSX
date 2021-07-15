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


#include "SinglePhasePoromechanicsSolverEmbeddedFractures.hpp"

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
#include "common/GEOS_RAJA_Interface.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

SinglePhasePoromechanicsSolverEmbeddedFractures::SinglePhasePoromechanicsSolverEmbeddedFractures( const std::string & name,
                                                                                                  Group * const parent ):
  SinglePhasePoromechanicsSolver( name, parent ),
  m_fracturesSolverName()
{
  registerWrapper( viewKeyStruct::fracturesSolverNameString(), &m_fracturesSolverName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of the fractures solver to use in the fractured poroelastic solver" );

  this->getWrapper< string >( viewKeyStruct::discretizationString() ).
    setInputFlag( InputFlags::FALSE );
}

SinglePhasePoromechanicsSolverEmbeddedFractures::~SinglePhasePoromechanicsSolverEmbeddedFractures()
{}

void SinglePhasePoromechanicsSolverEmbeddedFractures::postProcessInput()
{
  SinglePhasePoromechanicsSolver::postProcessInput();

  m_fracturesSolver  = &this->getParent().getGroup< SolidMechanicsEmbeddedFractures >( m_fracturesSolverName );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::registerDataOnMesh( dataRepository::Group & meshBodies )
{
  SinglePhasePoromechanicsSolver::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    ElementRegionManager & elemManager = meshLevel.getElemManager();
    elemManager.forElementRegions< SurfaceElementRegion >( [&] ( SurfaceElementRegion & region )
    {
      region.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
      {
        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dTraction_dPressureString() );
      } );
    } );
  } );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::initializePostInitialConditionsPreSubGroups()
{
  updateState( this->getGroupByPath< DomainPartition >( "/Problem/domain" ) );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::setupDofs( DomainPartition const & domain,
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

void SinglePhasePoromechanicsSolverEmbeddedFractures::setupSystem( DomainPartition & domain,
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

  // sets up the object that contains the deravative of the fluxes w.r.t. the aperture ( i.e., the normal jump )
  m_flowSolver->setUpDflux_dApertureMatrix( domain, dofManager, localMatrix );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::addCouplingNumNonzeros( DomainPartition & domain,
                                                                              DofManager & dofManager,
                                                                              arrayView1d< localIndex > const & rowLengths ) const
{
  // Add the nonzeros from coupling jump-displacement
  m_fracturesSolver->addCouplingNumNonzeros( domain, dofManager, rowLengths );

  // Add the number of nonzeros induced by coupling jump-displacement
  MeshLevel const & mesh                   = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  // Coupling jump - matrix pressure
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

  // Coupling jump (aperture) - fracture pressure due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

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

void SinglePhasePoromechanicsSolverEmbeddedFractures::addCouplingSparsityPattern( DomainPartition const & domain,
                                                                                  DofManager const & dofManager,
                                                                                  SparsityPatternView< globalIndex > const & pattern ) const
{
  m_fracturesSolver->addCouplingSparsityPattern( domain, dofManager, pattern );

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  // Coupling jump - matrix pressure
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

  // Coupling fracture pressure - jump (aperture) due to flux term
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_flowSolver->getDiscretization() );

  fluxApprox.forStencils< FaceElementStencil >( mesh, [&]( FaceElementStencil const & stencil )
  {
    for( localIndex iconn=0; iconn<stencil.size(); ++iconn )
    {
      localIndex const numFluxElems = stencil.stencilSize( iconn );
      typename FaceElementStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
      typename FaceElementStencil::IndexContainerViewConstType const & sei = stencil.getElementIndices();

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

void SinglePhasePoromechanicsSolverEmbeddedFractures::assembleSystem( real64 const time_n,
                                                                      real64 const dt,
                                                                      DomainPartition & domain,
                                                                      DofManager const & dofManager,
                                                                      CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                      arrayView1d< real64 > const & localRhs )
{

  updateState( domain );

  // assemble Kuu, Kuw, Kww, Kwu
  m_fracturesSolver->assembleSystem( time_n, dt,
                                     domain,
                                     dofManager,
                                     localMatrix,
                                     localRhs );

  // assemble flow matrices
  m_flowSolver->assembleSystem( time_n, dt,
                                domain,
                                dofManager,
                                localMatrix,
                                localRhs );

  // assemble mechanics-flow coupling blocks
  assembleCouplingTerms( domain,
                         dofManager,
                         localMatrix,
                         localRhs );

}

void SinglePhasePoromechanicsSolverEmbeddedFractures::assembleCouplingTerms( DomainPartition const & domain,
                                                                             DofManager const & dofManager,
                                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                                             arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Rock matrix poroelatic coupling
  SinglePhasePoromechanicsSolver::assembleCouplingTerms( domain, dofManager, localMatrix, localRhs );

  assembleFractureFlowResidualWrtJump( domain, dofManager, localMatrix, localRhs );

  assembleTractionBalanceResidualWrtPressure( domain, dofManager, localMatrix, localRhs );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::
  assembleTractionBalanceResidualWrtPressure( DomainPartition const & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{

  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const pressureDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString() );
  string const jumpDofKey     = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );

  // Get fracture's views
  ElementRegionManager const & elemManager = mesh.getElemManager();

  SurfaceElementRegion const & fractureRegion =
    elemManager.getRegion< SurfaceElementRegion >( m_fracturesSolver->getFractureRegionName() );

  EmbeddedSurfaceSubRegion const & embeddedSurfSubRegion = fractureRegion.getSubRegion< EmbeddedSurfaceSubRegion >( 0 );

  arrayView1d< globalIndex const > const jumpDofNumber = embeddedSurfSubRegion.getReference< globalIndex_array >( jumpDofKey );

  arrayView2d< real64 const > const nVec  = embeddedSurfSubRegion.getNormalVector();

  arrayView2d< real64 const > const tVec1 =  embeddedSurfSubRegion.getTangentVector1();

  arrayView2d< real64 const > const tVec2 =  embeddedSurfSubRegion.getTangentVector2();

  arrayView1d< real64 const > const fractureSurfaceArea = embeddedSurfSubRegion.getElementArea();

  arrayView1d< globalIndex const > const & fracturePresDofNumber =
    embeddedSurfSubRegion.getReference< array1d< globalIndex > >( pressureDofKey );

  arrayView1d< real64 const > const & dTraction_dPressure =
    embeddedSurfSubRegion.getReference< array1d< real64 > >( viewKeyStruct::dTraction_dPressureString() );


  // begin subregion loop
  forTargetSubRegionsComplete< CellElementSubRegion >( mesh, [&]( localIndex const,
                                                                  localIndex const,
                                                                  localIndex const,
                                                                  ElementRegionBase const & region,
                                                                  CellElementSubRegion const & elementSubRegion )
  {
    string const & solidName = m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex( region.getName() )];
    SolidBase const & solid = getConstitutiveModel< SolidBase >( elementSubRegion, solidName );

    real64 const biotCoefficient = solid.getReference< real64 >( "BiotCoefficient" );

    // Get FEM information for integration
    finiteElement::FiniteElementBase const &
    fe = elementSubRegion.getReference< finiteElement::FiniteElementBase >( m_solidSolver->getDiscretizationName() );
    localIndex const numQuadraturePoints = fe.getNumQuadraturePoints();

    // Get cell element subregion views

    arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

    SortedArrayView< localIndex const > const fracturedElements = elementSubRegion.fracturedElementsList();

    ArrayOfArraysView< localIndex const > const cellsToEmbeddedSurfaces = elementSubRegion.embeddedSurfacesList().toViewConst();

    arrayView1d< real64 const > const & cellVolume = elementSubRegion.getElementVolume();

    arrayView1d< real64 const > const & matrixPres =
      elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );

    arrayView1d< real64 const > const & matrixDeltaPres =
      elementSubRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );

    arrayView1d< globalIndex const > const matrixPresDofNumber =
      elementSubRegion.getReference< array1d< globalIndex > >( pressureDofKey );

    forAll< parallelDevicePolicy< 32 > >( fracturedElements.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      localIndex cellIndex     = fracturedElements[ei];
      localIndex embSurfIndex = cellsToEmbeddedSurfaces[cellIndex][0];

      // Initialize local variables
      real64 jumpEqnRowIndices[3], eqMatrix[3][6], Rfrac[3], Kwpm_gauss[3], Kwpm_elem[3];

      real64 const pm = matrixPres[cellIndex] + matrixDeltaPres[cellIndex];

      // Fill in equilibrium operator
      LvArray::tensorOps::fill< 3, 6 >( eqMatrix, 0 );
      real64 hInv = fractureSurfaceArea[embSurfIndex] / cellVolume[cellIndex];  // AreaFrac / cellVolume
      SolidMechanicsEFEMKernelsHelper::assembleEquilibriumOperator( eqMatrix,
                                                                    nVec[embSurfIndex],
                                                                    tVec1[embSurfIndex],
                                                                    tVec2[embSurfIndex],
                                                                    hInv );

      // 1. Assembly of element matrices
      LvArray::tensorOps::fill< 3 >( Kwpm_elem, 0 );
      LvArray::tensorOps::fill< 3 >( Kwpm_gauss, 0 );

      for( int i=0; i < 3; ++i )
      {
        Kwpm_gauss[0] += eqMatrix[0][i];
        Kwpm_gauss[1] += eqMatrix[1][i];
        Kwpm_gauss[2] += eqMatrix[2][i];
      }

      // dTdpf
      // signs for the mechanics are flipped so we add a minus.
      real64 dTdpf = -dTraction_dPressure[embSurfIndex] * fractureSurfaceArea[embSurfIndex];

      for( integer q=0; q<numQuadraturePoints; ++q )
      {
        const real64 detJq = detJ[cellIndex][q];

        // No neg coz the effective stress is total stress - porePressure
        // and all signs are flipped here.
        Kwpm_elem[0] += Kwpm_gauss[0] * biotCoefficient * detJq;
        Kwpm_elem[1] += Kwpm_gauss[1] * biotCoefficient * detJq;
        Kwpm_elem[2] += Kwpm_gauss[2] * biotCoefficient * detJq;
      }

      Rfrac[0] = Kwpm_elem[0] * pm;
      Rfrac[1] = Kwpm_elem[1] * pm;
      Rfrac[2] = Kwpm_elem[2] * pm;

      /// Add derivatives and residual contribution for porelastic case.

      // Row and column indexes.
      globalIndex matrixPressureColIndex = matrixPresDofNumber[cellIndex];
      // Dof number of jump enrichment
      for( int i= 0; i < 3; i++ )
      {
        jumpEqnRowIndices[i] = jumpDofNumber[embSurfIndex] + i - rankOffset;
      }

      for( localIndex i=0; i < 3; ++i )
      {
        if( jumpEqnRowIndices[i] >= 0 && jumpEqnRowIndices[i] < localMatrix.numRows() )
        {
          RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[jumpEqnRowIndices[i]], Rfrac[i] );

          // fill in matrix
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( jumpEqnRowIndices[i],
                                                                            &matrixPressureColIndex,
                                                                            &Kwpm_elem[i],
                                                                            1 );
        }
      }
      if( jumpEqnRowIndices[0] >= 0 && jumpEqnRowIndices[0] < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( jumpEqnRowIndices[0],
                                                                          &fracturePresDofNumber[embSurfIndex],
                                                                          &dTdpf,
                                                                          1 );
      }


    } );
  } );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::
  assembleFractureFlowResidualWrtJump( DomainPartition const & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( localRhs );

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const presDofKey = dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString ());
  string const jumpDofKey = dofManager.getKey( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  CRSMatrixView< real64 const, localIndex const > const &
  dFluxResidual_dAperture = m_flowSolver->getDerivativeFluxResidual_dAperture().toViewConst();

  forTargetSubRegionsComplete< EmbeddedSurfaceSubRegion >( mesh,
                                                           [&]( localIndex const,
                                                                localIndex const,
                                                                localIndex const,
                                                                ElementRegionBase const & region,
                                                                EmbeddedSurfaceSubRegion const & subRegion )
  {
    string const & fluidName = m_flowSolver->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
    SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, fluidName );

    arrayView1d< globalIndex const > const presDofNumber = subRegion.getReference< array1d< globalIndex > >( presDofKey );
    arrayView1d< globalIndex const > const jumpDofNumber = subRegion.getReference< array1d< globalIndex > >( jumpDofKey );

    arrayView2d< real64 const > const dens = fluid.density();

    arrayView1d< real64 const > const area = subRegion.getElementArea();

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_DEVICE ( localIndex ei )
    {
      globalIndex const elemDOF = presDofNumber[ei];
      // row number associated to the pressure dof
      globalIndex rowNumber = elemDOF - rankOffset;

      real64 const dAccumulationResidualdAperture = dens[ei][0] * area[ei];
      globalIndex const jumpDOF = jumpDofNumber[ei];

      if( rowNumber >= 0  && rowNumber < localMatrix.numRows() )
      {
        localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                          &jumpDOF,
                                                                          &dAccumulationResidualdAperture,
                                                                          1 );
      }

      // flux derivative
      localIndex const numColumns = dFluxResidual_dAperture.numNonZeros( ei );
      arraySlice1d< localIndex const > const & columns = dFluxResidual_dAperture.getColumns( ei );
      arraySlice1d< real64 const > const & values = dFluxResidual_dAperture.getEntries( ei );

      for( localIndex kfe2 = 0; kfe2 < numColumns; ++kfe2 )
      {
        real64 dRdAper = values[kfe2];
        localIndex const ei2 = columns[kfe2];

        globalIndex const jumpDOF2 = jumpDofNumber[ei2];

        if( rowNumber >= 0 && rowNumber < localMatrix.numRows() )
        {
          localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >( rowNumber,
                                                                            &jumpDOF2,
                                                                            &dRdAper,
                                                                            1 );
        }
      }
    } );
  } );

}

void SinglePhasePoromechanicsSolverEmbeddedFractures::applyBoundaryConditions( real64 const time_n,
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

void SinglePhasePoromechanicsSolverEmbeddedFractures::implicitStepSetup( real64 const & time_n,
                                                                         real64 const & dt,
                                                                         DomainPartition & domain )
{
  m_flowSolver->implicitStepSetup( time_n, dt, domain );
  m_fracturesSolver->implicitStepSetup( time_n, dt, domain );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::implicitStepComplete( real64 const & time_n,
                                                                            real64 const & dt,
                                                                            DomainPartition & domain )
{
  m_fracturesSolver->implicitStepComplete( time_n, dt, domain );
  m_flowSolver->implicitStepComplete( time_n, dt, domain );
}



void SinglePhasePoromechanicsSolverEmbeddedFractures::resetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->resetStateToBeginningOfStep( domain );
  m_fracturesSolver->resetStateToBeginningOfStep( domain );
}

real64 SinglePhasePoromechanicsSolverEmbeddedFractures::solverStep( real64 const & time_n,
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

real64 SinglePhasePoromechanicsSolverEmbeddedFractures::calculateResidualNorm( DomainPartition const & domain,
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

void SinglePhasePoromechanicsSolverEmbeddedFractures::applySystemSolution( DofManager const & dofManager,
                                                                           arrayView1d< real64 const > const & localSolution,
                                                                           real64 const scalingFactor,
                                                                           DomainPartition & domain )
{
  // update displacement and jump
  m_fracturesSolver->applySystemSolution( dofManager, localSolution, scalingFactor, domain );
  // update pressure field
  m_flowSolver->applySystemSolution( dofManager, localSolution, -scalingFactor, domain );

  updateState( domain );
}

void SinglePhasePoromechanicsSolverEmbeddedFractures::updateState( DomainPartition & domain )
{
  // update aperture to be equal to the normal displacement jump and traction on the fracture to include the pressure contribution
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = meshLevel.getElemManager();

  ConstitutiveManager const & constitutiveManager = domain.getConstitutiveManager();

  ContactRelationBase const &
  contactRelation = constitutiveManager.getGroup< ContactRelationBase >( m_fracturesSolver->getContactRelationName() );

  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >( [&]( EmbeddedSurfaceSubRegion & subRegion )
  {
    arrayView2d< real64 const > const dispJump =
      subRegion.getReference< array2d< real64 > >( SolidMechanicsEmbeddedFractures::viewKeyStruct::dispJumpString() );

    arrayView1d< real64 > const aperture = subRegion.getElementAperture();

    arrayView1d< real64 > const effectiveAperture =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::effectiveApertureString() );

    arrayView1d< real64 const > const volume = subRegion.getElementVolume();

    arrayView1d< real64 > const deltaVolume =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaVolumeString() );
    arrayView1d< real64 const > const area = subRegion.getElementArea().toViewConst();

    arrayView2d< real64 > const & fractureTraction =
      subRegion.getReference< array2d< real64 > >( SolidMechanicsEmbeddedFractures::viewKeyStruct::fractureTractionString() );

    arrayView1d< real64 >  const & dTdpf =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTraction_dPressureString() );

    arrayView1d< real64 const > const & pressure =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::pressureString() );

    arrayView1d< real64 const > const & deltaPressure =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::deltaPressureString() );

    forAll< serialPolicy >( subRegion.size(), [=, &contactRelation] ( localIndex const k )
    {
      aperture[k] = dispJump[k][0];   // the first component of the jump is the normal one.

      effectiveAperture[k] = contactRelation.effectiveAperture( aperture[k] );

      deltaVolume[k] = effectiveAperture[k] * area[k] - volume[k];

      contactRelation.addPressureToTraction( pressure[k] + deltaPressure[k], fractureTraction[k] );

      bool open = aperture[k] >= 0 ? true : false;
      contactRelation.dTraction_dPressure( dTdpf[k], open );
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhasePoromechanicsSolverEmbeddedFractures, std::string const &, Group * const )

} /* namespace geosx */
