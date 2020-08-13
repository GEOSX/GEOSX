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
 * SolidMechanicsEmbeddedFractures.cpp
 */

#include "SolidMechanicsEmbeddedFractures.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/solid/LinearElasticIsotropic.hpp"
#include "finiteElement/elementFormulations/FiniteElementBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/EmbeddedSurfaceRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"
#include "linearAlgebra/interfaces/BlasLapackLA.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

SolidMechanicsEmbeddedFractures::SolidMechanicsEmbeddedFractures(
  const std::string & name,
  Group * const parent ) :
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_solidSolver( nullptr )
{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Name of the solid mechanics solver in the rock matrix" );

  registerWrapper( viewKeyStruct::contactRelationNameString, &m_contactRelationName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of contact relation to enforce constraints on fracture boundary." );
}

SolidMechanicsEmbeddedFractures::~SolidMechanicsEmbeddedFractures()
{
  // TODO Auto-generated destructor stub
}

void
SolidMechanicsEmbeddedFractures::RegisterDataOnMesh(
  dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel =
      Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();
    {
      elemManager->forElementRegions< EmbeddedSurfaceRegion >(
        [&]( EmbeddedSurfaceRegion & region ) {
          region.forElementSubRegions< EmbeddedSurfaceSubRegion >(
            [&]( EmbeddedSurfaceSubRegion & subRegion ) {
              //subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::dispJumpString
              // )->setPlotLevel(PlotLevel::LEVEL_0);
              // subRegion->registerWrapper< array1d<real64> >( viewKeyStruct::deltaDispJumpString );
              subRegion
                .registerWrapper< array1d< R1Tensor > >( viewKeyStruct::dispJumpString )
                ->setPlotLevel( PlotLevel::LEVEL_0 );
              subRegion.registerWrapper< array1d< R1Tensor > >(
                viewKeyStruct::deltaDispJumpString );
            } );
        } );
    }
  }
}

void
SolidMechanicsEmbeddedFractures::ResetStateToBeginningOfStep(
  DomainPartition & domain )
{
  m_solidSolver->ResetStateToBeginningOfStep( domain );
}

void
SolidMechanicsEmbeddedFractures::ImplicitStepSetup( real64 const & time_n,
                                                    real64 const & dt,
                                                    DomainPartition & domain )
{
  m_solidSolver =
    this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  m_solidSolver->ImplicitStepSetup( time_n, dt, domain );
}

void
SolidMechanicsEmbeddedFractures::ImplicitStepComplete( real64 const & time_n,
                                                       real64 const & dt,
                                                       DomainPartition & domain )
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
}

real64
SolidMechanicsEmbeddedFractures::SolverStep( real64 const & time_n,
                                             real64 const & dt,
                                             int const cycleNumber,
                                             DomainPartition & domain )
{
  real64 dtReturn = dt;

  ImplicitStepSetup( time_n, dt, domain );

  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

  // currently the only method is implicit time integration
  dtReturn = this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // m_solidSolver->updateStress( domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dtReturn, domain );

  return dtReturn;
}

void
SolidMechanicsEmbeddedFractures::SetupDofs( DomainPartition const & domain,
                                            DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );

  MeshLevel const & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *meshLevel.getElemManager();

  array1d< string > regions;
  elemManager.forElementRegions< EmbeddedSurfaceRegion >(
    [&]( EmbeddedSurfaceRegion const & region ) {
      regions.emplace_back( region.getName() );
    } );

  dofManager.addField( viewKeyStruct::dispJumpString,
                       DofManager::Location::Elem,
                       3,
                       regions );

  dofManager.addCoupling( viewKeyStruct::dispJumpString,
                          viewKeyStruct::dispJumpString,
                          DofManager::Connector::Elem,
                          regions );
}

void
SolidMechanicsEmbeddedFractures::SetupSystem(
  DomainPartition & domain,
  DofManager & dofManager,
  CRSMatrix< real64, globalIndex > & localMatrix,
  array1d< real64 > & localRhs,
  array1d< real64 > & localSolution,
  bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  GEOSX_UNUSED_VAR( setSparsity );
  // By not calling dofManager.reorderByRank(), we keep separate dof numbering for each field,
  // which allows constructing separate sparsity patterns for off-diagonal blocks of the matrix.
  // Once the solver moves to monolithic matrix, we can remove this method and just use SolverBase::SetupSystem.
  // setup coupled DofManager
  dofManager.setMesh( domain, 0, 0 );
  SetupDofs( domain, dofManager );
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

  // Add the number of nonzeros induced by coupling
  AddCouplingNumNonzeros( domain, dofManager, rowLengths.toView() );

  // Create a new pattern with enough capacity for coupled matrix
  SparsityPattern< globalIndex > pattern;
  pattern.resizeFromRowCapacities< parallelHostPolicy >( patternDiag.numRows(),
                                                         patternDiag.numColumns(),
                                                         rowLengths.data() );

  // Copy the original nonzeros
  for( localIndex localRow = 0; localRow < patternDiag.numRows(); ++localRow )
  {
    globalIndex const * cols = patternDiag.getColumns( localRow ).dataIfContiguous();
    pattern.insertNonZeros( localRow,
                            cols,
                            cols + patternDiag.numNonZeros( localRow ) );
  }

  // Add the nonzeros from coupling
  AddCouplingSparsityPattern( domain, dofManager, pattern.toView() );

  // Finally, steal the pattern into a CRS matrix
  localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  localRhs.resize( localMatrix.numRows() );
  localSolution.resize( localMatrix.numRows() );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );

  ///
}

void
SolidMechanicsEmbeddedFractures::AssembleSystem(
  real64 const time,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->AssembleSystem( time, dt, domain, dofManager, localMatrix, localRhs );

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  ConstitutiveManager const * const constitutiveManager =
    domain.getConstitutiveManager();

  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & disp =
    nodeManager.totalDisplacement();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & dDisp =
    nodeManager.incrementalDisplacement();
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord =
    nodeManager.referencePosition();

  r1_array const uhattilde;

  string const dofKey = dofManager.getKey( keys::TotalDisplacement );
  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  arrayView1d< globalIndex const > const & globalDofNumber =
    nodeManager.getReference< globalIndex_array >( dofKey );

  //  ElementRegionManager::ConstitutiveRelationAccessor< ConstitutiveBase const >
  //  constitutiveRelations = elemManager.ConstructFullConstitutiveAccessor< ConstitutiveBase const >( constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< real64 > const density =
    elemManager.ConstructFullMaterialViewAccessor< real64 >( "density0",
                                                             constitutiveManager );
  globalIndex const rankOffset = dofManager.rankOffset();

  constexpr int dim = 3;
  static constexpr int maxNumUdof = dim * 8;  // this is hard-coded for now.

  // Initialise local matrices and vectors
  array1d< globalIndex > dispEqnRowIndices( maxNumUdof );
  array1d< globalIndex > jumpEqnRowIndices( 3 );
  array1d< globalIndex > dispColIndices( maxNumUdof );
  array1d< globalIndex > jumpColIndices( 3 );

  array2d< real64 > Kwu_elem( 3, maxNumUdof );
  array2d< real64 > Kuw_elem( maxNumUdof, 3 );
  array2d< real64 > Kww_elem( 3, 3 );
  array1d< real64 > R1( 3 );
  array1d< real64 > R0( maxNumUdof );
  array1d< real64 > tractionVec( 3 );
  array2d< real64 > dTdw( 3, 3 );

  // Equilibrium and compatibility operators for the element
  // number of strain components x number of jump enrichments. The comp operator is different
  // at each Gauss point.
  array2d< real64 > eqMatrix( 3, 6 );
  array2d< real64 > compMatrix( 6, 3 );
  array2d< real64 > strainMatrix( 6, maxNumUdof );

  // local storage of contribution of each gauss point
  array2d< real64 > Kwu_gauss( 3, maxNumUdof );
  array2d< real64 > Kuw_gauss( maxNumUdof, 3 );
  array2d< real64 > Kww_gauss( 3, 3 );

  // intermediate objects to do BDC, EDB, EDC
  array2d< real64 > matBD( maxNumUdof, 6 );
  array2d< real64 > matED( 3, 6 );

  array1d< R1Tensor > u_local( 8 );
  array1d< R1Tensor > du_local( 8 );

  array1d< real64 > u( maxNumUdof );
  array1d< real64 > w( 3 );

  array2d< real64 > dMatrix( 6, 6 );

  // begin region loop
  elemManager.forElementRegions< EmbeddedSurfaceRegion >(
    [&]( EmbeddedSurfaceRegion const & embeddedRegion ) -> void {
      // loop of embeddeSubregions
      embeddedRegion.forElementSubRegions< EmbeddedSurfaceSubRegion >(
        [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion ) -> void {
          localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();
          arrayView1d< localIndex const > const & embeddedSurfaceToRegion =
            embeddedSurfaceSubRegion.getSurfaceToRegionList();
          arrayView1d< localIndex const > const & embeddedSurfaceToSubRegion =
            embeddedSurfaceSubRegion.getSurfaceToSubRegionList();
          arrayView1d< localIndex const > const & embeddedSurfaceToCell =
            embeddedSurfaceSubRegion.getSurfaceToCellList();

          arrayView1d< globalIndex const > const & embeddedElementDofNumber =
            embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >(
              jumpDofKey );
          arrayView1d< R1Tensor const > const & w_global =
            embeddedSurfaceSubRegion.getReference< array1d< R1Tensor > >(
              viewKeyStruct::dispJumpString );
          arrayView1d< R1Tensor const > const & dw_global =
            embeddedSurfaceSubRegion.getReference< array1d< R1Tensor > >(
              viewKeyStruct::deltaDispJumpString );

          arrayView1d< real64 const > const & fractureSurfaceArea =
            embeddedSurfaceSubRegion.getElementArea();

          arrayView1d< integer const > const & ghostRank =
            embeddedSurfaceSubRegion.ghostRank();

          jumpEqnRowIndices.resize(
            embeddedSurfaceSubRegion.numOfJumpEnrichments() );
          jumpColIndices.resize( embeddedSurfaceSubRegion.numOfJumpEnrichments() );

          // loop over embedded surfaces
          for( localIndex k = 0; k < numEmbeddedElems; ++k )
          {
            if( ghostRank[k] < 0 )
            {
              // Get rock matrix element subregion
              CellElementSubRegion const * const elementSubRegion =
                Group::group_cast< CellElementSubRegion const * const >(
                  elemManager.GetRegion( embeddedSurfaceToRegion[k] )
                    ->GetSubRegion( embeddedSurfaceToSubRegion[k] ) );
              CellBlock::NodeMapType const & elemsToNodes =
                elementSubRegion->nodeList();
              // Get the number of nodes per element
              localIndex const numNodesPerElement = elemsToNodes.size( 1 );

              // Get finite element discretization info
              finiteElement::FiniteElementBase const & fe =
                elementSubRegion->getReference< finiteElement::FiniteElementBase >(
                  m_solidSolver->getDiscretizationName() );

              // Resize based on number of dof of the subregion
              int nUdof = numNodesPerElement * 3;
              dispEqnRowIndices.resize( nUdof );
              dispColIndices.resize( nUdof );
              Kwu_elem.resizeDimension< 1 >( nUdof );
              Kuw_elem.resizeDimension< 0 >( nUdof );
              R0.resize( nUdof );
              u.resize( nUdof );

              // Initialize
              Kwu_elem.setValues< serialPolicy >( 0 );
              Kuw_elem.setValues< serialPolicy >( 0 );
              Kww_elem.setValues< serialPolicy >( 0 );
              R0.setValues< serialPolicy >( 0 );
              R1.setValues< serialPolicy >( 0 );
              dTdw.setValues< serialPolicy >( 0 );
              eqMatrix.setValues< serialPolicy >( 0 );
              compMatrix.setValues< serialPolicy >( 0 );
              strainMatrix.setValues< serialPolicy >( 0 );

              u_local.resize( numNodesPerElement );
              du_local.resize( numNodesPerElement );

              // Get mechanical moduli tensor
              LinearElasticIsotropic const * constitutiveRelation =
                elementSubRegion->getConstitutiveModel< LinearElasticIsotropic >(
                  m_solidSolver->solidMaterialNames()[0] );
              LinearElasticIsotropic::KernelWrapper const & solidConstitutive =
                constitutiveRelation->createKernelUpdates();
              solidConstitutive.GetStiffness( embeddedSurfaceToCell[k], dMatrix );

              // Basis functions derivatives
              arrayView4d< real64 const > const & dNdX = elementSubRegion->dNdX();

              // transformation determinant
              arrayView2d< real64 const > const & detJ = elementSubRegion->detJ();

              // Fill in equilibrium operator
              arrayView1d< real64 const > const & cellVolume =
                elementSubRegion->getElementVolume();
              real64 hInv = fractureSurfaceArea[k] /
                cellVolume[embeddedSurfaceToCell[k]];  // AreaFrac / cellVolume
              AssembleEquilibriumOperator( eqMatrix,
                                           embeddedSurfaceSubRegion,
                                           k,
                                           hInv );

              // Dof index of nodal displacements and row indices
              for( localIndex a = 0; a < numNodesPerElement; ++a )
              {
                localIndex localNodeIndex =
                  elemsToNodes[embeddedSurfaceToCell[k]][a];

                for( int i = 0; i < dim; ++i )
                {
                  dispEqnRowIndices[static_cast< int >( a ) * dim + i] =
                    globalDofNumber[localNodeIndex] + i - rankOffset;
                  dispColIndices[static_cast< int >( a ) * dim + i] =
                    globalDofNumber[localNodeIndex] + i;
                }
              }

              for( localIndex i = 0; i < numNodesPerElement; ++i )
              {
                localIndex const nodeID =
                  elemsToNodes( embeddedSurfaceToCell[k], i );
                u_local[i] = disp[nodeID];
                du_local[i] = dDisp[nodeID];
              }

              // Dof number of jump enrichment
              for( int i = 0; i < embeddedSurfaceSubRegion.numOfJumpEnrichments();
                   i++ )
              {
                jumpEqnRowIndices[i] =
                  embeddedElementDofNumber[k] + i - rankOffset;
                jumpColIndices[i] = embeddedElementDofNumber[k] + i;
                w( i ) = w_global[k][i] + 0 * dw_global[k][i];
              }

              // copy values in the R1Tensor object to use in BlasLapack interface
              for( localIndex j = 0; j < numNodesPerElement; ++j )
              {
                for( int i = 0; i < dim; ++i )
                {
                  u( j * dim + i ) = u_local[j][i] + 0 * du_local[j][i];
                }
              }

              // 1. Assembly of element matrices

              // Resize local storages.
              Kwu_gauss.resizeDimension< 1 >( nUdof );
              Kuw_gauss.resizeDimension< 0 >( nUdof );
              matBD.resizeDimension< 0 >( nUdof );
              strainMatrix.resizeDimension< 1 >( nUdof );

              BlasLapackLA::matrixMatrixMultiply( eqMatrix, dMatrix, matED );

              // Compute traction
              ComputeTraction( constitutiveManager, w, tractionVec, dTdw );

              for( integer q = 0; q < fe.getNumQuadraturePoints(); ++q )
              {
                const realT detJq = detJ[embeddedSurfaceToCell[k]][q];
                AssembleCompatibilityOperator( compMatrix,
                                               embeddedSurfaceSubRegion,
                                               k,
                                               q,
                                               elemsToNodes,
                                               nodesCoord,
                                               embeddedSurfaceToCell,
                                               numNodesPerElement,
                                               dNdX );

                AssembleStrainOperator( strainMatrix,
                                        embeddedSurfaceToCell[k],
                                        q,
                                        numNodesPerElement,
                                        dNdX );
                // transp(B)D
                BlasLapackLA::matrixTMatrixMultiply( strainMatrix, dMatrix, matBD );
                // EDC
                BlasLapackLA::matrixMatrixMultiply( matED, compMatrix, Kww_gauss );
                // EDB
                BlasLapackLA::matrixMatrixMultiply( matED, strainMatrix, Kwu_gauss );
                // transp(B)DB
                BlasLapackLA::matrixMatrixMultiply( matBD, compMatrix, Kuw_gauss );

                // multiply by determinant
                BlasLapackLA::matrixScale( detJq, Kwu_gauss );
                BlasLapackLA::matrixScale( detJq, Kuw_gauss );
                BlasLapackLA::matrixScale( detJq, Kww_gauss );

                // Add Gauss point contribution to element matrix
                BlasLapackLA::matrixMatrixAdd( Kww_gauss, Kww_elem, -1 );
                BlasLapackLA::matrixMatrixAdd( Kwu_gauss, Kwu_elem, -1 );
                BlasLapackLA::matrixMatrixAdd( Kuw_gauss, Kuw_elem, -1 );
              }

              BlasLapackLA::matrixMatrixAdd( dTdw, Kww_elem, -1 );

              BlasLapackLA::matrixVectorMultiply( Kww_elem, w, R1 );
              BlasLapackLA::matrixVectorMultiply( Kwu_elem, u, R1, 1, 1 );

              BlasLapackLA::matrixVectorMultiply( Kuw_elem, w, R0, 1, 1 );

              BlasLapackLA::vectorVectorAdd( tractionVec, R1, 1 );

              // 2. Assembly into global system
              // fill in residual vector
              for( localIndex i = 0; i < dispEqnRowIndices.size(); ++i )
              {
                if( dispEqnRowIndices[i] >= 0 &&
                    dispEqnRowIndices[i] < localMatrix.numRows() )
                {
                  RAJA::atomicAdd< parallelDeviceAtomic >(
                    &localRhs[dispEqnRowIndices[i]],
                    R0[i] );

                  localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                    dispEqnRowIndices[i],
                    jumpColIndices.data(),
                    Kuw_elem[i],
                    embeddedSurfaceSubRegion.numOfJumpEnrichments() );
                }
              }

              for( localIndex i = 0; i < jumpEqnRowIndices.size(); ++i )
              {
                if( jumpEqnRowIndices[i] >= 0 &&
                    jumpEqnRowIndices[i] < localMatrix.numRows() )
                {
                  RAJA::atomicAdd< parallelDeviceAtomic >(
                    &localRhs[jumpEqnRowIndices[i]],
                    R1[i] );

                  // fill in matrix
                  localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                    jumpEqnRowIndices[i],
                    jumpColIndices.data(),
                    Kww_elem[i],
                    embeddedSurfaceSubRegion.numOfJumpEnrichments() );

                  localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
                    jumpEqnRowIndices[i],
                    dispColIndices.data(),
                    Kwu_elem[i],
                    numNodesPerElement * dim );
                }
              }
            }
          }   // loop over embedded surfaces
        } );  // subregion loop
    } );      // region loop
}

void
SolidMechanicsEmbeddedFractures::AddCouplingNumNonzeros(
  DomainPartition & domain,
  DofManager & dofManager,
  arrayView1d< localIndex > const & rowLengths ) const
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< globalIndex_array >( dispDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >(
    [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion ) {
      localIndex const numEmbeddedElems = embeddedSurfaceSubRegion.size();
      arrayView1d< localIndex const > const & embeddedSurfaceToRegion =
        embeddedSurfaceSubRegion.getSurfaceToRegionList();
      arrayView1d< localIndex const > const & embeddedSurfaceToSubRegion =
        embeddedSurfaceSubRegion.getSurfaceToSubRegionList();
      arrayView1d< localIndex const > const & embeddedSurfaceToCell =
        embeddedSurfaceSubRegion.getSurfaceToCellList();

      arrayView1d< globalIndex const > const & embeddedElementDofNumber =
        embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank =
        embeddedSurfaceSubRegion.ghostRank();

      for( localIndex k = 0; k < numEmbeddedElems; ++k )
      {
        CellBlock const * const subRegion =
          Group::group_cast< CellBlock const * const >(
            elemManager.GetRegion( embeddedSurfaceToRegion[k] )
              ->GetSubRegion( embeddedSurfaceToSubRegion[k] ) );

        if( ghostRank[k] < 0 )
        {
          localIndex const localRow = LvArray::integerConversion< localIndex >(
            embeddedElementDofNumber[k] - rankOffset );
          GEOSX_ASSERT_GE( localRow, 0 );
          GEOSX_ASSERT_GE(
            rowLengths.size(),
            localRow + embeddedSurfaceSubRegion.numOfJumpEnrichments() );

          for( localIndex i = 0;
               i < embeddedSurfaceSubRegion.numOfJumpEnrichments();
               ++i )
          {
            rowLengths[localRow + i] += 3 * subRegion->numNodesPerElement();
          }

          for( localIndex a = 0; a < subRegion->numNodesPerElement(); ++a )
          {
            const localIndex & node =
              subRegion->nodeList( embeddedSurfaceToCell[k], a );
            localIndex const localDispRow =
              LvArray::integerConversion< localIndex >( dispDofNumber[node] -
                                                        rankOffset );
            GEOSX_ASSERT_GE( localDispRow, 0 );
            GEOSX_ASSERT_GE( rowLengths.size(),
                             localDispRow + 3 * subRegion->numNodesPerElement() );

            for( int d = 0; d < 3; ++d )
            {
              rowLengths[localDispRow + d] +=
                embeddedSurfaceSubRegion.numOfJumpEnrichments();
            }
          }
        }
      }
    } );
}

void
SolidMechanicsEmbeddedFractures::AddCouplingSparsityPattern(
  DomainPartition const & domain,
  DofManager const & dofManager,
  SparsityPatternView< globalIndex > const & pattern ) const
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );
  string const dispDofKey = dofManager.getKey( keys::TotalDisplacement );

  arrayView1d< globalIndex const > const & dispDofNumber =
    nodeManager.getReference< globalIndex_array >( dispDofKey );

  globalIndex const rankOffset = dofManager.rankOffset();

  static constexpr int maxNumDispDof = 3 * 8;  // this is hard-coded for now.

  elemManager.forElementSubRegions< EmbeddedSurfaceSubRegion >(
    [&]( EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion ) {
      arrayView1d< localIndex const > const & embeddedSurfaceToRegion =
        embeddedSurfaceSubRegion.getSurfaceToRegionList();
      arrayView1d< localIndex const > const & embeddedSurfaceToSubRegion =
        embeddedSurfaceSubRegion.getSurfaceToSubRegionList();
      arrayView1d< localIndex const > const & embeddedSurfaceToCell =
        embeddedSurfaceSubRegion.getSurfaceToCellList();

      arrayView1d< globalIndex const > const & embeddedElementDofNumber =
        embeddedSurfaceSubRegion.getReference< array1d< globalIndex > >( jumpDofKey );

      // Insert the entries corresponding to jump-disp coupling
      // This will fill K_wu, and K_uw
      for( localIndex k = 0; k < embeddedSurfaceSubRegion.size(); ++k )
      {
        CellBlock const * const elemSubRegion =
          Group::group_cast< CellBlock const * const >(
            elemManager.GetRegion( embeddedSurfaceToRegion[k] )
              ->GetSubRegion( embeddedSurfaceToSubRegion[k] ) );

        // working arrays
        stackArray1d< globalIndex, maxNumDispDof > eqnRowIndicesDisp(
          3 * elemSubRegion->numNodesPerElement() );
        stackArray1d< globalIndex, 3 > eqnRowIndicesJump(
          embeddedSurfaceSubRegion.numOfJumpEnrichments() );
        stackArray1d< globalIndex, maxNumDispDof > dofColIndicesDisp(
          3 * elemSubRegion->numNodesPerElement() );
        stackArray1d< globalIndex, 3 > dofColIndicesJump(
          embeddedSurfaceSubRegion.numOfJumpEnrichments() );

        for( localIndex idof = 0;
             idof < embeddedSurfaceSubRegion.numOfJumpEnrichments();
             ++idof )
        {
          eqnRowIndicesJump[idof] =
            embeddedElementDofNumber[k] + idof - rankOffset;
          dofColIndicesJump[idof] = embeddedElementDofNumber[k] + idof;
        }

        for( localIndex a = 0; a < elemSubRegion->numNodesPerElement(); ++a )
        {
          const localIndex & node =
            elemSubRegion->nodeList( embeddedSurfaceToCell[k], a );
          for( localIndex idof = 0; idof < 3; ++idof )
          {
            eqnRowIndicesDisp[3 * a + idof] =
              dispDofNumber[node] + idof - rankOffset;
            dofColIndicesDisp[3 * a + idof] = dispDofNumber[node] + idof;
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesDisp.size(); ++i )
        {
          if( eqnRowIndicesDisp[i] >= 0 && eqnRowIndicesDisp[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesJump.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesDisp[i], dofColIndicesJump[j] );
            }
          }
        }

        for( localIndex i = 0; i < eqnRowIndicesJump.size(); ++i )
        {
          if( eqnRowIndicesJump[i] >= 0 && eqnRowIndicesJump[i] < pattern.numRows() )
          {
            for( localIndex j = 0; j < dofColIndicesDisp.size(); ++j )
            {
              pattern.insertNonZero( eqnRowIndicesJump[i], dofColIndicesDisp[j] );
            }
          }
        }
      }
    } );
}

void
SolidMechanicsEmbeddedFractures::AssembleEquilibriumOperator(
  array2d< real64 > & eqMatrix,
  EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
  const localIndex k,
  const real64 hInv )
{
  GEOSX_MARK_FUNCTION;
  // Normal and tangent unit vectors
  R1Tensor const nVec = embeddedSurfaceSubRegion.getNormalVector( k );
  R1Tensor const tVec1 = embeddedSurfaceSubRegion.getTangentVector1( k );
  R1Tensor const tVec2 = embeddedSurfaceSubRegion.getTangentVector2( k );

  BlasLapackLA::matrixScale( 0, eqMatrix );

  real64 nDn[3][3], t1DnSym[3][3], t2DnSym[3][3];

  // n dyadic n
  LvArray::tensorOps::AiBj< 3, 3 >( nDn, nVec, nVec );

  // sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::AiBj< 3, 3 >( t1DnSym, nVec, tVec1 );
  LvArray::tensorOps::plusAiBj< 3, 3 >( t1DnSym, tVec1, nVec );
  LvArray::tensorOps::scale< 3, 3 >( t1DnSym, 0.5 );

  LvArray::tensorOps::AiBj< 3, 3 >( t2DnSym, nVec, tVec2 );
  LvArray::tensorOps::plusAiBj< 3, 3 >( t2DnSym, tVec2, nVec );
  LvArray::tensorOps::scale< 3, 3 >( t2DnSym, 0.5 );

  int VoigtIndex;

  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      if( i == j )
      {
        VoigtIndex = 1;
      }
      else
      {
        VoigtIndex = 6 - i - j;
      }
      eqMatrix( 0, VoigtIndex ) += nDn[i][j];
      eqMatrix( 1, VoigtIndex ) += t1DnSym[i][j];
      eqMatrix( 2, VoigtIndex ) += t2DnSym[i][j];
    }
  }
  BlasLapackLA::matrixScale( -hInv, eqMatrix );
}

void
SolidMechanicsEmbeddedFractures::AssembleCompatibilityOperator(
  array2d< real64 > & compMatrix,
  EmbeddedSurfaceSubRegion const & embeddedSurfaceSubRegion,
  localIndex const k,
  localIndex const q,
  CellBlock::NodeMapType const & elemsToNodes,
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodesCoord,
  arrayView1d< localIndex const > const & embeddedSurfaceToCell,
  localIndex const numNodesPerElement,
  arrayView4d< real64 const > const & dNdX )
{
  GEOSX_MARK_FUNCTION;

  // Normal and tangent unit vectors
  R1Tensor const nVec = embeddedSurfaceSubRegion.getNormalVector( k );
  R1Tensor const tVec1 = embeddedSurfaceSubRegion.getTangentVector1( k );
  R1Tensor const tVec2 = embeddedSurfaceSubRegion.getTangentVector2( k );

  // Fill in compatibility operator

  // 1. construct mvector sum(dNdX(a) * H(a)) value for each Gauss point
  R1Tensor mVec;
  real64 heavisideFun;
  mVec = 0.0;
  for( integer a = 0; a < numNodesPerElement; ++a )
  {
    // Heaviside
    heavisideFun = embeddedSurfaceSubRegion.ComputeHeavisideFunction(
      nodesCoord[elemsToNodes[embeddedSurfaceToCell[k]][a]],
      k );
    // sum contribution of each node
    mVec[0] -= dNdX( embeddedSurfaceToCell[k], q, a, 0 ) * heavisideFun;
    mVec[1] -= dNdX( embeddedSurfaceToCell[k], q, a, 1 ) * heavisideFun;
    mVec[2] -= dNdX( embeddedSurfaceToCell[k], q, a, 2 ) * heavisideFun;
  }

  BlasLapackLA::matrixScale( 0, compMatrix );

  // 2. fill in the operator itself

  real64 nDmSym[3][3], t1DmSym[3][3], t2DmSym[3][3];

  // sym(n dyadic m)
  LvArray::tensorOps::AiBj< 3, 3 >( nDmSym, mVec, nVec );
  LvArray::tensorOps::plusAiBj< 3, 3 >( nDmSym, nVec, mVec );
  LvArray::tensorOps::scale< 3, 3 >( nDmSym, 0.5 );

  // sym(n dyadic t1) and sym (n dyadic t2)
  LvArray::tensorOps::AiBj< 3, 3 >( t1DmSym, mVec, tVec1 );
  LvArray::tensorOps::plusAiBj< 3, 3 >( t1DmSym, tVec1, mVec );
  LvArray::tensorOps::scale< 3, 3 >( t1DmSym, 0.5 );

  LvArray::tensorOps::AiBj< 3, 3 >( t2DmSym, mVec, tVec2 );
  LvArray::tensorOps::plusAiBj< 3, 3 >( t2DmSym, tVec2, mVec );
  LvArray::tensorOps::scale< 3, 3 >( t2DmSym, 0.5 );

  int VoigtIndex;

  for( int i = 0; i < 3; ++i )
  {
    for( int j = 0; j < 3; ++j )
    {
      if( i == j )
      {
        VoigtIndex = 1;
      }
      else
      {
        VoigtIndex = 6 - i - j;
      }
      compMatrix( VoigtIndex, 0 ) += nDmSym[i][j];
      compMatrix( VoigtIndex, 1 ) += t1DmSym[i][j];
      compMatrix( VoigtIndex, 2 ) += t2DmSym[i][j];
    }
  }
}

void
SolidMechanicsEmbeddedFractures::AssembleStrainOperator(
  array2d< real64 > & strainMatrix,
  localIndex const elIndex,
  localIndex const q,
  localIndex const numNodesPerElement,
  arrayView4d< real64 const > const & dNdX )
{
  GEOSX_MARK_FUNCTION;
  strainMatrix.setValues< serialPolicy >( 0 );  // make 0

  for( integer a = 0; a < numNodesPerElement; ++a )
  {
    strainMatrix( 0, a * 3 + 0 ) = dNdX( elIndex, q, a, 0 );
    strainMatrix( 1, a * 3 + 1 ) = dNdX( elIndex, q, a, 1 );
    strainMatrix( 2, a * 3 + 2 ) = dNdX( elIndex, q, a, 2 );

    strainMatrix( 3, a * 3 + 1 ) = dNdX( elIndex, q, a, 2 );
    strainMatrix( 3, a * 3 + 2 ) = dNdX( elIndex, q, a, 1 );

    strainMatrix( 4, a * 3 + 0 ) = dNdX( elIndex, q, a, 2 );
    strainMatrix( 4, a * 3 + 2 ) = dNdX( elIndex, q, a, 0 );

    strainMatrix( 5, a * 3 + 0 ) = dNdX( elIndex, q, a, 1 );
    strainMatrix( 5, a * 3 + 1 ) = dNdX( elIndex, q, a, 0 );
  }
}

void
SolidMechanicsEmbeddedFractures::ApplyBoundaryConditions(
  real64 const time,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplyBoundaryConditions( time,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );
}

real64
SolidMechanicsEmbeddedFractures::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // Matrix residual
  real64 const solidResidualNorm =
    m_solidSolver->CalculateResidualNorm( domain, dofManager, localRhs );
  // Fracture residual
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  string const jumpDofKey = dofManager.getKey( viewKeyStruct::dispJumpString );

  globalIndex const rankOffset = dofManager.rankOffset();

  RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

  // globalResidualNorm[0]: the sum of all the local sum(rhs^2).
  // globalResidualNorm[1]: max of max force of each rank. Basically max force globally
  real64 globalResidualNorm[2] = { 0, 0 };

  forTargetSubRegions< EmbeddedSurfaceSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, EmbeddedSurfaceSubRegion const & subRegion ) {
      GEOSX_UNUSED_VAR( targetIndex );

      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( jumpDofKey );
      arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

      for( localIndex k = 0; k < subRegion.size(); ++k )
      {
        if( ghostRank[k] < 0 )
        {
          localIndex const localRow =
            LvArray::integerConversion< localIndex >( dofNumber[k] - rankOffset );
          for( localIndex i = 0; i < subRegion.numOfJumpEnrichments(); ++i )
          {
            localSum += localRhs[localRow + i] * localRhs[localRow + i];
          }
        }
      }

      real64 const localResidualNorm[2] = { localSum.get(),
                                            m_solidSolver->getMaxForce() };

      int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
      int const numRanks = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
      array1d< real64 > globalValues( numRanks * 2 );

      // Everything is done on rank 0
      MpiWrapper::gather( localResidualNorm,
                          2,
                          globalValues.data(),
                          2,
                          0,
                          MPI_COMM_GEOSX );

      if( rank == 0 )
      {
        for( int r = 0; r < numRanks; ++r )
        {
          // sum/max across all ranks
          globalResidualNorm[0] += globalValues[r * 2];
          globalResidualNorm[1] =
            std::max( globalResidualNorm[1], globalValues[r * 2 + 1] );
        }
      }

      MpiWrapper::bcast( globalResidualNorm, 2, 0, MPI_COMM_GEOSX );
    } );

  real64 const fractureResidualNorm = sqrt( globalResidualNorm[0] ) /
    ( globalResidualNorm[1] + 1 );  // the + 1 is for the first
  // time-step when maxForce = 0;

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    char output[200] = { 0 };
    sprintf( output, "( RFracture ) = (%4.2e) ; ", fractureResidualNorm );
    std::cout << output;
  }

  return sqrt( solidResidualNorm * solidResidualNorm +
               fractureResidualNorm * fractureResidualNorm );
}

void
SolidMechanicsEmbeddedFractures::ApplySystemSolution(
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor,
  DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  m_solidSolver->ApplySystemSolution( dofManager,
                                      localSolution,
                                      scalingFactor,
                                      domain );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dispJumpString,
                               viewKeyStruct::deltaDispJumpString,
                               -scalingFactor );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dispJumpString,
                               viewKeyStruct::dispJumpString,
                               -scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::dispJumpString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaDispJumpString ) );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         domain.getMeshBody( 0 )->getMeshLevel( 0 ),
                                         domain.getNeighbors(),
                                         true );
}

void
SolidMechanicsEmbeddedFractures::ComputeTraction(
  ConstitutiveManager const * const constitutiveManager,
  array1d< real64 > const & dispJump,
  array1d< real64 > & tractionVector,
  array2d< real64 > & dTdw )
{
  // Compute traction vector on the fracture element
  ContactRelationBase const * const contactRelation =
    constitutiveManager->GetGroup< ContactRelationBase >( m_contactRelationName );

  // check if fracture is open
  bool open = dispJump[0] >= 0 ? true : false;

  if( open )
  {
    tractionVector[0] = 1e5;
    tractionVector[1] = 0.0;
    tractionVector[2] = 0.0;
    dTdw( 0, 0 ) = 0.0;
    dTdw( 0, 1 ) = 0.0;
    dTdw( 0, 2 ) = 0.0;
    dTdw( 1, 0 ) = 0.0;
    dTdw( 1, 1 ) = 0.0;
    dTdw( 1, 2 ) = 0.0;
    dTdw( 2, 0 ) = 0.0;
    dTdw( 2, 1 ) = 0.0;
    dTdw( 2, 2 ) = 0.0;
  }
  else
  {
    // Contact through penalty condition.
    tractionVector[0] = contactRelation->stiffness() * dispJump[0];
    tractionVector[1] = 0;
    tractionVector[2] = 0;
  }
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        SolidMechanicsEmbeddedFractures,
                        std::string const &,
                        Group * const )
} /* namespace geosx */
