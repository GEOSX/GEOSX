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
 * @file PoroelasticSolver.cpp
 *
 */

#include "PoroelasticSolver.hpp"

#include "../solidMechanics/SolidMechanicsPoroElasticKernel.hpp"
#include "common/DataLayouts.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/solid/PoroElastic.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "finiteElement/Kinematics.h"
#include "linearAlgebra/solvers/BlockPreconditioner.hpp"
#include "linearAlgebra/solvers/SeparateComponentPreconditioner.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBase.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsLagrangianFEM.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

PoroelasticSolver::PoroelasticSolver( const std::string & name, Group * const parent ) :
  SolverBase( name, parent ),
  m_solidSolverName(),
  m_flowSolverName(),
  m_couplingTypeOptionString(),
  m_couplingTypeOption()

{
  registerWrapper( viewKeyStruct::solidSolverNameString, &m_solidSolverName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of the solid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::fluidSolverNameString, &m_flowSolverName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Name of the fluid mechanics solver to use in the poroelastic solver" );

  registerWrapper( viewKeyStruct::couplingTypeOptionStringString,
                   &m_couplingTypeOptionString )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Coupling option: (FIM, SIM_FixedStress)" );

  m_linearSolverParameters.get().mgr.strategy = "Poroelastic";
  m_linearSolverParameters.get().mgr.separateComponents = true;
  m_linearSolverParameters.get().mgr.displacementFieldName =
    keys::TotalDisplacement;
  m_linearSolverParameters.get().dofsPerNode = 3;
}

void
PoroelasticSolver::RegisterDataOnMesh( dataRepository::Group * const MeshBodies )
{
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    ElementRegionManager * const elemManager =
      mesh.second->group_cast< MeshBody * >()->getMeshLevel( 0 )->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >(
      [&]( ElementSubRegionBase & elementSubRegion ) {
        elementSubRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::totalMeanStressString )
          ->setDescription( "Total Mean Stress" );
        elementSubRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::oldTotalMeanStressString )
          ->setDescription( "Total Mean Stress" );
      } );
  }
}

void
PoroelasticSolver::SetupDofs( DomainPartition const & domain,
                              DofManager & dofManager ) const
{
  GEOSX_MARK_FUNCTION;
  m_solidSolver->SetupDofs( domain, dofManager );
  m_flowSolver->SetupDofs( domain, dofManager );

  dofManager.addCoupling( keys::TotalDisplacement,
                          FlowSolverBase::viewKeyStruct::pressureString,
                          DofManager::Connector::Elem );
}

void
PoroelasticSolver::SetupSystem( DomainPartition & domain,
                                DofManager & dofManager,
                                CRSMatrix< real64, globalIndex > & localMatrix,
                                array1d< real64 > & localRhs,
                                array1d< real64 > & localSolution,
                                bool const setSparsity )
{
  if( m_precond )
  {
    m_precond->clear();
  }

  // setup monolithic coupled system
  SolverBase::SetupSystem( domain,
                           dofManager,
                           localMatrix,
                           localRhs,
                           localSolution,
                           setSparsity );

  if( !m_precond && m_linearSolverParameters.get().solverType != "direct" )
  {
    CreatePreconditioner();
  }
}

void
PoroelasticSolver::ImplicitStepSetup( real64 const & time_n,
                                      real64 const & dt,
                                      DomainPartition & domain )
{
  m_flowSolver->ImplicitStepSetup( time_n, dt, domain );
  m_solidSolver->ImplicitStepSetup( time_n, dt, domain );

  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

    forTargetSubRegions( mesh,
                         [&]( localIndex const, ElementSubRegionBase & subRegion ) {
                           arrayView1d< real64 const > const & totalMeanStress =
                             subRegion.getReference< array1d< real64 > >(
                               viewKeyStruct::totalMeanStressString );
                           arrayView1d< real64 > const & oldTotalMeanStress =
                             subRegion.getReference< array1d< real64 > >(
                               viewKeyStruct::oldTotalMeanStressString );

                           forAll< parallelDevicePolicy<> >(
                             subRegion.size(),
                             [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                               oldTotalMeanStress[ei] = totalMeanStress[ei];
                             } );
                         } );
  }
}

void
PoroelasticSolver::ImplicitStepComplete( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition & domain )
{
  m_solidSolver->ImplicitStepComplete( time_n, dt, domain );
  m_flowSolver->ImplicitStepComplete( time_n, dt, domain );
}

void
PoroelasticSolver::PostProcessInput()
{
  SolverBase::PostProcessInput();

  m_flowSolver = this->getParent()->GetGroup< SinglePhaseBase >( m_flowSolverName );
  m_solidSolver =
    this->getParent()->GetGroup< SolidMechanicsLagrangianFEM >( m_solidSolverName );

  GEOSX_ERROR_IF( m_flowSolver == nullptr,
                  "Flow solver not found or invalid type: " << m_flowSolverName );
  GEOSX_ERROR_IF( m_solidSolver == nullptr,
                  "Solid solver not found or invalid type: " << m_solidSolverName );

  string ctOption =
    this->getReference< string >( viewKeyStruct::couplingTypeOptionStringString );

  m_solidSolver->setEffectiveStress( 1 );

  if( ctOption == "SIM_FixedStress" )
  {
    m_couplingTypeOption = couplingTypeOption::SIM_FixedStress;

    // For this coupled solver the minimum number of Newton Iter should be 0 for both flow and solid solver,
    // otherwise it will never converge.
    m_flowSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
    m_solidSolver->getNonlinearSolverParameters().m_minIterNewton = 0;
  }
  else if( ctOption == "FIM" )
  {
    m_couplingTypeOption = couplingTypeOption::FIM;
  }
  else
  {
    GEOSX_ERROR( "invalid coupling type option: " + ctOption );
  }
}

void
PoroelasticSolver::InitializePostInitialConditions_PreSubGroups(
  Group * const problemManager )
{
  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    m_flowSolver->setPoroElasticCoupling();
    // Calculate initial total mean stress
    UpdateDeformationForCoupling(
      *problemManager->GetGroup< DomainPartition >( keys::domain ) );
  }
}

PoroelasticSolver::~PoroelasticSolver()
{
  // TODO Auto-generated destructor stub
}

void
PoroelasticSolver::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  m_flowSolver->ResetStateToBeginningOfStep( domain );
  m_solidSolver->ResetStateToBeginningOfStep( domain );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    arrayView1d< real64 const > const & oldTotalMeanStress =
      subRegion.getReference< array1d< real64 > >(
        viewKeyStruct::oldTotalMeanStressString );
    arrayView1d< real64 > const & totalMeanStress =
      subRegion.getReference< array1d< real64 > >(
        viewKeyStruct::totalMeanStressString );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                        totalMeanStress[ei] = oldTotalMeanStress[ei];
                                      } );
  } );
}

real64
PoroelasticSolver::SolverStep( real64 const & time_n,
                               real64 const & dt,
                               int const cycleNumber,
                               DomainPartition & domain )
{
  real64 dt_return = dt;
  if( m_couplingTypeOption == couplingTypeOption::SIM_FixedStress )
  {
    dt_return = SplitOperatorStep( time_n, dt, cycleNumber, domain );
  }
  else if( m_couplingTypeOption == couplingTypeOption::FIM )
  {
    SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

    ImplicitStepSetup( time_n, dt, domain );

    dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

    ImplicitStepComplete( time_n, dt_return, domain );
  }
  return dt_return;
}

void
PoroelasticSolver::UpdateDeformationForCoupling( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager & nodeManager = *mesh.getNodeManager();

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager.referencePosition();
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & u =
    nodeManager.totalDisplacement();

  forTargetSubRegionsComplete< CellElementSubRegion >(
    mesh,
    [&]( localIndex const,
         localIndex const,
         localIndex const,
         ElementRegionBase & elemRegion,
         CellElementSubRegion & elementSubRegion ) {
      string const & solidName =
        m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex(
          elemRegion.getName() )];
      SolidBase const & solid =
        GetConstitutiveModel< SolidBase >( elementSubRegion, solidName );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes =
        elementSubRegion.nodeList();

      arrayView1d< real64 > const & totalMeanStress =
        elementSubRegion.getReference< array1d< real64 > >(
          viewKeyStruct::totalMeanStressString );

      arrayView1d< real64 > const & oldTotalMeanStress =
        elementSubRegion.getReference< array1d< real64 > >(
          viewKeyStruct::oldTotalMeanStressString );

      arrayView1d< real64 const > const & pres =
        elementSubRegion.getReference< array1d< real64 > >(
          FlowSolverBase::viewKeyStruct::pressureString );

      arrayView1d< real64 const > const & dPres =
        elementSubRegion.getReference< array1d< real64 > >(
          FlowSolverBase::viewKeyStruct::deltaPressureString );

      arrayView1d< real64 > const & poro =
        elementSubRegion.getReference< array1d< real64 > >(
          SinglePhaseBase::viewKeyStruct::porosityString );

      arrayView1d< real64 const > const & poroOld =
        elementSubRegion.getReference< array1d< real64 > >(
          SinglePhaseBase::viewKeyStruct::porosityOldString );

      arrayView1d< real64 const > const & volume =
        elementSubRegion.getReference< array1d< real64 > >(
          CellBlock::viewKeyStruct::elementVolumeString );

      arrayView1d< real64 > const & dVol =
        elementSubRegion.getReference< array1d< real64 > >(
          SinglePhaseBase::viewKeyStruct::deltaVolumeString );

      arrayView1d< real64 const > const & bulkModulus =
        solid.getReference< array1d< real64 > >( "BulkModulus" );

      real64 const biotCoefficient =
        solid.getReference< real64 >( "BiotCoefficient" );

      arrayView3d< real64 const, solid::STRESS_USD > const & stress =
        solid.getStress();

      localIndex const numNodesPerElement = elemsToNodes.size( 1 );
      finiteElement::FiniteElementBase const & fe =
        elementSubRegion.getReference< finiteElement::FiniteElementBase >(
          m_solidSolver->getDiscretizationName() );
      localIndex const numQuadraturePoints = fe.getNumQuadraturePoints();

      // TODO: remove use of R1Tensor and use device policy
      forAll< parallelDevicePolicy< 32 > >(
        elementSubRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          real64 effectiveMeanStress = 0.0;
          for( localIndex q = 0; q < numQuadraturePoints; ++q )
          {
            effectiveMeanStress +=
              ( stress( ei, q, 0 ) + stress( ei, q, 1 ) + stress( ei, q, 2 ) );
          }
          effectiveMeanStress /= ( 3 * numQuadraturePoints );

          totalMeanStress[ei] =
            effectiveMeanStress - biotCoefficient * ( pres[ei] + dPres[ei] );

          poro[ei] = poroOld[ei] +
            ( biotCoefficient - poroOld[ei] ) / bulkModulus[ei] *
              ( totalMeanStress[ei] - oldTotalMeanStress[ei] + dPres[ei] );

          // update element volume
          R1Tensor Xlocal[ElementRegionManager::maxNumNodesPerElem];
          for( localIndex a = 0; a < numNodesPerElement; ++a )
          {
            Xlocal[a] = X[elemsToNodes[ei][a]];
            Xlocal[a] += u[elemsToNodes[ei][a]];
          }

          dVol[ei] = computationalGeometry::HexVolume( Xlocal ) - volume[ei];
        } );
    } );
}

void
PoroelasticSolver::AssembleSystem(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  // assemble J_SS
  //  m_solidSolver->AssembleSystem( time_n, dt,
  //                                 domain,
  //                                 dofManager,
  //                                 localMatrix,
  //                                 localRhs );

  m_solidSolver
    ->AssemblyLaunch< constitutive::PoroElasticBase,
                      SolidMechanicsLagrangianFEMKernels::QuasiStaticPoroElastic >(
      domain,
      dofManager,
      localMatrix,
      localRhs );

  // assemble J_FF
  m_flowSolver->AssembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

  // assemble J_SF
  AssembleCouplingTerms( domain, dofManager, localMatrix, localRhs );
}

void
PoroelasticSolver::AssembleCouplingTerms(
  DomainPartition const & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();

  string const uDofKey = dofManager.getKey( keys::TotalDisplacement );
  arrayView1d< globalIndex const > const & uDofNumber =
    nodeManager.getReference< globalIndex_array >( uDofKey );

  arrayView2d< real64 const, nodes::INCR_DISPLACEMENT_USD > const & incr_disp =
    nodeManager.incrementalDisplacement();

  globalIndex const rankOffset = dofManager.rankOffset();
  string const pDofKey =
    dofManager.getKey( FlowSolverBase::viewKeyStruct::pressureString );

  // begin subregion loop
  forTargetSubRegionsComplete< CellElementSubRegion >(
    mesh,
    [&]( localIndex const,
         localIndex const,
         localIndex const,
         ElementRegionBase const & region,
         CellElementSubRegion const & elementSubRegion ) {
      string const & fluidName =
        m_flowSolver
          ->fluidModelNames()[m_flowSolver->targetRegionIndex( region.getName() )];
      SingleFluidBase const & fluid =
        GetConstitutiveModel< SingleFluidBase >( elementSubRegion, fluidName );

      string const & solidName =
        m_solidSolver->solidMaterialNames()[m_solidSolver->targetRegionIndex(
          region.getName() )];
      SolidBase const & solid =
        GetConstitutiveModel< SolidBase >( elementSubRegion, solidName );

      arrayView4d< real64 const > const & dNdX = elementSubRegion.dNdX();

      arrayView2d< real64 const > const & detJ = elementSubRegion.detJ();

      arrayView1d< globalIndex const > const & pDofNumber =
        elementSubRegion.getReference< globalIndex_array >( pDofKey );

      arrayView2d< localIndex const, cells::NODE_MAP_USD > const & elemsToNodes =
        elementSubRegion.nodeList();
      localIndex const numNodesPerElement = elemsToNodes.size( 1 );

      finiteElement::FiniteElementBase const & fe =
        elementSubRegion.getReference< finiteElement::FiniteElementBase >(
          m_solidSolver->getDiscretizationName() );
      localIndex const numQuadraturePoints = fe.getNumQuadraturePoints();

      real64 const biotCoefficient =
        solid.getReference< real64 >( "BiotCoefficient" );

      arrayView2d< real64 const > const & density = fluid.density();

      int dim = 3;
      localIndex constexpr maxNumUDof =
        24;  // TODO: assuming linear HEX at most for the moment
      localIndex constexpr maxNumPDof =
        1;  // TODO: assuming piecewise constant (P0) only for the moment
      localIndex const nUDof = dim * numNodesPerElement;
      localIndex const nPDof = m_flowSolver->numDofPerCell();
      GEOSX_ERROR_IF_GT( nPDof, maxNumPDof );

      // TODO: remove use of R1Tensor and use device policy
      forAll< parallelDevicePolicy< 32 > >(
        elementSubRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const k ) {
          stackArray2d< real64, maxNumUDof * maxNumPDof > dRsdP( nUDof, nPDof );
          stackArray2d< real64, maxNumUDof * maxNumPDof > dRfdU( nPDof, nUDof );
          stackArray1d< real64, maxNumPDof > Rf( nPDof );

          for( integer q = 0; q < numQuadraturePoints; ++q )
          {
            const real64 detJq = detJ[k][q];

            for( integer a = 0; a < numNodesPerElement; ++a )
            {
              dRsdP( a * dim + 0, 0 ) += biotCoefficient * dNdX[k][q][a][0] * detJq;
              dRsdP( a * dim + 1, 0 ) += biotCoefficient * dNdX[k][q][a][1] * detJq;
              dRsdP( a * dim + 2, 0 ) += biotCoefficient * dNdX[k][q][a][2] * detJq;
              dRfdU( 0, a * dim + 0 ) +=
                density[k][0] * biotCoefficient * dNdX[k][q][a][0] * detJq;
              dRfdU( 0, a * dim + 1 ) +=
                density[k][0] * biotCoefficient * dNdX[k][q][a][1] * detJq;
              dRfdU( 0, a * dim + 2 ) +=
                density[k][0] * biotCoefficient * dNdX[k][q][a][2] * detJq;

              localIndex localNodeIndex = elemsToNodes[k][a];

              real64 Rf_tmp = dNdX[k][q][a][0] * incr_disp[localNodeIndex][0] +
                dNdX[k][q][a][1] * incr_disp[localNodeIndex][1] +
                dNdX[k][q][a][2] * incr_disp[localNodeIndex][2];
              Rf_tmp *= density[k][0] * biotCoefficient * detJq;
              Rf[0] += Rf_tmp;
            }
          }

          stackArray1d< globalIndex, maxNumUDof > elementULocalDofIndex( nUDof );
          stackArray1d< globalIndex, maxNumPDof > elementPLocalDofIndex( nPDof );

          // Get dof local to global mapping
          for( localIndex a = 0; a < numNodesPerElement; ++a )
          {
            for( int i = 0; i < dim; ++i )
            {
              elementULocalDofIndex[a * dim + i] =
                uDofNumber[elemsToNodes[k][a]] + i;
            }
          }
          for( localIndex i = 0; i < nPDof; ++i )
          {
            elementPLocalDofIndex[i] = pDofNumber[k] + i;
          }

          for( localIndex i = 0; i < nUDof; ++i )
          {
            localIndex const dof = LvArray::integerConversion< localIndex >(
              elementULocalDofIndex[i] - rankOffset );
            if( dof < 0 || dof >= localMatrix.numRows() ) continue;
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
              dof,
              elementPLocalDofIndex.data(),
              dRsdP[i].dataIfContiguous(),
              nPDof );
          }
          for( localIndex i = 0; i < nPDof; ++i )
          {
            localIndex const dof = LvArray::integerConversion< localIndex >(
              elementPLocalDofIndex[i] - rankOffset );
            if( dof < 0 || dof >= localMatrix.numRows() ) continue;
            localMatrix.addToRowBinarySearchUnsorted< parallelDeviceAtomic >(
              dof,
              elementULocalDofIndex.data(),
              dRfdU[i].dataIfContiguous(),
              nUDof );

            RAJA::atomicAdd< parallelDeviceAtomic >( &localRhs[dof], Rf[i] );
          }
        } );
    } );
}

void
PoroelasticSolver::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  m_solidSolver->ApplyBoundaryConditions( time_n,
                                          dt,
                                          domain,
                                          dofManager,
                                          localMatrix,
                                          localRhs );

  m_flowSolver->ApplyBoundaryConditions( time_n,
                                         dt,
                                         domain,
                                         dofManager,
                                         localMatrix,
                                         localRhs );
}

real64
PoroelasticSolver::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localRhs )
{
  // compute norm of momentum balance residual equations
  real64 const momementumResidualNorm =
    m_solidSolver->CalculateResidualNorm( domain, dofManager, localRhs );

  // compute norm of mass balance residual equations
  real64 const massResidualNorm =
    m_flowSolver->CalculateResidualNorm( domain, dofManager, localRhs );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    char output[200] = { 0 };
    sprintf( output,
             "    ( Rsolid, Rfluid ) = ( %4.2e, %4.2e )",
             momementumResidualNorm,
             massResidualNorm );
    std::cout << output << std::endl;
  }

  return sqrt( momementumResidualNorm * momementumResidualNorm +
               massResidualNorm * massResidualNorm );
}

void
PoroelasticSolver::CreatePreconditioner()
{
  if( m_linearSolverParameters.get().preconditionerType == "block" )
  {
    auto precond = std::make_unique< BlockPreconditioner< LAInterface > >(
      BlockShapeOption::UpperTriangular,
      SchurComplementOption::RowsumDiagonalProbing,
      BlockScalingOption::FrobeniusNorm );

    auto mechPrecond = LAInterface::createPreconditioner(
      m_solidSolver->getLinearSolverParameters() );
    precond->setupBlock(
      0,
      { { keys::TotalDisplacement, 0, 3 } },
      std::make_unique< SeparateComponentPreconditioner< LAInterface > >(
        3,
        std::move( mechPrecond ) ) );

    auto flowPrecond = LAInterface::createPreconditioner(
      m_flowSolver->getLinearSolverParameters() );
    precond->setupBlock( 1,
                         { { SinglePhaseBase::viewKeyStruct::pressureString, 0, 1 } },
                         std::move( flowPrecond ) );

    m_precond = std::move( precond );
  }
  else
  {
    //TODO: Revisit this part such that is coherent across physics solver
    //m_precond = LAInterface::createPreconditioner( m_linearSolverParameters.get() );
  }
}

void
PoroelasticSolver::SolveSystem( DofManager const & dofManager,
                                ParallelMatrix & matrix,
                                ParallelVector & rhs,
                                ParallelVector & solution )
{
  solution.zero();
  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void
PoroelasticSolver::ApplySystemSolution(
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor,
  DomainPartition & domain )
{
  // update displacement field
  m_solidSolver->ApplySystemSolution( dofManager,
                                      localSolution,
                                      scalingFactor,
                                      domain );
  // update pressure field
  m_flowSolver->ApplySystemSolution( dofManager,
                                     localSolution,
                                     -scalingFactor,
                                     domain );
}

real64
PoroelasticSolver::SplitOperatorStep( real64 const & time_n,
                                      real64 const & dt,
                                      integer const cycleNumber,
                                      DomainPartition & domain )
{
  real64 dtReturn = dt;
  real64 dtReturnTemporary;

  m_flowSolver->SetupSystem( domain,
                             m_flowSolver->getDofManager(),
                             m_flowSolver->getLocalMatrix(),
                             m_flowSolver->getLocalRhs(),
                             m_flowSolver->getLocalSolution() );

  m_solidSolver->SetupSystem( domain,
                              m_solidSolver->getDofManager(),
                              m_solidSolver->getLocalMatrix(),
                              m_solidSolver->getLocalRhs(),
                              m_solidSolver->getLocalSolution() );

  ImplicitStepSetup( time_n, dt, domain );

  int iter = 0;
  while( iter < m_nonlinearSolverParameters.m_maxIterNewton )
  {
    if( iter == 0 )
    {
      // reset the states of all child solvers if any of them has been reset
      ResetStateToBeginningOfStep( domain );
    }

    GEOSX_LOG_LEVEL_RANK_0( 1, "\tIteration: " << iter + 1 << ", FlowSolver: " );

    dtReturnTemporary =
      m_flowSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }

    if( m_flowSolver->getNonlinearSolverParameters().m_numNewtonIterations == 0 &&
        iter > 0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1,
                              "***** The iterative coupling has converged in "
                                << iter << " iterations! *****\n" );
      break;
    }

    GEOSX_LOG_LEVEL_RANK_0( 1,
                            "\tIteration: " << iter + 1 << ", MechanicsSolver: " );

    m_solidSolver->ResetStressToBeginningOfStep( domain );
    dtReturnTemporary =
      m_solidSolver->NonlinearImplicitStep( time_n, dtReturn, cycleNumber, domain );

    if( dtReturnTemporary < dtReturn )
    {
      iter = 0;
      dtReturn = dtReturnTemporary;
      continue;
    }
    if( m_solidSolver->getNonlinearSolverParameters().m_numNewtonIterations > 0 )
    {
      UpdateDeformationForCoupling( domain );
    }
    ++iter;
  }

  ImplicitStepComplete( time_n, dt, domain );

  return dtReturn;
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        PoroelasticSolver,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
