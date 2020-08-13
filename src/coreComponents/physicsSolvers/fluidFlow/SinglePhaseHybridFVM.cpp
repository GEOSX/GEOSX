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
 * @file SinglePhaseHybridFVM.cpp
 */

#include "SinglePhaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseHybridFVMKernels;
using namespace HybridFVMInnerProduct;

SinglePhaseHybridFVM::SinglePhaseHybridFVM( const std::string & name,
                                            Group * const parent ) :
  SinglePhaseBase( name, parent ),
  m_faceDofKey( "" ),
  m_areaRelTol( 1e-8 )
{
  // one cell-centered dof per cell
  m_numDofPerCell = 1;
}

void
SinglePhaseHybridFVM::RegisterDataOnMesh( Group * const MeshBodies )
{
  // 1) Register the cell-centered data
  SinglePhaseBase::RegisterDataOnMesh( MeshBodies );

  // 2) Register the face data
  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * const meshLevel =
      Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    FaceManager * const faceManager = meshLevel->getFaceManager();

    // primary variables: face pressures changes
    faceManager
      ->registerWrapper< array1d< real64 > >( viewKeyStruct::deltaFacePressureString )
      ->setPlotLevel( PlotLevel::LEVEL_0 )
      ->setRegisteringObjects( this->getName() )
      ->setDescription(
        "An array that holds the accumulated pressure updates at the faces." );
  }
}

void
SinglePhaseHybridFVM::InitializePostInitialConditions_PreSubGroups(
  Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel const & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  // in the flux kernel, we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  for( string const & regionName : targetRegionNames() )
  {
    m_regionFilter.insert( elemManager.GetRegions().getIndex( regionName ) );
  }
}

void
SinglePhaseHybridFVM::ImplicitStepSetup( real64 const & time_n,
                                         real64 const & dt,
                                         DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::ImplicitStepSetup( time_n, dt, domain );

  // setup the face fields
  MeshLevel & meshLevel =
    *domain.getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the accumulated pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >(
      viewKeyStruct::deltaFacePressureString );

  // zero out the face pressures
  dFacePres.setValues< parallelDevicePolicy<> >( 0.0 );
}

void
SinglePhaseHybridFVM::ImplicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  SinglePhaseBase::ImplicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel & meshLevel = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *meshLevel.getFaceManager();

  // get the face-based pressures
  arrayView1d< real64 > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >(
      viewKeyStruct::deltaFacePressureString );

  forAll< parallelDevicePolicy<> >( faceManager.size(),
                                    [=] GEOSX_HOST_DEVICE( localIndex const iface ) {
                                      facePres[iface] += dFacePres[iface];
                                      dFacePres[iface] = 0.0;
                                    } );
}

void
SinglePhaseHybridFVM::SetupDofs(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager & dofManager ) const
{
  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::pressureString,
                       DofManager::Location::Elem,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::pressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString,
                       DofManager::Location::Face,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::facePressureString,
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::facePressureString,
                          viewKeyStruct::pressureString,
                          DofManager::Connector::Elem,
                          true );
}

void
SinglePhaseHybridFVM::AssembleFluxTerms(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const dt,
  DomainPartition const & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  NodeManager const & nodeManager = *mesh.getNodeManager();
  FaceManager const & faceManager = *mesh.getFaceManager();

  // node data (for transmissibility computation)

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition =
    nodeManager.referencePosition();

  // face data

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

  // get the element dof numbers for the assembly
  string const & elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
    mesh.getElemManager()->ConstructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
  elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

  // get the face-centered pressures
  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >(
      viewKeyStruct::deltaFacePressureString );

  // get the face-centered depth
  arrayView1d< real64 const > const & faceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView< localIndex const > const & faceToNodes =
    faceManager.nodeList().toViewConst();

  arrayView2d< localIndex const > const & elemRegionList =
    faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList =
    faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = faceManager.elementList();

  // tolerance for transmissibility calculation
  real64 const lengthTolerance =
    domain.getMeshBody( 0 )->getGlobalLengthScale() * m_areaRelTol;

  forTargetSubRegionsComplete< CellElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex,
         localIndex const er,
         localIndex const esr,
         ElementRegionBase const &,
         auto const & subRegion ) {
      SingleFluidBase const & fluid =
        GetConstitutiveModel< SingleFluidBase >( subRegion,
                                                 m_fluidModelNames[targetIndex] );

      KernelLaunchSelector< FluxKernel >( subRegion.numFacesPerElement(),
                                          er,
                                          esr,
                                          subRegion,
                                          fluid,
                                          m_regionFilter.toViewConst(),
                                          nodePosition,
                                          elemRegionList,
                                          elemSubRegionList,
                                          elemList,
                                          faceToNodes,
                                          faceDofNumber,
                                          faceGhostRank,
                                          facePres,
                                          dFacePres,
                                          faceGravCoef,
                                          m_mobility.toViewConst(),
                                          m_dMobility_dPres.toViewConst(),
                                          elemDofNumber.toViewConst(),
                                          dofManager.rankOffset(),
                                          lengthTolerance,
                                          dt,
                                          localMatrix,
                                          localRhs );
    } );
}

void
SinglePhaseHybridFVM::ApplyBoundaryConditions(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_MARK_FUNCTION;
  // will implement boundary conditions later
}

real64
SinglePhaseHybridFVM::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();

  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point

  // get a view into local residual vector

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );

  globalIndex const rankOffset = dofManager.rankOffset();

  // local residual
  real64 localResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };
  real64 globalResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };

  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume

  real64 defaultViscosity = 0;  // for the normalization of the face residuals
  localIndex subRegionCounter = 0;

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & refPoro =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::referencePorosityString );
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & densOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );

      SinglePhaseBaseKernels::ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                                          parallelDeviceReduce >(
        localRhs,
        rankOffset,
        elemDofNumber,
        elemGhostRank,
        refPoro,
        volume,
        densOld,
        localResidualNorm );

      SingleFluidBase const & fluid =
        GetConstitutiveModel< SingleFluidBase >( subRegion,
                                                 m_fluidModelNames[targetIndex] );
      defaultViscosity += fluid.defaultViscosity();
      subRegionCounter++;
    } );

  arrayView1d< integer const > const & faceGhostRank =
    faceManager.getReference< array1d< integer > >(
      ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView2d< localIndex const > const & elemRegionList =
    faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList =
    faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = faceManager.elementList();

  defaultViscosity /= subRegionCounter;

  // 2. Compute the residual for the face-based constraints
  SinglePhaseHybridFVMKernels::ResidualNormKernel::Launch< parallelDevicePolicy<>,
                                                           parallelDeviceReduce >(
    localRhs,
    rankOffset,
    faceDofNumber.toViewConst(),
    faceGhostRank.toViewConst(),
    elemRegionList.toViewConst(),
    elemSubRegionList.toViewConst(),
    elemList.toViewConst(),
    m_volume.toViewConst(),
    defaultViscosity,
    &localResidualNorm[3] );

  // 3. Combine the two norms

  // compute global residual norm
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         4,
                         MPI_SUM,
                         MPI_COMM_GEOSX );

  real64 const elemResidualNorm = sqrt( globalResidualNorm[0] ) /
    ( ( globalResidualNorm[1] + m_fluxEstimate ) / ( globalResidualNorm[2] + 1 ) );
  real64 const faceResidualNorm = sqrt( globalResidualNorm[3] );

  real64 const residualNorm =
    ( elemResidualNorm > faceResidualNorm ) ? elemResidualNorm : faceResidualNorm;

  return residualNorm;
}

bool
SinglePhaseHybridFVM::CheckSystemSolution(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor )
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager const & faceManager = *mesh.getFaceManager();

  localIndex localCheck = 1;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString );

  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetSubRegions(
    mesh,
    [&]( localIndex const, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & elemDofNumber =
        subRegion.getReference< array1d< globalIndex > >( elemDofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & pres =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
      arrayView1d< real64 const > const & dPres =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::deltaPressureString );

      localIndex const subRegionSolutionCheck =
        SinglePhaseBaseKernels::SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                                             parallelDeviceReduce >(
          localSolution,
          rankOffset,
          elemDofNumber,
          elemGhostRank,
          pres,
          dPres,
          scalingFactor );

      if( subRegionSolutionCheck == 0 )
      {
        localCheck = 0;
      }
    } );

  arrayView1d< integer const > const & faceGhostRank =
    faceManager.getReference< array1d< integer > >(
      ObjectManagerBase::viewKeyStruct::ghostRankString );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >(
      viewKeyStruct::deltaFacePressureString );

  localIndex const faceSolutionCheck =
    SinglePhaseBaseKernels::SolutionCheckKernel::Launch< parallelDevicePolicy<>,
                                                         parallelDeviceReduce >(
      localSolution,
      rankOffset,
      faceDofNumber,
      faceGhostRank,
      facePres,
      dFacePres,
      scalingFactor );

  if( faceSolutionCheck == 0 )
  {
    localCheck = 0;
  }

  return MpiWrapper::Min( localCheck );
}

void
SinglePhaseHybridFVM::ApplySystemSolution(
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor,
  DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // here we apply the cell-centered update in the derived class
  // to avoid duplicating a synchronization point

  // 1. apply the cell-centered update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::pressureString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor );

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::facePressureString,
                               viewKeyStruct::deltaFacePressureString,
                               scalingFactor );

  // 3. synchronize

  // the tags in fieldNames have to match the tags used in NeighborCommunicator.cpp
  std::map< string, string_array > fieldNames;
  fieldNames["face"].emplace_back( string( viewKeyStruct::deltaFacePressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );

  CommunicationTools::SynchronizeFields( fieldNames,
                                         &mesh,
                                         domain.getNeighbors(),
                                         true );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      UpdateState( subRegion, targetIndex );
    } );
}

void
SinglePhaseHybridFVM::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseBase::ResetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  FaceManager & faceManager = *mesh.getFaceManager();

  // get the accumulated face pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >(
      viewKeyStruct::deltaFacePressureString );

  // zero out the face pressures
  dFacePres.setValues< parallelDevicePolicy<> >( 0.0 );
}

REGISTER_CATALOG_ENTRY( SolverBase,
                        SinglePhaseHybridFVM,
                        std::string const &,
                        Group * const )
} /* namespace geosx */
