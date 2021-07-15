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
 * @file SinglePhaseHybridFVM.cpp
 */

#include "SinglePhaseHybridFVM.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "finiteVolume/HybridMimeticDiscretization.hpp"
#include "finiteVolume/MimeticInnerProductDispatch.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mainInterface/ProblemManager.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseHybridFVMKernels;
using namespace mimeticInnerProduct;

SinglePhaseHybridFVM::SinglePhaseHybridFVM( const string & name,
                                            Group * const parent ):
  SinglePhaseBase( name, parent ),
  m_faceDofKey( "" ),
  m_areaRelTol( 1e-8 )
{

  // one cell-centered dof per cell
  m_numDofPerCell = 1;
  m_linearSolverParameters.get().mgr.strategy = LinearSolverParameters::MGR::StrategyType::singlePhaseHybridFVM;

}


void SinglePhaseHybridFVM::registerDataOnMesh( Group & meshBodies )
{

  // 1) Register the cell-centered data
  SinglePhaseBase::registerDataOnMesh( meshBodies );

  // 2) Register the face data
  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
    FaceManager & faceManager = meshLevel.getFaceManager();

    // primary variables: face pressures changes
    faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
  } );
}

void SinglePhaseHybridFVM::initializePreSubGroups()
{
  SinglePhaseBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( !fvManager.hasGroup< HybridMimeticDiscretization >( m_discretizationName ) )
  {
    GEOSX_ERROR( "The HybridMimeticDiscretization must be selected with SinglePhaseHybridFVM" );
  }
}

void SinglePhaseHybridFVM::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();
  FaceManager const & faceManager = mesh.getFaceManager();
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );

  // in the flux kernel, we need to make sure that we act only on the target regions
  // for that, we need the following region filter
  for( string const & regionName : targetRegionNames() )
  {
    m_regionFilter.insert( elemManager.getRegions().getIndex( regionName ) );
  }

  // check that multipliers are stricly larger than 0, which would work with SinglePhaseFVM, but not with SinglePhaseHybridFVM.
  // To deal with a 0 multiplier, we would just have to skip the corresponding face in the FluxKernel
  string const & coeffName = hmDiscretization.getReference< string >( HybridMimeticDiscretization::viewKeyStruct::coeffNameString() );
  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( coeffName + HybridMimeticDiscretization::viewKeyStruct::transMultiplierString() );

  RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );
  forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
  {
    minVal.min( transMultiplier[iface] );
  } );

  GEOSX_ERROR_IF_LE_MSG( minVal.get(), 0.0,
                         "The transmissibility multipliers used in SinglePhaseHybridFVM must strictly larger than 0.0" );
}

void SinglePhaseHybridFVM::implicitStepSetup( real64 const & time_n,
                                              real64 const & dt,
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // setup the cell-centered fields
  SinglePhaseBase::implicitStepSetup( time_n, dt, domain );

  // setup the face fields
  MeshLevel & meshLevel     = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = meshLevel.getFaceManager();

  // get the accumulated pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() );

  // zero out the face pressures
  dFacePres.zero();
}

void SinglePhaseHybridFVM::implicitStepComplete( real64 const & time_n,
                                                 real64 const & dt,
                                                 DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // increment the cell-centered fields
  SinglePhaseBase::implicitStepComplete( time_n, dt, domain );

  // increment the face fields
  MeshLevel & meshLevel     = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = meshLevel.getFaceManager();

  // get the face-based pressures
  arrayView1d< real64 > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString() );
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() );

  forAll< parallelDevicePolicy<> >( faceManager.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iface )
  {
    facePres[iface] += dFacePres[iface];
  } );
}

void SinglePhaseHybridFVM::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                      DofManager & dofManager ) const
{

  // setup the connectivity of elem fields
  // we need Connectivity::Face because of the two-point upwinding
  // in AssembleOneSidedMassFluxes
  dofManager.addField( viewKeyStruct::pressureString(),
                       DofManager::Location::Elem,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::pressureString(),
                          viewKeyStruct::pressureString(),
                          DofManager::Connector::Face );

  // setup the connectivity of face fields
  dofManager.addField( viewKeyStruct::facePressureString(),
                       DofManager::Location::Face,
                       1,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::facePressureString(),
                          viewKeyStruct::facePressureString(),
                          DofManager::Connector::Elem );

  // setup coupling between pressure and face pressure
  dofManager.addCoupling( viewKeyStruct::facePressureString(),
                          viewKeyStruct::pressureString(),
                          DofManager::Connector::Elem );
}

void SinglePhaseHybridFVM::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const dt,
                                              DomainPartition const & domain,
                                              DofManager const & dofManager,
                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh          = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  HybridMimeticDiscretization const & hmDiscretization = fvManager.getHybridMimeticDiscretization( m_discretizationName );
  MimeticInnerProductBase const & mimeticInnerProductBase =
    hmDiscretization.getReference< MimeticInnerProductBase >( HybridMimeticDiscretization::viewKeyStruct::innerProductString() );

  // node data (for transmissibility computation)

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition = nodeManager.referencePosition();

  // face data

  // get the face-based DOF numbers for the assembly
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString() );
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );
  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();

  // get the element dof numbers for the assembly
  string const & elemDofKey = dofManager.getKey( viewKeyStruct::pressureString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
    mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( elemDofKey );
  elemDofNumber.setName( getName() + "/accessors/" + elemDofKey );

  // get the face-centered pressures
  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString() );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() );

  // get the face-centered depth
  arrayView1d< real64 const > const & faceGravCoef =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

  // get the face-centered transMultiplier
  string const & coeffName = hmDiscretization.getReference< string >( HybridMimeticDiscretization::viewKeyStruct::coeffNameString() );
  arrayView1d< real64 const > const & transMultiplier =
    faceManager.getReference< array1d< real64 > >( coeffName + HybridMimeticDiscretization::viewKeyStruct::transMultiplierString() );

  // get the face-to-nodes connectivity for the transmissibility calculation
  ArrayOfArraysView< localIndex const > const & faceToNodes = faceManager.nodeList().toViewConst();

  arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

  // tolerance for transmissibility calculation
  real64 const lengthTolerance = domain.getMeshBody( 0 ).getGlobalLengthScale() * m_areaRelTol;

  forTargetSubRegionsComplete< CellElementSubRegion >( mesh,
                                                       [&]( localIndex const targetIndex,
                                                            localIndex const er,
                                                            localIndex const esr,
                                                            ElementRegionBase const &,
                                                            auto const & subRegion )
  {
    SingleFluidBase const & fluid =
      getConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    arrayView2d< real64 const > const & permeability =
      getConstitutiveModel< PermeabilityBase >( subRegion, m_permeabilityModelNames[targetIndex] ).permeability();

    mimeticInnerProductDispatch( mimeticInnerProductBase,
                                 [&] ( auto const mimeticInnerProduct )
    {
      using IP_TYPE = TYPEOFREF( mimeticInnerProduct );

      KernelLaunchSelector< IP_TYPE, FluxKernel >( subRegion.numFacesPerElement(),
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
                                                   permeability,
                                                   faceGravCoef,
                                                   transMultiplier,
                                                   m_mobility.toNestedViewConst(),
                                                   m_dMobility_dPres.toNestedViewConst(),
                                                   elemDofNumber.toNestedViewConst(),
                                                   dofManager.rankOffset(),
                                                   lengthTolerance,
                                                   dt,
                                                   localMatrix,
                                                   localRhs );
    } );
  } );
}

void SinglePhaseHybridFVM::applyBoundaryConditions( real64 const time_n,
                                                    real64 const dt,
                                                    DomainPartition & domain,
                                                    DofManager const & dofManager,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  SinglePhaseBase::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
}


real64 SinglePhaseHybridFVM::calculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager const & faceManager = mesh.getFaceManager();

  // here we compute the cell-centered residual norm in the derived class
  // to avoid duplicating a synchronization point

  // get a view into local residual vector

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString() );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  // local residual
  real64 localResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };
  real64 globalResidualNorm[4] = { 0.0, 0.0, 0.0, 0.0 };

  // 1. Compute the residual for the mass conservation equations

  // compute the norm of local residual scaled by cell pore volume

  real64 defaultViscosity = 0; // for the normalization of the face residuals
  localIndex subRegionCounter = 0;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase const & subRegion )
  {

    arrayView1d< globalIndex const > const & elemDofNumber = subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & refPoro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & densOld = subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString() );

    SinglePhaseBaseKernels::ResidualNormKernel::launch< parallelDevicePolicy<>,
                                                        parallelDeviceReduce >( localRhs,
                                                                                rankOffset,
                                                                                elemDofNumber,
                                                                                elemGhostRank,
                                                                                refPoro,
                                                                                volume,
                                                                                densOld,
                                                                                localResidualNorm );

    SingleFluidBase const & fluid = getConstitutiveModel< SingleFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    defaultViscosity += fluid.defaultViscosity();
    subRegionCounter++;
  } );

  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView2d< localIndex const > const & elemRegionList    = faceManager.elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList = faceManager.elementSubRegionList();
  arrayView2d< localIndex const > const & elemList          = faceManager.elementList();

  defaultViscosity /= subRegionCounter;

  // 2. Compute the residual for the face-based constraints
  SinglePhaseHybridFVMKernels::ResidualNormKernel::launch< parallelDevicePolicy<>,
                                                           parallelDeviceReduce >( localRhs,
                                                                                   rankOffset,
                                                                                   faceDofNumber.toNestedViewConst(),
                                                                                   faceGhostRank.toNestedViewConst(),
                                                                                   elemRegionList.toNestedViewConst(),
                                                                                   elemSubRegionList.toNestedViewConst(),
                                                                                   elemList.toNestedViewConst(),
                                                                                   m_volume.toNestedViewConst(),
                                                                                   defaultViscosity,
                                                                                   &localResidualNorm[3] );

  // 3. Combine the two norms

  // compute global residual norm
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         4,
                         MPI_SUM,
                         MPI_COMM_GEOSX );


  real64 const elemResidualNorm = sqrt( globalResidualNorm[0] )
                                  / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );
  real64 const faceResidualNorm = sqrt( globalResidualNorm[3] );

  real64 const residualNorm = ( elemResidualNorm > faceResidualNorm )
                            ? elemResidualNorm
                            : faceResidualNorm;

  return residualNorm;
}


bool SinglePhaseHybridFVM::checkSystemSolution( DomainPartition const & domain,
                                                DofManager const & dofManager,
                                                arrayView1d< real64 const > const & localSolution,
                                                real64 const scalingFactor )
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager const & faceManager = mesh.getFaceManager();

  localIndex localCheck = 1;

  string const elemDofKey = dofManager.getKey( viewKeyStruct::pressureString() );
  string const faceDofKey = dofManager.getKey( viewKeyStruct::facePressureString() );

  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & elemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( elemDofKey );
    arrayView1d< integer const > const & elemGhostRank =
      subRegion.ghostRank();

    arrayView1d< real64 const > const & pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

    localIndex const subRegionSolutionCheck =
      SinglePhaseBaseKernels::SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                                           parallelDeviceReduce >( localSolution,
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

  arrayView1d< integer const > const & faceGhostRank = faceManager.ghostRank();
  arrayView1d< globalIndex const > const & faceDofNumber =
    faceManager.getReference< array1d< globalIndex > >( faceDofKey );

  arrayView1d< real64 const > const & facePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::facePressureString() );
  arrayView1d< real64 const > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() );

  localIndex const faceSolutionCheck =
    SinglePhaseBaseKernels::SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                                         parallelDeviceReduce >( localSolution,
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

  return MpiWrapper::min( localCheck );
}


void SinglePhaseHybridFVM::applySystemSolution( DofManager const & dofManager,
                                                arrayView1d< real64 const > const & localSolution,
                                                real64 const scalingFactor,
                                                DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // here we apply the cell-centered update in the derived class
  // to avoid duplicating a synchronization point

  // 1. apply the cell-centered update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::pressureString(),
                               viewKeyStruct::deltaPressureString(),
                               scalingFactor );

  // 2. apply the face-based update

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::facePressureString(),
                               viewKeyStruct::deltaFacePressureString(),
                               scalingFactor );

  // 3. synchronize

  // the tags in fieldNames have to match the tags used in NeighborCommunicator.cpp
  std::map< string, string_array > fieldNames;
  fieldNames["face"].emplace_back( string( viewKeyStruct::deltaFacePressureString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    updateState( subRegion, targetIndex );
  } );
}


void SinglePhaseHybridFVM::resetStateToBeginningOfStep( DomainPartition & domain )
{
  // 1. Reset the cell-centered fields
  SinglePhaseBase::resetStateToBeginningOfStep( domain );

  // 2. Reset the face-based fields
  MeshLevel & mesh          = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = mesh.getFaceManager();

  // get the accumulated face pressure updates
  arrayView1d< real64 > const & dFacePres =
    faceManager.getReference< array1d< real64 > >( viewKeyStruct::deltaFacePressureString() );

  // zero out the face pressures
  dFacePres.zero();
}

REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseHybridFVM, string const &, Group * const )
} /* namespace geosx */
