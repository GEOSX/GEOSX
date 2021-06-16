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
 * @file ProppantTransport.cpp
 */

#include "ProppantTransport.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/slurryFluidSelector.hpp"
#include "constitutive/fluid/particleFluidSelector.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "mesh/SurfaceElementRegion.hpp"

#include "physicsSolvers/fluidFlow/ProppantTransportKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ProppantTransportKernels;

ProppantTransport::ProppantTransport( const string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{
  this->registerWrapper( viewKeyStruct::proppantNamesString(), &m_proppantModelNames ).setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of proppant constitutive object to use for this solver." );

  registerWrapper( viewKeyStruct::bridgingFactorString(), &m_bridgingFactor ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bridging factor used for bridging/screen-out calculation" );

  registerWrapper( viewKeyStruct::maxProppantConcentrationString(), &m_maxProppantConcentration ).setApplyDefaultValue( 0.6 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Maximum proppant concentration" );

  registerWrapper( viewKeyStruct::proppantDiameterString(), &m_proppantDiameter ).setApplyDefaultValue( 0.4e-3 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant diameter" );

  registerWrapper( viewKeyStruct::proppantDensityString(), &m_proppantDensity ).setApplyDefaultValue( 2500.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Proppant density" );

  registerWrapper( viewKeyStruct::criticalShieldsNumberString(), &m_criticalShieldsNumber ).setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Critical Shields number" );

  registerWrapper( viewKeyStruct::frictionCoefficientString(), &m_frictionCoefficient ).setApplyDefaultValue( 0.03 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Friction coefficient" );

  registerWrapper( viewKeyStruct::updateProppantPackingString(), &m_updateProppantPacking ).setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag that enables/disables proppant-packing update" );

}

void ProppantTransport::postProcessInput()
{
  FlowSolverBase::postProcessInput();
  checkModelNames( m_proppantModelNames, viewKeyStruct::proppantNamesString() );
}

void ProppantTransport::registerDataOnMesh( Group & meshBodies )
{
  FlowSolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() ).
        setDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() ).
        setDefaultValue( 0.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::componentConcentrationString() ).
        setDefaultValue( 0.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() ).
        setDefaultValue( 0.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::cellBasedFluxString() ).
        setDefaultValue( 0.0 ).
        reference().resizeDimension< 1 >( 3 );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantBoundaryString() ).
        setDefaultValue( 0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString() ).
        setDefaultValue( 0.0 );
    } );


    forTargetSubRegions< FaceElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 FaceElementSubRegion & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::componentConcentrationString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::oldComponentDensityString() );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::cellBasedFluxString() ).
        reference().resizeDimension< 1 >( 3 );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantBoundaryString() );
      subRegion.registerWrapper< array1d< integer > >( viewKeyStruct::isProppantMobileString() );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString() );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString() );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::poroMultiplierString() ).
        setDefaultValue( 1.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::transTMultiplierString() ).
        setDefaultValue( 1.0 ).reference().resizeDimension< 1 >( 3 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString() ).
        setDefaultValue( 0.0 );
    } );
  } );
}

void ProppantTransport::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager & cm = domain.getConstitutiveManager();

  // Validate proppant models in regions
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    validateModelMapping< SlurryFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
    validateModelMapping< ParticleFluidBase >( meshLevel.getElemManager(), m_proppantModelNames );
  }

  SlurryFluidBase const & fluid0 = cm.getConstitutiveRelation< SlurryFluidBase >( m_fluidModelNames[0] );

  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() ).resizeDimension< 1 >( NC );
    } );
  }


}

void ProppantTransport::resizeFractureFields( MeshLevel & mesh )
{
  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    forTargetSubRegions< FaceElementSubRegion >( mesh, [&]( localIndex const,
                                                            FaceElementSubRegion & subRegion )
    {
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString() ).resizeDimension< 1 >( NC );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString() ).resizeDimension< 1 >( NC );
    } );
  }
}

void ProppantTransport::updateFluidModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView2d< real64 const > const componentConc  = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() );
  arrayView2d< real64 const > const dComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    FluidUpdateKernel::launch( fluidWrapper,
                               pres,
                               dPres,
                               componentConc,
                               dComponentConc );
  } );
}

void ProppantTransport::updateComponentDensity( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView2d< real64 const > const componentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() );
  arrayView2d< real64 const > const dComponentConc = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );

  SlurryFluidBase & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();
    ComponentDensityUpdateKernel::launch( fluidWrapper,
                                          pres,
                                          dPres,
                                          componentConc,
                                          dComponentConc );
  } );
}


void ProppantTransport::updateProppantModel( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const proppantConc  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );
  arrayView1d< real64 const > const dProppantConc = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );

  SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView2d< real64 const > const fluidDens            = fluid.fluidDensity();
  arrayView2d< real64 const > const dFluidDens_dPres     = fluid.dFluidDensity_dPressure();
  arrayView3d< real64 const > const dFluidDens_dCompConc = fluid.dFluidDensity_dComponentConcentration();
  arrayView2d< real64 const > const fluidVisc            = fluid.fluidViscosity();
  arrayView2d< real64 const > const dFluidVisc_dPres     = fluid.dFluidViscosity_dPressure();
  arrayView3d< real64 const > const dFluidVisc_dCompConc = fluid.dFluidViscosity_dComponentConcentration();

  ParticleFluidBase & proppant = getConstitutiveModel< ParticleFluidBase >( dataGroup, m_proppantModelNames[targetIndex] );

  constitutiveUpdatePassThru( proppant, [&]( auto & castedProppant )
  {
    typename TYPEOFREF( castedProppant ) ::KernelWrapper proppantWrapper = castedProppant.createKernelWrapper();
    ProppantUpdateKernel::launch( proppantWrapper,
                                  proppantConc,
                                  dProppantConc,
                                  fluidDens,
                                  dFluidDens_dPres,
                                  dFluidDens_dCompConc,
                                  fluidVisc,
                                  dFluidVisc_dPres,
                                  dFluidVisc_dCompConc );
  } );
}

void ProppantTransport::updateProppantMobility( Group & dataGroup )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const conc = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );
  arrayView1d< real64 const > const aperture = dataGroup.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );
  arrayView1d< integer > const isProppantMobile = dataGroup.getReference< array1d< integer > >( viewKeyStruct::isProppantMobileString() );

  real64 const minAperture = m_minAperture;
  real64 const maxProppantConcentration = m_maxProppantConcentration;

  forAll< parallelDevicePolicy<> >( dataGroup.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
  {
    isProppantMobile[a] = aperture[a] > minAperture && conc[a] < maxProppantConcentration;
  } );

}

void ProppantTransport::updateState( Group & dataGroup, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  updateFluidModel( dataGroup, targetIndex );
  updateProppantModel( dataGroup, targetIndex );
}

void ProppantTransport::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::componentConcentrationString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

  resetViews( mesh );

  // We have to redo the below loop after fractures are generated

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    updateState( subRegion, targetIndex );

    SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );
    arrayView3d< real64 const > const componentDens = fluid.componentDensity();
    arrayView2d< real64 > const componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        componentDensOld[ei][c] = componentDens[ei][0][c];
      }
    } );
  } );

  m_minAperture = m_bridgingFactor * m_proppantDiameter;

  real64 const oneMinusMaxConc = (1.0 - m_maxProppantConcentration);
  m_proppantPackPermeability = m_proppantDiameter * m_proppantDiameter / 180.0 *
                               (oneMinusMaxConc * oneMinusMaxConc * oneMinusMaxConc) / (m_maxProppantConcentration * m_maxProppantConcentration);

}

// By now the proppant trasnport solver is coupled with single-phase flow solver. This SolverStep function is used for
// stand-alone proppant transport modeling which is not fully implemented and tested

real64 ProppantTransport::solverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  implicitStepSetup( time_n, dt, domain );

  if( cycleNumber == 0 )
  {
    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::getInstance();

    boundaryConditionManager.applyInitialConditions( domain );

  }

  preStepUpdate( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dtReturn= this->nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.

  implicitStepComplete( time_n, dtReturn, domain );

  postStepUpdate( time_n, dtReturn, domain );

  return dtReturn;
}


void ProppantTransport::preStepUpdate( real64 const & time,
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  FlowSolverBase::precomputeData( mesh );

  NodeManager const & nodeManager = mesh.getNodeManager();
  FaceManager const & faceManager = mesh.getFaceManager();

  if( time <= 0 )
  {
    forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {

      subRegion.calculateElementGeometricQuantities( nodeManager, faceManager );

      updateProppantMobility( subRegion );

      arrayView1d< real64 > const packVf = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() );
      arrayView1d< real64 > const poroMultiplier = subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultiplierString() );
      arrayView2d< real64 > const transTMultiplier = subRegion.getReference< array2d< real64 > >( viewKeyStruct::transTMultiplierString() );
      arrayView1d< integer > const isProppantMobile = subRegion.getReference< array1d< integer > >( viewKeyStruct::isProppantMobileString() );
      arrayView1d< real64 > const proppantLiftFlux = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString() );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        packVf[ei] = 0.0;
        poroMultiplier[ei] = 1.0;
        LvArray::tensorOps::fill< 3 >( transTMultiplier[ei], 1.0 );
        isProppantMobile[ei] = 1;
        proppantLiftFlux[ei] = 0.0;
      } );
    } );
  }

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );
    arrayView3d< real64 const > const componentDens = fluid.componentDensity();
    arrayView2d< real64 > const componentDensOld = subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString() );

    arrayView1d< real64 > const excessPackV = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString() );
    arrayView2d< real64 > const cellBasedFlux = subRegion.getReference< array2d< real64 > >( viewKeyStruct::cellBasedFluxString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        componentDensOld[ei][c] = componentDens[ei][0][c];
      }

      excessPackV[ei] = 0.0;
      LvArray::tensorOps::fill< 3 >( cellBasedFlux[ei], 0.0 );
    } );
  } );

  updateCellBasedFlux( time, domain );
}

void ProppantTransport::postStepUpdate( real64 const & time_n,
                                        real64 const & dt_return,
                                        DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    updateProppantMobility( subRegion );
  } );

  real64 const maxProppantConcentration = m_maxProppantConcentration;
  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView1d< real64 > const & packVf = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() );
    arrayView1d< real64 > const & proppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( proppantConc[ei] >= maxProppantConcentration || packVf[ei] >= 1.0 )
      {
        packVf[ei] = 1.0;
        proppantConc[ei] = maxProppantConcentration;
      }

    } );
  } );

  if( m_updateProppantPacking == 1 )
  {
    updateProppantPackVolume( time_n, dt_return, domain );
  }
}

void ProppantTransport::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  resetViews( mesh );

}

void ProppantTransport::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const proppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );
    arrayView1d< real64 const > const dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );

    arrayView2d< real64 > const componentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() );
    arrayView2d< real64 const > const dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );

    arrayView1d< real64 > const proppantLiftFlux =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      proppantConc[ei] += dProppantConc[ei];
      proppantLiftFlux[ei] = 0.0;

      for( localIndex c = 0; c < NC; ++c )
      {
        componentConc[ei][c] += dComponentConc[ei][c];
      }
    } );
  } );

}

void ProppantTransport::setupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                   DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::proppantConcentrationString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::proppantConcentrationString(),
                          viewKeyStruct::proppantConcentrationString(),
                          DofManager::Connector::Face );
}


void ProppantTransport::assembleSystem( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms( dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  assembleFluxTerms( time,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void ProppantTransport::assembleAccumulationTerms( real64 const dt,
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString() );
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

    arrayView2d< real64 const > const componentDensOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::oldComponentDensityString() );
    arrayView1d< real64 const > const proppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );
    arrayView1d< real64 const > const dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );
    arrayView1d< real64 const > const proppantPackVf =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() );
    arrayView1d< real64 const > const proppantLiftFlux =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantLiftFluxString() );

    SlurryFluidBase const & fluid = getConstitutiveModel< SlurryFluidBase >( subRegion, targetIndex );

    arrayView3d< real64 const > const componentDens = fluid.componentDensity();
    arrayView3d< real64 const > const dCompDens_dPres = fluid.dComponentDensity_dPressure();
    arrayView4d< real64 const > const dCompDens_dCompConc = fluid.dComponentDensity_dComponentConcentration();

    AccumulationKernel::launch( subRegion.size(),
                                m_numComponents,
                                m_numDofPerCell,
                                dofManager.rankOffset(),
                                dofNumber,
                                elemGhostRank,
                                proppantConc,
                                dProppantConc,
                                componentDensOld,
                                componentDens,
                                dCompDens_dPres,
                                dCompDens_dCompConc,
                                volume,
                                proppantPackVf,
                                proppantLiftFlux,
                                dt,
                                m_maxProppantConcentration,
                                localMatrix,
                                localRhs );
  } );
}


void ProppantTransport::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager.constructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const dofNumber = dofNumberAccessor.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const pres  = m_pressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dPres = m_deltaPressure.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const proppantConc   = m_proppantConcentration.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dProppantConc  = m_deltaProppantConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const gravCoef = m_gravCoef.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dens        = m_density.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dDens_dPres = m_dDensity_dPressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dDens_dProppantConc = m_dDensity_dProppantConcentration.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const dDens_dComponentConc = m_dDensity_dComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const visc        = m_viscosity.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dVisc_dPres = m_dViscosity_dPressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dVisc_dProppantConc = m_dViscosity_dProppantConcentration.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const dVisc_dComponentConc = m_dViscosity_dComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const componentDens = m_componentDensity.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const dComponentDens_dPres = m_dComponentDensity_dPressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView4d< real64 const > > const dComponentDens_dComponentConc = m_dComponentDensity_dComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const fluidDensity = m_fluidDensity.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dFluidDens_dPres = m_dFluidDensity_dPressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView3d< real64 const > > const dFluidDens_dComponentConc = m_dFluidDensity_dComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const settlingFactor = m_settlingFactor.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dSettlingFactor_dPres = m_dSettlingFactor_dPressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dSettlingFactor_dProppantConc = m_dSettlingFactor_dProppantConcentration.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dSettlingFactor_dComponentConc = m_dSettlingFactor_dComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const collisionFactor = m_collisionFactor.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dCollisionFactor_dProppantConc = m_dCollisionFactor_dProppantConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const isProppantMobile  = m_isProppantMobile.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const transTMultiplier  = m_transTMultiplier.toNestedViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const
  aperture = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const proppantPackVf  = m_proppantPackVolumeFraction.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const elemGhostRank = m_elemGhostRank.toNestedViewConst();

  fluxApprox.forStencils< SurfaceElementStencil >( mesh, [&]( auto const & stencil )
  {

    FluxKernel::launch( stencil,
                        m_numDofPerCell,
                        dt,
                        dofManager.rankOffset(),
                        transTMultiplier,
                        m_updateProppantPacking,
                        downVector,
                        dofNumber,
                        elemGhostRank,
                        pres,
                        dPres,
                        proppantConc,
                        dProppantConc,
                        componentDens,
                        dComponentDens_dPres,
                        dComponentDens_dComponentConc,
                        gravCoef,
                        dens,
                        dDens_dPres,
                        dDens_dProppantConc,
                        dDens_dComponentConc,
                        visc,
                        dVisc_dPres,
                        dVisc_dProppantConc,
                        dVisc_dComponentConc,
                        fluidDensity,
                        dFluidDens_dPres,
                        dFluidDens_dComponentConc,
                        settlingFactor,
                        dSettlingFactor_dPres,
                        dSettlingFactor_dProppantConc,
                        dSettlingFactor_dComponentConc,
                        collisionFactor,
                        dCollisionFactor_dProppantConc,
                        isProppantMobile,
                        proppantPackVf,
                        aperture.toNestedViewConst(),
                        localMatrix,
                        localRhs );
  } );
}

void ProppantTransport::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  //  Apply Dirichlet BC for proppant concentration

  fsManager.apply( time_n + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::proppantConcentrationString(),
                   [&]( FieldSpecificationBase const & fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group & subRegion,
                        string const & )
  {
    arrayView1d< globalIndex const > const
    dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const
    proppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantConcentrationString() );

    arrayView1d< real64 const > const
    dProppantConc = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );

    fs.applyBoundaryConditionToSystem< FieldSpecificationEqual,
                                       parallelDevicePolicy<> >( lset,
                                                                 time_n + dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 rankOffset,
                                                                 localMatrix,
                                                                 localRhs,
                                                                 [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      return proppantConc[a] + dProppantConc[a];
    } );
  } );

  //  Apply Dirichlet BC for component concentration

  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     viewKeyStruct::proppantConcentrationString(),
                     [&]( FieldSpecificationBase const &,
                          string const & setName,
                          SortedArrayView< localIndex const > const &,
                          Group & subRegion,
                          string const & )
    {

      string const & subRegionName = subRegion.getName();
      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0, "Conflicting proppant boundary conditions on set " << setName );
      bcStatusMap[subRegionName][setName].resize( NC );
      bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

    } );

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     viewKeyStruct::componentConcentrationString(),
                     [&] ( FieldSpecificationBase const & fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group & subRegion,
                           string const & )
    {

      string const & subRegionName = subRegion.getName();
      localIndex const comp = fs.getComponent();

      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0, "Proppant boundary condition not prescribed on set '" << setName << "'" );
      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      fs.applyFieldValue< FieldSpecificationEqual >( targetSet,
                                                     time_n + dt,
                                                     subRegion,
                                                     viewKeyStruct::bcComponentConcentrationString() );

    } );

    bool bcConsistent = true;
    for( auto const & bcStatusEntryOuter : bcStatusMap )
    {
      for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          bcConsistent &= bcStatusEntryInner.second[ic];
          GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic
                                                                                                      << " on region '" << bcStatusEntryOuter.first << "',"
                                                                                                      << " set '" << bcStatusEntryInner.first << "'" );
        }
      }
    }

    GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

    fsManager.apply( time_n + dt,
                     domain,
                     "ElementRegions",
                     viewKeyStruct::proppantConcentrationString(),
                     [&] ( FieldSpecificationBase const &,
                           string const &,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group & subRegion,
                           string const & )
    {
      arrayView1d< integer const > const ghostRank =
        subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
      arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const > const compConc =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::componentConcentrationString() );
      arrayView2d< real64 const > const deltaCompConc =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );
      arrayView2d< real64 const > const bcCompConc =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString() );

      forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
      {
        localIndex const ei = targetSet[a];
        if( ghostRank[ei] >= 0 )
          return;

        globalIndex const dofIndex = dofNumber[ei];
        localIndex const localRow = dofIndex - rankOffset;
        real64 rhsValue;

        for( localIndex ic = 0; ic < NC; ++ic )
        {
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompConc[ei][ic],
                                                      compConc[ei][ic] + deltaCompConc[ei][ic] );
          localRhs[localRow + ic + 1] = rhsValue;
        }
      } );
    } );
  }
}

real64
ProppantTransport::calculateResidualNorm( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  localIndex const NDOF = m_numDofPerCell;

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::proppantConcentrationString() );

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const volume = subRegion.getElementVolume();

    RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        localIndex const lid = dofNumber[ei] - rankOffset;
        for( localIndex idof = 0; idof < NDOF; ++idof )
        {
          real64 const val = localRhs[lid] / volume[ei];
          localSum += val * val;
        }
      }
    } );

    localResidualNorm += localSum.get();
  } );

  // compute global residual norm
  real64 const globalResidualNorm = MpiWrapper::sum( localResidualNorm, MPI_COMM_GEOSX );
  return sqrt( globalResidualNorm );
}

void ProppantTransport::applySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::proppantConcentrationString(),
                               viewKeyStruct::deltaProppantConcentrationString(),
                               scalingFactor,
                               { m_numDofPerCell, 0, 1 } );


  if( m_numDofPerCell > 1 )
  {
    dofManager.addVectorToField( localSolution,
                                 viewKeyStruct::proppantConcentrationString(),
                                 viewKeyStruct::deltaComponentConcentrationString(),
                                 scalingFactor,
                                 { m_numDofPerCell, 1, m_numDofPerCell } );
  }

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaProppantConcentrationString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaComponentConcentrationString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    //UpdateState( subRegion );
    updateComponentDensity( subRegion, targetIndex );
  } );

}

void ProppantTransport::solveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void ProppantTransport::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dProppantConc =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaProppantConcentrationString() );
    arrayView2d< real64 > const & dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      dProppantConc[ei] = 0.0;
      for( localIndex c = 0; c < NC; ++c )
      {
        dComponentConc[ei][c] = 0.0;
      }
    } );

    updateState( subRegion, targetIndex );
  } );

}

void ProppantTransport::resetViews( MeshLevel & mesh )
{
  FlowSolverBase::resetViews( mesh );
  ElementRegionManager & elemManager = mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::pressureString() );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString() );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaPressureString() );
  m_deltaPressure.setName( getName() + "/accessors/" + viewKeyStruct::deltaPressureString() );

  m_proppantConcentration.clear();
  m_proppantConcentration = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::proppantConcentrationString() );
  m_proppantConcentration.setName( getName() + "/accessors/" + viewKeyStruct::proppantConcentrationString() );

  m_deltaProppantConcentration.clear();
  m_deltaProppantConcentration = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaProppantConcentrationString() );
  m_deltaProppantConcentration.setName( getName() + "/accessors/" + viewKeyStruct::deltaProppantConcentrationString() );

  m_cellBasedFlux.clear();
  m_cellBasedFlux = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::cellBasedFluxString() );
  m_cellBasedFlux.setName( getName() + "/accessors/" + viewKeyStruct::cellBasedFluxString() );

  m_proppantPackVolumeFraction.clear();
  m_proppantPackVolumeFraction = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::proppantPackVolumeFractionString() );
  m_proppantPackVolumeFraction.setName( getName() + "/accessors/" + viewKeyStruct::proppantPackVolumeFractionString() );

  m_proppantExcessPackVolume.clear();
  m_proppantExcessPackVolume = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::proppantExcessPackVolumeString() );
  m_proppantExcessPackVolume.setName( getName() + "/accessors/" + viewKeyStruct::proppantExcessPackVolumeString() );

  m_proppantLiftFlux.clear();
  m_proppantLiftFlux = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::proppantLiftFluxString() );
  m_proppantLiftFlux.setName( getName() + "/accessors/" + viewKeyStruct::proppantLiftFluxString() );

  m_isProppantMobile.clear();
  m_isProppantMobile = elemManager.constructArrayViewAccessor< integer, 1 >( viewKeyStruct::isProppantMobileString() );
  m_isProppantMobile.setName( getName() + "/accessors/" + viewKeyStruct::isProppantMobileString() );

  m_isProppantBoundaryElement.clear();
  m_isProppantBoundaryElement = elemManager.constructArrayViewAccessor< integer, 1 >( viewKeyStruct::isProppantBoundaryString() );
  m_isProppantBoundaryElement.setName( getName() + "/accessors/" + viewKeyStruct::isProppantBoundaryString() );

  {
    using keys = SlurryFluidBase::viewKeyStruct;

    m_density.clear();
    m_density = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::densityString(),
                                                                             targetRegionNames(),
                                                                             fluidModelNames() );
    m_density.setName( getName() + "/accessors/" + keys::densityString() );

    m_dDensity_dPressure.clear();
    m_dDensity_dPressure = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dDens_dPresString(),
                                                                                        targetRegionNames(),
                                                                                        fluidModelNames() );
    m_dDensity_dPressure.setName( getName() + "/accessors/" + keys::dDens_dPresString() );

    m_dDensity_dProppantConcentration.clear();
    m_dDensity_dProppantConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dDens_dProppantConcString(),
                                                                                                     targetRegionNames(),
                                                                                                     fluidModelNames() );
    m_dDensity_dProppantConcentration.setName( getName() + "/accessors/" + keys::dDens_dProppantConcString() );

    m_dDensity_dComponentConcentration.clear();
    m_dDensity_dComponentConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dDens_dCompConcString(),
                                                                                                      targetRegionNames(),
                                                                                                      fluidModelNames() );
    m_dDensity_dComponentConcentration.setName( getName() + "/accessors/" + keys::dDens_dCompConcString() );

    m_componentDensity.clear();
    m_componentDensity = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::componentDensityString(),
                                                                                      targetRegionNames(),
                                                                                      fluidModelNames() );
    m_componentDensity.setName( getName() + "/accessors/" + keys::componentDensityString() );

    m_dComponentDensity_dPressure.clear();
    m_dComponentDensity_dPressure = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dCompDens_dPresString(),
                                                                                                 targetRegionNames(),
                                                                                                 fluidModelNames() );
    m_dComponentDensity_dPressure.setName( getName() + "/accessors/" + keys::dCompDens_dPresString() );

    m_dComponentDensity_dComponentConcentration.clear();
    m_dComponentDensity_dComponentConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dCompDens_dCompConcString(),
                                                                                                               targetRegionNames(),
                                                                                                               fluidModelNames() );
    m_dComponentDensity_dComponentConcentration.setName( getName() + "/accessors/" + keys::dCompDens_dCompConcString() );

    m_fluidDensity.clear();
    m_fluidDensity = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::fluidDensityString(),
                                                                                  targetRegionNames(),
                                                                                  fluidModelNames() );
    m_fluidDensity.setName( getName() + "/accessors/" + keys::fluidDensityString() );

    m_dFluidDensity_dPressure.clear();
    m_dFluidDensity_dPressure = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dFluidDens_dPresString(),
                                                                                             targetRegionNames(),
                                                                                             fluidModelNames() );
    m_dFluidDensity_dPressure.setName( getName() + "/accessors/" + keys::dFluidDens_dPresString() );

    m_dFluidDensity_dComponentConcentration.clear();
    m_dFluidDensity_dComponentConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dFluidDens_dCompConcString(),
                                                                                                           targetRegionNames(),
                                                                                                           fluidModelNames() );
    m_dFluidDensity_dComponentConcentration.setName( getName() + "/accessors/" + keys::dFluidDens_dCompConcString() );

    m_fluidViscosity.clear();
    m_fluidViscosity = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::fluidViscosityString(),
                                                                                    targetRegionNames(),
                                                                                    fluidModelNames() );
    m_fluidViscosity.setName( getName() + "/accessors/" + keys::fluidViscosityString() );

    m_viscosity.clear();
    m_viscosity = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::viscosityString(),
                                                                               targetRegionNames(),
                                                                               fluidModelNames() );
    m_viscosity.setName( getName() + "/accessors/" + keys::viscosityString() );

    m_dViscosity_dPressure.clear();
    m_dViscosity_dPressure = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dVisc_dPresString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dViscosity_dPressure.setName( getName() + "/accessors/" + keys::dVisc_dPresString() );

    m_dViscosity_dProppantConcentration.clear();
    m_dViscosity_dProppantConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dVisc_dProppantConcString(),
                                                                                                       targetRegionNames(),
                                                                                                       fluidModelNames() );
    m_dViscosity_dProppantConcentration.setName( getName() + "/accessors/" + keys::dVisc_dProppantConcString() );

    m_dViscosity_dComponentConcentration.clear();
    m_dViscosity_dComponentConcentration = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dVisc_dCompConcString(),
                                                                                                        targetRegionNames(),
                                                                                                        fluidModelNames() );
    m_dViscosity_dComponentConcentration.setName( getName() + "/accessors/" + keys::dVisc_dCompConcString() );
  }

  {
    using keys = ParticleFluidBase::viewKeyStruct;

    m_settlingFactor.clear();
    m_settlingFactor =
      elemManager.constructMaterialArrayViewAccessor< real64, 1 >( keys::settlingFactorString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_settlingFactor.setName( getName() + "/accessors/" + keys::settlingFactorString() );

    m_dSettlingFactor_dPressure.clear();
    m_dSettlingFactor_dPressure =
      elemManager.constructMaterialArrayViewAccessor< real64, 1 >( keys::dSettlingFactor_dPressureString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_dSettlingFactor_dPressure.setName( getName() + "/accessors/" + keys::dSettlingFactor_dPressureString() );

    m_dSettlingFactor_dProppantConcentration.clear();
    m_dSettlingFactor_dProppantConcentration =
      elemManager.constructMaterialArrayViewAccessor< real64, 1 >( keys::dSettlingFactor_dProppantConcentrationString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_dSettlingFactor_dProppantConcentration.setName( getName() + "/accessors/" + keys::dSettlingFactor_dProppantConcentrationString() );

    m_dSettlingFactor_dComponentConcentration.clear();
    m_dSettlingFactor_dComponentConcentration =
      elemManager.constructMaterialArrayViewAccessor< real64, 2 >( keys::dSettlingFactor_dComponentConcentrationString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_dSettlingFactor_dComponentConcentration.setName( getName() + "/accessors/" + keys::dSettlingFactor_dComponentConcentrationString() );

    m_collisionFactor.clear();
    m_collisionFactor =
      elemManager.constructMaterialArrayViewAccessor< real64, 1 >( keys::collisionFactorString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_collisionFactor.setName( getName() + "/accessors/" + keys::collisionFactorString() );

    m_dCollisionFactor_dProppantConcentration.clear();
    m_dCollisionFactor_dProppantConcentration =
      elemManager.constructMaterialArrayViewAccessor< real64, 1 >( keys::dCollisionFactor_dProppantConcentrationString(),
                                                                   targetRegionNames(),
                                                                   proppantModelNames() );
    m_dCollisionFactor_dProppantConcentration.setName( getName() + "/accessors/" + keys::dCollisionFactor_dProppantConcentrationString() );
  }

  m_transTMultiplier.clear();
  m_transTMultiplier = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::transTMultiplierString() );
  m_transTMultiplier.setName( getName() + "/accessors/" + viewKeyStruct::transTMultiplierString() );
}



void ProppantTransport::updateCellBasedFlux( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const pres             = m_pressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const gravCoef         = m_gravCoef.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dens             = m_density.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const visc             = m_viscosity.toNestedViewConst();
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const
  aperture = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const proppantPackVf   = m_proppantPackVolumeFraction.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const transTMultiplier = m_transTMultiplier.toNestedViewConst();

  ElementRegionManager::ElementViewAccessor< arrayView2d< real64 > > const & cellBasedFluxAccessor =
    mesh.getElemManager().constructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::cellBasedFluxString() );

  FluxKernel::ElementView< arrayView2d< real64 > > const & cellBasedFlux = cellBasedFluxAccessor.toNestedView();

  fluxApprox.forAllStencils( mesh, [&]( auto const & stencil )
  {
    FluxKernel::launchCellBasedFluxCalculation( stencil,
                                                transTMultiplier,
                                                downVector,
                                                pres,
                                                gravCoef,
                                                dens,
                                                visc,
                                                aperture.toNestedViewConst(),
                                                proppantPackVf,
                                                cellBasedFlux );
  } );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::cellBasedFluxString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

}

void ProppantTransport::updateProppantPackVolume( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;

  R1Tensor downVector = gravityVector();
  LvArray::tensorOps::normalize< 3 >( downVector );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  // For data modified through an accessor, we must create the view accessor
  // every time in order to ensure the data gets properly touched on device

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantConc =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantConcentrationString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantPackVf =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantExcessPackV =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantExcessPackVolumeString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 > > const proppantLiftFlux =
    elemManager.constructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::proppantLiftFluxString() );

  ElementRegionManager::ElementViewAccessor< arrayView1d< real64 const > > const
  aperture = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );

  fluxApprox.forAllStencils( mesh, [&]( auto const & stencil )
  {
    ProppantPackVolumeKernel::launchProppantPackVolumeCalculation( stencil,
                                                                   dt,
                                                                   m_proppantDensity,
                                                                   m_proppantDiameter,
                                                                   m_maxProppantConcentration,
                                                                   downVector,
                                                                   m_criticalShieldsNumber,
                                                                   m_frictionCoefficient,
                                                                   m_settlingFactor.toNestedViewConst(),
                                                                   m_density.toNestedViewConst(),
                                                                   m_fluidDensity.toNestedViewConst(),
                                                                   m_fluidViscosity.toNestedViewConst(),
                                                                   m_isProppantMobile.toNestedViewConst(),
                                                                   m_isProppantBoundaryElement.toNestedViewConst(),
                                                                   aperture.toNestedViewConst(),
                                                                   m_volume.toNestedViewConst(),
                                                                   m_elemGhostRank.toNestedViewConst(),
                                                                   m_cellBasedFlux.toNestedViewConst(),
                                                                   proppantConc.toNestedView(),
                                                                   proppantPackVf.toNestedView(),
                                                                   proppantExcessPackV.toNestedView(),
                                                                   proppantLiftFlux.toNestedView() );
  } );

  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString() ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantPackVolumeFractionString() ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantExcessPackVolumeString() ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantLiftFluxString() ) );

    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
  }

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    updateProppantMobility( subRegion );
  } );


  fluxApprox.forAllStencils( mesh, [&]( auto const & stencil )
  {
    ProppantPackVolumeKernel::launchProppantPackVolumeUpdate( stencil,
                                                              downVector,
                                                              m_maxProppantConcentration,
                                                              m_isProppantMobile.toNestedViewConst(),
                                                              m_proppantExcessPackVolume.toNestedView(),
                                                              proppantConc.toNestedView(),
                                                              proppantPackVf.toNestedView() );
  } );


  {
    std::map< string, string_array > fieldNames;
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantConcentrationString() ) );
    fieldNames["elems"].emplace_back( string( viewKeyStruct::proppantPackVolumeFractionString() ) );

    CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
  }

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    updateProppantMobility( subRegion );
  } );

  // update poroMultiplier and transTMultiplier

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const proppantPackVfNew =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::proppantPackVolumeFractionString() );
    arrayView1d< real64 const > const aperture0 =
      subRegion.getReference< array1d< real64 > >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );

    arrayView1d< real64 > const poroMultiplier =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultiplierString() );
    arrayView2d< real64 > const transTMultiplier =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::transTMultiplierString() );

    real64 const maxProppantConcentration = m_maxProppantConcentration;
    real64 const proppantPackPermeability = m_proppantPackPermeability;

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      poroMultiplier[ei] = 1.0 - maxProppantConcentration * proppantPackVfNew[ei];

      //K0 horizontal
      transTMultiplier[ei][0] = ( 1.0 - proppantPackVfNew[ei] )
                                + proppantPackVfNew[ei] * proppantPackPermeability * 12.0 / ( aperture0[ei] * aperture0[ei] );

      //K1 vertical
      transTMultiplier[ei][1] = 1.0 /
                                ( proppantPackVfNew[ei] * aperture0[ei] * aperture0[ei] / 12.0 / proppantPackPermeability
                                  + ( 1.0 - proppantPackVfNew[ei] ) );
    } );
  } );
}


REGISTER_CATALOG_ENTRY( SolverBase, ProppantTransport, string const &, Group * const )
} /* namespace geosx */
