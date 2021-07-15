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
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

FlowSolverBase::FlowSolverBase( string const & name,
                                Group * const parent ):
  SolverBase( name, parent ),
  m_fluidModelNames(),
  m_solidModelNames(),
  m_poroElasticFlag( 0 ),
  m_coupledWellsFlag( 0 ),
  m_numDofPerCell( 0 ),
  m_derivativeFluxResidual_dAperture(),
  m_fluxEstimate(),
  m_elemGhostRank(),
  m_volume(),
  m_gravCoef(),
  m_porosityRef(),
  m_elementArea(),
  m_elementAperture0(),
  m_elementAperture()
{
  this->registerWrapper( viewKeyStruct::discretizationString(), &m_discretizationName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of discretization object to use for this solver." );

  this->registerWrapper( viewKeyStruct::fluidNamesString(), &m_fluidModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of fluid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::solidNamesString(), &m_solidModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of solid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::permeabilityNamesString(), &m_permeabilityModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of permeability constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::inputFluxEstimateString(), &m_fluxEstimate ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Initial estimate of the input flux used only for residual scaling. This should be "
                    "essentially equivalent to the input flux * dt." );

  this->registerWrapper( viewKeyStruct::meanPermCoeffString(), &m_meanPermCoeff ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Coefficient to move between harmonic mean (1.0) and arithmetic mean (0.0) for the "
                    "calculation of permeability between elements." );

}

void FlowSolverBase::registerDataOnMesh( Group & meshBodies )
{
  SolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & mesh = meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const,
                                    ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).
        setApplyDefaultValue( 0.0 );
    } );

    ElementRegionManager & elemManager = mesh.getElemManager();

    elemManager.forElementSubRegionsComplete< SurfaceElementSubRegion >( [&]( localIndex const,
                                                                              localIndex const,
                                                                              ElementRegionBase & region,
                                                                              SurfaceElementSubRegion & subRegion )
    {
      SurfaceElementRegion & faceRegion = dynamicCast< SurfaceElementRegion & >( region );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString() ).
        setApplyDefaultValue( 1.0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).
        setApplyDefaultValue( 0.0 );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::aperture0String() ).
        setDefaultValue( faceRegion.getDefaultAperture() );

      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::effectiveApertureString() ).
        setApplyDefaultValue( faceRegion.getDefaultAperture() ).
        setPlotLevel( PlotLevel::LEVEL_0 );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString() ).setApplyDefaultValue( 0.0 );
  } );
}

void FlowSolverBase::postProcessInput()
{
  SolverBase::postProcessInput();
  checkModelNames( m_fluidModelNames, viewKeyStruct::fluidNamesString() );
  checkModelNames( m_solidModelNames, viewKeyStruct::solidNamesString() );
}

void FlowSolverBase::initializePreSubGroups()
{
  SolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  // Validate solid models in regions (fluid models are validated by derived classes)
  domain.getMeshBodies().forSubGroups< MeshBody >( [&] ( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );
    validateModelMapping( meshLevel.getElemManager(), m_solidModelNames );
  } );

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();

  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
    array1d< string > & stencilTargetRegions = fluxApprox.targetRegions();
    std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(), stencilTargetRegions.end() );
    for( auto const & targetRegion : targetRegionNames() )
    {
      stencilTargetRegionsSet.insert( targetRegion );
    }

    stencilTargetRegions.clear();
    for( auto const & targetRegion : stencilTargetRegionsSet )
    {
      stencilTargetRegions.emplace_back( targetRegion );
    }
  }
}

void FlowSolverBase::initializePostInitialConditionsPreSubGroups()
{
  SolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );;
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  resetViews( mesh );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  precomputeData( mesh );
}

void FlowSolverBase::precomputeData( MeshLevel & mesh )
{
  FaceManager & faceManager = mesh.getFaceManager();
  real64 const gravVector[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( gravityVector() );

  forTargetSubRegions( mesh, [&]( localIndex const,
                                  ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 const > const elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const gravityCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

    forAll< parallelHostPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      gravityCoef[ ei ] = LvArray::tensorOps::AiBi< 3 >( elemCenter[ ei ], gravVector );
    } );
  } );

  {
    arrayView2d< real64 const > const faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const gravityCoef =
      faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

    forAll< parallelHostPolicy >( faceManager.size(), [=] ( localIndex const kf )
    {
      gravityCoef[ kf ] = LvArray::tensorOps::AiBi< 3 >( faceCenter[ kf ], gravVector );
    } );
  }
}

FlowSolverBase::~FlowSolverBase() = default;

void FlowSolverBase::resetViews( MeshLevel & mesh )
{
  ElementRegionManager const & elemManager = mesh.getElemManager();

  m_elemGhostRank.clear();
  m_elemGhostRank = elemManager.constructArrayViewAccessor< integer, 1 >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
  m_elemGhostRank.setName( getName() + "/accessors/" + ObjectManagerBase::viewKeyStruct::ghostRankString() );

  m_volume.clear();
  m_volume = elemManager.constructArrayViewAccessor< real64, 1 >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );
  m_volume.setName( getName() + "/accessors/" + ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

  m_gravCoef.clear();
  m_gravCoef = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::gravityCoefString() );
  m_gravCoef.setName( getName() + "/accessors/" + viewKeyStruct::gravityCoefString() );

  m_porosityRef.clear();
  m_porosityRef = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::referencePorosityString() );
  m_porosityRef.setName( getName() + "/accessors/" + viewKeyStruct::referencePorosityString() );

  m_elementArea.clear();
  m_elementArea = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementAreaString() );
  m_elementArea.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::elementAreaString() );

  m_elementAperture.clear();
  m_elementAperture = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::elementApertureString() );
  m_elementAperture.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::elementApertureString() );

  m_elementAperture0.clear();
  m_elementAperture0 = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::aperture0String() );
  m_elementAperture0.setName( getName() + "/accessors/" + viewKeyStruct::aperture0String() );

  m_effectiveAperture.clear();
  m_effectiveAperture = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::effectiveApertureString() );
  m_effectiveAperture.setName( getName() + "/accessors/" + viewKeyStruct::effectiveApertureString() );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  m_elementSeparationCoefficient.clear();
  m_elementSeparationCoefficient = elemManager.constructArrayViewAccessor< real64, 1 >( FaceElementSubRegion::viewKeyStruct::separationCoeffString() );
  m_elementSeparationCoefficient.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::separationCoeffString() );

  m_element_dSeparationCoefficient_dAperture.clear();
  m_element_dSeparationCoefficient_dAperture = elemManager.constructArrayViewAccessor< real64, 1 >(
    FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
  m_element_dSeparationCoefficient_dAperture.setName( getName() + "/accessors/" + FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString() );
#endif
}

std::vector< string > FlowSolverBase::getConstitutiveRelations( string const & regionName ) const
{

  localIndex const regionIndex = this->targetRegionIndex( regionName );

  std::vector< string > rval{ m_solidModelNames[regionIndex], m_fluidModelNames[regionIndex] };

  return rval;
}

void FlowSolverBase::setUpDflux_dApertureMatrix( DomainPartition & GEOSX_UNUSED_PARAM( domain ),
                                                 DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
                                                 CRSMatrix< real64, globalIndex > & GEOSX_UNUSED_PARAM( localMatrix ) )
{
  GEOSX_ERROR( "FlowSolverBase::setUpDfluxDapertureMatrix. Should be overridden." );
}


} // namespace geosx
