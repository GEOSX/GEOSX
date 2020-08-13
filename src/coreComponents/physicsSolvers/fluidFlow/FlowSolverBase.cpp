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
 * @file FlowSolverBase.cpp
 */

#include "FlowSolverBase.hpp"

#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;

FlowSolverBase::FlowSolverBase( std::string const & name, Group * const parent ) :
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
  this
    ->registerWrapper( viewKeyStruct::discretizationString, &m_discretizationName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Name of discretization object to use for this solver." );

  this->registerWrapper( viewKeyStruct::fluidNamesString, &m_fluidModelNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Names of fluid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::solidNamesString, &m_solidModelNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Names of solid constitutive models for each region." );

  this->registerWrapper( viewKeyStruct::inputFluxEstimateString, &m_fluxEstimate )
    ->setApplyDefaultValue( 1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Initial estimate of the input flux used only for residual scaling. This "
      "should be "
      "essentially equivalent to the input flux * dt." );

  this->registerWrapper( viewKeyStruct::meanPermCoeffString, &m_meanPermCoeff )
    ->setApplyDefaultValue( 1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Coefficient to move between harmonic mean (1.0) and arithmetic mean "
      "(0.0) for the "
      "calculation of permeability between elements." );
}

void
FlowSolverBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  SolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & subgroup : MeshBodies->GetSubGroups() )
  {
    MeshBody & meshBody = *subgroup.second->group_cast< MeshBody * >();
    MeshLevel & mesh = *meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
      subRegion
        .registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString )
        ->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion
        .registerWrapper< array1d< R1Tensor > >( viewKeyStruct::permeabilityString )
        ->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion
        .registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )
        ->setApplyDefaultValue( 0.0 );
    } );

    ElementRegionManager * const elemManager = mesh.getElemManager();

    elemManager->forElementSubRegionsComplete< FaceElementSubRegion >(
      [&]( localIndex const,
           localIndex const,
           ElementRegionBase & region,
           FaceElementSubRegion & subRegion ) {
        FaceElementRegion & faceRegion = dynamicCast< FaceElementRegion & >( region );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::referencePorosityString )
          ->setApplyDefaultValue( 1.0 );

        subRegion
          .registerWrapper< array1d< R1Tensor > >( viewKeyStruct::permeabilityString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )
          ->setApplyDefaultValue( 0.0 );
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::aperture0String )
          ->setDefaultValue( faceRegion.getDefaultAperture() );
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::effectiveApertureString )
          ->setApplyDefaultValue(
            subRegion
              .getWrapper< array1d< real64 > >(
                FaceElementSubRegion::viewKeyStruct::elementApertureString )
              ->getDefaultValue() )
          ->setPlotLevel( PlotLevel::LEVEL_0 );
      } );

    FaceManager * const faceManager = mesh.getFaceManager();
    faceManager
      ->registerWrapper< array1d< real64 > >( viewKeyStruct::gravityCoefString )
      ->setApplyDefaultValue( 0.0 );
  }
}

void
FlowSolverBase::PostProcessInput()
{
  SolverBase::PostProcessInput();
  CheckModelNames( m_fluidModelNames, viewKeyStruct::fluidNamesString );
  CheckModelNames( m_solidModelNames, viewKeyStruct::solidNamesString );
}

void
FlowSolverBase::InitializePreSubGroups( Group * const rootGroup )
{
  SolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain =
    rootGroup->GetGroup< DomainPartition >( keys::domain );

  // Validate solid models in regions (fluid models are validated by derived classes)
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel =
      *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping( *meshLevel.getElemManager(), m_solidModelNames );
  }

  // fill stencil targetRegions
  NumericalMethodsManager & numericalMethodManager =
    domain->getNumericalMethodManager();

  FiniteVolumeManager & fvManager =
    numericalMethodManager.getFiniteVolumeManager();

  FluxApproximationBase & fluxApprox =
    fvManager.getFluxApproximation( m_discretizationName );
  array1d< string > & stencilTargetRegions = fluxApprox.targetRegions();
  std::set< string > stencilTargetRegionsSet( stencilTargetRegions.begin(),
                                              stencilTargetRegions.end() );
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

void
FlowSolverBase::InitializePostInitialConditions_PreSubGroups(
  Group * const rootGroup )
{
  SolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  ResetViews( mesh );

  // Precompute solver-specific constant data (e.g. gravity-depth)
  PrecomputeData( mesh );
}

void
FlowSolverBase::PrecomputeData( MeshLevel & mesh )
{
  FaceManager & faceManager = *mesh.getFaceManager();
  R1Tensor const gravVector = gravityVector();

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    arrayView2d< real64 const > const & elemCenter = subRegion.getElementCenter();

    arrayView1d< real64 > const & gravityCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    forAll< parallelHostPolicy >( subRegion.size(), [=]( localIndex const ei ) {
      gravityCoef[ei] = elemCenter( ei, 0 ) * gravVector[0] +
        elemCenter( ei, 1 ) * gravVector[1] + elemCenter( ei, 2 ) * gravVector[2];
    } );
  } );

  {
    arrayView2d< real64 const > const & faceCenter = faceManager.faceCenter();

    arrayView1d< real64 > const & gravityCoef =
      faceManager.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString );

    forAll< parallelHostPolicy >( faceManager.size(), [=]( localIndex const kf ) {
      // TODO change to LvArray::tensorOps::AiBi once gravVector is a c-array.
      gravityCoef[kf] = faceCenter[kf][0] * gravVector[0];
      gravityCoef[kf] += faceCenter[kf][1] * gravVector[1];
      gravityCoef[kf] += faceCenter[kf][2] * gravVector[2];
    } );
  }
}

FlowSolverBase::~FlowSolverBase() = default;

void
FlowSolverBase::ResetViews( MeshLevel & mesh )
{
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  m_elemGhostRank.clear();
  m_elemGhostRank = elemManager.ConstructArrayViewAccessor< integer, 1 >(
    ObjectManagerBase::viewKeyStruct::ghostRankString );
  m_elemGhostRank.setName( getName() + "/accessors/" +
                           ObjectManagerBase::viewKeyStruct::ghostRankString );

  m_volume.clear();
  m_volume = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    ElementSubRegionBase::viewKeyStruct::elementVolumeString );
  m_volume.setName( getName() + "/accessors/" +
                    ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  m_gravCoef.clear();
  m_gravCoef = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::gravityCoefString );
  m_gravCoef.setName( getName() + "/accessors/" + viewKeyStruct::gravityCoefString );

  m_porosityRef.clear();
  m_porosityRef = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::referencePorosityString );
  m_porosityRef.setName( getName() + "/accessors/" +
                         viewKeyStruct::referencePorosityString );

  m_elementArea.clear();
  m_elementArea = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    FaceElementSubRegion::viewKeyStruct::elementAreaString );
  m_elementArea.setName( getName() + "/accessors/" +
                         FaceElementSubRegion::viewKeyStruct::elementAreaString );

  m_elementAperture.clear();
  m_elementAperture = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    FaceElementSubRegion::viewKeyStruct::elementApertureString );
  m_elementAperture.setName(
    getName() + "/accessors/" +
    FaceElementSubRegion::viewKeyStruct::elementApertureString );

  m_elementAperture0.clear();
  m_elementAperture0 = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::aperture0String );
  m_elementAperture0.setName( getName() + "/accessors/" +
                              viewKeyStruct::aperture0String );

  m_effectiveAperture.clear();
  m_effectiveAperture = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::effectiveApertureString );
  m_effectiveAperture.setName( getName() + "/accessors/" +
                               viewKeyStruct::effectiveApertureString );

#ifdef GEOSX_USE_SEPARATION_COEFFICIENT
  m_elementSeparationCoefficient.clear();
  m_elementSeparationCoefficient =
    elemManager.ConstructArrayViewAccessor< real64, 1 >(
      FaceElementSubRegion::viewKeyStruct::separationCoeffString );
  m_elementSeparationCoefficient.setName(
    getName() + "/accessors/" +
    FaceElementSubRegion::viewKeyStruct::separationCoeffString );

  m_element_dSeparationCoefficient_dAperture.clear();
  m_element_dSeparationCoefficient_dAperture =
    elemManager.ConstructArrayViewAccessor< real64, 1 >(
      FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
  m_element_dSeparationCoefficient_dAperture.setName(
    getName() + "/accessors/" +
    FaceElementSubRegion::viewKeyStruct::dSeparationCoeffdAperString );
#endif
}

std::vector< string >
FlowSolverBase::getConstitutiveRelations(
  string const & regionName ) const
{
  localIndex const regionIndex = this->targetRegionIndex( regionName );

  std::vector< string > rval { m_solidModelNames[regionIndex],
                               m_fluidModelNames[regionIndex] };

  return rval;
}

}  // namespace geosx
