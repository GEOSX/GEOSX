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
 * @file CompositionalMultiphaseFlow.cpp
 */

#include "CompositionalMultiphaseFlow.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/relativePermeabilitySelector.hpp"
#include "dataRepository/Group.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseFlowKernels;

static constexpr real64 minDensForDivision = 1e-10;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          Group * const parent ) :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_capPressureFlag( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{
  //START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::temperatureString, &m_temperature )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useMassFlagString, &m_useMass )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNamesString, &m_relPermModelNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription(
      "Name of the relative permeability constitutive model to use" );

  this
    ->registerWrapper( viewKeyStruct::capPressureNamesString,
                       &m_capPressureModelNames )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Name of the capillary pressure constitutive model to use" );

  this
    ->registerWrapper( viewKeyStruct::maxCompFracChangeString, &m_maxCompFracChange )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setApplyDefaultValue( 1.0 )
    ->setDescription(
      "Maximum (absolute) change in a component fraction between two Newton "
      "iterations" );

  this
    ->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString,
                       &m_allowCompDensChopping )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setApplyDefaultValue( 1 )
    ->setDescription(
      "Flag indicating whether local (cell-wise) chopping of negative "
      "compositions is allowed" );

  m_linearSolverParameters.get().mgr.strategy = "CompositionalMultiphaseFlow";
}

void
CompositionalMultiphaseFlow::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();
  CheckModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString );
  m_capPressureFlag = CheckModelNames( m_capPressureModelNames,
                                       viewKeyStruct::capPressureNamesString,
                                       true );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange,
                         1.0,
                         "The maximum absolute change in component fraction "
                         "must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange,
                         0.0,
                         "The maximum absolute change in component fraction "
                         "must larger or equal to 0.0" );
}

void
CompositionalMultiphaseFlow::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel =
      *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions(
      meshLevel,
      [&]( localIndex const, ElementSubRegionBase & elementSubRegion ) {
        elementSubRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array1d< real64 > >(
          viewKeyStruct::bcPressureString );

        elementSubRegion
          .registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompDensityString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion
          .registerWrapper< array2d< real64 > >(
            viewKeyStruct::deltaGlobalCompDensityString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion
          .registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompFractionString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion
          .registerWrapper< array3d< real64 > >(
            viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion
          .registerWrapper< array2d< real64 > >(
            viewKeyStruct::phaseVolumeFractionString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion
          .registerWrapper< array2d< real64 > >(
            viewKeyStruct::dPhaseVolumeFraction_dPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion
          .registerWrapper< array3d< real64 > >(
            viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion
          .registerWrapper< array2d< real64 > >( viewKeyStruct::phaseMobilityString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion
          .registerWrapper< array2d< real64 > >(
            viewKeyStruct::dPhaseMobility_dPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion
          .registerWrapper< array3d< real64 > >(
            viewKeyStruct::dPhaseMobility_dGlobalCompDensityString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionOldString );
        elementSubRegion.registerWrapper< array2d< real64 > >(
          viewKeyStruct::phaseDensityOldString );
        elementSubRegion.registerWrapper< array3d< real64 > >(
          viewKeyStruct::phaseComponentFractionOldString );
        elementSubRegion.registerWrapper< array1d< real64 > >(
          viewKeyStruct::porosityOldString );
      } );
  }
}

namespace
{
template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void
CompareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG(
    lhs.numFluidPhases(),
    rhs.numFluidPhases(),
    "Mismatch in number of phases between constitutive models "
      << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOSX_ERROR_IF_NE_MSG( lhs.phaseNames()[ip],
                           rhs.phaseNames()[ip],
                           "Mismatch in phase names between constitutive models "
                             << lhs.getName() << " and " << rhs.getName() );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void
CompareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG(
    lhs.numFluidComponents(),
    rhs.numFluidComponents(),
    "Mismatch in number of components between constitutive models "
      << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOSX_ERROR_IF_NE_MSG(
      lhs.componentNames()[ic],
      rhs.componentNames()[ic],
      "Mismatch in component names between constitutive models "
        << lhs.getName() << " and " << rhs.getName() );
  }
}

}  // namespace

void
CompositionalMultiphaseFlow::ValidateConstitutiveModels(
  constitutive::ConstitutiveManager const & cm ) const
{
  MultiFluidBase const & fluid0 =
    *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  for( localIndex i = 1; i < m_fluidModelNames.size(); ++i )
  {
    MultiFluidBase const & fluid =
      *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[i] );
    CompareMultiphaseModels( fluid, fluid0 );
    CompareMulticomponentModels( fluid, fluid0 );
  }

  RelativePermeabilityBase const & relPerm0 =
    *cm.GetConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[0] );
  CompareMultiphaseModels( relPerm0, fluid0 );

  for( localIndex i = 1; i < m_relPermModelNames.size(); ++i )
  {
    RelativePermeabilityBase const & relPerm =
      *cm.GetConstitutiveRelation< RelativePermeabilityBase >(
        m_relPermModelNames[i] );
    CompareMultiphaseModels( relPerm, relPerm0 );
  }

  if( m_capPressureFlag )
  {
    CapillaryPressureBase const & capPres0 =
      *cm.GetConstitutiveRelation< CapillaryPressureBase >(
        m_capPressureModelNames[0] );
    CompareMultiphaseModels( capPres0, fluid0 );

    for( localIndex i = 1; i < m_capPressureModelNames.size(); ++i )
    {
      CapillaryPressureBase const & capPres =
        *cm.GetConstitutiveRelation< CapillaryPressureBase >(
          m_capPressureModelNames[i] );
      CompareMultiphaseModels( capPres, capPres0 );
    }
  }
}

void
CompositionalMultiphaseFlow::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * const domain =
    rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager const & cm = *domain->getConstitutiveManager();

  // 1. Set key dimensions of the problem
  MultiFluidBase const & fluid0 =
    *cm.GetConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  // 2. Validate various models against each other (must have same phases and components)
  ValidateConstitutiveModels( cm );

  // 3. Resize all fields as necessary, validate constitutive models in regions
  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel & meshLevel =
      *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    ResizeFields( meshLevel );

    ValidateModelMapping< MultiFluidBase >( *meshLevel.getElemManager(),
                                            m_fluidModelNames );
    ValidateModelMapping< RelativePermeabilityBase >( *meshLevel.getElemManager(),
                                                      m_relPermModelNames );
    if( m_capPressureFlag )
    {
      ValidateModelMapping< CapillaryPressureBase >( *meshLevel.getElemManager(),
                                                     m_capPressureModelNames );
    }
  }
}

void
CompositionalMultiphaseFlow::ResizeFields( MeshLevel & meshLevel ) const
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forTargetSubRegions( meshLevel, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString )
      .resizeDimension< 1 >( NC );
    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString )
      .resizeDimension< 1 >( NC );

    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString )
      .resizeDimension< 1 >( NC );
    subRegion
      .getReference< array3d< real64 > >(
        viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString )
      .resizeDimension< 1, 2 >( NC, NC );

    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array2d< real64 > >(
        viewKeyStruct::dPhaseVolumeFraction_dPressureString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array3d< real64 > >(
        viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString )
      .resizeDimension< 1, 2 >( NP, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array3d< real64 > >(
        viewKeyStruct::dPhaseMobility_dGlobalCompDensityString )
      .resizeDimension< 1, 2 >( NP, NC );

    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array2d< real64 > >( viewKeyStruct::phaseDensityOldString )
      .resizeDimension< 1 >( NP );
    subRegion
      .getReference< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString )
      .resizeDimension< 1, 2 >( NP, NC );
  } );
}

void
CompositionalMultiphaseFlow::UpdateComponentFraction( Group & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & compFrac = dataGroup.getReference< array2d< real64 > >(
    viewKeyStruct::globalCompFractionString );

  arrayView3d< real64 > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & compDens =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::deltaGlobalCompDensityString );

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    dataGroup.size(),
                                                    compDens,
                                                    dCompDens,
                                                    compFrac,
                                                    dCompFrac_dCompDens );
}

void
CompositionalMultiphaseFlow::UpdatePhaseVolumeFraction(
  Group & dataGroup,
  localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::phaseVolumeFractionString );

  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 > const & dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  // inputs

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  arrayView2d< real64 const > const & compDens =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::globalCompDensityString );

  arrayView2d< real64 const > const & dCompDens =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::deltaGlobalCompDensityString );

  MultiFluidBase const & fluid =
    GetConstitutiveModel< MultiFluidBase >( dataGroup,
                                            m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres =
    fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp =
    fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres =
    fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp =
    fluid.dPhaseDensity_dGlobalCompFraction();

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents,
                                                      m_numPhases,
                                                      dataGroup.size(),
                                                      compDens,
                                                      dCompDens,
                                                      dCompFrac_dCompDens,
                                                      phaseDens,
                                                      dPhaseDens_dPres,
                                                      dPhaseDens_dComp,
                                                      phaseFrac,
                                                      dPhaseFrac_dPres,
                                                      dPhaseFrac_dComp,
                                                      phaseVolFrac,
                                                      dPhaseVolFrac_dPres,
                                                      dPhaseVolFrac_dComp );
}

void
CompositionalMultiphaseFlow::UpdatePhaseMobility( Group & dataGroup,
                                                  localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseMob =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString );

  arrayView2d< real64 > const & dPhaseMob_dPres =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::dPhaseMobility_dPressureString );

  arrayView3d< real64 > const & dPhaseMob_dComp =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  // inputs

  arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  arrayView3d< real64 const > const & dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >(
      viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  MultiFluidBase const & fluid =
    GetConstitutiveModel< MultiFluidBase >( dataGroup,
                                            m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres =
    fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp =
    fluid.dPhaseDensity_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseVisc = fluid.phaseViscosity();
  arrayView3d< real64 const > const & dPhaseVisc_dPres =
    fluid.dPhaseViscosity_dPressure();
  arrayView4d< real64 const > const & dPhaseVisc_dComp =
    fluid.dPhaseViscosity_dGlobalCompFraction();

  RelativePermeabilityBase const & relperm =
    GetConstitutiveModel< RelativePermeabilityBase >(
      dataGroup,
      m_relPermModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseRelPerm = relperm.phaseRelPerm();
  arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac =
    relperm.dPhaseRelPerm_dPhaseVolFraction();

  KernelLaunchSelector2< PhaseMobilityKernel >( m_numComponents,
                                                m_numPhases,
                                                dataGroup.size(),
                                                dCompFrac_dCompDens,
                                                phaseDens,
                                                dPhaseDens_dPres,
                                                dPhaseDens_dComp,
                                                phaseVisc,
                                                dPhaseVisc_dPres,
                                                dPhaseVisc_dComp,
                                                phaseRelPerm,
                                                dPhaseRelPerm_dPhaseVolFrac,
                                                dPhaseVolFrac_dPres,
                                                dPhaseVolFrac_dComp,
                                                phaseMob,
                                                dPhaseMob_dPres,
                                                dPhaseMob_dComp );
}

void
CompositionalMultiphaseFlow::UpdateFluidModel( Group & dataGroup,
                                               localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView2d< real64 const > const compFrac =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::globalCompFractionString );

  MultiFluidBase & fluid =
    GetConstitutiveModel< MultiFluidBase >( dataGroup,
                                            m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid ) {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
      castedFluid.createKernelWrapper();

    // MultiFluid models are not thread-safe or device-capable yet
    FluidUpdateKernel::Launch< serialPolicy >( dataGroup.size(),
                                               fluidWrapper,
                                               pres,
                                               dPres,
                                               m_temperature,
                                               compFrac );
  } );
}

void
CompositionalMultiphaseFlow::UpdateSolidModel( Group & dataGroup,
                                               localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase & solid =
    GetConstitutiveModel< ConstitutiveBase >( dataGroup,
                                              m_solidModelNames[targetIndex] );

  arrayView1d< real64 const > const & pres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  solid.StateUpdateBatchPressure( pres, dPres );
}

void
CompositionalMultiphaseFlow::UpdateRelPermModel( Group & dataGroup,
                                                 localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< real64 const > const & phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >(
      viewKeyStruct::phaseVolumeFractionString );

  RelativePermeabilityBase & relPerm =
    GetConstitutiveModel< RelativePermeabilityBase >(
      dataGroup,
      m_relPermModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( relPerm, [&]( auto & castedRelPerm ) {
    typename TYPEOFREF( castedRelPerm )::KernelWrapper relPermWrapper =
      castedRelPerm.createKernelWrapper();

    RelativePermeabilityUpdateKernel::Launch< parallelDevicePolicy<> >(
      dataGroup.size(),
      relPermWrapper,
      phaseVolFrac );
  } );
}

void
CompositionalMultiphaseFlow::UpdateCapPressureModel(
  Group & dataGroup,
  localIndex const targetIndex ) const
{
  if( m_capPressureFlag )
  {
    arrayView2d< real64 const > const & phaseVolFrac =
      dataGroup.getReference< array2d< real64 > >(
        viewKeyStruct::phaseVolumeFractionString );

    CapillaryPressureBase & capPressure =
      GetConstitutiveModel< CapillaryPressureBase >(
        dataGroup,
        m_capPressureModelNames[targetIndex] );

    constitutive::constitutiveUpdatePassThru( capPressure, [&]( auto & castedCapPres ) {
      typename TYPEOFREF( castedCapPres )::KernelWrapper capPresWrapper =
        castedCapPres.createKernelWrapper();

      CapillaryPressureUpdateKernel::Launch< parallelDevicePolicy<> >(
        dataGroup.size(),
        capPresWrapper,
        phaseVolFrac );
    } );
  }
}

void
CompositionalMultiphaseFlow::UpdateState( Group & dataGroup,
                                          localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateComponentFraction( dataGroup );
  UpdateFluidModel( dataGroup, targetIndex );
  UpdatePhaseVolumeFraction( dataGroup, targetIndex );
  UpdateSolidModel( dataGroup, targetIndex );
  UpdateRelPermModel( dataGroup, targetIndex );
  UpdatePhaseMobility( dataGroup, targetIndex );
  UpdateCapPressureModel( dataGroup, targetIndex );
}

void
CompositionalMultiphaseFlow::InitializeFluidState( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      // 1. Assume global component fractions have been prescribed.
      // Initialize constitutive state to get fluid density.
      UpdateFluidModel( subRegion, targetIndex );

      // 2. Back-calculate global component densities from fractions and total fluid density
      // in order to initialize the primary solution variables
      MultiFluidBase const & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                fluidModelNames()[targetIndex] );
      arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

      arrayView2d< real64 const > const & compFrac =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::globalCompFractionString );
      arrayView2d< real64 > const & compDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::globalCompDensityString );

      forAll< parallelDevicePolicy<> >( subRegion.size(),
                                        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                          for( localIndex ic = 0; ic < NC; ++ic )
                                          {
                                            compDens[ei][ic] =
                                              totalDens[ei][0] * compFrac[ei][ic];
                                          }
                                        } );

      // 3. Update dependent state quantities
      UpdatePhaseVolumeFraction( subRegion, targetIndex );
      UpdateSolidModel( subRegion, targetIndex );
      UpdateRelPermModel( subRegion, targetIndex );
      UpdatePhaseMobility( subRegion, targetIndex );
      UpdateCapPressureModel( subRegion, targetIndex );
    } );
}

void
CompositionalMultiphaseFlow::InitializePostInitialConditions_PreSubGroups(
  Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::pressureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::globalCompDensityString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );

  // set mass fraction flag on fluid models
  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      MultiFluidBase & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                m_fluidModelNames[targetIndex] );
      fluid.setMassFlag( m_useMass );
    } );

  // Initialize primary variables from applied initial conditions
  InitializeFluidState( mesh );
}

real64
CompositionalMultiphaseFlow::SolverStep( real64 const & time_n,
                                         real64 const & dt,
                                         integer const cycleNumber,
                                         DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

  ImplicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void
CompositionalMultiphaseFlow::BackupFields( MeshLevel & mesh ) const
{
  // backup some fields used in time derivative approximation
  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & poroRef =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::referencePorosityString );
      arrayView2d< real64 const > const & phaseVolFrac =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionString );

      MultiFluidBase const & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                fluidModelNames()[targetIndex] );
      arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
      arrayView4d< real64 const > const & phaseCompFrac = fluid.phaseCompFraction();

      ConstitutiveBase const & solid =
        GetConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
      arrayView2d< real64 const > const & pvMult =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );

      arrayView2d< real64 > const & phaseDensOld =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseDensityOldString );
      arrayView2d< real64 > const & phaseVolFracOld =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionOldString );
      arrayView3d< real64 > const & phaseCompFracOld =
        subRegion.getReference< array3d< real64 > >(
          viewKeyStruct::phaseComponentFractionOldString );
      arrayView1d< real64 > const & poroOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );

      localIndex const NC = m_numComponents;
      localIndex const NP = m_numPhases;

      forAll< parallelDevicePolicy<> >(
        subRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          if( elemGhostRank[ei] >= 0 ) return;

          for( localIndex ip = 0; ip < NP; ++ip )
          {
            phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
            phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];

            for( localIndex ic = 0; ic < NC; ++ic )
            {
              phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
            }
          }
          poroOld[ei] = poroRef[ei] * pvMult[ei][0];
        } );
    } );
}

void
CompositionalMultiphaseFlow::ImplicitStepSetup(
  real64 const & GEOSX_UNUSED_PARAM( time_n ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  // bind the stored views to the current domain
  ResetViews( mesh );

  // set deltas to zero and recompute dependent quantities
  ResetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  BackupFields( mesh );
}

void
CompositionalMultiphaseFlow::SetupDofs( DomainPartition const & domain,
                                        DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::dofFieldString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager =
    domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager =
    numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox =
    fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::dofFieldString, fluxApprox );
}

void
CompositionalMultiphaseFlow::AssembleSystem(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  AssembleAccumulationTerms( domain, dofManager, localMatrix, localRhs );

  AssembleFluxTerms( dt, domain, dofManager, localMatrix, localRhs );

  AssembleVolumeBalanceTerms( domain, dofManager, localMatrix, localRhs );
}

void
CompositionalMultiphaseFlow::AssembleAccumulationTerms(
  DomainPartition const & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & porosityRef =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::referencePorosityString );

      arrayView2d< real64 const > const & phaseVolFrac =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionString );
      arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::dPhaseVolumeFraction_dPressureString );
      arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
        subRegion.getReference< array3d< real64 > >(
          viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );
      arrayView3d< real64 const > const & dCompFrac_dCompDens =
        subRegion.getReference< array3d< real64 > >(
          viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

      arrayView1d< real64 const > const & porosityOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );
      arrayView2d< real64 const > const & phaseVolFracOld =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionOldString );
      arrayView2d< real64 const > const & phaseDensOld =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseDensityOldString );
      arrayView3d< real64 const > const & phaseCompFracOld =
        subRegion.getReference< array3d< real64 > >(
          viewKeyStruct::phaseComponentFractionOldString );

      ConstitutiveBase const & solid =
        GetConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
      arrayView2d< real64 const > const & pvMult =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );
      arrayView2d< real64 const > const & dPvMult_dPres =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::dPVMult_dPresString );

      MultiFluidBase const & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                fluidModelNames()[targetIndex] );
      arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
      arrayView3d< real64 const > const & dPhaseDens_dPres =
        fluid.dPhaseDensity_dPressure();
      arrayView4d< real64 const > const & dPhaseDens_dComp =
        fluid.dPhaseDensity_dGlobalCompFraction();
      arrayView4d< real64 const > const & phaseCompFrac = fluid.phaseCompFraction();
      arrayView4d< real64 const > const & dPhaseCompFrac_dPres =
        fluid.dPhaseCompFraction_dPressure();
      arrayView5d< real64 const > const & dPhaseCompFrac_dComp =
        fluid.dPhaseCompFraction_dGlobalCompFraction();

      KernelLaunchSelector1< AccumulationKernel >( m_numComponents,
                                                   m_numPhases,
                                                   subRegion.size(),
                                                   dofManager.rankOffset(),
                                                   dofNumber,
                                                   elemGhostRank,
                                                   volume,
                                                   porosityOld,
                                                   porosityRef,
                                                   pvMult,
                                                   dPvMult_dPres,
                                                   dCompFrac_dCompDens,
                                                   phaseVolFracOld,
                                                   phaseVolFrac,
                                                   dPhaseVolFrac_dPres,
                                                   dPhaseVolFrac_dCompDens,
                                                   phaseDensOld,
                                                   phaseDens,
                                                   dPhaseDens_dPres,
                                                   dPhaseDens_dComp,
                                                   phaseCompFracOld,
                                                   phaseCompFrac,
                                                   dPhaseCompFrac_dPres,
                                                   dPhaseCompFrac_dComp,
                                                   localMatrix,
                                                   localRhs );
    } );
}

void
CompositionalMultiphaseFlow::AssembleFluxTerms(
  real64 const dt,
  DomainPartition const & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  /*
   * Force phase compositions to be moved to device.
   *
   * An issue with ElementViewAccessors is that if the outer arrays are already on device,
   * but an inner array gets touched and updated on host, capturing outer arrays in a device kernel
   * DOES NOT call move() on the inner array (see implementation of NewChaiBuffer::moveNested()).
   * Here we force the move by launching a dummy kernel.
   *
   * This is not a problem in normal solver execution, as these arrays get moved by AccumulationKernel.
   * But it fails unit tests, which test flux assembly separately.
   *
   * TODO: See if this can be fixed in NewChaiBuffer (I have not found a way - Sergey).
   *       Alternatively, stop using ElementViewAccessors altogether and just roll with
   *       accessors' outer arrays being moved on every jacobian assembly (maybe disable output).
   *       Or stop testing through the solver interface and test separate kernels instead.
   *       Finally, the problem should go away when fluid updates are executed on device.
   */
  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion ) {
      MultiFluidBase const & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                fluidModelNames()[targetIndex] );
      arrayView4d< real64 const > const & phaseCompFrac = fluid.phaseCompFraction();
      arrayView4d< real64 const > const & dPhaseCompFrac_dPres =
        fluid.dPhaseCompFraction_dPressure();
      arrayView5d< real64 const > const & dPhaseCompFrac_dComp =
        fluid.dPhaseCompFraction_dGlobalCompFraction();

      forAll< parallelDevicePolicy<> >(
        subRegion.size(),
        [phaseCompFrac,
         dPhaseCompFrac_dPres,
         dPhaseCompFrac_dComp] GEOSX_HOST_DEVICE( localIndex const ) {
          GEOSX_UNUSED_VAR( phaseCompFrac )
          GEOSX_UNUSED_VAR( dPhaseCompFrac_dPres )
          GEOSX_UNUSED_VAR( dPhaseCompFrac_dComp )
        } );
    } );

  NumericalMethodsManager const & numericalMethodManager =
    domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager =
    numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox =
    fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > > elemDofNumber =
    mesh.getElemManager()->ConstructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( getName() + "/accessors/" + dofKey );

  fluxApprox.forAllStencils( mesh, [&]( auto const & stencil ) {
    KernelLaunchSelector1< FluxKernel >(
      m_numComponents,
      m_numPhases,
      stencil,
      dofManager.rankOffset(),
      elemDofNumber.toViewConst(),
      m_elemGhostRank.toViewConst(),
      m_pressure.toViewConst(),
      m_deltaPressure.toViewConst(),
      m_gravCoef.toViewConst(),
      m_phaseMob.toViewConst(),
      m_dPhaseMob_dPres.toViewConst(),
      m_dPhaseMob_dCompDens.toViewConst(),
      m_dPhaseVolFrac_dPres.toViewConst(),
      m_dPhaseVolFrac_dCompDens.toViewConst(),
      m_dCompFrac_dCompDens.toViewConst(),
      m_phaseDens.toViewConst(),
      m_dPhaseDens_dPres.toViewConst(),
      m_dPhaseDens_dComp.toViewConst(),
      m_phaseCompFrac.toViewConst(),
      m_dPhaseCompFrac_dPres.toViewConst(),
      m_dPhaseCompFrac_dComp.toViewConst(),
      m_phaseCapPressure.toViewConst(),
      m_dPhaseCapPressure_dPhaseVolFrac.toViewConst(),
      m_capPressureFlag,
      dt,
      localMatrix.toViewConstSizes(),
      localRhs.toView() );
  } );
}

void
CompositionalMultiphaseFlow::AssembleVolumeBalanceTerms(
  DomainPartition const & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & porosityRef =
        subRegion.getReference< array1d< real64 > >(
          FlowSolverBase::viewKeyStruct::referencePorosityString );
      arrayView2d< real64 const > const & phaseVolFrac =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::phaseVolumeFractionString );
      arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::dPhaseVolumeFraction_dPressureString );
      arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
        subRegion.getReference< array3d< real64 > >(
          viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

      ConstitutiveBase const & solid =
        GetConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );
      arrayView2d< real64 const > const & pvMult =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );
      arrayView2d< real64 const > const & dPvMult_dPres =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::dPVMult_dPresString );

      KernelLaunchSelector2< VolumeBalanceKernel >( m_numComponents,
                                                    m_numPhases,
                                                    subRegion.size(),
                                                    dofManager.rankOffset(),
                                                    dofNumber,
                                                    elemGhostRank,
                                                    volume,
                                                    porosityRef,
                                                    pvMult,
                                                    dPvMult_dPres,
                                                    phaseVolFrac,
                                                    dPhaseVolFrac_dPres,
                                                    dPhaseVolFrac_dCompDens,
                                                    localMatrix.toViewConstSizes(),
                                                    localRhs.toView() );
    } );
}

void
CompositionalMultiphaseFlow::ApplyBoundaryConditions(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // apply pressure boundary conditions.
  ApplyDirichletBC( time_n,
                    dt,
                    dofManager,
                    domain,
                    localMatrix.toViewConstSizes(),
                    localRhs.toView() );

  // apply flux boundary conditions
  ApplySourceFluxBC( time_n,
                     dt,
                     dofManager,
                     domain,
                     localMatrix.toViewConstSizes(),
                     localRhs.toView() );
}

void
CompositionalMultiphaseFlow::ApplySourceFluxBC(
  real64 const time,
  real64 const dt,
  DofManager const & dofManager,
  DomainPartition & domain,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  fsManager.Apply(
    time + dt,
    &domain,
    "ElementRegions",
    FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
    [&]( FieldSpecificationBase const * const fs,
         string const &,
         SortedArrayView< localIndex const > const & lset,
         Group * const subRegion,
         string const & ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion->getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & ghostRank =
        subRegion->getReference< array1d< integer > >(
          ObjectManagerBase::viewKeyStruct::ghostRankString );

      SortedArray< localIndex > localSet;
      for( localIndex const a : lset )
      {
        if( ghostRank[a] < 0 )
        {
          localSet.insert( a );
        }
      }

      fs->ApplyBoundaryConditionToSystem< FieldSpecificationAdd,
                                          parallelDevicePolicy<> >(
        localSet.toViewConst(),
        time + dt,
        dt,
        subRegion,
        dofNumber,
        dofManager.rankOffset(),
        localMatrix,
        localRhs,
        [] GEOSX_HOST_DEVICE( localIndex const ) { return 0.0; } );
    } );
}

void
CompositionalMultiphaseFlow::ApplyDirichletBC(
  real64 const time,
  real64 const dt,
  DofManager const & dofManager,
  DomainPartition & domain,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  localIndex const NC = m_numComponents;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();

  map< string, map< string, array1d< bool > > > bcStatusMap;  // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  fsManager.Apply(
    time + dt,
    &domain,
    "ElementRegions",
    viewKeyStruct::pressureString,
    [&]( FieldSpecificationBase const * const fs,
         string const & setName,
         SortedArrayView< localIndex const > const & targetSet,
         Group * const subRegion,
         string const & ) {
      // 1.0. Check whether pressure has already been applied to this set
      string const & subRegionName = subRegion->getName();
      GEOSX_ERROR_IF(
        bcStatusMap[subRegionName].count( setName ) > 0,
        "Conflicting pressure boundary conditions on set " << setName );

      bcStatusMap[subRegionName][setName].resize( m_numComponents );
      bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

      // 1.1. Apply BC to set the field values
      fs->ApplyFieldValue< FieldSpecificationEqual, parallelHostPolicy >(
        targetSet,
        time + dt,
        subRegion,
        viewKeyStruct::bcPressureString );
    } );

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager.Apply(
    time + dt,
    &domain,
    "ElementRegions",
    viewKeyStruct::globalCompFractionString,
    [&]( FieldSpecificationBase const * const fs,
         string const & setName,
         SortedArrayView< localIndex const > const & targetSet,
         Group * const subRegion,
         string const & ) {
      // 2.0. Check pressure and record composition bc application
      string const & subRegionName = subRegion->getName();
      localIndex const comp = fs->GetComponent();
      GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0,
                      "Pressure boundary condition not prescribed on set '"
                        << setName << "'" );
      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp],
                      "Conflicting composition["
                        << comp << "] boundary conditions on set '" << setName
                        << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      // 2.1. Apply BC to set the field values
      fs->ApplyFieldValue< FieldSpecificationEqual, parallelHostPolicy >(
        targetSet,
        time + dt,
        subRegion,
        viewKeyStruct::globalCompFractionString );
    } );

  // 2.3 Check consistency between composition BC applied to sets
  bool bcConsistent = true;
  for( auto const & bcStatusEntryOuter : bcStatusMap )
  {
    for( auto const & bcStatusEntryInner : bcStatusEntryOuter.second )
    {
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        bcConsistent &= bcStatusEntryInner.second[ic];
        GEOSX_WARNING_IF(
          !bcConsistent,
          "Composition boundary condition not applied to component "
            << ic << " on region " << bcStatusEntryOuter.first << ", set "
            << bcStatusEntryInner.first );
      }
    }
  }
  GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager.Apply(
    time + dt,
    &domain,
    "ElementRegions",
    viewKeyStruct::pressureString,
    [&]( FieldSpecificationBase const * const,
         string const &,
         SortedArrayView< localIndex const > const & targetSet,
         Group * const subRegion,
         string const & ) {
      // TODO: hack! Find a better way to get the fluid
      Group const * const region = subRegion->getParent()->getParent();
      string const & fluidName =
        m_fluidModelNames[targetRegionIndex( region->getName() )];
      MultiFluidBase & fluid =
        GetConstitutiveModel< MultiFluidBase >( *subRegion, fluidName );

      arrayView1d< integer const > const ghostRank =
        subRegion->getReference< array1d< integer > >(
          ObjectManagerBase::viewKeyStruct::ghostRankString );
      arrayView1d< globalIndex const > const dofNumber =
        subRegion->getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const pres =
        subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
      arrayView1d< real64 const > const dPres =
        subRegion->getReference< array1d< real64 > >(
          viewKeyStruct::deltaPressureString );
      arrayView1d< real64 const > const bcPres =
        subRegion->getReference< array1d< real64 > >( viewKeyStruct::bcPressureString );
      arrayView2d< real64 const > const compFrac =
        subRegion->getReference< array2d< real64 > >(
          viewKeyStruct::globalCompFractionString );
      arrayView2d< real64 const > const compDens =
        subRegion->getReference< array2d< real64 > >(
          viewKeyStruct::globalCompDensityString );
      arrayView2d< real64 const > const dCompDens =
        subRegion->getReference< array2d< real64 > >(
          viewKeyStruct::deltaGlobalCompDensityString );
      arrayView2d< real64 const > const totalDens = fluid.totalDensity();

      constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid ) {
        typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
          castedFluid.createKernelWrapper();

        // MultiFluid models are not thread-safe or device-capable yet
        FluidUpdateKernel::Launch< serialPolicy >( targetSet,
                                                   fluidWrapper,
                                                   bcPres,
                                                   m_temperature,
                                                   compFrac );
      } );

      forAll< parallelDevicePolicy<> >(
        targetSet.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const a ) {
          localIndex const ei = targetSet[a];
          if( ghostRank[ei] >= 0 ) return;

          globalIndex const dofIndex = dofNumber[ei];
          localIndex const localRow = dofIndex - rankOffset;
          real64 rhsValue;

          // 3.1. Apply pressure value to the matrix/rhs
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcPres[ei],
                                                      pres[ei] + dPres[ei] );
          localRhs[localRow] = rhsValue;

          // 3.2. For each component, apply target global density value
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            FieldSpecificationEqual::SpecifyFieldValue(
              dofIndex + ic + 1,
              rankOffset,
              localMatrix,
              rhsValue,
              totalDens[ei][0] * compFrac[ei][ic],
              compDens[ei][ic] + dCompDens[ei][ic] );
            localRhs[localRow + ic + 1] = rhsValue;
          }
        } );
    } );
}

real64
CompositionalMultiphaseFlow::CalculateResidualNorm(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localRhs )
{
  localIndex const NDOF = m_numComponents + 1;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  real64 localResidualNorm = 0.0;

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion ) {
      MultiFluidBase const & fluid =
        GetConstitutiveModel< MultiFluidBase >( subRegion,
                                                m_fluidModelNames[targetIndex] );

      arrayView1d< globalIndex const > dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & refPoro =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::referencePorosityString );
      arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

      forAll< parallelDevicePolicy<> >(
        subRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          if( elemGhostRank[ei] < 0 )
          {
            localIndex const localRow = dofNumber[ei] - rankOffset;
            real64 const normalizer = totalDens[ei][0] * refPoro[ei] * volume[ei];

            for( localIndex idof = 0; idof < NDOF; ++idof )
            {
              real64 const val = localRhs[localRow + idof] / normalizer;
              localSum += val * val;
            }
          }
        } );

      localResidualNorm += localSum.get();
    } );

  // compute global residual norm
  real64 const residual = std::sqrt( MpiWrapper::Sum( localResidualNorm ) );

  if( getLogLevel() >= 1 && logger::internal::rank == 0 )
  {
    char output[200] = { 0 };
    sprintf( output, "    ( Rfluid ) = (%4.2e) ; ", residual );
    std::cout << output;
  }

  return residual;
}

void
CompositionalMultiphaseFlow::SolveSystem( DofManager const & dofManager,
                                          ParallelMatrix & matrix,
                                          ParallelVector & rhs,
                                          ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

real64
CompositionalMultiphaseFlow::ScalingForSystemSolution(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  real64 constexpr eps = minDensForDivision;
  real64 const maxCompFracChange = m_maxCompFracChange;

  localIndex const NC = m_numComponents;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );
  real64 scalingFactor = 1.0;

  forTargetSubRegions(
    mesh,
    [&]( localIndex const, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView2d< real64 const > const & compDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::globalCompDensityString );
      arrayView2d< real64 const > const & dCompDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::deltaGlobalCompDensityString );

      RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

      forAll< parallelDevicePolicy<> >(
        dofNumber.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          if( elemGhostRank[ei] < 0 )
          {
            real64 prevTotalDens = 0;
            for( localIndex ic = 0; ic < NC; ++ic )
            {
              prevTotalDens += compDens[ei][ic] + dCompDens[ei][ic];
            }

            // compute the change in component densities and component fractions
            for( localIndex ic = 0; ic < NC; ++ic )
            {
              localIndex const lid = dofNumber[ei] + ic + 1 - rankOffset;

              // compute scaling factor based on relative change in component densities
              real64 const absCompDensChange = fabs( localSolution[lid] );
              real64 const maxAbsCompDensChange =
                maxCompFracChange * prevTotalDens;

              // This actually checks the change in component fraction, using a lagged total density
              // Indeed we can rewrite the following check as:
              //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
              // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
              // because I found it more robust than using directly newTotalDens (which can vary also
              // wildly when the compDens change is large)
              if( absCompDensChange > maxAbsCompDensChange &&
                  absCompDensChange > eps )
              {
                minVal.min( maxAbsCompDensChange / absCompDensChange );
              }
            }
          }
        } );

      if( minVal.get() < scalingFactor )
      {
        scalingFactor = minVal.get();
      }
    } );

  return LvArray::math::max( MpiWrapper::Min( scalingFactor, MPI_COMM_GEOSX ),
                             m_minScalingFactor );
}

bool
CompositionalMultiphaseFlow::CheckSystemSolution(
  DomainPartition const & domain,
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor )
{
  real64 constexpr eps = minDensForDivision;

  localIndex const NC = m_numComponents;
  integer const allowCompDensChopping = m_allowCompDensChopping;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString );
  int localCheck = 1;

  forTargetSubRegions(
    mesh,
    [&]( localIndex const, ElementSubRegionBase const & subRegion ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

      arrayView1d< real64 const > const & pres =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
      arrayView1d< real64 const > const & dPres =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::deltaPressureString );
      arrayView2d< real64 const > const & compDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::globalCompDensityString );
      arrayView2d< real64 const > const & dCompDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::deltaGlobalCompDensityString );

      RAJA::ReduceMin< parallelDeviceReduce, integer > check( 1 );

      forAll< parallelDevicePolicy<> >(
        subRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          if( elemGhostRank[ei] < 0 )
          {
            localIndex const localRow = dofNumber[ei] - rankOffset;
            {
              real64 const newPres =
                pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
              check.min( newPres >= 0.0 );
            }

            // if component density chopping is not allowed, the time step fails if a component density is negative
            // otherwise, we just check that the total density is positive, and negative component densities
            // will be chopped (i.e., set to zero) in ApplySystemSolution)
            if( !allowCompDensChopping )
            {
              for( localIndex ic = 0; ic < NC; ++ic )
              {
                real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] +
                  scalingFactor * localSolution[localRow + ic + 1];
                check.min( newDens >= 0.0 );
              }
            }
            else
            {
              real64 totalDens = 0.0;
              for( localIndex ic = 0; ic < NC; ++ic )
              {
                real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] +
                  scalingFactor * localSolution[localRow + ic + 1];
                totalDens += ( newDens > 0.0 ) ? newDens : 0.0;
              }
              check.min( totalDens >= eps );
            }
          }
        } );

      localCheck = std::min( localCheck, check.get() );
    } );

  return MpiWrapper::Min( localCheck, MPI_COMM_GEOSX );
}

void
CompositionalMultiphaseFlow::ApplySystemSolution(
  DofManager const & dofManager,
  arrayView1d< real64 const > const & localSolution,
  real64 const scalingFactor,
  DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dofFieldString,
                               viewKeyStruct::deltaPressureString,
                               scalingFactor,
                               0,
                               1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dofFieldString,
                               viewKeyStruct::deltaGlobalCompDensityString,
                               scalingFactor,
                               1,
                               m_numDofPerCell );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    ChopNegativeDensities( domain );
  }

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ) );
  fieldNames["elems"].emplace_back(
    string( viewKeyStruct::deltaGlobalCompDensityString ) );
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
CompositionalMultiphaseFlow::ChopNegativeDensities( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;
  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    arrayView1d< integer const > const & ghostRank =
      subRegion.getReference< array1d< integer > >(
        ObjectManagerBase::viewKeyStruct::ghostRankString );

    arrayView2d< real64 const > const & compDens =
      subRegion.getReference< array2d< real64 > >(
        viewKeyStruct::globalCompDensityString );
    arrayView2d< real64 > const & dCompDens =
      subRegion.getReference< array2d< real64 > >(
        viewKeyStruct::deltaGlobalCompDensityString );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                        if( ghostRank[ei] < 0 )
                                        {
                                          for( localIndex ic = 0; ic < NC; ++ic )
                                          {
                                            real64 const newDens =
                                              compDens[ei][ic] + dCompDens[ei][ic];
                                            if( newDens < 0 )
                                            {
                                              dCompDens[ei][ic] = -compDens[ei][ic];
                                            }
                                          }
                                        }
                                      } );
  } );
}

void
CompositionalMultiphaseFlow::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      arrayView1d< real64 > const & dPres = subRegion.getReference< array1d< real64 > >(
        viewKeyStruct::deltaPressureString );
      arrayView2d< real64 > const & dCompDens =
        subRegion.getReference< array2d< real64 > >(
          viewKeyStruct::deltaGlobalCompDensityString );

      dPres.setValues< parallelDevicePolicy<> >( 0.0 );
      dCompDens.setValues< parallelDevicePolicy<> >( 0.0 );

      UpdateState( subRegion, targetIndex );
    } );
}

void
CompositionalMultiphaseFlow::ImplicitStepComplete(
  real64 const & GEOSX_UNUSED_PARAM( time ),
  real64 const & GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & domain )
{
  localIndex const NC = m_numComponents;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    arrayView1d< real64 const > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView2d< real64 const > const & dCompDens =
      subRegion.getReference< array2d< real64 > >(
        viewKeyStruct::deltaGlobalCompDensityString );

    arrayView1d< real64 > const & pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView2d< real64 > const & compDens = subRegion.getReference< array2d< real64 > >(
      viewKeyStruct::globalCompDensityString );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                        pres[ei] += dPres[ei];
                                        for( localIndex ic = 0; ic < NC; ++ic )
                                        {
                                          compDens[ei][ic] += dCompDens[ei][ic];
                                        }
                                      } );
  } );
}

void
CompositionalMultiphaseFlow::ResetViews( MeshLevel & mesh )
{
  FlowSolverBase::ResetViews( mesh );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::pressureString );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::deltaPressureString );
  m_deltaPressure.setName( getName() + "/accessors/" +
                           viewKeyStruct::deltaPressureString );

  m_dCompFrac_dCompDens.clear();
  m_dCompFrac_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >(
    viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );
  m_dCompFrac_dCompDens.setName(
    getName() + "/accessors/" +
    viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString );

  m_dPhaseVolFrac_dPres.clear();
  m_dPhaseVolFrac_dPres = elemManager.ConstructArrayViewAccessor< real64, 2 >(
    viewKeyStruct::dPhaseVolumeFraction_dPressureString );
  m_dPhaseVolFrac_dPres.setName(
    getName() + "/accessors/" +
    viewKeyStruct::dPhaseVolumeFraction_dPressureString );

  m_dPhaseVolFrac_dCompDens.clear();
  m_dPhaseVolFrac_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >(
    viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );
  m_dPhaseVolFrac_dCompDens.setName(
    getName() + "/accessors/" +
    viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString );

  m_phaseMob.clear();
  m_phaseMob = elemManager.ConstructArrayViewAccessor< real64, 2 >(
    viewKeyStruct::phaseMobilityString );
  m_phaseMob.setName( getName() + "/accessors/" +
                      viewKeyStruct::phaseMobilityString );

  m_dPhaseMob_dPres.clear();
  m_dPhaseMob_dPres = elemManager.ConstructArrayViewAccessor< real64, 2 >(
    viewKeyStruct::dPhaseMobility_dPressureString );
  m_dPhaseMob_dPres.setName( getName() + "/accessors/" +
                             viewKeyStruct::dPhaseMobility_dPressureString );

  m_dPhaseMob_dCompDens.clear();
  m_dPhaseMob_dCompDens = elemManager.ConstructArrayViewAccessor< real64, 3 >(
    viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );
  m_dPhaseMob_dCompDens.setName(
    getName() + "/accessors/" +
    viewKeyStruct::dPhaseMobility_dGlobalCompDensityString );

  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_phaseDens.clear();
    m_phaseDens = elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >(
      keys::phaseDensityString,
      targetRegionNames(),
      fluidModelNames() );
    m_phaseDens.setName( getName() + "/accessors/" + keys::phaseDensityString );

    m_dPhaseDens_dPres.clear();
    m_dPhaseDens_dPres =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >(
        keys::dPhaseDensity_dPressureString,
        targetRegionNames(),
        fluidModelNames() );
    m_dPhaseDens_dPres.setName( getName() + "/accessors/" +
                                keys::dPhaseDensity_dPressureString );

    m_dPhaseDens_dComp.clear();
    m_dPhaseDens_dComp =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >(
        keys::dPhaseDensity_dGlobalCompFractionString,
        targetRegionNames(),
        fluidModelNames() );
    m_dPhaseDens_dComp.setName( getName() + "/accessors/" +
                                keys::dPhaseDensity_dGlobalCompFractionString );

    m_phaseCompFrac.clear();
    m_phaseCompFrac = elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >(
      keys::phaseCompFractionString,
      targetRegionNames(),
      fluidModelNames() );
    m_phaseCompFrac.setName( getName() + "/accessors/" +
                             keys::phaseCompFractionString );

    m_dPhaseCompFrac_dPres.clear();
    m_dPhaseCompFrac_dPres =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >(
        keys::dPhaseCompFraction_dPressureString,
        targetRegionNames(),
        fluidModelNames() );
    m_dPhaseCompFrac_dPres.setName( getName() + "/accessors/" +
                                    keys::dPhaseCompFraction_dPressureString );

    m_dPhaseCompFrac_dComp.clear();
    m_dPhaseCompFrac_dComp =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 5 >(
        keys::dPhaseCompFraction_dGlobalCompFractionString,
        targetRegionNames(),
        fluidModelNames() );
    m_dPhaseCompFrac_dComp.setName(
      getName() + "/accessors/" +
      keys::dPhaseCompFraction_dGlobalCompFractionString );
  }
  if( m_capPressureFlag )
  {
    using keys = CapillaryPressureBase::viewKeyStruct;

    m_phaseCapPressure.clear();
    m_phaseCapPressure =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 3 >(
        keys::phaseCapPressureString,
        targetRegionNames(),
        capPresModelNames() );
    m_phaseCapPressure.setName( getName() + "/accessors/" +
                                keys::phaseCapPressureString );

    m_dPhaseCapPressure_dPhaseVolFrac.clear();
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager.ConstructMaterialArrayViewAccessor< real64, 4 >(
        keys::dPhaseCapPressure_dPhaseVolFractionString,
        targetRegionNames(),
        capPresModelNames() );
    m_dPhaseCapPressure_dPhaseVolFrac.setName(
      getName() + "/accessors/" + keys::dPhaseCapPressure_dPhaseVolFractionString );
  }
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase,
                        CompositionalMultiphaseFlow,
                        string const &,
                        Group * const )

}  // namespace geosx
