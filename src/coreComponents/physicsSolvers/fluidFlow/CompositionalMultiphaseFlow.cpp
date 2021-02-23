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
#include "managers/GeosxState.hpp"
#include "managers/ProblemManager.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseFlowKernels.hpp"
#include "physicsSolvers/fluidFlow/IsothermalCompositionalMultiphaseFlowKernels.hpp"
#include "constitutive/solid/SolidBase.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

static constexpr real64 minDensForDivision = 1e-10;

CompositionalMultiphaseFlow::CompositionalMultiphaseFlow( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_capPressureFlag( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{
//START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::uniformTemperatureString(), &m_uniformTemperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Uniform temperature for isothermal simulations." );

  this->registerWrapper( viewKeyStruct::isothermalFlagString(), &m_isothermalFlag ).
    setApplyDefaultValue( 1 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether the problem  is isothermal or not." );

  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNamesString(), &m_relPermModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Name of the relative permeability constitutive model to use" );

  this->registerWrapper( viewKeyStruct::capPressureNamesString(), &m_capPressureModelNames ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Name of the capillary pressure constitutive model to use" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

  m_linearSolverParameters.get().mgr.strategy = "CompositionalMultiphaseFlow";

}

void CompositionalMultiphaseFlow::postProcessInput()
{
  FlowSolverBase::postProcessInput();
  checkModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString() );
  m_capPressureFlag = checkModelNames( m_capPressureModelNames, viewKeyStruct::capPressureNamesString(), true );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );
}

void CompositionalMultiphaseFlow::registerDataOnMesh( Group & meshBodies )
{
  FlowSolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & meshLevel = meshBody.getMeshLevel( 0 );

    forTargetSubRegions( meshLevel, [&]( localIndex const, ElementSubRegionBase & elementSubRegion )
    {

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::bcPressureString() );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompDensityString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompFractionString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseMobilityString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      if( m_isothermalFlag == 0 )
      {
        elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::temperatureString() ).
          setPlotLevel( PlotLevel::LEVEL_0 );

        elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaTemperatureString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dTemperatureString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseInternalEnergyOldString() ).
          setRestartFlags( RestartFlags::NO_WRITE );

        elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::rockInternalEnergyOldString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
      }
    } );
  } );
}

namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void CompareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG( lhs.numFluidPhases(), rhs.numFluidPhases(),
                         "Mismatch in number of phases between constitutive models "
                         << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ip = 0; ip < lhs.numFluidPhases(); ++ip )
  {
    GEOSX_ERROR_IF_NE_MSG( lhs.phaseNames()[ip], rhs.phaseNames()[ip],
                           "Mismatch in phase names between constitutive models "
                           << lhs.getName() << " and " << rhs.getName() );
  }
}

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void CompareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
{
  GEOSX_ERROR_IF_NE_MSG( lhs.numFluidComponents(), rhs.numFluidComponents(),
                         "Mismatch in number of components between constitutive models "
                         << lhs.getName() << " and " << rhs.getName() );

  for( localIndex ic = 0; ic < lhs.numFluidComponents(); ++ic )
  {
    GEOSX_ERROR_IF_NE_MSG( lhs.componentNames()[ic], rhs.componentNames()[ic],
                           "Mismatch in component names between constitutive models "
                           << lhs.getName() << " and " << rhs.getName() );
  }
}

}

void CompositionalMultiphaseFlow::validateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const
{
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  for( localIndex i = 1; i < m_fluidModelNames.size(); ++i )
  {
    MultiFluidBase const & fluid = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[i] );
    CompareMultiphaseModels( fluid, fluid0 );
    CompareMulticomponentModels( fluid, fluid0 );
  }

  RelativePermeabilityBase const & relPerm0 = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[0] );
  CompareMultiphaseModels( relPerm0, fluid0 );

  for( localIndex i = 1; i < m_relPermModelNames.size(); ++i )
  {
    RelativePermeabilityBase const & relPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[i] );
    CompareMultiphaseModels( relPerm, relPerm0 );
  }

  if( m_capPressureFlag )
  {
    CapillaryPressureBase const & capPres0 = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[0] );
    CompareMultiphaseModels( capPres0, fluid0 );

    for( localIndex i = 1; i < m_capPressureModelNames.size(); ++i )
    {
      CapillaryPressureBase const & capPres = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[i] );
      CompareMultiphaseModels( capPres, capPres0 );
    }
  }
}


void CompositionalMultiphaseFlow::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Set key dimensions of the problem
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases     = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  if( m_isothermalFlag )
  {
    // Nc + pressure
    m_numDofPerCell = m_numComponents + 1;
  }
  else
  {
    // Nc + pressure + tempertature
    m_numDofPerCell = m_numComponents + 2;
  }

  // 2. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( cm );

  // 3. Check that the discretization is valid for this solver
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with CompositionalMultiphaseFlow" );
  }

  // 3. Resize all fields as necessary, validate constitutive models in regions
  for( auto & mesh : domain.getMeshBodies().getSubGroups() )
  {
    MeshLevel & meshLevel = dynamicCast< MeshBody * >( mesh.second )->getMeshLevel( 0 );
    resizeFields( meshLevel );

    validateModelMapping< MultiFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
    validateModelMapping< RelativePermeabilityBase >( meshLevel.getElemManager(), m_relPermModelNames );
    if( m_capPressureFlag )
    {
      validateModelMapping< CapillaryPressureBase >( meshLevel.getElemManager(), m_capPressureModelNames );
    }
  }
}

void CompositionalMultiphaseFlow::resizeFields( MeshLevel & meshLevel ) const
{
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  forTargetSubRegions( meshLevel, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() ).resizeDimension< 1 >( NC );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() ).resizeDimension< 1 >( NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() ).resizeDimension< 1 >( NC );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).resizeDimension< 1, 2 >( NC, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).resizeDimension< 1, 2 >( NP, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() ).resizeDimension< 1, 2 >( NP, NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseDensityOldString() ).resizeDimension< 1 >( NP );
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString() ).resizeDimension< 1, 2 >( NP, NC );
  } );

  if( m_isothermalFlag == 0 )
  {
    forTargetSubRegions( meshLevel, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() ).resizeDimension< 1 >( NP );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dTemperatureString() ).resizeDimension< 1 >( NP );
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseInternalEnergyOldString() ).resizeDimension< 1 >( NP );
    } );
  }
}

void CompositionalMultiphaseFlow::updateComponentFraction( Group & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & compFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );

  arrayView3d< real64 > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  // inputs

  arrayView2d< real64 const > const compDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );

  arrayView2d< real64 const > const dCompDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  CompositionalMultiphaseFlowKernels::KernelLaunchSelector1
  < CompositionalMultiphaseFlowKernels::ComponentFractionKernel >( m_numComponents,
                                                                   dataGroup.size(),
                                                                   compDens,
                                                                   dCompDens,
                                                                   compFrac,
                                                                   dCompFrac_dCompDens );
}

void CompositionalMultiphaseFlow::updatePhaseVolumeFraction( Group & dataGroup,
                                                             localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );

  arrayView2d< real64 > const dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );

  arrayView3d< real64 > const dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  // inputs

  arrayView3d< real64 const > const dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  arrayView2d< real64 const > const compDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );

  arrayView2d< real64 const > const dCompDens =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  if( m_isothermalFlag )
  {
    IsothermalCompositionalMultiphaseFlowKernels::KernelLaunchSelector2
    < IsothermalCompositionalMultiphaseFlowKernels::PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
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
  else
  {
    // Get thermal views
    arrayView3d< real64 const > const & dPhaseFrac_dTemp = fluid.dPhaseFraction_dTemperature();
    arrayView3d< real64 const > const & dPhaseDens_dTemp = fluid.dPhaseDensity_dTemperature();

    arrayView2d< real64 > const dPhaseVolFrac_dTemp =
      dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );

    CompositionalMultiphaseFlowKernels::KernelLaunchSelector2
    < CompositionalMultiphaseFlowKernels::PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
                                                                       dataGroup.size(),
                                                                       compDens,
                                                                       dCompDens,
                                                                       dCompFrac_dCompDens,
                                                                       phaseDens,
                                                                       dPhaseDens_dPres,
                                                                       dPhaseDens_dTemp,
                                                                       dPhaseDens_dComp,
                                                                       phaseFrac,
                                                                       dPhaseFrac_dPres,
                                                                       dPhaseFrac_dTemp,
                                                                       dPhaseFrac_dComp,
                                                                       phaseVolFrac,
                                                                       dPhaseVolFrac_dPres,
                                                                       dPhaseVolFrac_dTemp,
                                                                       dPhaseVolFrac_dComp );

  }
}

void CompositionalMultiphaseFlow::updatePhaseMobility( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const phaseMob =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseMobilityString() );

  arrayView2d< real64 > const dPhaseMob_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dPressureString() );

  arrayView3d< real64 > const dPhaseMob_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() );

  // inputs

  arrayView2d< real64 const > const dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );

  arrayView3d< real64 const > const dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  arrayView3d< real64 const > const dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseVisc = fluid.phaseViscosity();
  arrayView3d< real64 const > const & dPhaseVisc_dPres = fluid.dPhaseViscosity_dPressure();
  arrayView4d< real64 const > const & dPhaseVisc_dComp = fluid.dPhaseViscosity_dGlobalCompFraction();

  RelativePermeabilityBase const & relperm = getConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseRelPerm = relperm.phaseRelPerm();
  arrayView4d< real64 const > const & dPhaseRelPerm_dPhaseVolFrac = relperm.dPhaseRelPerm_dPhaseVolFraction();

  if( m_isothermalFlag )
  {
    IsothermalCompositionalMultiphaseFlowKernels::KernelLaunchSelector2
    < IsothermalCompositionalMultiphaseFlowKernels::PhaseMobilityKernel >( m_numComponents, m_numPhases,
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
  else
  {
    // get thermal views
    arrayView3d< real64 const > const & dPhaseDens_dTemp = fluid.dPhaseDensity_dTemperature();
    arrayView3d< real64 const > const & dPhaseVisc_dTemp = fluid.dPhaseViscosity_dTemperature();

    arrayView2d< real64 > const dPhaseMob_dTemp =
      dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseMobility_dTemperatureString() );

    arrayView2d< real64 const > const dPhaseVolFrac_dTemp =
      dataGroup.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );

    CompositionalMultiphaseFlowKernels::KernelLaunchSelector2
    < CompositionalMultiphaseFlowKernels::PhaseMobilityKernel >( m_numComponents, m_numPhases,
                                                                 dataGroup.size(),
                                                                 dCompFrac_dCompDens,
                                                                 phaseDens,
                                                                 dPhaseDens_dPres,
                                                                 dPhaseDens_dTemp,
                                                                 dPhaseDens_dComp,
                                                                 phaseVisc,
                                                                 dPhaseVisc_dPres,
                                                                 dPhaseVisc_dTemp,
                                                                 dPhaseVisc_dComp,
                                                                 phaseRelPerm,
                                                                 dPhaseRelPerm_dPhaseVolFrac,
                                                                 dPhaseVolFrac_dPres,
                                                                 dPhaseVolFrac_dTemp,
                                                                 dPhaseVolFrac_dComp,
                                                                 phaseMob,
                                                                 dPhaseMob_dPres,
                                                                 dPhaseMob_dTemp,
                                                                 dPhaseMob_dComp );
  }

}

void CompositionalMultiphaseFlow::updateFluidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
  arrayView2d< real64 const > const compFrac = dataGroup.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // MultiFluid models are not thread-safe or device-capable yet
    if( m_isothermalFlag )
    {
      IsothermalCompositionalMultiphaseFlowKernels::FluidUpdateKernel::launch< serialPolicy >( dataGroup.size(),
                                                                                               fluidWrapper,
                                                                                               pres,
                                                                                               dPres,
                                                                                               m_uniformTemperature,
                                                                                               compFrac );
    }
    else
    {
      arrayView1d< real64 const > const temperature  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::temperatureString() );
      arrayView1d< real64 const > const dTemperature = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaTemperatureString() );

      CompositionalMultiphaseFlowKernels::FluidUpdateKernel::launch< serialPolicy >( dataGroup.size(),
                                                                                     fluidWrapper,
                                                                                     pres,
                                                                                     dPres,
                                                                                     temperature,
                                                                                     dTemperature,
                                                                                     compFrac );
    }

  } );
}

void CompositionalMultiphaseFlow::updateSolidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase & solid = getConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidModelNames[targetIndex] );

  arrayView1d< real64 const > const pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  solid.stateUpdateBatchPressure( pres, dPres );
}

void CompositionalMultiphaseFlow::updateRelPermModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< real64 const > const phaseVolFrac =
    dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );

  RelativePermeabilityBase & relPerm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    CompositionalMultiphaseFlowKernels::RelativePermeabilityUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                                                            relPermWrapper,
                                                                                                            phaseVolFrac );
  } );
}

void CompositionalMultiphaseFlow::updateCapPressureModel( Group & dataGroup, localIndex const targetIndex ) const
{
  if( m_capPressureFlag )
  {
    arrayView2d< real64 const > const phaseVolFrac =
      dataGroup.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );

    CapillaryPressureBase & capPressure =
      getConstitutiveModel< CapillaryPressureBase >( dataGroup, m_capPressureModelNames[targetIndex] );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      CompositionalMultiphaseFlowKernels::CapillaryPressureUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                                                           capPresWrapper,
                                                                                                           phaseVolFrac );
    } );
  }
}

void CompositionalMultiphaseFlow::updateState( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  updateComponentFraction( dataGroup );
  updateFluidModel( dataGroup, targetIndex );
  updatePhaseVolumeFraction( dataGroup, targetIndex );
  updateSolidModel( dataGroup, targetIndex );
  updateRelPermModel( dataGroup, targetIndex );
  updatePhaseMobility( dataGroup, targetIndex );
  updateCapPressureModel( dataGroup, targetIndex );
}

void CompositionalMultiphaseFlow::initializeFluidState( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    // 1. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion, targetIndex );

    // 2. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const > const totalDens = fluid.totalDensity();

    arrayView2d< real64 const > const compFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
    arrayView2d< real64 > const
    compDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compDens[ei][ic] = totalDens[ei][0] * compFrac[ei][ic];
      }
    } );

    // 3. Update dependent state quantities
    updatePhaseVolumeFraction( subRegion, targetIndex );
    updateSolidModel( subRegion, targetIndex );
    updateRelPermModel( subRegion, targetIndex );
    updatePhaseMobility( subRegion, targetIndex );
    updateCapPressureModel( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseFlow::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = getGlobalState().getProblemManager().getDomainPartition();

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::pressureString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::globalCompDensityString() ) );

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames, mesh, domain.getNeighbors() );

  // set mass fraction flag on fluid models
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );

  // Initialize primary variables from applied initial conditions
  initializeFluidState( mesh );
}

real64 CompositionalMultiphaseFlow::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseFlow::backupFields( MeshLevel & mesh ) const
{
  // backup some fields used in time derivative approximation
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const poroRef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
    arrayView2d< real64 const > const phaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView3d< real64 const > const phaseDens = fluid.phaseDensity();
    arrayView4d< real64 const > const phaseCompFrac = fluid.phaseCompFraction();

    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
    arrayView2d< real64 const > const pvMult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );

    arrayView2d< real64 > const phaseDensOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseDensityOldString() );
    arrayView2d< real64 > const phaseVolFracOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString() );
    arrayView3d< real64 > const phaseCompFracOld =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString() );
    arrayView1d< real64 > const poroOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString() );

    localIndex const NC = m_numComponents;
    localIndex const NP = m_numPhases;

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
        return;

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

  //  backup thermal fields
  if( m_isothermalFlag == 0 )
  {
    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
      SolidBase const & solid = getConstitutiveModel< SolidBase >( subRegion, solidModelNames()[targetIndex] );

      arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();

      arrayView3d< real64 const > const phaseInternalEnergy = fluid.phaseInternalEnergy();
      arrayView2d< real64 const > const rockInternalEnergy = solid.internalEnergy();

      arrayView2d< real64 > const phaseInternalEnergyOld =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseInternalEnergyOldString() );
      arrayView1d< real64 > rockInternalEnergyOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::rockInternalEnergyOldString() );

      localIndex const NP = m_numPhases;

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] >= 0 )
          return;

        for( localIndex ip = 0; ip < NP; ++ip )
        {
          phaseInternalEnergyOld[ei][ip] = phaseInternalEnergy[ei][0][ip];
          rockInternalEnergyOld[ei] = rockInternalEnergy[ei][0];
        }

      } );
    } );
  }
}

void
CompositionalMultiphaseFlow::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // bind the stored views to the current domain
  resetViews( mesh );

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  backupFields( mesh );
}

void CompositionalMultiphaseFlow::setupDofs( DomainPartition const & domain,
                                             DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::dofFieldString(),
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );
  dofManager.addCoupling( viewKeyStruct::dofFieldString(), fluxApprox );
}

void CompositionalMultiphaseFlow::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                  real64 const dt,
                                                  DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  assembleAccumulationTerms( domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  assembleFluxTerms( dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );

  assembleVolumeBalanceTerms( domain,
                              dofManager,
                              localMatrix,
                              localRhs );
}

void CompositionalMultiphaseFlow::assembleAccumulationTerms( DomainPartition const & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & porosityRef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );

    arrayView2d< real64 const > const & phaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );
    arrayView3d< real64 const > const & dCompFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

    arrayView1d< real64 const > const & porosityOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString() );
    arrayView2d< real64 const > const & phaseVolFracOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionOldString() );
    arrayView2d< real64 const > const & phaseDensOld =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseDensityOldString() );
    arrayView3d< real64 const > const & phaseCompFracOld =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::phaseComponentFractionOldString() );

    SolidBase const & solid = getConstitutiveModel< SolidBase >( subRegion, solidModelNames()[targetIndex] );

    arrayView2d< real64 const > const & pvMult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );
    arrayView2d< real64 const > const & dPvMult_dPres =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString() );

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
    arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
    arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();
    arrayView4d< real64 const > const & phaseCompFrac = fluid.phaseCompFraction();
    arrayView4d< real64 const > const & dPhaseCompFrac_dPres = fluid.dPhaseCompFraction_dPressure();
    arrayView5d< real64 const > const & dPhaseCompFrac_dComp = fluid.dPhaseCompFraction_dGlobalCompFraction();

    if( m_isothermalFlag )
    {
      IsothermalCompositionalMultiphaseFlowKernels::KernelLaunchSelector1
      < IsothermalCompositionalMultiphaseFlowKernels::AccumulationKernel >( m_numComponents,
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
    }
    else
    {

      arrayView2d< real64 const > const & dPhaseVolFrac_dTemp =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );

      arrayView3d< real64 const > const & dPhaseDens_dTemp = fluid.dPhaseDensity_dTemperature();
      arrayView4d< real64 const > const & dPhaseCompFrac_dTemp = fluid.dPhaseCompFraction_dTemperature();

      arrayView2d< real64 const > const & phaseInternalEnergyOld =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseInternalEnergyOldString() );
      arrayView3d< real64 const > const & phaseInternalEnergy = fluid.phaseInternalEnergy();
      arrayView3d< real64 const > const & dPhaseInternalEnergy_dPres = fluid.dPhaseInternalEnergy_dPressure();
      arrayView3d< real64 const > const & dPhaseInternalEnergy_dTemp = fluid.dPhaseInternalEnergy_dTemperature();
      arrayView4d< real64 const > const & dPhaseInternalEnergy_dComp = fluid.dPhaseInternalEnergy_dGlobalCompFraction();

      arrayView1d< real64 const > const & rockInternalEnergyOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::rockInternalEnergyOldString() );
      arrayView2d< real64 const > const & rockInternalEnergy = solid.internalEnergy();
      arrayView2d< real64 const > const & dRockInternalEnergy_dTemp = solid.dInternalEnergy_dTemperature();
      arrayView2d< real64 const > const & rockDensity = solid.getDensity();

      CompositionalMultiphaseFlowKernels::KernelLaunchSelector1
      < CompositionalMultiphaseFlowKernels::AccumulationKernel >( m_numComponents,
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
                                                                  dPhaseVolFrac_dTemp,
                                                                  dPhaseVolFrac_dCompDens,
                                                                  phaseDensOld,
                                                                  phaseDens,
                                                                  dPhaseDens_dPres,
                                                                  dPhaseDens_dTemp,
                                                                  dPhaseDens_dComp,
                                                                  phaseCompFracOld,
                                                                  phaseCompFrac,
                                                                  dPhaseCompFrac_dPres,
                                                                  dPhaseCompFrac_dTemp,
                                                                  dPhaseCompFrac_dComp,
                                                                  phaseInternalEnergyOld,
                                                                  phaseInternalEnergy,
                                                                  dPhaseInternalEnergy_dPres,
                                                                  dPhaseInternalEnergy_dTemp,
                                                                  dPhaseInternalEnergy_dComp,
                                                                  rockInternalEnergyOld,
                                                                  rockInternalEnergy,
                                                                  dRockInternalEnergy_dTemp,
                                                                  rockDensity,
                                                                  localMatrix,
                                                                  localRhs );
    }

  } );
}

void CompositionalMultiphaseFlow::assembleFluxTerms( real64 const dt,
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

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
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView4d< real64 const > const & phaseCompFrac = fluid.phaseCompFraction();
    arrayView4d< real64 const > const & dPhaseCompFrac_dPres = fluid.dPhaseCompFraction_dPressure();
    arrayView5d< real64 const > const & dPhaseCompFrac_dComp = fluid.dPhaseCompFraction_dGlobalCompFraction();

    arrayView3d< real64 const > const & phaseMassDens = fluid.phaseMassDensity();
    arrayView3d< real64 const > const & dPhaseMassDens_dPres = fluid.dPhaseMassDensity_dPressure();
    arrayView4d< real64 const > const & dPhaseMassDens_dComp = fluid.dPhaseMassDensity_dGlobalCompFraction();

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [phaseCompFrac, dPhaseCompFrac_dPres, dPhaseCompFrac_dComp,
                                       phaseMassDens, dPhaseMassDens_dPres, dPhaseMassDens_dComp]
                                      GEOSX_HOST_DEVICE ( localIndex const )
    {
      GEOSX_UNUSED_VAR( phaseCompFrac )
      GEOSX_UNUSED_VAR( dPhaseCompFrac_dPres )
      GEOSX_UNUSED_VAR( dPhaseCompFrac_dComp )
      GEOSX_UNUSED_VAR( phaseMassDens )
      GEOSX_UNUSED_VAR( dPhaseMassDens_dPres )
      GEOSX_UNUSED_VAR( dPhaseMassDens_dComp )
    } );
  } );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( getName() + "/accessors/" + dofKey );

  fluxApprox.forAllStencils( mesh, [&] ( auto const & stencil )
  {
    if( m_isothermalFlag )
    {
      IsothermalCompositionalMultiphaseFlowKernels::KernelLaunchSelector1
      < IsothermalCompositionalMultiphaseFlowKernels::FluxKernel >( m_numComponents,
                                                                    m_numPhases,
                                                                    stencil,
                                                                    dofManager.rankOffset(),
                                                                    elemDofNumber.toNestedViewConst(),
                                                                    m_elemGhostRank.toNestedViewConst(),
                                                                    m_pressure.toNestedViewConst(),
                                                                    m_deltaPressure.toNestedViewConst(),
                                                                    m_gravCoef.toNestedViewConst(),
                                                                    m_phaseMob.toNestedViewConst(),
                                                                    m_dPhaseMob_dPres.toNestedViewConst(),
                                                                    m_dPhaseMob_dCompDens.toNestedViewConst(),
                                                                    m_dPhaseVolFrac_dPres.toNestedViewConst(),
                                                                    m_dPhaseVolFrac_dCompDens.toNestedViewConst(),
                                                                    m_dCompFrac_dCompDens.toNestedViewConst(),
                                                                    m_phaseMassDens.toNestedViewConst(),
                                                                    m_dPhaseMassDens_dPres.toNestedViewConst(),
                                                                    m_dPhaseMassDens_dComp.toNestedViewConst(),
                                                                    m_phaseCompFrac.toNestedViewConst(),
                                                                    m_dPhaseCompFrac_dPres.toNestedViewConst(),
                                                                    m_dPhaseCompFrac_dComp.toNestedViewConst(),
                                                                    m_phaseCapPressure.toNestedViewConst(),
                                                                    m_dPhaseCapPressure_dPhaseVolFrac.toNestedViewConst(),
                                                                    m_capPressureFlag,
                                                                    dt,
                                                                    localMatrix.toViewConstSizes(),
                                                                    localRhs.toView()  );
    }
    else
    {
      CompositionalMultiphaseFlowKernels::KernelLaunchSelector1
      < CompositionalMultiphaseFlowKernels::FluxKernel >( m_numComponents,
                                                          m_numPhases,
                                                          stencil,
                                                          dofManager.rankOffset(),
                                                          elemDofNumber.toNestedViewConst(),
                                                          m_elemGhostRank.toNestedViewConst(),
                                                          m_pressure.toNestedViewConst(),
                                                          m_deltaPressure.toNestedViewConst(),
                                                          m_gravCoef.toNestedViewConst(),
                                                          m_phaseMob.toNestedViewConst(),
                                                          m_dPhaseMob_dPres.toNestedViewConst(),
                                                          m_dPhaseMob_dTemp.toNestedViewConst(),
                                                          m_dPhaseMob_dCompDens.toNestedViewConst(),
                                                          m_phaseVolFrac.toNestedViewConst(),
                                                          m_dPhaseVolFrac_dPres.toNestedViewConst(),
                                                          m_dPhaseVolFrac_dTemp.toNestedViewConst(),
                                                          m_dPhaseVolFrac_dCompDens.toNestedViewConst(),
                                                          m_dCompFrac_dCompDens.toNestedViewConst(),
                                                          m_phaseMassDens.toNestedViewConst(),
                                                          m_dPhaseMassDens_dPres.toNestedViewConst(),
                                                          m_dPhaseMassDens_dTemp.toNestedViewConst(),
                                                          m_dPhaseMassDens_dComp.toNestedViewConst(),
                                                          m_phaseCompFrac.toNestedViewConst(),
                                                          m_dPhaseCompFrac_dPres.toNestedViewConst(),
                                                          m_dPhaseCompFrac_dTemp.toNestedViewConst(),
                                                          m_dPhaseCompFrac_dComp.toNestedViewConst(),
                                                          m_phaseCapPressure.toNestedViewConst(),
                                                          m_dPhaseCapPressure_dPhaseVolFrac.toNestedViewConst(),
                                                          m_temperature.toNestedViewConst(),
                                                          m_deltaTemperature.toNestedViewConst(),
                                                          m_phaseEnthalpy.toNestedViewConst(),
                                                          m_dPhaseEnthalpy_dPres.toNestedViewConst(),
                                                          m_dPhaseEnthalpy_dTemp.toNestedViewConst(),
                                                          m_dPhaseEnthalpy_dComp.toNestedViewConst(),
                                                          m_capPressureFlag,
                                                          dt,
                                                          localMatrix.toViewConstSizes(),
                                                          localRhs.toView()  );
    }
  } );
}

void CompositionalMultiphaseFlow::assembleVolumeBalanceTerms( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & porosityRef =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::referencePorosityString() );
    arrayView2d< real64 const > const & phaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
    arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );
    arrayView2d< real64 const > const & pvMult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );
    arrayView2d< real64 const > const & dPvMult_dPres =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString() );

    if( m_isothermalFlag )
    {
      IsothermalCompositionalMultiphaseFlowKernels::KernelLaunchSelector2
      < IsothermalCompositionalMultiphaseFlowKernels::VolumeBalanceKernel >( m_numComponents, m_numPhases,
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
    }
    else
    {
      arrayView2d< real64 const > const & dPhaseVolFrac_dTemp =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );

      CompositionalMultiphaseFlowKernels::KernelLaunchSelector2
      < CompositionalMultiphaseFlowKernels::VolumeBalanceKernel >( m_numComponents, m_numPhases,
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
                                                                   dPhaseVolFrac_dTemp,
                                                                   dPhaseVolFrac_dCompDens,
                                                                   localMatrix.toViewConstSizes(),
                                                                   localRhs.toView() );
    }
  } );
}

void CompositionalMultiphaseFlow::applyBoundaryConditions( real64 const time_n,
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // apply pressure boundary conditions.
  applyDirichletBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );

  // apply flux boundary conditions
  applySourceFluxBC( time_n, dt, dofManager, domain, localMatrix.toViewConstSizes(), localRhs.toView() );
}

void CompositionalMultiphaseFlow::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{

  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();

  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );

  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString(),
                   [&]( FieldSpecificationBase const & fs,
                        string const &,
                        SortedArrayView< localIndex const > const & lset,
                        Group & subRegion,
                        string const & )
  {

    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );

    SortedArray< localIndex > localSet;
    for( localIndex const a : lset )
    {
      if( ghostRank[a] < 0 )
      {
        localSet.insert( a );
      }
    }

    fs.applyBoundaryConditionToSystem< FieldSpecificationAdd,
                                       parallelDevicePolicy<> >( localSet.toViewConst(),
                                                                 time + dt,
                                                                 dt,
                                                                 subRegion,
                                                                 dofNumber,
                                                                 dofManager.rankOffset(),
                                                                 localMatrix,
                                                                 localRhs,
                                                                 [] GEOSX_HOST_DEVICE ( localIndex const )
    {
      return 0.0;
    } );

  } );
}


void CompositionalMultiphaseFlow::applyDirichletBC( real64 const time,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  localIndex const NC = m_numComponents;

  FieldSpecificationManager & fsManager = getGlobalState().getFieldSpecificationManager();

  map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

  // 1. apply pressure Dirichlet BCs
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString(),
                   [&]( FieldSpecificationBase const & fs,
                        string const & setName,
                        SortedArrayView< localIndex const > const & targetSet,
                        Group & subRegion,
                        string const & )
  {
    // 1.0. Check whether pressure has already been applied to this set
    string const & subRegionName = subRegion.getName();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) > 0,
                    "Conflicting pressure boundary conditions on set " << setName );

    bcStatusMap[subRegionName][setName].resize( m_numComponents );
    bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

    // 1.1. Apply BC to set the field values
    fs.applyFieldValue< FieldSpecificationEqual, parallelHostPolicy >( targetSet,
                                                                       time + dt,
                                                                       subRegion,
                                                                       viewKeyStruct::bcPressureString() );
  } );

  // 2. Apply composition BC (global component fraction) and store them for constitutive call
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::globalCompFractionString(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & subRegion,
                         string const & )
  {
    // 2.0. Check pressure and record composition bc application
    string const & subRegionName = subRegion.getName();
    localIndex const comp = fs.getComponent();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0,
                    "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp],
                    "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[subRegionName][setName][comp] = true;

    // 2.1. Apply BC to set the field values
    fs.applyFieldValue< FieldSpecificationEqual, parallelHostPolicy >( targetSet,
                                                                       time + dt,
                                                                       subRegion,
                                                                       viewKeyStruct::globalCompFractionString() );
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
        GEOSX_WARNING_IF( !bcConsistent, "Composition boundary condition not applied to component " << ic <<
                          " on region " << bcStatusEntryOuter.first << ", set " << bcStatusEntryInner.first );
      }
    }
  }
  GEOSX_ERROR_IF( !bcConsistent, "Inconsistent composition boundary conditions" );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );

  // 3. Call constitutive update, back-calculate target global component densities and apply to the system
  fsManager.apply( time + dt,
                   domain,
                   "ElementRegions",
                   viewKeyStruct::pressureString(),
                   [&] ( FieldSpecificationBase const &,
                         string const &,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & subRegion,
                         string const & )
  {
    // TODO: hack! Find a better way to get the fluid
    Group const & region = subRegion.getParent().getParent();
    string const & fluidName = m_fluidModelNames[ targetRegionIndex( region.getName() ) ];
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidName );

    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
    arrayView1d< globalIndex const > const dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< real64 const > const pres      = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const dPres     = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView1d< real64 const > const bcPres    = subRegion.getReference< array1d< real64 > >( viewKeyStruct::bcPressureString() );
    arrayView2d< real64 const > const compFrac  = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
    arrayView2d< real64 const > const compDens  = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const dCompDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );
    arrayView2d< real64 const > const totalDens = fluid.totalDensity();

    constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      // MultiFluid models are not thread-safe or device-capable yet

      IsothermalCompositionalMultiphaseFlowKernels::FluidUpdateKernel::launch< serialPolicy >( targetSet,
                                                                                               fluidWrapper,
                                                                                               bcPres,
                                                                                               m_uniformTemperature,
                                                                                               compFrac );
    } );

    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const ei = targetSet[a];
      if( ghostRank[ei] >= 0 )
        return;

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
        FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic + 1,
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

real64 CompositionalMultiphaseFlow::calculateResidualNorm( DomainPartition const & domain,
                                                           DofManager const & dofManager,
                                                           arrayView1d< real64 const > const & localRhs )
{

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );

  real64 residual;

  if( m_isothermalFlag )
  {
    localIndex const NDOF = m_numComponents + 1;

    real64 localResidualNorm = 0.0;

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
    {
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

      arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & refPoro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
      arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > localSum( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
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
    residual = std::sqrt( MpiWrapper::sum( localResidualNorm ) );

    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[200] = {0};
      sprintf( output, "    ( Rfluid ) = (%4.2e) ; ", residual );
      std::cout<<output;
    }
  }
  else
  {
    localIndex const numFlowDOF = m_numComponents + 1;

    real64 localFlowResidualNorm = 0.0;
    real64 localEnergyResidualNorm = 0.0;

    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
    {
      MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

      arrayView1d< globalIndex const > dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & refPoro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
      arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

      RAJA::ReduceSum< parallelDeviceReduce, real64 > localFlowSum( 0.0 );
      RAJA::ReduceSum< parallelDeviceReduce, real64 > localEnergySum( 0.0 );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        if( elemGhostRank[ei] < 0 )
        {
          localIndex const localRow = dofNumber[ei] - rankOffset;
          real64 const normalizer = totalDens[ei][0] * refPoro[ei] * volume[ei];

          for( localIndex idof = 0; idof < numFlowDOF; ++idof )
          {
            real64 const val = localRhs[localRow + idof] / normalizer;
            localFlowSum += val * val;
          }
          real64 const val = localRhs[localRow + numFlowDOF];
          localEnergySum += val*val;
        }
      } );
      localFlowResidualNorm += localFlowSum.get();
      localEnergyResidualNorm += localEnergySum.get();
    } );

    // compute global residual norm
    real64 const flowResidual   = std::sqrt( MpiWrapper::sum( localFlowResidualNorm ) );
    real64 const energyResidual = std::sqrt( MpiWrapper::sum( localEnergyResidualNorm ) );

    residual = ( flowResidual + energyResidual ) / 2;
    if( getLogLevel() >= 1 && logger::internal::rank==0 )
    {
      char output[200] = {0};
      sprintf( output, "    ( Rflow ) = (%4.2e) - ( Renergy ) = (%4.2e)  ; ", flowResidual, energyResidual );
      std::cout<<output;
    }
  }

  return residual;
}

void CompositionalMultiphaseFlow::solveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

real64 CompositionalMultiphaseFlow::scalingForSystemSolution( DomainPartition const & domain,
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

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );
  real64 scalingFactor = 1.0;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const > const & compDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const & dCompDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    RAJA::ReduceMin< parallelDeviceReduce, real64 > minVal( 1.0 );

    forAll< parallelDevicePolicy<> >( dofNumber.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
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
          real64 const maxAbsCompDensChange = maxCompFracChange * prevTotalDens;

          // This actually checks the change in component fraction, using a lagged total density
          // Indeed we can rewrite the following check as:
          //    | prevCompDens / prevTotalDens - newCompDens / prevTotalDens | > maxCompFracChange
          // Note that the total density in the second term is lagged (i.e, we use prevTotalDens)
          // because I found it more robust than using directly newTotalDens (which can vary also
          // wildly when the compDens change is large)
          if( absCompDensChange > maxAbsCompDensChange && absCompDensChange > eps )
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

  return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}


bool CompositionalMultiphaseFlow::checkSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor )
{
  real64 constexpr eps = minDensForDivision;

  localIndex const NC = m_numComponents;
  integer const allowCompDensChopping = m_allowCompDensChopping;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  globalIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::dofFieldString() );
  int localCheck = 1;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64 const > const & compDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const & dCompDens = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    RAJA::ReduceMin< parallelDeviceReduce, integer > check( 1 );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        localIndex const localRow = dofNumber[ei] - rankOffset;
        {
          real64 const newPres = pres[ei] + dPres[ei] + scalingFactor * localSolution[localRow];
          check.min( newPres >= 0.0 );
        }

        // if component density chopping is not allowed, the time step fails if a component density is negative
        // otherwise, we just check that the total density is positive, and negative component densities
        // will be chopped (i.e., set to zero) in applySystemSolution)
        if( !allowCompDensChopping )
        {
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            check.min( newDens >= 0.0 );
          }
        }
        else
        {
          real64 totalDens = 0.0;
          for( localIndex ic = 0; ic < NC; ++ic )
          {
            real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic] + scalingFactor * localSolution[localRow + ic + 1];
            totalDens += (newDens > 0.0) ? newDens : 0.0;
          }
          check.min( totalDens >= eps );
        }
      }
    } );

    localCheck = std::min( localCheck, check.get() );
  } );

  return MpiWrapper::min( localCheck, MPI_COMM_GEOSX );
}

void CompositionalMultiphaseFlow::applySystemSolution( DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution,
                                                       real64 const scalingFactor,
                                                       DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dofFieldString(),
                               viewKeyStruct::deltaPressureString(),
                               scalingFactor,
                               0, 1 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::dofFieldString(),
                               viewKeyStruct::deltaGlobalCompDensityString(),
                               scalingFactor,
                               1, m_numDofPerCell );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  std::map< string, string_array > fieldNames;

  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaGlobalCompDensityString() ) );

  if( m_isothermalFlag == 0 )
  {
    fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaTemperatureString() ) );
  }

  getGlobalState().getCommunicationTools().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    updateState( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseFlow::chopNegativeDensities( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  localIndex const NC = m_numComponents;
  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    arrayView2d< real64 const > const compDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 > const dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          real64 const newDens = compDens[ei][ic] + dCompDens[ei][ic];
          if( newDens < 0 )
          {
            dCompDens[ei][ic] = -compDens[ei][ic];
          }
        }
      }
    } );
  } );
}

void CompositionalMultiphaseFlow::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  if( m_isothermalFlag )
  {
    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & dPres =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
      arrayView2d< real64 > const & dCompDens =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );
      dPres.setValues< parallelDevicePolicy<> >( 0.0 );
      dCompDens.setValues< parallelDevicePolicy<> >( 0.0 );

      updateState( subRegion, targetIndex );
    } );
  }
  else
  {
    forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 > const & dPres =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
      arrayView2d< real64 > const & dCompDens =
        subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );
      arrayView1d< real64 > const & dTemp =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaTemperatureString() );

      dPres.setValues< parallelDevicePolicy<> >( 0.0 );
      dCompDens.setValues< parallelDevicePolicy<> >( 0.0 );
      dTemp.setValues< parallelDevicePolicy<> >( 0.0 );

      updateState( subRegion, targetIndex );
    } );
  }
}

void CompositionalMultiphaseFlow::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{
  localIndex const NC = m_numComponents;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64 const > const dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    arrayView1d< real64 > const pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView2d< real64 > const compDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        compDens[ei][ic] += dCompDens[ei][ic];
      }
    } );
  } );

  if( m_isothermalFlag == 0 )
  {
    forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
    {
      arrayView1d< real64 const > const dTemp =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaTemperatureString() );
      arrayView1d< real64 > const temp =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::temperatureString() );

      forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
      {
        temp[ei] += dTemp[ei];
      } );
    } );
  }
}

void CompositionalMultiphaseFlow::resetViews( MeshLevel & mesh )
{
  FlowSolverBase::resetViews( mesh );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::pressureString() );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString() );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaPressureString() );
  m_deltaPressure.setName( getName() + "/accessors/" + viewKeyStruct::deltaPressureString() );

  m_dCompFrac_dCompDens.clear();
  m_dCompFrac_dCompDens = elemManager.constructArrayViewAccessor< real64, 3 >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );
  m_dCompFrac_dCompDens.setName( getName() + "/accessors/" + viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  m_phaseVolFrac.clear();
  m_phaseVolFrac = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::phaseVolumeFractionString() );
  m_phaseVolFrac.setName( getName() + "/accessors/" + viewKeyStruct::phaseVolumeFractionString() );

  m_dPhaseVolFrac_dPres.clear();
  m_dPhaseVolFrac_dPres = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
  m_dPhaseVolFrac_dPres.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseVolumeFraction_dPressureString() );

  m_dPhaseVolFrac_dCompDens.clear();
  m_dPhaseVolFrac_dCompDens = elemManager.constructArrayViewAccessor< real64, 3 >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );
  m_dPhaseVolFrac_dCompDens.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  m_phaseMob.clear();
  m_phaseMob = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::phaseMobilityString() );
  m_phaseMob.setName( getName() + "/accessors/" + viewKeyStruct::phaseMobilityString() );

  m_dPhaseMob_dPres.clear();
  m_dPhaseMob_dPres = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::dPhaseMobility_dPressureString() );
  m_dPhaseMob_dPres.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseMobility_dPressureString() );

  m_dPhaseMob_dCompDens.clear();
  m_dPhaseMob_dCompDens = elemManager.constructArrayViewAccessor< real64, 3 >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() );
  m_dPhaseMob_dCompDens.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() );

  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_phaseMassDens.clear();
    m_phaseMassDens = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseMassDensityString(),
                                                                                   targetRegionNames(),
                                                                                   fluidModelNames() );
    m_phaseMassDens.setName( getName() + "/accessors/" + keys::phaseMassDensityString() );

    m_dPhaseMassDens_dPres.clear();
    m_dPhaseMassDens_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseMassDensity_dPressureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseMassDens_dPres.setName( getName() + "/accessors/" + keys::dPhaseMassDensity_dPressureString() );

    m_dPhaseMassDens_dComp.clear();
    m_dPhaseMassDens_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseMassDensity_dGlobalCompFractionString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseMassDens_dComp.setName( getName() + "/accessors/" + keys::dPhaseMassDensity_dGlobalCompFractionString() );

    m_phaseCompFrac.clear();
    m_phaseCompFrac = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::phaseCompFractionString(),
                                                                                   targetRegionNames(),
                                                                                   fluidModelNames() );
    m_phaseCompFrac.setName( getName() + "/accessors/" + keys::phaseCompFractionString() );

    m_dPhaseCompFrac_dPres.clear();
    m_dPhaseCompFrac_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseCompFraction_dPressureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseCompFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dPressureString() );

    m_dPhaseCompFrac_dComp.clear();
    m_dPhaseCompFrac_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 5 >( keys::dPhaseCompFraction_dGlobalCompFractionString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseCompFrac_dComp.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dGlobalCompFractionString() );
  }
  if( m_capPressureFlag )
  {
    using keys = CapillaryPressureBase::viewKeyStruct;

    m_phaseCapPressure.clear();
    m_phaseCapPressure = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseCapPressureString(),
                                                                                      targetRegionNames(),
                                                                                      capPresModelNames() );
    m_phaseCapPressure.setName( getName() + "/accessors/" + keys::phaseCapPressureString() );

    m_dPhaseCapPressure_dPhaseVolFrac.clear();
    m_dPhaseCapPressure_dPhaseVolFrac = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseCapPressure_dPhaseVolFractionString(),
                                                                                                     targetRegionNames(),
                                                                                                     capPresModelNames() );
    m_dPhaseCapPressure_dPhaseVolFrac.setName( getName() + "/accessors/" + keys::dPhaseCapPressure_dPhaseVolFractionString() );
  }
  if( m_isothermalFlag == 0 )
  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_temperature.clear();
    m_temperature = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::temperatureString() );
    m_temperature.setName( getName() + "/accessors/" + viewKeyStruct::temperatureString() );

    m_deltaTemperature.clear();
    m_deltaTemperature = elemManager.constructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaTemperatureString() );
    m_deltaTemperature.setName( getName() + "/accessors/" + viewKeyStruct::deltaTemperatureString() );

    m_dPhaseVolFrac_dTemp.clear();
    m_dPhaseVolFrac_dTemp = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );
    m_dPhaseVolFrac_dTemp.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseVolumeFraction_dTemperatureString() );

    m_dPhaseMob_dTemp.clear();
    m_dPhaseMob_dTemp = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::dPhaseMobility_dTemperatureString() );
    m_dPhaseMob_dTemp.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseMobility_dTemperatureString() );

    m_dPhaseMassDens_dTemp.clear();
    m_dPhaseMassDens_dTemp = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseMassDensity_dTemperatureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseMassDens_dTemp.setName( getName() + "/accessors/" + keys::dPhaseMassDensity_dTemperatureString() );

    m_dPhaseCompFrac_dTemp.clear();
    m_dPhaseCompFrac_dTemp = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseCompFraction_dTemperatureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseCompFrac_dTemp.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dTemperatureString() );

    m_dPhaseMob_dTemp.clear();
    m_dPhaseMob_dTemp = elemManager.constructArrayViewAccessor< real64, 2 >( viewKeyStruct::dPhaseMobility_dTemperatureString() );
    m_dPhaseMob_dTemp.setName( getName() + "/accessors/" + viewKeyStruct::dPhaseMobility_dTemperatureString() );


    m_phaseEnthalpy.clear();
    m_phaseEnthalpy = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseEnthalpyString(),
                                                                                   targetRegionNames(),
                                                                                   fluidModelNames() );
    m_phaseEnthalpy.setName( getName() + "/accessors/" + keys::phaseEnthalpyString() );

    m_dPhaseEnthalpy_dPres.clear();
    m_dPhaseEnthalpy_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseEnthalpy_dPressureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseEnthalpy_dPres.setName( getName() + "/accessors/" + keys::dPhaseEnthalpy_dPressureString() );

    m_dPhaseEnthalpy_dTemp.clear();
    m_dPhaseEnthalpy_dTemp = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseEnthalpy_dTemperatureString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseEnthalpy_dTemp.setName( getName() + "/accessors/" + keys::dPhaseEnthalpy_dTemperatureString() );

    m_dPhaseEnthalpy_dComp.clear();
    m_dPhaseEnthalpy_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseEnthalpy_dGlobalCompFractionString(),
                                                                                          targetRegionNames(),
                                                                                          fluidModelNames() );
    m_dPhaseEnthalpy_dComp.setName( getName() + "/accessors/" + keys::dPhaseEnthalpy_dGlobalCompFractionString() );



  }
}

//START_SPHINX_INCLUDE_01
REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseFlow, string const &, Group * const )

} // namespace geosx
