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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"
#include "mesh/DomainPartition.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "mesh/PerforationData.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/SinglePhaseWellKernels.hpp"
#include "physicsSolvers/fluidFlow/wells/WellControls.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseWellKernels;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                          Group * const parent )
  :
  WellSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_temperature( 0.0 ),
  m_useMass( false ),
  m_maxCompFracChange( 1.0 ),
  m_maxRelativePresChange( 0.2 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 ),
  m_targetPhaseIndex( -1 )
{
  this->registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );

  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::relPermNamesString(), &m_relPermModelNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Names of relative permeability constitutive models to use" );

  this->registerWrapper( viewKeyStruct::maxCompFracChangeString(), &m_maxCompFracChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (absolute) change in a component fraction between two Newton iterations" );

  this->registerWrapper( viewKeyStruct::maxRelativePresChangeString(), &m_maxRelativePresChange ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0 ).
    setDescription( "Maximum (relative) change in pressure between two Newton iterations (recommended with rate control)" );

  this->registerWrapper( viewKeyStruct::allowLocalCompDensChoppingString(), &m_allowCompDensChopping ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1 ).
    setDescription( "Flag indicating whether local (cell-wise) chopping of negative compositions is allowed" );

}

void CompositionalMultiphaseWell::postProcessInput()
{
  WellSolverBase::postProcessInput();
  checkModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString() );

  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );
  GEOSX_UNUSED_VAR( flowSolver );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );

}

void CompositionalMultiphaseWell::registerDataOnMesh( Group & meshBodies )
{
  WellSolverBase::registerDataOnMesh( meshBodies );

  MeshLevel & meshLevel = meshBodies.getGroup< MeshBody >( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompDensityString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mixtureConnRateString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::globalCompFractionString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    subRegion.registerWrapper< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalMassDensityString() ).setPlotLevel( PlotLevel::LEVEL_0 );
    subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    PerforationData & perforationData = *subRegion.getPerforationData();
    perforationData.registerWrapper< array2d< real64 > >( viewKeyStruct::compPerforationRateString() );
    perforationData.registerWrapper< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    perforationData.registerWrapper< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString() ).
      setRestartFlags( RestartFlags::NO_WRITE );

    WellControls & wellControls = getWellControls( subRegion );
    wellControls.registerWrapper< real64 >( viewKeyStruct::currentBHPString() );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentBHP_dPresString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentBHP_dCompDensString() ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );

    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dPresString() ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dRateString() ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );

    wellControls.registerWrapper< real64 >( viewKeyStruct::currentTotalVolRateString() );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dPresString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
    wellControls.registerWrapper< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRate_dCompDensString() ).
      setRestartFlags( RestartFlags::NO_WRITE ).
      setSizedFromParent( 0 );
    wellControls.registerWrapper< real64 >( viewKeyStruct::dCurrentTotalVolRate_dRateString() ).
      setRestartFlags( RestartFlags::NO_WRITE );
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

void CompositionalMultiphaseWell::validateConstitutiveModels( MeshLevel const & meshLevel, ConstitutiveManager const & cm ) const
{
  CompositionalMultiphaseBase const & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  arrayView1d< string const > const & flowTargetRegionNames = flowSolver.targetRegionNames();
  arrayView1d< string const > const & flowFluidModels = flowSolver.fluidModelNames();
  arrayView1d< string const > const & flowRelPermModels = flowSolver.relPermModelNames();

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // Make a set of unique reservoir region indices the well is connected to
    PerforationData const & perforationData = *subRegion.getPerforationData();
    arrayView1d< localIndex const > const & resElementRegion = perforationData.getMeshElements().m_toElementRegion;
    std::set< string > reservoirRegionNames;
    for( localIndex const ei : resElementRegion )
    {
      reservoirRegionNames.insert( meshLevel.getElemManager().getRegion( ei ).getName() );
    }

    // Check that each well model is compatible with all models in perforated reservoir regions
    MultiFluidBase const & wellFluid = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[targetIndex] );
    RelativePermeabilityBase const & wellRelPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[targetIndex] );
    for( localIndex resTargetIndex = 0; resTargetIndex < flowTargetRegionNames.size(); ++resTargetIndex )
    {
      if( reservoirRegionNames.count( flowTargetRegionNames[resTargetIndex] ) > 0 )
      {
        MultiFluidBase const & resFluid = cm.getConstitutiveRelation< MultiFluidBase >( flowFluidModels[resTargetIndex] );
        CompareMultiphaseModels( wellFluid, resFluid );
        CompareMulticomponentModels( wellFluid, resFluid );

        RelativePermeabilityBase const & resRelPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( flowRelPermModels[resTargetIndex] );
        CompareMultiphaseModels( wellRelPerm, resRelPerm );
      }
    }
  } );
}

void CompositionalMultiphaseWell::validateInjectionStreams( MeshLevel const & meshLevel ) const
{
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = getWellControls( subRegion );

    // check well injection stream for injectors
    if( wellControls.getType() == WellControls::Type::INJECTOR )
    {
      arrayView1d< real64 const > const & injection = wellControls.getInjectionStream();
      real64 compFracSum = 0;
      for( localIndex ic = 0; ic < m_numComponents; ++ic )
      {
        real64 const compFrac = injection[ic];
        GEOSX_ERROR_IF( ( compFrac < 0.0 ) || ( compFrac > 1.0 ),
                        "Invalid injection stream for well " << subRegion.getName() );
        compFracSum += compFrac;
      }
      GEOSX_ERROR_IF( ( compFracSum < 1.0 - std::numeric_limits< real64 >::epsilon() ) ||
                      ( compFracSum > 1.0 + std::numeric_limits< real64 >::epsilon() ),
                      "Invalid injection stream for well " << subRegion.getName() );
    }
  } );
}

void CompositionalMultiphaseWell::validateWellConstraints( MeshLevel const & meshLevel, MultiFluidBase const & fluid )
{
  // now that we know we are single-phase, we can check a few things in the constraints
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = getWellControls( subRegion );
    WellControls::Type const wellType = wellControls.getType();
    WellControls::Control const currentControl = wellControls.getControl();
    real64 const & targetTotalRate = wellControls.getTargetTotalRate();
    real64 const & targetPhaseRate = wellControls.getTargetPhaseRate();
    integer const useSurfaceConditions = wellControls.useSurfaceConditions();
    real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

    GEOSX_ERROR_IF( wellType == WellControls::Type::INJECTOR && currentControl == WellControls::Control::PHASEVOLRATE,
                    "Phase rate control is not available for injectors" );
    GEOSX_ERROR_IF( wellType == WellControls::Type::PRODUCER && currentControl == WellControls::Control::TOTALVOLRATE,
                    "Phase rate control is not available for producers" );

    GEOSX_ERROR_IF( wellType == WellControls::Type::INJECTOR && isZero( targetTotalRate ),
                    "Target total rate cannot be equal to zero for injectors" );
    GEOSX_ERROR_IF( wellType == WellControls::Type::INJECTOR && !isZero( targetPhaseRate ),
                    "Target phase rate cannot be used for injectors" );

    GEOSX_ERROR_IF( wellType == WellControls::Type::PRODUCER && isZero( targetPhaseRate ),
                    "Target phase rate cannot be equal to zero for producers" );
    GEOSX_ERROR_IF( wellType == WellControls::Type::PRODUCER && !isZero( targetTotalRate ),
                    "Target total rate cannot be used for producers" );

    GEOSX_ERROR_IF( useSurfaceConditions && surfaceTemp <= 0,
                    "Surface temperature must be set to a strictly positive value" );

    // Find target phase index for phase rate constraint
    for( localIndex ip = 0; ip < fluid.numFluidPhases(); ++ip )
    {
      if( fluid.phaseNames()[ip] == wellControls.getTargetPhaseName() )
      {
        m_targetPhaseIndex = ip;
      }
    }
    GEOSX_ERROR_IF( wellType == WellControls::Type::PRODUCER && m_targetPhaseIndex == -1,
                    "Phase " << wellControls.getTargetPhaseName() << " not found for well control " << wellControls.getName() );

  } );

}

void CompositionalMultiphaseWell::initializePreSubGroups()
{
  WellSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  validateConstitutiveModels( meshLevel, cm );

  // Set key dimensions (phases, components) from one of the fluids - they should all be compatible
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases     = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerWellElement = m_numComponents + 2; // 1 pressure + NC compositions + 1 connectionRate

  validateModelMapping< MultiFluidBase >( meshLevel.getElemManager(), m_fluidModelNames );
  validateModelMapping< RelativePermeabilityBase >( meshLevel.getElemManager(), m_relPermModelNames );
  validateInjectionStreams( meshLevel );
  validateWellConstraints( meshLevel, fluid0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    resizeFields( subRegion );
  } );
}

void CompositionalMultiphaseWell::resizeFields( WellElementSubRegion & subRegion )
{
  PerforationData & perforationData = *subRegion.getPerforationData();
  WellControls & wellControls = getWellControls( subRegion );

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() ).resizeDimension< 1 >( NC );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() ).resizeDimension< 1 >( NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() ).resizeDimension< 1 >( NC );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).resizeDimension< 1, 2 >( NC, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() ).resizeDimension< 1 >( NP );
  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).resizeDimension< 1 >( NP );
  subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).resizeDimension< 1, 2 >( NP, NC );

  subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() ).resizeDimension< 1 >( NC );

  perforationData.getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString() ).resizeDimension< 1 >( NC );
  perforationData.getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString() ).resizeDimension< 1, 2 >( 2, NC );
  perforationData.getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString() ).resizeDimension< 1, 2, 3 >( 2, NC, NC );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentBHP_dCompDensString() ).resizeDimension< 0 >( NC );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::currentPhaseVolRateString() ).resizeDimension< 0 >( NP );
  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dPresString() ).resizeDimension< 0 >( NP );
  wellControls.getReference< array2d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() ).resizeDimension< 0, 1 >( NP, NC );
  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentPhaseVolRate_dRateString() ).resizeDimension< 0 >( NP );

  wellControls.getReference< array1d< real64 > >( viewKeyStruct::dCurrentTotalVolRate_dCompDensString() ).resizeDimension< 0 >( NC );
}

void CompositionalMultiphaseWell::initializePostInitialConditionsPreSubGroups()
{
  WellSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );
}

void CompositionalMultiphaseWell::updateComponentFraction( WellElementSubRegion & subRegion ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs
  arrayView2d< real64 > const & compFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
  arrayView3d< real64 > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  // inputs
  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  CompositionalMultiphaseBaseKernels::KernelLaunchSelector1< CompositionalMultiphaseBaseKernels::ComponentFractionKernel
                                                             >( numFluidComponents(),
                                                                subRegion.size(),
                                                                compDens,
                                                                dCompDens,
                                                                compFrac,
                                                                dCompFrac_dCompDens );

}

void CompositionalMultiphaseWell::updateBHPForConstraint( WellElementSubRegion & subRegion,
                                                          localIndex const targetIndex )
{
  GEOSX_UNUSED_VAR( targetIndex );

  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex const NC = m_numComponents;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView1d< real64 > const & totalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString() );
  arrayView1d< real64 > const & dTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString() );
  arrayView2d< real64 > const & dTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() );

  arrayView1d< real64 const > const wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  real64 const & refGravCoef = wellControls.getReferenceGravityCoef();

  real64 & currentBHP =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentBHPString() );
  real64 & dCurrentBHP_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dPresString() );
  arrayView1d< real64 > const & dCurrentBHP_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentBHP_dCompDensString() );

  // bring everything back to host, capture the scalars by reference
  forAll< serialPolicy >( 1, [&NC,
                              pres,
                              dPres,
                              totalMassDens,
                              dTotalMassDens_dPres,
                              dTotalMassDens_dCompDens,
                              wellElemGravCoef,
                              &currentBHP,
                              &dCurrentBHP_dPres,
                              dCurrentBHP_dCompDens,
                              &iwelemRef,
                              &refGravCoef] ( localIndex const )
  {
    real64 const diffGravCoef = refGravCoef - wellElemGravCoef[iwelemRef];
    currentBHP = pres[iwelemRef] + dPres[iwelemRef] + totalMassDens[iwelemRef] * diffGravCoef;
    dCurrentBHP_dPres = 1 + dTotalMassDens_dPres[iwelemRef] * diffGravCoef;
    for( localIndex ic = 0; ic < NC; ++ic )
    {
      dCurrentBHP_dCompDens[ic] = dTotalMassDens_dCompDens[iwelemRef][ic] * diffGravCoef;
    }
  } );
}

void CompositionalMultiphaseWell::updateVolRatesForConstraint( WellElementSubRegion & subRegion,
                                                               localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  // the rank that owns the reference well element is responsible for the calculations below.
  if( !subRegion.isLocallyOwned() )
  {
    return;
  }

  localIndex constexpr maxNumComp = constitutive::MultiFluidBase::MAX_NUM_COMPONENTS;
  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;
  localIndex const iwelemRef = subRegion.getTopWellElementIndex();
  real64 const temp = m_temperature;

  // subRegion data

  arrayView1d< real64 const > const & pres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  arrayView1d< real64 const > const & connRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString() );
  arrayView1d< real64 const > const & dConnRate =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString() );

  arrayView2d< real64 const > const & compFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  // fluid data

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  // control data

  WellControls & wellControls = getWellControls( subRegion );

  integer const useSurfaceConditions = wellControls.useSurfaceConditions();
  real64 const & surfacePres = wellControls.getSurfacePressure();
  real64 const & surfaceTemp = wellControls.getSurfaceTemperature();

  arrayView1d< real64 > const & currentPhaseVolRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::currentPhaseVolRateString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dPres =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dPresString() );
  arrayView2d< real64 > const & dCurrentPhaseVolRate_dCompDens =
    wellControls.getReference< array2d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dCompDensString() );
  arrayView1d< real64 > const & dCurrentPhaseVolRate_dRate =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentPhaseVolRate_dRateString() );

  real64 & currentTotalVolRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::currentTotalVolRateString() );
  real64 & dCurrentTotalVolRate_dPres =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dPresString() );
  arrayView1d< real64 > const & dCurrentTotalVolRate_dCompDens =
    wellControls.getReference< array1d< real64 > >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dCompDensString() );
  real64 & dCurrentTotalVolRate_dRate =
    wellControls.getReference< real64 >( CompositionalMultiphaseWell::viewKeyStruct::dCurrentTotalVolRate_dRateString() );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    // bring everything back to host, capture the scalars by reference
    forAll< serialPolicy >( 1, [&NC,
                                &NP,
                                fluidWrapper,
                                pres,
                                dPres,
                                compDens,
                                dCompDens,
                                compFrac,
                                dCompFrac_dCompDens,
                                connRate,
                                dConnRate,
                                phaseDens,
                                dPhaseDens_dPres,
                                dPhaseDens_dComp,
                                phaseFrac,
                                dPhaseFrac_dPres,
                                dPhaseFrac_dComp,
                                &useSurfaceConditions,
                                &surfacePres,
                                &surfaceTemp,
                                &temp,
                                &currentTotalVolRate,
                                &dCurrentTotalVolRate_dPres,
                                dCurrentTotalVolRate_dCompDens,
                                &dCurrentTotalVolRate_dRate,
                                currentPhaseVolRate,
                                dCurrentPhaseVolRate_dPres,
                                dCurrentPhaseVolRate_dCompDens,
                                dCurrentPhaseVolRate_dRate,
                                &iwelemRef] ( localIndex const )
    {
      //    We need to evaluate the density as follows:
      //      - Surface conditions: using the surface pressure provided by the user
      //      - Reservoir conditions: using the pressure in the top element

      if( useSurfaceConditions )
      {
        // we need to compute the surface density
        fluidWrapper.update( iwelemRef, 0, surfacePres, surfaceTemp, compFrac[iwelemRef] );
      }
      else
      {
        real64 const refPres = pres[iwelemRef] + dPres[iwelemRef];
        fluidWrapper.update( iwelemRef, 0, refPres, temp, compFrac[iwelemRef] );
      }

      // total volume rate
      real64 const currentTotalRate = connRate[iwelemRef] + dConnRate[iwelemRef];
      real64 totalDens = 0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        totalDens += compDens[iwelemRef][ic] + dCompDens[iwelemRef][ic];
      }
      real64 const totalDensInv = 1.0 / totalDens;
      currentTotalVolRate = currentTotalRate * totalDensInv;
      dCurrentTotalVolRate_dPres = 0; // using the fact that totalDens = \sum_c compDens (independent of pressure)
      dCurrentTotalVolRate_dRate = totalDensInv;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dCurrentTotalVolRate_dCompDens[ic] = -currentTotalVolRate * totalDensInv; // using the fact that totalDens = \sum_c compDens
      }

      // phase volume rate
      stackArray1d< real64, maxNumComp > work( NC );
      for( localIndex ip = 0; ip < NP; ++ip )
      {
        real64 const phaseDensInv = 1.0 / phaseDens[iwelemRef][0][ip];
        real64 const phaseFracTimesPhaseDensInv = phaseFrac[iwelemRef][0][ip] * phaseDensInv;
        real64 const dPhaseFracTimesPhaseDensInv_dPres = dPhaseFrac_dPres[iwelemRef][0][ip] * phaseDensInv
                                                         - dPhaseDens_dPres[iwelemRef][0][ip] * phaseFracTimesPhaseDensInv * phaseDensInv;

        currentPhaseVolRate[ip] = currentTotalRate * phaseFracTimesPhaseDensInv;
        dCurrentPhaseVolRate_dPres[ip] = ( useSurfaceConditions ==  0 ) * currentTotalRate * dPhaseFracTimesPhaseDensInv_dPres;
        dCurrentPhaseVolRate_dRate[ip] = phaseFracTimesPhaseDensInv;
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          dCurrentPhaseVolRate_dCompDens[ip][ic] = -phaseFracTimesPhaseDensInv * dPhaseDens_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] += dPhaseFrac_dComp[iwelemRef][0][ip][ic] * phaseDensInv;
          dCurrentPhaseVolRate_dCompDens[ip][ic] *= currentTotalRate;
        }
        applyChainRuleInPlace( NC, dCompFrac_dCompDens[iwelemRef], dCurrentPhaseVolRate_dCompDens[ip], work.data() );
      }
    } );
  } );
}


void CompositionalMultiphaseWell::updateFluidModel( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & dPres = subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
  arrayView2d< real64 const > const & compFrac = subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    CompositionalMultiphaseBaseKernels::FluidUpdateKernel::launch< ExecPolicy >( subRegion.size(),
                                                                                 fluidWrapper,
                                                                                 pres,
                                                                                 dPres,
                                                                                 m_temperature,
                                                                                 compFrac );
  } );
}

void CompositionalMultiphaseWell::updatePhaseVolumeFraction( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64 > const & phaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );
  arrayView2d< real64 > const & dPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
  arrayView3d< real64 > const & dPhaseVolFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  // inputs

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );
  arrayView2d< real64 const > const & compDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
  arrayView2d< real64 const > const & dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  CompositionalMultiphaseBaseKernels::KernelLaunchSelector2< CompositionalMultiphaseBaseKernels::PhaseVolumeFractionKernel
                                                             >( numFluidComponents(), numFluidPhases(),
                                                                subRegion.size(),
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
                                                                dPhaseVolFrac_dCompDens );
}

void CompositionalMultiphaseWell::updateTotalMassDensity( WellElementSubRegion & subRegion, localIndex const targetIndex ) const
{
  // outputs

  arrayView1d< real64 > const & totalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString() );
  arrayView1d< real64 > const & dTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString() );
  arrayView2d< real64 > const & dTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() );

  // inputs

  arrayView2d< real64 const > const & phaseVolFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );
  arrayView2d< real64 const > const & dPhaseVolFrac_dPres =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
  arrayView3d< real64 const > const & dPhaseVolFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  arrayView3d< real64 const > const & dCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const > const & phaseMassDens = fluid.phaseMassDensity();
  arrayView3d< real64 const > const & dPhaseMassDens_dPres = fluid.dPhaseMassDensity_dPressure();
  arrayView4d< real64 const > const & dPhaseMassDens_dComp = fluid.dPhaseMassDensity_dGlobalCompFraction();

  TotalMassDensityKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                            numFluidComponents(),
                                                            numFluidPhases(),
                                                            phaseVolFrac,
                                                            dPhaseVolFrac_dPres,
                                                            dPhaseVolFrac_dCompDens,
                                                            dCompFrac_dCompDens,
                                                            phaseMassDens,
                                                            dPhaseMassDens_dPres,
                                                            dPhaseMassDens_dComp,
                                                            totalMassDens,
                                                            dTotalMassDens_dPres,
                                                            dTotalMassDens_dCompDens );

}

void CompositionalMultiphaseWell::updateSubRegionState( WellElementSubRegion & subRegion, localIndex const targetIndex )
{
  // update properties
  updateComponentFraction( subRegion );

  // update volumetric rates for the well constraints
  // note: this must be called before updateFluidModel
  updateVolRatesForConstraint( subRegion, targetIndex );

  // update densities, phase fractions, phase volume fractions
  updateFluidModel( subRegion, targetIndex );
  updatePhaseVolumeFraction( subRegion, targetIndex );
  updateTotalMassDensity( subRegion, targetIndex );

  // update the current BHP pressure
  updateBHPForConstraint( subRegion, targetIndex );

  // update perforation rates
  computePerforationRates( subRegion, targetIndex );
}

void CompositionalMultiphaseWell::initializeWells( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  localIndex const NC = m_numComponents;
  localIndex const NP = m_numPhases;

  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion & subRegion )
  {

    WellControls const & wellControls = getWellControls( subRegion );
    PerforationData const & perforationData = *subRegion.getPerforationData();

    // get well primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView2d< real64 > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString() );

    // get the info stored on well elements
    arrayView2d< real64 > const & wellElemCompFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

    // get the element region, subregion, index
    arrayView1d< localIndex const > const & resElementRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
    arrayView1d< localIndex const > const & resElementSubRegion =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
    arrayView1d< localIndex const > const & resElementIndex =
      perforationData.getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

    // 1) Loop over all perforations to compute an average mixture density and component fraction
    // 2) Initialize the reference pressure
    // 3) Estimate the pressures in the well elements using the average density
    PresCompFracInitializationKernel::launch< parallelDevicePolicy<> >( perforationData.size(),
                                                                        subRegion.size(),
                                                                        NC,
                                                                        NP,
                                                                        perforationData.getNumPerforationsGlobal(),
                                                                        wellControls,
                                                                        m_resPres.toNestedViewConst(),
                                                                        m_resCompDens.toNestedViewConst(),
                                                                        m_resPhaseVolFrac.toNestedViewConst(),
                                                                        m_resPhaseMassDens.toNestedViewConst(),
                                                                        resElementRegion,
                                                                        resElementSubRegion,
                                                                        resElementIndex,
                                                                        wellElemGravCoef,
                                                                        wellElemPressure,
                                                                        wellElemCompFrac );

    // get well secondary variables on well elements
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView3d< real64 const > const & wellElemPhaseDens = fluid.phaseDensity();
    arrayView2d< real64 const > const & wellElemTotalDens = fluid.totalDensity();

    // 4) Back calculate component densities
    constitutive::constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      typename TYPEOFREF( castedFluid ) ::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      CompositionalMultiphaseBaseKernels::FluidUpdateKernel::launch< serialPolicy >( subRegion.size(),
                                                                                     fluidWrapper,
                                                                                     wellElemPressure,
                                                                                     m_temperature,
                                                                                     wellElemCompFrac );
    } );

    CompDensInitializationKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                    NC,
                                                                    wellElemCompFrac,
                                                                    wellElemTotalDens,
                                                                    wellElemCompDens );

    // 5) Recompute the pressure-dependent properties
    // Note: I am leaving that here because I would like to use the perforationRates (computed in UpdateState)
    //       to better initialize the rates
    updateSubRegionState( subRegion, targetIndex );

    // 6) Estimate the well rates
    // TODO: initialize rates using perforation rates
    CompositionalMultiphaseWellKernels::RateInitializationKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                                                                    m_targetPhaseIndex,
                                                                                                    wellControls,
                                                                                                    wellElemPhaseDens,
                                                                                                    wellElemTotalDens,
                                                                                                    connRate );
  } );
}

void CompositionalMultiphaseWell::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                     real64 const dt,
                                                     DomainPartition const & domain,
                                                     DofManager const & dofManager,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  // saved current dt for residual normalization
  m_currentDt = dt;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // loop over the wells
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    WellControls const & wellControls = getWellControls( subRegion );

    // get a reference to the degree-of-freedom numbers
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString() );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString() );

    // get the info stored on well elements
    arrayView2d< real64 const > const & wellElemCompFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
    arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

    FluxKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                  dofManager.rankOffset(),
                                                  numFluidComponents(),
                                                  numDofPerResElement(),
                                                  wellControls,
                                                  wellElemDofNumber,
                                                  nextWellElemIndex,
                                                  connRate,
                                                  dConnRate,
                                                  wellElemCompFrac,
                                                  dWellElemCompFrac_dCompDens,
                                                  dt,
                                                  localMatrix,
                                                  localRhs );
  } );
}

void CompositionalMultiphaseWell::assembleVolumeBalanceTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                                              real64 const GEOSX_UNUSED_PARAM( dt ),
                                                              DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degrees of freedom and ghosting info
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get the properties on the well element
    arrayView2d< real64 const > const & wellElemPhaseVolFrac =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const > const & dWellElemPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
    arrayView3d< real64 const > const & dWellElemPhaseVolFrac_dComp =
      subRegion.getReference< array3d< real64 > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

    arrayView1d< real64 const > const & wellElemVolume =
      subRegion.getReference< array1d< real64 > >( ElementSubRegionBase::viewKeyStruct::elementVolumeString() );

    VolumeBalanceKernel::launch< parallelDevicePolicy<> >( subRegion.size(),
                                                           numFluidComponents(),
                                                           numFluidPhases(),
                                                           numDofPerWellElement(),
                                                           dofManager.rankOffset(),
                                                           wellElemDofNumber,
                                                           wellElemGhostRank,
                                                           wellElemPhaseVolFrac,
                                                           dWellElemPhaseVolFrac_dPres,
                                                           dWellElemPhaseVolFrac_dComp,
                                                           wellElemVolume,
                                                           localMatrix,
                                                           localRhs );
  } );
}


real64
CompositionalMultiphaseWell::calculateResidualNorm( DomainPartition const & domain,
                                                    DofManager const & dofManager,
                                                    arrayView1d< real64 const > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 localResidualNorm = 0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const targetIndex,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    arrayView3d< real64 const > const & phaseDens = fluid.phaseDensity();
    arrayView2d< real64 const > const & totalDens = fluid.totalDensity();

    WellControls const & wellControls = getWellControls( subRegion );

    ResidualNormKernel::launch< parallelDevicePolicy<>,
                                parallelDeviceReduce >( localRhs,
                                                        dofManager.rankOffset(),
                                                        subRegion.isLocallyOwned(),
                                                        subRegion.getTopWellElementIndex(),
                                                        numDofPerWellElement(),
                                                        m_targetPhaseIndex,
                                                        wellControls,
                                                        wellElemDofNumber,
                                                        wellElemGhostRank,
                                                        phaseDens,
                                                        totalDens,
                                                        m_currentDt,
                                                        &localResidualNorm );
  } );
  return sqrt( MpiWrapper::sum( localResidualNorm, MPI_COMM_GEOSX ) );
}

real64
CompositionalMultiphaseWell::scalingForSystemSolution( DomainPartition const & domain,
                                                       DofManager const & dofManager,
                                                       arrayView1d< real64 const > const & localSolution )
{
  GEOSX_MARK_FUNCTION;

  // check if we want to rescale the Newton update
  if( m_maxCompFracChange >= 1.0 && m_maxRelativePresChange >= 1.0 )
  {
    // no rescaling wanted, we just return 1.0;
    return 1.0;
  }

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  real64 scalingFactor = 1.0;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers on well elements and ghosting info
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    real64 const subRegionScalingFactor =
      SolutionScalingKernel::launch< parallelDevicePolicy<>,
                                     parallelDeviceReduce >( localSolution,
                                                             dofManager.rankOffset(),
                                                             numFluidComponents(),
                                                             wellElemDofNumber,
                                                             wellElemGhostRank,
                                                             wellElemPressure,
                                                             dWellElemPressure,
                                                             wellElemCompDens,
                                                             dWellElemCompDens,
                                                             m_maxRelativePresChange,
                                                             m_maxCompFracChange );


    if( subRegionScalingFactor < scalingFactor )
    {
      scalingFactor = subRegionScalingFactor;
    }
  } );

  return LvArray::math::max( MpiWrapper::min( scalingFactor, MPI_COMM_GEOSX ), m_minScalingFactor );
}

bool
CompositionalMultiphaseWell::checkSystemSolution( DomainPartition const & domain,
                                                  DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  int localCheck = 1;
  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {
    // get the degree of freedom numbers on well elements and ghosting info
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    // get a reference to the primary variables on well elements
    arrayView1d< real64 const > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    localIndex const subRegionSolutionCheck =
      SolutionCheckKernel::launch< parallelDevicePolicy<>,
                                   parallelDeviceReduce >( localSolution,
                                                           dofManager.rankOffset(),
                                                           numFluidComponents(),
                                                           wellElemDofNumber,
                                                           wellElemGhostRank,
                                                           wellElemPressure,
                                                           dWellElemPressure,
                                                           wellElemCompDens,
                                                           dWellElemCompDens,
                                                           m_allowCompDensChopping,
                                                           scalingFactor );

    if( subRegionSolutionCheck == 0 )
    {
      localCheck = 0;
    }
  } );
  return MpiWrapper::min( localCheck );
}

void CompositionalMultiphaseWell::computePerforationRates( WellElementSubRegion & subRegion,
                                                           localIndex const GEOSX_UNUSED_PARAM( targetIndex ) )
{
  GEOSX_MARK_FUNCTION;

  PerforationData * const perforationData = subRegion.getPerforationData();

  // get depth
  arrayView1d< real64 const > const & wellElemGravCoef =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );

  // get well primary variables on well elements
  arrayView1d< real64 const > const & wellElemPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const & dWellElemPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  arrayView2d< real64 const > const & wellElemCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
  arrayView2d< real64 const > const & dWellElemCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

  arrayView1d< real64 const > const & wellElemTotalMassDens =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString() );
  arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString() );
  arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() );

  arrayView2d< real64 const > const & wellElemCompFrac =
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompFractionString() );
  arrayView3d< real64 const > const & dWellElemCompFrac_dCompDens =
    subRegion.getReference< array3d< real64 > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  // get well variables on perforations
  arrayView1d< real64 const > const & perfGravCoef =
    perforationData->getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );
  arrayView1d< localIndex const > const & perfWellElemIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::wellElementIndexString() );
  arrayView1d< real64 const > const & perfTrans =
    perforationData->getReference< array1d< real64 > >( PerforationData::viewKeyStruct::wellTransmissibilityString() );

  arrayView2d< real64 > const & compPerfRate =
    perforationData->getReference< array2d< real64 > >( viewKeyStruct::compPerforationRateString() );
  arrayView3d< real64 > const & dCompPerfRate_dPres =
    perforationData->getReference< array3d< real64 > >( viewKeyStruct::dCompPerforationRate_dPresString() );
  arrayView4d< real64 > const & dCompPerfRate_dComp =
    perforationData->getReference< array4d< real64 > >( viewKeyStruct::dCompPerforationRate_dCompString() );

  // get the element region, subregion, index
  arrayView1d< localIndex const > const & resElementRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementRegionString() );
  arrayView1d< localIndex const > const & resElementSubRegion =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementSubregionString() );
  arrayView1d< localIndex const > const & resElementIndex =
    perforationData->getReference< array1d< localIndex > >( PerforationData::viewKeyStruct::reservoirElementIndexString() );

  PerforationKernel::launch< parallelDevicePolicy<> >( perforationData->size(),
                                                       numFluidComponents(),
                                                       numFluidPhases(),
                                                       m_resPres.toNestedViewConst(),
                                                       m_deltaResPres.toNestedViewConst(),
                                                       m_dResPhaseVolFrac_dPres.toNestedViewConst(),
                                                       m_dResPhaseVolFrac_dCompDens.toNestedViewConst(),
                                                       m_dResCompFrac_dCompDens.toNestedViewConst(),
                                                       m_resPhaseDens.toNestedViewConst(),
                                                       m_dResPhaseDens_dPres.toNestedViewConst(),
                                                       m_dResPhaseDens_dComp.toNestedViewConst(),
                                                       m_resPhaseVisc.toNestedViewConst(),
                                                       m_dResPhaseVisc_dPres.toNestedViewConst(),
                                                       m_dResPhaseVisc_dComp.toNestedViewConst(),
                                                       m_resPhaseCompFrac.toNestedViewConst(),
                                                       m_dResPhaseCompFrac_dPres.toNestedViewConst(),
                                                       m_dResPhaseCompFrac_dComp.toNestedViewConst(),
                                                       m_resPhaseRelPerm.toNestedViewConst(),
                                                       m_dResPhaseRelPerm_dPhaseVolFrac.toNestedViewConst(),
                                                       wellElemGravCoef,
                                                       wellElemPres,
                                                       dWellElemPres,
                                                       wellElemCompDens,
                                                       dWellElemCompDens,
                                                       wellElemTotalMassDens,
                                                       dWellElemTotalMassDens_dPres,
                                                       dWellElemTotalMassDens_dCompDens,
                                                       wellElemCompFrac,
                                                       dWellElemCompFrac_dCompDens,
                                                       perfGravCoef,
                                                       perfWellElemIndex,
                                                       perfTrans,
                                                       resElementRegion,
                                                       resElementSubRegion,
                                                       resElementIndex,
                                                       compPerfRate,
                                                       dCompPerfRate_dPres,
                                                       dCompPerfRate_dComp );

}


void
CompositionalMultiphaseWell::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  // update all the fields using the global damping coefficients
  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               viewKeyStruct::deltaPressureString(),
                               scalingFactor,
                               { m_numDofPerWellElement, 0, 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               viewKeyStruct::deltaGlobalCompDensityString(),
                               scalingFactor,
                               { m_numDofPerWellElement, 1, m_numDofPerWellElement - 1 } );

  dofManager.addVectorToField( localSolution,
                               wellElementDofName(),
                               viewKeyStruct::deltaMixtureConnRateString(),
                               scalingFactor,
                               { m_numDofPerWellElement, m_numDofPerWellElement - 1, m_numDofPerWellElement } );

  // if component density chopping is allowed, some component densities may be negative after the update
  // these negative component densities are set to zero in this function
  if( m_allowCompDensChopping )
  {
    chopNegativeDensities( domain );
  }

  // synchronize
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaGlobalCompDensityString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaMixtureConnRateString() ) );
  CommunicationTools::getInstance().synchronizeFields( fieldNames,
                                                       domain.getMeshBody( 0 ).getMeshLevel( 0 ),
                                                       domain.getNeighbors(),
                                                       true );
}

void CompositionalMultiphaseWell::chopNegativeDensities( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {
    arrayView1d< integer const > const & wellElemGhostRank = subRegion.ghostRank();

    arrayView2d< real64 const > const & wellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 > const & dWellElemCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      if( wellElemGhostRank[iwelem] < 0 )
      {
        for( localIndex ic = 0; ic < NC; ++ic )
        {
          // get the latest component density (i.e., after the update of dWellElemCompDens)
          real64 const newDens = wellElemCompDens[iwelem][ic] + dWellElemCompDens[iwelem][ic];
          // we allowed for some densities to be slightly negative in CheckSystemSolution
          // if the new density is negative, chop back to zero
          if( newDens < 0 )
          {
            dWellElemCompDens[iwelem][ic] = -wellElemCompDens[iwelem][ic];
          }
        }
      }
    } );
  } );
}


void CompositionalMultiphaseWell::resetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64 > const & dWellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );
    arrayView1d< real64 > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      // extract solution and apply to dP
      dWellElemPressure[iwelem] = 0;
      dConnRate[iwelem] = 0;
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        dWellElemGlobalCompDensity[iwelem][ic] = 0;
      }
    } );
  } );

  // call constitutive models
  updateState( domain );
}


void CompositionalMultiphaseWell::resetViews( DomainPartition & domain )
{
  WellSolverBase::resetViews( domain );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  ElementRegionManager & elemManager = mesh.getElemManager();

  CompositionalMultiphaseBase & flowSolver = getParent().getGroup< CompositionalMultiphaseBase >( getFlowSolverName() );

  {
    using keys = CompositionalMultiphaseBase::viewKeyStruct;

    m_resPres.clear();
    m_resPres = elemManager.constructArrayViewAccessor< real64, 1 >( keys::pressureString() );
    m_resPres.setName( getName() + "/accessors/" + keys::pressureString() );

    m_deltaResPres.clear();
    m_deltaResPres = elemManager.constructArrayViewAccessor< real64, 1 >( keys::deltaPressureString() );
    m_deltaResPres.setName( getName() + "/accessors/" + keys::deltaPressureString() );

    m_resCompDens.clear();
    m_resCompDens = elemManager.constructArrayViewAccessor< real64, 2 >( keys::globalCompDensityString() );
    m_resCompDens.setName( getName() + "/accessors/" + keys::globalCompDensityString() );

    m_dResCompFrac_dCompDens.clear();
    m_dResCompFrac_dCompDens = elemManager.constructArrayViewAccessor< real64, 3 >( keys::dGlobalCompFraction_dGlobalCompDensityString() );
    m_dResCompFrac_dCompDens.setName( getName() + "/accessors/" + keys::dGlobalCompFraction_dGlobalCompDensityString() );

    m_resPhaseVolFrac.clear();
    m_resPhaseVolFrac = elemManager.constructArrayViewAccessor< real64, 2 >( keys::phaseVolumeFractionString() );
    m_resPhaseVolFrac.setName( getName() + "/accessors/" + keys::phaseVolumeFractionString() );

    m_dResPhaseVolFrac_dPres.clear();
    m_dResPhaseVolFrac_dPres = elemManager.constructArrayViewAccessor< real64, 2 >( keys::dPhaseVolumeFraction_dPressureString() );
    m_dResPhaseVolFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dPressureString() );

    m_dResPhaseVolFrac_dCompDens.clear();
    m_dResPhaseVolFrac_dCompDens = elemManager.constructArrayViewAccessor< real64, 3 >( keys::dPhaseVolumeFraction_dGlobalCompDensityString() );
    m_dResPhaseVolFrac_dCompDens.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dGlobalCompDensityString() );

  }
  {
    using keys = MultiFluidBase::viewKeyStruct;

    m_resPhaseDens.clear();
    m_resPhaseDens = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseDensityString(),
                                                                                  flowSolver.targetRegionNames(),
                                                                                  flowSolver.fluidModelNames() );
    m_resPhaseDens.setName( getName() + "/accessors/" + keys::phaseDensityString() );

    m_dResPhaseDens_dPres.clear();
    m_dResPhaseDens_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseDensity_dPressureString(),
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseDens_dPres.setName( getName() + "/accessors/" + keys::dPhaseDensity_dPressureString() );

    m_dResPhaseDens_dComp.clear();
    m_dResPhaseDens_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseDensity_dGlobalCompFractionString(),
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseDens_dComp.setName( getName() + "/accessors/" + keys::dPhaseDensity_dGlobalCompFractionString() );


    m_resPhaseMassDens.clear();
    m_resPhaseMassDens = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseMassDensityString(),
                                                                                      flowSolver.targetRegionNames(),
                                                                                      flowSolver.fluidModelNames() );
    m_resPhaseMassDens.setName( getName() + "/accessors/" + keys::phaseMassDensityString() );

    m_resPhaseVisc.clear();
    m_resPhaseVisc = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseViscosityString(),
                                                                                  flowSolver.targetRegionNames(),
                                                                                  flowSolver.fluidModelNames() );
    m_resPhaseVisc.setName( getName() + "/accessors/" + keys::phaseViscosityString() );

    m_dResPhaseVisc_dPres.clear();
    m_dResPhaseVisc_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::dPhaseViscosity_dPressureString(),
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dPres.setName( getName() + "/accessors/" + keys::dPhaseViscosity_dPressureString() );

    m_dResPhaseVisc_dComp.clear();
    m_dResPhaseVisc_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseViscosity_dGlobalCompFractionString(),
                                                                                         flowSolver.targetRegionNames(),
                                                                                         flowSolver.fluidModelNames() );
    m_dResPhaseVisc_dComp.setName( getName() + "/accessors/" + keys::dPhaseViscosity_dGlobalCompFractionString() );

    m_resPhaseCompFrac.clear();
    m_resPhaseCompFrac = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::phaseCompFractionString(),
                                                                                      flowSolver.targetRegionNames(),
                                                                                      flowSolver.fluidModelNames() );
    m_resPhaseCompFrac.setName( getName() + "/accessors/" + keys::phaseCompFractionString() );

    m_dResPhaseCompFrac_dPres.clear();
    m_dResPhaseCompFrac_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseCompFraction_dPressureString(),
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dPressureString() );

    m_dResPhaseCompFrac_dComp.clear();
    m_dResPhaseCompFrac_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 5 >( keys::dPhaseCompFraction_dGlobalCompFractionString(),
                                                                                             flowSolver.targetRegionNames(),
                                                                                             flowSolver.fluidModelNames() );
    m_dResPhaseCompFrac_dComp.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dGlobalCompFractionString() );

  }
  {
    using keys = RelativePermeabilityBase::viewKeyStruct;

    m_resPhaseRelPerm.clear();
    m_resPhaseRelPerm = elemManager.constructMaterialArrayViewAccessor< real64, 3 >( keys::phaseRelPermString(),
                                                                                     flowSolver.targetRegionNames(),
                                                                                     flowSolver.relPermModelNames() );
    m_resPhaseRelPerm.setName( getName() + "/accessors/" + keys::phaseRelPermString() );

    m_dResPhaseRelPerm_dPhaseVolFrac.clear();
    m_dResPhaseRelPerm_dPhaseVolFrac = elemManager.constructMaterialArrayViewAccessor< real64, 4 >( keys::dPhaseRelPerm_dPhaseVolFractionString(),
                                                                                                    flowSolver.targetRegionNames(),
                                                                                                    flowSolver.relPermModelNames() );
    m_dResPhaseRelPerm_dPhaseVolFrac.setName( getName() + "/accessors/" + keys::dPhaseRelPerm_dPhaseVolFractionString() );

  }
}


void CompositionalMultiphaseWell::formPressureRelations( DomainPartition const & domain,
                                                         DofManager const & dofManager,
                                                         CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                         arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion const & subRegion )
  {

    WellControls & wellControls = getWellControls( subRegion );

    // get the degrees of freedom, depth info, next welem index
    string const wellDofKey = dofManager.getKey( wellElementDofName() );
    arrayView1d< globalIndex const > const & wellElemDofNumber =
      subRegion.getReference< array1d< globalIndex > >( wellDofKey );
    arrayView1d< real64 const > const & wellElemGravCoef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::gravityCoefString() );
    arrayView1d< localIndex const > const & nextWellElemIndex =
      subRegion.getReference< array1d< localIndex > >( WellElementSubRegion::viewKeyStruct::nextWellElementIndexString() );

    // get primary variables on well elements
    arrayView1d< real64 const > const & wellElemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dWellElemPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

    // get total mass density on well elements (for potential calculations)
    arrayView1d< real64 const > const & wellElemTotalMassDens =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalMassDensityString() );
    arrayView1d< real64 const > const & dWellElemTotalMassDens_dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::dTotalMassDensity_dPressureString() );
    arrayView2d< real64 const > const & dWellElemTotalMassDens_dCompDens =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::dTotalMassDensity_dGlobalCompDensityString() );


    localIndex const controlHasSwitched =
      PressureRelationKernel::launch< parallelDevicePolicy<>,
                                      parallelDeviceReduce >( subRegion.size(),
                                                              dofManager.rankOffset(),
                                                              subRegion.isLocallyOwned(),
                                                              subRegion.getTopWellElementIndex(),
                                                              numFluidComponents(),
                                                              m_targetPhaseIndex,
                                                              numDofPerResElement(),
                                                              wellControls,
                                                              wellElemDofNumber,
                                                              wellElemGravCoef,
                                                              nextWellElemIndex,
                                                              wellElemPres,
                                                              dWellElemPres,
                                                              wellElemTotalMassDens,
                                                              dWellElemTotalMassDens_dPres,
                                                              dWellElemTotalMassDens_dCompDens,
                                                              localMatrix,
                                                              localRhs );

    if( controlHasSwitched == 1 )
    {
      // TODO: move the switch logic into wellControls
      // TODO: implement a more general switch when more then two constraints per well type are allowed

      if( wellControls.getControl() == WellControls::Control::BHP )
      {
        WellControls::Type const wellType = wellControls.getType();
        if( wellType == WellControls::Type::PRODUCER )
        {
          wellControls.switchToPhaseRateControl( wellControls.getTargetPhaseRate() );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from BHP constraint to phase volumetric rate constraint" );
        }
        else
        {
          wellControls.switchToTotalRateControl( wellControls.getTargetTotalRate() );
          GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                                << " from BHP constraint to total volumetric rate constraint" );
        }
      }
      else
      {
        wellControls.switchToBHPControl( wellControls.getTargetBHP() );
        GEOSX_LOG_LEVEL_RANK_0( 1, "Control switch for well " << subRegion.getName()
                                                              << " from rate constraint to BHP constraint" );
      }
    }
  } );
}


void CompositionalMultiphaseWell::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{
  MeshLevel & meshLevel = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  localIndex const NC = m_numComponents;

  forTargetSubRegions< WellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                               WellElementSubRegion & subRegion )
  {

    // get a reference to the primary variables on well elements
    arrayView1d< real64 > const & wellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const & dWellElemPressure =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

    arrayView2d< real64 > const & wellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const > const & dWellElemGlobalCompDensity =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaGlobalCompDensityString() );

    arrayView1d< real64 > const & connRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::mixtureConnRateString() );
    arrayView1d< real64 const > const & dConnRate =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaMixtureConnRateString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const iwelem )
    {
      wellElemPressure[iwelem] += dWellElemPressure[iwelem];
      for( localIndex ic = 0; ic < NC; ++ic )
      {
        wellElemGlobalCompDensity[iwelem][ic] += dWellElemGlobalCompDensity[iwelem][ic];
      }
      connRate[iwelem] += dConnRate[iwelem];
    } );
  } );
}

REGISTER_CATALOG_ENTRY( SolverBase, CompositionalMultiphaseWell, string const &, Group * const )
}// namespace geosx
