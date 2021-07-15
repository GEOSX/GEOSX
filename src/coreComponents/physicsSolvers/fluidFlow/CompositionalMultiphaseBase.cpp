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
 * @file CompositionalMultiphaseBase.cpp
 */

#include "CompositionalMultiphaseBase.hpp"

#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/capillaryPressure/capillaryPressureSelector.hpp"
#include "constitutive/fluid/multiFluidSelector.hpp"
#include "constitutive/relativePermeability/relativePermeabilitySelector.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "physicsSolvers/fluidFlow/CompositionalMultiphaseBaseKernels.hpp"

#if defined( __INTEL_COMPILER )
#pragma GCC optimize "O0"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace CompositionalMultiphaseBaseKernels;

CompositionalMultiphaseBase::CompositionalMultiphaseBase( const string & name,
                                                          Group * const parent )
  :
  FlowSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 ),
  m_computeCFLNumbers( 0 ),
  m_capPressureFlag( 0 ),
  m_maxCompFracChange( 1.0 ),
  m_minScalingFactor( 0.01 ),
  m_allowCompDensChopping( 1 )
{
//START_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::temperatureString(), &m_temperature ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Temperature" );
//END_SPHINX_INCLUDE_00
  this->registerWrapper( viewKeyStruct::useMassFlagString(), &m_useMass ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Use mass formulation instead of molar" );

  this->registerWrapper( viewKeyStruct::computeCFLNumbersString(), &m_computeCFLNumbers ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Flag indicating whether CFL numbers are computed or not" );

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

}

void CompositionalMultiphaseBase::postProcessInput()
{
  FlowSolverBase::postProcessInput();
  checkModelNames( m_relPermModelNames, viewKeyStruct::relPermNamesString() );
  m_capPressureFlag = checkModelNames( m_capPressureModelNames, viewKeyStruct::capPressureNamesString(), true );

  GEOSX_ERROR_IF_GT_MSG( m_maxCompFracChange, 1.0,
                         "The maximum absolute change in component fraction must smaller or equal to 1.0" );
  GEOSX_ERROR_IF_LT_MSG( m_maxCompFracChange, 0.0,
                         "The maximum absolute change in component fraction must larger or equal to 0.0" );
}

void CompositionalMultiphaseBase::registerDataOnMesh( Group & meshBodies )
{
  FlowSolverBase::registerDataOnMesh( meshBodies );

  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    MeshLevel & mesh = meshBody.getMeshLevel( 0 );

    forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & elementSubRegion )
    {

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString() ).setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::bcPressureString() );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array3d< real64, compflow::LAYOUT_COMP_DC > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityString() ).
        setPlotLevel( PlotLevel::LEVEL_0 );

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseMobility_dPressureString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      elementSubRegion.registerWrapper< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() ).
        setRestartFlags( RestartFlags::NO_WRITE );

      if( m_computeCFLNumbers )
      {
        elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseOutfluxString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
        elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::componentOutfluxString() ).
          setRestartFlags( RestartFlags::NO_WRITE );
        elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::phaseCFLNumberString() ).
          setPlotLevel( PlotLevel::LEVEL_0 ).
          setRestartFlags( RestartFlags::NO_WRITE );
        elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::componentCFLNumberString() )
          .setPlotLevel( PlotLevel::LEVEL_0 )
          .setRestartFlags( RestartFlags::NO_WRITE );
      }

      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionOldString() );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::totalDensityOldString() );
      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseDensityOldString() );
      elementSubRegion.registerWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityOldString() );
      elementSubRegion.registerWrapper< array3d< real64, compflow::LAYOUT_PHASE_COMP > >( viewKeyStruct::phaseComponentFractionOldString() );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityString() );
      elementSubRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString() );
    } );

    FaceManager & faceManager = mesh.getFaceManager();
    {
      faceManager.registerWrapper< array1d< real64 > >( viewKeyStruct::facePressureString() ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the pressures at the faces." );
    }

  } );
}

namespace
{

template< typename MODEL1_TYPE, typename MODEL2_TYPE >
void compareMultiphaseModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
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
void compareMulticomponentModels( MODEL1_TYPE const & lhs, MODEL2_TYPE const & rhs )
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

void CompositionalMultiphaseBase::validateConstitutiveModels( constitutive::ConstitutiveManager const & cm ) const
{
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );

  for( localIndex i = 1; i < m_fluidModelNames.size(); ++i )
  {
    MultiFluidBase const & fluid = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[i] );
    compareMultiphaseModels( fluid, fluid0 );
    compareMulticomponentModels( fluid, fluid0 );
  }

  RelativePermeabilityBase const & relPerm0 = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[0] );
  compareMultiphaseModels( relPerm0, fluid0 );

  for( localIndex i = 1; i < m_relPermModelNames.size(); ++i )
  {
    RelativePermeabilityBase const & relPerm = cm.getConstitutiveRelation< RelativePermeabilityBase >( m_relPermModelNames[i] );
    compareMultiphaseModels( relPerm, relPerm0 );
  }

  if( m_capPressureFlag )
  {
    CapillaryPressureBase const & capPres0 = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[0] );
    compareMultiphaseModels( capPres0, fluid0 );

    for( localIndex i = 1; i < m_capPressureModelNames.size(); ++i )
    {
      CapillaryPressureBase const & capPres = cm.getConstitutiveRelation< CapillaryPressureBase >( m_capPressureModelNames[i] );
      compareMultiphaseModels( capPres, capPres0 );
    }
  }
}


void CompositionalMultiphaseBase::initializePreSubGroups()
{
  FlowSolverBase::initializePreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  ConstitutiveManager const & cm = domain.getConstitutiveManager();

  // 1. Set key dimensions of the problem
  MultiFluidBase const & fluid0 = cm.getConstitutiveRelation< MultiFluidBase >( m_fluidModelNames[0] );
  m_numPhases     = fluid0.numFluidPhases();
  m_numComponents = fluid0.numFluidComponents();
  m_numDofPerCell = m_numComponents + 1;

  // 2. Validate various models against each other (must have same phases and components)
  validateConstitutiveModels( cm );

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

void CompositionalMultiphaseBase::resizeFields( MeshLevel & meshLevel ) const
{
  forTargetSubRegions( meshLevel, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );

    subRegion.getWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() ).
      setDimLabels( 1, fluid.componentNames() ).
      reference().resizeDimension< 1 >( m_numComponents );
    subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() ).
      resizeDimension< 1 >( m_numComponents );

    subRegion.getWrapper< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() ).
      setDimLabels( 1, fluid.componentNames() ).
      reference().resizeDimension< 1 >( m_numComponents );
    subRegion.getReference< array3d< real64, compflow::LAYOUT_COMP_DC > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() ).
      resizeDimension< 1, 2 >( m_numComponents, m_numComponents );

    subRegion.getWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() ).
      setDimLabels( 1, fluid.phaseNames() ).
      reference().resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() ).
      resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() ).
      resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

    subRegion.getWrapper< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityString() ).
      setDimLabels( 1, fluid.phaseNames() ).
      reference().resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseMobility_dPressureString() ).
      resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseMobility_dGlobalCompDensityString() ).
      resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

    subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionOldString() ).
      resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityOldString() ).
      resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseDensityOldString() ).
      resizeDimension< 1 >( m_numPhases );
    subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_COMP > >( viewKeyStruct::phaseComponentFractionOldString() ).
      resizeDimension< 1, 2 >( m_numPhases, m_numComponents );

    if( m_computeCFLNumbers )
    {
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseOutfluxString() ).
        resizeDimension< 1 >( m_numPhases );
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::componentOutfluxString() ).
        resizeDimension< 1 >( m_numComponents );
    }
  } );
}

void CompositionalMultiphaseBase::updateComponentFraction( Group & dataGroup ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64, compflow::USD_COMP > const & compFrac =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() );

  arrayView3d< real64, compflow::USD_COMP_DC > const & dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64, compflow::LAYOUT_COMP_DC > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  // inputs

  arrayView2d< real64 const, compflow::USD_COMP > const compDens =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );

  arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );

  KernelLaunchSelector1< ComponentFractionKernel >( m_numComponents,
                                                    dataGroup.size(),
                                                    compDens,
                                                    dCompDens,
                                                    compFrac,
                                                    dCompFrac_dCompDens );
}

void CompositionalMultiphaseBase::updatePhaseVolumeFraction( Group & dataGroup,
                                                             localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // outputs

  arrayView2d< real64, compflow::USD_PHASE > const phaseVolFrac =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );

  arrayView2d< real64, compflow::USD_PHASE > const dPhaseVolFrac_dPres =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );

  arrayView3d< real64, compflow::USD_PHASE_DC > const dPhaseVolFrac_dComp =
    dataGroup.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

  // inputs

  arrayView3d< real64 const, compflow::USD_COMP_DC > const dCompFrac_dCompDens =
    dataGroup.getReference< array3d< real64, compflow::LAYOUT_COMP_DC > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

  arrayView2d< real64 const, compflow::USD_COMP > const compDens =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );

  arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );

  MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseFrac = fluid.phaseFraction();
  arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseFrac_dPres = fluid.dPhaseFraction_dPressure();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseFrac_dComp = fluid.dPhaseFraction_dGlobalCompFraction();

  arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
  arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
  arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();

  KernelLaunchSelector2< PhaseVolumeFractionKernel >( m_numComponents, m_numPhases,
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

void CompositionalMultiphaseBase::updateFluidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const pres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
  arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() );

  MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( dataGroup, m_fluidModelNames[targetIndex] );

  constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
  {
    using FluidType = TYPEOFREF( castedFluid );
    using ExecPolicy = typename FluidType::exec_policy;
    typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

    FluidUpdateKernel::launch< ExecPolicy >( dataGroup.size(),
                                             fluidWrapper,
                                             pres,
                                             dPres,
                                             m_temperature,
                                             compFrac );
  } );
}

void CompositionalMultiphaseBase::updateSolidModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  ConstitutiveBase & solid = getConstitutiveModel< ConstitutiveBase >( dataGroup, m_solidModelNames[targetIndex] );

  arrayView1d< real64 const > const pres  = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
  arrayView1d< real64 const > const dPres = dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );

  solid.stateUpdateBatchPressure( pres, dPres );
}

void CompositionalMultiphaseBase::updateRelPermModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
    dataGroup.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );

  RelativePermeabilityBase & relPerm =
    getConstitutiveModel< RelativePermeabilityBase >( dataGroup, m_relPermModelNames[targetIndex] );

  constitutive::constitutiveUpdatePassThru( relPerm, [&] ( auto & castedRelPerm )
  {
    typename TYPEOFREF( castedRelPerm ) ::KernelWrapper relPermWrapper = castedRelPerm.createKernelWrapper();

    RelativePermeabilityUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                        relPermWrapper,
                                                                        phaseVolFrac );
  } );
}

void CompositionalMultiphaseBase::updateCapPressureModel( Group & dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  if( m_capPressureFlag )
  {
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      dataGroup.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );

    CapillaryPressureBase & capPressure =
      getConstitutiveModel< CapillaryPressureBase >( dataGroup, m_capPressureModelNames[targetIndex] );

    constitutive::constitutiveUpdatePassThru( capPressure, [&] ( auto & castedCapPres )
    {
      typename TYPEOFREF( castedCapPres ) ::KernelWrapper capPresWrapper = castedCapPres.createKernelWrapper();

      CapillaryPressureUpdateKernel::launch< parallelDevicePolicy<> >( dataGroup.size(),
                                                                       capPresWrapper,
                                                                       phaseVolFrac );
    } );
  }
}

void CompositionalMultiphaseBase::updateState( Group & dataGroup, localIndex const targetIndex ) const
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

void CompositionalMultiphaseBase::initializeFluidState( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    // 1. Assume global component fractions have been prescribed.
    // Initialize constitutive state to get fluid density.
    updateFluidModel( subRegion, targetIndex );

    // 2. Back-calculate global component densities from fractions and total fluid density
    // in order to initialize the primary solution variables
    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() );
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex ic = 0; ic < numComp; ++ic )
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

void CompositionalMultiphaseBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::initializePostInitialConditionsPreSubGroups();

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::pressureString() ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::globalCompDensityString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), false );

  // set mass fraction flag on fluid models
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    MultiFluidBase & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, m_fluidModelNames[targetIndex] );
    fluid.setMassFlag( m_useMass );
  } );

  // Initialize primary variables from applied initial conditions
  initializeFluidState( mesh );


  // TODO this can be removed once we fully move to porosity model
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  ElementSubRegionBase & subRegion )
  {
    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );
    arrayView1d< real64 const > const poroRef = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
    arrayView1d< real64 > const poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString() );

    bool poroInit = false;
    if( solid.hasWrapper( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() ) )
    {
      arrayView2d< real64 const > const
      pvmult = solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );
      if( pvmult.size() == poro.size() )
      {
        forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
        {
          poro[ei] = poroRef[ei] * pvmult[ei][0];
        } );
        poroInit = true;
      }
    }
    if( !poroInit )
    {
      poro.setValues< parallelDevicePolicy<> >( poroRef );
    }
  } );
}

real64 CompositionalMultiphaseBase::solverStep( real64 const & time_n,
                                                real64 const & dt,
                                                integer const cycleNumber,
                                                DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  // Only build the sparsity pattern once
  // TODO: this should be triggered by a topology change indicator
  static bool systemSetupDone = false;
  if( !systemSetupDone )
  {
    setupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );
    systemSetupDone = true;
  }

  implicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dt_return = nonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  implicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseBase::backupFields( MeshLevel & mesh ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;
  integer const numPhase = m_numPhases;

  // backup some fields used in time derivative approximation
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const poroRef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const phaseVolFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseMob =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityString() );

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const phaseDens = fluid.phaseDensity();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const phaseCompFrac = fluid.phaseCompFraction();

    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
    arrayView1d< real64 const > const & poro = subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString() );

    arrayView1d< real64 > const totalDensOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::totalDensityOldString() );
    arrayView2d< real64, compflow::USD_PHASE > const phaseDensOld =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseDensityOldString() );
    arrayView2d< real64, compflow::USD_PHASE > const phaseVolFracOld =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionOldString() );
    arrayView2d< real64, compflow::USD_PHASE > const phaseMobOld =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseMobilityOldString() );
    arrayView3d< real64, compflow::USD_PHASE_COMP > const phaseCompFracOld =
      subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_COMP > >( viewKeyStruct::phaseComponentFractionOldString() );
    arrayView1d< real64 > const poroOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( elemGhostRank[ei] >= 0 )
        return;

      for( localIndex ip = 0; ip < numPhase; ++ip )
      {
        phaseDensOld[ei][ip] = phaseDens[ei][0][ip];
        phaseVolFracOld[ei][ip] = phaseVolFrac[ei][ip];
        phaseMobOld[ei][ip] = phaseMob[ei][ip];

        for( localIndex ic = 0; ic < numComp; ++ic )
        {
          phaseCompFracOld[ei][ip][ic] = phaseCompFrac[ei][0][ip][ic];
        }
      }
      totalDensOld[ei] = totalDens[ei][0];
      poroOld[ei] = poro[ei];
    } );
  } );
}

void
CompositionalMultiphaseBase::implicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                                real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  // bind the stored views to the current domain
  static bool viewsSet = false;
  if( !viewsSet )
  {
    resetViews( mesh );
    viewsSet = true;
  }

  // set deltas to zero and recompute dependent quantities
  resetStateToBeginningOfStep( domain );

  // backup fields used in time derivative approximation
  backupFields( mesh );
}

void CompositionalMultiphaseBase::assembleSystem( real64 const GEOSX_UNUSED_PARAM( time_n ),
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

void CompositionalMultiphaseBase::assembleAccumulationTerms( DomainPartition & domain,
                                                             DofManager const & dofManager,
                                                             CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                             arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & porosityRef =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString() );

    arrayView1d< real64 > const porosity =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString() );

    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
    arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );
    arrayView3d< real64 const, compflow::USD_COMP_DC > const & dCompFrac_dCompDens =
      subRegion.getReference< array3d< real64, compflow::LAYOUT_COMP_DC > >( viewKeyStruct::dGlobalCompFraction_dGlobalCompDensityString() );

    arrayView1d< real64 const > const & porosityOld =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFracOld =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionOldString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseDensOld =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseDensityOldString() );
    arrayView3d< real64 const, compflow::USD_PHASE_COMP > const & phaseCompFracOld =
      subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_COMP > >( viewKeyStruct::phaseComponentFractionOldString() );

    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
    arrayView2d< real64 const > const & pvMult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );
    arrayView2d< real64 const > const & dPvMult_dPres =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString() );

    MultiFluidBase const & fluid = getConstitutiveModel< MultiFluidBase >( subRegion, fluidModelNames()[targetIndex] );
    arrayView3d< real64 const, multifluid::USD_PHASE > const & phaseDens = fluid.phaseDensity();
    arrayView3d< real64 const, multifluid::USD_PHASE > const & dPhaseDens_dPres = fluid.dPhaseDensity_dPressure();
    arrayView4d< real64 const, multifluid::USD_PHASE_DC > const & dPhaseDens_dComp = fluid.dPhaseDensity_dGlobalCompFraction();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & phaseCompFrac = fluid.phaseCompFraction();
    arrayView4d< real64 const, multifluid::USD_PHASE_COMP > const & dPhaseCompFrac_dPres = fluid.dPhaseCompFraction_dPressure();
    arrayView5d< real64 const, multifluid::USD_PHASE_COMP_DC > const & dPhaseCompFrac_dComp = fluid.dPhaseCompFraction_dGlobalCompFraction();

    KernelLaunchSelector1< AccumulationKernel >( m_numComponents,
                                                 m_numPhases,
                                                 subRegion.size(),
                                                 dofManager.rankOffset(),
                                                 dofNumber,
                                                 elemGhostRank,
                                                 volume,
                                                 porosity,
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

void CompositionalMultiphaseBase::assembleVolumeBalanceTerms( DomainPartition const & domain,
                                                              DofManager const & dofManager,
                                                              CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                              arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();

    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
    arrayView1d< real64 const > const & porosityRef =
      subRegion.getReference< array1d< real64 > >( FlowSolverBase::viewKeyStruct::referencePorosityString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & phaseVolFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::phaseVolumeFractionString() );
    arrayView2d< real64 const, compflow::USD_PHASE > const & dPhaseVolFrac_dPres =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_PHASE > >( viewKeyStruct::dPhaseVolumeFraction_dPressureString() );
    arrayView3d< real64 const, compflow::USD_PHASE_DC > const & dPhaseVolFrac_dCompDens =
      subRegion.getReference< array3d< real64, compflow::LAYOUT_PHASE_DC > >( viewKeyStruct::dPhaseVolumeFraction_dGlobalCompDensityString() );

    ConstitutiveBase const & solid = getConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );
    arrayView2d< real64 const > const & pvMult =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString() );
    arrayView2d< real64 const > const & dPvMult_dPres =
      solid.getReference< array2d< real64 > >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString() );

    KernelLaunchSelector2< VolumeBalanceKernel >( m_numComponents, m_numPhases,
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

void CompositionalMultiphaseBase::applyBoundaryConditions( real64 const time_n,
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

void CompositionalMultiphaseBase::applySourceFluxBC( real64 const time,
                                                     real64 const dt,
                                                     DofManager const & dofManager,
                                                     DomainPartition & domain,
                                                     CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                     arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

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

    // TODO: can we avoid having to rebuild the set and moving host->device?
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


void CompositionalMultiphaseBase::applyDirichletBC( real64 const time,
                                                    real64 const dt,
                                                    DofManager const & dofManager,
                                                    DomainPartition & domain,
                                                    CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                    arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  integer const numComp = m_numComponents;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

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
    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
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
    integer const comp = fs.getComponent();
    GEOSX_ERROR_IF( bcStatusMap[subRegionName].count( setName ) == 0,
                    "Pressure boundary condition not prescribed on set '" << setName << "'" );
    GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp],
                    "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
    bcStatusMap[subRegionName][setName][comp] = true;

    // 2.1. Apply BC to set the field values
    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
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
  string const dofKey = dofManager.getKey( viewKeyStruct::elemDofFieldString() );

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

    arrayView1d< real64 const > const bcPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::bcPressureString() );
    arrayView2d< real64 const, compflow::USD_COMP > const compFrac =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompFractionString() );

    constitutiveUpdatePassThru( fluid, [&] ( auto & castedFluid )
    {
      using FluidType = TYPEOFREF( castedFluid );
      using ExecPolicy = typename FluidType::exec_policy;
      typename FluidType::KernelWrapper fluidWrapper = castedFluid.createKernelWrapper();

      FluidUpdateKernel::launch< ExecPolicy >( targetSet,
                                               fluidWrapper,
                                               bcPres,
                                               m_temperature,
                                               compFrac );
    } );

    arrayView1d< integer const > const ghostRank =
      subRegion.getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString() );
    arrayView1d< globalIndex const > const dofNumber =
      subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< real64 const > const pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );
    arrayView2d< real64 const, multifluid::USD_FLUID > const totalDens = fluid.totalDensity();

    forAll< parallelDevicePolicy<> >( targetSet.size(), [=] GEOSX_HOST_DEVICE ( localIndex const a )
    {
      localIndex const ei = targetSet[a];
      if( ghostRank[ei] >= 0 )
      {
        return;
      }

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
      for( localIndex ic = 0; ic < numComp; ++ic )
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

void CompositionalMultiphaseBase::solveSystem( DofManager const & dofManager,
                                               ParallelMatrix & matrix,
                                               ParallelVector & rhs,
                                               ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::solveSystem( dofManager, matrix, rhs, solution );
}

void CompositionalMultiphaseBase::chopNegativeDensities( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  integer const numComp = m_numComponents;
  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< integer const > const ghostRank = subRegion.ghostRank();

    arrayView2d< real64 const, compflow::USD_COMP > const compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );
    arrayView2d< real64, compflow::USD_COMP > const dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      if( ghostRank[ei] < 0 )
      {
        for( localIndex ic = 0; ic < numComp; ++ic )
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

void CompositionalMultiphaseBase::resetStateToBeginningOfStep( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64, compflow::USD_COMP > const & dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );

    dPres.zero();
    dCompDens.zero();

    updateState( subRegion, targetIndex );
  } );
}

void CompositionalMultiphaseBase::implicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                                                        DomainPartition & domain )
{
  integer const numComp = m_numComponents;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView1d< real64 const > const dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString() );
    arrayView2d< real64 const, compflow::USD_COMP > const dCompDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::deltaGlobalCompDensityString() );

    arrayView1d< real64 > const pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString() );
    arrayView2d< real64, compflow::USD_COMP > const compDens =
      subRegion.getReference< array2d< real64, compflow::LAYOUT_COMP > >( viewKeyStruct::globalCompDensityString() );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      for( localIndex ic = 0; ic < numComp; ++ic )
      {
        compDens[ei][ic] += dCompDens[ei][ic];
      }
    } );
  } );
}

void CompositionalMultiphaseBase::resetViews( MeshLevel & mesh )
{
  FlowSolverBase::resetViews( mesh );
  ElementRegionManager const & elemManager = mesh.getElemManager();

  {
    using keys = viewKeyStruct;
    using namespace compflow;

    m_pressure.clear();
    m_pressure = elemManager.constructArrayViewAccessor< real64, 1 >( keys::pressureString() );
    m_pressure.setName( getName() + "/accessors/" + keys::pressureString() );

    m_deltaPressure.clear();
    m_deltaPressure = elemManager.constructArrayViewAccessor< real64, 1 >( keys::deltaPressureString() );
    m_deltaPressure.setName( getName() + "/accessors/" + keys::deltaPressureString() );

    m_dCompFrac_dCompDens.clear();
    m_dCompFrac_dCompDens =
      elemManager.constructArrayViewAccessor< real64, 3, LAYOUT_COMP_DC >( keys::dGlobalCompFraction_dGlobalCompDensityString() );
    m_dCompFrac_dCompDens.setName( getName() + "/accessors/" + keys::dGlobalCompFraction_dGlobalCompDensityString() );

    m_dPhaseVolFrac_dPres.clear();
    m_dPhaseVolFrac_dPres =
      elemManager.constructArrayViewAccessor< real64, 2, LAYOUT_PHASE >( keys::dPhaseVolumeFraction_dPressureString() );
    m_dPhaseVolFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dPressureString() );

    m_dPhaseVolFrac_dCompDens.clear();
    m_dPhaseVolFrac_dCompDens =
      elemManager.constructArrayViewAccessor< real64, 3, LAYOUT_PHASE_DC >( keys::dPhaseVolumeFraction_dGlobalCompDensityString() );
    m_dPhaseVolFrac_dCompDens.setName( getName() + "/accessors/" + keys::dPhaseVolumeFraction_dGlobalCompDensityString() );

    m_phaseMob.clear();
    m_phaseMob =
      elemManager.constructArrayViewAccessor< real64, 2, LAYOUT_PHASE >( keys::phaseMobilityString() );
    m_phaseMob.setName( getName() + "/accessors/" + keys::phaseMobilityString() );

    m_dPhaseMob_dPres.clear();
    m_dPhaseMob_dPres =
      elemManager.constructArrayViewAccessor< real64, 2, LAYOUT_PHASE >( keys::dPhaseMobility_dPressureString() );
    m_dPhaseMob_dPres.setName( getName() + "/accessors/" + keys::dPhaseMobility_dPressureString() );

    m_dPhaseMob_dCompDens.clear();
    m_dPhaseMob_dCompDens =
      elemManager.constructArrayViewAccessor< real64, 3, LAYOUT_PHASE_DC >( keys::dPhaseMobility_dGlobalCompDensityString() );
    m_dPhaseMob_dCompDens.setName( getName() + "/accessors/" + keys::dPhaseMobility_dGlobalCompDensityString() );

    m_phaseMobOld.clear();
    m_phaseMobOld =
      elemManager.constructArrayViewAccessor< real64, 2, LAYOUT_PHASE >( keys::phaseMobilityOldString() );
    m_phaseMobOld.setName( getName() + "/accessors/" + keys::phaseMobilityOldString() );

    m_totalDensOld.clear();
    m_totalDensOld = elemManager.constructArrayViewAccessor< real64, 1 >( keys::totalDensityOldString() );
    m_totalDensOld.setName( getName() + "/accessors/" + keys::totalDensityOldString() );
  }

  {
    using keys = MultiFluidBase::viewKeyStruct;
    using namespace constitutive::multifluid;

    m_phaseVisc.clear();
    m_phaseVisc = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_PHASE >( keys::phaseViscosityString(),
                                                                                             targetRegionNames(),
                                                                                             fluidModelNames() );
    m_phaseVisc.setName( getName() + "/accessors/" + keys::phaseViscosityString() );
    m_phaseDens.clear();
    m_phaseDens = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_PHASE >( keys::phaseDensityString(),
                                                                                             targetRegionNames(),
                                                                                             fluidModelNames() );
    m_phaseDens.setName( getName() + "/accessors/" + keys::phaseDensityString() );

    m_dPhaseDens_dPres.clear();
    m_dPhaseDens_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_PHASE >( keys::dPhaseDensity_dPressureString(),
                                                                                                    targetRegionNames(),
                                                                                                    fluidModelNames() );
    m_dPhaseDens_dPres.setName( getName() + "/accessors/" + keys::dPhaseDensity_dPressureString() );

    m_dPhaseDens_dComp.clear();
    m_dPhaseDens_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4, LAYOUT_PHASE_DC >( keys::dPhaseDensity_dGlobalCompFractionString(),
                                                                                                       targetRegionNames(),
                                                                                                       fluidModelNames() );
    m_dPhaseDens_dComp.setName( getName() + "/accessors/" + keys::dPhaseDensity_dGlobalCompFractionString() );

    m_phaseMassDens.clear();
    m_phaseMassDens = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_PHASE >( keys::phaseMassDensityString(),
                                                                                                 targetRegionNames(),
                                                                                                 fluidModelNames() );
    m_phaseMassDens.setName( getName() + "/accessors/" + keys::phaseMassDensityString() );

    m_dPhaseMassDens_dPres.clear();
    m_dPhaseMassDens_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_PHASE >( keys::dPhaseMassDensity_dPressureString(),
                                                                                                        targetRegionNames(),
                                                                                                        fluidModelNames() );
    m_dPhaseMassDens_dPres.setName( getName() + "/accessors/" + keys::dPhaseMassDensity_dPressureString() );

    m_dPhaseMassDens_dComp.clear();
    m_dPhaseMassDens_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 4, LAYOUT_PHASE_DC >( keys::dPhaseMassDensity_dGlobalCompFractionString(),
                                                                                                           targetRegionNames(),
                                                                                                           fluidModelNames() );
    m_dPhaseMassDens_dComp.setName( getName() + "/accessors/" + keys::dPhaseMassDensity_dGlobalCompFractionString() );

    m_phaseCompFrac.clear();
    m_phaseCompFrac = elemManager.constructMaterialArrayViewAccessor< real64, 4, LAYOUT_PHASE_COMP >( keys::phaseCompFractionString(),
                                                                                                      targetRegionNames(),
                                                                                                      fluidModelNames() );
    m_phaseCompFrac.setName( getName() + "/accessors/" + keys::phaseCompFractionString() );

    m_dPhaseCompFrac_dPres.clear();
    m_dPhaseCompFrac_dPres = elemManager.constructMaterialArrayViewAccessor< real64, 4, LAYOUT_PHASE_COMP >( keys::dPhaseCompFraction_dPressureString(),
                                                                                                             targetRegionNames(),
                                                                                                             fluidModelNames() );
    m_dPhaseCompFrac_dPres.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dPressureString() );

    m_dPhaseCompFrac_dComp.clear();
    m_dPhaseCompFrac_dComp = elemManager.constructMaterialArrayViewAccessor< real64, 5, LAYOUT_PHASE_COMP_DC >( keys::dPhaseCompFraction_dGlobalCompFractionString(),
                                                                                                                targetRegionNames(),
                                                                                                                fluidModelNames() );
    m_dPhaseCompFrac_dComp.setName( getName() + "/accessors/" + keys::dPhaseCompFraction_dGlobalCompFractionString() );
  }
  {
    using keys = RelativePermeabilityBase::viewKeyStruct;
    using namespace constitutive::relperm;

    m_phaseRelPerm.clear();
    m_phaseRelPerm = elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_RELPERM >( keys::phaseRelPermString(),
                                                                                                  targetRegionNames(),
                                                                                                  relPermModelNames() );
    m_phaseRelPerm.setName( getName() + "/accessors/" + keys::phaseRelPermString() );
  }
  if( m_capPressureFlag )
  {
    using keys = CapillaryPressureBase::viewKeyStruct;
    using namespace constitutive::cappres;

    m_phaseCapPressure.clear();
    m_phaseCapPressure =
      elemManager.constructMaterialArrayViewAccessor< real64, 3, LAYOUT_CAPPRES >( keys::phaseCapPressureString(),
                                                                                   targetRegionNames(),
                                                                                   capPresModelNames() );
    m_phaseCapPressure.setName( getName() + "/accessors/" + keys::phaseCapPressureString() );

    m_dPhaseCapPressure_dPhaseVolFrac.clear();
    m_dPhaseCapPressure_dPhaseVolFrac =
      elemManager.constructMaterialArrayViewAccessor< real64, 4, LAYOUT_CAPPRES_DS >( keys::dPhaseCapPressure_dPhaseVolFractionString(),
                                                                                      targetRegionNames(),
                                                                                      capPresModelNames() );
    m_dPhaseCapPressure_dPhaseVolFrac.setName( getName() + "/accessors/" + keys::dPhaseCapPressure_dPhaseVolFractionString() );
  }
}

}// namespace geosx
