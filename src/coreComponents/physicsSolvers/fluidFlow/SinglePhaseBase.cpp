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
 * @file SinglePhaseBase.cpp
 */

#include "SinglePhaseBase.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"

namespace geosx
{
using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;

SinglePhaseBase::SinglePhaseBase( const std::string & name, Group * const parent ) :
  FlowSolverBase( name, parent )
{
  m_numDofPerCell = 1;
}

void
SinglePhaseBase::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel * meshLevel =
      Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forElementSubRegions< CellElementSubRegion >(
      [&]( CellElementSubRegion & subRegion ) {
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::deltaVolumeString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mobilityString );

        subRegion
          .registerWrapper< array1d< real64 > >(
            viewKeyStruct::dMobility_dPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::porosityString )
          ->setPlotLevel( PlotLevel::LEVEL_1 );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )
          ->setRestartFlags( RestartFlags::NO_WRITE );
      } );

    elemManager->forElementSubRegions< FaceElementSubRegion >(
      [&]( FaceElementSubRegion & subRegion ) {
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )
          ->setPlotLevel( PlotLevel::LEVEL_0 );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::deltaVolumeString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::mobilityString );

        subRegion
          .registerWrapper< array1d< real64 > >(
            viewKeyStruct::dMobility_dPressureString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::porosityString )
          ->setDefaultValue( 1.0 );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::porosityOldString )
          ->setDefaultValue( 1.0 )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )
          ->setRestartFlags( RestartFlags::NO_WRITE );

        subRegion
          .registerWrapper< array1d< R1Tensor > >( viewKeyStruct::transTMultString )
          ->setDefaultValue( { 1.0, 1.0, 1.0 } );
        subRegion
          .registerWrapper< array1d< real64 > >( viewKeyStruct::poroMultString )
          ->setDefaultValue( 1.0 );
      } );

    FaceManager * const faceManager = meshLevel->getFaceManager();
    {
      faceManager
        ->registerWrapper< array1d< real64 > >( viewKeyStruct::facePressureString )
        ->setPlotLevel( PlotLevel::LEVEL_0 )
        ->setRegisteringObjects( this->getName() )
        ->setDescription( "An array that holds the pressures at the faces." );
    }
  }
}

void
SinglePhaseBase::ValidateFluidModels( DomainPartition const & domain ) const
{
  for( auto & mesh : domain.getMeshBodies()->GetSubGroups() )
  {
    MeshLevel const & meshLevel =
      *Group::group_cast< MeshBody const * >( mesh.second )->getMeshLevel( 0 );
    ValidateModelMapping< SingleFluidBase >( *meshLevel.getElemManager(),
                                             m_fluidModelNames );
  }
}

SinglePhaseBase::FluidPropViews
SinglePhaseBase::getFluidProperties(
  ConstitutiveBase const & fluid ) const
{
  SingleFluidBase const & singleFluid =
    dynamicCast< SingleFluidBase const & >( fluid );
  return { singleFluid.density(),
           singleFluid.dDensity_dPressure(),
           singleFluid.viscosity(),
           singleFluid.dViscosity_dPressure(),
           singleFluid.defaultDensity(),
           singleFluid.defaultViscosity() };
}

arrayView1d< real64 const > const &
SinglePhaseBase::getPoreVolumeMult(
  ElementSubRegionBase const & subRegion ) const
{
  return subRegion.getReference< array1d< real64 > >( viewKeyStruct::poroMultString );
}

void
SinglePhaseBase::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  ValidateFluidModels( domain );
}

void
SinglePhaseBase::UpdateFluidModel( Group & dataGroup,
                                   localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  SingleFluidBase & fluid =
    GetConstitutiveModel< SingleFluidBase >( dataGroup,
                                             m_fluidModelNames[targetIndex] );

  constitutiveUpdatePassThru( fluid, [&]( auto & castedFluid ) {
    typename TYPEOFREF( castedFluid )::KernelWrapper fluidWrapper =
      castedFluid.createKernelWrapper();
    FluidUpdateKernel::Launch( fluidWrapper, pres, dPres );
  } );
}

void
SinglePhaseBase::UpdateSolidModel( Group & dataGroup,
                                   localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  arrayView1d< real64 const > const & pres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  ConstitutiveBase & solid =
    GetConstitutiveModel< ConstitutiveBase >( dataGroup,
                                              m_solidModelNames[targetIndex] );
  solid.StateUpdateBatchPressure( pres, dPres );
}

void
SinglePhaseBase::UpdateMobility( Group & dataGroup,
                                 localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  // output

  arrayView1d< real64 > const & mob =
    dataGroup.getReference< array1d< real64 > >( viewKeyStruct::mobilityString );

  arrayView1d< real64 > const & dMob_dPres = dataGroup.getReference< array1d< real64 > >(
    viewKeyStruct::dMobility_dPressureString );

  // input

  ConstitutiveBase & fluid =
    GetConstitutiveModel( dataGroup, m_fluidModelNames[targetIndex] );
  FluidPropViews fluidProps = getFluidProperties( fluid );

  SinglePhaseBaseKernels::MobilityKernel::Launch< parallelDevicePolicy<> >(
    dataGroup.size(),
    fluidProps.dens,
    fluidProps.dDens_dPres,
    fluidProps.visc,
    fluidProps.dVisc_dPres,
    mob,
    dMob_dPres );
}

void
SinglePhaseBase::UpdateState( Group & dataGroup,
                              localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateFluidModel( dataGroup, targetIndex );
  UpdateSolidModel( dataGroup, targetIndex );
  UpdateMobility( dataGroup, targetIndex );
}

void
SinglePhaseBase::InitializePostInitialConditions_PreSubGroups(
  Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::pressureString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain->getNeighbors() );

  ResetViews( mesh );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in ApplySystemSolution and ImplicitStepComplete, respectively
  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      ConstitutiveBase const & fluid =
        GetConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );

      real64 const defaultDensity = getFluidProperties( fluid ).defaultDensity;
      subRegion.getWrapper< array1d< real64 > >( viewKeyStruct::densityOldString )
        ->setDefaultValue( defaultDensity );

      UpdateState( subRegion, targetIndex );

      ConstitutiveBase const & solid =
        GetConstitutiveModel( subRegion, m_solidModelNames[targetIndex] );
      arrayView1d< real64 const > const & poroRef =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::referencePorosityString );
      arrayView2d< real64 const > const & pvmult =
        solid.getReference< array2d< real64 > >(
          ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );

      arrayView1d< real64 > const & poro =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );

      if( pvmult.size() == poro.size() )
      {
        forAll< parallelDevicePolicy<> >( subRegion.size(),
                                          [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                            poro[ei] = poroRef[ei] * pvmult[ei][0];
                                          } );
      }
      else
      {
        poro.setValues< parallelDevicePolicy<> >( poroRef );
      }
    } );

  mesh.getElemManager()->forElementRegions< FaceElementRegion >(
    targetRegionNames(),
    [&]( localIndex const targetIndex, FaceElementRegion & region ) {
      region.forElementSubRegions< FaceElementSubRegion >(
        [&]( FaceElementSubRegion & subRegion ) {
          ConstitutiveBase & fluid =
            GetConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
          real64 const defaultDensity = getFluidProperties( fluid ).defaultDensity;

          subRegion
            .getWrapper< real64_array >( viewKeyStruct::effectiveApertureString )
            ->setApplyDefaultValue( region.getDefaultAperture() );

          subRegion
            .getWrapper< real64_array >(
              FaceElementSubRegion::viewKeyStruct::creationMassString )
            ->setApplyDefaultValue( defaultDensity * region.getDefaultAperture() );
        } );
    } );

  BackupFields( mesh );
}

real64
SinglePhaseBase::SolverStep( real64 const & time_n,
                             real64 const & dt,
                             const int cycleNumber,
                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return;

  // setup dof numbers and linear system
  SetupSystem( domain, m_dofManager, m_localMatrix, m_localRhs, m_localSolution );

  ImplicitStepSetup( time_n, dt, domain );

  // currently the only method is implicit time integration
  dt_return = NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void
SinglePhaseBase::SetupSystem( DomainPartition & domain,
                              DofManager & dofManager,
                              CRSMatrix< real64, globalIndex > & localMatrix,
                              array1d< real64 > & localRhs,
                              array1d< real64 > & localSolution,
                              bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  ResetViews( *( domain.getMeshBody( 0 )->getMeshLevel( 0 ) ) );

  SolverBase::SetupSystem( domain,
                           dofManager,
                           localMatrix,
                           localRhs,
                           localSolution,
                           setSparsity );
}

void
SinglePhaseBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                    real64 const & GEOSX_UNUSED_PARAM( dt ),
                                    DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  ResetViews( mesh );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      arrayView1d< real64 > const & dPres = subRegion.getReference< array1d< real64 > >(
        viewKeyStruct::deltaPressureString );
      arrayView1d< real64 > const & dVol =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );

      dPres.setValues< parallelDevicePolicy<> >( 0.0 );
      dVol.setValues< parallelDevicePolicy<> >( 0.0 );

      // This should fix NaN density in newly created fracture elements
      UpdateState( subRegion, targetIndex );
    } );

  forTargetSubRegions< FaceElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, FaceElementSubRegion & subRegion ) {
      arrayView1d< real64 const > const & aper =
        subRegion.getReference< array1d< real64 > >(
          viewKeyStruct::effectiveApertureString );
      arrayView1d< real64 > const & aper0 =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::aperture0String );

      aper0.setValues< parallelDevicePolicy<> >( aper );

      // UpdateMobility( &subRegion );
      UpdateState( subRegion, targetIndex );
    } );

  BackupFields( mesh );
}

void
SinglePhaseBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion ) {
    arrayView1d< real64 const > const & dPres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
    arrayView1d< real64 const > const & dVol =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );

    arrayView1d< real64 > const & pres =
      subRegion.getReference< array1d< real64 > >( viewKeyStruct::pressureString );
    arrayView1d< real64 > const & vol = subRegion.getReference< array1d< real64 > >(
      CellBlock::viewKeyStruct::elementVolumeString );

    forAll< parallelDevicePolicy<> >( subRegion.size(),
                                      [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                        pres[ei] += dPres[ei];
                                        vol[ei] += dVol[ei];
                                      } );
  } );

  forTargetSubRegions< FaceElementSubRegion >(
    mesh,
    [&]( localIndex const, FaceElementSubRegion & subRegion ) {
      arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
      arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
      arrayView1d< real64 const > const & densOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
      arrayView1d< real64 > const & creationMass =
        subRegion.getReference< real64_array >(
          FaceElementSubRegion::viewKeyStruct::creationMassString );

      forAll< parallelDevicePolicy<> >(
        subRegion.size(),
        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
          if( elemGhostRank[ei] < 0 )
          {
            if( volume[ei] * densOld[ei] > 1.1 * creationMass[ei] )
            {
              creationMass[ei] *= 0.75;
              if( creationMass[ei] < 1.0e-20 )
              {
                creationMass[ei] = 0.0;
              }
            }
          }
        } );
    } );
}

void
SinglePhaseBase::AssembleSystem(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  if( m_poroElasticFlag )
  {
    AssembleAccumulationTerms< true, parallelDevicePolicy<> >( domain,
                                                               dofManager,
                                                               localMatrix,
                                                               localRhs );
  }
  else
  {
    AssembleAccumulationTerms< false, parallelDevicePolicy<> >( domain,
                                                                dofManager,
                                                                localMatrix,
                                                                localRhs );
  }

  AssembleFluxTerms( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

template< bool ISPORO, typename POLICY >
void
SinglePhaseBase::AccumulationLaunch(
  localIndex const targetIndex,
  CellElementSubRegion & subRegion,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dofNumber =
    subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const & deltaPressure =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );
  arrayView1d< real64 const > const & densityOld =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
  arrayView1d< real64 > const & porosity =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
  arrayView1d< real64 const > const & porosityOld =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );
  arrayView1d< real64 const > const & porosityRef =
    subRegion.getReference< array1d< real64 > >(
      viewKeyStruct::referencePorosityString );
  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const & deltaVolume =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );

  arrayView1d< real64 const > const & totalMeanStressOld = ISPORO
    ? subRegion.getReference< array1d< real64 > >( "oldTotalMeanStress" )
    : porosityOld;
  arrayView1d< real64 const > const & totalMeanStress = ISPORO
    ? subRegion.getReference< array1d< real64 > >( "totalMeanStress" )
    : porosityOld;

  ConstitutiveBase const & fluid =
    GetConstitutiveModel( subRegion, fluidModelNames()[targetIndex] );
  FluidPropViews fluidProps = getFluidProperties( fluid );
  arrayView2d< real64 const > const & density = fluidProps.dens;
  arrayView2d< real64 const > const & dDens_dPres = fluidProps.dDens_dPres;

  ConstitutiveBase const & solid =
    GetConstitutiveModel( subRegion, solidModelNames()[targetIndex] );
  arrayView2d< real64 const > const & pvMult = solid.getReference< array2d< real64 > >(
    ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString );
  arrayView2d< real64 const > const & dPvMult_dPres =
    solid.getReference< array2d< real64 > >(
      ConstitutiveBase::viewKeyStruct::dPVMult_dPresString );
  arrayView1d< real64 const > const & bulkModulus =
    ISPORO ? solid.getReference< array1d< real64 > >( "BulkModulus" ) : porosityOld;
  real64 const biotCoefficient =
    ISPORO ? solid.getReference< real64 >( "BiotCoefficient" ) : 0.0;

  using Kernel = AccumulationKernel< CellElementSubRegion >;

  Kernel::template Launch< ISPORO, POLICY >( subRegion.size(),
                                             rankOffset,
                                             dofNumber,
                                             ghostRank,
                                             deltaPressure,
                                             densityOld,
                                             porosity,
                                             porosityOld,
                                             porosityRef,
                                             volume,
                                             deltaVolume,
                                             density,
                                             dDens_dPres,
                                             pvMult,
                                             dPvMult_dPres,
                                             totalMeanStressOld,
                                             totalMeanStress,
                                             bulkModulus,
                                             biotCoefficient,
                                             localMatrix,
                                             localRhs );
}

template< bool ISPORO, typename POLICY >
void
SinglePhaseBase::AccumulationLaunch(
  localIndex const targetIndex,
  FaceElementSubRegion const & subRegion,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );
  globalIndex const rankOffset = dofManager.rankOffset();
  arrayView1d< globalIndex const > const & dofNumber =
    subRegion.getReference< array1d< globalIndex > >( dofKey );
  arrayView1d< integer const > const & ghostRank = subRegion.ghostRank();

  arrayView1d< real64 const > const & densityOld =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
  arrayView1d< real64 const > const & volume = subRegion.getElementVolume();
  arrayView1d< real64 const > const & deltaVolume =
    subRegion.getReference< array1d< real64 > >( viewKeyStruct::deltaVolumeString );
  arrayView1d< real64 const > const & poroMult = getPoreVolumeMult( subRegion );

  ConstitutiveBase const & fluid =
    GetConstitutiveModel( subRegion, fluidModelNames()[targetIndex] );
  FluidPropViews const fluidProps = getFluidProperties( fluid );
  arrayView2d< real64 const > const & density = fluidProps.dens;
  arrayView2d< real64 const > const & dDens_dPres = fluidProps.dDens_dPres;

#if !defined( ALLOW_CREATION_MASS )
  static_assert( true, "must have ALLOW_CREATION_MASS defined" );
#endif

#if ALLOW_CREATION_MASS
  arrayView1d< real64 const > const & creationMass =
    subRegion.getReference< real64_array >(
      FaceElementSubRegion::viewKeyStruct::creationMassString );
#endif

  using Kernel = AccumulationKernel< FaceElementSubRegion >;

  Kernel::template Launch< ISPORO, POLICY >( subRegion.size(),
                                             rankOffset,
                                             dofNumber,
                                             ghostRank,
                                             densityOld,
                                             volume,
                                             deltaVolume,
                                             density,
                                             dDens_dPres,
                                             poroMult,
#if ALLOW_CREATION_MASS
                                             creationMass,
#endif
                                             localMatrix,
                                             localRhs );
}

template< bool ISPORO, typename POLICY >
void
SinglePhaseBase::AssembleAccumulationTerms(
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions< CellElementSubRegion, FaceElementSubRegion >(
    mesh,
    [&]( localIndex const targetIndex, auto & subRegion ) {
      AccumulationLaunch< ISPORO, POLICY >( targetIndex,
                                            subRegion,
                                            dofManager,
                                            localMatrix,
                                            localRhs );
    } );
}

void
SinglePhaseBase::ApplyBoundaryConditions(
  real64 time_n,
  real64 dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  ApplySourceFluxBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
  ApplyDiricletBC( time_n, dt, domain, dofManager, localMatrix, localRhs );
}

void
SinglePhaseBase::ApplyDiricletBC(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  fsManager.Apply(
    time_n + dt,
    &domain,
    "ElementRegions",
    viewKeyStruct::pressureString,
    [&]( FieldSpecificationBase const * const fs,
         string const &,
         SortedArrayView< localIndex const > const & lset,
         Group * subRegion,
         string const & ) {
      arrayView1d< globalIndex const > const & dofNumber =
        subRegion->getReference< array1d< globalIndex > >( dofKey );

      arrayView1d< real64 const > const & pres =
        subRegion->getReference< array1d< real64 > >( viewKeyStruct::pressureString );

      arrayView1d< real64 const > const & dPres =
        subRegion->getReference< array1d< real64 > >(
          viewKeyStruct::deltaPressureString );

      // call the application of the boundary condition to alter the matrix and rhs
      fs->ApplyBoundaryConditionToSystem< FieldSpecificationEqual,
                                          parallelDevicePolicy<> >(
        lset,
        time_n + dt,
        subRegion,
        dofNumber,
        dofManager.rankOffset(),
        localMatrix,
        localRhs,
        [=] GEOSX_HOST_DEVICE( localIndex const a ) { return pres[a] + dPres[a]; } );
    } );
}

void
SinglePhaseBase::ApplySourceFluxBC(
  real64 const time_n,
  real64 const dt,
  DomainPartition & domain,
  DofManager const & dofManager,
  CRSMatrixView< real64, globalIndex const > const & localMatrix,
  arrayView1d< real64 > const & localRhs ) const
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::pressureString );

  fsManager.Apply(
    time_n + dt,
    &domain,
    "ElementRegions",
    FieldSpecificationBase::viewKeyStruct::fluxBoundaryConditionString,
    [&]( FieldSpecificationBase const * const fs,
         string const &,
         SortedArrayView< localIndex const > const & lset,
         Group * subRegion,
         string const & ) -> void {
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
        time_n + dt,
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
SinglePhaseBase::SolveSystem( DofManager const & dofManager,
                              ParallelMatrix & matrix,
                              ParallelVector & rhs,
                              ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void
SinglePhaseBase::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      arrayView1d< real64 > const & dPres = subRegion.getReference< array1d< real64 > >(
        viewKeyStruct::deltaPressureString );

      dPres.setValues< parallelDevicePolicy<> >( 0.0 );

      UpdateState( subRegion, targetIndex );
    } );
}

void
SinglePhaseBase::BackupFields( MeshLevel & mesh ) const
{
  forTargetSubRegions(
    mesh,
    [&]( localIndex const targetIndex, ElementSubRegionBase & subRegion ) {
      ConstitutiveBase const & fluid =
        GetConstitutiveModel( subRegion, m_fluidModelNames[targetIndex] );
      arrayView2d< real64 const > const & dens = getFluidProperties( fluid ).dens;

      arrayView1d< real64 > const & poro =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityString );
      arrayView1d< real64 > const & densOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::densityOldString );
      arrayView1d< real64 > const & poroOld =
        subRegion.getReference< array1d< real64 > >( viewKeyStruct::porosityOldString );

      forAll< parallelDevicePolicy<> >( subRegion.size(),
                                        [=] GEOSX_HOST_DEVICE( localIndex const ei ) {
                                          densOld[ei] = dens[ei][0];
                                          poroOld[ei] = poro[ei];
                                        } );
    } );
}

void
SinglePhaseBase::ResetViews( MeshLevel & mesh )
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

  m_volume.clear();
  m_volume = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    ElementSubRegionBase::viewKeyStruct::elementVolumeString );
  m_volume.setName( string( "accessors/" ) +
                    ElementSubRegionBase::viewKeyStruct::elementVolumeString );

  m_deltaVolume.clear();
  m_deltaVolume = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::deltaVolumeString );
  m_deltaVolume.setName( getName() + "/accessors/" +
                         viewKeyStruct::deltaVolumeString );

  m_mobility.clear();
  m_mobility = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::mobilityString );
  m_mobility.setName( getName() + "/accessors/" + viewKeyStruct::mobilityString );

  m_dMobility_dPres.clear();
  m_dMobility_dPres = elemManager.ConstructArrayViewAccessor< real64, 1 >(
    viewKeyStruct::dMobility_dPressureString );
  m_dMobility_dPres.setName( getName() + "/accessors/" +
                             viewKeyStruct::dMobility_dPressureString );

  ResetViewsPrivate( elemManager );
}

void
SinglePhaseBase::ResetViewsPrivate( ElementRegionManager const & elemManager )
{
  m_density.clear();
  m_density = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SingleFluidBase::viewKeyStruct::densityString,
    targetRegionNames(),
    fluidModelNames() );
  m_density.setName( getName() + "/accessors/" +
                     SingleFluidBase::viewKeyStruct::densityString );

  m_dDens_dPres.clear();
  m_dDens_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SingleFluidBase::viewKeyStruct::dDens_dPresString,
    targetRegionNames(),
    fluidModelNames() );
  m_dDens_dPres.setName( getName() + "/accessors/" +
                         SingleFluidBase::viewKeyStruct::dDens_dPresString );

  m_viscosity.clear();
  m_viscosity = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SingleFluidBase::viewKeyStruct::viscosityString,
    targetRegionNames(),
    fluidModelNames() );
  m_viscosity.setName( getName() + "/accessors/" +
                       SingleFluidBase::viewKeyStruct::viscosityString );

  m_dVisc_dPres.clear();
  m_dVisc_dPres = elemManager.ConstructMaterialArrayViewAccessor< real64, 2 >(
    SingleFluidBase::viewKeyStruct::dVisc_dPresString,
    targetRegionNames(),
    fluidModelNames() );
  m_dVisc_dPres.setName( getName() + "/accessors/" +
                         SingleFluidBase::viewKeyStruct::dVisc_dPresString );

  m_transTMultiplier.clear();
  m_transTMultiplier = elemManager.ConstructArrayViewAccessor< R1Tensor, 1 >(
    viewKeyStruct::transTMultString );
  m_transTMultiplier.setName( getName() + "/accessors/" +
                              viewKeyStruct::transTMultString );
}

} /* namespace geosx */
