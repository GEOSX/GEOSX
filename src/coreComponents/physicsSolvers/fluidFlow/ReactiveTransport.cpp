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

#include "ReactiveTransport.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "physicsSolvers/fluidFlow/ReactiveTransportKernels.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace ReactiveTransportKernels;

ReactiveTransport::ReactiveTransport( const std::string & name,
                                      Group * const parent ):
  FlowSolverBase( name, parent )
{

  this->registerWrapper( viewKeyStruct::viscosityString, &m_viscosity )->setApplyDefaultValue( 0.001 )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Fluid viscosity" );

  this->registerWrapper( GeochemicalModel::viewKeyStruct::reactiveFluidNamesString, &m_reactiveFluidNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of chemical system constitutive objects to use each target region." );


}

void ReactiveTransport::PostProcessInput()
{
  FlowSolverBase::PostProcessInput();

}

void ReactiveTransport::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {
    MeshLevel & meshLevel = *Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions< CellElementSubRegion >( meshLevel, [&]( localIndex const,
                                                                 CellElementSubRegion & subRegion )
    {

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString )->setDefaultValue( 0.0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString )->setDefaultValue( 0.0 );

    } );
  }

}

void ReactiveTransport::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  ConstitutiveManager & cm = *domain->getConstitutiveManager();

  ReactiveFluidBase const * reactiveFluid = cm.GetConstitutiveRelation< ReactiveFluidBase >( m_reactiveFluidNames[0] );

  m_numComponents = reactiveFluid->numBasisSpecies();
  m_numDofPerCell = m_numComponents;

  MeshLevel & meshLevel = *domain->getMeshBody( 0 )->getMeshLevel( 0 );

  ResizeFields( meshLevel );

}

void ReactiveTransport::ResizeFields( MeshLevel & mesh )
{

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh,
                       [&]
                         ( localIndex const,
                         ElementSubRegionBase & subRegion )
  {

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString ).resizeDimension< 1 >( NC );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString ).resizeDimension< 1 >( NC );

  } );

}

void ReactiveTransport::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition & domain = *rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaComponentConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors(), true );

  ResetViews( mesh );

}

real64 ReactiveTransport::SolverStep( real64 const & time_n,
                                      real64 const & dt,
                                      const int cycleNumber,
                                      DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  ImplicitStepSetup( time_n, dt, domain );

  if( cycleNumber == 0 )
  {
    FieldSpecificationManager const & boundaryConditionManager = FieldSpecificationManager::get();

    boundaryConditionManager.ApplyInitialConditions( &domain );

  }

  PreStepUpdate( time_n, dt, domain );

  // currently the only method is implicit time integration
  real64 const dtReturn= this->NonlinearImplicitStep( time_n, dt, cycleNumber, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.

  ImplicitStepComplete( time_n, dtReturn, domain );

  PostStepUpdate( time_n, dtReturn, domain );

  return dtReturn;
}


void ReactiveTransport::PreStepUpdate( real64 const & GEOSX_UNUSED_PARAM( time ),
                                       real64 const & GEOSX_UNUSED_PARAM( dt ),
                                       DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  FlowSolverBase::PrecomputeData( mesh );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView2d< real64 > const deltaConc = subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        deltaConc[ei][c] = 0.0;
      }
    } );
  } );

}

void ReactiveTransport::PostStepUpdate( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                        real64 const & GEOSX_UNUSED_PARAM( dt_return ),
                                        DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void ReactiveTransport::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const & GEOSX_UNUSED_PARAM( dt ),
                                           DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ResetViews( mesh );

}

void ReactiveTransport::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                                              real64 const & GEOSX_UNUSED_PARAM( dt ),
                                              DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {

    arrayView2d< real64 > const componentConc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::totalConcentrationString );
    arrayView2d< real64 const > const dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {

      for( localIndex c = 0; c < NC; ++c )
      {
        componentConc[ei][c] += dComponentConc[ei][c];
      }
    } );
  } );

}

void ReactiveTransport::SetupDofs( DomainPartition const & domain,
                                   DofManager & dofManager ) const
{

  dofManager.addField( viewKeyStruct::reactiveTransportModelString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( viewKeyStruct::reactiveTransportModelString, fluxApprox );

}


void ReactiveTransport::AssembleSystem( real64 const time,
                                        real64 const dt,
                                        DomainPartition & domain,
                                        DofManager const & dofManager,
                                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                        arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  AssembleAccumulationTerms( dt,
                             domain,
                             dofManager,
                             localMatrix,
                             localRhs );

  AssembleFluxTerms( time,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );
}

void ReactiveTransport::AssembleAccumulationTerms( real64 const dt,
                                                   DomainPartition const & domain,
                                                   DofManager const & dofManager,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase const & subRegion )
  {
    string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume = subRegion.getElementVolume();

    arrayView2d< real64 const > const componentConc =
      subRegion.getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::totalConcentrationString );

    arrayView2d< real64 const > const deltaComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    arrayView1d< real64 const > const & porosity = subRegion.getReference< array1d< real64 > >( viewKeyStruct::referencePorosityString );

    AccumulationKernel::Launch( subRegion.size(),
                                m_numComponents,
                                m_numDofPerCell,
                                dofManager.rankOffset(),
                                dofNumber,
                                elemGhostRank,
                                componentConc,
                                deltaComponentConc,
                                porosity,
                                volume,
                                dt,
                                localMatrix,
                                localRhs );
  } );
}


void ReactiveTransport::AssembleFluxTerms( real64 const GEOSX_UNUSED_PARAM( time_n ),
                                           real64 const dt,
                                           DomainPartition const & domain,
                                           DofManager const & dofManager,
                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                           arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );
  ElementRegionManager const & elemManager = *mesh.getElemManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  dofNumberAccessor = elemManager.ConstructViewAccessor< array1d< globalIndex >, arrayView1d< globalIndex const > >( dofKey );

  FluxKernel::ElementViewConst< arrayView1d< globalIndex const > > const dofNumber = dofNumberAccessor.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const pres  = m_pressure.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView1d< real64 const > > const dPres = m_deltaPressure.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const componentConc   = m_componentConcentration.toNestedViewConst();
  FluxKernel::ElementViewConst< arrayView2d< real64 const > > const dComponentConc  = m_deltaComponentConcentration.toNestedViewConst();

  FluxKernel::ElementViewConst< arrayView1d< integer const > > const elemGhostRank = m_elemGhostRank.toNestedViewConst();

  fluxApprox.forStencils< CellElementStencilTPFA >( mesh, [&]( auto const & stencil )
  {

    FluxKernel::Launch( stencil,
                        m_numDofPerCell,
                        dt,
                        dofManager.rankOffset(),
                        dofNumber,
                        elemGhostRank,
                        pres,
                        dPres,
                        componentConc,
                        dComponentConc,
                        m_viscosity,
                        localMatrix,
                        localRhs );
  } );
}

void ReactiveTransport::ApplyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::get();
  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );
  globalIndex const rankOffset = dofManager.rankOffset();

  //  Apply Dirichlet BC for component concentration

  localIndex const NC = m_numComponents;

  if( NC > 0 )
  {
    map< string, map< string, array1d< bool > > > bcStatusMap; // map to check consistent application of BC

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::totalConcentrationString,
                     [&]( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( fs ),
                          string const & setName,
                          SortedArrayView< localIndex const > const & GEOSX_UNUSED_PARAM( targetSet ),
                          Group * const subRegion,
                          string const & )
    {

      string const & subRegionName = subRegion->getName();
      bcStatusMap[subRegionName][setName].resize( NC );
      bcStatusMap[subRegionName][setName].setValues< serialPolicy >( false );

    } );

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::totalConcentrationString,
                     [&] ( FieldSpecificationBase const * const fs,
                           string const & setName,
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * const subRegion,
                           string const & )
    {

      string const & subRegionName = subRegion->getName();
      localIndex const comp = fs->GetComponent();

      GEOSX_ERROR_IF( bcStatusMap[subRegionName][setName][comp], "Conflicting composition[" << comp << "] boundary conditions on set '" << setName << "'" );
      bcStatusMap[subRegionName][setName][comp] = true;

      fs->ApplyFieldValue< FieldSpecificationEqual >( targetSet,
                                                      time_n + dt,
                                                      subRegion,
                                                      viewKeyStruct::bcComponentConcentrationString );

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

    fsManager.Apply( time_n + dt,
                     &domain,
                     "ElementRegions",
                     GeochemicalModel::viewKeyStruct::totalConcentrationString,
                     [&] ( FieldSpecificationBase const * const GEOSX_UNUSED_PARAM( bc ),
                           string const & GEOSX_UNUSED_PARAM( setName ),
                           SortedArrayView< localIndex const > const & targetSet,
                           Group * const subRegion,
                           string const & )
    {
      arrayView1d< integer const > const ghostRank =
        subRegion->getReference< array1d< integer > >( ObjectManagerBase::viewKeyStruct::ghostRankString );
      arrayView1d< globalIndex const > const dofNumber = subRegion->getReference< array1d< globalIndex > >( dofKey );

      arrayView2d< real64 const > const compConc =
        subRegion->getReference< array2d< real64 > >( GeochemicalModel::viewKeyStruct::totalConcentrationString );
      arrayView2d< real64 const > const deltaCompConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );
      arrayView2d< real64 const > const bcCompConc =
        subRegion->getReference< array2d< real64 > >( viewKeyStruct::bcComponentConcentrationString );

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
          FieldSpecificationEqual::SpecifyFieldValue( dofIndex + ic,
                                                      rankOffset,
                                                      localMatrix,
                                                      rhsValue,
                                                      bcCompConc[ei][ic],
                                                      compConc[ei][ic] + deltaCompConc[ei][ic] );
          localRhs[localRow + ic] = rhsValue;
        }
      } );
    } );
  }
}

real64
ReactiveTransport::CalculateResidualNorm( DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NDOF = m_numDofPerCell;

  localIndex const rankOffset = dofManager.rankOffset();
  string const dofKey = dofManager.getKey( viewKeyStruct::reactiveTransportModelString );

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
  real64 const globalResidualNorm = MpiWrapper::Sum( localResidualNorm, MPI_COMM_GEOSX );
  return sqrt( globalResidualNorm );
}

void ReactiveTransport::ApplySystemSolution( DofManager const & dofManager,
                                             arrayView1d< real64 const > const & localSolution,
                                             real64 const scalingFactor,
                                             DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::reactiveTransportModelString,
                               viewKeyStruct::deltaComponentConcentrationString,
                               scalingFactor,
                               0, m_numDofPerCell );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaComponentConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors(), true );


}

void ReactiveTransport::SolveSystem( DofManager const & dofManager,
                                     ParallelMatrix & matrix,
                                     ParallelVector & rhs,
                                     ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );
}

void ReactiveTransport::ResetStateToBeginningOfStep( DomainPartition & domain )
{
  MeshLevel & mesh = *domain.getMeshBody( 0 )->getMeshLevel( 0 );

  localIndex const NC = m_numComponents;

  forTargetSubRegions( mesh, [&]( localIndex const, ElementSubRegionBase & subRegion )
  {
    arrayView2d< real64 > const & dComponentConc =
      subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaComponentConcentrationString );

    forAll< parallelDevicePolicy<> >( subRegion.size(), [=] GEOSX_HOST_DEVICE ( localIndex const ei )
    {
      for( localIndex c = 0; c < NC; ++c )
      {
        dComponentConc[ei][c] = 0.0;
      }
    } );

  } );

}

void ReactiveTransport::ResetViews( MeshLevel & mesh )
{
  FlowSolverBase::ResetViews( mesh );
  ElementRegionManager & elemManager = *mesh.getElemManager();

  m_pressure.clear();
  m_pressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::pressureString );
  m_pressure.setName( getName() + "/accessors/" + viewKeyStruct::pressureString );

  m_deltaPressure.clear();
  m_deltaPressure = elemManager.ConstructArrayViewAccessor< real64, 1 >( viewKeyStruct::deltaPressureString );
  m_deltaPressure.setName( getName() + "/accessors/" + viewKeyStruct::deltaPressureString );

  m_componentConcentration.clear();
  m_componentConcentration = elemManager.ConstructArrayViewAccessor< real64, 2 >( GeochemicalModel::viewKeyStruct::totalConcentrationString );
  m_componentConcentration.setName( getName() + "/accessors/" + GeochemicalModel::viewKeyStruct::totalConcentrationString );

  m_deltaComponentConcentration.clear();
  m_deltaComponentConcentration = elemManager.ConstructArrayViewAccessor< real64, 2 >( viewKeyStruct::deltaComponentConcentrationString );
  m_deltaComponentConcentration.setName( getName() + "/accessors/" + viewKeyStruct::deltaComponentConcentrationString );

}

REGISTER_CATALOG_ENTRY( SolverBase, ReactiveTransport, std::string const &, Group * const )
} /* namespace geosx */
