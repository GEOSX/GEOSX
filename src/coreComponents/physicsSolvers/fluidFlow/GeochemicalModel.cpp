/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file GeochemicalModel.cpp
 */

#include "GeochemicalModel.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/ReactiveFluidBase.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "mpiCommunications/CommunicationTools.hpp"
#include "mpiCommunications/NeighborCommunicator.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"


/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;

GeochemicalModel::GeochemicalModel( const std::string & name,
                                    Group * const parent ):
  FlowSolverBase( name, parent )
{

  this->registerWrapper( viewKeyStruct::reactiveFluidNamesString, &m_reactiveFluidNames )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Name of chemical system constitutive objects to use each target region." );

  this->registerWrapper( viewKeyStruct::outputSpeciesFileNameString, &m_outputSpeciesFileName )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Output species to file" );

}

void GeochemicalModel::RegisterDataOnMesh( Group * const MeshBodies )
{
  FlowSolverBase::RegisterDataOnMesh( MeshBodies );

  for( auto & mesh : MeshBodies->GetSubGroups() )
  {

    MeshLevel * meshLevel = Group::group_cast< MeshBody * >( mesh.second )->getMeshLevel( 0 );

    forTargetSubRegions( *meshLevel,
                         [&]
                           ( localIndex const GEOSX_UNUSED_PARAM( targetIndex ),
                           ElementSubRegionBase & subRegion )
    {
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::pressureString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaPressureString );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::temperatureString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array1d< real64 > >( viewKeyStruct::deltaTemperatureString );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::concentrationString )->setPlotLevel( PlotLevel::LEVEL_0 );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::concentrationOutString )->setPlotLevel( PlotLevel::LEVEL_0 );

      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::deltaConcentrationString );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::totalConcentrationString )->setPlotLevel( PlotLevel::LEVEL_0 );
      subRegion.registerWrapper< array2d< real64 > >( viewKeyStruct::concentrationNewString );

    } );

  }
}

void GeochemicalModel::InitializePreSubGroups( Group * const rootGroup )
{
  FlowSolverBase::InitializePreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );


  MeshLevel & meshLevel = *(domain->getMeshBody( 0 )->getMeshLevel( 0 ));

  m_numBasisSpecies.resize( m_reactiveFluidNames.size() );
  m_numDependentSpecies.resize( m_reactiveFluidNames.size() );
  this->forTargetSubRegions( meshLevel,
                             [&]
                               ( localIndex const targetRegionIndex,
                               ElementSubRegionBase const & subRegion )
  {
    string const & fluidName = m_reactiveFluidNames[targetRegionIndex];
    ReactiveFluidBase const & reactiveFluid = *(subRegion.getConstitutiveModel< ReactiveFluidBase >( fluidName ) );
    m_numBasisSpecies[targetRegionIndex]     = reactiveFluid.numBasisSpecies();
    m_numDependentSpecies[targetRegionIndex]  = reactiveFluid.numDependentSpecies();

    GEOSX_ERROR_IF( m_numDofPerCell!=0 || m_numDofPerCell==m_numBasisSpecies[targetRegionIndex],
                    "FlowSolverBase::m_numDofPerCell is set inconsistently. "
                    "Implementation not yet capable of having a different number of dof per cell across regions." );
    m_numDofPerCell = m_numBasisSpecies[targetRegionIndex];
  } );


  ResizeFields( &meshLevel );
}

void GeochemicalModel::ResizeFields( MeshLevel * const meshLevel )
{

  forTargetSubRegions( *meshLevel,
                       [&]
                         ( localIndex const targetRegionIndex,
                         ElementSubRegionBase & subRegion )
  {
    localIndex const NC = m_numBasisSpecies[targetRegionIndex];
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::concentrationString ).resizeDimension< 1 >( NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::concentrationOutString ).resizeDimension< 1 >( NC );

    subRegion.getReference< array2d< real64 > >( viewKeyStruct::deltaConcentrationString ).resizeDimension< 1 >( NC );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::totalConcentrationString ).resizeDimension< 1 >( NC );
    subRegion.getReference< array2d< real64 > >( viewKeyStruct::concentrationNewString ).resizeDimension< 1 >( NC );

  } );

}

void GeochemicalModel::UpdateReactiveFluidModel( Group * const dataGroup, localIndex const targetIndex )const
{
  GEOSX_MARK_FUNCTION;

  ReactiveFluidBase & reactiveFluid = GetConstitutiveModel< ReactiveFluidBase >( *dataGroup, m_reactiveFluidNames[targetIndex] );

  arrayView1d< real64 const > const & pres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::pressureString );
  arrayView1d< real64 const > const & dPres = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaPressureString );

  arrayView1d< real64 const > const & temp = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::temperatureString );
  arrayView1d< real64 const > const & dTemp = dataGroup->getReference< array1d< real64 > >( viewKeyStruct::deltaTemperatureString );

  arrayView2d< real64 const > const & conc = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::concentrationString );
  arrayView2d< real64 const > const & dConc = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::deltaConcentrationString );
  arrayView2d< real64 > & concNew = dataGroup->getReference< array2d< real64 > >( viewKeyStruct::concentrationNewString );

  forAll< serialPolicy >( dataGroup->size(), [&] ( localIndex const a )
  {
    for( localIndex ic = 0; ic < m_numBasisSpecies[targetIndex]; ++ic )
    {
      concNew[a][ic] = conc[a][ic] + dConc[a][ic];
    }
    reactiveFluid.PointUpdate( pres[a] + dPres[a], temp[a] + dTemp[a], concNew[a], a );
  } );

}

void GeochemicalModel::UpdateState( Group * dataGroup, localIndex const targetIndex ) const
{
  GEOSX_MARK_FUNCTION;

  UpdateReactiveFluidModel( dataGroup, targetIndex );

}

void GeochemicalModel::InitializePostInitialConditions_PreSubGroups( Group * const rootGroup )
{
  GEOSX_MARK_FUNCTION;

  FlowSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup< DomainPartition >( keys::domain );
  MeshLevel * mesh = domain->getMeshBody( 0 )->getMeshLevel( 0 );

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::pressureString ));
  fieldNames["elems"].emplace_back( string( viewKeyStruct::temperatureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::concentrationString ) );

  std::vector< NeighborCommunicator > & comms =
    domain->getNeighbors();

  CommunicationTools::SynchronizeFields( fieldNames, mesh, comms );

  ResetViews( *mesh );

}

real64 GeochemicalModel::SolverStep( real64 const & time_n,
                                     real64 const & dt,
                                     const int cycleNumber,
                                     DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  real64 dt_return = dt;

  ImplicitStepSetup( time_n,
                     dt,
                     domain );


  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;

}


void GeochemicalModel::ImplicitStepSetup( real64 const & time_n,
                                          real64 const & dt,
                                          DomainPartition & domain )
{
  GEOSX_UNUSED_VAR( time_n, dt );

  MeshLevel & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));
  ResetViews( mesh );

  /* The loop below could be moved to SolverStep after ImplicitStepSetup */

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {

    arrayView1d< real64 > const & dPres   = m_deltaPressure[er][esr];
    arrayView1d< real64 > const & dTemp   = m_deltaTemperature[er][esr];
    arrayView2d< real64 > const & dConc   = m_deltaConcentration[er][esr];
    arrayView2d< real64 > const & conc   = m_concentration[er][esr];
    arrayView2d< real64 > const & totalConc   = m_totalConcentration[er][esr];


    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      dPres[ei] = 0.0;
      dTemp[ei] = 0.0;
      for( localIndex ic = 0; ic < m_numBasisSpecies[targetRegionIndex]; ++ic )
      {
        dConc[ei][ic] = 0.0;
        conc[ei][ic] = log10( totalConc[ei][ic] );
      }

    } );

    UpdateState( &subRegion, targetRegionIndex );

  } );

  // setup dof numbers and linear system
  SetupSystem( domain,
               m_dofManager,
               m_localMatrix,
               m_localRhs,
               m_localSolution );

}

void GeochemicalModel::ImplicitStepComplete( real64 const & time_n,
                                             real64 const & GEOSX_UNUSED_PARAM( dt ),
                                             DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase & subRegion )
  {


    ReactiveFluidBase const & reactiveFluid  = GetConstitutiveModel< ReactiveFluidBase >( subRegion, m_reactiveFluidNames[targetRegionIndex] );

    arrayView1d< bool > const isHplus = reactiveFluid.IsHplus();

    arrayView1d< real64 > const & pres = m_pressure[er][esr];
    arrayView1d< real64 > const & temp = m_temperature[er][esr];

    arrayView2d< real64 > const & conc = m_concentration[er][esr];

    arrayView2d< real64 > const & concOut = m_concentrationOut[er][esr];

    arrayView1d< real64 const > const & dPres = m_deltaPressure[er][esr];
    arrayView1d< real64 const > const & dTemp = m_deltaTemperature[er][esr];
    arrayView2d< real64 const > const & dConc = m_deltaConcentration[er][esr];

    forAll< serialPolicy >( subRegion.size(),
                            [=]
                              ( localIndex const ei )
    {
      pres[ei] += dPres[ei];
      temp[ei] += dTemp[ei];
      for( localIndex ic = 0; ic < m_numBasisSpecies[targetRegionIndex]; ++ic )
      {
        conc[ei][ic] += dConc[ei][ic];

        if( isHplus[ic] )
          concOut[ei][ic] = -conc[ei][ic];
        else
          concOut[ei][ic] = pow( 10.0, conc[ei][ic] );

      }
    } );

  } );


  //  CalculateTotalConcentration( domain);

  WriteSpeciesToFile( &domain, time_n );

}

void GeochemicalModel::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                                  DofManager & dofManager ) const
{
  dofManager.addField( viewKeyStruct::geochemicalModelString,
                       DofManager::Location::Elem,
                       m_numDofPerCell,
                       targetRegionNames() );

  dofManager.addCoupling( viewKeyStruct::geochemicalModelString,
                          viewKeyStruct::geochemicalModelString,
                          DofManager::Connector::None );

}


void GeochemicalModel::AssembleSystem( real64 const time,
                                       real64 const dt,
                                       DomainPartition & domain,
                                       DofManager const & dofManager,
                                       CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                       arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  bool isInitialization = time <= 0.0 ? 1 : 0;

  isInitialization = 0;

  AssembleAccumulationTerms( domain, dofManager, localMatrix, localRhs, isInitialization );

  AssembleFluxTerms( time,
                     dt,
                     domain,
                     dofManager,
                     localMatrix,
                     localRhs );



  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK( "After GeochemicalModel::AssembleSystem" );
    GEOSX_LOG_RANK_0( "\nJacobian:\n" );
    std::cout<<localMatrix.toViewConst()<<std::endl;
    GEOSX_LOG_RANK_0( "\nResidual:\n" );
    std::cout<<localRhs.toViewConst()<<std::endl;
  }

}

void GeochemicalModel::AssembleAccumulationTerms( DomainPartition & domain,
                                                  DofManager const & dofManager,
                                                  CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                  arrayView1d< real64 > const & localRhs,
                                                  bool isInitialization )
{

  GEOSX_MARK_FUNCTION;


  MeshLevel const & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));
  globalIndex const rankOffset = dofManager.rankOffset();

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase const & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase const & subRegion )
  {


    ReactiveFluidBase const & reactiveFluid  = GetConstitutiveModel< ReactiveFluidBase >( subRegion,
                                                                                          m_reactiveFluidNames[targetRegionIndex] );

    arrayView2d< real64 > const stochMatrix = reactiveFluid.StochMatrix();
    arrayView1d< bool > const isHplus = reactiveFluid.IsHplus();


    string const dofKey = dofManager.getKey( viewKeyStruct::geochemicalModelString );
    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );

    arrayView1d< integer const >     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView2d< real64 const > const & conc          = m_concentration[er][esr];
    arrayView2d< real64 const > const & dConc         = m_deltaConcentration[er][esr];
    arrayView2d< real64 const > const & totalConc     = m_totalConcentration[er][esr];

    arrayView2d< real64 const > const & dependentConc     = m_dependentConc[er][esr];
    arrayView3d< real64 const > const & dDependentConc_dConc     = m_dDependentConc_dConc[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {
        stackArray1d< globalIndex, ReactiveFluidBase::MAX_NUM_SPECIES > localAccumDOF( m_numDofPerCell );
        stackArray1d< real64, ReactiveFluidBase::MAX_NUM_SPECIES > localAccum( m_numDofPerCell );
        stackArray2d< real64, ReactiveFluidBase::MAX_NUM_SPECIES * ReactiveFluidBase::MAX_NUM_SPECIES > localAccumJacobian( m_numDofPerCell, m_numDofPerCell );

        globalIndex const elemDOF = dofNumber[ei];

        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          localAccumDOF[idof] = elemDOF + idof;
        }

//        localAccumJacobian = 0.0;

        for( localIndex ic = 0; ic < m_numBasisSpecies[targetRegionIndex]; ++ic )
        {

          real64 concBasis = pow( 10.0, conc[ei][ic]+dConc[ei][ic] );

          localAccum[ic] = concBasis - totalConc[ei][ic];

          localAccumJacobian[ic][ic] = log( 10.0 ) * concBasis;

          if( isHplus[ic] && isInitialization )
          {
            localAccum[ic] = 0.0;
            localAccumJacobian[ic][ic] = 1.0;
            continue;
          }

          for( localIndex id = 0; id < m_numDependentSpecies[targetRegionIndex]; ++id )
          {
            real64 concDependent = pow( 10.0, dependentConc[ei][id] );
            localAccum[ic] -= stochMatrix[ic][id] * concDependent;

            for( localIndex idc = 0; idc < m_numBasisSpecies[targetRegionIndex]; ++idc )
            {

              localAccumJacobian[ic][idc] -= stochMatrix[ic][id] * log( 10.0 ) * concDependent * dDependentConc_dConc[ei][id][idc];

            }
          }
        }


        // add contribution to global residual and jacobian
//        rhs->add( localAccumDOF.data(),
//                  localAccum.data(),
//                  m_numDofPerCell );
//
//        matrix->add( localAccumDOF.data(),
//                     localAccumDOF.data(),
//                     localAccumJacobian.data(),
//                     m_numDofPerCell, m_numDofPerCell );

        for( localIndex i = 0; i < m_numDofPerCell; ++i )
        {
          localIndex const localRow = localAccumDOF[i] - rankOffset;
          localRhs[localRow] += localAccum[i];
          localMatrix.addToRow< serialAtomic >( localRow,
                                                localAccumDOF.data(),
                                                localAccumJacobian[i],
                                                m_numDofPerCell );
        }

      }

    } );
  } );

}


void GeochemicalModel::CalculateTotalConcentration( DomainPartition & domain )
{

  GEOSX_MARK_FUNCTION;


  MeshLevel const & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));
  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase const & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase const & subRegion )
  {

    ReactiveFluidBase const & reactiveFluid  = GetConstitutiveModel< ReactiveFluidBase >( subRegion,
                                                                                          m_reactiveFluidNames[targetRegionIndex] );

    arrayView2d< real64 > const stochMatrix = reactiveFluid.StochMatrix();
    arrayView1d< bool > const isHplus = reactiveFluid.IsHplus();

    arrayView1d< integer const >     const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView2d< real64 const > const & conc          = m_concentration[er][esr];
    arrayView2d< real64 > const & totalConc     = m_totalConcentration[er][esr];

    arrayView2d< real64 const > const & dependentConc     = m_dependentConc[er][esr];

    forAll< serialPolicy >( subRegion.size(), [=] ( localIndex const ei )
    {
      if( elemGhostRank[ei] < 0 )
      {

        for( localIndex ic = 0; ic < m_numBasisSpecies[targetRegionIndex]; ++ic )
        {


          real64 concBasis = pow( 10.0, conc[ei][ic] );

          real64 totalConcNew = concBasis;

          for( localIndex id = 0; id < m_numDependentSpecies[targetRegionIndex]; ++id )
          {
            real64 concDependent = pow( 10.0, dependentConc[ei][id] );
            totalConcNew -= stochMatrix[ic][id] * concDependent;

          }

          totalConc[ei][ic] = totalConcNew;

        }

      }

    } );
  } );

}

void GeochemicalModel::AssembleFluxTerms( real64 const time_n,
                                          real64 const dt,
                                          DomainPartition const & domain,
                                          DofManager const & dofManager,
                                          CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                          arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( time_n );
  GEOSX_UNUSED_VAR( dt );
  GEOSX_UNUSED_VAR( domain );
  GEOSX_UNUSED_VAR( dofManager );
  GEOSX_UNUSED_VAR( localMatrix );
  GEOSX_UNUSED_VAR( localRhs );
}

void GeochemicalModel::ApplyBoundaryConditions( real64 const time,
                                                real64 const dt,
                                                DomainPartition & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs )
{
  GEOSX_UNUSED_VAR( time );
  GEOSX_UNUSED_VAR( dt );
  GEOSX_UNUSED_VAR( domain );
  GEOSX_UNUSED_VAR( dofManager );
  GEOSX_UNUSED_VAR( localMatrix );
  GEOSX_UNUSED_VAR( localRhs );

}


real64
GeochemicalModel::
  CalculateResidualNorm( DomainPartition const & domain,
                         DofManager const & dofManager,
                         arrayView1d< real64 const > const & localRhs )
{

  MeshLevel const & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));

  string const dofKey = dofManager.getKey( viewKeyStruct::geochemicalModelString );
  globalIndex const rankOffset = dofManager.rankOffset();

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = 0.0;

  forTargetSubRegionsComplete( mesh,
                               [&]
                                 ( localIndex const,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase const & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase const & subRegion )
  {

    arrayView1d< globalIndex const > const & dofNumber = subRegion.getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    localIndex const subRegionSize = subRegion.size();
    for( localIndex a = 0; a < subRegionSize; ++a )
    {

      if( elemGhostRank[a] < 0 )
      {
        for( localIndex idof = 0; idof < m_numDofPerCell; ++idof )
        {
          localIndex const lid = dofNumber[a] - rankOffset;
          real64 const val = localRhs[lid];
          localResidualNorm += val * val;
        }
      }
    }

  } );

  // compute global residual norm
  real64 globalResidualNorm;
  MPI_Allreduce( &localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX );

  return sqrt( globalResidualNorm );
}


void GeochemicalModel::ApplySystemSolution( DofManager const & dofManager,
                                            arrayView1d< real64 const > const & localSolution,
                                            real64 const scalingFactor,
                                            DomainPartition & domain )
{

  MeshLevel & mesh = *(domain.getMeshBody( 0 )->getMeshLevel( 0 ));

  dofManager.addVectorToField( localSolution,
                               viewKeyStruct::geochemicalModelString,
                               viewKeyStruct::deltaConcentrationString,
                               scalingFactor,
                               0, m_numDofPerCell );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaPressureString ));
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaTemperatureString ) );
  fieldNames["elems"].emplace_back( string( viewKeyStruct::deltaConcentrationString ) );

  CommunicationTools::SynchronizeFields( fieldNames, &mesh, domain.getNeighbors() );

  this->forTargetSubRegions( mesh, [&] ( localIndex const targetRegionIndex,
                                         ElementSubRegionBase & subRegion )
  {
    UpdateState( &subRegion, targetRegionIndex );
  } );

}

void GeochemicalModel::SolveSystem( DofManager const & dofManager,
                                    ParallelMatrix & matrix,
                                    ParallelVector & rhs,
                                    ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  rhs.scale( -1.0 );
  solution.zero();

  SolverBase::SolveSystem( dofManager, matrix, rhs, solution );

  if( getLogLevel() >= 2 )
  {
    GEOSX_LOG_RANK( "After GeochemicalModel::SolveSystem" );
    GEOSX_LOG_RANK( "\nsolution\n" );
    std::cout<<solution<<std::endl;
  }

}

void GeochemicalModel::ResetStateToBeginningOfStep( DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{}

void GeochemicalModel::ResetViews( MeshLevel & mesh )
{
  FlowSolverBase::ResetViews( mesh );

  ElementRegionManager * const elemManager = mesh.getElemManager();

  m_pressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::pressureString );
  m_deltaPressure =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaPressureString );

  m_temperature =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::temperatureString );
  m_deltaTemperature =
    elemManager->ConstructViewAccessor< array1d< real64 >, arrayView1d< real64 > >( viewKeyStruct::deltaTemperatureString );

  m_concentration =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::concentrationString );
  m_deltaConcentration =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::deltaConcentrationString );
  m_totalConcentration =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::totalConcentrationString );

  m_concentrationNew =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::concentrationNewString );

  m_concentrationOut =
    elemManager->ConstructViewAccessor< array2d< real64 >, arrayView2d< real64 > >( viewKeyStruct::concentrationOutString );

  m_dependentConc =
    elemManager->ConstructMaterialArrayViewAccessor< real64, 2 >( ReactiveFluidBase::viewKeyStruct::dependentConcString,
                                                                  targetRegionNames(),
                                                                  m_reactiveFluidNames );

  m_dDependentConc_dConc =
    elemManager->ConstructMaterialArrayViewAccessor< real64, 3 >( ReactiveFluidBase::viewKeyStruct::dDependentConc_dConcString,
                                                                  targetRegionNames(),
                                                                  m_reactiveFluidNames );

}

void GeochemicalModel::WriteSpeciesToFile( DomainPartition * const domain, real64 const & time )
{

  GEOSX_MARK_FUNCTION;

  if( m_outputSpeciesFileName.empty())
    return;

  std::ofstream os( m_outputSpeciesFileName );
  GEOSX_ERROR_IF( !os.is_open(), "Cannot open the species-output file" );

  os << "==============================================" << std::endl;

  os << "time  =" << time << std::endl;

  os << "==============================================" << std::endl;


  MeshLevel const * const mesh = domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  forTargetSubRegionsComplete( *mesh,
                               [&]
                                 ( localIndex const targetRegionIndex,
                                 localIndex const er,
                                 localIndex const esr,
                                 ElementRegionBase const & GEOSX_UNUSED_PARAM( region ),
                                 ElementSubRegionBase const & subRegion )
  {
    arrayView1d< integer const > const & elemGhostRank = m_elemGhostRank[er][esr];

    arrayView2d< real64 const > const & conc          = m_concentration[er][esr];

    arrayView2d< real64 const > const & totalConc          = m_totalConcentration[er][esr];

    arrayView2d< real64 const > const & dependentConc     = m_dependentConc[er][esr];

    ReactiveFluidBase const & reactiveFluid  = *(subRegion.getConstitutiveModel< ReactiveFluidBase >( m_reactiveFluidNames[targetRegionIndex] ));
    const string_array & basisSpeciesNames = reactiveFluid.basisiSpeciesNames();
    const string_array & dependentSpeciesNames = reactiveFluid.dependentSpeciesNames();

    for( localIndex ei = 0; ei < subRegion.size(); ++ei )
    {

      if( elemGhostRank[ei] < 0 )
      {

        array1d< localIndex > indices;
        array1d< real64 > speciesConc;
        localIndex count  = 0;

        for( localIndex ic = 0; ic < m_numBasisSpecies[targetRegionIndex]; ++ic )
        {

          indices.emplace_back( count++ );

          speciesConc.emplace_back( conc[ei][ic] );

        }

        for( localIndex ic = 0; ic < m_numDependentSpecies[targetRegionIndex]; ++ic )
        {

          indices.emplace_back( count++ );

          speciesConc.emplace_back( dependentConc[ei][ic] );

        }

        std::sort( indices.begin(), indices.end(), [&]( localIndex i, localIndex j ){return speciesConc[i] > speciesConc[j];} );

        os << "   --- Distribution of Aqueous Solute Species ---" << std::endl;

        os << std::endl;

        os << "Species                   Molality            Log Molality" << std::endl;

        os << std::endl;

        for( localIndex ic = 0; ic < indices.size(); ic++ )
        {

          localIndex idx = indices[ic];
          real64 spC, spLogC;
          string spName;

          if( idx < m_numBasisSpecies[targetRegionIndex] )
          {
            spName = basisSpeciesNames[idx];
            spLogC = conc[ei][idx];
          }
          else
          {
            idx -= m_numBasisSpecies[targetRegionIndex];
            spName = dependentSpeciesNames[idx];
            spLogC = dependentConc[ei][idx];
          }

          auto found = spName.find( "(g)" );
          if( found != std::string::npos )
            continue;

          spC = pow( 10.0, spLogC );

          if( fabs( spLogC ) < 1e-64 || spC < 1e-40 )
            continue;

          os <<  std::left << std::setw( 25 ) << spName << std::setw( 10 ) << std::scientific << std::setprecision( 4 )<< std::right << spC << std::fixed << std::setw( 20 ) << spLogC << std::endl;

        }

        os << std::endl;
        os << std::endl;


        os << "            --- Total Conc ---" << std::endl;


        os << std::endl;
        os << std::endl;


        for( localIndex ic = 0; ic < indices.size(); ic++ )
        {

          localIndex idx = indices[ic];
          real64 spC;
          string spName;

          if( idx < m_numBasisSpecies[targetRegionIndex] )
          {
            spName = basisSpeciesNames[idx];
            spC = totalConc[ei][idx];
          }
          else
          {
            continue;
          }

          os <<  std::left << std::setw( 25 ) << spName << std::setw( 10 ) << std::scientific << std::setprecision( 4 )<< std::right << spC << std::endl;

        }

        os << std::endl;
        os << std::endl;


        os << std::endl;
        os << std::endl;

        os << "            --- Gas Fugacities ---" << std::endl;

        os << std::endl;

        os << "Gas                      Log Fugacity           Fugacity" << std::endl;

        os << std::endl;

        for( localIndex ic = 0; ic < indices.size(); ic++ )
        {

          localIndex idx = indices[ic];
          real64 spC, spLogC;
          string spName;

          if( idx < m_numBasisSpecies[targetRegionIndex] )
          {
            spName = basisSpeciesNames[idx];
            spLogC = conc[ei][idx];
          }
          else
          {
            idx -= m_numBasisSpecies[targetRegionIndex];
            spName = dependentSpeciesNames[idx];
            spLogC = dependentConc[ei][idx];
          }

          auto found = spName.find( "(g)" );
          if( found == std::string::npos )
            continue;

          spC = pow( 10.0, spLogC );

          os <<  std::left << std::setw( 20 ) << spName << std::setw( 15 ) << std::fixed << std::setprecision( 5 )<< std::right << spLogC << std::scientific << std::setw( 23 ) << spC << std::endl;

        }

        os << std::endl;
        os << std::endl;
        os << std::endl;

        os << "           --- Ionic Strength ---" << std::endl;

        os << std::endl;

        os <<  "Ionic Strength = " <<  std::fixed << std::setprecision( 4 ) << dependentConc[ei][m_numDependentSpecies[targetRegionIndex]] << std::endl;

      }
    }
  } );

  os.close();

}

REGISTER_CATALOG_ENTRY( SolverBase, GeochemicalModel, std::string const &, Group * const )
} /* namespace geosx */
