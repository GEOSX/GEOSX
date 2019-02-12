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
 * @file CompositionalMultiphaseWell.cpp
 */

#include "CompositionalMultiphaseWell.hpp"

#include "managers/FieldSpecification/FieldSpecificationManager.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/Fluid/MultiFluidBase.hpp"
#include "constitutive/RelPerm/RelativePermeabilityBase.hpp"
#include "dataRepository/ManagedGroup.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "wells/WellManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/NeighborCommunicator.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;

CompositionalMultiphaseWell::CompositionalMultiphaseWell( const string & name,
                                                                      ManagedGroup * const parent )
  :
  WellSolverBase( name, parent ),
  m_numPhases( 0 ),
  m_numComponents( 0 )
{
  this->RegisterViewWrapper( viewKeyStruct::temperatureString, &m_temperature, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Temperature");

  this->RegisterViewWrapper( viewKeyStruct::useMassFlagString, &m_useMass, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Use mass formulation instead of molar");

  this->RegisterViewWrapper( viewKeyStruct::relPermNameString,  &m_relPermName,  false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("Name of the relative permeability constitutive model to use");

  this->RegisterViewWrapper( viewKeyStruct::relPermIndexString, &m_relPermIndex, false );
  
}

localIndex CompositionalMultiphaseWell::numFluidComponents() const
{
  return m_numComponents;
}

localIndex CompositionalMultiphaseWell::numFluidPhases() const
{
  return m_numPhases;
}


void CompositionalMultiphaseWell::InitializePreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePreSubGroups( rootGroup );
}

MultiFluidBase * CompositionalMultiphaseWell::GetFluidModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

MultiFluidBase const * CompositionalMultiphaseWell::GetFluidModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  MultiFluidBase const * const fluid = constitutiveModels->GetGroup<MultiFluidBase>( m_fluidName );
  GEOS_ERROR_IF( fluid == nullptr, "Target group does not contain the fluid model" );

  return fluid;
}

RelativePermeabilityBase * CompositionalMultiphaseWell::GetRelPermModel( ManagedGroup * dataGroup ) const
{
  ManagedGroup * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

RelativePermeabilityBase const * CompositionalMultiphaseWell::GetRelPermModel( ManagedGroup const * dataGroup ) const
{
  ManagedGroup const * const constitutiveModels = dataGroup->GetGroup( ConstitutiveManager::groupKeyStruct::constitutiveModelsString );
  GEOS_ERROR_IF( constitutiveModels == nullptr, "Target group does not contain constitutive models" );

  RelativePermeabilityBase const * const relPerm = constitutiveModels->GetGroup<RelativePermeabilityBase>( m_relPermName );
  GEOS_ERROR_IF( relPerm == nullptr, "Target group does not contain the relative permeability model" );

  return relPerm;
}

void CompositionalMultiphaseWell::UpdateComponentFraction( ManagedGroup * const dataGroup )
{
  
}

void CompositionalMultiphaseWell::UpdateComponentFractionAll( DomainPartition * const domain )
{
  
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFraction( ManagedGroup * const dataGroup )
{
  
}

void CompositionalMultiphaseWell::UpdatePhaseVolumeFractionAll( DomainPartition * const domain )
{
  
}

void CompositionalMultiphaseWell::UpdateFluidModel( ManagedGroup * const dataGroup )
{
  
}

void CompositionalMultiphaseWell::UpdateFluidModelAll( DomainPartition * const domain )
{

}
  
void CompositionalMultiphaseWell::UpdateRelPermModel( ManagedGroup * dataGroup )
{
  
}

void CompositionalMultiphaseWell::UpdateRelPermModelAll( DomainPartition * const domain )
{
  
}

void CompositionalMultiphaseWell::UpdateState( ManagedGroup * const dataGroup )
{
  UpdateComponentFraction( dataGroup );
  UpdateFluidModel( dataGroup );
  UpdatePhaseVolumeFraction( dataGroup );
  UpdateRelPermModel( dataGroup );
}

void CompositionalMultiphaseWell::UpdateStateAll( DomainPartition * const domain )
{
  UpdateComponentFractionAll( domain );
  UpdateFluidModelAll( domain );
  UpdatePhaseVolumeFractionAll( domain );
  UpdateRelPermModelAll( domain );
}

void CompositionalMultiphaseWell::InitializeWellState( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    // do something
  });
}

void CompositionalMultiphaseWell::InitializePostInitialConditions_PreSubGroups( ManagedGroup * const rootGroup )
{
  WellSolverBase::InitializePostInitialConditions_PreSubGroups( rootGroup );
  
}

real64 CompositionalMultiphaseWell::SolverStep( real64 const & time_n,
                                                      real64 const & dt,
                                                      integer const cycleNumber,
                                                      DomainPartition * const domain )
{
  real64 dt_return = dt;

  ImplicitStepSetup( time_n, dt, domain, getLinearSystemRepository() );

  // currently the only method is implicit time integration
  dt_return= this->NonlinearImplicitStep( time_n,
                                          dt,
                                          cycleNumber,
                                          domain,
                                          getLinearSystemRepository() );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt_return, domain );

  return dt_return;
}

void CompositionalMultiphaseWell::BackupFields( DomainPartition * const domain )
{
  
}

void
CompositionalMultiphaseWell::ImplicitStepSetup( real64 const & time_n, real64 const & dt,
                                                DomainPartition * const domain,
                                                EpetraBlockSystem * const blockSystem )
{

}

void CompositionalMultiphaseWell::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                                localIndex & numLocalRows,
                                                                globalIndex & numGlobalRows,
                                                                localIndex offset )
{
  
}

void CompositionalMultiphaseWell::SetSparsityPattern( DomainPartition const * const domain,
                                                      Epetra_FECrsGraph * const sparsity )
{
  
}

void CompositionalMultiphaseWell::SetupSystem( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem )
{
  
}

void CompositionalMultiphaseWell::AssembleSystem( DomainPartition * const domain,
                                                  EpetraBlockSystem * const blockSystem,
                                                  real64 const time_n, real64 const dt )
{ 
  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix( BlockIDs::compositionalBlock,
                                                                BlockIDs::compositionalBlock );
  Epetra_FEVector * const residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  // TODO: make sure to remove that if wells write into compositionalBlock
  jacobian->Scale(0.0);
  residual->Scale(0.0);

  CheckWellControlSwitch( domain );
  
  AssembleAccumulationTerms( domain, jacobian, residual, time_n, dt );
  AssembleFluxTerms( domain, jacobian, residual, time_n, dt );
  AssembleVolumeBalanceTerms( domain, jacobian, residual, time_n, dt );
  AssembleSourceTerms( domain, blockSystem, time_n, dt );
  
  jacobian->GlobalAssemble( true );
  residual->GlobalAssemble();

  if( verboseLevel() >= 3 )
  {
    GEOS_LOG_RANK("After CompositionalMultiphaseWell::AssembleSystem");
    GEOS_LOG_RANK("\nJacobian:\n" << *jacobian);
    GEOS_LOG_RANK("\nResidual:\n" << *residual);
  }

}

void CompositionalMultiphaseWell::AssembleAccumulationTerms( DomainPartition * const domain,
                                                             Epetra_FECrsMatrix * const jacobian,
                                                             Epetra_FEVector * const residual,
                                                             real64 const time_n,
                                                             real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    // loop over the segments
  });  
}

void CompositionalMultiphaseWell::AssembleFluxTerms( DomainPartition * const domain,
                                                     Epetra_FECrsMatrix * const jacobian,
                                                     Epetra_FEVector * const residual,
                                                     real64 const time_n,
                                                     real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    // loop over the connections
  });    
}

void CompositionalMultiphaseWell::AssembleVolumeBalanceTerms( DomainPartition * const domain,
                                                              Epetra_FECrsMatrix * const jacobian,
                                                              Epetra_FEVector * const residual,
                                                              real64 const time_n,
                                                              real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    // loop over the segments
  });    
}


void CompositionalMultiphaseWell::AssembleSourceTerms( DomainPartition * const domain,
                                                       EpetraBlockSystem * const blockSystem,
                                                       real64 const time_n,
                                                       real64 const dt )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    //  -- Compute the rates at each perforation
    //  -- Form the control equations
    //  -- Add to residual and Jacobian matrix
  });
}

void CompositionalMultiphaseWell::CheckWellControlSwitch( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    // check if the well control needs to be switched
  });
}


real64
CompositionalMultiphaseWell::CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                                                    DomainPartition * const domain )
{

  return 0.0;
}


void CompositionalMultiphaseWell::SolveSystem( EpetraBlockSystem * const blockSystem,
                                               SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
    solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );

  Epetra_FEVector * const
    residual = blockSystem->GetResidualVector( BlockIDs::compositionalBlock );

  residual->Scale(-1.0);
  solution->Scale(0.0);

  if( verboseLevel() >= 2 )
  {
    GEOS_LOG_RANK("\nSolution:\n" << *solution);
  }
}

bool
CompositionalMultiphaseWell::CheckSystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::compositionalBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::compositionalBlock );


  return false;
}

void
CompositionalMultiphaseWell::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                                  real64 const scalingFactor,
                                                  DomainPartition * const domain )
{
  
}

void
CompositionalMultiphaseWell::ApplyBoundaryConditions( DomainPartition * const domain,
                                                      systemSolverInterface::EpetraBlockSystem * const blockSystem,
                                                      real64 const time_n,
                                                      real64 const dt )
{

}

void CompositionalMultiphaseWell::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  WellManager * const wellManager = domain->getWellManager();

  wellManager->forSubGroups<CompositionalMultiphaseWell>( [&] ( CompositionalMultiphaseWell * well ) -> void
  {
    //   -- set dWellPres = 0;
    //   -- update the other well variables
  });
}

void CompositionalMultiphaseWell::ImplicitStepComplete( real64 const & time,
                                                              real64 const & dt,
                                                              DomainPartition * const domain )
{
  
}


REGISTER_CATALOG_ENTRY(SolverBase, CompositionalMultiphaseWell, string const &, ManagedGroup * const)
}// namespace geosx
