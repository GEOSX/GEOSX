/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2018, Lawrence Livermore National Security, LLC.
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
 * @file SinglePhaseFlow.cpp
 */

#include "SinglePhaseFlow.hpp"

#include "ArrayView.hpp"
#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "managers/BoundaryConditions/BoundaryConditionManager.hpp"
#include "managers/DomainPartition.hpp"
#include "managers/NumericalMethodsManager.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "MPI_Communications/CommunicationTools.hpp"
#include "systemSolverInterface/LinearSolverWrapper.hpp"
#include "systemSolverInterface/EpetraBlockSystem.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace systemSolverInterface;
using namespace multidimensionalArray;


SinglePhaseFlow::SinglePhaseFlow( const std::string& name,
                                  ManagedGroup * const parent ):
  FlowSolverBase(name, parent)
{
  m_numDofPerCell = 1;

  // set the blockID for the block system interface
  getLinearSystemRepository()->SetBlockID( BlockIDs::fluidPressureBlock, this->getName() );
}


void SinglePhaseFlow::FillDocumentationNode(  )
{
  FlowSolverBase::FillDocumentationNode();

  cxx_utilities::DocumentationNode * const docNode = this->getDocumentationNode();

  docNode->setName(SinglePhaseFlow::CatalogName());
  docNode->setSchemaType("Node");
  docNode->setShortDescription("A single phase flow solver");
}

void SinglePhaseFlow::FillOtherDocumentationNodes( dataRepository::ManagedGroup * const rootGroup )
{
  FlowSolverBase::FillOtherDocumentationNodes( rootGroup );

  DomainPartition * domain  = rootGroup->GetGroup<DomainPartition>( keys::domain );

  for( auto & mesh : domain->getMeshBodies()->GetSubGroups() )
  {
    MeshLevel * meshLevel = ManagedGroup::group_cast<MeshBody*>(mesh.second)->getMeshLevel(0);

    ElementRegionManager * const elemManager = meshLevel->getElemManager();

    elemManager->forCellBlocks( [&]( CellBlockSubRegion * const cellBlock ) -> void
      {
        cxx_utilities::DocumentationNode * const docNode = cellBlock->getDocumentationNode();

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.pressure.Key(),
                                    viewKeysSinglePhaseFlow.pressure.Key(),
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Fluid pressure",
                                    "Fluid pressure",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    0 );

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.deltaPressure.Key(),
                                    viewKeysSinglePhaseFlow.deltaPressure.Key(),
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Change in fluid pressure",
                                    "Change in fluid pressure",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    1 );

        docNode->AllocateChildNode( viewKeyStruct::deltaVolumeString,
                                    viewKeyStruct::deltaVolumeString,
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Change in fluid volume",
                                    "Change in fluid volume",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    1 );

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.density.Key(),
                                    viewKeysSinglePhaseFlow.density.Key(),
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Fluid density",
                                    "Fluid density",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    3 );

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.porosity.Key(),
                                    viewKeysSinglePhaseFlow.porosity.Key(),
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Porosity",
                                    "Porosity",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    3 );

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.oldPorosity.Key(),
                                    viewKeysSinglePhaseFlow.oldPorosity.Key(),
                                    -1,
                                    "real64_array",
                                    "real64_array",
                                    "Porosity old value",
                                    "Porosity old value",
                                    "",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    3 );

        docNode->AllocateChildNode( viewKeysSinglePhaseFlow.blockLocalDofNumber.Key(),
                                    viewKeysSinglePhaseFlow.blockLocalDofNumber.Key(),
                                    -1,
                                    "globalIndex_array",
                                    "globalIndex_array",
                                    "DOF index",
                                    "DOF index",
                                    "0",
                                    elemManager->getName(),
                                    1,
                                    0,
                                    3 );

      });

    {
      FaceManager * const faceManager = meshLevel->getFaceManager();
      cxx_utilities::DocumentationNode * const docNode = faceManager->getDocumentationNode();

      docNode->AllocateChildNode( viewKeysSinglePhaseFlow.facePressure.Key(),
                                  viewKeysSinglePhaseFlow.facePressure.Key(),
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Fluid pressure",
                                  "Fluid pressure",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeysSinglePhaseFlow.density.Key(),
                                  viewKeysSinglePhaseFlow.density.Key(),
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Fluid density",
                                  "Fluid density",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  3 );

      docNode->AllocateChildNode( viewKeysSinglePhaseFlow.viscosity.Key(),
                                  viewKeysSinglePhaseFlow.viscosity.Key(),
                                  -1,
                                  "real64_array",
                                  "real64_array",
                                  "Fluid viscosity",
                                  "Fluid viscosity",
                                  "",
                                  faceManager->getName(),
                                  1,
                                  0,
                                  3 );
    }
  }
}

void SinglePhaseFlow::InitializePreSubGroups(ManagedGroup * const rootGroup)
{
  FlowSolverBase::InitializePreSubGroups(rootGroup);
}

void SinglePhaseFlow::FinalInitializationPreSubGroups( ManagedGroup * const rootGroup )
{
  FlowSolverBase::FinalInitializationPreSubGroups( rootGroup );

  DomainPartition * domain = rootGroup->GetGroup<DomainPartition>(keys::domain);

  //TODO this is a hack until the sets are fixed to include ghosts!!
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysSinglePhaseFlow.pressure.Key());
  CommunicationTools::SynchronizeFields(fieldNames,
                              domain->getMeshBody(0)->getMeshLevel(0),
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  // Moved the following part from ImplicitStepSetup to here since it only needs to be initialized once
  // They will be updated in AyyplySystemSolution and ImplicitStepComplete, respectively
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  auto pres = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::pressureString);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elemManager->
                          ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  auto densOld = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::densityString);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> >
  pvmult = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                            constitutiveManager );

  auto poro = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::porosityString);
  auto poroOld = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::oldPorosityString);
  auto poroRef = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.referencePorosity.Key());

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  dens = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::
                                                                        viewKeyStruct::
                                                                        densityString,
                                                                        constitutiveManager);

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei)->void
  {
    constitutiveRelations[er][esr][m_fluidIndex]->StateUpdatePointPressure(pres[er][esr][ei], ei, 0); // fluid
    constitutiveRelations[er][esr][m_solidIndex]->StateUpdatePointPressure(pres[er][esr][ei], ei, 0); // solid

    densOld[er][esr][ei] = dens[er][esr][m_fluidIndex][ei][0];
    poro[er][esr][ei] = poroRef[er][esr][ei] * pvmult[er][esr][m_solidIndex][ei][0];
    poroOld[er][esr][ei] = poro[er][esr][ei];
  });
}



real64 SinglePhaseFlow::SolverStep( real64 const& time_n,
                                    real64 const& dt,
                                    const int cycleNumber,
                                    DomainPartition * domain )
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


void SinglePhaseFlow::ImplicitStepSetup( real64 const& time_n,
                                         real64 const& dt,
                                         DomainPartition * const domain,
                                         EpetraBlockSystem * const blockSystem )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  auto dPres = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());
  auto dVol      = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaVolume.Key());

  auto densOld = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.density.Key());
  auto poro = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.porosity.Key());
  auto poroOld = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.oldPorosity.Key());

  ElementRegionManager::MaterialViewAccessor< array2d<real64> >
  dens = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::densityString,
                                                                        constitutiveManager );

  //***** loop over all elements and initialize the derivative arrays *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei)->void
  {
    dPres[er][esr][ei] = 0.0;
    dVol[er][esr][ei] = 0.0;
    densOld[er][esr][ei] = dens[er][esr][m_fluidIndex][ei][0];
    poroOld[er][esr][ei] = poro[er][esr][ei];
  });

  // setup dof numbers and linear system
  SetupSystem( domain, blockSystem );
}

void SinglePhaseFlow::ImplicitStepComplete( real64 const & time_n,
                                            real64 const & dt,
                                            DomainPartition * const domain)
{

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto pres = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.pressure.Key());
  auto dPres = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());

  auto volume    = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);
  auto dVol      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaVolumeString);

  //***** Loop over all elements and update pressure and 'previous' values *****
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei)->void
  {
    pres[er][esr][ei] += dPres[er][esr][ei];
    volume[er][esr][ei] += dVol[er][esr][ei];
  });
}

void SinglePhaseFlow::SetNumRowsAndTrilinosIndices( MeshLevel * const meshLevel,
                                                    localIndex & numLocalRows,
                                                    globalIndex & numGlobalRows,
                                                    localIndex offset )
{

  ElementRegionManager * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<globalIndex_array>
  blockLocalDofNumber = elementRegionManager->
                  ConstructViewAccessor<globalIndex_array>( viewKeysSinglePhaseFlow.blockLocalDofNumber.Key(),
                                                            string() );

  ElementRegionManager::ElementViewAccessor< integer_array >
  ghostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  int numMpiProcesses;
  MPI_Comm_size( MPI_COMM_GEOSX, &numMpiProcesses );

  int thisMpiProcess = 0;
  MPI_Comm_rank( MPI_COMM_GEOSX, &thisMpiProcess );

  localIndex numLocalRowsToSend = numLocalRows;
  array1d<localIndex> gather(numMpiProcesses);

  // communicate the number of local rows to each process
  m_linearSolverWrapper.m_epetraComm.GatherAll( &numLocalRowsToSend,
                                                &gather.front(),
                                                1 );

  GEOS_ERROR_IF( numLocalRows != numLocalRowsToSend, "number of local rows inconsistent" );

  // find the first local row on this partition, and find the number of total global rows.
  localIndex firstLocalRow = 0;
  numGlobalRows = 0;

  for( integer p=0 ; p<numMpiProcesses ; ++p)
  {
    numGlobalRows += gather[p];
    if(p<thisMpiProcess)
      firstLocalRow += gather[p];
  }

  // create trilinos dof indexing, setting initial values to -1 to indicate unset values.
  for( localIndex er=0 ; er<ghostRank.size() ; ++er )
  {
    for( localIndex esr=0 ; esr<ghostRank[er].size() ; ++esr )
    {
      blockLocalDofNumber[er][esr] = -1;
    }
  }

  // loop over all elements and set the dof number if the element is not a ghost
  raja::ReduceSum< reducePolicy, localIndex  > localCount(0);
  forAllElemsInMesh<RAJA::seq_exec>( meshLevel, [=]( localIndex const er,
                                                     localIndex const esr,
                                                     localIndex const k) mutable ->void
  {
    if( ghostRank[er][esr][k] < 0 )
    {
      blockLocalDofNumber[er][esr][k] = firstLocalRow+localCount+offset;
      localCount += 1;
    }
    else
    {
      blockLocalDofNumber[er][esr][k] = -1;
    }
  });

  GEOS_ERROR_IF(localCount != numLocalRows, "Number of DOF assigned does not match numLocalRows" );
}


void SinglePhaseFlow::SetupSystem ( DomainPartition * const domain,
                                         EpetraBlockSystem * const blockSystem )
{
  // assume that there is only a single MeshLevel for now
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  // for this solver, the dof are on the cell center, and the row corrosponds to an element
  localIndex numGhostRows  = 0;
  localIndex numLocalRows  = 0;
  globalIndex numGlobalRows = 0;

  // get the number of local elements, and ghost elements...i.e. local rows and ghost rows
  elementRegionManager->forCellBlocks( [&]( CellBlockSubRegion * const subRegion )
    {
      localIndex subRegionGhosts = subRegion->GetNumberOfGhosts();
      numGhostRows += subRegionGhosts;
      numLocalRows += subRegion->size() - subRegionGhosts;
    });

  SetNumRowsAndTrilinosIndices( mesh,
                                numLocalRows,
                                numGlobalRows,
                                0 );

  //TODO element sync doesn't work yet
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());
  CommunicationTools::
  SynchronizeFields(fieldNames,
                    mesh,
                    domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );


  // construct row map, and set a pointer to the row map
  Epetra_Map * const
  rowMap = blockSystem->
           SetRowMap( BlockIDs::fluidPressureBlock,
                      std::make_unique<Epetra_Map>( numGlobalRows,
                                                    numLocalRows,
                                                    0,
                                                    m_linearSolverWrapper.m_epetraComm ) );

  // construct sparisty matrix, set a pointer to the sparsity pattern matrix
  Epetra_FECrsGraph * const
  sparsity = blockSystem->SetSparsity( BlockIDs::fluidPressureBlock,
                                       BlockIDs::fluidPressureBlock,
                                       std::make_unique<Epetra_FECrsGraph>(Copy,*rowMap,0) );



  // set the sparsity patter
  SetSparsityPattern( domain, sparsity );

  // assemble the global sparsity matrix
  sparsity->GlobalAssemble();
  sparsity->OptimizeStorage();

  // construct system matrix
  blockSystem->SetMatrix( BlockIDs::fluidPressureBlock,
                          BlockIDs::fluidPressureBlock,
                          std::make_unique<Epetra_FECrsMatrix>(Copy,*sparsity) );

  // construct solution vector
  blockSystem->SetSolutionVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

  // construct residual vector
  blockSystem->SetResidualVector( BlockIDs::fluidPressureBlock,
                                  std::make_unique<Epetra_FEVector>(*rowMap) );

}

void SinglePhaseFlow::SetSparsityPattern( DomainPartition const * const domain,
                                               Epetra_FECrsGraph * const sparsity )
{
  MeshLevel const * const meshLevel = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager const * const elementRegionManager = meshLevel->getElemManager();
  ElementRegionManager::ElementViewAccessor<globalIndex_array const>
  blockLocalDofNumber = elementRegionManager->
                  ConstructViewAccessor<globalIndex_array>( viewKeysSinglePhaseFlow.blockLocalDofNumber.Key() );

  ElementRegionManager::ElementViewAccessor< integer_array const >
  elemGhostRank = elementRegionManager->
              ConstructViewAccessor<integer_array>( ObjectManagerBase::viewKeyStruct::ghostRankString );

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  globalIndex_array elementLocalDofIndexRow(1);
  globalIndex_array elementLocalDofIndexCol(1);

  //**** loop over all faces. Fill in sparsity for all pairs of DOF/elem that are connected by face
  constexpr localIndex numElems = 2;
  stencilCollection.forAll<RAJA::seq_exec>([&] (StencilCollection<CellDescriptor, real64>::Accessor stencil) -> void
  {
    elementLocalDofIndexRow.resize(numElems);
    stencil.forConnected([&] (CellDescriptor const & cell, localIndex const i) -> void
    {
      elementLocalDofIndexRow[i] = blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
    });

    localIndex const stencilSize = stencil.size();
    elementLocalDofIndexCol.resize(stencilSize);
    stencil.forAll([&] (CellDescriptor const & cell, real64 w, localIndex const i) -> void
    {
      elementLocalDofIndexCol[i] = blockLocalDofNumber[cell.region][cell.subRegion][cell.index];
    });

    sparsity->InsertGlobalIndices(integer_conversion<int>(numElems),
                                  elementLocalDofIndexRow.data(),
                                  integer_conversion<int>(stencilSize),
                                  elementLocalDofIndexCol.data());
  });

  // loop over all elements and add all locals just in case the above connector loop missed some
  forAllElemsInMesh<RAJA::seq_exec>(meshLevel, [&] (localIndex const er,
                                    localIndex const esr,
                                    localIndex const k) -> void
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      elementLocalDofIndexRow[0] = blockLocalDofNumber[er][esr][k];

      sparsity->InsertGlobalIndices( 1,
                                     elementLocalDofIndexRow.data(),
                                     1,
                                     elementLocalDofIndexRow.data());
    }
  });

  // add additional connectivity resulting from boundary stencils
  fluxApprox->forBoundaryStencils([&] (FluxApproximationBase::BoundaryStencil const & boundaryStencilCollection) -> void
  {
    boundaryStencilCollection.forAll<RAJA::seq_exec>([=] (StencilCollection<PointDescriptor, real64>::Accessor stencil) mutable -> void
    {
      elementLocalDofIndexRow.resize(1);
      stencil.forConnected([&] (PointDescriptor const & point, localIndex const i) -> void
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          elementLocalDofIndexRow[0] = blockLocalDofNumber[c.region][c.subRegion][c.index];
        }
      });

      localIndex const stencilSize = stencil.size();
      elementLocalDofIndexCol.resize(stencilSize);
      integer counter = 0;
      stencil.forAll([&] (PointDescriptor const & point, real64 w, localIndex i) -> void
      {
        if (point.tag == PointDescriptor::Tag::CELL)
        {
          CellDescriptor const & c = point.cellIndex;
          elementLocalDofIndexCol[counter++] = blockLocalDofNumber[c.region][c.subRegion][c.index];
        }
      });

      sparsity->InsertGlobalIndices(1,
                                    elementLocalDofIndexRow.data(),
                                    integer_conversion<int>(counter),
                                    elementLocalDofIndexCol.data());
    });
  });
}



void SinglePhaseFlow::AssembleSystem(DomainPartition * const  domain,
                                     EpetraBlockSystem * const blockSystem,
                                     real64 const time_n,
                                     real64 const dt)
{
  //***** extract data required for assembly of system *****
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  ElementRegionManager * const elemManager = mesh->getElemManager();

  NumericalMethodsManager const * numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager const * fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * fluxApprox = fvManager->getFluxApproximation(m_discretizationName);
  FluxApproximationBase::CellStencil const & stencilCollection = fluxApprox->getStencil();

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(BlockIDs::fluidPressureBlock,
                                                               BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(BlockIDs::fluidPressureBlock);

  jacobian->Scale(0.0);
  residual->Scale(0.0);

  auto
  elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                     viewKeyStruct::
                                                                     ghostRankString );

  auto blockLocalDofNumber = elemManager->
                             ConstructViewAccessor<globalIndex_array>(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());

  auto pres      = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.pressure.Key());
  auto dPres     = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());
  auto densOld   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.density.Key());
  auto poro      = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.porosity.Key());
  auto poroOld   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.oldPorosity.Key());
  auto poroRef   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.referencePorosity.Key());
  auto gravDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.gravityDepth.Key());
  auto volume    = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);
  auto dVol      = elemManager->ConstructViewAccessor<real64_array>(viewKeyStruct::deltaVolumeString);
  auto totalMeanStress = elemManager->ConstructViewAccessor<real64_array>("totalMeanStress");
  auto oldTotalMeanStress = elemManager->ConstructViewAccessor<real64_array>("oldTotalMeanStress");

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  dens = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::densityString,
                                                                        constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  dDens_dPres = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::dDens_dPresString,
                                                                               constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  visc = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::viscosityString,
                                                                        constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  dVisc_dPres = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::dVisc_dPresString,
                                                                               constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  pvmult = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::poreVolumeMultiplierString,
                                                                          constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  dPVMult_dPres = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::dPVMult_dPresString,
                                                                                 constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
  bulkModulus = elemManager->
                ConstructMaterialViewAccessor< array2d<real64> >( "BulkModulus",
                                                                  constitutiveManager);

  ElementRegionManager::MaterialViewAccessor<real64> const
  biotCoefficient = elemManager->ConstructMaterialViewAccessor<real64>( "BiotCoefficient", constitutiveManager);

  //***** Loop over all elements and assemble the change in volume/density terms *****
//  forAllElemsInMesh(mesh, [=] (localIndex const er,
//                               localIndex const esr,
//                               localIndex const ei) -> void

  for( localIndex er=0 ; er<elemManager->numRegions() ; ++er )
  {
    ElementRegion const * const elemRegion = elemManager->GetRegion(er);
    for( localIndex esr=0 ; esr<elemRegion->numSubRegions() ; ++esr )
    {
      CellBlockSubRegion const * const cellBlockSubRegion = elemRegion->GetSubRegion(esr);

      for( localIndex ei=0 ; ei<cellBlockSubRegion->size() ; ++ei )
      {
        if (elemGhostRank[er][esr][ei]<0)
        {
          globalIndex const elemDOF = blockLocalDofNumber[er][esr][ei];

          real64 const densNew = dens[er][esr][m_fluidIndex][ei][0];
          real64 const volNew = volume[er][esr][ei] + dVol[er][esr][ei];

          if (m_poroElasticFlag == 0)
          {
            poro[er][esr][ei] = poroRef[er][esr][ei] * pvmult[er][esr][m_solidIndex][ei][0];
          }
          else
          {
            poro[er][esr][ei] = poroOld[er][esr][ei] + (biotCoefficient[er][esr][m_solidIndex] - poroOld[er][esr][ei]) / bulkModulus[er][esr][m_solidIndex][ei][0]
                                                     * (totalMeanStress[er][esr][ei] - oldTotalMeanStress[er][esr][ei] + dPres[er][esr][ei]);
          }

          // Residual contribution is mass conservation in the cell
          real64 const localAccum = poro[er][esr][ei]    * densNew              * volNew
                                  - poroOld[er][esr][ei] * densOld[er][esr][ei] * volume[er][esr][ei];

          // Derivative of residual wrt to pressure in the cell
          real64 localAccumJacobian;
          if (m_poroElasticFlag == 0)
          {
            localAccumJacobian = dPVMult_dPres[er][esr][m_solidIndex][ei][0] * poroRef[er][esr][ei] * densNew           * volNew
                               + dDens_dPres[er][esr][m_fluidIndex][ei][0]                          * poro[er][esr][ei] * volNew;
          }
          else
          {
            localAccumJacobian = (biotCoefficient[er][esr][m_solidIndex] - poroOld[er][esr][ei]) / bulkModulus[er][esr][m_solidIndex][ei][0] * densNew           * volNew
                               + dDens_dPres[er][esr][m_fluidIndex][ei][0]                                                                * poro[er][esr][ei] * volNew;
          }

          // add contribution to global residual and dRdP
          residual->SumIntoGlobalValues(1, &elemDOF, &localAccum);
          jacobian->SumIntoGlobalValues(1, &elemDOF, 1, &elemDOF, &localAccumJacobian);
        }
      }
    }
  }//);


  constexpr localIndex numElems = 2;
  globalIndex eqnRowIndices[numElems] = { -1, -1 };
  globalIndex_array dofColIndices;
  real64 localFlux[numElems] = { 0.0, 0.0 };
  array2d<real64> localFluxJacobian;

  // temporary working arrays
  real64 densWeight[numElems] = { 0.5, 0.5 };
  real64 mobility[numElems] = { 0.0, 0.0 };
  real64 dMobility_dP[numElems] = { 0.0, 0.0 };
  real64_array dDensMean_dP, dFlux_dP;

  stencilCollection.forAll<RAJA::seq_exec>([=] (StencilCollection<CellDescriptor, real64>::Accessor stencil) mutable -> void
  {
    localIndex const stencilSize = stencil.size();

    // resize and clear local working arrays
    dDensMean_dP.resize(stencilSize); // doesn't need to be that large, but it's convenient
    dFlux_dP.resize(stencilSize);

    // clear working arrays
    dDensMean_dP = 0.0;

    // resize local matrices and vectors
    dofColIndices.resize(stencilSize);
    localFluxJacobian.resize(numElems, stencilSize);

    // calculate quantities on primary connected cells
    real64 densMean = 0.0;
    stencil.forConnected([&] (auto const & cell,
                              localIndex i) -> void
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      eqnRowIndices[i] = blockLocalDofNumber[er][esr][ei];

      // density
      real64 const density   = dens[er][esr][m_fluidIndex][ei][0];
      real64 const dDens_dP  = dDens_dPres[er][esr][m_fluidIndex][ei][0];

      // viscosity
      real64 const viscosity = visc[er][esr][m_fluidIndex][ei][0];
      real64 const dVisc_dP  = dVisc_dPres[er][esr][m_fluidIndex][ei][0];

      // mobility
      mobility[i]  = density / viscosity;
      dMobility_dP[i]  = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

      // average density
      densMean += densWeight[i] * density;
      dDensMean_dP[i] = densWeight[i] * dDens_dP;
    });

    //***** calculation of flux *****

    // compute potential difference MPFA-style
    real64 potDif = 0.0;
    stencil.forAll([&] (CellDescriptor cell, real64 w, localIndex i) -> void
    {
      localIndex const er  = cell.region;
      localIndex const esr = cell.subRegion;
      localIndex const ei  = cell.index;

      dofColIndices[i] = blockLocalDofNumber[er][esr][ei];

      real64 const gravD    = gravDepth[er][esr][ei];
      real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
      real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

      potDif += w * (pres[er][esr][ei] + dPres[er][esr][ei] - gravTerm);
      dFlux_dP[i] = w * (1.0 - dGrav_dP);
    });

    // upwinding of fluid properties (make this an option?)
    localIndex const k_up = (potDif >= 0) ? 0 : 1;

    // compute the final flux and derivatives
    real64 const flux = mobility[k_up] * potDif;
    for (localIndex ke = 0; ke < stencilSize; ++ke)
      dFlux_dP[ke] *= mobility[k_up];
    dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

    //***** end flux terms *****

    // populate local flux vector and derivatives
    localFlux[0] =  dt * flux;
    localFlux[1] = -localFlux[0];

    for (localIndex ke = 0; ke < stencilSize; ++ke)
    {
      localFluxJacobian[0][ke] =   dt * dFlux_dP[ke];
      localFluxJacobian[1][ke] = - dt * dFlux_dP[ke];
    }

    // Add to global residual/jacobian
    jacobian->SumIntoGlobalValues(2, eqnRowIndices,
                                  integer_conversion<int>(stencilSize), dofColIndices.data(),
                                  localFluxJacobian.data());

    residual->SumIntoGlobalValues(2, eqnRowIndices, localFlux);
  });

  jacobian->GlobalAssemble(true);
  residual->GlobalAssemble();

  if( verboseLevel() >= 2 )
  {
    jacobian->Print(std::cout);
    residual->Print(std::cout);
  }

}

void SinglePhaseFlow::ApplyBoundaryConditions( DomainPartition * const domain,
                                               EpetraBlockSystem * const blockSystem,
                                               real64 const time_n,
                                               real64 const dt )
{
  // apply pressure boundary conditions.
  ApplyDirichletBC_implicit(domain, time_n, dt, blockSystem);
  ApplyFaceDirichletBC_implicit(domain, time_n, dt, blockSystem);

  if (verboseLevel() >= 2)
  {
    blockSystem->GetMatrix( BlockIDs::fluidPressureBlock, BlockIDs::fluidPressureBlock )->Print( std::cout );
    blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock )->Print( std::cout );
  }

}

/**
 * This function currently applies Dirichlet boundary conditions on the elements/zones as they
 * hold the DOF. Futher work will need to be done to apply a Dirichlet bc to the connectors (faces)
 */
void SinglePhaseFlow::ApplyDirichletBC_implicit( DomainPartition * domain,
                                                 real64 const time_n, real64 const dt,
                                                 EpetraBlockSystem * const blockSystem )
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();

  // call the BoundaryConditionManager::ApplyBoundaryCondition function that will check to see
  // if the boundary condition should be applied to this subregion
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "ElementRegions",
                                     viewKeysSinglePhaseFlow.pressure.Key(),
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const &,
                                          set<localIndex> const & lset,
                                          ManagedGroup * subRegion,
                                          string const & ) -> void
  {
    auto & dofMap = subRegion->getReference<array1d<globalIndex>>( viewKeysSinglePhaseFlow.blockLocalDofNumber );
    auto & pres   = subRegion->getReference<array1d<real64>>( viewKeysSinglePhaseFlow.pressure );
    auto & dPres  = subRegion->getReference<array1d<real64>>( viewKeysSinglePhaseFlow.deltaPressure );

    // call the application of the boundary condition to alter the matrix and rhs
    bc->ApplyBoundaryConditionToSystem<BcEqual>( lset,
                                                 time_n + dt,
                                                 subRegion,
                                                 dofMap,
                                                 1,
                                                 blockSystem,
                                                 BlockIDs::fluidPressureBlock,
                                                 [&] (localIndex const a) -> real64
    {
      return pres[a] + dPres[a];
    });
  });
}


void SinglePhaseFlow::ApplyFaceDirichletBC_implicit(DomainPartition * domain,
                                                    real64 const time_n, real64 const dt,
                                                    EpetraBlockSystem * const blockSystem)
{
  BoundaryConditionManager * bcManager = BoundaryConditionManager::get();
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();
  FaceManager * const faceManager = mesh->getFaceManager();

  array2d<localIndex> const & elemRegionList     = faceManager->elementRegionList();
  array2d<localIndex> const & elemSubRegionList  = faceManager->elementSubRegionList();
  array2d<localIndex> const & elemList           = faceManager->elementList();

  integer_array const & faceGhostRank = faceManager->getReference<integer_array>(ObjectManagerBase::
                                                                                 viewKeyStruct::
                                                                                 ghostRankString);

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  NumericalMethodsManager * const numericalMethodManager = domain->
    getParent()->GetGroup<NumericalMethodsManager>(keys::numericalMethodsManager);

  FiniteVolumeManager * const fvManager = numericalMethodManager->
    GetGroup<FiniteVolumeManager>(keys::finiteVolumeManager);

  FluxApproximationBase const * const fluxApprox = fvManager->getFluxApproximation(m_discretizationName);

  Epetra_FECrsMatrix * const jacobian = blockSystem->GetMatrix(BlockIDs::fluidPressureBlock,
                                                               BlockIDs::fluidPressureBlock);
  Epetra_FEVector * const residual = blockSystem->GetResidualVector(BlockIDs::fluidPressureBlock);

  auto
    elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                       viewKeyStruct::
                                                                       ghostRankString );

  auto blockLocalDofNumber = elemManager->
    ConstructViewAccessor<globalIndex_array>(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());

  auto pres      = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.pressure.Key());
  auto dPres     = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());
  auto densOld   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.density.Key());
  auto poroOld   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.porosity.Key());
  auto poroRef   = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.referencePorosity.Key());
  auto gravDepth = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.gravityDepth.Key());
  auto volume    = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
    dens = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::densityString,
                                                                          constitutiveManager );

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
    dDens_dPres = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::dDens_dPresString,
                                                                                 constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
    visc = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::viscosityString,
                                                                          constitutiveManager);

  ElementRegionManager::MaterialViewAccessor< array2d<real64> > const
    dVisc_dPres = elemManager->ConstructMaterialViewAccessor< array2d<real64> >( ConstitutiveBase::viewKeyStruct::dVisc_dPresString,
                                                                                 constitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
    constitutiveRelations = elemManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  // use ArrayView to make capture by value easy in lambdas
  ArrayView<real64, 1, localIndex> presFace = faceManager->getReference<real64_array>( viewKeysSinglePhaseFlow.facePressure );
  ArrayView<real64, 1, localIndex> densFace = faceManager->getReference<real64_array>( viewKeysSinglePhaseFlow.density);
  ArrayView<real64, 1, localIndex> viscFace = faceManager->getReference<real64_array>( viewKeysSinglePhaseFlow.viscosity );
  ArrayView<real64, 1, localIndex> gravDepthFace = faceManager->getReference<real64_array>( viewKeysSinglePhaseFlow.gravityDepth );

  dataRepository::ManagedGroup const * sets = faceManager->GetGroup(dataRepository::keys::sets);

  // first, evaluate BC to get primary field values (pressure)
//  bcManager->ApplyBoundaryCondition(faceManager, viewKeys.facePressure.Key(), time_n + dt);
  bcManager->ApplyBoundaryCondition( time_n + dt,
                                     domain,
                                     "faceManager",
                                     viewKeysSinglePhaseFlow.facePressure.Key(),
                                     [&]( BoundaryConditionBase const * const bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * const targetGroup,
                                          string const fieldName )->void
  {
    bc->ApplyBoundaryConditionToField<BcEqual>(targetSet,time_n + dt, targetGroup, fieldName);
  });


  // call constitutive models to get dependent quantities needed for flux (density, viscosity)
  bcManager->ApplyBoundaryCondition(time_n + dt,
                                    domain,
                                    "faceManager",
                                    viewKeysSinglePhaseFlow.facePressure.Key(),
                                    [&] ( BoundaryConditionBase const * bc,
                                          string const &,
                                          set<localIndex> const & targetSet,
                                          ManagedGroup * const,
                                          string const & ) -> void
  {
    for (auto kf : targetSet)
    {
      if (faceGhostRank[kf] >= 0)
        continue;

      // since we don't have models on faces yet, we take them from an adjacent cell
      integer const ke = (elemRegionList[kf][0] >= 0) ? 0 : 1;
      localIndex const er  = elemRegionList[kf][ke];
      localIndex const esr = elemSubRegionList[kf][ke];
      localIndex const ei  = elemList[kf][ke];

      real64 dummy; // don't need derivatives on faces
      constitutiveRelations[er][esr][m_fluidIndex]->FluidDensityCompute(presFace[kf], ei, densFace[kf], dummy);
      constitutiveRelations[er][esr][m_fluidIndex]->FluidViscosityCompute(presFace[kf], ei, viscFace[kf], dummy);
    }
  });

  // *** assembly loop ***

  constexpr localIndex numElems = 2;
  globalIndex_array dofColIndices;
  real64_array localFluxJacobian;

  // temporary working arrays
  real64 densWeight[numElems] = { 0.5, 0.5 };
  real64 mobility[numElems] = {0.0,0.0};
  real64 dMobility_dP[numElems] = {0.0,0.0};
  real64_array dDensMean_dP, dFlux_dP;


  bcManager->ApplyBoundaryCondition(time_n + dt,
                                    domain,
                                    "faceManager",
                                    viewKeysSinglePhaseFlow.facePressure.Key(),
                                    [&]( BoundaryConditionBase const * bc,
                                         string const & setName,
                                         set<localIndex> const &,
                                         ManagedGroup * const,
                                         string const & ) -> void
  {
    if (!sets->hasView(setName) || !fluxApprox->hasBoundaryStencil(setName))
      return;

    FluxApproximationBase::BoundaryStencil const & stencilCollection = fluxApprox->getBoundaryStencil(setName);

    stencilCollection.forAll([=] (StencilCollection<PointDescriptor, real64>::Accessor stencil) mutable -> void
    {
      localIndex const stencilSize = stencil.size();

      // resize and clear local working arrays
      dDensMean_dP.resize(stencilSize); // doesn't need to be that large, but it's convenient
      dFlux_dP.resize(stencilSize);

      // clear working arrays
      dDensMean_dP = 0.0;

      // resize local matrices and vectors
      dofColIndices.resize(stencilSize);
      localFluxJacobian.resize(stencilSize);

      // calculate quantities on primary connected points
      real64 densMean = 0.0;
      globalIndex eqnRowIndex = -1;
      localIndex cell_order = -1;
      stencil.forConnected([&] (PointDescriptor const & point, localIndex i) -> void
      {
        real64 density = 0, dDens_dP = 0;
        real64 viscosity = 0, dVisc_dP = 0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er  = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei  = point.cellIndex.index;

            eqnRowIndex = blockLocalDofNumber[er][esr][ei];

            density   = dens[er][esr][m_fluidIndex][ei][0];
            dDens_dP  = dDens_dPres[er][esr][m_fluidIndex][ei][0];

            viscosity = visc[er][esr][m_fluidIndex][ei][0];
            dVisc_dP  = dVisc_dPres[er][esr][m_fluidIndex][ei][0];

            cell_order = i; // mark position of the cell in connection for sign consistency later
            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            localIndex const kf = point.faceIndex;

            density = densFace[kf];
            dDens_dP = 0.0;

            viscosity = viscFace[kf];
            dVisc_dP = 0.0;

            break;
          }
          default:
            GEOS_ERROR("Unsupported point type in stencil");
        }

        // mobility
        mobility[i]  = density / viscosity;
        dMobility_dP[i]  = dDens_dP / viscosity - mobility[i] / viscosity * dVisc_dP;

        // average density
        densMean += densWeight[i] * density;
        dDensMean_dP[i] = densWeight[i] * dDens_dP;
      });

      //***** calculation of flux *****

      // compute potential difference MPFA-style
      real64 potDif = 0.0;
      dofColIndices = -1;
      stencil.forAll([&] (PointDescriptor point, real64 w, localIndex i) -> void
      {
        real64 pressure = 0.0, gravD = 0.0;
        switch (point.tag)
        {
          case PointDescriptor::Tag::CELL:
          {
            localIndex const er = point.cellIndex.region;
            localIndex const esr = point.cellIndex.subRegion;
            localIndex const ei = point.cellIndex.index;

            dofColIndices[i] = blockLocalDofNumber[er][esr][ei];
            pressure = pres[er][esr][ei] + dPres[er][esr][ei];
            gravD = gravDepth[er][esr][ei];

            break;
          }
          case PointDescriptor::Tag::FACE:
          {
            localIndex const kf = point.faceIndex;

            pressure = presFace[kf];
            gravD = gravDepthFace[kf];

            break;
          }
          default:
          GEOS_ERROR("Unsupported point type in stencil");
        }

        real64 const gravTerm = m_gravityFlag ? densMean * gravD : 0.0;
        real64 const dGrav_dP = m_gravityFlag ? dDensMean_dP[i] * gravD : 0.0;

        potDif += w * (pressure + gravTerm);
        dFlux_dP[i] = w * (1.0 + dGrav_dP);
      });

      // upwinding of fluid properties (make this an option?)
      localIndex const k_up = (potDif >= 0) ? 0 : 1;

      // compute the final flux and derivatives
      real64 const flux = mobility[k_up] * potDif;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
        dFlux_dP[ke] *= mobility[k_up];
      dFlux_dP[k_up] += dMobility_dP[k_up] * potDif;

      //***** end flux terms *****

      // populate local flux vector and derivatives
      integer sign = (cell_order == 0 ? 1 : -1);
      real64 const localFlux =  dt * flux * sign;

      integer counter = 0;
      for (localIndex ke = 0; ke < stencilSize; ++ke)
      {
        // compress arrays, skipping face derivatives
        if (dofColIndices[ke] >= 0)
        {
          dofColIndices[counter] = dofColIndices[ke];
          localFluxJacobian[counter] = dt * dFlux_dP[ke] * sign;
          ++counter;
        }
      }

      // Add to global residual/jacobian
      jacobian->SumIntoGlobalValues(1, &eqnRowIndex,
                                    counter, dofColIndices.data(),
                                    localFluxJacobian.data());

      residual->SumIntoGlobalValues(1, &eqnRowIndex, &localFlux);
    });
  });
}


real64
SinglePhaseFlow::
CalculateResidualNorm( EpetraBlockSystem const * const blockSystem,
                       DomainPartition * const domain )
{

  Epetra_FEVector const * const residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  Epetra_Map      const * const rowMap   = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );

  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elemManager = mesh->getElemManager();

  auto elemGhostRank = elemManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                          viewKeyStruct::
                                                                          ghostRankString );

  auto blockLocalDofNumber = elemManager->
      ConstructViewAccessor<globalIndex_array>(viewKeysSinglePhaseFlow.blockLocalDofNumber.Key());

  auto refPoro = elemManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.referencePorosity.Key());
  auto volume  = elemManager->ConstructViewAccessor<real64_array>(CellBlock::viewKeyStruct::elementVolumeString);

  // get a view into local residual vector
  int localSizeInt;
  double* localResidual = nullptr;
  residual->ExtractView(&localResidual, &localSizeInt);

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm = sumOverElemsInMesh(mesh, [&] (localIndex const er,
                                                           localIndex const esr,
                                                           localIndex const k) -> real64
  {
    if (elemGhostRank[er][esr][k] < 0)
    {
      int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][k]));
      real64 const val = localResidual[lid] / (refPoro[er][esr][k] * volume[er][esr][k]);
      return val * val;
    }
    return 0.0;
  });

  // compute global residual norm
  realT globalResidualNorm;
  MPI_Allreduce(&localResidualNorm, &globalResidualNorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_GEOSX);

  return sqrt(globalResidualNorm);
}

void SinglePhaseFlow::ApplySystemSolution( EpetraBlockSystem const * const blockSystem,
                                           real64 const scalingFactor,
                                           DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);

  Epetra_Map const * const rowMap        = blockSystem->GetRowMap( BlockIDs::fluidPressureBlock );
  Epetra_FEVector const * const solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  int dummy;
  double* local_solution = nullptr;
  solution->ExtractView(&local_solution,&dummy);

  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  ElementRegionManager::ElementViewAccessor<globalIndex_array>
  blockLocalDofNumber = elementRegionManager->
                        ConstructViewAccessor<globalIndex_array>( viewKeysSinglePhaseFlow.blockLocalDofNumber.Key() );

  auto pres  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.pressure.Key());
  auto dPres = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());

  auto elemGhostRank = elementRegionManager->ConstructViewAccessor<integer_array>( ObjectManagerBase::
                                                                                   viewKeyStruct::
                                                                                   ghostRankString );

  // loop over all elements to update incremental pressure
  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    if( elemGhostRank[er][esr][ei]<0 )
    {
      // extract solution and apply to dP
      int const lid = rowMap->LID(integer_conversion<int>(blockLocalDofNumber[er][esr][ei]));
      dPres[er][esr][ei] += scalingFactor * local_solution[lid];
    }
  });


  // TODO Sync dP once element field syncing is reimplemented.
  std::map<string, string_array > fieldNames;
  fieldNames["elems"].push_back(viewKeysSinglePhaseFlow.deltaPressure.Key());
  CommunicationTools::SynchronizeFields(fieldNames,
                              mesh,
                              domain->getReference< array1d<NeighborCommunicator> >( domain->viewKeys.neighbors ) );

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elementRegionManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei)->void
  {
    real64 const new_pres = pres[er][esr][ei] + dPres[er][esr][ei];
    constitutiveRelations[er][esr][m_fluidIndex]->StateUpdatePointPressure( new_pres, ei, 0 );
    constitutiveRelations[er][esr][m_solidIndex]->StateUpdatePointPressure( new_pres, ei, 0 );
  });
}

void SinglePhaseFlow::SolveSystem( EpetraBlockSystem * const blockSystem,
                                        SystemSolverParameters const * const params )
{
  Epetra_FEVector * const
  solution = blockSystem->GetSolutionVector( BlockIDs::fluidPressureBlock );

  Epetra_FEVector * const
  residual = blockSystem->GetResidualVector( BlockIDs::fluidPressureBlock );
  residual->Scale(-1.0);

  solution->Scale(0.0);

  m_linearSolverWrapper.SolveSingleBlockSystem( blockSystem,
                                                params,
                                                BlockIDs::fluidPressureBlock );

  if( verboseLevel() >= 2 )
  {
    solution->Print(std::cout);
  }

}

void SinglePhaseFlow::ResetStateToBeginningOfStep( DomainPartition * const domain )
{
  MeshLevel * const mesh = domain->getMeshBodies()->GetGroup<MeshBody>(0)->getMeshLevel(0);
  ElementRegionManager * const elementRegionManager = mesh->getElemManager();

  auto pres  = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.pressure.Key());
  auto dPres = elementRegionManager->ConstructViewAccessor<real64_array>(viewKeysSinglePhaseFlow.deltaPressure.Key());

  ConstitutiveManager * const
  constitutiveManager = domain->GetGroup<ConstitutiveManager>(keys::ConstitutiveManager);

  ElementRegionManager::ConstitutiveRelationAccessor<ConstitutiveBase>
  constitutiveRelations = elementRegionManager->ConstructConstitutiveAccessor<ConstitutiveBase>(constitutiveManager);

  forAllElemsInMesh( mesh, [&]( localIndex const er,
                                localIndex const esr,
                                localIndex const ei )->void
  {
    dPres[er][esr][ei] = 0.0;
    constitutiveRelations[er][esr][m_fluidIndex]->StateUpdatePointPressure( pres[er][esr][ei], ei, 0 );
    constitutiveRelations[er][esr][m_solidIndex]->StateUpdatePointPressure( pres[er][esr][ei], ei, 0 );
  });

}


REGISTER_CATALOG_ENTRY( SolverBase, SinglePhaseFlow, std::string const &, ManagedGroup * const )
} /* namespace geosx */
