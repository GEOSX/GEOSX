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
 * @file SinglePhaseFVM.cpp
 */

#include "SinglePhaseFVM.hpp"

#include "mesh/mpiCommunications/CommunicationTools.hpp"
#include "common/TimingMacros.hpp"
#include "constitutive/fluid/singleFluidSelector.hpp"
#include "constitutive/solid/RockBase.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"
#include "discretizationMethods/NumericalMethodsManager.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "finiteVolume/BoundaryStencil.hpp"
#include "finiteVolume/FiniteVolumeManager.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseBaseKernels.hpp"
#include "physicsSolvers/fluidFlow/SinglePhaseFVMKernels.hpp"
#include "physicsSolvers/multiphysics/SinglePhasePoroelasticFluxKernels.hpp"
#include "physicsSolvers/solidMechanics/SolidMechanicsEmbeddedFractures.hpp"

/**
 * @namespace the geosx namespace that encapsulates the majority of the code
 */
namespace geosx
{

using namespace dataRepository;
using namespace constitutive;
using namespace SinglePhaseBaseKernels;
using namespace SinglePhaseFVMKernels;
using namespace SinglePhasePoroelasticFluxKernels;

template< typename BASE >
SinglePhaseFVM< BASE >::SinglePhaseFVM( const string & name,
                                        Group * const parent ):
  BASE( name, parent )
{
  m_numDofPerCell = 1;
}

template< typename BASE >
void SinglePhaseFVM< BASE >::initializePreSubGroups()
{
  BASE::initializePreSubGroups();

  DomainPartition & domain = this->template getGroupByPath< DomainPartition >( "/Problem/domain" );
  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();

  if( !fvManager.hasGroup< FluxApproximationBase >( m_discretizationName ) )
  {
    GEOSX_ERROR( "A discretization deriving from FluxApproximationBase must be selected with SinglePhaseFVM" );
  }
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupDofs( DomainPartition const & domain,
                                        DofManager & dofManager ) const
{
  dofManager.addField( BASE::viewKeyStruct::pressureString(),
                       DofManager::Location::Elem,
                       1,
                       targetRegionNames() );

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  dofManager.addCoupling( BASE::viewKeyStruct::pressureString(), fluxApprox );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::setupSystem( DomainPartition & domain,
                                          DofManager & dofManager,
                                          CRSMatrix< real64, globalIndex > & localMatrix,
                                          array1d< real64 > & localRhs,
                                          array1d< real64 > & localSolution,
                                          bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;
  BASE::setupSystem( domain,
                     dofManager,
                     localMatrix,
                     localRhs,
                     localSolution,
                     setSparsity );

}

template< typename BASE >
real64 SinglePhaseFVM< BASE >::calculateResidualNorm( DomainPartition const & domain,
                                                      DofManager const & dofManager,
                                                      arrayView1d< real64 const > const & localRhs )
{
  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  string const dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  globalIndex const rankOffset = dofManager.rankOffset();

  // compute the norm of local residual scaled by cell pore volume
  real64 localResidualNorm[3] = { 0.0, 0.0, 0.0 };
  forTargetSubRegions( mesh, [&]( localIndex const targetIndex,
                                  auto const & subRegion )
  {
    arrayView1d< globalIndex const > const & dofNumber = subRegion.template getReference< array1d< globalIndex > >( dofKey );
    arrayView1d< integer const > const & elemGhostRank = subRegion.ghostRank();
    arrayView1d< real64 const > const & volume         = subRegion.getElementVolume();
    arrayView1d< real64 const > const & densOld        = subRegion.template getReference< array1d< real64 > >( BASE::viewKeyStruct::densityOldString() );

    RockBase const & solidModel = this->template getConstitutiveModel< RockBase >( subRegion,
                                                                                   m_solidModelNames[targetIndex] );
    arrayView2d< real64 const > const & poroOld = solidModel.getOldPorosity();

    ResidualNormKernel::launch< parallelDevicePolicy<>, parallelDeviceReduce >( localRhs,
                                                                                rankOffset,
                                                                                dofNumber,
                                                                                elemGhostRank,
                                                                                volume,
                                                                                densOld,
                                                                                poroOld,
                                                                                localResidualNorm );
  } );

  // compute global residual norm
  real64 globalResidualNorm[3] = {0, 0, 0};
  MpiWrapper::allReduce( localResidualNorm,
                         globalResidualNorm,
                         3,
                         MPI_SUM,
                         MPI_COMM_GEOSX );


  real64 const residual = sqrt( globalResidualNorm[0] ) / ( ( globalResidualNorm[1] + m_fluxEstimate ) / (globalResidualNorm[2]+1) );
  return residual;
}


template< typename BASE >
void SinglePhaseFVM< BASE >::applySystemSolution( DofManager const & dofManager,
                                                  arrayView1d< real64 const > const & localSolution,
                                                  real64 const scalingFactor,
                                                  DomainPartition & domain )
{
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  dofManager.addVectorToField( localSolution,
                               BASE::viewKeyStruct::pressureString(),
                               BASE::viewKeyStruct::deltaPressureString(),
                               scalingFactor );

  std::map< string, string_array > fieldNames;
  fieldNames["elems"].emplace_back( string( BASE::viewKeyStruct::deltaPressureString() ) );

  CommunicationTools::getInstance().synchronizeFields( fieldNames, mesh, domain.getNeighbors(), true );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::assembleFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                real64 const dt,
                                                DomainPartition & domain,
                                                DofManager const & dofManager,
                                                CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  fluxApprox.forAllStencils( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    FluxKernel::launch( stencilWrapper,
                        dt,
                        dofManager.rankOffset(),
                        elemDofNumber.toNestedViewConst(),
                        m_elemGhostRank.toNestedViewConst(),
                        m_pressure.toNestedViewConst(),
                        m_deltaPressure.toNestedViewConst(),
                        m_gravCoef.toNestedViewConst(),
                        m_density.toNestedViewConst(),
                        m_dDens_dPres.toNestedViewConst(),
                        m_mobility.toNestedViewConst(),
                        m_dMobility_dPres.toNestedViewConst(),
                        m_permeability.toNestedViewConst(),
                        m_dPerm_dPressure.toNestedViewConst(),
                        m_transTMultiplier.toNestedViewConst(),
                        this->gravityVector(),
                        localMatrix,
                        localRhs );

  } );
}


template< typename BASE >
void SinglePhaseFVM< BASE >::assemblePoroelasticFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs,
                                                           string const & jumpDofKey )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & pressureDofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  pressureDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( pressureDofKey );
  pressureDofNumber.setName( this->getName() + "/accessors/" + pressureDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  jumpDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( jumpDofKey );
  jumpDofNumber.setName( this->getName() + "/accessors/" + jumpDofKey );

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dPerm_dAper =
    mesh.getElemManager().constructMaterialArrayViewAccessor< real64, 3 >( PermeabilityBase::viewKeyStruct::dPerm_dApertureString(),
                                                                           targetRegionNames(),
                                                                           m_permeabilityModelNames,
                                                                           true );

  fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, EmbeddedSurfaceToCellStencil >( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    EmbeddedSurfaceFluxKernel::launch( stencilWrapper,
                                       dt,
                                       dofManager.rankOffset(),
                                       pressureDofNumber.toNestedViewConst(),
                                       jumpDofNumber.toNestedViewConst(),
                                       m_elemGhostRank.toNestedViewConst(),
                                       m_pressure.toNestedViewConst(),
                                       m_deltaPressure.toNestedViewConst(),
                                       m_gravCoef.toNestedViewConst(),
                                       m_density.toNestedViewConst(),
                                       m_dDens_dPres.toNestedViewConst(),
                                       m_mobility.toNestedViewConst(),
                                       m_dMobility_dPres.toNestedViewConst(),
                                       m_permeability.toNestedViewConst(),
                                       m_dPerm_dPressure.toNestedViewConst(),
                                       dPerm_dAper.toNestedViewConst(),
                                       m_transTMultiplier.toNestedViewConst(),
                                       this->gravityVector(),
                                       localMatrix,
                                       localRhs );
  } );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::assemblePoroelasticFluxTerms( real64 const GEOSX_UNUSED_PARAM ( time_n ),
                                                           real64 const dt,
                                                           DomainPartition & domain,
                                                           DofManager const & dofManager,
                                                           CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                           arrayView1d< real64 > const & localRhs,
                                                           CRSMatrixView< real64, localIndex const > const & dR_dAper )
{
  GEOSX_MARK_FUNCTION;

  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );

  NumericalMethodsManager & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  ElementRegionManager::ElementViewAccessor< arrayView3d< real64 const > > dPerm_dAper =
    mesh.getElemManager().constructMaterialArrayViewAccessor< real64, 3 >( PermeabilityBase::viewKeyStruct::dPerm_dApertureString(),
                                                                           targetRegionNames(),
                                                                           m_permeabilityModelNames );

  fluxApprox.forStencils< CellElementStencilTPFA, SurfaceElementStencil, FaceElementToCellStencil >( mesh, [&]( auto & stencil )
  {
    typename TYPEOFREF( stencil ) ::StencilWrapper stencilWrapper = stencil.createStencilWrapper();

    FaceElementFluxKernel::launch( stencilWrapper,
                                   dt,
                                   dofManager.rankOffset(),
                                   elemDofNumber.toNestedViewConst(),
                                   m_elemGhostRank.toNestedViewConst(),
                                   m_pressure.toNestedViewConst(),
                                   m_deltaPressure.toNestedViewConst(),
                                   m_gravCoef.toNestedViewConst(),
                                   m_density.toNestedViewConst(),
                                   m_dDens_dPres.toNestedViewConst(),
                                   m_mobility.toNestedViewConst(),
                                   m_dMobility_dPres.toNestedViewConst(),
                                   m_permeability.toNestedViewConst(),
                                   m_dPerm_dPressure.toNestedViewConst(),
                                   dPerm_dAper.toNestedViewConst(),
                                   m_transTMultiplier.toNestedViewConst(),
                                   this->gravityVector(),
                                   localMatrix,
                                   localRhs,
                                   dR_dAper );
  } );
}

template< typename BASE >
void
SinglePhaseFVM< BASE >::applyBoundaryConditions( real64 const time_n,
                                                 real64 const dt,
                                                 DomainPartition & domain,
                                                 DofManager const & dofManager,
                                                 CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                 arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  BASE::applyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );
  applyFaceDirichletBC( time_n, dt, dofManager, domain, localMatrix, localRhs );
}

template< typename BASE >
void SinglePhaseFVM< BASE >::applyFaceDirichletBC( real64 const time_n,
                                                   real64 const dt,
                                                   DofManager const & dofManager,
                                                   DomainPartition & domain,
                                                   CRSMatrixView< real64, globalIndex const > const & localMatrix,
                                                   arrayView1d< real64 > const & localRhs )
{
  GEOSX_MARK_FUNCTION;

  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  MeshLevel & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  FaceManager & faceManager = mesh.getFaceManager();

  ConstitutiveManager & constitutiveManager = domain.getConstitutiveManager();

  NumericalMethodsManager const & numericalMethodManager = domain.getNumericalMethodManager();
  FiniteVolumeManager const & fvManager = numericalMethodManager.getFiniteVolumeManager();
  FluxApproximationBase const & fluxApprox = fvManager.getFluxApproximation( m_discretizationName );

  // make a list of region indices to be included
  map< localIndex, localIndex > regionFluidMap;
  forTargetRegionsComplete( mesh, [&]( localIndex const targetIndex, localIndex const er, ElementRegionBase & )
  {
    localIndex const modelIndex = constitutiveManager.getSubGroups().getIndex( m_fluidModelNames[targetIndex] );
    regionFluidMap.emplace( er, modelIndex );
  } );

  arrayView1d< real64 const > const presFace =
    faceManager.getReference< array1d< real64 > >( BASE::viewKeyStruct::facePressureString() );

  arrayView1d< real64 const > const gravCoefFace =
    faceManager.getReference< array1d< real64 > >( BASE::viewKeyStruct::gravityCoefString() );

  string const & dofKey = dofManager.getKey( BASE::viewKeyStruct::pressureString() );
  ElementRegionManager::ElementViewAccessor< arrayView1d< globalIndex const > >
  elemDofNumber = mesh.getElemManager().constructArrayViewAccessor< globalIndex, 1 >( dofKey );
  elemDofNumber.setName( this->getName() + "/accessors/" + dofKey );

  // Take BCs defined for "pressure" field and apply values to "facePressure"
  fsManager.apply( time_n + dt,
                   domain,
                   "faceManager",
                   BASE::viewKeyStruct::pressureString(),
                   [&] ( FieldSpecificationBase const & fs,
                         string const & setName,
                         SortedArrayView< localIndex const > const & targetSet,
                         Group & targetGroup,
                         string const & )
  {
    BoundaryStencil const & stencil = fluxApprox.getStencil< BoundaryStencil >( mesh, setName );
    if( stencil.size() == 0 )
    {
      return;
    }

    // first, evaluate BC to get primary field values (pressure)
    fs.applyFieldValue< FieldSpecificationEqual, parallelDevicePolicy<> >( targetSet,
                                                                           time_n + dt,
                                                                           targetGroup,
                                                                           BASE::viewKeyStruct::facePressureString() );

    // Now run the actual kernel
    BoundaryStencil::IndexContainerViewConstType const & seri = stencil.getElementRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sesri = stencil.getElementSubRegionIndices();
    BoundaryStencil::IndexContainerViewConstType const & sefi = stencil.getElementIndices();
    BoundaryStencil::WeightContainerViewConstType const & trans = stencil.getWeights();

    // TODO: currently we just use model from the first cell in this stencil
    //       since it's not clear how to create fluid kernel wrappers for arbitrary models.
    //       Can we just use cell properties for an approximate flux computation?
    //       Then we can forget about capturing the fluid model.
    SingleFluidBase & fluidBase = constitutiveManager.getConstitutiveRelation< SingleFluidBase >( regionFluidMap[seri( 0, 0 )] );

    constitutiveUpdatePassThru( fluidBase, [&]( auto & fluid )
    {
      // create the fluid compute wrapper suitable for capturing in a kernel lambda
      typename TYPEOFREF( fluid ) ::KernelWrapper fluidWrapper = fluid.createKernelWrapper();

      FaceDirichletBCKernel::launch( seri, sesri, sefi, trans,
                                     m_elemGhostRank.toNestedViewConst(),
                                     elemDofNumber.toNestedViewConst(),
                                     dofManager.rankOffset(),
                                     m_pressure.toNestedViewConst(),
                                     m_deltaPressure.toNestedViewConst(),
                                     m_gravCoef.toNestedViewConst(),
                                     m_density.toNestedViewConst(),
                                     m_dDens_dPres.toNestedViewConst(),
                                     m_mobility.toNestedViewConst(),
                                     m_dMobility_dPres.toNestedViewConst(),
                                     presFace,
                                     gravCoefFace,
                                     fluidWrapper,
                                     dt,
                                     localMatrix,
                                     localRhs );
    } );
  } );
}

namespace
{
typedef SinglePhaseFVM< SinglePhaseBase > NoProppant;
typedef SinglePhaseFVM< SinglePhaseProppantBase > Proppant;
REGISTER_CATALOG_ENTRY( SolverBase, NoProppant, string const &, Group * const )
REGISTER_CATALOG_ENTRY( SolverBase, Proppant, string const &, Group * const )
}
} /* namespace geosx */
