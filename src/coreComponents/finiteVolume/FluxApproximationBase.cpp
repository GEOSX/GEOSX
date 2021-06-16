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
 * @file FluxApproximationBase.cpp
 *
 */

#include "FluxApproximationBase.hpp"

#include "fieldSpecification/FieldSpecificationManager.hpp"
#include "mesh/mpiCommunications/CommunicationTools.hpp"

namespace geosx
{

using namespace dataRepository;

FluxApproximationBase::FluxApproximationBase( string const & name, Group * const parent )
  : Group( name, parent ),
  m_lengthScale( 1.0 )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( viewKeyStruct::fieldNameString(), &m_fieldName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of primary solution field" );

  registerWrapper( viewKeyStruct::coeffNameString(), &m_coeffName ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Name of coefficient field" );

  registerWrapper( viewKeyStruct::targetRegionsString(), &m_targetRegions ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of regions to build the stencil for" );

  registerWrapper( viewKeyStruct::coefficientModelNamesString(), &m_coefficientModelNames ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "List of constitutive models that contain the coefficient used to build the stencil" );

  registerWrapper( viewKeyStruct::areaRelativeToleranceString(), &m_areaRelTol ).
    setInputFlag( InputFlags::OPTIONAL ).
    setApplyDefaultValue( 1.0e-8 ).
    setDescription( "Relative tolerance for area calculations." );
}

FluxApproximationBase::CatalogInterface::CatalogType &
FluxApproximationBase::getCatalog()
{
  static FluxApproximationBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

void FluxApproximationBase::registerDataOnMesh( Group & meshBodies )
{
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();
  meshBodies.forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    meshBody.forSubGroups< MeshLevel >( [&]( MeshLevel & mesh )
    {
      // Group structure: mesh1/finiteVolumeStencils/myTPFA

      Group & stencilParentGroup = mesh.hasGroup( groupKeyStruct::stencilMeshGroupString() ) ?
                                   mesh.getGroup( groupKeyStruct::stencilMeshGroupString() ) :
                                   mesh.registerGroup( groupKeyStruct::stencilMeshGroupString() );

      Group & stencilGroup = stencilParentGroup.hasGroup( getName() ) ?
                             stencilParentGroup.getGroup( getName() ) :
                             stencilParentGroup.registerGroup( getName() );

      registerCellStencil( stencilGroup );
      registerFractureStencil( stencilGroup );
      // For each face-based boundary condition on target field, create a boundary stencil
      fsManager.apply( 0.0,
                       dynamicCast< DomainPartition & >( meshBodies.getParent() ), // TODO: Apply() should take a MeshLevel directly
                       "faceManager",
                       m_fieldName,
                       [&] ( FieldSpecificationBase const &,
                             string const & setName,
                             SortedArrayView< localIndex const > const &,
                             Group const &,
                             string const & )
      {
        registerBoundaryStencil( stencilGroup, setName );
      } );

      FaceManager & faceManager = mesh.getFaceManager();
      faceManager.registerWrapper< array1d< real64 > >( m_coeffName + viewKeyStruct::transMultiplierString() ).
        setApplyDefaultValue( 1.0 ).
        setPlotLevel( PlotLevel::LEVEL_0 ).
        setRegisteringObjects( this->getName() ).
        setDescription( "An array that holds the transmissibility multipliers" );

    } );
  } );
}

void FluxApproximationBase::initializePostInitialConditionsPreSubGroups()
{
  GEOSX_MARK_FUNCTION;

  DomainPartition & domain = this->getGroupByPath< DomainPartition >( "/Problem/domain" );
  FieldSpecificationManager & fsManager = FieldSpecificationManager::getInstance();

  domain.getMeshBodies().forSubGroups< MeshBody >( [&]( MeshBody & meshBody )
  {
    m_lengthScale = meshBody.getGlobalLengthScale();

    meshBody.forSubGroups< MeshLevel >( [&]( MeshLevel & mesh )
    {
      // Compute the main cell-based stencil
      computeCellStencil( mesh );

      // For each face-based boundary condition on target field, create a boundary stencil
      fsManager.apply( 0.0,
                       domain,
                       "faceManager",
                       m_fieldName,
                       [&] ( FieldSpecificationBase const &,
                             string const & setName,
                             SortedArrayView< localIndex const > const & faceSet,
                             Group const &,
                             string const & )
      {
        computeBoundaryStencil( mesh, setName, faceSet );
      } );
    } );
  } );
}

} //namespace geosx
