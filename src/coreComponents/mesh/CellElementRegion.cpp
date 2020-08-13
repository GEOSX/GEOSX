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
 * @file CellElementRegion.cpp
 */

#include "CellElementRegion.hpp"
#include "AggregateElementSubRegion.hpp"
#include "common/TimingMacros.hpp"
#include "LvArray/src/SparsityPattern.hpp"
#include "metis.h"

namespace geosx
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ) :
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceCellBlockNames, &m_cellBlockNames )
    ->setInputFlag( InputFlags::OPTIONAL );

  registerWrapper( viewKeyStruct::coarseningRatioString, &m_coarseningRatio )
    ->setInputFlag( InputFlags::OPTIONAL );
}

CellElementRegion::~CellElementRegion()
{}

void
CellElementRegion::GenerateMesh( Group * const cellBlocks )
{
  Group * const elementSubRegions =
    this->GetGroup( viewKeyStruct::elementSubRegions );

  for( string const & cellBlockName : this->m_cellBlockNames )
  {
    CellElementSubRegion * const subRegion =
      elementSubRegions->RegisterGroup< CellElementSubRegion >( cellBlockName );
    CellBlock * const source =
      cellBlocks->GetGroup< CellBlock >( subRegion->getName() );
    GEOSX_ERROR_IF(
      source == nullptr,
      "Cell block named " + subRegion->getName() + " does not exist" );
    subRegion->CopyFromCellBlock( source );
  }
}

void
CellElementRegion::GenerateAggregates(
  FaceManager const * const faceManager,
  NodeManager const * const GEOSX_UNUSED_PARAM( nodeManager ) )
{
  GEOSX_MARK_FUNCTION;

  if( m_coarseningRatio <= 0. )
  {
    return;
  }
  Group * elementSubRegions = this->GetGroup( viewKeyStruct::elementSubRegions );
  localIndex regionIndex = getIndexInParent();
  AggregateElementSubRegion * const aggregateSubRegion =
    elementSubRegions->RegisterGroup< AggregateElementSubRegion >( "coarse" );

  arrayView2d< localIndex const > const & elemRegionList =
    faceManager->elementRegionList();
  arrayView2d< localIndex const > const & elemSubRegionList =
    faceManager->elementSubRegionList();
  arrayView2d< localIndex const > const & elemList = faceManager->elementList();

  // Counting the total number of cell and number of vertices
  localIndex nbCellElements = 0;
  this->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >(
    [&]( auto & elementSubRegion ) { nbCellElements += elementSubRegion.size(); } );

  // Number of aggregate computation
  localIndex nbAggregates = LvArray::integerConversion< localIndex >(
    int( nbCellElements * m_coarseningRatio ) );
  GEOSX_LOG_RANK_0( "Generating " << nbAggregates << " aggregates on region "
                                  << this->getName() );

  // METIS variable declarations
  using idx_t = ::idx_t;
  idx_t options[METIS_NOPTIONS];       // Contains the METIS options
  METIS_SetDefaultOptions( options );  // ... That are set by default
  idx_t nnodes = LvArray::integerConversion< idx_t >(
    nbCellElements );                // Number of connectivity graph nodes
  idx_t nconst = 1;                  // Number of balancy constraints
  idx_t objval;                      // Total communication volume
  array1d< idx_t > parts( nnodes );  // Map element index -> aggregate index
  idx_t nparts = LvArray::integerConversion< idx_t >(
    nbAggregates );  // Number of aggregates to be generated

  // Compute the connectivity graph
  SparsityPattern< idx_t, idx_t > graph(
    LvArray::integerConversion< idx_t >( nbCellElements ),
    LvArray::integerConversion< idx_t >( nbCellElements ) );
  localIndex nbConnections = 0;
  array1d< localIndex > offsetSubRegions( this->GetSubRegions().size() );
  for( localIndex subRegionIndex = 1; subRegionIndex < offsetSubRegions.size();
       subRegionIndex++ )
  {
    offsetSubRegions[subRegionIndex] = offsetSubRegions[subRegionIndex - 1] +
      this->GetSubRegion( subRegionIndex )->size();
  }
  for( localIndex kf = 0; kf < faceManager->size(); ++kf )
  {
    if( elemRegionList[kf][0] == regionIndex &&
        elemRegionList[kf][1] == regionIndex && elemRegionList[kf][0] )
    {
      localIndex const esr0 = elemSubRegionList[kf][0];
      idx_t const ei0 = LvArray::integerConversion< idx_t >(
        elemList[kf][0] + offsetSubRegions[esr0] );
      localIndex const esr1 = elemSubRegionList[kf][1];
      idx_t const ei1 = LvArray::integerConversion< idx_t >(
        elemList[kf][1] + offsetSubRegions[esr1] );
      graph.insertNonZero( ei0, ei1 );
      graph.insertNonZero( ei1, ei0 );
      nbConnections++;
    }
  }

  // METIS partitionning
  idx_t * offsets = const_cast< idx_t * >( graph.getOffsets() );
  idx_t * columns = const_cast< idx_t * >( &graph.getColumns( 0 )[0] );
  METIS_PartGraphRecursive( &nnodes,
                            &nconst,
                            offsets,
                            columns,
                            nullptr,
                            nullptr,
                            nullptr,
                            &nparts,
                            nullptr,
                            nullptr,
                            options,
                            &objval,
                            parts.data() );

  // Compute Aggregate barycenters
  array1d< R1Tensor > aggregateBarycenters( nparts );
  array1d< real64 > aggregateVolumes( nparts );
  array1d< real64 > normalizeVolumes( nbCellElements );

  // First, compute the volume of each aggregates
  this->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >(
    [&]( ElementSubRegionBase & elementSubRegion ) {
      arrayView1d< integer const > const & ghostRank = elementSubRegion.ghostRank();
      localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
      for( localIndex cellIndex = 0; cellIndex < elementSubRegion.size();
           cellIndex++ )
      {
        if( ghostRank[cellIndex] >= 0 ) continue;
        aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]] +=
          elementSubRegion.getElementVolume()[cellIndex];
      }
    } );

  // Second, compute the normalized volume of each fine elements
  this->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >(
    [&]( ElementSubRegionBase & elementSubRegion ) {
      arrayView1d< integer const > const & ghostRank = elementSubRegion.ghostRank();
      localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
      for( localIndex cellIndex = 0; cellIndex < elementSubRegion.size();
           cellIndex++ )
      {
        if( ghostRank[cellIndex] >= 0 ) continue;
        normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]] =
          elementSubRegion.getElementVolume()[cellIndex] /
          aggregateVolumes[parts[cellIndex + offsetSubRegions[subRegionIndex]]];
      }
    } );

  // Third, normalize the centers
  this->forElementSubRegions< CellElementSubRegion, FaceElementSubRegion >(
    [&]( ElementSubRegionBase & elementSubRegion ) {
      arrayView1d< integer const > const & ghostRank = elementSubRegion.ghostRank();
      localIndex const subRegionIndex = elementSubRegion.getIndexInParent();
      for( localIndex cellIndex = 0; cellIndex < elementSubRegion.size();
           cellIndex++ )
      {
        if( ghostRank[cellIndex] >= 0 ) continue;

        // TODO Change the rest of this to
        // LvArray::tensorOps::scaledAdd< 3 >( aggregateBarycenters[ parts[ cellIndex + offsetSubRegions[ subRegionIndex ]
        // ] ],
        //                                     elementSubRegion.getElementCenter()[ cellIndex ],
        //                                     normalizeVolumes[ cellIndex + offsetSubRegions[ subRegionIndex ] ] )
        real64 const center[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3(
          normalizeVolumes[cellIndex + offsetSubRegions[subRegionIndex]] *
          elementSubRegion.getElementCenter()[cellIndex] );
        aggregateBarycenters[parts[cellIndex + offsetSubRegions[subRegionIndex]]][0] +=
          center[0];
        aggregateBarycenters[parts[cellIndex + offsetSubRegions[subRegionIndex]]][1] +=
          center[1];
        aggregateBarycenters[parts[cellIndex + offsetSubRegions[subRegionIndex]]][2] +=
          center[2];
      }
    } );

  // Convert from metis to GEOSX types
  array1d< localIndex > partsGEOS( parts.size() );
  for( localIndex fineCellIndex = 0; fineCellIndex < partsGEOS.size();
       fineCellIndex++ )
  {
    partsGEOS[fineCellIndex] =
      LvArray::integerConversion< localIndex >( parts[fineCellIndex] );
  }
  aggregateSubRegion->CreateFromFineToCoarseMap( nbAggregates,
                                                 partsGEOS,
                                                 aggregateBarycenters );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        CellElementRegion,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
