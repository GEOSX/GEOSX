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
 * @file MetisPartitioner.cpp
 */

#include "MetisPartitioner.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshUtils.hpp"

#include <metis.h>

namespace geosx
{
namespace multiscale
{

MetisPartitioner::MetisPartitioner( LinearSolverParameters::Multiscale::Coarsening params )
  : m_params( std::move( params ) )
{}

static void callMetis( SparsityPatternView< idx_t const, idx_t const > const & graph,
                       LinearSolverParameters::Multiscale::Coarsening::Metis const & params,
                       localIndex const numPart,
                       arrayView1d< localIndex > const & partition )
{
  GEOSX_MARK_FUNCTION;

  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions( options );
  options[METIS_OPTION_SEED] = 2020;
  options[METIS_OPTION_OBJTYPE] = 0;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_UFACTOR] = params.ufactor;

  idx_t nnodes = graph.numRows();
  idx_t nconst = 1;
  idx_t objval = 0;
  idx_t nparts = LvArray::integerConversion< idx_t >( numPart );

  // Must cast away constness to comply with METIS API
  idx_t * offsets = const_cast< idx_t * >( graph.getOffsets() );
  idx_t * columns = const_cast< idx_t * >( graph.getColumns() );

  int result;
  switch( params.method )
  {
    case LinearSolverParameters::Multiscale::Coarsening::Metis::Method::kway:
    {
      result =
        METIS_PartGraphKway( &nnodes, &nconst, offsets, columns, nullptr, nullptr, nullptr,
                             &nparts, nullptr, nullptr, options, &objval, partition.data() );
      break;
    }
    case LinearSolverParameters::Multiscale::Coarsening::Metis::Method::recursive:
    {
      result =
        METIS_PartGraphRecursive( &nnodes, &nconst, offsets, columns, nullptr, nullptr, nullptr,
                                  &nparts, nullptr, nullptr, options, &objval, partition.data() );
      break;
    }
    default:
    {
      GEOSX_THROW( "Invalid METIS partitioning method", std::runtime_error );
    }
  }
  GEOSX_THROW_IF_NE_MSG( result, METIS_OK, "METIS call returned an error", std::runtime_error );
}

localIndex MetisPartitioner::generate( MultiscaleMeshLevel const & mesh,
                                       arrayView1d< localIndex > const & partition )
{
  GEOSX_MARK_FUNCTION;

  localIndex const numCells = mesh.cellManager().numOwnedObjects();
  localIndex const numPart = std::max( localIndex( 1 ), localIndex( real64( numCells ) / m_params.ratio ) ); // rounded down

  if( numPart == 1 )
  {
    partition.zero();
    return numPart;
  }

  MultiscaleMeshObjectManager::MapViewConst const cellToNode = mesh.cellManager().toDualRelation().toViewConst();
  MultiscaleMeshObjectManager::MapViewConst const nodeToCell = mesh.nodeManager().toDualRelation().toViewConst();
  arrayView1d< integer const > const cellGhostRank = mesh.cellManager().ghostRank().toViewConst();

  // Count exact length of each row (METIS does not accept holes in CRS graph)
  array1d< localIndex > rowCounts( numCells );
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView(),
                                           minCommonNodes = m_params.metis.minCommonNodes] ( localIndex const ic )
  {
    integer numUniqueNeighborCells = 0;
    meshUtils::forUniqueNeighbors< 128 >( ic, cellToNode, nodeToCell, [&]( localIndex const nbrIdx, localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != ic && cellGhostRank[nbrIdx] < 0 )
      {
        ++numUniqueNeighborCells;
      }
    } );
    rowCounts[ic] = numUniqueNeighborCells;
  } );

  SparsityPattern< idx_t, idx_t > graph;
  graph.resizeFromRowCapacities< parallelHostPolicy >( numCells, numCells, rowCounts.data() );

  // Fill the graph
  forAll< parallelHostPolicy >( numCells, [=, rowCounts = rowCounts.toView(),
                                           minCommonNodes = m_params.metis.minCommonNodes,
                                           graph = graph.toView()] ( localIndex const ic )
  {
    meshUtils::forUniqueNeighbors< 128 >( ic, cellToNode, nodeToCell, [&]( localIndex const nbrIdx, localIndex const numCommonNodes )
    {
      if( numCommonNodes >= minCommonNodes && nbrIdx != ic && cellGhostRank[nbrIdx] < 0 )
      {
        graph.insertNonZero( LvArray::integerConversion< idx_t >( ic ),
                             LvArray::integerConversion< idx_t >( nbrIdx ) );
      }
    } );
  } );

  callMetis( graph.toViewConst(), m_params.metis, numPart, partition );
  return numPart;
}

} // namespace multiscale
} // namespace geosx
