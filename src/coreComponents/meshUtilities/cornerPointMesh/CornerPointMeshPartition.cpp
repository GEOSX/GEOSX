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
 * @file CornerPointMeshPartition.cpp
 */

#include "CornerPointMeshPartition.hpp"

#include "mpiCommunications/MpiWrapper.hpp"

namespace geosx
{

namespace cornerPointMesh
{

CornerPointMeshPartition::CornerPointMeshPartition( string const & name )
  :
  m_mainRank( 0 ),
  m_myRank( -1 ),
  m_size( 0 ),
  m_meshName( name )
{}

void CornerPointMeshPartition::setupMPIPartition( CornerPointMeshDimensions & dims )
{
  m_myRank = MpiWrapper::commRank( MPI_COMM_GEOSX );
  m_size = MpiWrapper::commSize( MPI_COMM_GEOSX );

  if( m_size == 1 ) // serial run
  {
    dims.definePartitionBoundaries( 0, 0, dims.nX()-1, dims.nY()-1 );
    dims.definePartitionOverlaps( 0, 0, 0, 0 );
  }
  else // parallel run
  {
    broadcastDomainDimensions( dims );
    scatterPartitionBoundaries( dims );
  }
}

void CornerPointMeshPartition::broadcastDomainDimensions( CornerPointMeshDimensions & dims ) const
{
  localIndex nCellsInEachDirection[3];

  // the main rank has read the top of the file and knows the domain dimensions
  if( m_myRank == m_mainRank )
  {
    nCellsInEachDirection[0] = dims.nX();
    nCellsInEachDirection[1] = dims.nY();
    nCellsInEachDirection[2] = dims.nZ();
  }
  localIndex const sendRecvSize = 3;
  MpiWrapper::bcast( nCellsInEachDirection, sendRecvSize, m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();

  // the other ranks have received the data and can set these values in CornerPointMeshData
  if( m_myRank != m_mainRank )
  {
    dims.defineDomainDimensions( nCellsInEachDirection[0],
                                 nCellsInEachDirection[1],
                                 nCellsInEachDirection[2] );
  }
}

void CornerPointMeshPartition::scatterPartitionBoundaries( CornerPointMeshDimensions & dims )
{
  localIndex const sizePerRank = 4;
  array1d< localIndex > boundaries( sizePerRank*m_size );
  array1d< localIndex > overlaps( sizePerRank*m_size );
  array1d< localIndex > neighbors( sizePerRank*m_size );

  // Step 0: the main rank partitions the mesh (temporary partitioning strategy, see below)
  if( m_myRank == m_mainRank )
  {
    // Temporary code to facilitate debugging: 1D partitioning along the x-direction.
    // Ultimately here we will partition in both X and Y, taking into account active cells in each partition.
    // Once I am sure that everything works on large test case, I will rewrite this entirely
    localIndex const spacing = floor( static_cast< real64 >( dims.nX() )
                                      / static_cast< real64 >( m_size ) );
    for( localIndex rank = 0; rank < m_size; ++rank )
    {
      boundaries( rank*sizePerRank )   = rank * spacing;
      boundaries( rank*sizePerRank+1 ) = 0;
      boundaries( rank*sizePerRank+2 ) = ( rank < m_size-1 ) ? ( rank+1 ) * spacing - 1 : dims.nX()-1;
      boundaries( rank*sizePerRank+3 ) = dims.nY()-1;

      overlaps( rank*sizePerRank ) = ( boundaries( rank*sizePerRank ) == 0 ) ? 0 : 1;
      overlaps( rank*sizePerRank+1 ) = ( boundaries( rank*sizePerRank+1 ) == 0 ) ? 0 : 1;
      overlaps( rank*sizePerRank+2 ) = ( boundaries( rank*sizePerRank+2 ) == dims.nX()-1 ) ? 0 : 1;
      overlaps( rank*sizePerRank+3 ) = ( boundaries( rank*sizePerRank+3 ) == dims.nY()-1 ) ? 0 : 1;

      neighbors( rank*sizePerRank )   = ( rank > 0 ) ? rank - 1 : -1;
      neighbors( rank*sizePerRank+1 ) = -1;
      neighbors( rank*sizePerRank+2 ) = ( rank < m_size-1 ) ? rank + 1 : -1;
      neighbors( rank*sizePerRank+3 ) = -1;
    }
  }

  // Step 1: communicate partition boundaries and set them
  array1d< localIndex > myBoundaries( sizePerRank );
  MpiWrapper::scatter( boundaries.toViewConst(), myBoundaries.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  dims.definePartitionBoundaries( myBoundaries( 0 ), myBoundaries( 1 ),
                                  myBoundaries( 2 ), myBoundaries( 3 ) );

  // Step 2: communicate overlap between partitions and set them
  array1d< localIndex > myOverlaps( sizePerRank );
  MpiWrapper::scatter( overlaps.toViewConst(), myOverlaps.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  dims.definePartitionOverlaps( myOverlaps( 0 ), myOverlaps( 1 ),
                                myOverlaps( 2 ), myOverlaps( 3 ) );

  // Step 3: communicate neighbor lists
  array1d< localIndex > myNeighbors( sizePerRank );
  MpiWrapper::scatter( neighbors.toViewConst(), myNeighbors.toView(), m_mainRank, MPI_COMM_GEOSX );
  MpiWrapper::barrier();
  for( localIndex i = 0; i < sizePerRank; ++i )
  {
    if( myNeighbors( i ) != -1 )
    {
      m_neighborsList.insert( myNeighbors( i ) );
    }
  }
}

REGISTER_CATALOG_ENTRY( CornerPointMeshPartition, CornerPointMeshPartition, string const & )

} // namespace cornerPointMesh

} // namespace geosx
