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

#include "SpatialPartition.hpp"
#include "codingUtilities/Utilities.hpp"
#include "LvArray/src/genericTensorOps.hpp"

#include <cmath>

namespace geosx
{

namespace
{

// Modulo
// returns a positive value regardless of the sign of numerator
real64 Mod( real64 num, real64 denom )
{
  if( fabs( denom )<fabs( num )*1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor( num/denom );
}

// MapValueToRange
// returns a periodic value in the range [min, max)
real64 MapValueToRange( real64 value, real64 min, real64 max )
{
  return Mod( value-min, max-min )+min;
}

}

SpatialPartition::SpatialPartition():
  PartitionBase(),
  m_Partitions(),
  m_Periodic( nsdof ),
  m_coords( nsdof ),
  m_min{ 0.0 },
  m_max{ 0.0 },
  m_blockSize{ 1.0 },
  m_gridSize{ 0.0 },
  m_gridMin{ 0.0 },
  m_gridMax{ 0.0 }
{
  m_size = 0;
  m_rank = 0;
  m_numColors = 8,
  setPartitions( 1, 1, 1 );
}

SpatialPartition::~SpatialPartition()
{}

void SpatialPartition::setPartitions( unsigned int xPartitions,
                                      unsigned int yPartitions,
                                      unsigned int zPartitions )
{
  m_Partitions.resize( 3 );
  m_Partitions( 0 ) = xPartitions;
  m_Partitions( 1 ) = yPartitions;
  m_Partitions( 2 ) = zPartitions;
  m_size = 1;
  for( int i = 0; i < nsdof; i++ )
  {
    m_size *= m_Partitions( i );
  }
  setContactGhostRange( 0.0 );
}

int SpatialPartition::getColor()
{
  int color = 0;

  if( isOdd( m_coords[0] ) )
  {
    color += 1;
  }

  if( isOdd( m_coords[1] ) )
  {
    color += 2;
  }

  if( isOdd( m_coords[2] ) )
  {
    color += 4;
  }

  m_numColors = 8;

  return color;
}

void SpatialPartition::addNeighbors( const unsigned int idim,
                                     MPI_Comm & cartcomm,
                                     int * ncoords )
{

  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != this->m_coords( i ))
      {
        me = false;
        break;
      }
    }
    if( !me )
    {
      int const rank = MpiWrapper::cartRank( cartcomm, ncoords );
      m_neighbors.push_back( NeighborCommunicator( rank ) );
    }
  }
  else
  {
    const int dim = this->m_Partitions( LvArray::integerConversion< localIndex >( idim ) );
    const bool periodic = this->m_Periodic( LvArray::integerConversion< localIndex >( idim ) );
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] = this->m_coords( LvArray::integerConversion< localIndex >( idim ) ) + i;
      bool ok = true;
      if( periodic )
      {
        if( ncoords[idim] < 0 )
          ncoords[idim] = dim - 1;
        else if( ncoords[idim] >= dim )
          ncoords[idim] = 0;
      }
      else
      {
        ok = ncoords[idim] >= 0 && ncoords[idim] < dim;
      }
      if( ok )
      {
        addNeighbors( idim + 1, cartcomm, ncoords );
      }
    }
  }
}

void SpatialPartition::setSizes( real64 const ( &min )[ 3 ],
                                 real64 const ( &max )[ 3 ] )
{

  {
    //get size of problem and decomposition
    m_size = MpiWrapper::commSize( MPI_COMM_GEOSX );

    //check to make sure our dimensions agree
    {
      int check = 1;
      for( int i = 0; i < nsdof; i++ )
      {
        check *= this->m_Partitions( i );
      }
      GEOSX_ERROR_IF_NE( check, m_size );
    }

    //get communicator, rank, and coordinates
    MPI_Comm cartcomm;
    {
      int reorder = 0;
      MpiWrapper::cartCreate( MPI_COMM_GEOSX, nsdof, m_Partitions.data(), m_Periodic.data(), reorder, &cartcomm );
    }
    m_rank = MpiWrapper::commRank( cartcomm );
    MpiWrapper::cartCoords( cartcomm, m_rank, nsdof, m_coords.data());

    //add neighbors
    {
      int ncoords[nsdof];
      m_neighbors.clear();
      addNeighbors( 0, cartcomm, ncoords );
    }

    MpiWrapper::commFree( cartcomm );
  }

  // global values
  LvArray::tensorOps::copy< 3 >( m_gridMin, min );
  LvArray::tensorOps::copy< 3 >( m_gridMax, max );
  LvArray::tensorOps::copy< 3 >( m_gridSize, max );
  LvArray::tensorOps::subtract< 3 >( m_gridSize, min );

  // block values
  LvArray::tensorOps::copy< 3 >( m_blockSize, m_gridSize );

  LvArray::tensorOps::copy< 3 >( m_min, min );
  for( int i = 0; i < nsdof; ++i )
  {
    const int nloc = m_Partitions( i ) - 1;
    const localIndex nlocl = static_cast< localIndex >(nloc);
    if( m_PartitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize[ i ] /= m_Partitions( i );
      m_min[ i ] += m_coords( i ) * m_blockSize[ i ];
      m_max[ i ] = min[ i ] + (m_coords( i ) + 1) * m_blockSize[ i ];

      m_PartitionLocations[i].resize( nlocl );
      for( localIndex j = 0; j < m_PartitionLocations[ i ].size(); ++j )
      {
        m_PartitionLocations[ i ][ j ] = (j+1) * m_blockSize[ i ];
      }
    }
    else if( nlocl == m_PartitionLocations[i].size() )
    {
      const int parIndex = m_coords[i];
      if( parIndex == 0 )
      {
        m_min[i] = min[i];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
      else if( parIndex == nloc )
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_PartitionLocations[i][parIndex-1];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
    }
    else
    {
      GEOSX_ERROR( "SpatialPartition::setSizes(): number of partition locations does not equal number of partitions - 1\n" );
    }
  }
}

bool SpatialPartition::isCoordInPartition( const real64 & coord, const int dir )
{
  bool rval = true;
  const int i = dir;
  if( m_Periodic( i ))
  {
    if( m_Partitions( i ) != 1 )
    {
      real64 localCenter = MapValueToRange( coord, m_gridMin[ i ], m_gridMax[ i ] );
      rval = rval && localCenter >= m_min[ i ] && localCenter < m_max[ i ];
    }

  }
  else
  {
    rval = rval && (m_Partitions[ i ] == 1 || (coord >= m_min[ i ] && coord < m_max[ i ]));
  }

  return rval;
}

void SpatialPartition::setContactGhostRange( const real64 bufferSize )
{
  LvArray::tensorOps::copy< 3 >( m_contactGhostMin, m_min );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMin, -bufferSize );

  LvArray::tensorOps::copy< 3 >( m_contactGhostMax, m_max );
  LvArray::tensorOps::addScalar< 3 >( m_contactGhostMax, bufferSize );
}

}
