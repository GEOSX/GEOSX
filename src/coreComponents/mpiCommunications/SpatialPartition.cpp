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
 * @file SpatialPartition.cpp
 */

#include "mpiCommunications/SpatialPartition.hpp"

#include "codingUtilities/Utilities.hpp"

//#include "Common/intrinsic_typedefs.h"
#include <cmath>

namespace
{
// Modulo
// returns a positive value regardless of the sign of numerator
realT
Mod( realT num, realT denom )
{
  if( fabs( denom ) < fabs( num ) * 1.0e-14 )
  {
    return num;
  }

  return num - denom * std::floor( num / denom );
}

// MapValueToRange
// returns a periodic value in the range [min, max)
realT
MapValueToRange( realT value, realT min, realT max )
{
  return Mod( value - min, max - min ) + min;
}

}  // namespace

namespace geosx
{
using namespace dataRepository;

SpatialPartition::SpatialPartition() :
  PartitionBase(),
  m_Partitions(),
  m_Periodic( nsdof ),
  m_coords( nsdof ),
  m_min( 0.0 ),
  m_max( 0.0 ),
  m_blockSize( 1 ),
  m_gridSize( 0.0 ),
  m_gridMin( 0.0 ),
  m_gridMax( 0.0 )
{
  m_size = 0;
  m_rank = 0;
  m_numColors = 8, setPartitions( 1, 1, 1 );
}

SpatialPartition::~SpatialPartition()
{}

//void SpatialPartition::ReadXML( xmlWrapper::xmlNode const & targetNode )
//{
//  int xpar  = targetNode.attribute("xpar").as_int(1);
//
//}

void
SpatialPartition::InitializePostSubGroups( Group * const )
{
  //get size of problem and decomposition
  m_size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

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
    MpiWrapper::Cart_create( MPI_COMM_GEOSX,
                             nsdof,
                             m_Partitions.data(),
                             m_Periodic.data(),
                             reorder,
                             &cartcomm );
  }
  m_rank = MpiWrapper::Comm_rank( cartcomm );
  MpiWrapper::Cart_coords( cartcomm, m_rank, nsdof, m_coords.data() );

  m_color = GetColor();

  //add neighbors
  {
    int ncoords[nsdof];
    m_neighbors.clear();
    AddNeighbors( 0, cartcomm, ncoords );
  }

  MpiWrapper::Comm_free( cartcomm );

  //initialize cached requests and status
  m_mpiRequest.resize( 2 * m_neighbors.size() );
  m_mpiStatus.resize( 2 * m_neighbors.size() );
}

void
SpatialPartition::InitializeMetis()
{
  //get size of problem and decomposition
  m_size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  m_rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

  //check to make sure our dimensions agree
  GEOSX_ERROR_IF_NE( m_sizeMetis, m_size );

  //initialize cached requests and status
  m_mpiRequest.resize( 100 );
  m_mpiStatus.resize( 100 );
}

int
SpatialPartition::GetColor()
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

void
SpatialPartition::AddNeighbors( const unsigned int idim,
                                MPI_Comm & cartcomm,
                                int * ncoords )
{
  if( idim == nsdof )
  {
    bool me = true;
    for( int i = 0; i < nsdof; i++ )
    {
      if( ncoords[i] != this->m_coords( i ) )
      {
        me = false;
        break;
      }
    }
    if( !me )
    {
      m_neighbors.push_back( NeighborCommunicator() );
      int rank;
      rank = MpiWrapper::Cart_rank( cartcomm, ncoords );
      m_neighbors.back().SetNeighborRank( rank );
      //      m_neighbors.back().Initialize( rank, this->m_rank, this->m_size );

      //      array1d<int> nbrcoords(nsdof);
      //      for(unsigned int i =0; i < nsdof; ++i) nbrcoords[i] = ncoords[i];
      //      neighborCommPtrIndx[nbrcoords] = m_neighbors.size()-1;
    }
  }
  else
  {
    const int dim =
      this->m_Partitions( LvArray::integerConversion< localIndex >( idim ) );
    const bool periodic =
      this->m_Periodic( LvArray::integerConversion< localIndex >( idim ) );
    for( int i = -1; i < 2; i++ )
    {
      ncoords[idim] =
        this->m_coords( LvArray::integerConversion< localIndex >( idim ) ) + i;
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
        AddNeighbors( idim + 1, cartcomm, ncoords );
      }
    }
  }
}

void
SpatialPartition::AddNeighborsMetis( SortedArray< globalIndex > & neighborList )
{
  for( globalIndex const gid : neighborList )
  {
    m_neighbors.push_back( NeighborCommunicator() );
    m_neighbors.back().SetNeighborRank( LvArray::integerConversion< int >( gid ) );

    //    m_neighbors.back().Initialize( LvArray::integerConversion<int>(gid), this->m_rank, this->m_size );
  }
}

//void SpatialPartition::GetPartitionBoundingBox( R1Tensor & xmin, R1Tensor & xmax )
//{
//  xmin = m_xBoundingBoxMin;
//  xmax = m_xBoundingBoxMax;
//}

void
SpatialPartition::setSizes( const R1Tensor & min, const R1Tensor & max )
{
  {
    //get size of problem and decomposition
    m_size = MpiWrapper::Comm_size( MPI_COMM_GEOSX );

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
      MpiWrapper::Cart_create( MPI_COMM_GEOSX,
                               nsdof,
                               m_Partitions.data(),
                               m_Periodic.data(),
                               reorder,
                               &cartcomm );
    }
    m_rank = MpiWrapper::Comm_rank( cartcomm );
    MpiWrapper::Cart_coords( cartcomm, m_rank, nsdof, m_coords.data() );

    m_color = GetColor();

    //add neighbors
    {
      int ncoords[nsdof];
      m_neighbors.clear();
      AddNeighbors( 0, cartcomm, ncoords );
    }

    MpiWrapper::Comm_free( cartcomm );
  }

  // global values
  m_gridMin = min;
  m_gridMax = max;
  m_gridSize = max;
  m_gridSize -= min;

  // block values
  m_blockSize = m_gridSize;

  m_min = min;
  for( int i = 0; i < nsdof; ++i )
  {
    const int nloc = m_Partitions( i ) - 1;
    const localIndex nlocl = static_cast< localIndex >( nloc );
    if( m_PartitionLocations[i].empty() )
    {
      // the default "even" spacing
      m_blockSize( i ) /= m_Partitions( i );
      m_min( i ) += m_coords( i ) * m_blockSize( i );
      m_max( i ) = min( i ) + ( m_coords( i ) + 1 ) * m_blockSize( i );

      m_PartitionLocations[i].resize( nlocl );
      for( localIndex j = 0; j < m_PartitionLocations[i].size(); ++j )
      {
        m_PartitionLocations[i][j] = ( j + 1 ) * m_blockSize( i );
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
        m_min[i] = m_PartitionLocations[i][parIndex - 1];
        m_max[i] = max[i];
      }
      else
      {
        m_min[i] = m_PartitionLocations[i][parIndex - 1];
        m_max[i] = m_PartitionLocations[i][parIndex];
      }
    }
    else
    {
      GEOSX_ERROR(
        "SpatialPartition::setSizes(): number of partition locations does not "
        "equal number of partitions - 1\n" );
    }
  }
}

//void SpatialPartition::setGlobalDomainSizes( const R1Tensor & min, const R1Tensor & max )
//{
//  // global values
//  // without updating partition sizes.  We need this in mesh generator when we
//  // have extension zones.
//  m_gridMin = min;
//  m_gridMax = max;
//  m_gridSize = max;
//  m_gridSize -= min;
//}

void
SpatialPartition::SetPartitionGeometricalBoundary( R1Tensor & min,
                                                   R1Tensor & max )
{
  // We need this in mesh generator when we have extension zones.
  m_min = min;
  m_max = max;
}

bool
SpatialPartition::IsCoordInPartition( const realT & coord, const int dir )
{
  bool rval = true;
  const int i = dir;
  if( m_Periodic( i ) )
  {
    if( m_Partitions( i ) != 1 )
    {
      realT localCenter = MapValueToRange( coord, m_gridMin( i ), m_gridMax( i ) );
      rval = rval && localCenter >= m_min( i ) && localCenter < m_max( i );
    }
  }
  else
  {
    rval =
      rval && ( m_Partitions( i ) == 1 || ( coord >= m_min( i ) && coord < m_max( i ) ) );
  }

  return rval;
}

bool
SpatialPartition::IsCoordInPartition( const R1Tensor & elemCenter )
{
  bool rval = true;
  for( int i = 0; i < nsdof; i++ )
  {
    if( m_Periodic( i ) )
    {
      if( m_Partitions( i ) != 1 )
      {
        realT localCenter =
          MapValueToRange( elemCenter( i ), m_gridMin( i ), m_gridMax( i ) );
        rval = rval && localCenter >= m_min( i ) && localCenter < m_max( i );
      }
    }
    else
    {
      rval = rval &&
        ( m_Partitions( i ) == 1 ||
          ( elemCenter( i ) >= m_min( i ) && elemCenter( i ) < m_max( i ) ) );
    }
  }
  return rval;
}

bool
SpatialPartition::IsCoordInPartition( const R1Tensor & elemCenter,
                                      const int numDistPartition )
{
  bool rval = true;
  R1Tensor m_xBoundingBoxMinTemp, m_xBoundingBoxMaxTemp;
  for( unsigned int i = 0; i < nsdof; i++ )
  {
    m_xBoundingBoxMinTemp( i ) = m_min( i ) - numDistPartition * m_blockSize( i );
    m_xBoundingBoxMaxTemp( i ) = m_max( i ) + numDistPartition * m_blockSize( i );
  }

  for( int i = 0; i < nsdof; i++ )
  {
    if( m_Periodic( i ) )
    {
      if( m_Partitions( i ) != 1 )
      {
        realT localCenter =
          MapValueToRange( elemCenter( i ), m_gridMin( i ), m_gridMax( i ) );
        rval = rval && localCenter >= m_xBoundingBoxMinTemp( i ) &&
          localCenter <= m_xBoundingBoxMaxTemp( i );
      }
    }
    else
    {
      rval = rval &&
        ( m_Partitions( i ) == 1 ||
          ( elemCenter( i ) >= m_xBoundingBoxMinTemp( i ) &&
            elemCenter( i ) <= m_xBoundingBoxMaxTemp( i ) ) );
    }
  }
  return rval;
}

bool
SpatialPartition::IsCoordInPartitionClosed( const R1Tensor & elemCenter )
// A variant with intervals closed at both ends
{
  bool rval = true;
  for( int i = 0; i < nsdof; i++ )
  {
    if( m_Periodic( i ) )
    {
      if( m_Partitions( i ) != 1 )
      {
        realT localCenter =
          MapValueToRange( elemCenter( i ), m_gridMin( i ), m_gridMax( i ) );
        rval = rval && localCenter >= m_min( i ) && localCenter < m_max( i );
      }
    }
    else
    {
      rval = rval &&
        ( m_Partitions( i ) == 1 ||
          ( elemCenter( i ) >= m_min( i ) && elemCenter( i ) <= m_max( i ) ) );
    }
  }
  return rval;
}

bool
SpatialPartition::IsCoordInPartitionBoundingBox( const R1Tensor & elemCenter )

{
  bool rval = true;
  for( int i = 0; i < nsdof; i++ )
  {
    if( m_Periodic( i ) )
    {
      if( m_Partitions( i ) != 1 )
      {
        realT localCenter =
          MapValueToRange( elemCenter( i ), m_gridMin( i ), m_gridMax( i ) );
        rval = rval && localCenter >= m_xBoundingBoxMin( i ) &&
          localCenter <= m_xBoundingBoxMax( i );
      }
    }
    else
    {
      rval = rval &&
        ( m_Partitions( i ) == 1 ||
          ( elemCenter( i ) >= m_xBoundingBoxMin( i ) &&
            elemCenter( i ) <= m_xBoundingBoxMax( i ) ) );
    }
  }
  return rval;
}

void
SpatialPartition::SetContactGhostRange( const realT bufferSize )
{
  m_contactGhostMin = m_min;
  m_contactGhostMin -= bufferSize;

  m_contactGhostMax = m_max;
  m_contactGhostMax += bufferSize;
}

bool
SpatialPartition::IsCoordInContactGhostRange( const R1Tensor & elemCenter )
{
  bool rval = true;
  for( int i = 0; i < nsdof; i++ )
  {
    if( m_Periodic( i ) )
    {
      if( m_Partitions( i ) != 1 )
      {
        realT minBuffer = m_min( i ) - m_contactGhostMin( i );
        realT maxBuffer = m_contactGhostMax( i ) - m_max( i );
        realT localCenterA = MapValueToRange( elemCenter( i ),
                                              m_gridMin( i ) - minBuffer,
                                              m_gridMax( i ) - minBuffer );
        realT localCenterB = MapValueToRange( elemCenter( i ),
                                              m_gridMin( i ) + maxBuffer,
                                              m_gridMax( i ) + maxBuffer );

        rval = rval &&
          ( m_Partitions( i ) == 1 ||
            ( localCenterA >= m_contactGhostMin( i ) &&
              localCenterB < m_contactGhostMax( i ) ) );
      }
    }
    else
    {
      rval = rval &&
        ( m_Partitions( i ) == 1 ||
          ( elemCenter( i ) >= m_contactGhostMin( i ) &&
            elemCenter( i ) < m_contactGhostMax( i ) ) );
    }
  }
  return rval;
}

}  // namespace geosx
