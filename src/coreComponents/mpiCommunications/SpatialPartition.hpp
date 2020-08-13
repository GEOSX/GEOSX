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

#ifndef GEOSX_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOSX_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_

#include <map>

#include "mpiCommunications/PartitionBase.hpp"

constexpr int nsdof = 3;
namespace geosx
{
/**
 * @brief Concrete (cartesian?) partitioning.
 */
class SpatialPartition : public PartitionBase
{
public:
  SpatialPartition();
  virtual ~SpatialPartition();

  //  void ReadXML( xmlWrapper::xmlNode const & targetNode );

  /**
 * @brief Initialise Metis partitioning.
 *
 * @note Unused
 */
  virtual void
  InitializeMetis();

  /**
   * @brief Adds some neighbors to the neighbor communicators.
   * @param neighborList The neighbors to add.
   */
  void
  AddNeighborsMetis( SortedArray< globalIndex > & neighborList );
  virtual bool
  IsCoordInPartition( const realT & coord, const int dir );
  virtual bool
  IsCoordInPartition( const R1Tensor & elemCenter );
  virtual bool
  IsCoordInPartition( const R1Tensor & elemCenter,
                      const int numDistPartition );
  /**
   * @brief Variant of IsCoordInPartition with intervals closed at both ends
   * @param elemCenter The point coordinates.
   * @return The predicate result.
   *
   * @note Unused.
   */
  virtual bool
  IsCoordInPartitionClosed( const R1Tensor & elemCenter );

  /**
   * @brief Checks if the point located inside the current partition bouding box.
   * @param elemCenter The element center coordinates.
   * @return The predicate result.
   *
   * @note Unused.
   */
  virtual bool
  IsCoordInPartitionBoundingBox( const R1Tensor & elemCenter );

  virtual bool
  IsCoordInContactGhostRange( const R1Tensor & elemCenter );

  void
  setSizes( const R1Tensor & min, const R1Tensor & max );
  //  void setGlobalDomainSizes( const R1Tensor & min, const R1Tensor & max );
  //  void getSizes( R1Tensor & min, R1Tensor & max ) const;
  //  void getPartitionSizes( R1Tensor & min, R1Tensor & max ) const;
  //  void getPartitionGeometricalBoundary( R1Tensor & min, R1Tensor & max ) const;
  //  void UpdatePartitionBoundingBox(NodeManager& nodeManager);
  //  void GetPartitionBoundingBox( R1Tensor & xmin, R1Tensor & xmax );
  /**
   * @brief Defines the boundaries of the partition
   * @param min The minimum.
   * @param max The maximum.
   */
  void
  SetPartitionGeometricalBoundary( R1Tensor & min, R1Tensor & max );

  void
  setPartitions( unsigned int xPartitions,
                 unsigned int yPartitions,
                 unsigned int zPartitions )
  {
    m_Partitions.resize( 3 );
    m_Partitions( 0 ) = xPartitions;
    m_Partitions( 1 ) = yPartitions;
    m_Partitions( 2 ) = zPartitions;
    m_size = 1;
    for( int i = 0; i < nsdof; i++ ) m_size *= m_Partitions( i );
    SetContactGhostRange( 0.0 );
  }

  /**
   * @brief Defines periodicity along the three (x, y, z) axis. An argument of 1 indicates periodicity.
   * @param xPeriodic Periodicity in x.
   * @param yPeriodic Periodicity in y.
   * @param zPeriodic Periodicity in z.
   */
  void
  setPeriodic( unsigned int xPeriodic,
               unsigned int yPeriodic,
               unsigned int zPeriodic )
  {
    m_Periodic( 0 ) = xPeriodic;
    m_Periodic( 1 ) = yPeriodic;
    m_Periodic( 2 ) = zPeriodic;
  }

  //  void setRadialPeriodic( unsigned int rPeriodic )
  //  {
  //    m_Periodic( 1 ) = rPeriodic;
  //  }

  //  virtual void ReadXMLInput(TICPP::HierarchicalDataNsode& hdn);

  virtual void
  SetContactGhostRange( const realT bufferSize );

  //  virtual void ResetSinglePartitionGlobalToLocalMap(PhysicalDomainT& domain);

  //  virtual void SetPeriodicDomainBoundaryObjects(PhysicalDomainT& domain);
  //  virtual void CorrectReferencePositionsForPeriodicBoundaries(
  //      PhysicalDomainT& domain);

  //  void CreateSinglePartitionGhostObjects(PhysicalDomainT& domain,
  //      const bool contactActive, const int elementGhostingDepth);
  //  void SetSinglePartitionGhostArrays(PhysicalDomainT& domain);

  int
  GetColor();

  //  const array1d< integer > & GetPartitions() const
  //  {
  //    return m_Partitions;
  //  }
  //
  //  const array1d< integer > & GetCoords() const
  //  {
  //    return m_coords;
  //  }
  //
  //  const R1Tensor & xMin() const
  //  {
  //    return m_min;
  //  }
  //  const R1Tensor & xMax() const
  //  {
  //    return m_max;
  //  }
  //
  //  realT xMin( const int i ) const
  //  {
  //    return m_min[i];
  //  }
  //  realT xMax( const int i ) const
  //  {
  //    return m_max[i];
  //  }

protected:
  void
  InitializePostSubGroups( dataRepository::Group * const );

public:
  /// number of partitions
  array1d< int > m_Partitions;
  /**
   * @brief Boolean like array of length 3 (space dimensions).
   *
   * 1 means periodic.
   */
  array1d< int > m_Periodic;
  /// ijk partition indexes
  array1d< int > m_coords;

  /// Minimum extent of partition dimensions (excluding ghost objects)
  R1Tensor m_min;
  /// Maximum extent of partition dimensions (excluding ghost objects)
  R1Tensor m_max;

  /**
   * @brief Bounding box (minimum value) along the x direction.
   * @note Probably unused.
   */
  R1Tensor m_xBoundingBoxMin;
  /**
   * @brief Bounding box (maximum value) along the x direction.
   * @note Probably unused.
   */
  R1Tensor m_xBoundingBoxMax;

  /// Locations of partition boundaries
  array1d< real64 > m_PartitionLocations[3];

  /// Length of partition dimensions (excluding ghost objects).
  R1Tensor m_blockSize;

  /// Total length of problem dimensions (excluding ghost objects).
  R1Tensor m_gridSize;
  /// Minimum extent of problem dimensions (excluding ghost objects).
  R1Tensor m_gridMin;
  /// Maximum extent of problem dimensions (excluding ghost objects).
  R1Tensor m_gridMax;

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void
  AddNeighbors( const unsigned int idim, MPI_Comm & cartcomm, int * ncoords );
};
}  // namespace geosx
#endif /* GEOSX_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
