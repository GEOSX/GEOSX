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

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_


#include "PartitionBase.hpp"

#include <map>

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

  ~SpatialPartition() override;

  bool isCoordInPartition( const real64 & coord, const int dir ) override;

  void setSizes( real64 const ( &min )[ 3 ],
                 real64 const ( &max )[ 3 ] ) override;

  void setPartitions( unsigned int xPartitions,
                      unsigned int yPartitions,
                      unsigned int zPartitions ) override;

  int getColor() override;

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

private:

  /**
   * @brief Recursively builds neighbors if an MPI cartesian topology is used (i.e. not metis).
   * @param idim Dimension index in the cartesian.
   * @param cartcomm Communicator with cartesian structure.
   * @param ncoords Cartesian coordinates of a process (assumed to be of length 3).
   *
   * @note Rough copy/paste of DomainPartition::AddNeighbors
   */
  void addNeighbors( const unsigned int idim,
                     MPI_Comm & cartcomm,
                     int * ncoords );

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  void setContactGhostRange( const real64 bufferSize );

  /// Minimum extent of partition dimensions (excluding ghost objects)
  real64 m_min[3];
  /// Maximum extent of partition dimensions (excluding ghost objects)
  real64 m_max[3];

  /// Locations of partition boundaries
  array1d< real64 > m_PartitionLocations[3];

  /// Length of partition dimensions (excluding ghost objects).
  real64 m_blockSize[3];

  /// Total length of problem dimensions (excluding ghost objects).
  real64 m_gridSize[3];
  /// Minimum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMin[3];
  /// Maximum extent of problem dimensions (excluding ghost objects).
  real64 m_gridMax[3];

  /**
   * @brief Ghost position (min).
   */
  real64 m_contactGhostMin[3];

  /**
   * @brief Ghost position (max).
   */
  real64 m_contactGhostMax[3];
};

}
#endif /* GEOSX_MESH_MPICOMMUNICATIONS_SPATIALPARTITION_HPP_ */
