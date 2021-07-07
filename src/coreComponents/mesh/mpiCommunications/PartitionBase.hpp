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

#ifndef GEOSX_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_
#define GEOSX_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @brief Base class for partitioning.
 */
class PartitionBase
{
public:

  /**
   * @brief Virtual empty destructor for C++ inheritance reasons
   */
  virtual ~PartitionBase();

  /**
   * @brief Checks if the point located inside the current partition in the given direction dir.
   * @param coord The point coordinates.
   * @param dir The considered direction.
   * @return The predicate result.
   */
  virtual bool isCoordInPartition( const real64 & coord, const int dir ) = 0;

  /**
   * @brief Defines the dimensions of the grid.
   * @param min Global minimum spatial dimensions.
   * @param max Global maximum spatial dimensions.
   */
  virtual void setSizes( real64 const ( &min )[ 3 ],
                         real64 const ( &max )[ 3 ] ) = 0;

  /**
   * @brief Defines the number of partitions along the three (x, y, z) axis.
   * @param xPartitions Number of partitions along x.
   * @param yPartitions Number of partitions along y.
   * @param zPartitions Number of partitions along z.
   */
  virtual void setPartitions( unsigned int xPartitions,
                              unsigned int yPartitions,
                              unsigned int zPartitions ) = 0;

  /**
   * @brief Computes an associated color.
   * @return The color
   *
   * @note The other Color member function.
   */
  virtual int getColor() = 0;

  /**
   * @brief Returns the number of colors.
   * @return The number of associated colors.
   */
  int numColor() const
  { return m_numColors; }

protected:
  /**
   * @brief Preventing dummy default constructor.
   */
  PartitionBase() = default;

  /**
   * @brief Builds from the size of partitions and the current rank of the partition
   * @param numPartitions Size of the partitions.
   * @param thisPartition The rank of the build partition.
   */
  PartitionBase( const unsigned int numPartitions,
                 const unsigned int thisPartition );

  /**
   * @brief Array of neighbor communicators.
   */
  std::vector< NeighborCommunicator > m_neighbors;

  /// Size of the group associated with the MPI communicator
  int m_size;
  /// MPI rank of the current partition
  int m_rank;

  /**
   * @brief Number of colors
   */
  int m_numColors;
};

}

#endif /* GEOSX_MESH_MPICOMMUNICATIONS_PARTITIONBASE_HPP_ */
