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
 * @file PartitionBase.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_
#define GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_

#include "common/DataTypes.hpp"

#include "mpiCommunications/NeighborCommunicator.hpp"

// class oBinStream;
// class iBinStream;

namespace geosx
{
namespace dataRepository
{
class Group;
}

class DomainPartition;
class ObjectManagerBase;

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
   * @brief Defines the domain
   * @param domain The new domain
   *
   * Actually unused function, consider removing.
   */
  void
  SetDomain( DomainPartition * domain );

  /**
   * @brief Checks if the point located inside the current partition.
   * @param elemCenter The point coordinates.
   * @return The predicate result.
   */
  virtual bool
  IsCoordInPartition( const R1Tensor & elemCenter ) = 0;
  /**
   * @brief Checks if the point located inside the current partition, taking in to account the distant partition.
   * @param elemCenter The point coordinates.
   * @param numDistPartition The number of distant partitions.
   * @return The predicate result.
   */
  virtual bool
  IsCoordInPartition( const R1Tensor & elemCenter,
                      const int numDistPartition ) = 0;
  /**
   * @brief Checks if the point located inside the current partition in the given direction dir.
   * @param coord The point coordinates.
   * @param dir The considered direction.
   * @return The predicate result.
   */
  virtual bool
  IsCoordInPartition( const realT & coord, const int dir ) = 0;

  /**
   * @brief Defines the dimensions of the grid.
   * @param min Global minimum spatial dimensions.
   * @param max Global maximum spatial dimensions.
   */
  virtual void
  setSizes( const R1Tensor & min, const R1Tensor & max ) = 0;

  /**
   * @brief Defines the number of partitions along the three (x, y, z) axis.
   * @param xPartitions Number of partitions along x.
   * @param yPartitions Number of partitions along y.
   * @param zPartitions Number of partitions along z.
   */
  virtual void
  setPartitions( unsigned int xPartitions,
                 unsigned int yPartitions,
                 unsigned int zPartitions ) = 0;

  /**
   * @brief Checks if the point (as an element center) is in contact of a ghost.
   * @param elemCenter The position of an element center.
   * @return The predicate result.
   */
  virtual bool
  IsCoordInContactGhostRange( const R1Tensor & elemCenter ) = 0;

  //  virtual void ReadXML( xmlWrapper::xmlNode const & targetNode ) = 0;

  //virtual void AssignGlobalIndices( DomainPartition * domain );

  //  virtual void FindMatchedBoundaryIndices( string const & key,
  //                                           const ObjectManagerBase& object );

  //  virtual void SetUpNeighborLists( DomainPartition * domain,
  //                                   const bool contactActive );

  //  void SetRankOfNeighborNeighbors();

  //  virtual void ResetNeighborLists( PhysicalDomainT& domain,
  //                                   const int elementGhostingDepth );

  //  virtual void ModifyGhostsAndNeighborLists( const ModifiedObjectLists& modifiedObjects );

  //  template< typename T >
  //  void SendReceive( const array1d< array1d< T > > & sendArray, array1d< array1d< T > > & recvArray );

  //  void SynchronizeFields( const std::map<std::string, string_array >& fieldNames,
  //                          const CommRegistry::commID commID = CommRegistry::genericComm01 );

  //  void SetOwnedByRank( const std::map< std::string, globalIndex_array > & localBoundaryGlobalIndices,
  //                       std::map< std::string, std::map< globalIndex, int > > & boundaryOwnership );

  //  void SetGhostArrays( DomainPartition * domain );

  //  localIndex_array GetFaceSendIndices();

  /**
   * @brief Defines a distance/buffer below which we are considered in the contact zone ghosts.
   * @param bufferSize The distance.
   */
  virtual void
  SetContactGhostRange( const realT bufferSize ) = 0;
  //
  //  void SetBufferSizes( const std::map<string, string_array >& fieldNames,
  //                       const CommRegistry::commID commID  );
  //
  //  int NumberOfNeighbors( ) {return LvArray::integerConversion<int>(m_neighbors.size());}

  /// Size of the group associated with the MPI communicator
  int m_size;
  /// Metis size, unused
  int m_sizeMetis;
  /// MPI rank of the current partition
  int m_rank;

  /**
   * @brief Computes an associated color.
   * @return The color
   *
   * @note The other Color member function.
   */
  virtual int
  GetColor() = 0;

  /**
   * @brief Returns the associated color.
   * @return The color.
   *
   * @note The other GetColor member function.
   */
  int
  Color() const
  {
    return m_color;
  }
  /**
   * @brief Returns the number of colors.
   * @return The number of associated colors.
   */
  int
  NumColor() const
  {
    return m_numColors;
  }

  //  void DeleteExcessNeighbors();
  //  void GraphBasedColoring();

protected:
  /**
   * @brief Preventing dummy default constructor.
   */
  PartitionBase();
  /**
   * @brief Builds from the size of partitions and the current rank of the partition
   * @param numPartitions Size of the partitions.
   * @param thisPartiton The rank of the build partition.
   */
  PartitionBase( const unsigned int numPartitions,
                 const unsigned int thisPartiton );

  /**
   * @brief Called by Initialize() after to initializing sub-Groups.
   * @param group A group that is passed in to the initialization functions
   *              in order to facilitate the initialization.
   */
  virtual void
  InitializePostSubGroups( dataRepository::Group * const group ) = 0;

  /**
   * @brief Array of neighbor communicators.
   */
  std::vector< NeighborCommunicator > m_neighbors;

  /**
   * @brief Array of mpi_requests
   */
  array1d< MPI_Request > m_mpiRequest;
  /**
   * @brief Array of mpi statuses
   */
  array1d< MPI_Status > m_mpiStatus;

  /**
   * @brief Ghost position (min).
   */
  R1Tensor m_contactGhostMin;
  /**
   * @brief Ghost position (max).
   */
  R1Tensor m_contactGhostMax;

  /**
   * @brief Associated color
   */
  int m_color;
  /**
   * @brief Number of colors
   */
  int m_numColors;

  /**
   * @brief Reference to the associated domain.
   */
  DomainPartition * const m_domain;

public:
  /// Unused parameter
  realT m_t1;
  /// Unused parameter
  realT m_t2;
  /// Unused parameter
  realT m_t3;
  /// Unused parameter
  realT m_t4;
  /// Unused parameter
  bool m_hasLocalGhosts;
  /// Unused parameter
  std::map< std::string, localIndex_array > m_localGhosts;
  /// Unused parameter
  std::map< std::string, localIndex_array > m_elementRegionsLocalGhosts;
  /// Unused parameter
  std::map< std::string, localIndex_array > m_localGhostSources;
  /// Unused parameter
  std::map< std::string, localIndex_array > m_elementRegionsLocalGhostSources;
  /// Unused parameter
  int m_ghostDepth;

private:
  //  virtual void AssignGlobalIndices( ObjectDataStructureBaseT& object, const ObjectDataStructureBaseT&
  // compositionObject );

  //  void CommunicateRequiredObjectIndices();
};

}  // namespace geosx

#endif /* GEOSX_MPICOMMUNICATIONS_PARTITIONBASE_HPP_ */
