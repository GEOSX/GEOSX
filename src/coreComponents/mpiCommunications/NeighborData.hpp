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
 * @file NeighborData.hpp
 */

#ifndef GEOSX_MPICOMMUNICATIONS_NEIGHBORDATA_HPP_
#define GEOSX_MPICOMMUNICATIONS_NEIGHBORDATA_HPP_

#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

namespace geosx
{
/**
 * @class NeighborData
 * @brief Class that holds the data associated with a specific neighbor and object class.
 */
class NeighborData : public dataRepository::Group
{
public:
  /**
   * @brief Create a NeighborData object with name @p name and parent @p parent.
   */
  NeighborData( std::string const & name, Group * const parent ) :
    dataRepository::Group( name, parent )
  {
    this
      ->registerWrapper( "matchedPartitionBoundaryObjects",
                         &m_matchedPartitionBoundary )
      ->setSizedFromParent( 0 );

    this->registerWrapper( "ghostsToSend", &m_ghostsToSend )->setSizedFromParent( 0 );

    this->registerWrapper( "ghostsToReceive", &m_ghostsToReceive )
      ->setSizedFromParent( 0 );

    this->registerWrapper( "adjacencyList", &m_adjacencyList )->setSizedFromParent( 0 );
  }

  /**
   * @brief @return True iff there is communication of the associated objects with the neighbor.
   */
  bool
  communicationExists() const
  {
    return m_ghostsToSend.size() && m_ghostsToReceive.size();
  }

  /**
   * @brief @return An array containing the indices of the objects on the domain boundary with the neighbor.
   */
  array1d< localIndex > &
  matchedPartitionBoundary()
  {
    return m_matchedPartitionBoundary;
  }

  /// @copydoc matchedPartitionBoundary()
  arrayView1d< localIndex const > const &
  matchedPartitionBoundary() const
  {
    return m_matchedPartitionBoundary;
  }

  /**
   * @brief @return An array containing the indices of the objects to send to the neighbor.
   */
  array1d< localIndex > &
  ghostsToSend()
  {
    return m_ghostsToSend;
  }

  /// @copydoc ghostsToSend()
  arrayView1d< localIndex const > const &
  ghostsToSend() const
  {
    return m_ghostsToSend;
  }

  /**
   * @brief @return An array containing the indices of the objects to receive from the neighbor.
   */
  array1d< localIndex > &
  ghostsToReceive()
  {
    return m_ghostsToReceive;
  }

  /// @copydoc ghostsToReceive()
  arrayView1d< localIndex const > const &
  ghostsToReceive() const
  {
    return m_ghostsToReceive;
  }

  /**
   * @brief @return An array containing the indices of the adjacent objects (Improve this).
   */
  array1d< localIndex > &
  adjacencyList()
  {
    return m_adjacencyList;
  }

  /// @copydoc adjacencyList()
  arrayView1d< localIndex const > const &
  adjacencyList() const
  {
    return m_adjacencyList;
  }

  /**
   * @brief @return An array containing the globalIndex of any objects the neighbor
   *                mistakenly thinks we own along with the rank that actually does own them.
   */
  array1d< std::pair< globalIndex, int > > &
  nonLocalGhosts()
  {
    return m_nonLocalGhosts;
  }

  /// @copydoc nonLocalGhosts()
  arrayView1d< std::pair< globalIndex, int > const > const &
  nonLocalGhosts() const
  {
    return m_nonLocalGhosts;
  }

private:
  /// Array containing the indices of the objects on the bomain boundary with the neighbor.
  array1d< localIndex > m_matchedPartitionBoundary;

  /// Array containing the indices of the objects to send to the neighbor.
  array1d< localIndex > m_ghostsToSend;

  /// Array containing the indices of the objects to receive from the neighbor.
  array1d< localIndex > m_ghostsToReceive;

  /// Array containing the indices of the adjacent objects (Improve this).
  array1d< localIndex > m_adjacencyList;

  /// An array containing the globalIndex of any objects the neighbor mistakenly
  /// thinks we own along with the rank that actually does own them.
  array1d< std::pair< globalIndex, int > > m_nonLocalGhosts;
};

} /* namespace geosx */

#endif /* GEOSX_MPICOMMUNICATIONS_NEIGHBORDATA_HPP_ */
