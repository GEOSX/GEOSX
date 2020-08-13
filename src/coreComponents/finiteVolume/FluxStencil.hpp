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
 * @file FluxStencil.hpp
 *
 */

#ifndef GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_

#include <common/DataTypes.hpp>
#include "codingUtilities/Utilities.hpp"

namespace geosx
{
/**
 * @class FluxStencil
 *
 * Class representing a spatial discretization stencil.
 */
template< typename INDEX, typename WEIGHT >
class FluxStencil
{
public:
  /**
   * @brief Number of points the flux is between (normally 2).
   */
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /**
   * @brief Maximum number of points in a stencil (required to use static arrays in kernels).
   */
  static localIndex constexpr MAX_STENCIL_SIZE = 18;

  /// Alias for INDEX
  using index_type = INDEX;
  /// Alias for WEIGHT
  using weight_type = WEIGHT;

  explicit FluxStencil();

  /**
   * @brief Constructor.
   * @param[in] numConn number of connections
   * @param[in] avgStencilSize average stencil size
   */
  explicit FluxStencil( localIndex const numConn, localIndex const avgStencilSize );

  /**
   * @brief Return the size of the stencil collection (i.e. number of connections).
   * @return the size of the stencil collection (i.e. number of connections)
   */
  localIndex
  numConnections() const;

  /**
   * @brief Resize the collection.
   * @param[in] numConn number of connections
   * @param[in] avgStencilSize average stencil size
   */
  void
  reserve( localIndex const numConn, localIndex const avgStencilSize );

  /**
   * @brief Add data for one connection.
   * @param[in] numPts number of points to be added
   * @param[in] indices the INDEX array
   * @param[in] weights the WEIGHT array
   * @param[in] connectorIndex the connector index
   */
  void
  add( localIndex const numPts,
       INDEX const * const indices,
       WEIGHT const * const weights,
       localIndex const connectorIndex );

  /**
   * @brief Zero out connections.
   * @param[in] connectorIndex the connector index
   * @return true if the stencil is zeroed out
   */
  bool
  zero( localIndex const connectorIndex );

  /**
   * @brief Called after adding connections is done to compress the data.
   */
  void
  compress();

  /**
   * @struct Entry
   * @brief Structure containing the index and the weight of the single edge in the stencil graph.
   */
  struct Entry
  {
    INDEX index;    ///< edge index
    WEIGHT weight;  ///< edge weight
  };

  /**
   * @brief Return the connections.
   * @return the connections
   */
  ArrayOfArraysView< Entry const, true >
  getConnections() const
  {
    return m_connections.toViewConst();
  }

private:
  ArrayOfArrays< Entry > m_connections;
  map< localIndex, localIndex > m_connectorIndices;
};

template< typename INDEX, typename WEIGHT >
FluxStencil< INDEX, WEIGHT >::FluxStencil() :
  FluxStencil( 0, 0 )
{}

template< typename INDEX, typename WEIGHT >
FluxStencil< INDEX, WEIGHT >::FluxStencil( localIndex const numConn,
                                           localIndex const avgStencilSize ) :
  m_connections()
{
  reserve( numConn, avgStencilSize );
}

template< typename INDEX, typename WEIGHT >
localIndex
FluxStencil< INDEX, WEIGHT >::numConnections() const
{
  return m_connections.size();
}

template< typename INDEX, typename WEIGHT >
void
FluxStencil< INDEX, WEIGHT >::reserve( localIndex const numConn,
                                       localIndex const avgStencilSize )
{
  m_connections.reserve( numConn );
  m_connections.reserveValues( numConn * avgStencilSize );
}

template< typename INDEX, typename WEIGHT >
void
FluxStencil< INDEX, WEIGHT >::add( localIndex const numPts,
                                   INDEX const * const indices,
                                   WEIGHT const * const weights,
                                   localIndex const connectorIndex )
{
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  stackArray1d< Entry, MAX_STENCIL_SIZE > entries( numPts );
  for( localIndex i = 0; i < numPts; ++i )
  {
    entries[i] = { indices[i], weights[i] };
  }

  m_connections.appendArray( entries.begin(), entries.end() );
  m_connectorIndices[connectorIndex] = m_connections.size() - 1;
}

template< typename INDEX, typename WEIGHT >
bool
FluxStencil< INDEX, WEIGHT >::zero( localIndex const connectorIndex )
{
  return executeOnMapValue(
    m_connectorIndices,
    connectorIndex,
    [&]( localIndex const connectionListIndex ) {
      Entry * const entries = m_connections[connectionListIndex];
      for( localIndex i = 0; i < m_connections.sizeOfArray( connectionListIndex );
           ++i )
      {
        entries[i].weight = 0;  // TODO remove entries altogether?
      }
      //    m_connections.resizeArray( connectionListIndex, 0 );
    } );
}

template< typename INDEX, typename WEIGHT >
void
FluxStencil< INDEX, WEIGHT >::compress()
{
  // nothing for the moment
}

}  // namespace geosx

#endif  //GEOSX_FINITEVOLUME_FLUXSTENCIL_HPP_
