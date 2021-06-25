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
 * @file WaveSolverBase.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{

class WaveSolverBase : public SolverBase
{
public:
  WaveSolverBase( const std::string & name,
                  Group * const parent );

  virtual ~WaveSolverBase() override;

  WaveSolverBase() = delete;
  WaveSolverBase( WaveSolverBase const & ) = delete;
  WaveSolverBase( WaveSolverBase && ) = default;

  WaveSolverBase & operator=( WaveSolverBase const & ) = delete;
  WaveSolverBase & operator=( WaveSolverBase && ) = delete;

  virtual void initializePreSubGroups() override;

  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * sourceCoordinatesString() { return "sourceCoordinates"; }

    static constexpr char const * timeSourceFrequencyString() { return "timeSourceFrequency"; }

    static constexpr char const * receiverCoordinatesString() { return "receiverCoordinates"; }

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    static constexpr char const * outputSismoTraceString() { return "outputSismoTrace"; }

  };

protected:

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  virtual void applyFreeSurfaceBC( real64 const time, DomainPartition & domain ) = 0;

  /**
   * @brief Compute the value of a Ricker (a Gaussian function)
   * @param time_n time to evaluate the Ricker
   * @param f0 central frequency of the Ricker
   * @param order order of the ricker
   * @return the value of a Ricker evaluated a time_n with f0
   */
  virtual
  real64 evaluateRicker( real64 const & time_n, real64 const & f0, localIndex order );

  /**
   * @brief Convert a mesh element point coordinate into a coorinate on the reference element
   * @param coords coordinate of the point
   * @param coordsOnRefElem to contain the coordinate computed in the reference element
   * @param indexElement index of the element containing the coords
   * @param faceNodes array of face of the element
   * @param elemsToNodes map to obtaint global nodes from element index
   * @param X array of mesh nodes coordinates
   * @return true if coords is inside the element num index
   */
  template< typename FE_TYPE >
  bool computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                             real64 ( &coordsOnRefElem )[3],
                                             localIndex const & indexElement,
                                             array1d< array1d< localIndex > > const & faceNodes,
                                             arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                             arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X );

  /**
   * @brief Locate sources and receivers position in the mesh elements, evaluate the basis functions at each point and save them to the
   * corresponding elements nodes.
   * @param mesh mesh of the computational domain
   */
  virtual void precomputeSourceAndReceiverTerm( MeshLevel & mesh ) = 0;

  /**
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param time_n the time of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  virtual void addSourceToRightHandSide( real64 const & time_n, arrayView1d< real64 > const rhs ) = 0;

  /**
   * @brief Compute the pressure at each receiver coordinate in one time step
   * @param iseismo index number of the seismo trace
   * @param val_np1 the array to save the value at the receiver position
   */
  virtual void computeSeismoTrace( localIndex const iseismo, arrayView1d< real64 > const pressure_np1 ) = 0;

  /**
   * @brief Save the sismo trace in file
   * @param iseismo index number of the seismo trace
   * @param val_pressure value of the pressure for iseismo
   * @param filename name of the output file
   */
  virtual void saveSeismo( localIndex iseismo, real64 val_pressure, char *filename );

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;

  /// Coordinates of the receivers in the mesh
  array2d< real64 > m_receiverCoordinates;

  /// Flag that indicates the order of the Ricker to be used, order 2 by default
  localIndex m_rickerOrder;

  /// Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise
  localIndex m_outputSismoTrace;



};


template< typename FE_TYPE >
bool WaveSolverBase::computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
                                                           real64 (& coordsOnRefElem)[3],
                                                           localIndex const & indexElement,
                                                           array1d< array1d< localIndex > > const & faceNodes,
                                                           arrayView2d< localIndex const, cells::NODE_MAP_USD > const elemsToNodes,
                                                           arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X )
{
  if( computationalGeometry::IsPointInsidePolyhedron( X, faceNodes, coords ) )
  {
    constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
    real64 xLocal[numNodesPerElem][3];
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      for( localIndex i=0; i<3; ++i )
      {
        xLocal[a][i] = X( elemsToNodes( indexElement, a ), i );
      }
    }

    /// coordsOnRefElem = invJ*(coords-coordsNode_0)
    localIndex q=0;

    real64 invJ[3][3]={{0}};
    FE_TYPE::invJacobianTransformation( q, xLocal, invJ );

    real64 coordsRef[3]={0};
    for( localIndex i=0; i<3; ++i )
    {
      coordsRef[i] = coords[i] - xLocal[q][i];
    }

    for( localIndex i=0; i<3; ++i )
    {
      // Init at (-1,-1,-1) as the origin of the referential elem
      coordsOnRefElem[i] =-1.0;
      for( localIndex j=0; j<3; ++j )
      {
        coordsOnRefElem[i] += invJ[i][j]*coordsRef[j];
      }
    }
    return true;
  }
  else
  {
    return false;
  }
}

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_WAVESOLVERBASE_HPP_ */
