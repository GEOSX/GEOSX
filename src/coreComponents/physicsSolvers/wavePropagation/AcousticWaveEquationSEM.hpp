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
 * @file AcousticWaveEquationSEM.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_
#define GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{

class AcousticWaveEquationSEM : public SolverBase
{
public:

  using EXEC_POLICY = parallelDevicePolicy<32>;
  using OMP_EXEC_POLICY = parallelHostPolicy;
  
  AcousticWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~AcousticWaveEquationSEM() override;

  AcousticWaveEquationSEM() = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM( AcousticWaveEquationSEM && ) = default;

  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM const & ) = delete;
  AcousticWaveEquationSEM & operator=( AcousticWaveEquationSEM && ) = delete;


  static string catalogName() { return "AcousticSEM"; }

  virtual void initializePreSubGroups() override;

  virtual void registerDataOnMesh( Group & meshBodies ) override final;


  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/
  virtual
  real64 solverStep( real64 const & time_n,
                     real64 const & dt,
                     integer const cycleNumber,
                     DomainPartition & domain ) override;

  virtual
  real64 explicitStep( real64 const & time_n,
                       real64 const & dt,
                       integer const cycleNumber,
                       DomainPartition & domain ) override;

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
   * @brief Multiply the precomputed term by the Ricker and add to the right-hand side
   * @param time_n the time of evaluation of the source
   * @param rhs the right hand side vector to be computed
   */
  void addSourceToRightHandSide( real64 const & time_n, arrayView1d< real64 > const rhs );

  /**@}*/


  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    static constexpr char const * sourceCoordinatesString() { return "sourceCoordinates"; }
    static constexpr char const * sourceNodeIdsString() { return "sourceNodeIds"; }
    static constexpr char const * sourceConstantsString() { return "sourceConstants"; }
    static constexpr char const * sourceIsLocalString() { return "sourceIsLocal"; }

    static constexpr char const * timeSourceFrequencyString() { return "timeSourceFrequency"; }

    static constexpr char const * receiverCoordinatesString() { return "receiverCoordinates"; }
    static constexpr char const * receiverNodeIdsString() { return "receiverNodeIds"; }
    static constexpr char const * receiverConstantsString() {return "receiverConstants"; }
    static constexpr char const * receiverIsLocalString() { return "receiverIsLocal"; }

    static constexpr char const * pressureNp1AtReceiversString() { return "pressureNp1AtReceivers"; }

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    static constexpr char const * outputSismoTraceString() { return "outputSismoTrace"; }
    static constexpr char const * dtSismoTraceString() { return "dtSismoTrace"; }
    static constexpr char const * nSampleSismoTraceString() { return "nSampleSismoTrace"; }
    static constexpr char const * indexSismoTraceString() { return "indexSismoTrace"; }


  } waveEquationViewKeys;


protected:

  virtual void postProcessInput() override final;

  virtual void initializePostInitialConditionsPreSubGroups() override final;

private:


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
  void precomputeSourceAndReceiverTerm( MeshLevel & mesh );

  /**
   * @brief Apply free surface condition to the face define in the geometry box from the xml
   * @param time the time to apply the BC
   * @param domain the partition domain
   */
  void applyFreeSurfaceBC( real64 const time, DomainPartition & domain );

  /**
   * @brief Compute the pressure at each receiver coordinate in one time step
   * @param num_timeStep the cycle number of timestep
   * @param pressure_np1 the array to save the pressure value at the receiver position
   */
  void computeSismoTrace( real64 const time_n, real64 const dt, localIndex iSismoTrace, arrayView1d< real64 > const pressure_np1, arrayView1d< real64 > const pressure_n );

  /**
   * @brief Save the sismo trace in file
   * @param pressure_receivers array of pressure values at the receivers locations
   * @param filename name of the output file
   */
  void saveSismo( arrayView2d< real64 const > const pressure_receivers );

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;

  /// Flag that indicates whether the source is local or not to the MPI rank
  array1d< localIndex > m_sourceIsLocal;

  /// Central frequency for the Ricker time source
  real64 m_timeSourceFrequency;

  /// Coordinates of the receivers in the mesh
  array2d< real64 > m_receiverCoordinates;

  /// Indices of the element nodes (in the right order) for each receiver point
  array2d< localIndex > m_receiverNodeIds;

  /// Basis function evaluated at the receiver for the nodes listed in m_receiverNodeIds
  array2d< real64 > m_receiverConstants;

  /// Flag that indicates whether the receiver is local or not to the MPI rank
  array1d< localIndex > m_receiverIsLocal;

  /// Pressure_np1 at the receiver location for each time step for each receiver
  array2d< real64 > m_pressureNp1AtReceivers;


  /// Flag that indicates the order of the Ricker to be used, order 2 by default
  localIndex m_rickerOrder;

  /// Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise
  localIndex m_outputSismoTrace;

  /// Time step size to compute the sismo trace
  real64 m_dtSismoTrace;

  /// Number of sismo trace to be coputed
  localIndex m_nSampleSismoTrace;

  /// Index of the sismo trace
  localIndex m_indexSismoTrace;


};


template< typename FE_TYPE >
bool AcousticWaveEquationSEM::computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
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


namespace extrinsicMeshData
{

EXTRINSIC_MESH_DATA_TRAIT( Pressure_nm1,
                           "pressure_nm1",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Scalar pressure at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_n,
                           "pressure_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Scalar pressure at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Pressure_np1,
                           "pressure_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Scalar pressure at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS,
                           "rhs",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "RHS" );

EXTRINSIC_MESH_DATA_TRAIT( MassVector,
                           "massVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal of the Mass Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector,
                           "dampingVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal of the Damping Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocity,
                           "mediumVelocity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Medium velocity of the cell" );

EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector,
                           "stiffnessVector",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Stiffness vector contains R_h*Pressure_n." );

EXTRINSIC_MESH_DATA_TRAIT( FreeSurfaceFaceIndicator,
                           "freeSurfaceFaceIndicator",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Free surface indicator, 1 if a face is on free surface 0 otherwise." );

EXTRINSIC_MESH_DATA_TRAIT( FreeSurfaceNodeIndicator,
                           "freeSurfaceNodeIndicator",
                           array1d< localIndex >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Free surface indicator, 1 if a node is on free surface 0 otherwise." );


}


} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_WAVEPROPAGATION_ACOUSTICWAVEEQUATIONSEM_HPP_ */
