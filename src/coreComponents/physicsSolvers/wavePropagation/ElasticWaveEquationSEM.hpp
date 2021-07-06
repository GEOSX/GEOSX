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
 * @file ElasticWaveEquationSEM.hpp
 */

#ifndef SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_
#define SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASTICWAVEEQUATIONSEM_HPP_

#include "mesh/ExtrinsicMeshData.hpp"
#include "physicsSolvers/SolverBase.hpp"


namespace geosx
{

class ElasticWaveEquationSEM : public SolverBase
{
public:
  ElasticWaveEquationSEM( const std::string & name,
                           Group * const parent );

  virtual ~ElasticWaveEquationSEM() override;

  ElasticWaveEquationSEM() = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM( ElasticWaveEquationSEM && ) = default;

  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM const & ) = delete;
  ElasticWaveEquationSEM & operator=( ElasticWaveEquationSEM && ) = delete;


  static string catalogName() { return "ElasticSEM"; }

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

    static constexpr char const * rickerOrderString() { return "rickerOrder"; }
    
    static constexpr char const * displacementNp1AtReceiversString() { return "displacementNp1AtReceivers"; }

    static constexpr char const * outputSismoTraceString() { return "outputSismoTrace";}

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

  /// Locates the source term and precomputes the constant part of the source term
  /// And locate receivers and pre_evaluate the basis functions at each receiver coordinate
  void precomputeSourceAndReceiverTerm( MeshLevel & mesh );

  /// Multiply the precomputed term by the ricker and add to the right-hand side
  void addSourceToRightHandSide( real64 const & time, arrayView1d< real64 > const rhs_x, arrayView1d< real64 > const rhs_y, arrayView1d< real64 > const rhs_z);

  /// Apply free surface condition to the face define in the geometry box from the xml
  void applyFreeSurfaceBC( real64 const time, DomainPartition & domain );

  /// Apply absorbing boundary condition to the face define in the geometry box from the xml
  void applyABC( real64 const time, DomainPartition & domain );

  /// Compute the pressure at each receiver coordinate in one time step
  void computeSismoTrace( localIndex const num_timestep, arrayView1d< real64 > const displacementx_np1, arrayView1d< real64 > const displacementy_np1, arrayView1d< real64 > const displacementz_np1 );

  /// save the sismo trace in file
  void saveSismo( localIndex isismo, real64 val_displacement, char *filename );

  /// Coordinates of the sources in the mesh
  array2d< real64 > m_sourceCoordinates;

  /// Indices of the nodes (in the right order) for each source point
  array2d< localIndex > m_sourceNodeIds;

  /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants;

    /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants_x;

    /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants_y;

    /// Constant part of the source for the nodes listed in m_sourceNodeIds
  array2d< real64 > m_sourceConstants_z;

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

  /// Displacement_np1 at the receiver location for each time step for each receiver
  array2d< real64 > m_displacementNp1AtReceivers;

    /// Flag that indicates the order of the Ricker to be used, order 2 by default
  localIndex m_rickerOrder;

  /// Flag that indicates if we write the sismo trace in a file .txt, 0 no output, 1 otherwise
  localIndex m_outputSismoTrace;

};

template< typename FE_TYPE >
bool ElasticWaveEquationSEM::computeCoordinatesOnReferenceElement( real64 const (&coords)[3],
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

EXTRINSIC_MESH_DATA_TRAIT( Displacementx_nm1,
                           "displacementx_nm1",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "x-component of displacement at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementy_nm1,
                           "displacementy_nm1",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "y-component of displacement at time n-1." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementz_nm1,
                           "displacementz_nm1",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "z-component of displacement at time n-1." );


EXTRINSIC_MESH_DATA_TRAIT( Displacementx_n,
                           "displacementx_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "x-component of displacement at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementy_n,
                           "displacementy_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "y-component of displacement at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementz_n,
                           "displacementz_n",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "z-component of displacement at time n." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementx_np1,
                           "displacementx_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "x-component of displacement at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementy_np1,
                           "displacementy_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "y-component of displacement at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( Displacementz_np1,
                           "displacementz_np1",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "z-component of displacement at time n+1." );

EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS_x,
                           "rhs_x",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "RHS" );
                          
EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS_y,
                           "rhs_y",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "RHS" );
          
EXTRINSIC_MESH_DATA_TRAIT( ForcingRHS_z,
                           "rhs_z",
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
                           "Diagonal Mass Matrix." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector_x,
                           "dampingVector_x",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix in x-direction." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector_y,
                           "dampingVector_y",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix in y-direction." );

EXTRINSIC_MESH_DATA_TRAIT( DampingVector_z,
                           "dampingVector_z",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Diagonal Damping Matrix in z-direction." );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocityVp,
                           "mediumVelocityVp",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "P-waves speed in the cell" );

EXTRINSIC_MESH_DATA_TRAIT( MediumVelocityVs,
                           "mediumVelocityVs",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "S-waves speed in the cell" );

EXTRINSIC_MESH_DATA_TRAIT( MediumDensity,
                           "mediumDensity",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Medium density of the cell" );

EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector_x,
                           "stiffnessVector_x",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "x-component of stiffness vector." );
                          
EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector_y,
                           "stiffnessVector_y",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "y-component of stiffness vector." );
          
EXTRINSIC_MESH_DATA_TRAIT( StiffnessVector_z,
                           "stiffnessVector_z",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "z-component of stiffness vector." );

EXTRINSIC_MESH_DATA_TRAIT( LameCoefficientLambda,
                           "lambda",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "First coefficient of Lame." );

EXTRINSIC_MESH_DATA_TRAIT( LameCoefficientMu,
                           "mu",
                           array1d< real64 >,
                           0,
                           NOPLOT,
                           WRITE_AND_READ,
                           "Second coefficient of Lame." );                        

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

#endif /* SRC_CORECOMPONENTS_PHYSICSSOLVERS_WAVEPROPAGATION_ELASSTICWAVEEQUATIONSEM_HPP_ */
