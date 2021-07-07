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
 * @file SurfaceGenerator.hpp
 */
#ifndef GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_
#define GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_

#include "mesh/mpiCommunications/NeighborCommunicator.hpp"
#include "physicsSolvers/SolverBase.hpp"
#include "mesh/DomainPartition.hpp"

namespace geosx
{

struct ModifiedObjectLists
{
  std::set< localIndex > newNodes;
  std::set< localIndex > newEdges;
  std::set< localIndex > newFaces;
  std::set< localIndex > modifiedNodes;
  std::set< localIndex > modifiedEdges;
  std::set< localIndex > modifiedFaces;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > newElements;
  map< std::pair< localIndex, localIndex >, std::set< localIndex > > modifiedElements;

  void clearNewFromModified();

  void insert( ModifiedObjectLists const & lists );
};


class SpatialPartition;

class NodeManager;
class EdgeManager;
class FaceManager;
class ElementRegionManager;
class ElementRegionBase;

/**
 * @class SurfaceGenerator
 *
 * This solver manages the mesh topology splitting methods.
 *
 */
class SurfaceGenerator : public SolverBase
{
public:
  SurfaceGenerator( const string & name,
                    Group * const parent );
  ~SurfaceGenerator() override;


  static string catalogName() { return "SurfaceGenerator"; }

  virtual void registerDataOnMesh( Group & MeshBody ) override final;

  /**
   * @defgroup Solver Interface Functions
   *
   * These functions provide the primary interface that is required for derived classes
   */
  /**@{*/

  virtual bool execute( real64 const time_n,
                        real64 const dt,
                        integer const cycleNumber,
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        DomainPartition & domain ) override
  {
    solverStep( time_n, dt, cycleNumber, domain );
    return false;
  }

  virtual real64 solverStep( real64 const & time_n,
                             real64 const & dt,
                             integer const cycleNumber,
                             DomainPartition & domain ) override;

  /**@}*/


  int separationDriver( DomainPartition & domain,
                        MeshLevel & mesh,
                        std::vector< NeighborCommunicator > & neighbors,
                        int const tileColor,
                        int const numTileColors,
                        const bool prefrac,
                        const real64 time_np1 );

  /**
   * @brief Function to generate new global indices of a simple object (node, edge, face)
   * @param[in/out] object A reference to the object that needs new global indices
   * @param[in] indexList the list of local indices that need new global indices
   */
  void assignNewGlobalIndicesSerial( ObjectManagerBase & object,
                                     std::set< localIndex > const & indexList );


  /**
   * @brief Function to generate new global indices for elements
   * @param[in/out] elementManager A reference to the ElementRegionManager that needs new global indices
   * @param[in] indexList the list of local indices that need new global indices
   */
  void
  assignNewGlobalIndicesSerial( ElementRegionManager & elementManager,
                                map< std::pair< localIndex, localIndex >, std::set< localIndex > > const & indexList );

  // SortedArray< localIndex > & getSurfaceElementsRupturedThisSolve() { return m_faceElemsRupturedThisSolve; }

  inline string const getFractureRegionName() const { return m_fractureRegionName; }

protected:

  virtual void initializePostInitialConditionsPreSubGroups() override final;
  virtual void postRestartInitialization() override final;

private:


  /**
   * @brief Function to identify which faces are ready for rupture
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param partition
   * @param prefrac
   */
  void identifyRupturedFaces( DomainPartition & domain,
                              NodeManager & nodeManager,
                              EdgeManager & edgeManager,
                              FaceManager & faceManager,
                              ElementRegionManager & elementManager,
                              const bool prefrac );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param vecTipNorm
   * @param vecTip
   * @return
   */
  real64 calculateEdgeSif ( DomainPartition & domain,
                            const localIndex edgeID,
                            localIndex & trailFaceID,
                            NodeManager & nodeManager,
                            EdgeManager & edgeManager,
                            FaceManager & faceManager,
                            ElementRegionManager & elementManager,
                            real64 ( &vecTipNorm )[ 3 ],
                            real64 ( &vecTip )[ 3 ] );

  /**
   * @brief
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @return
   */
  void calculateNodeAndFaceSif ( DomainPartition & domain,
                                 NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager,
                                 ElementRegionManager & elementManager );

  /**
   * @brief Function to calculate f_disconnect and f_u.
   * @param edgeID
   * @param edgeLength
   * @param nodeIndices
   * @param nodeManager
   * @param edgeManager
   * @param elementManager
   * @param vecTipNorm
   * @param fNode
   * @param GdivBeta
   * @param threeNodesPinched
   * @param calculatef_u. True: calculate f_u; False: calculate f_disconnect.
   */
  int calculateElementForcesOnEdge ( DomainPartition & domain,
                                     const localIndex edgeID,
                                     real64 edgeLength,
                                     localIndex_array & nodeIndices,
                                     NodeManager & nodeManager,
                                     EdgeManager & edgeManager,
                                     ElementRegionManager & elementManager,
                                     real64 ( &vecTipNorm )[3],
                                     real64 ( &fNode )[3],
                                     real64 & GdivBeta,
                                     bool threeNodesPinched,
                                     bool calculatef_u );

  /**
   * @brief
   * @param edgeID
   * @param trailFaceID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param vecTipNorm
   * @param vecTip
   * @param modifiedObjects
   * @param edgeMode
   */
  void markRuptureFaceFromEdge ( const localIndex edgeID,
                                 localIndex & trailFaceID,
                                 NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager,
                                 ElementRegionManager & elementManager,
                                 real64 ( &vecTipNorm )[ 3 ],
                                 real64 ( &vecTip )[ 3 ],
                                 ModifiedObjectLists & modifiedObjects,
                                 const int edgeMode );

  /**
   * @brief
   * @param nodeIndex
   * @param nodeManager
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param modifiedObjects
   */
  void markRuptureFaceFromNode ( const localIndex nodeIndex,
                                 NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager,
                                 ElementRegionManager & elementManager,
                                 ModifiedObjectLists & modifiedObjects );

  /**
   *
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   */
  void postUpdateRuptureStates( NodeManager & nodeManager,
                                EdgeManager & edgeManager,
                                FaceManager & faceManager,
                                ElementRegionManager & elementManager,
                                std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                                std::vector< std::set< localIndex > > & edgesToRupturedFaces );

  /**
   *
   * @param elementManager
   * @param faceManager
   * @param iFace
   */
  int checkOrphanElement( ElementRegionManager & elementManager,
                          FaceManager & faceManager,
                          localIndex iFace );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param faceManager
   * @param edgeManager
   * @param prefrac
   * @return
   */
  int checkEdgeSplitability( const localIndex edgeID,
                             NodeManager & nodeManager,
                             FaceManager & faceManager,
                             EdgeManager & edgeManager,
                             const bool prefrac );


//  void UpdatePathCheckingArrays();

  /**
   * @brief check and split node in mesh
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param elementManager
   * @param modifiedObjects
   * @param prefrac
   * @return
   */
  bool processNode( const localIndex nodeID,
                    real64 const time,
                    NodeManager & nodeManager,
                    EdgeManager & edgeManager,
                    FaceManager & faceManager,
                    ElementRegionManager & elemManager,
                    std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                    std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                    ElementRegionManager & elementManager,
                    ModifiedObjectLists & modifiedObjects,
                    const bool prefrac );

  /**
   * @brief Find a fracture path for surface generation
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elemManager
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool findFracturePlanes( const localIndex nodeID,
                           const NodeManager & nodeManager,
                           const EdgeManager & edgeManager,
                           const FaceManager & faceManager,
                           ElementRegionManager & elemManager,
                           const std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                           const std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                           std::set< localIndex > & separationPathFaces,
                           map< localIndex, int > & edgeLocations,
                           map< localIndex, int > & faceLocations,
                           map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations );


  /**
   * @brief given a fracture path, split the mesh according to the path.
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @param elementManager
   * @param modifiedObjects
   * @param nodesToRupturedFaces
   * @param edgesToRupturedFaces
   * @param separationPathFaces
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   */
  void performFracture( const localIndex nodeID,
                        real64 const time_np1,
                        NodeManager & nodeManager,
                        EdgeManager & edgeManager,
                        FaceManager & faceManager,
                        ElementRegionManager & elementManager,
                        ModifiedObjectLists & modifiedObjects,
                        std::vector< std::set< localIndex > > & nodesToRupturedFaces,
                        std::vector< std::set< localIndex > > & edgesToRupturedFaces,
                        const std::set< localIndex > & separationPathFaces,
                        const map< localIndex, int > & edgeLocations,
                        const map< localIndex, int > & faceLocations,
                        const map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations );

  void mapConsistencyCheck( const localIndex nodeID,
                            NodeManager const & nodeManager,
                            EdgeManager const & edgeManager,
                            FaceManager const & faceManager,
                            ElementRegionManager const & elementManager,
                            const map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodesToElements
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool setLocations( const std::set< localIndex > & separationPathFaces,
                     ElementRegionManager & elemManager,
                     const FaceManager & faceManager,
                     const std::set< std::pair< CellElementSubRegion *, localIndex > > & nodesToElements,
                     const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                     map< localIndex, int > & edgeLocations,
                     map< localIndex, int > & faceLocations,
                     map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations );

  /**
   * @brief function to set which side of the fracture plane all objects are on
   * @param side
   * @param elem
   * @param separationPathFaces
   * @param elemManager
   * @param faceManager
   * @param nodesToElements
   * @param localFacesToEdges
   * @param edgeLocations
   * @param faceLocations
   * @param elemLocations
   * @return
   */
  bool setElemLocations( const int side,
                         const std::pair< CellElementSubRegion *, localIndex > & elem,
                         const std::set< localIndex > & separationPathFaces,
                         ElementRegionManager & elemManager,
                         const FaceManager & faceManager,
                         const std::set< std::pair< CellElementSubRegion *, localIndex > > & nodesToElements,
                         const map< localIndex, std::pair< localIndex, localIndex > > & localFacesToEdges,
                         map< localIndex, int > & edgeLocations,
                         map< localIndex, int > & faceLocations,
                         map< std::pair< CellElementSubRegion *, localIndex >, int > & elemLocations );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 calculateKinkAngle ( const localIndex nodeID,
                              const NodeManager & nodeManager,
                              EdgeManager & edgeManager,
                              FaceManager & faceManager );

  /**
   *
   * @param faceManager
   * @param edgeManager
   * @param nodeManager
   * @param modifiedObjects
   * @param prefrac
   */
  void calculateKinkAngles ( FaceManager & faceManager,
                             EdgeManager & edgeManager,
                             NodeManager & nodeManager,
                             ModifiedObjectLists & modifiedObjects,
                             const bool prefrac );

  /**
   *
   * @param ModifiedObjectLists
   */
  void synchronizeTipSets ( FaceManager & faceManager,
                            EdgeManager & edgeManager,
                            NodeManager & nodeManager,
                            ModifiedObjectLists & receivedObjects );


//  void setDegreeFromCrackTip( NodeManager & nodeManager,
//                              FaceManager & faceManager );

  /**
   *
   * @param edgeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 minimumToughnessOnEdge( const localIndex edgeID,
                                 const NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager );

  /**
   *
   * @param nodeID
   * @param nodeManager
   * @param edgeManager
   * @param faceManager
   * @return
   */
  real64 minimumToughnessOnNode( const localIndex nodeID,
                                 const NodeManager & nodeManager,
                                 EdgeManager & edgeManager,
                                 FaceManager & faceManager );


  real64 calculateRuptureRate( SurfaceElementRegion & faceElementRegion,
                               EdgeManager const & edgeManager );

  /**
   * @struct viewKeyStruct holds char strings and viewKeys for fast lookup
   */
  struct viewKeyStruct : SolverBase::viewKeyStruct
  {
    constexpr static char const * failCriterionString() { return "failCriterion"; }
    constexpr static char const * solidMaterialNameString() { return "solidMaterialNames"; }
    constexpr static char const * fExternalString() { return "fExternal"; }
    constexpr static char const * tipNodesString() { return "tipNodes"; }
    constexpr static char const * tipEdgesString() { return "tipEdges"; }
    constexpr static char const * tipFacesString() { return "tipFaces"; }
    constexpr static char const * trailingFacesString() { return "trailingFaces"; }
    constexpr static char const * fractureRegionNameString() { return "fractureRegion"; }
    constexpr static char const * mpiCommOrderString() { return "mpiCommOrder"; }

    //TODO: rock toughness should be a material parameter, and we need to make rock toughness to KIC a constitutive
    // relation.
    constexpr static char const * rockToughnessString() { return "rockToughness"; }

//    //TODO: Once the node-based SIF criterion becomes mature and robust, remove the edge-based criterion.
    constexpr static char const * nodeBasedSIFString() { return "nodeBasedSIF"; }

  };


private:

  constexpr static real64 m_nonRuptureTime = 1e9;

  /// choice of failure criterion
  integer m_failCriterion=1;

  // solid solver name
  array1d< string > m_solidMaterialNames;

  localIndex m_solidMaterialFullIndex;

  int m_nodeBasedSIF;

  real64 m_rockToughness;

  // Flag for consistent communication ordering
  int m_mpiCommOrder;

  /// set of separable faces
  SortedArray< localIndex > m_separableFaceSet;

  /// copy of the original node->face mapping prior to any separation
  ArrayOfSets< localIndex > m_originalNodetoFaces;

  /// copy of the original node->edge mapping prior to any separation
  ArrayOfSets< localIndex > m_originalNodetoEdges;

  /// copy of the original face->edge mapping prior to any separation
  ArrayOfArrays< localIndex >  m_originalFaceToEdges;

  /// collection of faces that have been used for separation of each node
  array1d< SortedArray< localIndex > > m_usedFacesForNode;

  /// copy of the original face->elemRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemRegion;

  /// copy of the original face->elemSubRegion mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemSubRegion;

  /// copy of the original face->elemIndex mapping prior to any separation
  array2d< localIndex > m_originalFacesToElemIndex;

  /// name of the element region to place all new fractures
  string m_fractureRegionName;

  SortedArray< localIndex > m_tipNodes;

  SortedArray< localIndex > m_tipEdges;

  SortedArray< localIndex > m_tipFaces;

  SortedArray< localIndex > m_trailingFaces;

  SortedArray< localIndex > m_faceElemsRupturedThisSolve;

};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_SURFACEGENERATION_SURFACEGENERATOR_HPP_ */
