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
 * @file EmbeddedSurfaceSubRegion.hpp
 */

#ifndef GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_
#define GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_

#include "SurfaceElementSubRegion.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "EmbeddedSurfaceNodeManager.hpp"
#include "CellElementSubRegion.hpp"
#include "simpleGeometricObjects/BoundedPlane.hpp"

namespace geosx
{

/**
 * @brief Struct defining an embedded element which has at least on node which is a ghost on this rank
 * @struct surfaceWithGhostNodes
 */
struct surfaceWithGhostNodes
{
  /// local index of the surface element
  localIndex surfaceIndex;
  /// index of the parent edge of each node
  std::vector< globalIndex > parentEdgeIndex;
  ///number of nodes of the element
  localIndex numOfNodes;

  /**
   * @brief Constructor
   */
  surfaceWithGhostNodes():
    surfaceIndex(),
    parentEdgeIndex(),
    numOfNodes( 0 ){}
  /**
   * @brief insert a new node
   * @param edgeIndex global index of the parent edge
   */
  void insert ( globalIndex const & edgeIndex );
};

/**
 * @class EmbeddedSurfaceSubRegion
 *
 * The EmbeddedSurfaceSubRegion class contains the functionality to support the concept of an embedded
 * surface element. It consists of a 2D surface that cuts a 3D matrix cell.
 */
class EmbeddedSurfaceSubRegion : public SurfaceElementSubRegion
{
public:

  /// Embedded surface element to faces map type
  using FaceMapType = FixedOneToManyRelation;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  static const string catalogName()
  { return "EmbeddedSurfaceSubRegion"; }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual const string getCatalogName() const override
  {
    return EmbeddedSurfaceSubRegion::catalogName();
  }

  ///@}

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor.
   * @param name the group name
   * @param parent the parent group
   */
  EmbeddedSurfaceSubRegion( string const & name,
                            dataRepository::Group * const parent );

  /// @brief Destructor
  virtual ~EmbeddedSurfaceSubRegion() override;

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  virtual void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                                    FaceManager const & facemanager ) override;

  /**
   * @brief Function to compute the geometric quantities of a specific embedded surface element.
   * @param intersectionPoints array containing the nodes defining the embedded surface elements
   * @param k index of the face element
   */
  void CalculateElementGeometricQuantities( arrayView2d< real64 const > const intersectionPoints,
                                            localIndex k );

  /**
   * @brief Function to add a new embedded surface element.
   * @param cellIndex cell element index
   * @param regionIndex cell element region index
   * @param subRegionIndex cell element subregion index
   * @param nodeManager the nodemanager group
   * @param embSurfNodeManager the embSurfNodeManager group
   * @param edgeManager the edgemanager group
   * @param cellToEdges cellElement to edges map
   * @param fracture pointer to the bounded plane which is defining the embedded surface element
   * @return boolean defining whether the embedded element was added or not
   */
  bool addNewEmbeddedSurface( localIndex const cellIndex,
                              localIndex const regionIndex,
                              localIndex const subRegionIndex,
                              NodeManager const & nodeManager,
                              EmbeddedSurfaceNodeManager & embSurfNodeManager,
                              EdgeManager const & edgeManager,
                              FixedOneToManyRelation const & cellToEdges,
                              BoundedPlane const * fracture );

  /**
   * @brief inherit ghost rank from cell elements.
   * @param cellGhostRank cell element ghost ranks
   */
  void inheritGhostRank( array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank );

  /**
   * @brief Given the coordinates of a node, it computes the Heaviside function iside a cut element with respect to the fracture element.
   * @param nodeCoord coordinate of the node
   * @param k embedded surface cell index
   * @return value of the Heaviside
   */
  real64 computeHeavisideFunction( ArraySlice< real64 const, 1, nodes::REFERENCE_POSITION_USD - 1 > const nodeCoord,
                                   localIndex const k ) const;


  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  ///@}

  /**
   * @brief Struct containing the keys to all embedded surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : SurfaceElementSubRegion::viewKeyStruct
  {
    /// @return Embedded surface element normal vector string
    static constexpr char const * normalVectorString()      { return "normalVector"; }

    /// @return Tangent vector 1 string
    static constexpr char const * t1VectorString()          { return "tangentVector1"; }

    /// @return Tangent vector 2 string
    static constexpr char const * t2VectorString()          { return "tangentVector2"; }

    /// @return Connectivity index string
    static constexpr char const * connectivityIndexString() { return "connectivityIndex"; }

    /// @return Displacement jump string
    static constexpr char const * dispJumpString()          { return "displacementJump"; }

    /// @return Delta displacement jump string
    static constexpr char const * deltaDispJumpString()     { return "deltaDisplacementJump"; }

    /// @return Fracture traction string
    static constexpr char const * fractureTractionString()  { return "fractureTraction"; }

    /// @return Fracture traction derivative w.r.t. jump string
    static constexpr char const * dTraction_dJumpString()   { return "dTraction_dJump"; }

    static constexpr char const * dTraction_dPressureString()   { return "dTraction_dPressure"; }

    /// @return surfaces with ghost nodes list string
    static constexpr char const * surfaceWithGhostNodesString() { return "surfaceWithGhostNodes"; }

    /// Displacement jump key
    dataRepository::ViewKey dispJump        = { dispJumpString() };

    /// Delta displacement jump key
    dataRepository::ViewKey deltaDispJump   = { deltaDispJumpString() };

    /// traction vector key
    dataRepository::ViewKey tractionVector  = { fractureTractionString() };

    /// dTraction_dJump key
    dataRepository::ViewKey dTraction_dJump = { dTraction_dJumpString() };

    /// dTraction_dPressure key
    dataRepository::ViewKey dTraction_dPressure = { dTraction_dPressureString() };

  }
  /// viewKey struct for the EmbeddedSurfaceSubRegion class
  viewKeys;

  virtual void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  virtual string getElementTypeString() const override final { return "Embedded"; }

  /**
   * @name Properties Getters
   * @brief Getters to embedded surface elements properties.
   */
  ///@{

  /**
   * @brief Get number of jump enrichments.
   * @return a reference to the number of jump enrichments
   */
  localIndex & numOfJumpEnrichments()       {return m_numOfJumpEnrichments;}

  /**
   * @brief Get number of jump enrichments.
   * @return  a constant reference to the number of jump enrichments
   */
  localIndex const & numOfJumpEnrichments() const {return m_numOfJumpEnrichments;}

  /**
   * @brief Get normal vectors.
   * @return an array of normal vectors.
   */
  array2d< real64 > & getNormalVector() { return m_normalVector; }

  /**
   * @copydoc getNormalVector()
   */
  arrayView2d< real64 const > getNormalVector() const { return m_normalVector; }

  /**
   * @brief Get normal vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the normal vector of a specific embedded surface element
   */
  arraySlice1d< real64 > getNormalVector( localIndex k ) { return m_normalVector[k]; }

  /**
   * @copydoc getNormalVector( localIndex k )
   */
  arraySlice1d< real64 const > getNormalVector( localIndex k ) const { return m_normalVector[k]; }

  /**
   * @brief Get an array of the first tangent vector of the embedded surface elements.
   * @return an array of the first tangent vector of the embedded surface elements
   */
  array2d< real64 > & getTangentVector1() { return m_tangentVector1; }

  /**
   * @copydoc getTangentVector1()
   */
  arrayView2d< real64 const > getTangentVector1() const { return m_tangentVector1; }

  /**
   * @brief Get the first tangent vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the first tangent vector of a specific embedded surface element
   */
  arraySlice1d< real64 > getTangentVector1( localIndex k ) { return m_tangentVector1[k];}

  /**
   * @copydoc getTangentVector1( localIndex k )
   */
  arraySlice1d< real64 const > getTangentVector1( localIndex k ) const { return m_tangentVector1[k]; }

  /**
   * @brief Get an array of the second tangent vector of the embedded surface elements.
   * @return an array of the second tangent vector of the embedded surface elements
   */
  array2d< real64 > & getTangentVector2() { return m_tangentVector2; }

  /**
   * @copydoc getTangentVector2()
   */
  arrayView2d< real64 const > getTangentVector2() const { return m_tangentVector2; }

  /**
   * @brief Get the second tangent vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the second tangent vector of a specific embedded surface element
   */
  arraySlice1d< real64 > getTangentVector2( localIndex k ) { return m_tangentVector2[k];}

  /**
   * @copydoc getTangentVector2( localIndex k )
   */
  arraySlice1d< real64 const > getTangentVector2( localIndex k ) const { return m_tangentVector2[k];}


  /**
   * @brief Get the connectivity index of the  embedded surface element.
   * @return the connectivity index
   */
  array1d< real64 > & getConnectivityIndex()   { return m_connectivityIndex;}

  /**
   * @copydoc getConnectivityIndex()
   */
  array1d< real64 > const & getConnectivityIndex() const { return m_connectivityIndex;}


  /**
   * @brief Get a mutable displacement jump array.
   * @return the displacement jump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the displacement jump does not exist
   */
  array2d< real64 > & displacementJump()
  { return getReference< array2d< real64 > >( viewKeys.dispJump ); }

  /**
   * @brief Provide an immutable arrayView to the displacement jump array.
   * @return immutable arrayView of the displacement jump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the displacement jump does not exist
   */
  arrayView2d< real64 const > displacementJump() const
  {return getReference< array2d< real64 > >( viewKeys.dispJump ); }

  /**
   * @brief Get a mutable incremental displacement jump array.
   * @return the incremental displacement jump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the incremental displacement jump does not exist
   */
  array2d< real64 > & incrementalDisplacementJump()
  { return getReference< array2d< real64 > >( viewKeys.deltaDispJump ); }

  /**
   * @brief Provide an immutable arrayView to the incremental displacement jump array.
   * @return immutable arrayView of the incremental displacement jump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the incremental displacement jump does not exist
   */
  arrayView2d< real64 const > incrementalDisplacementJump() const
  { return getReference< array2d< real64 > >( viewKeys.deltaDispJump ); }

  /**
   * @brief Get a mutable traction array.
   * @return the traction array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the traction does not exist
   */
  array2d< real64 > & tractionVector()
  { return getReference< array2d< real64 > >( viewKeys.tractionVector ); }

  /**
   * @brief Provide an immutable arrayView to the traction array.
   * @return immutable arrayView of the traction array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the traction does not exist
   */
  arrayView2d< real64 const > tractionVector() const
  {return getReference< array2d< real64 > >( viewKeys.tractionVector ); }

  /**
   * @brief Get a mutable dTraction_dJump array.
   * @return the dTraction_dJump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the dTraction_dJump does not exist
   */
  array3d< real64 > & dTraction_dJump()
  { return getReference< array3d< real64 > >( viewKeys.dTraction_dJump ); }

  /**
   * @brief Provide an immutable arrayView to the dTraction_dJump array.
   * @return immutable arrayView of the dTraction_dJump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the dTraction_dJump does not exist
   */
  arrayView3d< real64 const > dTraction_dJump() const
  { return getReference< array3d< real64 > >( viewKeys.dTraction_dJump ); }

  /**
   * @brief Get a mutable dTraction_dPressure array.
   * @return the dTraction_dJump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the dTraction_dJump does not exist
   */
  array1d< real64 > & dTraction_dPressure()
  { return getReference< array1d< real64 > >( viewKeys.dTraction_dPressure ); }

  /**
   * @brief Provide an immutable arrayView to the dTraction_dJump array.
   * @return immutable arrayView of the dTraction_dJump array if it exists, or an error is thrown if it does not exist
   * @note An error is thrown if the dTraction_dJump does not exist
   */
  arrayView1d< real64 const > dTraction_dPressure() const
  { return getReference< array1d< real64 > >( viewKeys.dTraction_dPressure ); }

  /**
   * @brief accessor to the m_surfaceWithGhostNodes list
   * @return the list of surfaces with at least one ghost node.
   */
  std::vector< struct surfaceWithGhostNodes > surfaceWithGhostNodes() { return m_surfaceWithGhostNodes; }

  ///@}

private:

  /**
   * @brief Pack element-to-node and element-to-face maps
   * @tparam the flag for the bufferOps::Pack function
   * @param buffer the buffer used in the bufferOps::Pack function
   * @param packList the packList used in the bufferOps::Pack function
   * @return the pack size
   */
  template< bool DOPACK >
  localIndex packUpDownMapsPrivate( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList ) const;

  /// normal vector to the embedded surface element
  array2d< real64 > m_normalVector;

  // tangential direction 1
  array2d< real64 > m_tangentVector1;

  // tangential direction 2
  array2d< real64 > m_tangentVector2;

  /// The number of jump enrichments
  localIndex m_numOfJumpEnrichments;

  /// The CI of the cells
  array1d< real64 > m_connectivityIndex;

  /// Surfaces with ghost nodes
  std::vector< struct surfaceWithGhostNodes > m_surfaceWithGhostNodes;
};


} /* namespace geosx */

#endif /* GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_ */
