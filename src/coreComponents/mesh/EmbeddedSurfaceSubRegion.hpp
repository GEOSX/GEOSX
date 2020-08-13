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
 * @file EmbeddedSurfaceSubRegion.hpp
 */

#ifndef GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_
#define GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "CellElementSubRegion.hpp"
#include "meshUtilities/SimpleGeometricObjects/BoundedPlane.hpp"

namespace geosx
{
/**
 * @class EmbeddedSurfaceSubRegion
 *
 * The EmbeddedSurfaceSubRegion class contains the functionality to support the concept of an embedded
 * surface element. It consists of a 2D surface that cuts a 3D matrix cell.
 */
class EmbeddedSurfaceSubRegion : public ElementSubRegionBase
{
public:
  /// Embedded surface element to nodes map type
  using NodeMapType =
    InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;

  /// Embedded surface element to edges map type
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

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
  static const string
  CatalogName()
  {
    return "EmbeddedSurfaceSubRegion";
  }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual const string
  getCatalogName() const override
  {
    return EmbeddedSurfaceSubRegion::CatalogName();
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

  virtual void
  CalculateElementGeometricQuantities(
    NodeManager const & nodeManager,
    FaceManager const & facemanager ) override;

  /**
   * @brief Function to compute the geometric quantities of a specific embedded surface element.
   * @param index index of the face element
   */
  void
  CalculateElementGeometricQuantities( localIndex const index );

  /**
   * @brief Function to compute the geometric quantities of a specific embedded surface element.
   * @param intersectionPoints array containing the nodes defining the embedded surface elements
   * @param k index of the face element
   */
  void
  CalculateElementGeometricQuantities( array1d< R1Tensor > const intersectionPoints,
                                       localIndex k );

  /**
   * @brief Function to add a new embedded surface element.
   * @param cellIndex index of the cell element cut by the new embedded surface element
   * @param normalVector unit normal vector to the embedded surface
   */
  void
  AddNewEmbeddedSurface( localIndex const cellIndex, R1Tensor normalVector );

  /**
   * @brief Function to add a new embedded surface element.
   * @param cellIndex cell element index
   * @param regionIndex cell element region index
   * @param subRegionIndex cell element subregion index
   * @param nodeManager the nodemanager group
   * @param edgeManager the edgemanager group
   * @param cellToEdges cellElement to edges map
   * @param fracture pointer to the bounded plane which is defining the embedded surface element
   * @return boolean defining whether the embedded element was added or not
   */
  bool
  AddNewEmbeddedSurface( localIndex const cellIndex,
                         localIndex const regionIndex,
                         localIndex const subRegionIndex,
                         NodeManager const & nodeManager,
                         EdgeManager const & edgeManager,
                         FixedOneToManyRelation const & cellToEdges,
                         BoundedPlane const * fracture );

  /**
   * @brief inherit ghost rank from cell elements.
   * @param cellGhostRank cell element ghost ranks
   */
  void
  inheritGhostRank(
    array1d< array1d< arrayView1d< integer const > > > const & cellGhostRank );

  /**
   * @brief Given the coordinates of a node, it computes the Heaviside function iside a cut element with respect to the fracture element.
   * @param nodeCoord coordinate of the node
   * @param k embedded surface cell index
   * @return value of the Heaviside
   */
  real64
  ComputeHeavisideFunction(
    ArraySlice< real64 const, 1, nodes::REFERENCE_POSITION_USD - 1 > const nodeCoord,
    localIndex const k ) const;

  /**
   * @brief Get list of intersection points (nodes of the embedded plane), offset and connectivity list
   * @param nodeManager the node manager
   * @param edgeManager the edged manager
   * @param elemManager the element region manager
   * @param intersectionPoints coordinates of the intersection points
   * @param connectivityList connectivity list
   * @param offSet offset
   */
  void
  getIntersectionPoints( NodeManager const & nodeManager,
                         EdgeManager const & edgeManager,
                         ElementRegionManager const & elemManager,
                         array1d< R1Tensor > & intersectionPoints,
                         array1d< localIndex > & connectivityList,
                         array1d< int > & offSet ) const;
  /**
   * @brief Compute intersection points (nodes of the embedded plane), offset and connectivity list
   * @param nodeManager the node manager
   * @param edgeManager the edged manager
   * @param elemManager the element region manager
   * @param intersectionPoints coordinates of the intersection points
   * @param connectivityList connectivity list
   * @param offSet offset
   * @param k embedded surface cell index
   */
  void
  ComputeIntersectionPoints( NodeManager const & nodeManager,
                             EdgeManager const & edgeManager,
                             ElementRegionManager const & elemManager,
                             array1d< R1Tensor > & intersectionPoints,
                             array1d< localIndex > & connectivityList,
                             array1d< int > & offSet,
                             localIndex const k ) const;
  ///@}

  /**
   * @brief Struct containing the keys to all embedded surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    /// Embedded surface element aperture string
    static constexpr auto elementApertureString = "elementAperture";

    /// Embedded surface element surface are string
    static constexpr auto elementAreaString = "elementArea";

    /// Embedded surface element cell list string
    static constexpr auto cellListString = "fractureElementsToCellIndices";

    /// Embedded surface element region list string
    static constexpr auto regionListString = "fractureElementsToRegionIndex";

    /// Embedded surface element subregion list string
    static constexpr auto subregionListString =
      "fractureElementsToSubRegionIndex";

    /// Embedded surface element normal vector string
    static constexpr auto normalVectorString = "normalVector";
  };

  virtual void
  setupRelatedObjectsInRelations( MeshLevel const * const mesh ) override;

  virtual string
  GetElementTypeString() const override
  {
    return "Embedded";
  }

  /**
   * @name Relation Accessors
   * @brief Accessor function for the various inter-object relations
   */
  ///@{

  /**
   * @brief Get the embedded surface element to nodes map (background grid nodes).
   * @return the embedded surface element to node map
   */
  NodeMapType const &
  nodeList() const
  {
    return m_toNodesRelation;
  }
  /**
   * @copydoc nodeList() const
   */
  NodeMapType &
  nodeList()
  {
    return m_toNodesRelation;
  }

  /**
   * @brief Get the embedded surface element to region map (background grid nodes).
   * @return the embedded surface element to region map
   */
  arrayView1d< localIndex > const &
  getSurfaceToRegionList()
  {
    return m_embeddedSurfaceToRegion;
  }

  /**
   * @copydoc getSurfaceToRegionList()
   */
  arrayView1d< localIndex const > const &
  getSurfaceToRegionList() const
  {
    return m_embeddedSurfaceToRegion;
  }

  /**
   * @brief Get the embedded surface element to subregion map (of cell elemtns being cut).
   * @return the embedded surface element to subregion map
   */
  arrayView1d< localIndex > const &
  getSurfaceToSubRegionList()
  {
    return m_embeddedSurfaceToSubRegion;
  }

  /**
   * @copydoc getSurfaceToSubRegionList()
   */
  arrayView1d< localIndex const > const &
  getSurfaceToSubRegionList() const
  {
    return m_embeddedSurfaceToSubRegion;
  }

  /**
   * @brief Get the embedded surface element to cell element map
   * @return the embedded surface element to cell element map
   */
  arrayView1d< localIndex > const &
  getSurfaceToCellList()
  {
    return m_embeddedSurfaceToCell;
  }

  /**
   * @copydoc getSurfaceToCellList()
   */
  arrayView1d< localIndex const > const &
  getSurfaceToCellList() const
  {
    return m_embeddedSurfaceToCell;
  }
  ///@}

  /**
   * @name Properties Getters
   * @brief Getters to embedded surface elements properties.
   */
  ///@{

  /**
   * @brief Get number of jump enrichments.
   * @return a reference to the number of jump enrichments
   */
  localIndex &
  numOfJumpEnrichments()
  {
    return m_numOfJumpEnrichments;
  }

  /**
   * @brief Get number of jump enrichments.
   * @return  a constant reference to the number of jump enrichments
   */
  localIndex const &
  numOfJumpEnrichments() const
  {
    return m_numOfJumpEnrichments;
  }

  /**
   * @brief Get face element aperture.
   * @return the aperture of the embedded surface elements
   */
  arrayView1d< real64 > const &
  getElementAperture()
  {
    return m_elementAperture;
  }

  /**
   * @copydoc getElementAperture()
   */
  arrayView1d< real64 const > const &
  getElementAperture() const
  {
    return m_elementAperture;
  }

  /**
   * @brief Get the embedded surface elements surface area.
   * @return the surface area of the embedded surface elements
   */
  arrayView1d< real64 > const &
  getElementArea()
  {
    return m_elementArea;
  }

  /**
   * @copydoc getElementArea()
   */
  arrayView1d< real64 const > const &
  getElementArea() const
  {
    return m_elementArea;
  }

  /**
   * @brief Get normal vectors.
   * @return an array of normal vectors.
   */
  array1d< R1Tensor > &
  getNormalVector()
  {
    return m_normalVector;
  }

  /**
   * @copydoc getNormalVector()
   */
  array1d< R1Tensor > const &
  getNormalVector() const
  {
    return m_normalVector;
  }

  /**
   * @brief Get normal vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the normal vector of a specific embedded surface element
   */
  R1Tensor &
  getNormalVector( localIndex k )
  {
    return m_normalVector[k];
  }

  /**
   * @copydoc getNormalVector( localIndex k )
   */
  R1Tensor const &
  getNormalVector( localIndex k ) const
  {
    return m_normalVector[k];
  }

  /**
   * @brief Get an array of the first tangent vector of the embedded surface elements.
   * @return an array of the first tangent vector of the embedded surface elements
   */
  array1d< R1Tensor > &
  getTangentVector1()
  {
    return m_tangentVector1;
  }

  /**
   * @copydoc getTangentVector1()
   */
  array1d< R1Tensor > const &
  getTangentVector1() const
  {
    return m_tangentVector1;
  }

  /**
   * @brief Get the first tangent vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the first tangent vector of a specific embedded surface element
   */
  R1Tensor &
  getTangentVector1( localIndex k )
  {
    return m_tangentVector1[k];
  }

  /**
   * @copydoc getTangentVector1( localIndex k )
   */
  R1Tensor const &
  getTangentVector1( localIndex k ) const
  {
    return m_tangentVector1[k];
  }

  /**
   * @brief Get an array of the second tangent vector of the embedded surface elements.
   * @return an array of the second tangent vector of the embedded surface elements
   */
  array1d< R1Tensor > &
  getTangentVector2()
  {
    return m_tangentVector2;
  }

  /**
   * @copydoc getTangentVector2()
   */
  array1d< R1Tensor > const &
  getTangentVector2() const
  {
    return m_tangentVector2;
  }

  /**
   * @brief Get the second tangent vector of a specific embedded surface element.
   * @param k index of the embedded surface element
   * @return the second tangent vector of a specific embedded surface element
   */
  R1Tensor &
  getTangentVector2( localIndex k )
  {
    return m_tangentVector2[k];
  }

  /**
   * @copydoc getTangentVector2( localIndex k )
   */
  R1Tensor const &
  getTangentVector2( localIndex k ) const
  {
    return m_tangentVector2[k];
  }
  ///@}

private:
  /// normal vector to the embedded surface element
  array1d< R1Tensor > m_normalVector;

  // tangential direction 1
  array1d< R1Tensor > m_tangentVector1;

  // tangential direction 2
  array1d< R1Tensor > m_tangentVector2;

  /// list of regions
  array1d< localIndex > m_embeddedSurfaceToRegion;

  /// list of subregions
  array1d< localIndex > m_embeddedSurfaceToSubRegion;

  /// list of elements cut by the embedded surface elem
  array1d< localIndex > m_embeddedSurfaceToCell;

  /// list of nodes
  NodeMapType m_toNodesRelation;  // Not used for now. Will need for Flow?

  /// list of edges (if necessary)
  EdgeMapType m_toEdgesRelation;  // Not used for now. Will need for Flow?

  /// The member level field for the element center
  array1d< real64 > m_elementAperture;

  /// The member level field for the element center
  array1d< real64 > m_elementArea;

  localIndex m_numOfJumpEnrichments;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_EMBEDDEDSURFACESUBREGION_HPP_ */
