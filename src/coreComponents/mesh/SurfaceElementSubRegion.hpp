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
 * @file SurfaceElementSubRegion.hpp
 */

#ifndef GEOSX_MESH_SURFACEELEMENTSUBREGION_HPP_
#define GEOSX_MESH_SURFACEELEMENTSUBREGION_HPP_

#include "ElementSubRegionBase.hpp"
#include "InterObjectRelation.hpp"
#include "ToElementRelation.hpp"
#include "EdgeManager.hpp"
#include "CellElementSubRegion.hpp"

namespace geosx
{

/**
 * @class SurfaceElementSubRegion
 *
 * The SurfaceElementSubRegion class contains the functionality to support the concept of a
 * surface element that can be either and embedded surface element or a face element.
 */
class SurfaceElementSubRegion : public ElementSubRegionBase
{
public:

  /// Surface element to nodes map type
  using NodeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /// Surface element to edges map type
  using EdgeMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @name Static factory catalog functions
   */
  ///@{

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  static const string catalogName()
  { return "SurfaceElementSubRegion"; }

  /**
   * @brief Get catalog name.
   * @return the catalog name
   */
  virtual const string getCatalogName() const override
  {
    return SurfaceElementSubRegion::catalogName();
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
  SurfaceElementSubRegion( string const & name,
                           dataRepository::Group * const parent );


  /// @brief Destructor
  virtual ~SurfaceElementSubRegion() override;

  ///@}

  /**
   * @name Relation getters
   * @brief Getter functions for the various inter-object relations
   */
  ///@{

  /**
   * @brief Get the face element to nodes map.
   * @return the face element to node map
   */
  NodeMapType const & nodeList() const
  {
    return m_toNodesRelation;
  }

  /**
   * @copydoc nodeList() const
   */
  NodeMapType & nodeList()
  {
    return m_toNodesRelation;
  }

  /**
   * @brief Get the local index of the a-th node of the k-th element.
   * @param[in] k the index of the element
   * @param[in] a the index of the node in the element
   * @return a reference to the local index of the node
   */
  localIndex & nodeList( localIndex const k, localIndex a ) { return m_toNodesRelation( k, a ); }

  /**
   * @copydoc nodeList( localIndex const k, localIndex a )
   */
  localIndex const & nodeList( localIndex const k, localIndex a ) const { return m_toNodesRelation( k, a ); }

  /**
   * @brief Get the surface element to edges map.
   * @return The face element to edge map
   */
  EdgeMapType const & edgeList() const
  {
    return m_toEdgesRelation;
  }

  /**
   * @copydoc edgeList() const
   */
  EdgeMapType & edgeList()
  {
    return m_toEdgesRelation;
  }

  /**
   * @brief Get the surface element to cells map.
   * @return The surface element to cells map
   */
  FixedToManyElementRelation & getToCellRelation()
  {
    return m_surfaceElementsToCells;
  }

  /**
   * @copydoc getToCellRelation()
   */
  FixedToManyElementRelation const & getToCellRelation() const
  {
    return m_surfaceElementsToCells;
  }

  ///@}


  /**
   * @name Properties Getters
   * @brief Getters to surface element properties.
   */
  ///@{

  /**
   * @brief Get face element aperture.
   * @return the aperture of the face elements
   */
  arrayView1d< real64 > getElementAperture() { return m_elementAperture; }

  /**
   * @copydoc getElementAperture()
   */
  arrayView1d< real64 const > getElementAperture() const { return m_elementAperture; }

  /**
   * @brief Get face element surface area.
   * @return the surface area of the face element
   */
  arrayView1d< real64 > getElementArea() { return m_elementArea; }

  /**
   * @copydoc getElementArea()
   */
  arrayView1d< real64 const > getElementArea() const { return m_elementArea; }

  ///@}

  /**
   * @brief Struct containing the keys to all embedded surface element views.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ElementSubRegionBase::viewKeyStruct
  {
    /// @return Face element to cell regions map string.
    static constexpr char const * surfaceElementsToCellRegionsString() { return "fractureElementsToCellRegions"; }

    /// @return Face element to cell subregions map string.
    static constexpr char const * surfaceElementsToCellSubRegionsString() { return "fractureElementsToCellSubRegions"; }

    /// @return Face element to cell indices map string.
    static constexpr char const * surfaceElementsToCellIndexString() { return "fractureElementsToCellIndices"; }


    /// @return Embedded surface element aperture string
    static constexpr char const * elementApertureString() { return "elementAperture"; }

    /// @return Embedded surface element surface are string
    static constexpr char const * elementAreaString() { return "elementArea"; }

    /// @return Mass creation string.
    constexpr static char const * creationMassString() { return "creationMass"; }
  };

protected:

  /// Map between the surface elements and the cells
  FixedToManyElementRelation m_surfaceElementsToCells;

  /// Unmapped surface elements to nodes map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInToNodes;

  /// list of nodes
  NodeMapType m_toNodesRelation;

  /// list of edges
  EdgeMapType m_toEdgesRelation;

  /// Member level field for the element center
  array1d< real64 > m_elementAperture;

  /// Member level field for the element center
  array1d< real64 > m_elementArea;

};

} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_MESH_SURFACEELEMENTSUBREGION_HPP_ */
