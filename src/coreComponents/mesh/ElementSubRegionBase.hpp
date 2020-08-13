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
 * @file ElementSubRegionBase.hpp
 */

#ifndef GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_
#define GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_

#include "managers/ObjectManagerBase.hpp"
namespace geosx
{
class NodeManager;
class FaceManager;
class MeshLevel;
class DomainPartition;
namespace constitutive
{
class ConstitutiveBase;
}

/**
 * @class ElementSubRegionBase
 * Abstract class for a collection of mesh elements that
 * will be derived and specialized for cell elements,
 * face (fracture) elements, embedded surface (fracture)
 * elements, well elements, etc
 */
class ElementSubRegionBase : public ObjectManagerBase
{
public:
  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  ElementSubRegionBase( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Destructor.
   */
  ~ElementSubRegionBase();

  ///@}

  /**
   * @name Geometry computation / Connectivity
   */
  ///@{

  /**
   * @brief Calculate the geometric quantities for each element in the subregion.
   * @param[in] nodeManager the nodeManager (for geometrical info and connectivity involving nodes)
   * @param[in] faceManager the faceManager (for geometrical info and connectivity involving faces)
   *
   */
  virtual void
  CalculateElementGeometricQuantities(
    NodeManager const & nodeManager,
    FaceManager const & faceManager ) = 0;

  /**
   * @brief Link the connectivity maps of the subregion to the managers storing the mesh information.
   * @param[in] mesh the meshLevel object (single level only)
   *
   * In the derived classes, this function is used to passe a pointer to the nodeManager,
   * faceManager, and (if needed) edgeManager to, respectively, the node list, face list, and
   * edge list of the subregion.
   */
  virtual void
  setupRelatedObjectsInRelations( MeshLevel const * const mesh ) = 0;

  /**
   * @brief Call ObjectManagerBase::FixUpDownMaps for the connectivity maps needed by
   *        the derived class (i.e., element-to-node map, element-to-face map, etc)
   * @param[in] clearIfUnmapped clearIfUnmapped
   */
  virtual void
  FixUpDownMaps( bool const clearIfUnmapped )
  {
    GEOSX_UNUSED_VAR( clearIfUnmapped );
  }

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the number of nodes per element.
   * @return number of nodes per element
   */
  localIndex const &
  numNodesPerElement() const
  {
    return m_numNodesPerElement;
  }

  /**
   * @brief Get the number of nodes per element.
   * @param[in] k cell index (not used)
   * @return number of nodes per element
   */
  virtual localIndex
  numNodesPerElement( localIndex const k ) const
  {
    GEOSX_UNUSED_VAR( k );
    return m_numNodesPerElement;
  }

  /**
   * @brief Set the number of nodes per element.
   * @param[in] numNodes the number of nodes per element
   */
  void
  setNumNodesPerElement( localIndex numNodes )
  {
    m_numNodesPerElement = numNodes;
  }

  /**
   * @brief Get the number of independent nodes per element.
   * @return numNodes the number of independent nodes per element
   *
   * Currently, the number of independent nodes per element is always
   * equal to the number of nodes per element, except for the case
   * of a triangular prism
   */
  localIndex const &
  numIndependentNodesPerElement() const
  {
    return m_numIndependentNodesPerElement;
  }

  /**
   * @brief Set the number of independent nodes per element.
   * @param[in] numNodes the number of independent nodes per element
   */
  void
  setNumIndependentNodesPerElement( localIndex const numNodes )
  {
    m_numIndependentNodesPerElement = numNodes;
  }

  /**
   * @brief Get the number of edges per element.
   * @return the number of edges per element
   */
  localIndex const &
  numEdgesPerElement() const
  {
    return m_numEdgesPerElement;
  }

  /**
   * @brief Set the number of edges per element.
   * @param[in] numEdges the number of edges per element
   */
  void
  setNumEdgesPerElement( localIndex const numEdges )
  {
    m_numEdgesPerElement = numEdges;
  }

  /**
   * @brief Get the number of faces per element.
   * @return number of faces per element
   */
  localIndex const &
  numFacesPerElement() const
  {
    return m_numFacesPerElement;
  }

  /**
   * @brief Set the number of faces per element.
   * @param[in] numFaces the number of faces per element
   */
  void
  setNumFacesPerElement( localIndex const numFaces )
  {
    m_numFacesPerElement = numFaces;
  }

  /**
   * @brief Get the center of each element in this subregion.
   * @return an arrayView1d of const element centers
   */
  arrayView2d< real64 const > const &
  getElementCenter() const
  {
    return m_elementCenter;
  }

  /**
   * @copydoc getElementCenter() const
   */
  arrayView2d< real64 > const &
  getElementCenter()
  {
    return m_elementCenter;
  }

  /**
   * @brief Get the volume of each element in this subregion.
   * @return an arrayView1d of const element volumes
   */
  arrayView1d< real64 const > const &
  getElementVolume() const
  {
    return m_elementVolume;
  }

  /**
   * @brief Get the group in which the constitutive models of this subregion are registered.
   * @return a pointer to the const group in which the constitutive models are registered
   */
  dataRepository::Group const *
  GetConstitutiveModels() const
  {
    return &m_constitutiveModels;
  }

  /**
   * @copydoc GetConstitutiveModels() const
   */
  dataRepository::Group *
  GetConstitutiveModels()
  {
    return &m_constitutiveModels;
  }

  /**
   * @brief Get a pointer to the constitutive model.
   * @tparam T The type of the constitutive model.
   * @param name The name of the constitutive model.
   * @return A pointer to the constitutive model.
   */
  template< typename T = constitutive::ConstitutiveBase >
  T const *
  getConstitutiveModel( string const & name ) const
  {
    return m_constitutiveModels.GetGroup< T >( name );
  }

  /**
   * @copydoc getConstitutiveModel( string const & ) const
   */
  template< typename T = constitutive::ConstitutiveBase >
  T *
  getConstitutiveModel( string const & name )
  {
    return m_constitutiveModels.GetGroup< T >( name );
  }

  /**
   * @brief Get the type of element in this subregion.
   * @return a string specifying the type of element in this subregion
   *
   * See class FiniteElementBase for possible element type.
   */
  virtual string
  GetElementTypeString() const
  {
    return m_elementTypeString;
  }

  /**
   * @brief Set the type of element in this subregion.
   * @param[in] elementType a string specifying the element type
   */
  virtual void
  SetElementType( string const & elementType );

  /**
   * @brief Get the VTK ordering for this subregion.
   * @return the VTK node ordering
   */
  std::vector< int >
  getVTKNodeOrdering() const;

  ///@}

  /**
   * @brief A struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : ObjectManagerBase::viewKeyStruct
  {
    /// String key for the number of nodes per element in this subregion.
    static constexpr auto numNodesPerElementString = "numNodesPerElement";
    /// String key for the element-to-node relation
    static constexpr auto nodeListString = "nodeList";
    /// String key for the number of edges per element in this subregion.
    static constexpr auto numEdgesPerElementString = "numEdgesPerElement";
    /// String key for the element-to-edge relation
    static constexpr auto edgeListString = "edgeList";
    /// String key for the number of faces per element in this subregion.
    static constexpr auto numFacesPerElementString = "numFacesPerElement";
    /// String key for the element-to-face relation
    static constexpr auto faceListString = "faceList";
    /// String key for the member level field for the element center.
    static constexpr auto elementCenterString = "elementCenter";
    /// String key for the member level field for the element volume.
    static constexpr auto elementVolumeString = "elementVolume";
  };

  /**
   * @brief struct to serve as a container for group strings and keys
   * @struct groupKeyStruct
   */
  struct groupKeyStruct : public ObjectManagerBase::groupKeyStruct
  {
    /// String key for the group in which the constitutive models of this subregion are registered.
    static constexpr auto constitutiveModelsString = "ConstitutiveModels";
  };

private:
  /// Group in which the constitutive models of this subregion are registered
  dataRepository::Group m_constitutiveModels;

protected:
  /// Number of nodes per element in this subregion.
  localIndex m_numNodesPerElement;

  /// Number of independent nodes per element in this subregion.
  localIndex m_numIndependentNodesPerElement;

  /// Number of edges per element in this subregion.
  localIndex m_numEdgesPerElement;

  /// Number of faces per element in this subregion.
  localIndex m_numFacesPerElement;

  /// Member level field for the element center.
  array2d< real64 > m_elementCenter;

  /// Member level field for the element volume.
  array1d< real64 > m_elementVolume;

  /// Type of element in this subregion.
  string m_elementTypeString;

  /// Type of element in this subregion.
};

} /* namespace geosx */

#endif /* GEOSX_MESH_ELEMENTSUBREGIONBASE_HPP_ */
