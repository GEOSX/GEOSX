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


#ifndef GEOSX_MESH_CELLELEMENTSUBREGION_HPP_
#define GEOSX_MESH_CELLELEMENTSUBREGION_HPP_

#include "mesh/generators/CellBlockABC.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/FaceManager.hpp"
#include "mesh/utilities/ComputationalGeometry.hpp"
#include "ElementSubRegionBase.hpp"


namespace geosx
{

class MeshLevel;

/**
 * @class CellElementSubRegion
 * Class deriving from CellBlock further specializing the element subregion
 * for a cell element. This is the class used in the physics solvers to
 * represent a collection of mesh cell elements
 */
class CellElementSubRegion : public ElementSubRegionBase
{
public:

  /// Alias for the type of the element-to-node map
  using NodeMapType = InterObjectRelation< array2d< localIndex, cells::NODE_MAP_PERMUTATION > >;
  /// Alias for the type of the element-to-edge map
  using EdgeMapType = FixedOneToManyRelation;
  /// Alias for the type of the element-to-face map
  using FaceMapType = FixedOneToManyRelation;
  /// Type of map between cell blocks and embedded elements
  using EmbSurfMapType = InterObjectRelation< ArrayOfArrays< localIndex > >;

  /**
   * @brief Const getter for the catalog name.
   * @return the name of this type in the catalog
   */
  static const string catalogName()
  { return "CellElementSubRegion"; } // FIXME I don't know what I am doing

  /**
   * @copydoc catalogName()
   */
  virtual const string getCatalogName() const override final
  { return CellElementSubRegion::catalogName(); }

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for this class.
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  CellElementSubRegion( string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~CellElementSubRegion() override;

  ///@}

  /**
   * @name Helpers for CellElementSubRegion construction
   */
  ///@{

  /**
   * @brief Fill the CellElementSubRegion by copying those of the source CellBlock
   * @param source the CellBlock whose properties (connectivity info) will be copied
   */
  void copyFromCellBlock( CellBlockABC & source );

  ///@}

  /**
   * @brief Add fractured element to list and relative entries to the map.
   * @param cellElemIndex cell element index
   * @param embSurfIndex embedded surface element index
   */
  void addFracturedElement( localIndex const cellElemIndex,
                            localIndex const embSurfIndex );

  /**
   * @name Overriding packing / Unpacking functions
   */
  ///@{

  virtual void viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const override;

  virtual localIndex packUpDownMapsSize( arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex packUpDownMaps( buffer_unit_type * & buffer,
                                     arrayView1d< localIndex const > const & packList ) const override;

  virtual localIndex unpackUpDownMaps( buffer_unit_type const * & buffer,
                                       array1d< localIndex > & packList,
                                       bool const overwriteUpMaps,
                                       bool const overwriteDownMaps ) override;

  /**
   * @brief Computes the pack size of the set of fractured elements and the element-to-embeddedSurfaces map of the elements in the @
   * packList.
   * @param packList The element we want packed.
   * @param  embeddedSurfacesLocalToGlobal LocaltoGlobal map of the embedded surfaces.
   * @return The packed size.
   */
  localIndex packFracturedElementsSize( arrayView1d< localIndex const > const & packList,
                                        arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;
  /**
   * @brief Packs the set of fractured elements and the element-to-embeddedSurfaces map of the elements in the @ packList.
   * @param buffer The buffer that will receive the packed data.
   * @param packList The element we want packed.
   * @param embeddedSurfacesLocalToGlobal LocaltoGlobal map of the embedded surfaces.
   * @return The packed size.
   */
  localIndex packFracturedElements( buffer_unit_type * & buffer,
                                    arrayView1d< localIndex const > const & packList,
                                    arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;

  /**
   * @brief Unpacks the set of fractured elemetn and the element-to-embeddedSurfaces map from @p buffer.
   * @param buffer The buffer containing the packed data.
   * @param packList The (un)packed element.
   * @param embeddedSurfacesGlobalToLocal GlobalToLocal map of the embedded surfaces.
   * @return The unpacked size.
   */
  localIndex unpackFracturedElements( buffer_unit_type const * & buffer,
                                      localIndex_array & packList,
                                      unordered_map< globalIndex, localIndex > const & embeddedSurfacesGlobalToLocal );

  virtual void fixUpDownMaps( bool const clearIfUnmapped ) final override;

  ///@}

  /**
   * @name Miscellaneous
   */
  ///@{

  /**
   * @brief Helper function to apply a lambda function over all constructive groups
   * @tparam LAMBDA the type of the lambda function
   * @param lambda the lambda function
   */
  template< typename LAMBDA >
  void forMaterials( LAMBDA lambda )
  {

    for( auto & constitutiveGroup : m_constitutiveGrouping )
    {
      lambda( constitutiveGroup );
    }
  }

  ///@}

  /**
   * @brief struct to serve as a container for variable strings and keys
   * @struct viewKeyStruct
   */
  struct viewKeyStruct : public ElementSubRegionBase::viewKeyStruct
  {
    /// @return String key for the constitutive point volume fraction
    static constexpr char const * constitutivePointVolumeFractionString() { return "ConstitutivePointVolumeFraction"; }
    /// @return String key for the derivatives of the shape functions with respect to the reference configuration
    static constexpr char const * dNdXString() { return "dNdX"; }
    /// @return String key for the derivative of the jacobian.
    static constexpr char const * detJString() { return "detJ"; }
    /// @return String key for the constitutive grouping
    static constexpr char const * constitutiveGroupingString() { return "ConstitutiveGrouping"; }
    /// @return String key for the constitutive map
    static constexpr char const * constitutiveMapString() { return "ConstitutiveMap"; }
    /// @return String key to embSurfMap
    static constexpr char const * toEmbSurfString() { return "ToEmbeddedSurfaces"; }
    /// @return String key to fracturedCells
    static constexpr char const * fracturedCellsString() { return "fracturedCells"; }

    /// ViewKey for the constitutive grouping
    dataRepository::ViewKey constitutiveGrouping  = { constitutiveGroupingString() };
    /// ViewKey for the constitutive map
    dataRepository::ViewKey constitutiveMap       = { constitutiveMapString() };
  }
  /// viewKey struct for the CellElementSubRegion class
  m_CellBlockSubRegionViewKeys;

  virtual viewKeyStruct & viewKeys() override { return m_CellBlockSubRegionViewKeys; }
  virtual viewKeyStruct const & viewKeys() const override { return m_CellBlockSubRegionViewKeys; }

  /**
   * @brief Get the local indices of the nodes in a face of the element.
   * @param elementIndex the local index of the target element
   * @param localFaceIndex the local index of the target face in the element  (this will be 0-numFacesInElement)
   * @param nodeIndices a reference to the array of node indices of the face
   */
  void getFaceNodes( localIndex const elementIndex,
                     localIndex const localFaceIndex,
                     array1d< localIndex > & nodeIndices ) const;

  /**
   * @brief Get the element-to-node map.
   * @return a reference to the element-to-node map
   */
  NodeMapType & nodeList() { return m_toNodesRelation; }

  /**
   * @copydoc nodeList()
   */
  NodeMapType const & nodeList() const { return m_toNodesRelation; }

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
   * @brief Get the element-to-edge map.
   * @return a reference to element-to-edge map
   */
  FixedOneToManyRelation & edgeList() { return m_toEdgesRelation; }

  /**
   * @copydoc edgeList()
   */
  FixedOneToManyRelation const & edgeList() const { return m_toEdgesRelation; }

  /**
   * @brief Get the element-to-face map.
   * @return a reference to the element to face map
   */
  FixedOneToManyRelation & faceList() { return m_toFacesRelation; }

  /**
   * @copydoc faceList()
   */
  FixedOneToManyRelation const & faceList() const { return m_toFacesRelation; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  array4d< real64 > & dNdX()
  { return m_dNdX; }

  /**
   * @brief @return The array of shape function derivatives.
   */
  arrayView4d< real64 const > dNdX() const
  { return m_dNdX; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  array2d< real64 > & detJ()
  { return m_detJ; }

  /**
   * @brief @return The array of jacobian determinantes.
   */
  arrayView2d< real64 const > detJ() const
  { return m_detJ; }

  /**
   * @brief @return The sorted array of local fractured elements.
   */
  SortedArray< localIndex > & fracturedElementsList()
  { return m_fracturedCells; }

  /**
   * @brief @return The sorted array view of local fractured elements.
   */
  SortedArrayView< localIndex const > const fracturedElementsList() const
  { return m_fracturedCells.toViewConst(); }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType & embeddedSurfacesList() { return m_toEmbeddedSurfaces; }

  /**
   * @brief @return The map to the embedded surfaces
   */
  EmbSurfMapType const & embeddedSurfacesList() const { return m_toEmbeddedSurfaces; }

  /**
   * @brief Compute the center of each element in the subregion.
   * @param[in] X an arrayView of (const) node positions
   */
  void calculateElementCenters( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    arrayView2d< real64 > const & elementCenters = m_elementCenter;
    localIndex nNodes = numNodesPerElement();

    forAll< parallelHostPolicy >( size(), [=]( localIndex const k )
    {
      LvArray::tensorOps::copy< 3 >( elementCenters[ k ], X[ m_toNodesRelation( k, 0 ) ] );
      for( localIndex a = 1; a < nNodes; ++a )
      {
        LvArray::tensorOps::add< 3 >( elementCenters[ k ], X[ m_toNodesRelation( k, a ) ] );
      }

      LvArray::tensorOps::scale< 3 >( elementCenters[ k ], 1.0 / nNodes );
    } );
  }

  void calculateElementGeometricQuantities( NodeManager const & nodeManager,
                                            FaceManager const & faceManager ) override;

private:

  /// Map used for constitutive grouping
  map< string, localIndex_array > m_constitutiveGrouping;

  /// Array of constitutive point volume fraction
  array3d< real64 > m_constitutivePointVolumeFraction;

  /// Element-to-node relation
  NodeMapType m_toNodesRelation;

  /// Element-to-edge relation
  EdgeMapType m_toEdgesRelation;

  /// Element-to-face relation
  FaceMapType m_toFacesRelation;

  /// Name of the properties registered from an external mesh
  string_array m_externalPropertyNames;

  /**
   * @brief Compute the volume of the k-th element in the subregion.
   * @param[in] k the index of the element in the subregion
   * @param[in] X an arrayView of (const) node positions
   */
  inline void calculateCellVolumesKernel( localIndex const k,
                                          arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X ) const
  {
    LvArray::tensorOps::fill< 3 >( m_elementCenter[ k ], 0 );

    real64 Xlocal[10][3];

    for( localIndex a = 0; a < m_numNodesPerElement; ++a )
    {
      LvArray::tensorOps::copy< 3 >( Xlocal[ a ], X[ m_toNodesRelation( k, a ) ] );
      LvArray::tensorOps::add< 3 >( m_elementCenter[ k ], X[ m_toNodesRelation( k, a ) ] );
    }
    LvArray::tensorOps::scale< 3 >( m_elementCenter[ k ], 1.0 / m_numNodesPerElement );

    if( m_numNodesPerElement == 8 )
    {
      m_elementVolume[k] = computationalGeometry::HexVolume( Xlocal );
    }
    else if( m_numNodesPerElement == 4 )
    {
      m_elementVolume[k] = computationalGeometry::TetVolume( Xlocal );
    }
    else if( m_numNodesPerElement == 6 )
    {
      m_elementVolume[k] = computationalGeometry::WedgeVolume( Xlocal );
    }
    else if( m_numNodesPerElement == 5 )
    {
      m_elementVolume[k] = computationalGeometry::PyramidVolume( Xlocal );
    }
    else
    {
      GEOSX_ERROR( "GEOX does not support cells with " << m_numNodesPerElement << " nodes" );
    }
  }

  /// The array of shape function derivaties.
  array4d< real64 > m_dNdX;

  /// The array of jacobian determinantes.
  array2d< real64 > m_detJ;

  /// Map of unmapped global indices in the element-to-node map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInNodelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInEdgelist;

  /// Map of unmapped global indices in the element-to-face map
  map< localIndex, array1d< globalIndex > > m_unmappedGlobalIndicesInFacelist;

  /// List of the fractured elements on this rank
  SortedArray< localIndex > m_fracturedCells;

  /// Map from local Cell Elements to Embedded Surfaces
  EmbSurfMapType m_toEmbeddedSurfaces;

  /**
   * @brief Defines the (unique) element type of this cell element region.
   * @param elementType The string defining the element type (like C3D8...)
   */
  void setElementType( string const & elementType ) final;

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

  /**
   * @brief Links the managers to their mappings.
   * @param mesh Holds the node/edge/face managers.
   *
   * Defines links the element to nodes, edges and faces to their respective node/edge/face managers.
   */
  void setupRelatedObjectsInRelations( MeshLevel const & mesh ) override;

  template< bool DOPACK >
  localIndex packFracturedElementsPrivate( buffer_unit_type * & buffer,
                                           arrayView1d< localIndex const > const & packList,
                                           arrayView1d< globalIndex const > const & embeddedSurfacesLocalToGlobal ) const;
};

} /* namespace geosx */

#endif /* GEOSX_MESH_CELLELEMENTSUBREGION_HPP_ */
