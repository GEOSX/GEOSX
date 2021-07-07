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
 * @file InternalMeshGenerator.hpp
 */

#ifndef GEOSX_MESH_GENERATORS_INTERNALMESHGENERATOR_HPP
#define GEOSX_MESH_GENERATORS_INTERNALMESHGENERATOR_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "MeshGeneratorBase.hpp"

namespace geosx
{

class SpatialPartition;

/**
 * @class InternalMeshGenerator
 * @brief The InternalMeshGenerator class is a class handling GEOSX generated meshes.
 */
class InternalMeshGenerator : public MeshGeneratorBase
{
public:

  /**
   * @brief Main constructor for InternalMeshGenerator.
   * @param[in] name of the InternalMeshGenerator
   * @param[in] parent point to the parent Group of the InternalMeshGenerator
   */
  InternalMeshGenerator( const string & name, Group * const parent );

  virtual ~InternalMeshGenerator() override = default;

  /**
   * @brief Return the name of the InternalMeshGenerator in object Catalog.
   * @return string that contains the key name to InternalMeshGenerator in the Catalog
   */
  static string catalogName() { return "InternalMesh"; }

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group * createChild( string const & childKey, string const & childName ) override;

  virtual void generateMesh( DomainPartition & domain ) override;

  /**
   * @return Whether or not a Cartesian mesh is being generated.
   */
  virtual inline bool isCartesian() const
  {
    return true;
  }

  /**
   * @brief Reduce the number of nodes in a block coordinate direction for
   * @param partition The partitioning object
   * @param numNodes The number of nodes in each coordinate direction.
   */
  virtual void reduceNumNodesForPeriodicBoundary( SpatialPartition & partition,
                                                  integer (& numNodes) [3] )
  {
    GEOSX_UNUSED_VAR( partition, numNodes );
  };

  /**
   * @brief Alter the directional indices for when the ending index should be
   *   set to the beginning of the index as is the case with periodic
   *   boundaries.
   * @param partition The partitioning object
   * @param index The indices to be evaluated for periodic indexing.
   *   merging.
   */
  virtual void
  setNodeGlobalIndicesOnPeriodicBoundary( SpatialPartition & partition,
                                          int (& index)[3] )
  {
    GEOSX_UNUSED_VAR( partition, index );
  }

  /**
   * @brief Alter connectivity to adhere to a specific periodic type boundary.
   * @param globalIJK The global ijk indices.
   * @param numNodesInDir The number of elements in each direction globally.
   * @param firstElemIndexInPartition The first index in the partition in each direction.
   * @param nodeOfBox The connectivity.
   */
  virtual void setConnectivityForPeriodicBoundaries( int ( & globalIJK )[3],
                                                     integer const ( &numNodesInDir )[3],
                                                     int const ( &firstElemIndexInPartition )[3],
                                                     localIndex ( & nodeOfBox )[8] )
  {
    GEOSX_UNUSED_VAR( globalIJK, numNodesInDir, firstElemIndexInPartition, nodeOfBox );
  }

  /**
   * @brief Alter connectivity to adhere to a periodic type boundary.
   * @param component The direction of the periodicity.
   * @param globalIJK The global ijk indices.
   * @param numNodesInDir The number of elements in each direction globally.
   * @param firstElemIndexInPartition The first index in the partition in each direction.
   * @param nodeOfBox The connectivity.
   */
  void setConnectivityForPeriodicBoundary( int const component,
                                           int const (&globalIJK)[3],
                                           integer const (&numNodesInDir)[3],
                                           int const (&firstElemIndexInPartition)[3],
                                           localIndex ( &nodeOfBox )[8] );

  /**
   * @brief Performs a coordinate transformation of all nodes.
   * @param X The nodes coordinates.
   * @param nodeSets The name to node sets mapping.
   */
  virtual void coordinateTransformation( array2d< real64, nodes::REFERENCE_POSITION_PERM > & X, std::map< string, SortedArray< localIndex > > & nodeSets )
  {
    GEOSX_UNUSED_VAR( X );
    GEOSX_UNUSED_VAR( nodeSets );
  }


protected:

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    constexpr static char const * xCoordsString() { return "xCoords"; }
    constexpr static char const * yCoordsString() { return "yCoords"; }
    constexpr static char const * zCoordsString() { return "zCoords"; }
    constexpr static char const * xElemsString() { return "nx"; }
    constexpr static char const * yElemsString() { return "ny"; }
    constexpr static char const * zElemsString() { return "nz"; }
    constexpr static char const * xBiasString() { return "xBias"; }
    constexpr static char const * yBiasString() { return "yBias"; }
    constexpr static char const * zBiasString() { return "zBias"; }
    constexpr static char const * cellBlockNamesString() { return "cellBlockNames"; }
    constexpr static char const * elementTypesString() { return "elementTypes"; }
    constexpr static char const * trianglePatternString() { return "trianglePattern"; }
    constexpr static char const * meshTypeString() { return "meshType"; }
    constexpr static char const * positionToleranceString() { return "positionTolerance"; }
  };
  /// @endcond

  void postProcessInput() override;

  /**
   * @brief Get the label mapping of element vertices indexes onto node indexes for a type of element.
   * @param[in] elementType the string identifier of the element type
   * @param[in] index ndim-sptialized Element index.
   * @param[in] iEle the index of Element begin processed
   * @param[out] nodeIDInBox array to map element vertices index to node indexes
   * @param[in] size the number of node on the element
   *
   */
  void getElemToNodesRelationInBox( const string & elementType,
                                    int const (&index)[3],
                                    int const & iEle,
                                    int (& nodeIDInBox)[8],
                                    int const size );

  /// Mesh number of dimension
  int m_dim;

  /// Minimum extent of mesh dimensions
  real64 m_min[3];

  /// Maximum extent of mesh dimensions
  real64 m_max[3];

  /// Position tolerance for adding nodes to nodesets
  real64 m_coordinatePrecision;

  /// Array of vertex coordinates
  array1d< real64 > m_vertices[3];

  /// Ndim x nElem spatialized for element indexes
  array1d< integer > m_nElems[3];

  /// Ndim x nElem spatialized array of element bias
  array1d< real64 > m_nElemBias[3];

  /// Coordinates of the nodes in each direction which may be calculated based
  /// on some geometry specifications and "set" for later use when generating
  /// the mesh.
  array1d< real64 > m_setCoords[3];

private:

  /// String array of region names
  array1d< string > m_regionNames;

  /// Ndim x nBlock spatialized array of first element index in the cellBlock
  array1d< integer > m_firstElemIndexForBlock[3];

  /// Ndim x nBlock spatialized array of last element index in the cellBlock
  array1d< integer > m_lastElemIndexForBlock[3];

  /// Array of number of elements per direction
  int m_numElemsTotal[3];

  /// String array listing the element type present
  array1d< string > m_elementType;

  /// Array of number of element per box
  array1d< integer > m_numElePerBox;

  /**
   * @brief Member variable for triangle pattern seletion.
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  int m_trianglePattern;

  /// Node perturbation amplitude value
  real64 m_fPerturb = 0.0;

  /// Random seed for generation of the node perturbation field
  int m_randSeed = 0;

  /// Skew angle in radians for skewed mesh generation
  real64 m_skewAngle = 0;

  /// Skew center for skew mesh generation
  real64 m_skewCenter[3] = { 0, 0, 0 };



  /**
   * @brief Convert ndim node spatialized index to node global index.
   * @param[in] node ndim spatialized array index
   */
  inline globalIndex nodeGlobalIndex( int const index[3] )
  {
    return index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
  }

  /**
   * @brief Convert ndim element spatialized index to element global index.
   * @param[in] element ndim spatialized array index
   */
  inline globalIndex elemGlobalIndex( int const index[3] )
  {
    return index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
  }

  /**
   * @brief Construct the node position for a spatially indexed node.
   * @tparam OUT_VECTOR type of output vector X
   * @param[in] a ndim spatial index for the considered node
   * @param[in] trianglePattern triangle pattern identifier
   * @param[out] X the node coordinates
   *
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  template< typename OUT_VECTOR >
  inline void getNodePosition( int const (&a)[3], int trianglePattern, OUT_VECTOR && X )
  {
    real64 xInterval( 0 );

    int xPosIndex = 0;
    if( trianglePattern == 1 )
    {
      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block=0; block<m_nElems[0].size(); ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[0][block];
      }
      xPosIndex = endingIndex;
    }

    for( int i=0; i<3; ++i )
    {
      if( m_setCoords[i].size()>0 )
      {
        X[i] = m_setCoords[i][a[i]];
      }
      else
      {

        int startingIndex = 0;
        int endingIndex = 0;
        int block = 0;
        for( block=0; block<m_nElems[i].size(); ++block )
        {
          startingIndex = endingIndex;
          endingIndex = startingIndex + m_nElems[i][block];
          if( a[i]>=startingIndex && a[i]<=endingIndex )
          {
            break;
          }
        }
        real64 min = m_vertices[i][block];
        real64 max = m_vertices[i][block+1];


        X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

        // First check if m_nElemBias contains values
        // Otherwise the next test will cause a segfault when looking for "block"
        if( m_nElemBias[i].size()>0 )
        {
          // Verify that the bias is non-zero and applied to more than one block:
          if( ( !isZero( m_nElemBias[i][block] ) ) && (m_nElems[i][block]>1))
          {
            GEOSX_ERROR_IF( fabs( m_nElemBias[i][block] ) >= 1, "Mesh bias must between -1 and 1!" );

            real64 len = max -  min;
            real64 xmean = len / m_nElems[i][block];
            real64 x0 = xmean * double( a[i] - startingIndex );
            real64 chi = m_nElemBias[i][block]/(xmean/len - 1.0);
            real64 dx = -x0*chi + x0*x0*chi/len;
            X[i] += dx;
          }
        }

        // This is for creating regular triangle pattern
        if( i==0 ) xInterval = (max-min) / m_nElems[i][block];
        if( trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex )
          X[0] -= xInterval * 0.5;
      }
    }
  }

  /**
   * @brief
   * @tparam OUT_VECTOR type of output vector X
   * @param[in] k the ijk-index of the element
   * @param[out] X the element center coordinates
   */
  template< typename OUT_VECTOR >
  inline void getElemCenterPosition( int const k[3], OUT_VECTOR && X )
  {
    for( int i=0; i<3; ++i )
    {
      X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }
  }

public:


};

//ENUM_STRINGS( InternalMeshGenerator::MeshType, "Cartesian", "Cylindrical", "CylindricalSquareBoundary" )

} /* namespace geosx */

#endif /* GEOSX_MESH_GENERATORS_INTERNALMESHGENERATOR_HPP */
