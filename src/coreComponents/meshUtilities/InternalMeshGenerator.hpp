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
 * @file InternalMeshGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_INTERNALMESHGENERATOR_HPP
#define GEOSX_MESHUTILITIES_INTERNALMESHGENERATOR_HPP

#include <stdlib.h> /* srand, rand */
#include <time.h>   /* time */

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "MeshGeneratorBase.hpp"

namespace geosx
{
namespace dataRepository
{
/*
 * @name keys for InternalMesh object
 */
///@{
namespace keys
{
/// key for x coordinates
string const xCoords = "xCoords";
/// key for y coordinates
string const yCoords = "yCoords";
/// key for z coordinates
string const zCoords = "zCoords";
/// key for number of element in x-direction
string const xElems = "nx";
/// key for number of element in y-direction
string const yElems = "ny";
/// key for number of elemnt in z-direction
string const zElems = "nz";
/// key for x bias
string const xBias = "xBias";
/// key for y bias
string const yBias = "yBias";
/// key for z bias
string const zBias = "zBias";
/// key for cellBlock names
string const cellBlockNames = "cellBlockNames";
/// key for element type
string const elementTypes = "elementTypes";
/// key for triangle pattern identifier
string const trianglePattern = "trianglePattern";
}  // namespace keys
///@}

}  // namespace dataRepository

class NodeManager;
class DomainPartition;
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
  InternalMeshGenerator( const std::string & name, Group * const parent );

  virtual ~InternalMeshGenerator() override;

  /**
   * @brief Return the name of the InternalMeshGenerator in object Catalog.
   * @return string that contains the key name to InternalMeshGenerator in the Catalog
   */
  static string
  CatalogName()
  {
    return "InternalMesh";
  }

  //  void ProcessInputFile( xmlWrapper::xmlNode const & targetNode ) override;
  //
  //

  virtual void
  GenerateElementRegions( DomainPartition & domain ) override;

  /**
   * @brief Create a new geometric object (box, plane, etc) as a child of this group.
   * @param childKey the catalog key of the new geometric object to create
   * @param childName the name of the new geometric object in the repository
   * @return the group child
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  virtual void
  GenerateMesh( DomainPartition * const domain ) override;

  // virtual void GenerateNodesets( xmlWrapper::xmlNode const & targetNode,
  //                                NodeManager * nodeManager ) override;

  virtual void
  GetElemToNodesRelationInBox( const std::string & elementType,
                               const int index[],
                               const int & iEle,
                               int nodeIDInBox[],
                               const int size ) override;

  virtual void
  RemapMesh( dataRepository::Group * const domain ) override;

  //  int m_delayMeshDeformation;

protected:
  void
  PostProcessInput() override final;

private:
  /// Mesh number of dimension
  int m_dim;
  /// Array of vertex coordinates
  array1d< real64 > m_vertices[3];
  /// Ndim x nElem spatialized for element indexes
  integer_array m_nElems[3];
  /// Ndim x nElem spatialized array of element scaling factors
  array1d< real64 > m_nElemScaling[3];

  //bool m_useBias = false;
  /// Ndim x nElem spatialized array of element bias
  array1d< real64 > m_nElemBias[3];

  /// String array of region names
  string_array m_regionNames;

  /// Minimum extent of mesh dimensions
  realT m_min[3];
  /// Maximum extent of mesh dimensions
  realT m_max[3];

  //int m_numElems[3];
  /// Ndim x nBlock spatialized array of first elemnt index in the cellBlock
  integer_array m_firstElemIndexForBlock[3];
  /// Ndim x nBlock spatialized array of last elemnt index in the cellBlock
  integer_array m_lastElemIndexForBlock[3];

  //  realT m_wExtensionMin[3];
  //  realT m_wExtensionMax[3];
  //  int m_nExtensionLayersMin[3];
  //  int m_nExtensionLayersMax[3];
  //  realT m_extendedMin[3];
  //  realT m_extendedMax[3]; // This is the domain size after we apply n layers
  // of elements which are of the same size as the core elements.  We will move
  // these nodes to where they should be later when we finish the meshing.

  /// Array of number of elements per direction
  int m_numElemsTotal[3];

  //  realT m_commonRatioMin[3];
  //  realT m_commonRatioMax[3];

  // String array listing the element type present
  string_array m_elementType;

  /// Array of number of element per box
  array1d< integer > m_numElePerBox;

  /**
   * @brief Member variable for triangle pattern seletion.
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  int m_trianglePattern;

  /// Node perturbation amplitude value
  realT m_fPerturb = 0.0;
  /// Random seed for generation of the node perturbation field
  int m_randSeed;

  /**
   * @brief Knob to map onto a radial mesh.
   * @note if 0 mesh is not radial, if positive mesh is, it larger than 1
   */
  int m_mapToRadial = 0;
  ///@cond DO_NOT_DOCUMENT
  /// axis index for cartesian to radial coordinates mapping
  // internal temp var
  int m_meshAxis;
  realT m_meshTheta;
  realT m_meshPhi;
  realT m_meshRout;
  realT m_meshRact;
  /// @endcond

  /// skew angle in radians for skewed mesh generation
  realT m_skewAngle = 0;
  /// skew center for skew mesh generation
  R1Tensor m_skewCenter = { 0, 0, 0 };

  ///@cond DO_NOT_DOCUMENT
  //unused
  std::string m_meshDx, m_meshDy, m_meshDz;
  ///@endcond

  /**
 * @brief Convert ndim node spatialized index to node global index.
 * @param[in] node ndim spatialized array index
 */
  inline globalIndex
  NodeGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0] * ( m_numElemsTotal[1] + 1 ) * ( m_numElemsTotal[2] + 1 ) +
      index[1] * ( m_numElemsTotal[2] + 1 ) + index[2];
    return rval;
  }

  /**
 * @brief Convert ndim element spatialized index to element global index.
 * @param[in] element ndim spatialized array index
 */
  inline globalIndex
  ElemGlobalIndex( const int index[3] )
  {
    globalIndex rval = 0;

    rval = index[0] * m_numElemsTotal[1] * m_numElemsTotal[2] +
      index[1] * m_numElemsTotal[2] + index[2];
    return rval;
  }

  /**
   * @brief Construct the node position for a spatially indexed node.
   * @param[in] a ndim spatial index for the considered node
   * @param[in] trianglePattern triangle pattern identifier
   * @return coordinate of the input node
   *
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  inline R1Tensor
  NodePosition( const int a[3], int trianglePattern )
  {
    R1Tensor X;
    realT xInterval( 0 );

    int xPosIndex = 0;
    if( trianglePattern == 1 )
    {
      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block = 0; block < m_nElems[0].size(); ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[0][block];
      }
      xPosIndex = endingIndex;
    }

    for( int i = 0; i < 3; ++i )
    {
      int startingIndex = 0;
      int endingIndex = 0;
      int block = 0;
      for( block = 0; block < m_nElems[i].size(); ++block )
      {
        startingIndex = endingIndex;
        endingIndex = startingIndex + m_nElems[i][block];
        if( a[i] >= startingIndex && a[i] <= endingIndex )
        {
          break;
        }
      }
      realT min = m_vertices[i][block];
      realT max = m_vertices[i][block + 1];

      X[i] =
        min + ( max - min ) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

      // First check if m_nElemBias contains values
      // Otherwise the next test will cause a segfault when looking for "block"
      if( m_nElemBias[i].size() > 0 )
      {
        // Verify that the bias is non-zero and applied to more than one block:
        if( ( !isZero( m_nElemBias[i][block] ) ) && ( m_nElems[i][block] > 1 ) )
        {
          GEOSX_ERROR_IF( fabs( m_nElemBias[i][block] ) >= 1,
                          "Mesh bias must between -1 and 1!" );

          realT len = max - min;
          realT xmean = len / m_nElems[i][block];
          realT x0 = xmean * double( a[i] - startingIndex );
          realT chi = m_nElemBias[i][block] / ( xmean / len - 1.0 );
          realT dx = -x0 * chi + x0 * x0 * chi / len;
          X[i] += dx;
        }
      }

      // This is for creating regular triangle pattern
      if( i == 0 ) xInterval = ( max - min ) / m_nElems[i][block];
      if( trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 &&
          a[0] != xPosIndex )
        X[0] -= xInterval * 0.5;
    }

    return X;
  }

  /**
   * @brief
   * @param[in]
   * @return an array of the element center coordinates
   */
  inline R1Tensor
  ElemCenterPosition( const int k[3] )
  {
    R1Tensor X;

    for( int i = 0; i < 3; ++i )
    {
      X[i] =
        m_min[i] + ( m_max[i] - m_min[i] ) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
    }

    return X;
  }

public:
  /**
   * @brief Check if the mesh is a radial mesh.
   * @return true if the Internal mesh is radial, false else
   */
  inline bool
  isRadial()
  {
    bool rval = ( m_mapToRadial > 0 );
    return rval;
  }
};
}  // namespace geosx

#endif /* GEOSX_MESHUTILITIES_INTERNALMESHGENERATOR_HPP */
