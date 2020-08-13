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
 * @file VTMMeshGenerator.hpp
 *
 */

#ifndef GEOSX_MESHUTILITIES_VTMMESHGENERATOR_HPP
#define GEOSX_MESHUTILITIES_VTMMESHGENERATOR_HPP

#include "dataRepository/Group.hpp"
#include "codingUtilities/Utilities.hpp"
#include "MeshGeneratorBase.hpp"

#include "fileIO/vtm/VtmFile.hpp"

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
///key for file path access
string const filePath = "file";
}  // namespace keys
///@}
}  // namespace dataRepository

class NodeManager;
class DomainPartition;

/**
 *  @class VTMMeshGenerator
 *  @brief The VTMMeshGenerator class provides a class implementation of VTK genrated meshes (.vtm).
 */
class VTMMeshGenerator : public MeshGeneratorBase
{
public:
  /**
   * @brief Main constructor for MeshGenerator base class.
   * @param[in] name of the VTMMeshGenerator object
   * @param[in] parent the parent Group pointer for the MeshGenerator object
   */
  VTMMeshGenerator( const std::string & name, Group * const parent );

  virtual ~VTMMeshGenerator() override;

  /**
   * @brief Return the name of the VTMMeshGenerator in object Catalog
   * @return string that contains the key name to VTMMeshGenerator in the Catalog
   */
  static string
  CatalogName()
  {
    return "MeshFile";
  }

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
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  void
  PostProcessInput() override final;

private:
  /// Contains the path to the VTM file
  string m_fileName;

  /// VTM file data structure from VTK API
  VtmFile m_vtmFile;

  /**
   * @brief Convert ndim node spatialized index to node global index.
   * @param[in] node ndim spatialized array index
   */
  inline globalIndex
  NodeGlobalIndex( const int GEOSX_UNUSED_PARAM( index )[3] )
  {
    globalIndex rval = 0;
    /*

       rval = index[0]*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) + index[1]*(m_numElemsTotal[2]+1) + index[2];
     */
    return rval;
  }

  /**
   * @brief Convert ndim element spatialized index to element global index.
   * @param[in] element ndim spatialized array index
   */
  inline globalIndex
  ElemGlobalIndex( const int GEOSX_UNUSED_PARAM( index )[3] )
  {
    globalIndex rval = 0;
    /*

       rval = index[0]*m_numElemsTotal[1]*m_numElemsTotal[2] + index[1]*m_numElemsTotal[2] + index[2];
     */
    return rval;
  }

  /**
   * @brief Construct the node position for a spatially indexed node.
   * @param[in] a ndim spatial index for the considered node
   * @param[in] trianglePattern triangle pattern identifier
   *
   * @note In pattern 0, half nodes have 4 edges and the other half have 8; for Pattern 1, every node has 6.
   */
  inline R1Tensor
  NodePosition( const int GEOSX_UNUSED_PARAM( a )[3],
                int GEOSX_UNUSED_PARAM( trianglePattern ) )
  {
    R1Tensor X;
    /*
       realT xInterval(0);

       int xPosIndex = 0;
       if (trianglePattern == 1)
       {
       int startingIndex = 0;
       int endingIndex = 0;
       int block = 0;
       for( block=0 ; block<m_nElems[0].size() ; ++block )
       {
       startingIndex = endingIndex;
       endingIndex = startingIndex + m_nElems[0][block];
       }
       xPosIndex = endingIndex;
       }

       for( int i=0 ; i<3 ; ++i )
       {

       int startingIndex = 0;
       int endingIndex = 0;
       int block = 0;
       for( block=0 ; block<m_nElems[i].size() ; ++block )
       {
       startingIndex = endingIndex;
       endingIndex = startingIndex + m_nElems[i][block];
       if( a[i]>=startingIndex && a[i]<=endingIndex )
       {
        break;
       }
       }
       realT min = m_vertices[i][block];
       realT max = m_vertices[i][block+1];


       X[i] = min + (max-min) * ( double( a[i] - startingIndex ) / m_nElems[i][block] );

       if (( !isZero(m_nElemBias[i][block]) ) & (m_nElems[i][block]>1))
       {
       if (fabs(m_nElemBias[i][block]) >= 1)
       {
        GEOSX_ERROR("Mesh bias must between -1 and 1!");
       }

       realT len = max -  min;
       realT xmean = len / m_nElems[i][block];
       realT x0 = xmean * double( a[i] - startingIndex );
       realT chi = m_nElemBias[i][block]/(xmean/len - 1.0);
       realT dx = -x0*chi + x0*x0*chi/len;
       X[i] += dx;
       }

       // This is for creating regular triangle pattern
       if (i==0) xInterval = (max-min) / m_nElems[i][block];
       if (trianglePattern == 1 && i == 1 && a[1] % 2 == 1 && a[0] != 0 && a[0] != xPosIndex)
       X[0] -= xInterval * 0.5;
       }

     */
    return X;
  }

  /**
   * @brief
   * @param[in]
   * @return an array of the element center coordinates
   */
  inline R1Tensor
  ElemCenterPosition( const int GEOSX_UNUSED_PARAM( k )[3] )
  {
    R1Tensor X;

    /*
       for( int i=0 ; i<3 ; ++i )
       {
       X[i] = m_min[i] + (m_max[i]-m_min[i]) * ( ( k[i] + 0.5 ) / m_numElemsTotal[i] );
       }
     */

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
    /*
       bool rval = (m_mapToRadial > 0);
       return rval;
     */
    return false;
  }
};
}  // namespace geosx

#endif /* GEOSX_MESHUTILITIES_VTMMESHGENERATOR_HPP */
