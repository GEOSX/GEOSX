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
 * @file InternalMeshGenerator.cpp
 */

#include "InternalMeshGenerator.hpp"


#include "codingUtilities/StringUtilities.hpp"
#include "common/DataTypes.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"


#include <math.h>
#include <algorithm>

namespace geosx
{
using namespace dataRepository;

InternalMeshGenerator::InternalMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent ),
  m_dim( 3 ),
  m_min(),
  m_max(),
  m_coordinatePrecision( 1e-10 ),
  m_vertices{},
  m_nElems{},
  m_nElemBias{},
  m_setCoords{}
{
  registerWrapper( viewKeyStruct::xCoordsString(), &(m_vertices[0]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "x-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::yCoordsString(), &(m_vertices[1]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "y-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::zCoordsString(), &(m_vertices[2]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "z-coordinates of each mesh block vertex" );

  registerWrapper( viewKeyStruct::xElemsString(), &(m_nElems[0]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the x-direction within each mesh block" );

  registerWrapper( viewKeyStruct::yElemsString(), &(m_nElems[1]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the y-direction within each mesh block" );

  registerWrapper( viewKeyStruct::zElemsString(), &(m_nElems[2]) ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Number of elements in the z-direction within each mesh block" );

  registerWrapper( viewKeyStruct::xBiasString(), &(m_nElemBias[0]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the x-direction within each mesh block (dx_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::yBiasString(), &(m_nElemBias[1]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the y-direction within each mesh block (dy_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::zBiasString(), &(m_nElemBias[2]) ).
    setApplyDefaultValue( 1.0 ).
    setSizedFromParent( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Bias of element sizes in the z-direction within each mesh block (dz_left=(1+b)*L/N, dz_right=(1-b)*L/N)" );

  registerWrapper( viewKeyStruct::cellBlockNamesString(), &m_regionNames ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Names of each mesh block" );

  registerWrapper( viewKeyStruct::elementTypesString(), &m_elementType ).
    setInputFlag( InputFlags::REQUIRED ).
    setSizedFromParent( 0 ).
    setDescription( "Element types of each mesh block" );

  registerWrapper( viewKeyStruct::trianglePatternString(), &m_trianglePattern ).
    setApplyDefaultValue( 0 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setDescription( "Pattern by which to decompose the hex mesh into prisms (more explanation required)" );

  registerWrapper( viewKeyStruct::positionToleranceString(), &m_coordinatePrecision ).
    setApplyDefaultValue( 1e-10 ).
    setInputFlag( InputFlags::OPTIONAL ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "A position tolerance to verify if a node belong to a nodeset" );
}

void InternalMeshGenerator::postProcessInput()
{
  if( m_elementType[0] == "C3D8" || m_elementType[0] == "C3D4" || m_elementType[0] == "C3D6" )
  {
    m_dim = 3;
  }
  else if( m_elementType[0] == "CPE4" || m_elementType[0] == "STRI" )
  {
    m_dim = 2;
  }
  else
  {
    GEOSX_ERROR( "InternalMeshGenerator: incorrect element type!" );
  }

  {
    // Check for vertex/element matching
    bool failFlag = false;
    for( int i=0; i<m_dim; ++i )
    {
      failFlag += ( m_nElems[i].size() != m_vertices[i].size()-1 );
    }
    if( failFlag )
    {
      GEOSX_ERROR( "vertex/element mismatch InternalMeshGenerator::ReadXMLPost()" );
    }

    // If specified, check to make sure bias values have the correct length
    for( int i=0; i<m_dim; ++i )
    {
      if( m_nElemBias[i].size() > 0 )
      {
        failFlag += ( m_nElems[i].size() != m_nElemBias[i].size() );
      }
    }
    if( failFlag )
    {
      GEOSX_ERROR( "element/bias mismatch InternalMeshGenerator::ReadXMLPost()" );
    }
  }

  m_numElePerBox.resize( m_nElems[0].size() * m_nElems[1].size() * m_nElems[2].size());

  if( LvArray::integerConversion< long >( m_elementType.size()) != m_numElePerBox.size())
  {
    if( m_elementType.size() == 1 )
    {
      m_elementType.resize( m_numElePerBox.size());
      for( localIndex i=1; i< m_elementType.size(); ++i )
      {
        m_elementType[i] = m_elementType[0];
      }
    }
    else
    {
      GEOSX_ERROR( "InternalMeshGenerator: The number of element types is inconsistent with the number of total block." );
    }
  }

  for( localIndex i = 0; i < LvArray::integerConversion< localIndex >( m_elementType.size() ); ++i )
  {
    if( m_elementType[i] == "C3D8" )
    {
      m_numElePerBox[i] = 1;
      m_dim = 3;
    }
    else if( m_elementType[i] == "C3D4" )
    {
      m_numElePerBox[i] = 6;
      m_dim = 3;
    }
    else if( m_elementType[i] == "C3D6" )
    {
      m_numElePerBox[i] = 2;
      m_dim = 3;
    }
    else if( m_elementType[i] == "CPE4" )
    {
      m_numElePerBox[i] = 1;
      m_dim = 2;
    }
    else if( m_elementType[i] == "STRI" )
    {
      m_numElePerBox[i] = 2;
      m_dim = 2;
    }
  }

  {
    localIndex numBlocks = 1;
    for( int i=0; i<m_dim; ++i )
    {
      numBlocks *= m_nElems[i].size();
    }
    if( numBlocks != m_regionNames.size() )
    {
      if( m_regionNames.size() == 1 )
      {
        m_regionNames.resize( numBlocks );
        for( localIndex i=1; i< m_elementType.size(); ++i )
        {
          m_regionNames[i] = m_regionNames[0];
        }
      }
      else
      {
        GEOSX_ERROR( "Incorrect number of regionLayout entries specified in InternalMeshGenerator::ReadXML()" );
      }
    }
  }

  for( int i=0; i<3; ++i )
  {
    m_min[i] = m_vertices[i].front();
    m_max[i] = m_vertices[i].back();
  }

  for( int dir=0; dir<3; ++dir )
  {
    m_firstElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_lastElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_firstElemIndexForBlock[dir][0] = 0;
    m_lastElemIndexForBlock[dir][0] = m_nElems[dir][0]-1;
    for( int block=1; block<m_nElems[dir].size(); ++block )
    {
      m_firstElemIndexForBlock[dir][block] = m_lastElemIndexForBlock[dir][block-1] + 1;
      m_lastElemIndexForBlock[dir][block] = m_firstElemIndexForBlock[dir][block] + m_nElems[dir][block]-1;
    }
  }

  m_fPerturb = 0.0;
}

Group * InternalMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                                            string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/**
 * @param partition
 * @param domain
 */
void InternalMeshGenerator::generateMesh( DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;

  Group & meshBodies = domain.getGroup( string( "MeshBodies" ));
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );
  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  // Make sure that the node manager fields are initialized

  CellBlockManager & elementManager = domain.getGroup< CellBlockManager >( keys::cellManager );
  Group & nodeSets = nodeManager.sets();

  SpatialPartition & partition = dynamic_cast< SpatialPartition & >(domain.getReference< PartitionBase >( keys::partitionManager ) );

//  bool isRadialWithOneThetaPartition = false;

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & cellBlockName : m_regionNames )
  {
    CellBlock & cellBlock = elementManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( cellBlockName );
    string elementType = m_elementType[aa++];
    cellBlock.setElementType( elementType );
  }

  SortedArray< localIndex > & xnegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "xneg" ) ).reference();
  SortedArray< localIndex > & xposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "xpos" ) ).reference();
  SortedArray< localIndex > & ynegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "yneg" ) ).reference();
  SortedArray< localIndex > & yposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "ypos" ) ).reference();
  SortedArray< localIndex > & znegNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "zneg" ) ).reference();
  SortedArray< localIndex > & zposNodes = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "zpos" ) ).reference();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();

  // Partition based on even spacing to get load balance
  // Partition geometrical boundaries will be corrected in the end.
  {
    m_min[0] = m_vertices[0].front();
    m_min[1] = m_vertices[1].front();
    m_min[2] = m_vertices[2].front();

    m_max[0] = m_vertices[0].back();
    m_max[1] = m_vertices[1].back();
    m_max[2] = m_vertices[2].back();

    partition.setSizes( m_min, m_max );

    real64 size[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( m_max );
    LvArray::tensorOps::subtract< 3 >( size, m_min );
    meshBody.setGlobalLengthScale( LvArray::tensorOps::l2Norm< 3 >( size ) );
  }

  // Find elemCenters for even uniform element sizes
  array1d< array1d< real64 > > elemCenterCoords( 3 );
  for( int i = 0; i < 3; ++i )
  {
    m_numElemsTotal[i] = 0;
    for( int block = 0; block < m_nElems[i].size(); ++block )
    {
      m_numElemsTotal[i] += m_nElems[i][block];
    }

    elemCenterCoords[i].resize( m_numElemsTotal[i] );
    array1d< real64 > elemCenterCoordsLocal( m_numElemsTotal[i] );
    for( int k = 0; k < m_numElemsTotal[i]; ++k )
    {
      elemCenterCoordsLocal[k] = m_min[i] + ( m_max[i] - m_min[i] ) * ( k + 0.5 ) / m_numElemsTotal[i];
    }
    MpiWrapper::allReduce( elemCenterCoordsLocal.data(),
                           elemCenterCoords[i].data(),
                           m_numElemsTotal[i],
                           MPI_MAX,
                           MPI_COMM_GEOSX );
  }

  // Find starting/ending index
  // Get the first and last indices in this partition each direction
  int firstElemIndexInPartition[3] = { -1, -1, -1 };
  int lastElemIndexInPartition[3] = { -2, -2, -2 };

  for( int i = 0; i < 3; ++i )
  {
    //    firstElemIndexInPartition[i] = -1;
    //    lastElemIndexInPartition[i] = -2;
    for( int k = 0; k < m_numElemsTotal[i]; ++k )
    {
      if( partition.isCoordInPartition( elemCenterCoords[i][k], i ) )
      {
        firstElemIndexInPartition[i] = k;
        break;
      }
    }

    if( firstElemIndexInPartition[i] > -1 )
    {
      for( int k = firstElemIndexInPartition[i]; k < m_numElemsTotal[i]; ++k )
      {
        if( partition.isCoordInPartition( elemCenterCoords[i][k], i ) )
        {
          lastElemIndexInPartition[i] = k;
        }
      }
    }
  }

  // Calculate number of elements in this partition from each region, and the
  // total number of nodes

  std::map< string, int > numElemsInRegions;
  std::map< string, string > elemTypeInRegions;

  array1d< integer > firstElemIndexForBlockInPartition[3];
  array1d< integer > lastElemIndexForBlockInPartition[3];

  for( int dir = 0; dir < 3; ++dir )
  {
    firstElemIndexForBlockInPartition[dir] = m_firstElemIndexForBlock[dir];
    lastElemIndexForBlockInPartition[dir] = m_lastElemIndexForBlock[dir];

    for( int block = 0; block < m_nElems[dir].size(); ++block )
    {
      if( firstElemIndexForBlockInPartition[dir][block] > lastElemIndexInPartition[dir] ||
          lastElemIndexForBlockInPartition[dir][block] < firstElemIndexInPartition[dir] )
      {
        firstElemIndexForBlockInPartition[dir][block] = -1;
        lastElemIndexForBlockInPartition[dir][block] = -2;
      }
      else
      {
        if( firstElemIndexForBlockInPartition[dir][block] < firstElemIndexInPartition[dir] )
        {
          firstElemIndexForBlockInPartition[dir][block] = firstElemIndexInPartition[dir];
        }
        if( lastElemIndexForBlockInPartition[dir][block] > lastElemIndexInPartition[dir] )
        {
          lastElemIndexForBlockInPartition[dir][block] = lastElemIndexInPartition[dir];
        }
      }
    }
  }

  // TODO This needs to be rewritten for dimensions lower than 3.
  localIndex regionOffset = 0;
  for( int iblock = 0; iblock < m_nElems[0].size(); ++iblock )
  {
    for( int jblock = 0; jblock < m_nElems[1].size(); ++jblock )
    {
      for( int kblock = 0; kblock < m_nElems[2].size(); ++kblock, ++regionOffset )
      {
        numElemsInRegions[ m_regionNames[ regionOffset ] ] = 0;
        elemTypeInRegions[ m_regionNames[ regionOffset ] ] = "";
      }
    }
  }

  regionOffset = 0;
  {
    localIndex iR = 0;
    for( int iblock = 0; iblock < m_nElems[0].size(); ++iblock )
    {
      for( int jblock = 0; jblock < m_nElems[1].size(); ++jblock )
      {
        for( int kblock = 0; kblock < m_nElems[2].size(); ++kblock, ++regionOffset, ++iR )
        {
          int numElemsInRegion = 1;
          numElemsInRegion *= lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1;

          if( m_dim > 1 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1;
          }
          if( m_dim > 2 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1;
          }

          numElemsInRegion *= m_numElePerBox[iR];
          numElemsInRegions[ m_regionNames[ regionOffset ] ] += numElemsInRegion;
          elemTypeInRegions[ m_regionNames[ regionOffset ] ] = m_elementType[iR];
        }
      }
    }
  }

  localIndex numNodes = 1;
  integer numNodesInDir[3] = { 1, 1, 1 };

  for( int i = 0; i < m_dim; ++i )
  {
    numNodesInDir[i] = lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
  }
  reduceNumNodesForPeriodicBoundary( partition, numNodesInDir );
  numNodes = numNodesInDir[0] * numNodesInDir[1] * numNodesInDir[2];

  nodeManager.resize( numNodes );
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();

  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  {
    localIndex localNodeIndex = 0;
    for( int i = 0; i < numNodesInDir[0]; ++i )
    {
      for( int j = 0; j < numNodesInDir[1]; ++j )
      {
        for( int k = 0; k < numNodesInDir[2]; ++k )
        {
          int globalIJK[3] = { i, j, k };

          for( int a = 0; a < m_dim; ++a )
          {
            globalIJK[a] += firstElemIndexInPartition[a];
          }

          getNodePosition( globalIJK, m_trianglePattern, X[localNodeIndex] );

          // Alter global node map for radial mesh
          setNodeGlobalIndicesOnPeriodicBoundary( partition, globalIJK );

          nodeLocalToGlobal[localNodeIndex] = nodeGlobalIndex( globalIJK );

          // Cartesian-specific nodesets
          if( isCartesian() )
          {
            if( isEqual( X( localNodeIndex, 0 ), m_min[0], m_coordinatePrecision ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 0 ), m_max[0], m_coordinatePrecision ) )
            {
              xposNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_min[1], m_coordinatePrecision ) )
            {
              ynegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_max[1], m_coordinatePrecision ) )
            {
              yposNodes.insert( localNodeIndex );
            }
          }

          // General nodesets
          if( isEqual( X( localNodeIndex, 2 ), m_min[2], m_coordinatePrecision ) )
          {
            znegNodes.insert( localNodeIndex );
          }
          if( isEqual( X( localNodeIndex, 2 ), m_max[2], m_coordinatePrecision ) )
          {
            zposNodes.insert( localNodeIndex );
          }
          allNodes.insert( localNodeIndex );

          ++localNodeIndex;
        }
      }
    }
  }

  {
    array1d< integer > numElements;
    array1d< string > elementRegionNames;
    array1d< string > elementTypes;
    std::map< string, localIndex > localElemIndexInRegion;

    for( std::map< string, int >::iterator iterNumElemsInRegion = numElemsInRegions.begin();
         iterNumElemsInRegion != numElemsInRegions.end(); ++iterNumElemsInRegion )
    {

      numElements.emplace_back( iterNumElemsInRegion->second );
      elementRegionNames.emplace_back( iterNumElemsInRegion->first );
      elementTypes.emplace_back( elemTypeInRegions[iterNumElemsInRegion->first] );

      localElemIndexInRegion[iterNumElemsInRegion->first] = 0;
    }

    elementManager.resize( numElements, elementRegionNames, elementTypes );

    // Assign global numbers to elements
    regionOffset = 0;
    SortedArray< string > processedRegionNames;
    localIndex iR = 0;

    // Reset the number of nodes in each dimension in case of periodic BCs so the element firstNodeIndex 
    //  calculation is correct? Not actually needed in parallel since we still have ghost nodes in that case and
    //  the count has not been altered due to periodicity.
    if( std::any_of( partition.m_Periodic.begin(), partition.m_Periodic.end(), [](int & dimPeriodic) { return dimPeriodic == 1; } ) )
    {
      for( int i = 0; i < m_dim; ++i )
      {
        numNodesInDir[i] = lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
      }
      numNodes = numNodesInDir[0] * numNodesInDir[1] * numNodesInDir[2];
    }

    for( int iblock = 0; iblock < m_nElems[0].size(); ++iblock )
    {
      for( int jblock = 0; jblock < m_nElems[1].size(); ++jblock )
      {
        for( int kblock = 0; kblock < m_nElems[2].size(); ++kblock, ++regionOffset, ++iR )
        {
          CellBlock & elemRegion =  elementManager.getRegion( m_regionNames[ regionOffset ] );
          int const numNodesPerElem = LvArray::integerConversion< int >( elemRegion.numNodesPerElement());
          integer nodeIDInBox[ 8 ];

          arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodes = elemRegion.nodeList();
          arrayView1d< globalIndex > const & elemLocalToGlobal = elemRegion.localToGlobalMap();

          int numElemsInDirForBlock[3] =
          { lastElemIndexForBlockInPartition[0][iblock] - firstElemIndexForBlockInPartition[0][iblock] + 1,
            lastElemIndexForBlockInPartition[1][jblock] - firstElemIndexForBlockInPartition[1][jblock] + 1,
            lastElemIndexForBlockInPartition[2][kblock] - firstElemIndexForBlockInPartition[2][kblock] + 1 };

          for( int i = 0; i < numElemsInDirForBlock[0]; ++i )
          {
            for( int j = 0; j < numElemsInDirForBlock[1]; ++j )
            {
              for( int k = 0; k < numElemsInDirForBlock[2]; ++k )
              {
                int globalIJK[3] =
                { i + firstElemIndexForBlockInPartition[0][iblock],
                  j + firstElemIndexForBlockInPartition[1][jblock],
                  k + firstElemIndexForBlockInPartition[2][kblock] };

                const localIndex firstNodeIndex = numNodesInDir[1] * numNodesInDir[2] * ( globalIJK[0] - firstElemIndexInPartition[0] )
                                                  + numNodesInDir[2] * ( globalIJK[1] - firstElemIndexInPartition[1] )
                                                  + ( globalIJK[2] - firstElemIndexInPartition[2] );
                localIndex nodeOfBox[8];

                if( m_elementType[iR] == "CPE4" || m_elementType[iR] == "STRI" )
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;
                }
                else
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                  nodeOfBox[4] = firstNodeIndex + 1;
                  nodeOfBox[5] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] + numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[7] = numNodesInDir[2] + firstNodeIndex + 1;

                  //               7___________________ 6
                  //               /                   /|
                  //              /                   / |
                  //             /                   /  |
                  //           4/__________________5/   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |                   |   |
                  //            |   3               |   /2        z
                  //            |                   |  /          |   y
                  //            |                   | /           |  /
                  //            |___________________|/            | /
                  //            0                   1             |/____ x

                }


                // Fix local connectivity for single theta (y) partition (radial meshes only)

                setConnectivityForPeriodicBoundaries( globalIJK,
                                                      numNodesInDir,
                                                      firstElemIndexInPartition,
                                                      nodeOfBox );

                for( int iEle = 0; iEle < m_numElePerBox[iR]; ++iEle )
                {
                  localIndex & localElemIndex = localElemIndexInRegion[ m_regionNames[ regionOffset ] ];
                  elemLocalToGlobal[localElemIndex] = elemGlobalIndex( globalIJK ) * m_numElePerBox[iR] + iEle;

                  getElemToNodesRelationInBox( m_elementType[iR],
                                               globalIJK,
                                               iEle,
                                               nodeIDInBox,
                                               numNodesPerElem );

                  for( localIndex iN = 0; iN < numNodesPerElem; ++iN )
                  {
                    elemsToNodes[localElemIndex][iN] = nodeOfBox[nodeIDInBox[iN]];
                  }
                  ++localElemIndex;
                }
              }
            }
          }
        }
      }
    }
  }

  // Node perturbation
  if( m_fPerturb > 0 )
  {

    for( localIndex iN = 0; iN != nodeManager.size(); ++iN )
    {

      for( int i = 0; i < m_dim; ++i )
      {
        if( X[iN][i] > m_min[i] && X[iN][i] < m_max[i] )
        {
          // This ensures that the perturbation pattern is unaffected by domain
          srand( LvArray::integerConversion< int >( nodeLocalToGlobal[iN] ) + m_randSeed + i );

          X[iN][i] += ( ( m_max[i] - m_min[i] ) / m_numElemsTotal[i] ) * ( ( rand() * 1.0 ) / RAND_MAX - 0.5 ) * 2 * m_fPerturb;
        }
      }
    }
  }

  if( std::fabs( m_skewAngle ) > 0.0 )
  {
    for( localIndex iN = 0; iN != nodeManager.size(); ++iN )
    {
      X[iN][0] -= ( X[iN][1] - m_skewCenter[1] ) * std::tan( m_skewAngle );
    }
  }

  coordinateTransformation( nodeManager );

  GEOSX_LOG_RANK_0( "Total number of nodes:"<<(m_numElemsTotal[0]+1)*(m_numElemsTotal[1]+1)*(m_numElemsTotal[2]+1) );
  GEOSX_LOG_RANK_0( "Total number of elems:"<<m_numElemsTotal[0]*m_numElemsTotal[1]*m_numElemsTotal[2] );

  GEOSX_LOG_RANK( "Total number of nodes:"<<nodeManager.size() );

}

/**
 * @param elementType
 * @param index
 * @param iEle
 * @param nodeIDInBox
 * @param node_size
 */
void InternalMeshGenerator::getElemToNodesRelationInBox( const string & elementType,
                                                         int const (&index)[3],
                                                         int const & iEle,
                                                         int (& nodeIDInBox)[8],
                                                         int const node_size )

{
  if( elementType == "C3D8" )
  {
    nodeIDInBox[0] = 0;
    nodeIDInBox[1] = 1;
    nodeIDInBox[2] = 3;
    nodeIDInBox[3] = 2;
    nodeIDInBox[4] = 4;
    nodeIDInBox[5] = 5;
    nodeIDInBox[6] = 7;
    nodeIDInBox[7] = 6;
  }
  else if( ( elementType == "C3D6" ) && ( m_trianglePattern == 0 ) )
  {
    if( ( index[0] + index[1] ) % 2 == 1 )
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 6;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
    }
    else
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 7;
        nodeIDInBox[3] = 3;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
    }
  }
  else if( ( elementType == "C3D6" ) && ( m_trianglePattern == 1 ) )
  {
    if( index[1] % 2 == 0 )
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
      else
      {
        nodeIDInBox[0] = 6;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
    }
    else
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 4;
        nodeIDInBox[3] = 0;
        nodeIDInBox[4] = 7;
        nodeIDInBox[5] = 3;
        nodeIDInBox[6] = 7;
        nodeIDInBox[7] = 3;
      }
      else
      {
        nodeIDInBox[0] = 5;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 7;
        nodeIDInBox[3] = 3;
        nodeIDInBox[4] = 6;
        nodeIDInBox[5] = 2;
        nodeIDInBox[6] = 6;
        nodeIDInBox[7] = 2;
      }
    }
  }
  else if( elementType == "CPE4" )
  {
    nodeIDInBox[0] = 0;
    nodeIDInBox[1] = 1;
    nodeIDInBox[2] = 3;
    nodeIDInBox[3] = 2;
  }
  else if( ( elementType == "STRI" ) && ( m_trianglePattern == 0 ) )
  {
    if( ( index[0] + index[1] ) % 2 == 1 )
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 2;
      }
    }
    else
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 1;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }
    }
  }
  else if( ( elementType == "STRI" ) && ( m_trianglePattern == 1 ) )
  {
    if( index[1] % 2 == 0 )
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 2;
      }
      else
      {
        nodeIDInBox[0] = 2;
        nodeIDInBox[1] = 3;
        nodeIDInBox[2] = 0;
      }
    }
    else
    {
      if( iEle % 2 == 0 )
      {
        nodeIDInBox[0] = 0;
        nodeIDInBox[1] = 1;
        nodeIDInBox[2] = 3;
      }
      else
      {
        nodeIDInBox[0] = 1;
        nodeIDInBox[1] = 2;
        nodeIDInBox[2] = 3;
      }
    }
  }
  else if( elementType == "C3D4" )
  {
    int mapBoxTet[8][6][4] =
    {
      {
        { 0, 3, 7, 6 },
        { 0, 7, 4, 6 },
        { 0, 5, 1, 6 },
        { 0, 4, 5, 6 },
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 }
      },
      {
        { 1, 5, 6, 7 },
        { 1, 6, 2, 7 },
        { 0, 4, 1, 7 },
        { 1, 4, 5, 7 },
        { 0, 1, 3, 7 },
        { 1, 2, 3, 7 }
      },
      {
        { 0, 3, 4, 2 },
        { 3, 7, 4, 2 },
        { 0, 4, 1, 2 },
        { 1, 4, 5, 2 },
        { 4, 7, 6, 2 },
        { 4, 6, 5, 2 }
      },
      {
        { 1, 5, 2, 3 },
        { 2, 5, 6, 3 },
        { 0, 5, 1, 3 },
        { 0, 4, 5, 3 },
        { 4, 7, 5, 3 },
        { 5, 7, 6, 3 }
      },
      {
        { 0, 3, 4, 5 },
        { 3, 7, 4, 5 },
        { 3, 6, 7, 5 },
        { 3, 2, 6, 5 },
        { 3, 0, 1, 5 },
        { 1, 2, 3, 5 }
      },
      {
        { 1, 5, 2, 4 },
        { 2, 5, 6, 4 },
        { 2, 6, 7, 4 },
        { 2, 7, 3, 4 },
        { 0, 2, 3, 4 },
        { 0, 1, 2, 4 }
      },
      {
        { 0, 7, 4, 1 },
        { 0, 3, 7, 1 },
        { 2, 7, 3, 1 },
        { 2, 6, 7, 1 },
        { 4, 7, 5, 1 },
        { 5, 7, 6, 1 }
      },
      {
        { 1, 5, 6, 0 },
        { 1, 6, 2, 0 },
        { 2, 6, 3, 0 },
        { 3, 6, 7, 0 },
        { 4, 6, 5, 0 },
        { 4, 7, 6, 0 }
      }
    };

    int mapBoxType[2][2][2];
    mapBoxType[0][0][0] = 0;
    mapBoxType[1][0][0] = 1;
    mapBoxType[0][0][1] = 2;
    mapBoxType[1][0][1] = 3;
    mapBoxType[0][1][0] = 4;
    mapBoxType[1][1][0] = 5;
    mapBoxType[0][1][1] = 6;
    mapBoxType[1][1][1] = 7;

    int boxType = mapBoxType[index[0] % 2][index[1] % 2][index[2] % 2];
    for( int i = 0; i < node_size; ++i )
    {
      nodeIDInBox[i] = mapBoxTet[boxType][iEle][i];
    }
  }
}

void
InternalMeshGenerator::
  setConnectivityForPeriodicBoundary( int const component,
                                      int const (&globalIJK)[3],
                                      integer const (&numNodesInDir)[3],
                                      int const (&firstElemIndexInPartition)[3],
                                      localIndex (& nodeOfBox)[8] )
{
  // Condition is:
  // 1) element is last index in component direction
  // 2) first local element in component partition is zero
  if( ( globalIJK[component] == m_numElemsTotal[component] - 1 ) &&
      ( firstElemIndexInPartition[component] == 0) )
  {
    // Last set of nodes
    int modGlobalIJK[3] = { globalIJK[0], globalIJK[1], globalIJK[2] };
    modGlobalIJK[component] = 0;
    const localIndex firstNodeIndex = numNodesInDir[1] * numNodesInDir[2] * ( modGlobalIJK[0] - firstElemIndexInPartition[0] )
                                      + numNodesInDir[2] * ( modGlobalIJK[1] - 0 )
                                      + ( modGlobalIJK[2] - firstElemIndexInPartition[2] );

    nodeOfBox[3] = firstNodeIndex;
    nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
    nodeOfBox[7] = firstNodeIndex + 1;
    nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex + 1;
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase, InternalMeshGenerator, string const &, Group * const )
} /* namespace geosx */
