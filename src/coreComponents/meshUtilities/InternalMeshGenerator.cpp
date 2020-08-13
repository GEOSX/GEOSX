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
 * @file InternalMeshGenerator.cpp
 */

#include "InternalMeshGenerator.hpp"

#include "managers/DomainPartition.hpp"

#include "codingUtilities/StringUtilities.hpp"
#include <math.h>
#include <algorithm>

#include "mpiCommunications/PartitionBase.hpp"
#include "mpiCommunications/SpatialPartition.hpp"
#include "common/DataTypes.hpp"

#include "mesh/MeshBody.hpp"

#include "common/TimingMacros.hpp"

namespace geosx
{
using namespace dataRepository;

InternalMeshGenerator::InternalMeshGenerator( string const & name,
                                              Group * const parent ) :
  MeshGeneratorBase( name, parent ),
  //    m_vertices({this->registerWrapper<real64_array>(keys::xCoords).reference(),
  //                this->registerWrapper<real64_array>(keys::yCoords).reference(),
  //                this->registerWrapper<real64_array>(keys::zCoords).reference()
  // }),
  m_dim( 0 ),
  m_min(),
  m_max()
{
  /*
     for( int i=0 ; i<3 ; ++i )
     {
     m_wExtensionMin[i] = 0;
     m_wExtensionMax[i] = 0;
     m_nExtensionLayersMin[i] = 0;
     m_nExtensionLayersMax[i] = 0;
     m_commonRatioMin[i] = 1.5;
     m_commonRatioMax[i] = 1.5;
     }
   */
  m_dim = 3;

  registerWrapper( keys::xCoords, &( m_vertices[0] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "x-coordinates of each mesh block vertex" );

  registerWrapper( keys::yCoords, &( m_vertices[1] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "y-coordinates of each mesh block vertex" );

  registerWrapper( keys::zCoords, &( m_vertices[2] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "z-coordinates of each mesh block vertex" );

  registerWrapper( keys::xElems, &( m_nElems[0] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription(
      "number of elements in the x-direction within each mesh block" );

  registerWrapper( keys::yElems, &( m_nElems[1] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription(
      "number of elements in the y-direction within each mesh block" );

  registerWrapper( keys::zElems, &( m_nElems[2] ) )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription(
      "number of elements in the z-direction within each mesh block" );

  registerWrapper( keys::xBias, &( m_nElemBias[0] ) )
    ->setApplyDefaultValue( 1.0 )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "bias of element sizes in the x-direction within each mesh block "
      "(dx_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( keys::yBias, &( m_nElemBias[1] ) )
    ->setApplyDefaultValue( 1.0 )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "bias of element sizes in the y-direction within each mesh block "
      "(dy_left=(1+b)*L/N, dx_right=(1-b)*L/N)" );

  registerWrapper( keys::zBias, &( m_nElemBias[2] ) )
    ->setApplyDefaultValue( 1.0 )
    ->setSizedFromParent( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "bias of element sizes in the z-direction within each mesh block "
      "(dz_left=(1+b)*L/N, dz_right=(1-b)*L/N)" );

  registerWrapper( keys::cellBlockNames, &m_regionNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "names of each mesh block" );

  registerWrapper( keys::elementTypes, &m_elementType )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "element types of each mesh block" );

  registerWrapper( keys::trianglePattern, &m_trianglePattern )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "pattern by which to decompose the hex mesh into prisms (more "
      "explanation required)" );
}

InternalMeshGenerator::~InternalMeshGenerator()
{
  // TODO Auto-generated destructor stub
}

/**
 * @param domain
 */
void
InternalMeshGenerator::GenerateElementRegions(
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  //  lvector numElements;
  //
  //  for( string_array::size_type r=0 ; r<m_regionNames.size() ; ++r )
  //  {
  //    numElements.emplace_back( 0 );
  //  }
  //
  //  domain.m_feElementManager->resize( numElements, m_regionNames,
  // m_elementType );
}

void
InternalMeshGenerator::PostProcessInput()
{
  if( m_elementType[0] == "C3D8" || m_elementType[0] == "C3D4" ||
      m_elementType[0] == "C3D6" )
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
    for( int i = 0; i < m_dim; ++i )
    {
      failFlag += ( m_nElems[i].size() != m_vertices[i].size() - 1 );
    }
    if( failFlag )
    {
      GEOSX_ERROR(
        "vertex/element mismatch InternalMeshGenerator::ReadXMLPost()" );
    }

    // If specified, check to make sure bias values have the correct length
    for( int i = 0; i < m_dim; ++i )
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

  m_numElePerBox.resize( m_nElems[0].size() * m_nElems[1].size() *
                         m_nElems[2].size() );

  if( LvArray::integerConversion< long >( m_elementType.size() ) !=
      m_numElePerBox.size() )
  {
    if( m_elementType.size() == 1 )
    {
      m_elementType.resize( m_numElePerBox.size() );
      for( localIndex i = 1; i < m_elementType.size(); ++i )
      {
        m_elementType[i] = m_elementType[0];
      }
    }
    else
    {
      GEOSX_ERROR(
        "InternalMeshGenerator: The number of element types is inconsistent "
        "with the number of total block." );
    }
  }

  for( localIndex i = 0;
       i < LvArray::integerConversion< localIndex >( m_elementType.size() );
       ++i )
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
    for( int i = 0; i < m_dim; ++i )
    {
      numBlocks *= m_nElems[i].size();
    }
    if( numBlocks != m_regionNames.size() )
    {
      if( m_regionNames.size() == 1 )
      {
        m_regionNames.resize( numBlocks );
        for( localIndex i = 1; i < m_elementType.size(); ++i )
        {
          m_regionNames[i] = m_regionNames[0];
        }
      }
      else
      {
        GEOSX_ERROR(
          "Incorrect number of regionLayout entries specified in "
          "InternalMeshGenerator::ReadXML()" );
      }
    }
  }

  for( int i = 0; i < 3; ++i )
  {
    m_min[i] = m_vertices[i].front();
    m_max[i] = m_vertices[i].back();
  }

  for( int dir = 0; dir < 3; ++dir )
  {
    m_firstElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_lastElemIndexForBlock[dir].resize( m_nElems[dir].size() );
    m_firstElemIndexForBlock[dir][0] = 0;
    m_lastElemIndexForBlock[dir][0] = m_nElems[dir][0] - 1;
    for( int block = 1; block < m_nElems[dir].size(); ++block )
    {
      m_firstElemIndexForBlock[dir][block] =
        m_lastElemIndexForBlock[dir][block - 1] + 1;
      m_lastElemIndexForBlock[dir][block] =
        m_firstElemIndexForBlock[dir][block] + m_nElems[dir][block] - 1;
    }
  }

  m_fPerturb = 0.0;

  //    m_fPerturb = hdn.GetAttributeOrDefault<realT>("perturbationFactor", 0.0);
  //    m_randSeed = hdn.GetAttributeOrDefault<int>("perturbationSeed",
  // time(NULL));
  //    srand(m_randSeed);
  //
  //    m_mapToRadial = hdn.GetAttributeOrDefault<int>("mapToRadial", 0);
  //
  //    m_skewAngle = hdn.GetAttributeOrDefault<realT>("skewAngle", 0.0);
  //    m_skewAngle *= 3.14159265/180;
  //    R1Tensor zeroVector;
  //    zeroVector *= 0.0;
  //    m_skewCenter = hdn.GetAttributeOrDefault<R1Tensor>("skewCenter",
  // zeroVector);
  //
  //
  //    // Mesh deformation
  //    m_delayMeshDeformation =
  // hdn.GetAttributeOrDefault<int>("delayMeshDeformation", 0);
  //    m_meshDx = hdn.GetAttributeString("dxTable");
  //    m_meshDy = hdn.GetAttributeString("dyTable");
  //    m_meshDz = hdn.GetAttributeString("dzTable");
}

Group *
InternalMeshGenerator::CreateChild(
  string const & GEOSX_UNUSED_PARAM( childKey ),
  string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/**
 * @param partition
 * @param domain
 */
void
InternalMeshGenerator::GenerateMesh( DomainPartition * const domain )
{
  GEOSX_MARK_FUNCTION;

  // This cannot find groupkeys:
  // Group * const meshBodies = domain->GetGroup(domain->groupKeys.meshBodies);
  Group * const meshBodies = domain->GetGroup( std::string( "MeshBodies" ) );
  MeshBody * const meshBody = meshBodies->RegisterGroup< MeshBody >( this->getName() );
  MeshLevel * const meshLevel0 =
    meshBody->RegisterGroup< MeshLevel >( std::string( "Level0" ) );

  // special case
  //  bool isRadialWithOneThetaPartition = (m_mapToRadial > 0) &&
  // (partition.GetPartitions()[1]==1);

  NodeManager * nodeManager = meshLevel0->getNodeManager();

  // Make sure that the node manager fields are initialized

  CellBlockManager * elementManager =
    domain->GetGroup< CellBlockManager >( keys::cellManager );
  Group & nodeSets = nodeManager->sets();

  PartitionBase & partition =
    domain->getReference< PartitionBase >( keys::partitionManager );

  bool isRadialWithOneThetaPartition = false;

  // This should probably handled elsewhere:
  int aa = 0;
  for( auto & cellBlockName : m_regionNames )
  {
    CellBlock * cellBlock = elementManager->GetGroup( keys::cellBlocks )
                              ->RegisterGroup< CellBlock >( cellBlockName );
    string elementType = m_elementType[aa++];
    cellBlock->SetElementType( elementType );
  }

  SortedArray< localIndex > & xnegNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "xneg" ) )
      ->reference();
  SortedArray< localIndex > & xposNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "xpos" ) )
      ->reference();
  SortedArray< localIndex > & ynegNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "yneg" ) )
      ->reference();
  SortedArray< localIndex > & yposNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "ypos" ) )
      ->reference();
  SortedArray< localIndex > & znegNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "zneg" ) )
      ->reference();
  SortedArray< localIndex > & zposNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "zpos" ) )
      ->reference();
  SortedArray< localIndex > & allNodes =
    nodeSets.registerWrapper< SortedArray< localIndex > >( std::string( "all" ) )
      ->reference();

  // partition based on even spacing to get load balance
  // Partition geometrical boundaries will be corrected in the end.
  {
    m_min[0] = m_vertices[0].front();
    m_min[1] = m_vertices[1].front();
    m_min[2] = m_vertices[2].front();

    m_max[0] = m_vertices[0].back();
    m_max[1] = m_vertices[1].back();
    m_max[2] = m_vertices[2].back();

    R1Tensor temp1( m_min );
    R1Tensor temp2( m_max );

    partition.setSizes( temp1, temp2 );
    temp2 -= temp1;
    meshBody->setGlobalLengthScale( std::fabs( temp2.L2_Norm() ) );
  }

  // find elemCenters for even uniform element sizes
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
      elemCenterCoordsLocal[k] =
        m_min[i] + ( m_max[i] - m_min[i] ) * ( k + 0.5 ) / m_numElemsTotal[i];
    }
    MpiWrapper::allReduce( elemCenterCoordsLocal.data(),
                           elemCenterCoords[i].data(),
                           m_numElemsTotal[i],
                           MPI_MAX,
                           MPI_COMM_GEOSX );
  }

  // find starting/ending index

  // get the first and last indices in this partition each direction
  int firstElemIndexInPartition[3] = { -1, -1, -1 };
  int lastElemIndexInPartition[3] = { -2, -2, -2 };

  for( int i = 0; i < 3; ++i )
  {
    //    firstElemIndexInPartition[i] = -1;
    //    lastElemIndexInPartition[i] = -2;
    for( int k = 0; k < m_numElemsTotal[i]; ++k )
    {
      if( partition.IsCoordInPartition( elemCenterCoords[i][k], i ) )
      {
        firstElemIndexInPartition[i] = k;
        break;
      }
    }

    if( firstElemIndexInPartition[i] > -1 )
    {
      for( int k = firstElemIndexInPartition[i]; k < m_numElemsTotal[i]; ++k )
      {
        if( partition.IsCoordInPartition( elemCenterCoords[i][k], i ) )
        {
          lastElemIndexInPartition[i] = k;
        }
      }
    }
  }

  // calculate number of elements in this partition from each region, and the
  // total number of nodes

  std::map< std::string, int > numElemsInRegions;
  std::map< std::string, std::string > elemTypeInRegions;

  integer_array firstElemIndexForBlockInPartition[3];
  integer_array lastElemIndexForBlockInPartition[3];

  for( int dir = 0; dir < 3; ++dir )
  {
    firstElemIndexForBlockInPartition[dir] = m_firstElemIndexForBlock[dir];
    lastElemIndexForBlockInPartition[dir] = m_lastElemIndexForBlock[dir];

    for( int block = 0; block < m_nElems[dir].size(); ++block )
    {
      if( firstElemIndexForBlockInPartition[dir][block] >
            lastElemIndexInPartition[dir] ||
          lastElemIndexForBlockInPartition[dir][block] <
            firstElemIndexInPartition[dir] )
      {
        firstElemIndexForBlockInPartition[dir][block] = -1;
        lastElemIndexForBlockInPartition[dir][block] = -2;
      }
      else
      {
        if( firstElemIndexForBlockInPartition[dir][block] <
            firstElemIndexInPartition[dir] )
        {
          firstElemIndexForBlockInPartition[dir][block] =
            firstElemIndexInPartition[dir];
        }
        if( lastElemIndexForBlockInPartition[dir][block] >
            lastElemIndexInPartition[dir] )
        {
          lastElemIndexForBlockInPartition[dir][block] =
            lastElemIndexInPartition[dir];
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
        numElemsInRegions[m_regionNames[regionOffset]] = 0;
        elemTypeInRegions[m_regionNames[regionOffset]] = "";
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
        for( int kblock = 0; kblock < m_nElems[2].size();
             ++kblock, ++regionOffset, ++iR )
        {
          int numElemsInRegion = 1;
          numElemsInRegion *= lastElemIndexForBlockInPartition[0][iblock] -
            firstElemIndexForBlockInPartition[0][iblock] + 1;

          if( m_dim > 1 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[1][jblock] -
              firstElemIndexForBlockInPartition[1][jblock] + 1;
          }
          if( m_dim > 2 )
          {
            numElemsInRegion *= lastElemIndexForBlockInPartition[2][kblock] -
              firstElemIndexForBlockInPartition[2][kblock] + 1;
          }

          numElemsInRegion *= m_numElePerBox[iR];
          numElemsInRegions[m_regionNames[regionOffset]] += numElemsInRegion;
          elemTypeInRegions[m_regionNames[regionOffset]] = m_elementType[iR];
        }
      }
    }
  }

  localIndex numNodes = 1;

  //  int numElemsInDir[3] = {1,1,1};
  localIndex numNodesInDir[3] = { 1, 1, 1 };

  for( int i = 0; i < m_dim; ++i )
  {
    //    numElemsInDir[i] = lastElemIndexInPartition[i] -
    // firstElemIndexInPartition[i] + 1;
    numNodesInDir[i] =
      lastElemIndexInPartition[i] - firstElemIndexInPartition[i] + 2;
    if( isRadialWithOneThetaPartition && i == 1 )
    {
      numNodesInDir[1] -= 1;
    }
    numNodes *= numNodesInDir[i];
  }

  nodeManager->resize( numNodes );
  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();

  arrayView1d< globalIndex > const & nodeLocalToGlobal =
    nodeManager->localToGlobalMap();

  {
    localIndex localNodeIndex = 0;
    for( int i = 0; i < numNodesInDir[0]; ++i )
    {
      for( int j = 0; j < numNodesInDir[1]; ++j )
      {
        for( int k = 0; k < numNodesInDir[2]; ++k )
        {
          int index[3] = { i, j, k };
          for( int a = 0; a < m_dim; ++a )
          {
            index[a] += firstElemIndexInPartition[a];
          }

          R1Tensor const pos = NodePosition( index, m_trianglePattern );
          for( int a = 0; a < 3; ++a )
          {
            X( localNodeIndex, a ) = pos[a];
          }

          // alter global node map for radial mesh
          if( m_mapToRadial > 0 )
          {
            if( isEqual( X( localNodeIndex, 1 ), m_max[1], 1e-10 ) )
            {
              index[1] = 0;
            }
          }

          nodeLocalToGlobal[localNodeIndex] = NodeGlobalIndex( index );

          // cartesian-specific nodesets
          if( m_mapToRadial == 0 )
          {
            if( isEqual( X( localNodeIndex, 0 ), m_min[0], 1e-10 ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 0 ), m_max[0], 1e-10 ) )
            {
              xposNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_min[1], 1e-10 ) )
            {
              ynegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 1 ), m_max[1], 1e-10 ) )
            {
              yposNodes.insert( localNodeIndex );
            }
          }
          else
          {
            // radial-specific nodesets
            if( isEqual( X( localNodeIndex, 0 ), m_min[0], 1e-10 ) )
            {
              xnegNodes.insert( localNodeIndex );
            }
            if( isEqual( X( localNodeIndex, 0 ), m_max[0], 1e-10 ) )
            {
              xposNodes.insert( localNodeIndex );
            }
          }

          // general nodesets
          if( isEqual( X( localNodeIndex, 2 ), m_min[2], 1e-10 ) )
          {
            znegNodes.insert( localNodeIndex );
          }
          if( isEqual( X( localNodeIndex, 2 ), m_max[2], 1e-10 ) )
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
    integer_array numElements;
    string_array elementRegionNames;
    string_array elementTypes;
    std::map< std::string, localIndex > localElemIndexInRegion;

    for( std::map< std::string, int >::iterator iterNumElemsInRegion =
           numElemsInRegions.begin();
         iterNumElemsInRegion != numElemsInRegions.end();
         ++iterNumElemsInRegion )
    {
      numElements.emplace_back( iterNumElemsInRegion->second );
      elementRegionNames.emplace_back( iterNumElemsInRegion->first );
      elementTypes.emplace_back( elemTypeInRegions[iterNumElemsInRegion->first] );

      localElemIndexInRegion[iterNumElemsInRegion->first] = 0;
    }

    elementManager->resize( numElements, elementRegionNames, elementTypes );

    // assign global numbers to elements
    regionOffset = 0;
    SortedArray< std::string > processedRegionNames;
    localIndex iR = 0;

    for( int iblock = 0; iblock < m_nElems[0].size(); ++iblock )
    {
      for( int jblock = 0; jblock < m_nElems[1].size(); ++jblock )
      {
        for( int kblock = 0; kblock < m_nElems[2].size();
             ++kblock, ++regionOffset, ++iR )
        {
          //          ElementRegionT& elemRegion =
          // domain->m_feElementManager->m_ElementRegions[*iterRegion];

          CellBlock * elemRegion =
            elementManager->GetRegion( m_regionNames[regionOffset] );
          int const numNodesPerElem =
            LvArray::integerConversion< int >( elemRegion->numNodesPerElement() );
          integer_array nodeIDInBox( 8 );

          arrayView2d< localIndex, cells::NODE_MAP_USD > elemsToNodes =
            elemRegion->nodeList();
          arrayView1d< globalIndex > const & elemLocalToGlobal =
            elemRegion->localToGlobalMap();

          int numElemsInDirForRegion[3] = {
            lastElemIndexForBlockInPartition[0][iblock] -
              firstElemIndexForBlockInPartition[0][iblock] + 1,
            lastElemIndexForBlockInPartition[1][jblock] -
              firstElemIndexForBlockInPartition[1][jblock] + 1,
            lastElemIndexForBlockInPartition[2][kblock] -
              firstElemIndexForBlockInPartition[2][kblock] + 1 };

          for( int i = 0; i < numElemsInDirForRegion[0]; ++i )
          {
            for( int j = 0; j < numElemsInDirForRegion[1]; ++j )
            {
              for( int k = 0; k < numElemsInDirForRegion[2]; ++k )
              {
                int index[3] = { i + firstElemIndexForBlockInPartition[0][iblock],
                                 j + firstElemIndexForBlockInPartition[1][jblock],
                                 k + firstElemIndexForBlockInPartition[2][kblock] };

                const localIndex firstNodeIndex = numNodesInDir[1] *
                    numNodesInDir[2] * ( index[0] - firstElemIndexInPartition[0] ) +
                  numNodesInDir[2] * ( index[1] - firstElemIndexInPartition[1] ) +
                  ( index[2] - firstElemIndexInPartition[2] );
                localIndex nodeOfBox[8];

                if( m_elementType[iR] == "CPE4" || m_elementType[iR] == "STRI" )
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] =
                    numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] +
                    numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;
                }
                else
                {
                  nodeOfBox[0] = firstNodeIndex;
                  nodeOfBox[1] =
                    numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] +
                    numNodesInDir[2] + firstNodeIndex;
                  nodeOfBox[3] = numNodesInDir[2] + firstNodeIndex;

                  nodeOfBox[4] = firstNodeIndex + 1;
                  nodeOfBox[5] =
                    numNodesInDir[1] * numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] +
                    numNodesInDir[2] + firstNodeIndex + 1;
                  nodeOfBox[7] = numNodesInDir[2] + firstNodeIndex + 1;

                  //               7___________________ 6
                  //               /                   /|
                  //              /                   / |f
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
                // fix local connectivity for single theta (y) partition (radial meshes only)
                if( isRadialWithOneThetaPartition )
                {
                  if( j == numElemsInDirForRegion[1] - 1 &&
                      jblock == m_nElems[1].size() - 1 )
                  {  // last set of elements
                    index[1] = -1;
                    const localIndex firstNodeIndexR = numNodesInDir[1] *
                        numNodesInDir[2] *
                        ( index[0] - firstElemIndexInPartition[0] ) +
                      numNodesInDir[2] *
                        ( index[1] - firstElemIndexInPartition[1] ) +
                      ( index[2] - firstElemIndexInPartition[2] );
                    nodeOfBox[2] = numNodesInDir[1] * numNodesInDir[2] +
                      numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[3] = numNodesInDir[2] + firstNodeIndexR;
                    nodeOfBox[6] = numNodesInDir[1] * numNodesInDir[2] +
                      numNodesInDir[2] + firstNodeIndexR + 1;
                    nodeOfBox[7] = numNodesInDir[2] + firstNodeIndexR + 1;
                  }
                }

                for( int iEle = 0; iEle < m_numElePerBox[iR]; ++iEle )
                {
                  localIndex & localElemIndex =
                    localElemIndexInRegion[m_regionNames[regionOffset]];
                  elemLocalToGlobal[localElemIndex] =
                    ElemGlobalIndex( index ) * m_numElePerBox[iR] + iEle;

                  GetElemToNodesRelationInBox( m_elementType[iR],
                                               index,
                                               iEle,
                                               nodeIDInBox.data(),
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

#if 0
  // Correct partition geometrical boundary.
  {
    R1Tensor pMin, pMax;
    partition.getPartitionGeometricalBoundary( pMin, pMax );
    for( int i = 0; i < m_dim; ++i )
    {
      realT xMinByNumElems = ( pMin[i] - m_min[i] ) / ( m_max[i] - m_min[i] ) * m_numElemsTotal[i];
      realT xMaxByNumElems = ( pMax[i] - m_min[i] ) / ( m_max[i] - m_min[i] ) * m_numElemsTotal[i];

      localIndex iBlockMin( 0 ), iBlockMax( 0 );
      while( ( xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 || xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 + 1.0 )
             && iBlockMin < m_nElems[i].size() - 1 )
      {
        ++iBlockMin;
      }
      while( ( xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 || xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 + 1.0 )
             && iBlockMax < m_nElems[i].size() - 1 )
      {
        ++iBlockMax;
      }
      //sanity check

      if( xMinByNumElems < m_firstElemIndexForBlock[i][iBlockMin] * 1.0 ||
          xMinByNumElems > m_lastElemIndexForBlock[i][iBlockMin] * 1.0 + 1.0 ||
          iBlockMin >= m_nElems[i].size() )
      {
        std::cout
          << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer"
          << std::endl;
      }
      if( xMaxByNumElems < m_firstElemIndexForBlock[i][iBlockMax] * 1.0 ||
          xMaxByNumElems > m_lastElemIndexForBlock[i][iBlockMax] * 1.0 + 1.0 ||
          iBlockMax >= m_nElems[i].size() )
      {
        std::cout
          << "WARNING: Had some trouble in correcting partition geometric boundaries.  If this affects the results, contact a developer"
          << std::endl;
      }

      pMin[i] = m_vertices[i][iBlockMin] + ( xMinByNumElems - m_firstElemIndexForBlock[i][iBlockMin] ) *
                ( m_vertices[i][iBlockMin + 1] - m_vertices[i][iBlockMin] ) / m_nElems[i][iBlockMin];
      pMax[i] = m_vertices[i][iBlockMax] + ( xMaxByNumElems - m_firstElemIndexForBlock[i][iBlockMax] ) *
                ( m_vertices[i][iBlockMax + 1] - m_vertices[i][iBlockMax] ) / m_nElems[i][iBlockMax];
    }

    partition.SetPartitionGeometricalBoundary( pMin, pMax );

  }
#endif

  /*
     {// Move nodes in the extension layers.

     for (localIndex iN = 0; iN != nodeManager->DataLengths(); ++iN)
     {
     for (int i=0; i<m_dim; ++i)
     {
     if ( X[iN][i] < m_min[i])
     {
     int eLayer = (int) ((m_min[i] -X[iN][i]) / ((m_max[i] - m_min[i]) /
        m_numElems[i]) + 0.5);
     X[iN][i] = m_min[i] - ((m_max[i] - m_min[i]) / m_numElems[i]) *
        m_commonRatioMin[i] * (1- pow(m_commonRatioMin[i], eLayer)) / (1 -
        m_commonRatioMin[i]);
     }
     else if (X[iN][i] > m_max[i])
     {
     int eLayer = (int) ((X[iN][i] - m_max[i] ) / ((m_max[i] - m_min[i]) /
        m_numElems[i]) + 0.5);
     X[iN][i] = m_max[i] + ((m_max[i] - m_min[i]) / m_numElems[i]) *
        m_commonRatioMax[i] * (1- pow(m_commonRatioMax[i], eLayer)) / (1 -
        m_commonRatioMax[i]);
     }

     }
     }
     }
   */

  // Node perturbation
  if( m_fPerturb > 0 )
  {
    for( localIndex iN = 0; iN != nodeManager->size(); ++iN )
    {
      for( int i = 0; i < m_dim; ++i )
      {
        if( X[iN][i] > m_min[i] && X[iN][i] < m_max[i] )
        {
          srand( LvArray::integerConversion< int >( nodeLocalToGlobal[iN] ) +
                 m_randSeed + i );  // This
          // ensures
          // that
          // the
          // perturbation
          // pattern
          // is
          // unaffected
          // by
          // domain
          X[iN][i] += ( ( m_max[i] - m_min[i] ) / m_numElemsTotal[i] ) *
            ( ( rand() * 1.0 ) / RAND_MAX - 0.5 ) * 2 * m_fPerturb;
        }
      }
    }
  }

  if( std::fabs( m_skewAngle ) > 0.0 )
  {
    for( localIndex iN = 0; iN != nodeManager->size(); ++iN )
    {
      X[iN][0] -= ( X[iN][1] - m_skewCenter[1] ) * std::tan( m_skewAngle );
    }
  }

  if( m_mapToRadial > 0 )
  {
    // Map to radial mesh
    for( localIndex iN = 0; iN != nodeManager->size(); ++iN )
    {
      m_meshTheta = X[iN][1] * M_PI / 180.0;
      m_meshAxis = static_cast< int >( round( m_meshTheta * 2.0 / M_PI ) );
      m_meshPhi = fabs( m_meshTheta - m_meshAxis * M_PI / 2.0 );
      m_meshRout = m_max[0] / cos( m_meshPhi );

      if( m_mapToRadial > 1 )
      {
        m_meshRact = ( ( m_meshRout - m_min[0] ) / ( m_max[0] - m_min[0] ) ) *
            ( X[iN][0] - m_min[0] ) +
          m_min[0];
      }
      else
      {
        m_meshRact = X[iN][0];
      }

      X[iN][0] = m_meshRact * cos( m_meshTheta );
      X[iN][1] = m_meshRact * sin( m_meshTheta );

      // add mapped values to nodesets
      if( m_mapToRadial > 1 )
      {
        if( isEqual( X[iN][0], -1 * m_max[0], 1e-6 ) )
        {
          xnegNodes.insert( iN );
        }
        if( isEqual( X[iN][0], m_max[0], 1e-6 ) )
        {
          xposNodes.insert( iN );
        }
        if( isEqual( X[iN][1], -1 * m_max[0], 1e-6 ) )
        {
          ynegNodes.insert( iN );
        }
        if( isEqual( X[iN][1], m_max[0], 1e-6 ) )
        {
          yposNodes.insert( iN );
        }
      }
    }
  }

  if( m_delayMeshDeformation == 0 )
  {
    RemapMesh( domain );
  }
}

/**
 * @param elementType
 * @param index
 * @param iEle
 * @param nodeIDInBox
 * @param node_size
 */
void
InternalMeshGenerator::GetElemToNodesRelationInBox( const std::string & elementType,
                                                    const int index[],
                                                    const int & iEle,
                                                    int nodeIDInBox[],
                                                    const int node_size )

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
    int mapBoxTet[8][6][4] = { { { 0, 3, 7, 6 },
                                 { 0, 7, 4, 6 },
                                 { 0, 5, 1, 6 },
                                 { 0, 4, 5, 6 },
                                 { 0, 1, 2, 6 },
                                 { 0, 2, 3, 6 } },
                               { { 1, 5, 6, 7 },
                                 { 1, 6, 2, 7 },
                                 { 0, 4, 1, 7 },
                                 { 1, 4, 5, 7 },
                                 { 0, 1, 3, 7 },
                                 { 1, 2, 3, 7 } },
                               { { 0, 3, 4, 2 },
                                 { 3, 7, 4, 2 },
                                 { 0, 4, 1, 2 },
                                 { 1, 4, 5, 2 },
                                 { 4, 7, 6, 2 },
                                 { 4, 6, 5, 2 } },
                               { { 1, 5, 2, 3 },
                                 { 2, 5, 6, 3 },
                                 { 0, 5, 1, 3 },
                                 { 0, 4, 5, 3 },
                                 { 4, 7, 5, 3 },
                                 { 5, 7, 6, 3 } },
                               { { 0, 3, 4, 5 },
                                 { 3, 7, 4, 5 },
                                 { 3, 6, 7, 5 },
                                 { 3, 2, 6, 5 },
                                 { 3, 0, 1, 5 },
                                 { 1, 2, 3, 5 } },
                               { { 1, 5, 2, 4 },
                                 { 2, 5, 6, 4 },
                                 { 2, 6, 7, 4 },
                                 { 2, 7, 3, 4 },
                                 { 0, 2, 3, 4 },
                                 { 0, 1, 2, 4 } },
                               { { 0, 7, 4, 1 },
                                 { 0, 3, 7, 1 },
                                 { 2, 7, 3, 1 },
                                 { 2, 6, 7, 1 },
                                 { 4, 7, 5, 1 },
                                 { 5, 7, 6, 1 } },
                               { { 1, 5, 6, 0 },
                                 { 1, 6, 2, 0 },
                                 { 2, 6, 3, 0 },
                                 { 3, 6, 7, 0 },
                                 { 4, 6, 5, 0 },
                                 { 4, 7, 6, 0 } } };

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
InternalMeshGenerator::RemapMesh(
  dataRepository::Group * const GEOSX_UNUSED_PARAM( domain ) )
{
  //  // Node mapping
  //  if (!m_meshDx.empty())
  //  {
  //    const Table3D* tableDx =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDx);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dx=tableDx->Lookup(X[iN]);
  //      X[iN][0] += dx;
  //    }
  //  }
  //
  //  if (!m_meshDy.empty())
  //  {
  //    const Table3D* tableDy =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDy);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dy=tableDy->Lookup(X[iN]);
  //      X[iN][1] += dy;
  //    }
  //  }
  //
  //  if (!m_meshDz.empty())
  //  {
  //    const Table3D* tableDz =
  // stlMapLookupPointer(TableManager::Instance().Tables<3>(), m_meshDz);
  //
  //    for (localIndex iN=0; iN!=nodeManager->DataLengths(); ++iN)
  //    {
  //      realT dz=tableDz->Lookup(X[iN]);
  //      X[iN][2] += dz;
  //    }
  //  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase,
                        InternalMeshGenerator,
                        std::string const &,
                        Group * const )
}  // namespace geosx
