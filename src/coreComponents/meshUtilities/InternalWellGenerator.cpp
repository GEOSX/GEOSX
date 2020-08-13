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

/*
 * @file InternalWellGenerator.cpp
 *
 */

#include "InternalWellGenerator.hpp"

#include "mpiCommunications/CommunicationTools.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "mesh/WellElementRegion.hpp"
#include "mesh/WellElementSubRegion.hpp"
#include "meshUtilities/PerforationData.hpp"
#include "meshUtilities/Perforation.hpp"

namespace geosx
{
using namespace dataRepository;

InternalWellGenerator::InternalWellGenerator( string const & name,
                                              Group * const parent ) :
  MeshGeneratorBase( name, parent ),
  m_numElemsPerSegment( 0 ),
  m_radius( 0 ),
  m_wellRegionName( "" ),
  m_wellControlsName( "" ),
  m_meshBodyName( "" ),
  m_numElems( 0 ),
  m_numNodesPerElem( 2 ),
  m_numNodes( 0 ),
  m_numPerforations( 0 ),
  m_nDims( 3 ),
  m_polylineHeadNodeId( -1 )
{
  registerWrapper( keys::nodeCoords, &m_inputPolyNodeCoords )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Physical coordinates of the well polyline nodes" );

  registerWrapper( keys::segmentConn, &m_segmentToPolyNodeMap )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Connectivity of the polyline segments" );

  registerWrapper( keys::radius, &m_radius )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Radius of the well" );

  registerWrapper( keys::nElems, &m_numElemsPerSegment )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Number of well elements per polyline segment" );

  registerWrapper( keys::wellRegionName, &m_wellRegionName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Name of the well element region" );

  registerWrapper( keys::wellControlsName, &m_wellControlsName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription(
      "Name of the set of constraints associated with this well" );

  registerWrapper( keys::meshBodyName, &m_meshBodyName )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setSizedFromParent( 0 )
    ->setDescription( "Name of the reservoir mesh associated with this well" );
}

InternalWellGenerator::~InternalWellGenerator()
{
  // TODO Auto-generated destructor stub
}

void
InternalWellGenerator::PostProcessInput()
{
  GEOSX_ERROR_IF(
    getName().find( "well" ) == std::string::npos,
    "Currently, the well generator must contain the word well in its name " );

  GEOSX_ERROR_IF( m_inputPolyNodeCoords.size( 1 ) != m_nDims,
                  "Invalid number of physical coordinates in "
                    << keys::nodeCoords << " for well " << getName() );

  GEOSX_ERROR_IF(
    m_segmentToPolyNodeMap.size( 1 ) != 2,
    "Invalid size in " << keys::segmentConn << " for well " << getName() );

  GEOSX_ERROR_IF(
    m_inputPolyNodeCoords.size( 0 ) - 1 != m_segmentToPolyNodeMap.size( 0 ),
    "Incompatible sizes of " << keys::nodeCoords << " and " << keys::segmentConn
                             << " in well " << getName() );

  GEOSX_ERROR_IF( m_radius <= 0,
                  "Invalid " << keys::radius << " in well " << getName() );

  GEOSX_ERROR_IF( m_wellRegionName.empty(),
                  "Invalid well region name in well " << getName() );

  GEOSX_ERROR_IF( m_meshBodyName.empty(),
                  "Invalid mesh name in well " << getName() );

  GEOSX_ERROR_IF( m_wellControlsName.empty(),
                  "Invalid well constraint name in well " << getName() );

  // convert the 2D array to an 1D array of R1Tensor
  m_polyNodeCoords.resize( m_inputPolyNodeCoords.size( 0 ) );
  for( globalIndex inode = 0; inode < m_inputPolyNodeCoords.size( 0 ); ++inode )
  {
    R1Tensor coords = { m_inputPolyNodeCoords[inode][0],
                        m_inputPolyNodeCoords[inode][1],
                        m_inputPolyNodeCoords[inode][2] };
    m_polyNodeCoords[inode] = coords;
  }

  // TODO: add more checks here
  // TODO: check that the connectivity of the well is valid
  // TODO: check that with no branching we can go from top to bottom and touch all the elements
}

Group *
InternalWellGenerator::CreateChild( string const & childKey,
                                    string const & childName )
{
  if( childKey == keys::perforation )
  {
    ++m_numPerforations;

    // keep track of the perforations that have been added
    m_perforationList.emplace_back( childName );

    return RegisterGroup< Perforation >( childName );
  }
  else
  {
    GEOSX_ERROR( "Unrecognized node: " << childKey );
  }
  return nullptr;
}

void
InternalWellGenerator::ExpandObjectCatalogs()
{
  CreateChild( keys::perforation, keys::perforation );
}

void
InternalWellGenerator::GenerateMesh( DomainPartition * const domain )
{
  // count the number of well elements to create
  m_numElems = m_numElemsPerSegment * m_segmentToPolyNodeMap.size( 0 );
  m_numNodes = m_numElems + 1;

  // resize the well element, node, and perforation arrays
  m_elemCenterCoords.resize( m_numElems );
  m_nextElemId.resize( m_numElems );
  m_prevElemId.resize( m_numElems );
  m_elemToNodesMap.resizeDimension< 0 >( m_numElems );
  m_elemToNodesMap.resizeDimension< 1 >( m_numNodesPerElem );
  m_elemVolume.resize( m_numElems );

  m_nodeDistFromHead.resize( m_numNodes );
  m_nodeCoords.resize( m_numNodes );

  m_perfCoords.resize( m_numPerforations );
  m_perfDistFromHead.resize( m_numPerforations );
  m_perfTransmissibility.resize( m_numPerforations );
  m_perfElemId.resize( m_numPerforations );

  // construct a reverse map from the polyline nodes to the segments
  ConstructPolylineNodeToSegmentMap();

  // detect the head polyline node based on depth
  FindPolylineHeadNodeIndex();

  // compute the location and distance from top of the well elements
  DiscretizePolyline();

  // map the perforations to the well elements
  ConnectPerforationsToWellElements();

  // merge perforations to make sure that no well element is shared between two MPI domains
  // TODO: instead of merging perforations, split the well elements and do not change the physical location of the
  // perforation
  int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  if( mpiSize > 1 )
  {
    MergePerforations();
  }

  // get the element (sub) region to populate and save the well generator and constraints names
  MeshLevel * const meshLevel =
    domain->getMeshBodies()->GetGroup< MeshBody >( 0 )->getMeshLevel( 0 );

  ElementRegionManager * const elemManager = meshLevel->getElemManager();
  WellElementRegion * const wellRegion =
    elemManager
      ->GetGroup( ElementRegionManager::groupKeyStruct::elementRegionsGroup )
      ->GetGroup< WellElementRegion >( this->m_wellRegionName );

  GEOSX_ERROR_IF( wellRegion == nullptr,
                  "Well region " << this->m_wellRegionName
                                 << " not found in well " << getName() );

  wellRegion->SetWellGeneratorName( this->getName() );
  wellRegion->SetWellControlsName( m_wellControlsName );
}

void
InternalWellGenerator::ConstructPolylineNodeToSegmentMap()
{
  m_polyNodeToSegmentMap.resize( m_polyNodeCoords.size() );

  // loop over the segments
  for( globalIndex iseg = 0; iseg < m_segmentToPolyNodeMap.size( 0 ); ++iseg )
  {
    globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[iseg][0];
    globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[iseg][1];

    R1Tensor v = m_polyNodeCoords[ipolyNode_a];
    v -= m_polyNodeCoords[ipolyNode_b];

    // various checks and warnings on the segment and element length
    GEOSX_ERROR_IF(
      v.L2_Norm() < 1e-2,
      "Error in the topology of well "
        << getName()
        << ": we detected a polyline segment measuring less than 1 cm" );

    if( v.L2_Norm() < 1 )
    {
      GEOSX_LOG_RANK_0(
        "Warning: we detected a segment measuring less than 1 m in the "
        "topology of well "
        << getName() );
    }

    if( v.L2_Norm() / static_cast< real64 >( m_numElemsPerSegment ) < 1e-2 )
    {
      GEOSX_LOG_RANK_0(
        "Warning: the chosen number of elements per polyline segment ("
        << m_numElemsPerSegment << ") leads to well elements measuring less than 1 cm in the topology of well "
        << getName() );
    }

    // map the polyline node ids to the polyline segment ids
    m_polyNodeToSegmentMap[ipolyNode_a].insert( iseg );
    m_polyNodeToSegmentMap[ipolyNode_b].insert( iseg );
  }
}

void
InternalWellGenerator::FindPolylineHeadNodeIndex()
{
  // we assume that the first *polyline segment* in the XML input is the well head segment
  globalIndex const wellHeadSegId = 0;

  // get the corresponding node indices
  globalIndex const ipolyNode_a = m_segmentToPolyNodeMap[wellHeadSegId][0];
  globalIndex const ipolyNode_b = m_segmentToPolyNodeMap[wellHeadSegId][1];

  // we determine which node is the head node based on depth
  // therefore here we throw an error if the well head segment is horizontal
  GEOSX_ERROR_IF(
    !( m_polyNodeCoords[ipolyNode_a][2] < m_polyNodeCoords[ipolyNode_b][2] ) &&
      !( m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2] ),
    "The head polyline segment cannot be horizontal in well "
      << getName() << " since we use depth to determine which of its nodes is to head node of the well" );

  // detect the top node, assuming z oriented upwards
  m_polylineHeadNodeId =
    ( m_polyNodeCoords[ipolyNode_a][2] > m_polyNodeCoords[ipolyNode_b][2] )
    ? ipolyNode_a
    : ipolyNode_b;

  real64 const headZcoord = m_polyNodeCoords[m_polylineHeadNodeId][2];
  for( globalIndex inode = 0; inode < m_polyNodeCoords.size(); ++inode )
  {
    if( inode == m_polylineHeadNodeId )
    {
      continue;
    }
    real64 const currentZcoord = m_polyNodeCoords[inode][2];

    GEOSX_ERROR_IF(
      !( currentZcoord < headZcoord ),
      "Error in the topology of well "
        << getName()
        << " since we found a well node that is above the head node" );
  }
}

void
InternalWellGenerator::DiscretizePolyline()
{
  // initialize well elements and node ids
  globalIndex ipolyNodeTop = m_polylineHeadNodeId;
  globalIndex iwelemCurrent = 0;
  globalIndex isegCurrent = 0;

  // set the location of the first well node and distance from well head
  m_nodeCoords[0] = m_polyNodeCoords[ipolyNodeTop];
  m_nodeDistFromHead[0] = 0.0;

  // note: this part of the code does not support well branching
  // TODO: check that there is only one branch
  // TODO: read wells with branching (already supported elsewhere in the code)

  // go through the well from top to bottom
  for( globalIndex is = 0; is < m_segmentToPolyNodeMap.size( 0 ); ++is )
  {
    GEOSX_ERROR_IF( isegCurrent == -1,
                    "Invalid segmentToNode map in well " << getName() );

    globalIndex const ipolyNodeBottom =
      ( ipolyNodeTop == m_segmentToPolyNodeMap[isegCurrent][0] )
      ? m_segmentToPolyNodeMap[isegCurrent][1]
      : m_segmentToPolyNodeMap[isegCurrent][0];

    R1Tensor vPoly = m_polyNodeCoords[ipolyNodeBottom];
    vPoly -= m_polyNodeCoords[ipolyNodeTop];

    // add the well elements and well nodes corresponding to this polyline segment
    for( localIndex iw = 0; iw < m_numElemsPerSegment; ++iw )
    {
      // 1) set the element location, connectivity
      real64 const scaleCenter =
        ( iw + 0.5 ) / static_cast< real64 >( m_numElemsPerSegment );
      m_elemCenterCoords[iwelemCurrent] = vPoly;
      m_elemCenterCoords[iwelemCurrent] *= scaleCenter;
      m_elemCenterCoords[iwelemCurrent] += m_polyNodeCoords[ipolyNodeTop];

      GEOSX_ERROR_IF( iwelemCurrent >= m_numElems,
                      "Invalid well topology in well " << getName() );

      globalIndex const iwelemTop = iwelemCurrent - 1;
      globalIndex const iwelemBottom = iwelemCurrent + 1;
      m_nextElemId[iwelemCurrent] = iwelemTop;
      m_prevElemId[iwelemCurrent].resize( 1 );
      m_prevElemId[iwelemCurrent][0] =
        ( iwelemBottom < m_numElems ) ? iwelemBottom : -1;

      // 2) set the node properties
      globalIndex const iwellNodeTop = iwelemCurrent;
      globalIndex const iwellNodeBottom = iwelemCurrent + 1;
      m_elemToNodesMap[iwelemCurrent][NodeLocation::TOP] = iwellNodeTop;
      m_elemToNodesMap[iwelemCurrent][NodeLocation::BOTTOM] = iwellNodeBottom;

      real64 const scaleBottom =
        ( iw + 1.0 ) / static_cast< real64 >( m_numElemsPerSegment );
      m_nodeCoords[iwellNodeBottom] = vPoly;
      m_nodeCoords[iwellNodeBottom] *= scaleBottom;
      m_nodeCoords[iwellNodeBottom] += m_polyNodeCoords[ipolyNodeTop];

      // 3) increment the distance from the well head to the bottom of current element
      R1Tensor vWellElem = m_nodeCoords[iwellNodeBottom];
      vWellElem -= m_nodeCoords[iwellNodeTop];
      m_nodeDistFromHead[iwellNodeBottom] = vWellElem.L2_Norm();
      m_nodeDistFromHead[iwellNodeBottom] += m_nodeDistFromHead[iwellNodeTop];

      // 4) set element volume
      m_elemVolume[iwelemCurrent] =
        vWellElem.L2_Norm() * M_PI * m_radius * m_radius;

      // 4) increment the element counter
      ++iwelemCurrent;
    }

    // then consider the next polyline segment
    ipolyNodeTop = ipolyNodeBottom;
    isegCurrent = GetNextSegmentIndex( isegCurrent, ipolyNodeTop );
  }

  // set the previous index for the bottom segment
  m_prevElemId[iwelemCurrent - 1].resize( 1 );
  m_prevElemId[iwelemCurrent - 1][0] = -1;
}

void
InternalWellGenerator::ConnectPerforationsToWellElements()
{
  // assign a well element to each perforation
  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    // get the perforation and its properties
    Perforation const * const perf =
      this->GetGroup< Perforation >( m_perforationList[iperf] );
    m_perfDistFromHead[iperf] = perf->GetDistanceFromWellHead();
    m_perfTransmissibility[iperf] = perf->GetWellTransmissibility();

    // search in all the elements of this well between head and bottom
    globalIndex iwelemTop = 0;
    globalIndex iwelemBottom = m_numElems - 1;

    // check the validity of the perforation before starting
    real64 const wellLength =
      m_nodeDistFromHead[m_elemToNodesMap[iwelemBottom][NodeLocation::BOTTOM]];

    GEOSX_ERROR_IF( m_perfDistFromHead[iperf] > wellLength,
                    "Distance from perforation "
                      << perf->getName()
                      << " to head is larger than well length for well "
                      << getName() );

    // start binary search
    const globalIndex maxNumSteps = m_numElems + 1;
    globalIndex currentNumSteps = 0;
    while( iwelemTop < iwelemBottom )
    {
      globalIndex iwelemMid = static_cast< globalIndex >(
        floor( static_cast< real64 >( iwelemTop + iwelemBottom ) / 2.0 ) );
      real64 const headToBottomDist =
        m_nodeDistFromHead[m_elemToNodesMap[iwelemMid][NodeLocation::BOTTOM]];

      if( headToBottomDist < m_perfDistFromHead[iperf] )
      {
        iwelemTop = iwelemMid + 1;
      }
      else
      {
        iwelemBottom = iwelemMid;
      }

      GEOSX_ERROR_IF( currentNumSteps > maxNumSteps,
                      "Perforation "
                        << perf->getName()
                        << " cannot be mapped to a well element in well "
                        << getName() );

      currentNumSteps++;
    }

    // set the index of the matched element
    globalIndex iwelemMatched = iwelemTop;
    GEOSX_ERROR_IF( iwelemMatched >= m_numElems,
                    "Invalid topology in well " << getName() );

    m_perfElemId[iperf] = iwelemMatched;

    // compute the physical location of the perforation
    globalIndex const inodeTop =
      m_elemToNodesMap[iwelemMatched][NodeLocation::TOP];
    globalIndex const inodeBottom =
      m_elemToNodesMap[iwelemMatched][NodeLocation::BOTTOM];
    real64 const elemLength =
      m_nodeDistFromHead[inodeBottom] - m_nodeDistFromHead[inodeTop];
    real64 const topToPerfDist =
      m_perfDistFromHead[iperf] - m_nodeDistFromHead[inodeTop];

    GEOSX_ERROR_IF( ( elemLength <= 0 ) || ( topToPerfDist < 0 ),
                    "Invalid topology in well " << getName() );

    R1Tensor v = m_nodeCoords[inodeBottom];
    v -= m_nodeCoords[inodeTop];
    m_perfCoords[iperf] = v;
    m_perfCoords[iperf] *= topToPerfDist / elemLength;
    m_perfCoords[iperf] += m_nodeCoords[inodeTop];
  }
}

globalIndex
InternalWellGenerator::GetNextSegmentIndex(
  globalIndex topSegId,
  globalIndex currentPolyNodeId ) const
{
  globalIndex nextSegId = -1;

  // get the index of the two segments sharing this node
  for( auto is : m_polyNodeToSegmentMap[currentPolyNodeId] )
  {
    // if this is not equal to the current segId, we found the next segId
    if( is != topSegId )
    {
      nextSegId = is;
      break;
    }
  }

  return nextSegId;
}

void
InternalWellGenerator::MergePerforations()
{
  array1d< array1d< localIndex > > elemToPerfMap;
  elemToPerfMap.resize( m_numElems );

  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    elemToPerfMap[m_perfElemId[iperf]].emplace_back( iperf );
  }

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    // collect the indices of the elems with more that one perforation
    if( elemToPerfMap[iwelem].size() > 1 )
    {
      // find the perforation with the largest Peaceman index and keep its location
      globalIndex iperfMaxTransmissibility = elemToPerfMap[iwelem][0];
      real64 maxTransmissibility =
        m_perfTransmissibility[iperfMaxTransmissibility];
      for( localIndex ip = 1; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( m_perfTransmissibility[elemToPerfMap[iwelem][ip]] > maxTransmissibility )
        {
          iperfMaxTransmissibility = elemToPerfMap[iwelem][ip];
          maxTransmissibility = m_perfTransmissibility[iperfMaxTransmissibility];
        }
      }

      // assign the coordinates of the perf with the largest trans to the other perfs on this elem
      for( localIndex ip = 0; ip < elemToPerfMap[iwelem].size(); ++ip )
      {
        if( elemToPerfMap[iwelem][ip] == iperfMaxTransmissibility )
        {
          continue;
        }

        GEOSX_LOG_RANK_0( "Moving perforation #"
                          << elemToPerfMap[iwelem][ip] << " of well " << getName()
                          << " from " << m_perfCoords[elemToPerfMap[iwelem][ip]]
                          << " to " << m_perfCoords[iperfMaxTransmissibility]
                          << " to make sure that no well element is shared "
                             "between two MPI ranks" );
        m_perfCoords[elemToPerfMap[iwelem][ip]] =
          m_perfCoords[iperfMaxTransmissibility];
      }
    }

    // in passing, check that there is always a perforation in the last well element (otherwise, the problem is not well
    // posed)
    for( localIndex iwelemPrev = 0; iwelemPrev < m_prevElemId[iwelem].size();
         ++iwelemPrev )
    {
      GEOSX_ERROR_IF( m_prevElemId[iwelem][iwelemPrev] == -1 &&
                        elemToPerfMap[iwelem].size() == 0,
                      "The bottom element of well "
                        << getName() << " does not have a perforation"
                        << " This is needed to have a well-posed problem" );
    }
  }
}

void
InternalWellGenerator::DebugWellGeometry() const
{
  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 0 )
  {
    return;
  }

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "InternalWellGenerator = " << getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX )
            << std::endl;
  std::cout << "Number of well elements = " << m_numElems << std::endl;

  for( globalIndex iwelem = 0; iwelem < m_numElems; ++iwelem )
  {
    std::cout << "m_elemCenterCoords[" << iwelem
              << "] = " << m_elemCenterCoords[iwelem] << std::endl;
    std::cout << "m_nextElemId[" << iwelem << "] = " << m_nextElemId[iwelem]
              << std::endl;
    std::cout << "m_prevElemId[" << iwelem << "] = " << m_prevElemId[iwelem][0]
              << std::endl;
    for( globalIndex inode = 0; inode < m_numNodesPerElem; ++inode )
    {
      std::cout << "m_elemToNodesMap[" << iwelem << "][" << inode
                << "] = " << m_elemToNodesMap[iwelem][inode] << std::endl;
    }
  }

  std::cout << "Number of well nodes = " << m_numNodes << std::endl;

  for( globalIndex inode = 0; inode < m_numNodes; ++inode )
  {
    std::cout << "m_nodeCoords[" << inode << "] = " << m_nodeCoords[inode]
              << std::endl;
    std::cout << "m_nodeDistFromHead[" << inode
              << "] = " << m_nodeDistFromHead[inode] << std::endl;
  }

  std::cout << "Number of perforations = " << m_numPerforations << std::endl;

  for( globalIndex iperf = 0; iperf < m_numPerforations; ++iperf )
  {
    std::cout << "m_perfCoords[" << iperf << "] = " << m_perfCoords[iperf]
              << std::endl;
    std::cout << "m_perfTransmissibility[" << iperf
              << "] = " << m_perfTransmissibility[iperf] << std::endl;
    std::cout << "m_perfElemId[" << iperf << "] = " << m_perfElemId[iperf]
              << std::endl;
  }
}

REGISTER_CATALOG_ENTRY( MeshGeneratorBase,
                        InternalWellGenerator,
                        std::string const &,
                        Group * const )
}  // namespace geosx
