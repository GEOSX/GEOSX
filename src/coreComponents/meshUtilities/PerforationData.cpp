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
 * @file PerforationData.cpp
 */

#include "PerforationData.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "managers/DomainPartition.hpp"
#include "mesh/MeshForLoopInterface.hpp"
#include "meshUtilities/Perforation.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{
using namespace dataRepository;

PerforationData::PerforationData( string const & name, Group * const parent ) :
  ObjectManagerBase( name, parent ),
  m_numPerforationsGlobal( 0 )
{
  registerWrapper( viewKeyStruct::numPerforationsGlobalString,
                   &m_numPerforationsGlobal );
  registerWrapper( viewKeyStruct::reservoirElementRegionString,
                   &m_toMeshElements.m_toElementRegion );
  registerWrapper( viewKeyStruct::reservoirElementSubregionString,
                   &m_toMeshElements.m_toElementSubRegion );
  registerWrapper( viewKeyStruct::reservoirElementIndexString,
                   &m_toMeshElements.m_toElementIndex );
  registerWrapper( viewKeyStruct::wellElementIndexString, &m_wellElementIndex );
  registerWrapper( viewKeyStruct::locationString, &m_location );
  registerWrapper( viewKeyStruct::wellTransmissibilityString,
                   &m_wellTransmissibility );
}

PerforationData::~PerforationData()
{}

void
PerforationData::ComputeWellTransmissibility(
  MeshLevel const & mesh,
  WellElementSubRegion const * const wellElemSubRegion,
  string const & permeabilityKey )
{
  // get the permeability in the domain
  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > const perm =
    mesh.getElemManager()
      ->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor const > >(
        permeabilityKey );

  // for all the local perforations on this well
  for( localIndex iperf = 0; iperf < size(); ++iperf )
  {
    // if the well transmissibility has been read from the XML
    // then skip the computation of the well transmissibility carried out below
    if( m_wellTransmissibility[iperf] > 0 )
    {
      continue;
    }

    // get the indices of the reservoir element
    localIndex const er = m_toMeshElements.m_toElementRegion[iperf];
    localIndex const esr = m_toMeshElements.m_toElementSubRegion[iperf];
    localIndex const ei = m_toMeshElements.m_toElementIndex[iperf];

    real64 dx = 0;
    real64 dy = 0;
    real64 dz = 0;

    // get an approximate dx, dy, dz for the reservoir element
    // this is done by computing a bounding box
    GetReservoirElementDimensions( mesh, er, esr, ei, dx, dy, dz );

    real64 d1 = 0;
    real64 d2 = 0;
    real64 h = 0;
    real64 k1 = 0;
    real64 k2 = 0;

    // compute the vector perforation - well elem center
    // this vector will be used to decide whether this is a vectical well or not
    localIndex const wellElemIndex = m_wellElementIndex[iperf];
    R1Tensor vecWellElemCenterToPerf =
      wellElemSubRegion->getElementCenter()[wellElemIndex];
    vecWellElemCenterToPerf -= m_location[iperf];

    // check if this is a vertical well or a horizontal well
    // assign d1, d2, h, k1, and k2 accordingly
    DecideWellDirection( vecWellElemCenterToPerf,
                         dx,
                         dy,
                         dz,
                         perm[er][esr][ei],
                         d1,
                         d2,
                         h,
                         k1,
                         k2 );

    real64 const k21 = k1 > 0 ? k2 / k1 : 0;
    real64 const k12 = k2 > 0 ? k1 / k2 : 0;

    // compute the equivalent radius
    real64 const num = 0.28 * sqrt( d1 * d1 * sqrt( k21 ) + d2 * d2 * sqrt( k12 ) );
    real64 const den = std::pow( k12, 0.25 ) + std::pow( k21, 0.25 );
    real64 const rEq = ( den > 0 ) ? num / den : 0.0;

    real64 const kh = h * sqrt( k1 * k2 );

    arrayView1d< real64 const > const & wellElemRadius =
      wellElemSubRegion->getReference< array1d< real64 > >(
        WellElementSubRegion::viewKeyStruct::radiusString );

    GEOSX_ERROR_IF( rEq < wellElemRadius[wellElemIndex],
                    "The equivalent radius r_eq = "
                      << rEq << " is smaller than the well radius (r = "
                      << wellElemRadius[wellElemIndex] << ") in " << getName() );

    // compute the well Peaceman index
    m_wellTransmissibility[iperf] =
      2 * M_PI * kh / std::log( rEq / wellElemRadius[wellElemIndex] );

    GEOSX_ERROR_IF(
      m_wellTransmissibility[iperf] <= 0,
      "The well index is negative or equal to zero in " << getName() );
  }
}

void
PerforationData::GetReservoirElementDimensions( MeshLevel const & mesh,
                                                localIndex const er,
                                                localIndex const esr,
                                                localIndex const ei,
                                                real64 & dx,
                                                real64 & dy,
                                                real64 & dz ) const
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager = mesh.getNodeManager();
  CellElementRegion const * const region =
    Group::group_cast< CellElementRegion const * >( elemManager->GetRegion( er ) );
  CellBlock const * const subRegion =
    Group::group_cast< CellElementSubRegion const * >( region->GetSubRegion( esr ) );

  // compute the bounding box of the element
  real64 boxDims[3];
  computationalGeometry::GetBoundingBox( ei,
                                         subRegion->nodeList(),
                                         nodeManager->referencePosition(),
                                         boxDims );

  // dx and dz from bounding box
  dx = boxDims[0];
  dy = boxDims[1];

  // dz is computed as vol / (dx * dy)
  dz = subRegion->getElementVolume()[ei];
  dz /= dx * dy;

  GEOSX_ERROR_IF(
    dx <= 0 || dy <= 0 || dz <= 0,
    "The reservoir element dimensions (dx, dy, and dz) should be positive in "
      << getName() );
}

void
PerforationData::DecideWellDirection( R1Tensor const & vecWellElemCenterToPerf,
                                      real64 const & dx,
                                      real64 const & dy,
                                      real64 const & dz,
                                      R1Tensor const & perm,
                                      real64 & d1,
                                      real64 & d2,
                                      real64 & h,
                                      real64 & k1,
                                      real64 & k2 ) const
{
  // vertical well (approximately along the z-direction)
  if( fabs( vecWellElemCenterToPerf[2] ) > fabs( vecWellElemCenterToPerf[0] ) &&
      fabs( vecWellElemCenterToPerf[2] ) > fabs( vecWellElemCenterToPerf[1] ) )
  {
    d1 = dx;
    d2 = dy;
    h = dz;
    k1 = perm[0];
    k2 = perm[1];
  }
  // well approximately along the y-direction
  else if( fabs( vecWellElemCenterToPerf[1] ) > fabs( vecWellElemCenterToPerf[0] ) &&
           fabs( vecWellElemCenterToPerf[1] ) > fabs( vecWellElemCenterToPerf[2] ) )
  {
    d1 = dx;
    d2 = dz;
    h = dy;
    k1 = perm[0];
    k2 = perm[2];
  }
  // well approximately along the x-direction
  else
  {
    d1 = dy;
    d2 = dz;
    h = dx;
    k1 = perm[1];
    k2 = perm[2];
  }
}

void
PerforationData::ConnectToMeshElements( MeshLevel const & mesh,
                                        InternalWellGenerator const & wellGeometry )
{
  ElementRegionManager const * const elemManager = mesh.getElemManager();
  NodeManager const * const nodeManager = mesh.getNodeManager();
  ElementRegionManager::ElementViewAccessor< arrayView1d< R1Tensor const > > elemCenter =
    elemManager
      ->ConstructViewAccessor< array1d< R1Tensor >, arrayView1d< R1Tensor const > >(
        ElementSubRegionBase::viewKeyStruct::elementCenterString );

  arrayView1d< R1Tensor const > const & perfCoordsGlobal =
    wellGeometry.GetPerfCoords();
  arrayView1d< real64 const > const & perfTransGlobal =
    wellGeometry.GetPerfTransmissibility();

  resize( perfCoordsGlobal.size() );
  localIndex iperfLocal = 0;

  // loop over all the perforations
  for( globalIndex iperfGlobal = 0; iperfGlobal < perfCoordsGlobal.size();
       ++iperfGlobal )
  {
    R1Tensor const & coords = perfCoordsGlobal[iperfGlobal];

    // TODO actually trace coords
    // TODO what if a fracture element is located

    // find the closest reservoir element
    auto ret = minLocOverElemsInMesh(
      &mesh,
      [&]( localIndex const er, localIndex const esr, localIndex const ei ) -> real64 {
        R1Tensor v = coords;
        v -= elemCenter[er][esr][ei];
        return v.L2_Norm();
      } );

    // save the region, subregion and index
    localIndex const er = std::get< 0 >( ret.second );
    localIndex const esr = std::get< 1 >( ret.second );
    localIndex const ei = std::get< 2 >( ret.second );

    // a ghostRank check may be inserted if needed
    // (not needed in the current implementation)

    // check if perforation location is indeed inside the element
    CellBlock const * const cellBlock =
      elemManager->GetRegion( er )->GetSubRegion< CellElementSubRegion >( esr );
    array1d< array1d< localIndex > > faceNodes( cellBlock->numFacesPerElement() );

    for( localIndex kf = 0; kf < cellBlock->numFacesPerElement(); ++kf )
    {
      cellBlock->GetFaceNodes( ei, kf, faceNodes[kf] );
    }

    if( !computationalGeometry::IsPointInsidePolyhedron(
          nodeManager->referencePosition(),
          faceNodes,
          coords ) )
    {
      continue;
    }

    // TODO: what happens when the boundary is at the boundary of the MPI domain??

    // now construct the local data

    // store the indices of the mesh element
    m_toMeshElements.m_toElementRegion[iperfLocal] = er;
    m_toMeshElements.m_toElementSubRegion[iperfLocal] = esr;
    m_toMeshElements.m_toElementIndex[iperfLocal] = ei;

    // construct the local transmissibility and location maps
    m_wellTransmissibility[iperfLocal] = perfTransGlobal[iperfGlobal];
    m_location[iperfLocal] = coords;
    m_localToGlobalMap[iperfLocal++] = iperfGlobal;
  }

  // set the size based on the number of perforations matched with local reservoir elements
  resize( iperfLocal );
  ConstructGlobalToLocalMap();
}

void
PerforationData::ConnectToWellElements(
  InternalWellGenerator const & wellGeometry,
  unordered_map< globalIndex, localIndex > const & globalToLocalWellElemMap,
  globalIndex elemOffsetGlobal )
{
  arrayView1d< globalIndex const > const & perfElemIndexGlobal =
    wellGeometry.GetPerfElemIndex();

  for( localIndex iperfLocal = 0; iperfLocal < size(); ++iperfLocal )
  {
    globalIndex const iwelemGlobal =
      perfElemIndexGlobal[m_localToGlobalMap[iperfLocal]];
    globalIndex const ielemGlobal = elemOffsetGlobal + iwelemGlobal;
    GEOSX_ASSERT( globalToLocalWellElemMap.count( ielemGlobal ) > 0 );
    m_wellElementIndex[iperfLocal] = globalToLocalWellElemMap.at( ielemGlobal );
  }

  //DebugLocalPerforations();
}

void
PerforationData::DebugLocalPerforations() const
{
  if( size() == 0 )
  {
    return;
  }

  if( MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) != 1 )
  {
    return;
  }

  std::cout << std::endl;
  std::cout << "++++++++++++++++++++++++++" << std::endl;
  std::cout << "PerforationData = " << getName() << " of "
            << getParent()->getName() << std::endl;
  std::cout << "MPI rank = " << MpiWrapper::Comm_rank( MPI_COMM_GEOSX )
            << std::endl;
  std::cout << "Number of local perforations = " << size() << std::endl;

  for( localIndex iperf = 0; iperf < size(); ++iperf )
  {
    std::cout << "m_toMeshElements.m_toElementRegion[" << iperf
              << "] = " << m_toMeshElements.m_toElementRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementSubRegion[" << iperf
              << "] = " << m_toMeshElements.m_toElementSubRegion[iperf]
              << std::endl;
    std::cout << "m_toMeshElements.m_toElementIndex[" << iperf
              << "] = " << m_toMeshElements.m_toElementIndex[iperf] << std::endl;

    std::cout << "m_wellElementIndexLocal[" << iperf
              << "] = " << m_wellElementIndex[iperf] << std::endl;
    std::cout << "m_location[" << iperf << "] = " << m_location[iperf]
              << std::endl;
    std::cout << "m_wellTransmissibility[" << iperf
              << "] = " << m_wellTransmissibility[iperf] << std::endl;
  }
}

}  //namespace geosx
