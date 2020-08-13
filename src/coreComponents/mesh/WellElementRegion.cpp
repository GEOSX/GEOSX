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
 * @file WellElementRegion.cpp
 */

#include "WellElementRegion.hpp"

#include "mpiCommunications/MpiWrapper.hpp"
#include "mesh/WellElementSubRegion.hpp"

namespace geosx
{
using namespace dataRepository;

WellElementRegion::WellElementRegion( string const & name, Group * const parent ) :
  ElementRegionBase( name, parent ),
  m_subRegionName( name + "uniqueSubRegion" ),
  m_wellControlsName( "" ),
  m_wellGeneratorName( "" )
{
  registerWrapper( viewKeyStruct::wellControlsString, &m_wellControlsName );
  registerWrapper( viewKeyStruct::wellGeneratorString, &m_wellGeneratorName );

  this->GetGroup( viewKeyStruct::elementSubRegions )
    ->RegisterGroup< WellElementSubRegion >( m_subRegionName );
}

WellElementRegion::~WellElementRegion()
{}

void
WellElementRegion::GenerateWell( MeshLevel & mesh,
                                 InternalWellGenerator const & wellGeometry,
                                 globalIndex nodeOffsetGlobal,
                                 globalIndex elemOffsetGlobal )
{
  // get the (unique) subregion
  WellElementSubRegion * const subRegion =
    this->GetGroup( ElementRegionBase::viewKeyStruct::elementSubRegions )
      ->GetGroup< WellElementSubRegion >( m_subRegionName );

  GEOSX_ERROR_IF( subRegion == nullptr,
                  "Well subRegion " << this->m_subRegionName
                                    << " not found in well region " << getName() );
  subRegion->SetWellControlsName( m_wellControlsName );

  PerforationData * const perforationData = subRegion->GetPerforationData();
  perforationData->SetNumPerforationsGlobal( wellGeometry.GetNumPerforations() );

  globalIndex const numElemsGlobal = wellGeometry.GetNumElements();
  globalIndex const numPerforationsGlobal = wellGeometry.GetNumPerforations();

  // 1) select the local perforations based on connectivity to the local reservoir elements
  subRegion->ConnectPerforationsToMeshElements( mesh, wellGeometry );

  globalIndex const matchedPerforations =
    MpiWrapper::Sum( perforationData->size() );
  GEOSX_ERROR_IF(
    matchedPerforations != numPerforationsGlobal,
    "Invalid mapping perforation-to-element in well " << this->getName() );

  // 2) classify well elements based on connectivity to local mesh partition
  array1d< integer > elemStatusGlobal;
  elemStatusGlobal.resizeDefault( numElemsGlobal,
                                  WellElementSubRegion::WellElemStatus::UNOWNED );

  arrayView1d< globalIndex const > const & perfElemIdGlobal =
    wellGeometry.GetPerfElemIndex();

  for( localIndex iperfGlobal = 0; iperfGlobal < numPerforationsGlobal;
       ++iperfGlobal )
  {
    globalIndex const iwelemGlobal = perfElemIdGlobal[iperfGlobal];

    if( perforationData->globalToLocalMap().count( iperfGlobal ) > 0 )
    {
      elemStatusGlobal[iwelemGlobal] |=
        WellElementSubRegion::WellElemStatus::LOCAL;
    }
    else
    {
      elemStatusGlobal[iwelemGlobal] |=
        WellElementSubRegion::WellElemStatus::REMOTE;
    }
  }

  // 3) select the local well elements and mark boundary nodes (for ghosting)
  subRegion->Generate( mesh,
                       wellGeometry,
                       elemStatusGlobal,
                       nodeOffsetGlobal,
                       elemOffsetGlobal );

  // 4) find out which rank is the owner of the top segment
  localIndex const refElemIdLocal = subRegion->GetTopWellElementIndex();

  array1d< localIndex > allRankTopElem;
  MpiWrapper::allGather( refElemIdLocal, allRankTopElem );
  int topRank = -1;
  for( int irank = 0; irank < allRankTopElem.size(); ++irank )
  {
    if( allRankTopElem[irank] >= 0 )
    {
      GEOSX_ASSERT( topRank < 0 );
      topRank = irank;
    }
  }
  GEOSX_ASSERT( topRank >= 0 );
  subRegion->SetTopRank( topRank );

  // 5) construct the local perforation to well element map
  perforationData->ConnectToWellElements( wellGeometry,
                                          subRegion->globalToLocalMap(),
                                          elemOffsetGlobal );
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase,
                        WellElementRegion,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
