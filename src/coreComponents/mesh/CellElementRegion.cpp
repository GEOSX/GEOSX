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

#include "CellElementRegion.hpp"
#include "mesh/generators/CellBlockABC.hpp"

namespace geosx
{
using namespace dataRepository;

CellElementRegion::CellElementRegion( string const & name, Group * const parent ):
  ElementRegionBase( name, parent )
{
  registerWrapper( viewKeyStruct::sourceCellBlockNamesString(), &m_cellBlockNames ).
    setInputFlag( InputFlags::OPTIONAL );

  registerWrapper( viewKeyStruct::coarseningRatioString(), &m_coarseningRatio ).
    setInputFlag( InputFlags::OPTIONAL );
}

CellElementRegion::~CellElementRegion()
{}

void CellElementRegion::generateMesh( Group & cellBlocks )
{
  Group & elementSubRegions = this->getGroup( viewKeyStruct::elementSubRegions() );

  for( string const & cellBlockName : this->m_cellBlockNames )
  {
    CellElementSubRegion & subRegion = elementSubRegions.registerGroup< CellElementSubRegion >( cellBlockName );
    CellBlockABC & source = cellBlocks.getGroup< CellBlockABC >( subRegion.getName() );
    subRegion.copyFromCellBlock( source );
  }
}

REGISTER_CATALOG_ENTRY( ObjectManagerBase, CellElementRegion, string const &, Group * const )

} /* namespace geosx */
