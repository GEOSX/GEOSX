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
 * @file MultiscaleMeshObjectManager.cpp
 */

#include "MultiscaleMeshObjectManager.hpp"

namespace geosx
{

using namespace dataRepository;

namespace multiscale
{

MultiscaleMeshObjectManager::MultiscaleMeshObjectManager( string const & name,
                                                          dataRepository::Group * const parent )
  : ObjectManagerBase( name, parent ),
  m_numOwnedObjects( 0 )
{
  registerWrapper( viewKeyStruct::dualObjectString(), &m_toDualRelation );
}

void MultiscaleMeshObjectManager::viewPackingExclusionList( SortedArray< localIndex > & exclusionList ) const
{
  exclusionList.insert( getWrapperIndex( viewKeyStruct::dualObjectString() ) );
}

void MultiscaleMeshObjectManager::setNumOwnedObjects( localIndex const n )
{
  GEOSX_ASSERT_GE( size(), n );
  m_numOwnedObjects = n;
}

} // namespace multiscale
} // namespace geosx
