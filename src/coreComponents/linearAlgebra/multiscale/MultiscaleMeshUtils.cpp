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
 * @file MultiscaleMeshUtils.cpp
 */

#include "MultiscaleMeshUtils.hpp"

#include "linearAlgebra/multiscale/MultiscaleMeshData.hpp"
#include "linearAlgebra/multiscale/MultiscaleMeshObjectManager.hpp"

namespace geosx
{
namespace multiscale
{
namespace meshUtils
{

void copySets( ObjectManagerBase const & srcManager,
               string const & mapKey,
               ObjectManagerBase & dstManager )
{
  arrayView1d< localIndex const > const map = srcManager.getReference< array1d< localIndex > >( mapKey );
  srcManager.sets().forWrappers< SortedArray< localIndex > >( [&]( dataRepository::Wrapper< SortedArray< localIndex > > const & setWrapper )
  {
    if( setWrapper.getName() != "all" ) // no use for "all" in MS mesh
    {
      SortedArrayView< localIndex const > const srcSet = setWrapper.referenceAsView();
      SortedArray< localIndex > & dstSet = dstManager.createSet( setWrapper.getName() );
      meshUtils::filterSet( srcSet, map, dstSet );
    }
  } );
}

} // namespace meshUtils
} // namespace multiscale
} // namespace geosx
