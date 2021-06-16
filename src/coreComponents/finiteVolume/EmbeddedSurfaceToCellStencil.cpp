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
 * @file EmbeddedSurfaceToCellStencil.cpp
 */

#include "EmbeddedSurfaceToCellStencil.hpp"

namespace geosx
{

EmbeddedSurfaceToCellStencil::EmbeddedSurfaceToCellStencil():
  StencilBase< EmbeddedSurfaceToCellStencil_Traits, EmbeddedSurfaceToCellStencil >()
{}

void EmbeddedSurfaceToCellStencil::move( LvArray::MemorySpace const space )
{
  StencilBase< EmbeddedSurfaceToCellStencil_Traits, EmbeddedSurfaceToCellStencil >::move( space );
}

void EmbeddedSurfaceToCellStencil::add( localIndex const numPts,
                                        localIndex const * const elementRegionIndices,
                                        localIndex const * const elementSubRegionIndices,
                                        localIndex const * const elementIndices,
                                        real64 const * const weights,
                                        localIndex const connectorIndex )
{
  GEOSX_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in TPFA stencil should be 2" );

  localIndex const oldSize = m_elementRegionIndices.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_elementRegionIndices.resize( newSize, numPts );
  m_elementSubRegionIndices.resize( newSize, numPts );
  m_elementIndices.resize( newSize, numPts );
  m_weights.resize( newSize, numPts );

  for( localIndex a=0; a<numPts; ++a )
  {
    m_elementRegionIndices( oldSize, a ) = elementRegionIndices[a];
    m_elementSubRegionIndices( oldSize, a ) = elementSubRegionIndices[a];
    m_elementIndices( oldSize, a ) = elementIndices[a];
    m_weights( oldSize, a ) = weights[a];
  }
  m_connectorIndices[connectorIndex] = oldSize;
}


} /* namespace geosx */
