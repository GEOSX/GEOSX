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
 * @file CellElementStencilMPFA.cpp
 */

#include "CellElementStencilMPFA.hpp"

namespace geosx
{
CellElementStencilMPFA::CellElementStencilMPFA() :
  StencilBase< CellElementStencilMPFA_Traits, CellElementStencilMPFA >()
{}

void
CellElementStencilMPFA::reserve( localIndex const size )
{
  m_elementRegionIndices.reserve( size * 9 );
  m_elementSubRegionIndices.reserve( size * 9 );
  m_elementIndices.reserve( size * 9 );
  m_weights.reserve( size * 9 );
}

void
CellElementStencilMPFA::add( localIndex const numPts,
                             localIndex const * const elementRegionIndices,
                             localIndex const * const elementSubRegionIndices,
                             localIndex const * const elementIndices,
                             real64 const * const weights,
                             localIndex const connectorIndex )
{
  GEOSX_ERROR_IF( numPts >= MAX_STENCIL_SIZE, "Maximum stencil size exceeded" );

  m_elementRegionIndices.appendArray( elementRegionIndices,
                                      elementRegionIndices + numPts );
  m_elementSubRegionIndices.appendArray( elementSubRegionIndices,
                                         elementSubRegionIndices + numPts );
  m_elementIndices.appendArray( elementIndices, elementIndices + numPts );
  m_weights.appendArray( weights, weights + numPts );

  m_connectorIndices[connectorIndex] = m_elementRegionIndices.size() - 1;
}

} /* namespace geosx */
