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
 * @file FaceElementToCellStencil.cpp
 */

#include "FaceElementToCellStencil.hpp"

namespace geosx
{

FaceElementToCellStencil::FaceElementToCellStencil():
  StencilBase< FaceElementToCellStencil_Traits, FaceElementToCellStencil >()
{
  m_faceNormal.resize( 0, 3 );
  m_cellToFaceVec.resize( 0, 3 );
}

void FaceElementToCellStencil::move( LvArray::MemorySpace const space )
{
  StencilBase< FaceElementToCellStencil_Traits, FaceElementToCellStencil >::move( space );
}

void FaceElementToCellStencil::add( localIndex const numPts,
                                    localIndex const * const elementRegionIndices,
                                    localIndex const * const elementSubRegionIndices,
                                    localIndex const * const elementIndices,
                                    real64 const * const weights,
                                    localIndex const connectorIndex )
{
  GEOSX_ERROR_IF_NE_MSG( numPts, 2, "Number of cells in faceToCell stencil should be 2" );

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

void FaceElementToCellStencil::addVectors( real64 const & transMultiplier,
                                           real64 const (&faceNormal)[3],
                                           real64 const (&cellToFaceVec)[3] )
{
  localIndex const oldSize = m_faceNormal.size( 0 );
  localIndex const newSize = oldSize + 1;
  m_faceNormal.resize( newSize );
  m_cellToFaceVec.resize( newSize );
  m_transMultiplier.resize( newSize );

  m_transMultiplier[oldSize] = transMultiplier;
  LvArray::tensorOps::copy< 3 >( m_faceNormal[oldSize], faceNormal );
  LvArray::tensorOps::copy< 3 >( m_cellToFaceVec[oldSize], cellToFaceVec );
}


} /* namespace geosx */
