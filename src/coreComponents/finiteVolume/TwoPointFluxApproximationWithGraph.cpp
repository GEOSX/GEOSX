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
 * @file TwoPointFluxApproximationWithGraph.cpp
 *
 */
#include "TwoPointFluxApproximationWithGraph.hpp"
#include "finiteVolume/FluxApproximationBase.hpp"
#include "finiteVolume/CellElementStencilTPFA.hpp"
#include "mesh/FaceElementRegion.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "mesh/GraphFromText.hpp"

namespace geosx
{

using namespace dataRepository;

TwoPointFluxApproximationWithGraph::TwoPointFluxApproximationWithGraph( std::string const & name,
                                                      Group * const parent )
  : TwoPointFluxApproximation( name, parent )
{
    registerWrapper( viewKeyStruct::graphString, &m_graphString )->
    setInputFlag( InputFlags::OPTIONAL )->
    setDescription( "Name of the graph to link to the graph." );


}

void TwoPointFluxApproximationWithGraph::recoverGraph() const
{
}


void TwoPointFluxApproximationWithGraph::computeCellStencil( MeshLevel & mesh) const
{
  const Group * tmp = this->GetGroupByPath( m_graphString );
  GEOSX_ERROR_IF( tmp == nullptr, "Can't find group with path " << m_graphString );
  const GraphFromText * graph = Group::group_cast< const GraphFromText * >( tmp );
  CellElementStencilTPFA & stencil = getStencil< CellElementStencilTPFA >( mesh, viewKeyStruct::cellStencilString );

  std::vector<Edge*> edges = graph->getEdges();
  for(long unsigned int i=0; i<edges.size(); i++)
  {
    stackArray1d< localIndex, 2 > regionIndex( 2 );
    stackArray1d< localIndex, 2 > subRegionIndex( 2 );
    stackArray1d< localIndex, 2 > elementIndex( 2 );
    stackArray1d< real64, 2 > stencilWeights( 2 );
    regionIndex[0]=LvArray::integerConversion< localIndex >(0);
    subRegionIndex[0]=LvArray::integerConversion< localIndex >(0);
    elementIndex[0]=LvArray::integerConversion< localIndex >(edges[i]->getN1()->getIndice());
    stencilWeights[0]=0.00000000000002;

    regionIndex[1]=LvArray::integerConversion< localIndex >(0);
    subRegionIndex[1]=LvArray::integerConversion< localIndex >(0);
    elementIndex[1]=LvArray::integerConversion< localIndex >(edges[i]->getN2()->getIndice());
    stencilWeights[1]=-0.00000000000002;
    localIndex faceIndex = LvArray::integerConversion< localIndex >(edges[i]->getIndice());
    std::cout<< regionIndex[0]<< " " << subRegionIndex[0] << " "
    << elementIndex[0] << " " << stencilWeights[0] << " " << regionIndex[1]<< " " << subRegionIndex[1] << " "
    << elementIndex[1] << " " << stencilWeights[1] << " " << faceIndex << "\n";
    stencil.add( 2,
                 regionIndex.data(),
                 subRegionIndex.data(),
                 elementIndex.data(),
                 stencilWeights.data(),
                 faceIndex);
  }

}



REGISTER_CATALOG_ENTRY( FluxApproximationBase, TwoPointFluxApproximationWithGraph, std::string const &, Group * const )

}
