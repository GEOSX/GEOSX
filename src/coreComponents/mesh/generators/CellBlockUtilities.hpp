/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2020-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_CELLBLOCKUTILITIES_HPP
#define GEOSX_CELLBLOCKUTILITIES_HPP

#include "common/DataTypes.hpp"

namespace geosx
{

/**
 * @brief Free function that generates face to edges, edge to faces and edge to nodes mappings.
 * @param[in] numNodes The number of nodes.
 * @param[in] faceToNodeMap Face to node mappings as an input.
 * @param[out] faceToEdgeMap Face to edges will be resized and filled.
 * @param[out] edgeToFaceMap Ege to faces will be resized and filled.
 * @param[out] edgeToNodeMap Edge to nodes will be resized and filled.
 * @return The number of edges.
 */
localIndex buildEdgeMaps( localIndex numNodes,
                          ArrayOfArraysView< localIndex const > const & faceToNodeMap,
                          ArrayOfArrays< localIndex > & faceToEdgeMap,
                          ArrayOfSets< localIndex > & edgeToFaceMap,
                          array2d< localIndex > & edgeToNodeMap );

/**
 * @brief Get the local indices of the nodes in a face of the element.
 * @param[in] elementType Type of the element (like C3D8, etc...)
 * @param[in] iElement the local index of the target element
 * @param[in] iFace the local index of the target face in the element  (this will be 0-numFacesInElement)
 * @param[in] elementToNodes Element to nodes mapping.
 * @param[out] nodeIndices A reference to the array of node indices of the face. Gets resized at the proper size.
 */
void getFaceNodes( string const & elementType,
                   localIndex const iElement,
                   localIndex const iFace,
                   array2d< localIndex, cells::NODE_MAP_PERMUTATION > const & elementToNodes,
                   array1d< localIndex > & nodeIndices );

}

#endif // include guard
