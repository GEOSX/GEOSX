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
 * @file MultiscaleMeshData.hpp
 */
#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHDATA_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHDATA_HPP

#include "mesh/ExtrinsicMeshData.hpp"

namespace geosx
{
namespace multiscale
{
namespace meshData
{

EXTRINSIC_MESH_DATA_TRAIT( OrigElementRegion,
                           "origElementRegion",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Region index of the original geosx element" );

EXTRINSIC_MESH_DATA_TRAIT( OrigElementSubRegion,
                           "origElementSubRegion",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Subregion index of the original geosx element" );

EXTRINSIC_MESH_DATA_TRAIT( OrigElementIndex,
                           "origElementIndex",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Index of the original geosx element" );

EXTRINSIC_MESH_DATA_TRAIT( OrigNodeIndex,
                           "origNodeIndex",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Index of the original geosx node" );

EXTRINSIC_MESH_DATA_TRAIT( CoarseCellLocalIndex,
                           "coarseCellLocalIndex",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Local index of the coarse scale cell" );

EXTRINSIC_MESH_DATA_TRAIT( CoarseCellGlobalIndex,
                           "coarseCellGlobalIndex",
                           array1d< globalIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Global index of the coarse scale cell" );

EXTRINSIC_MESH_DATA_TRAIT( FineCellLocalIndices,
                           "fineCellLocalIndices",
                           ArrayOfArrays< localIndex >,
                           -1,
                           NOPLOT,
                           NO_WRITE,
                           "Local indices of the fine scale cells that make up a coarse cell" );

EXTRINSIC_MESH_DATA_TRAIT( CoarseNodeLocalIndex,
                           "coarseNodeLocalIndex",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Local index of the coarse scale node" );

EXTRINSIC_MESH_DATA_TRAIT( CoarseNodeGlobalIndex,
                           "coarseNodeGlobalIndex",
                           array1d< globalIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Global index of the coarse scale node" );

EXTRINSIC_MESH_DATA_TRAIT( FineNodeLocalIndex,
                           "fineNodeLocalIndex",
                           array1d< localIndex >,
                           -1,
                           LEVEL_0,
                           NO_WRITE,
                           "Local index of the fine scale node corresponding to a coarse node" );

EXTRINSIC_MESH_DATA_TRAIT( NodeToCoarseSubdomain,
                           "nodeToCoarseSubdomain",
                           ArrayOfSets< globalIndex >,
                           -1,
                           NOPLOT,
                           NO_WRITE,
                           "Map of node to coarse subdomain global index (including boundaries)" );

} // namespace meshData
} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALEMESHDATA_HPP
