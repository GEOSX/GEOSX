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

#ifndef GEOSX_MESH_BUFFEROPS_HPP_
#define GEOSX_MESH_BUFFEROPS_HPP_

#include "common/DataTypes.hpp"

class OrderedVariableToManyElementRelation;
class FixedToManyElementRelation;
class ElementRegionManager;

namespace geosx
{
namespace bufferOps
{
/**
 * @brief Pack (or size) relationship information for specific indices from
 *        an OrderedVariableToManyElementRelation into a buffer.
 * @tparam DO_PACKING Whether to pack or just return the size of the
 *                    data to be packed.
 * @param buffer A properly-sized buffer to contain the packed data.
 *               If DO_PACKING is false this may be a null pointer.
 * @param var A relationship mapping the indices in the packList to
 *            a list of element information (region, subregion, and
 *            global index). Where each index is related to a variable
 *            number of regions.
 * @param packList A list of indices for which to pack the related
 *                 element information.
 * @param elementRegionManager The element region manager for the
 *                             related regions.
 * @return The size (in bytes) of the data
 */
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type *& buffer,
      OrderedVariableToManyElementRelation const & var,
      arrayView1d< localIndex const > const & packList,
      ElementRegionManager const * const elementRegionManager );
/**
 * @brief Pack (or size) relationship information for specifc indices from
 *        an FixedToManyElementRelation into a buffer.
 * @tparam DO_PACKING Whether to pack or just return the size of
 *                    the data to be packed.
 * @param buffer A properly-sized buffer to contain the packed data.
 *               If DO_PACKING is false this may be a null pointer.
 * @param var A relationship mapping the indices in the packList to
 *            a list of element information (region, subregion, and
 *            global index). Where each index is related to the same number
 *            of regions.
 * @param packList A list of indices for which to pack the related
 *                 element information.
 * @param elementRegionManager The element region manager for the
 *                             related regions.
 * @return The size (in bytes) of the data
 */
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type *& buffer,
      FixedToManyElementRelation const & var,
      arrayView1d< localIndex const > const & packList,
      ElementRegionManager const * const elementRegionManager );
/**
 * @brief Unpack expected indices into an OrderedVariableToManyElementRelation.
 * @param buffer A buffer containing packed data.
 * @param var A relationship to unpack into, maps the indices in the
 *            packList to a list of element information (region, subregion,
 *            and index). Where each index is related to the same number
 *            of regions.
 * @param packList A list of indices for which to unpack the related
 *                 element information.
 * @param elementRegionManager
 * @param clearFlag If false, include prexisting relationship entries
 *                  from the var, else only use the unpacked entries.
 * @return The size (in bytes) of the data
 */
localIndex
Unpack( buffer_unit_type const *& buffer,
        OrderedVariableToManyElementRelation & var,
        arrayView1d< localIndex const > const & packList,
        ElementRegionManager const * const elementRegionManager,
        bool const clearFlag );
/**
 * @brief Unpack expected indices into a FixedToManyElementRelation.
 * @param buffer A buffer containing packed data.
 * @param var A relationship to unpack inot, maps the indices in the
 *            packList to a list of element information (region, subregion,
 *            and index). Where each index is related to the same number
 *            of regions.
 * @param packList A list of indices for which to pack the related
 *                 element information.
 * @param elementRegionManager
 * @param clearFlag If false, include prexisting relationship entries
 *                  from the var, else only use the unpacked entries.
 * @return The size (in bytes) of the data
 */
localIndex
Unpack( buffer_unit_type const *& buffer,
        FixedToManyElementRelation & var,
        arrayView1d< localIndex const > const & packList,
        ElementRegionManager const * const elementRegionManager,
        bool const clearFlag );
}  // namespace bufferOps
}  // namespace geosx
#endif
