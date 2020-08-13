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

#ifndef GEOSX_TESTFLOWKERNELHELPERS_HPP
#define GEOSX_TESTFLOWKERNELHELPERS_HPP

#include "common/DataTypes.hpp"

namespace geosx
{
namespace detail
{
template< typename T, int NDIM >
void
setArrayElement( ArrayView< T, NDIM > const & arr,
                 localIndex const dstIndex,
                 localIndex const srcIndex,
                 T const * const data )
{
  localIndex const elemSize =
    arr.size() / arr.size( 0 );  // assume singleParameterResizeIndex == 0
  T const * const dataPtr = data + srcIndex * elemSize;

  localIndex i = 0;
  LvArray::forValuesInSlice( arr[dstIndex], [&i, dataPtr]( T & value ) {
    value = dataPtr[i];
    ++i;
  } );
}

}  // namespace detail

/**
 * @brief Specializations of this struct provide a way to create data structures
 *        that flow kernels want (element/material accessors)
 * @tparam FULL whether to create "full" accessors (i.e. nested arrays of Array)
 *              or for a single region (i.e. just an Array)
 */
template< bool FULL >
struct AccessorHelper
{};

template<>
struct AccessorHelper< false >
{
  template< int NDIM, typename T >
  using ElementAccessor = Array< T, NDIM >;

  template< int NDIM, typename T >
  using MaterialAccessor = Array< T, NDIM >;

  template< int NDIM, typename T, typename... DIMS >
  static ElementAccessor< NDIM, T >
  makeElementAccessor(
    T const * const data,
    localIndex const stencilSize,
    arraySlice1d< localIndex const > const & GEOSX_UNUSED_PARAM( stencilRegIndices ),
    arraySlice1d< localIndex const > const & GEOSX_UNUSED_PARAM( stencilSubRegIndices ),
    arraySlice1d< localIndex const > const & stencilElemIndices,
    DIMS... otherDims )
  {
    localIndex numElems = 0;
    for( int i = 0; i < stencilSize; ++i )
    {
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    ElementAccessor< NDIM, T > acc( numElems, otherDims... );

    for( int i = 0; i < stencilSize; ++i )
    {
      detail::setArrayElement( acc, stencilElemIndices[i], i, data );
    }

    return acc;
  }

  template< int NDIM, typename T, typename... DIMS >
  static MaterialAccessor< NDIM, T >
  makeMaterialAccessor(
    T const * const data,
    localIndex const stencilSize,
    arraySlice1d< localIndex const > const & GEOSX_UNUSED_PARAM( stencilRegIndices ),
    arraySlice1d< localIndex const > const & GEOSX_UNUSED_PARAM( stencilSubRegIndices ),
    arraySlice1d< localIndex const > const & stencilElemIndices,
    localIndex GEOSX_UNUSED_PARAM( matIndex ),
    DIMS... otherDims )
  {
    localIndex numElems = 0;
    for( int i = 0; i < stencilSize; ++i )
    {
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    MaterialAccessor< NDIM, T > acc( numElems, 1, otherDims... );

    for( int i = 0; i < stencilSize; ++i )
    {
      detail::setArrayElement( acc, stencilElemIndices[i], i, data );
    }

    return acc;
  }
};

template<>
struct AccessorHelper< true >
{
  template< int NDIM, typename T >
  using ElementAccessor =
    ElementRegionManager::ElementViewAccessor< Array< T, NDIM > >;

  template< int NDIM, typename T >
  using MaterialAccessor =
    ElementRegionManager::MaterialViewAccessor< Array< T, NDIM > >;

  template< int NDIM, typename T, typename... DIMS >
  static ElementAccessor< NDIM, T >
  makeElementAccessor(
    T const * const data,
    localIndex const stencilSize,
    arraySlice1d< localIndex const > const & stencilRegIndices,
    arraySlice1d< localIndex const > const & stencilSubRegIndices,
    arraySlice1d< localIndex const > const & stencilElemIndices,
    DIMS... otherDims )
  {
    localIndex numRegions = 0, numSubRegions = 0, numElems = 0;
    for( int i = 0; i < stencilSize; ++i )
    {
      numRegions = std::max( numRegions, stencilRegIndices[i] + 1 );
      numSubRegions = std::max( numSubRegions, stencilSubRegIndices[i] + 1 );
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    ElementAccessor< NDIM, T > acc;
    acc.resize( numRegions );
    for( localIndex kr = 0; kr < numRegions; ++kr )
    {
      acc[kr].resize( numSubRegions );
      for( localIndex ksr = 0; ksr < numSubRegions; ++ksr )
      {
        acc[kr][ksr].resize( numElems, otherDims... );
      }
    }

    for( int i = 0; i < stencilSize; ++i )
    {
      detail::setArrayElement( acc[stencilRegIndices[i]][stencilSubRegIndices[i]],
                               stencilElemIndices[i],
                               i,
                               data );
    }

    return acc;
  }

  template< int NDIM, typename T, typename... DIMS >
  static ElementAccessor< NDIM, T >
  makeMaterialAccessor(
    T const * const data,
    localIndex const stencilSize,
    arraySlice1d< localIndex const > const & stencilRegIndices,
    arraySlice1d< localIndex const > const & stencilSubRegIndices,
    arraySlice1d< localIndex const > const & stencilElemIndices,
    localIndex matIndex,
    DIMS... otherDims )
  {
    localIndex numRegions = 0, numSubRegions = 0, numElems = 0;
    for( int i = 0; i < stencilSize; ++i )
    {
      numRegions = std::max( numRegions, stencilRegIndices[i] + 1 );
      numSubRegions = std::max( numSubRegions, stencilSubRegIndices[i] + 1 );
      numElems = std::max( numElems, stencilElemIndices[i] + 1 );
    }

    ElementAccessor< NDIM, T > acc;
    acc.resize( numRegions );
    for( localIndex kr = 0; kr < numRegions; ++kr )
    {
      acc[kr].resize( numSubRegions );
      for( localIndex ksr = 0; ksr < numSubRegions; ++ksr )
      {
        acc[kr][ksr].resize( matIndex + 1 );
        acc[kr][ksr][matIndex].resize( numElems, 1, otherDims... );
      }
    }

    for( int i = 0; i < stencilSize; ++i )
    {
      detail::setArrayElement(
        acc[stencilRegIndices[i]][stencilSubRegIndices[i]][matIndex],
        stencilElemIndices[i],
        i,
        data );
    }

    return acc;
  }
};

}  // namespace geosx

#endif  //GEOSX_TESTFLOWKERNELHELPERS_HPP
