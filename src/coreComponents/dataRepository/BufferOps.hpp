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

#ifndef GEOSX_DATAREPOSITORY_BUFFEROPS_HPP_
#define GEOSX_DATAREPOSITORY_BUFFEROPS_HPP_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/traits.hpp"
#include "LvArray/src/limits.hpp"

#include <type_traits>

namespace geosx
{
/**
 * @brief Forward declaration of InterObjectRelation
 */
template< typename T >
class InterObjectRelation;

namespace bufferOps
{
/**
 * @brief Forward declaration of is_packable
 */
template< typename T >
struct is_packable_helper;

template< typename T >
constexpr bool is_noncontainer_type_packable = std::is_trivial< T >::value ||
  std::is_arithmetic< T >::value || traits::is_tensorT< T > || traits::is_string< T >;

template< typename T >
constexpr bool is_container = !is_noncontainer_type_packable< T >;

template< typename >
constexpr bool is_packable_array = false;

template< typename T, int NDIM, typename PERMUTATION >
constexpr bool is_packable_array< Array< T, NDIM, PERMUTATION > > =
  is_packable_helper< T >::value;

template< typename T, int NDIM, int USD >
constexpr bool is_packable_array< ArrayView< T, NDIM, USD > > =
  is_packable_helper< T >::value;

template< typename T, int NDIM, int USD >
constexpr bool is_packable_array< ArraySlice< T, NDIM, USD > > =
  is_packable_helper< T >::value;

template< typename T >
constexpr bool is_packable_array< ArrayOfArrays< T > > =
  is_packable_helper< T >::value;

template< typename >
constexpr bool is_packable_set = false;

template< typename T >
constexpr bool is_packable_set< SortedArray< T > > = is_packable_helper< T >::value;

template< typename >
constexpr bool is_packable_map = false;

template< typename T_KEY, typename T_VAL, typename SORTED >
constexpr bool is_packable_map< mapBase< T_KEY, T_VAL, SORTED > > =
  is_packable_helper< T_KEY >::value && is_packable_helper< T_VAL >::value;

template< typename T >
struct is_packable_helper
{
  static constexpr bool value = is_noncontainer_type_packable< T > ||
    is_packable_array< T > || is_packable_map< T > || is_packable_set< T >;
};

template< typename T >
constexpr bool is_packable = is_packable_helper< T >::value;

template< typename T >
constexpr bool is_packable_by_index = is_packable_array< T >;

template< typename >
constexpr bool is_map_packable_by_index = false;

template< typename T_KEY, typename T_VAL, typename SORTED >
constexpr bool is_map_packable_by_index< mapBase< T_KEY, T_VAL, SORTED > > =
  is_packable< T_KEY > && is_packable_by_index< T_VAL >;

template< typename T >
constexpr bool can_memcpy_helper =
  std::is_arithmetic< T >::value || traits::is_tensorT< T >;

template< typename T >
constexpr bool can_memcpy =
  can_memcpy_helper< std::remove_const_t< std::remove_pointer_t< T > > >;

//------------------------------------------------------------------------------
// Pack(buffer,var)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Pack(
  buffer_unit_type *& buffer,
  T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type *& buffer, const string & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type *& buffer, SortedArray< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Pack(
  buffer_unit_type *& buffer,
  T const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< is_packable< T >, localIndex >::type
Pack(
  buffer_unit_type *& buffer,
  ArrayView< T, NDIM, USD > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type *& buffer, ArrayOfArrays< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type *& buffer, ArrayOfSets< T > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Pack(
  buffer_unit_type *& buffer,
  MAP_TYPE const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type *& buffer,
      std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type *& buffer, InterObjectRelation< T > const & var );

//------------------------------------------------------------------------------
// fallthrough-implementation
template< bool DO_PACKING, typename T >
typename std::enable_if< !is_packable< T >, localIndex >::type
Pack(
  buffer_unit_type *& GEOSX_UNUSED_PARAM( buffer ),
  T const & GEOSX_UNUSED_PARAM( var ) )
{
  GEOSX_ERROR( "Trying to pack data type (" << typeid( T ).name()
                                            << ") but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
// PackArray(buffer,var,length)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackPointer(
  buffer_unit_type *& buffer,
  T const * const GEOSX_RESTRICT var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackPointer(
  buffer_unit_type *& buffer,
  T const * const GEOSX_RESTRICT var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
PackArray(
  buffer_unit_type *& buffer,
  arraySlice1d< T, USD > const & var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
PackArray(
  buffer_unit_type *& buffer,
  arraySlice1d< T, USD > const & var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
// PackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_indices >
typename std::enable_if< is_packable< T >, localIndex >::type
PackByIndex(
  buffer_unit_type *& buffer,
  ArrayView< T, NDIM, USD > const & var,
  const T_indices & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_indices >
localIndex
PackByIndex( buffer_unit_type *& buffer,
             ArrayOfArrays< T > const & var,
             T_indices const & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
PackByIndex( buffer_unit_type *& buffer,
             MAP_TYPE const & var,
             T_INDICES const & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !is_packable_by_index< T > && !is_map_packable_by_index< T >,
                         localIndex >::type
PackByIndex( buffer_unit_type *& GEOSX_UNUSED_PARAM( buffer ),
             T const & GEOSX_UNUSED_PARAM( var ),
             T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to pack data type ("
               << typeid( T ).name() << ") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
// Unpack(buffer,var)
//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
Unpack(
  buffer_unit_type const *& buffer,
  T & var );

//------------------------------------------------------------------------------
inline localIndex
Unpack( buffer_unit_type const *& buffer, string & var );

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< traits::is_tensorT< T >, localIndex >::type
Unpack(
  buffer_unit_type const *& buffer,
  T & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const *& buffer, SortedArray< T > & var );

//------------------------------------------------------------------------------
template< typename T, int NDIM, typename PERMUTATION >
typename std::enable_if< is_packable< T >, localIndex >::type
Unpack(
  buffer_unit_type const *& buffer,
  Array< T, NDIM, PERMUTATION > & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const *& buffer, ArrayOfArrays< T > & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const *& buffer, ArrayOfSets< T > & var );

//------------------------------------------------------------------------------
template< typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Unpack(
  buffer_unit_type const *& buffer,
  MAP_TYPE & map );

//------------------------------------------------------------------------------
template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( buffer_unit_type const *& buffer,
        std::pair< T_FIRST, T_SECOND > & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const *& buffer, InterObjectRelation< T > & var );

//------------------------------------------------------------------------------
template< typename T >
typename std::enable_if< !is_packable< T >, localIndex >::type
Unpack(
  buffer_unit_type const *& GEOSX_UNUSED_PARAM( buffer ),
  T & GEOSX_UNUSED_PARAM( var ) )
{
  GEOSX_ERROR( "Trying to unpack data type (" << typeid( T ).name()
                                              << ") but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
// UnpackArray(buffer,var,expectedLength)
//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackPointer(
  buffer_unit_type const *& buffer,
  T * const GEOSX_RESTRICT var,
  INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackPointer( buffer_unit_type const *& buffer,
               T * const GEOSX_RESTRICT var,
               INDEX_TYPE const expectedLength );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< std::is_trivial< T >::value, localIndex >::type
UnpackArray(
  buffer_unit_type const *& buffer,
  arraySlice1d< T, USD > const & var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE, int USD >
typename std::enable_if< !std::is_trivial< T >::value, localIndex >::type
UnpackArray(
  buffer_unit_type const *& buffer,
  arraySlice1d< T, USD > const & var,
  INDEX_TYPE const length );

//------------------------------------------------------------------------------
// UnpackByIndex(buffer,var,indices)
//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const *& buffer,
               ArrayView< T, NDIM, USD > & var,
               const T_indices & indices );

//------------------------------------------------------------------------------
template< typename T, typename T_indices >
localIndex
UnpackByIndex( buffer_unit_type const *& buffer,
               ArrayOfArrays< T > & var,
               T_indices const & indices );

//------------------------------------------------------------------------------
template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
UnpackByIndex( buffer_unit_type const *& buffer,
               MAP_TYPE & map,
               T_INDICES const & indices );

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
typename std::enable_if< !is_packable_by_index< T > && !is_map_packable_by_index< T >,
                         localIndex >::type
UnpackByIndex( buffer_unit_type const *& GEOSX_UNUSED_PARAM( buffer ),
               T & GEOSX_UNUSED_PARAM( var ),
               T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("
               << typeid( T ).name() << ") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type *& buffer,
      T const * const GEOSX_RESTRICT var,
      arraySlice1d< INDEX_TYPE const > const & indices,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const *& buffer,
        T * const GEOSX_RESTRICT var,
        arraySlice1d< INDEX_TYPE const > const & indices,
        INDEX_TYPE & length );

//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD >
localIndex
Pack( buffer_unit_type *& buffer,
      SortedArray< localIndex > const & var,
      SortedArray< globalIndex > const & unmappedGlobalIndices,
      arraySlice1d< globalIndex const, USD > const & localToGlobal );

//------------------------------------------------------------------------------
template< typename SORTED >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  SortedArray< localIndex > & var,
  SortedArray< globalIndex > & unmappedGlobalIndices,
  mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap,
  bool const clearExistingSet );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Pack( buffer_unit_type *& buffer, ArrayOfSets< T > const & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
Unpack( buffer_unit_type const *& buffer, ArrayOfSets< T > & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD >
localIndex
Pack( buffer_unit_type *& buffer,
      arraySlice1d< localIndex const, USD > const & var,
      globalIndex const * const unmappedGlobalIndices,
      localIndex const length,
      arraySlice1d< globalIndex const > const & localToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  localIndex_array & var,
  array1d< globalIndex > & unmappedGlobalIndices,
  mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
inline localIndex
UnpackSyncList(
  buffer_unit_type const *& buffer,
  localIndex_array & var,
  std::unordered_map< globalIndex, localIndex > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< typename SORTED, int USD >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arraySlice1d< localIndex, USD > & var,
  array1d< globalIndex > & unmappedGlobalIndices,
  localIndex const expectedLength,
  mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING >
localIndex
Pack( buffer_unit_type *& buffer,
      arrayView1d< localIndex const > const & var,
      arrayView1d< localIndex const > const & indices,
      arrayView1d< globalIndex const > const & localToGlobalMap,
      arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1 >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arrayView1d< localIndex > const & var,
  array1d< localIndex > const & indices,
  mapBase< globalIndex, localIndex, SORTED0 > const & globalToLocalMap,
  mapBase< globalIndex, localIndex, SORTED1 > const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED >
localIndex
Pack(
  buffer_unit_type *& buffer,
  arrayView1d< arrayView1d< localIndex const > const > const & var,
  mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
  arrayView1d< localIndex const > const & indices,
  arrayView1d< globalIndex const > const & localToGlobalMap,
  arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arrayView1d< localIndex_array > & var,
  array1d< localIndex > & indices,
  mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
  mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
  mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED >
localIndex
Pack(
  buffer_unit_type *& buffer,
  arrayView1d< SortedArray< localIndex > const > const & var,
  mapBase< localIndex, SortedArray< globalIndex >, SORTED > const & unmappedGlobalIndices,
  arrayView1d< localIndex const > const & indices,
  arrayView1d< globalIndex const > const & localToGlobalMap,
  arrayView1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2 >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arrayView1d< SortedArray< localIndex > > & var,
  localIndex_array & indices,
  mapBase< localIndex, SortedArray< globalIndex >, SORTED0 > & unmappedGlobalIndices,
  mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
  mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap,
  bool const clearFlag );

//------------------------------------------------------------------------------
template< bool DO_PACKING, int USD0, int USD1 >
localIndex
Pack( buffer_unit_type *& buffer,
      arrayView2d< localIndex const, USD0 > const & var,
      arrayView1d< localIndex > const & indices,
      arraySlice1d< globalIndex const, USD1 > const & localToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED, int USD >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arrayView2d< localIndex, USD > const & var,
  array1d< localIndex > & indices,
  mapBase< globalIndex, localIndex, SORTED > const & globalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename SORTED, int USD0 >
localIndex
Pack(
  buffer_unit_type *& buffer,
  arrayView2d< localIndex const, USD0 > const & var,
  mapBase< localIndex, array1d< globalIndex >, SORTED > const & unmappedGlobalIndices,
  arrayView1d< localIndex const > const & indices,
  arraySlice1d< globalIndex const > const & localToGlobalMap,
  arraySlice1d< globalIndex const > const & relatedObjectLocalToGlobalMap );

//------------------------------------------------------------------------------
template< typename SORTED0, typename SORTED1, typename SORTED2, int USD >
inline localIndex
Unpack(
  buffer_unit_type const *& buffer,
  arrayView2d< localIndex, USD > const & var,
  localIndex_array & indices,
  mapBase< localIndex, array1d< globalIndex >, SORTED0 > & unmappedGlobalIndices,
  mapBase< globalIndex, localIndex, SORTED1 > const & globalToLocalMap,
  mapBase< globalIndex, localIndex, SORTED2 > const & relatedObjectGlobalToLocalMap );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Pack(
  buffer_unit_type *& buffer,
  MAP_TYPE const & var );

//------------------------------------------------------------------------------
template< typename MAP_TYPE >
typename std::enable_if< is_packable_map< MAP_TYPE >, localIndex >::type
Unpack(
  buffer_unit_type const *& buffer,
  MAP_TYPE & map );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
Pack( buffer_unit_type *& buffer, MAP_TYPE const & var, T_INDICES const & packIndices );

//------------------------------------------------------------------------------
template< typename MAP_TYPE, typename T_INDICES >
typename std::enable_if< is_map_packable_by_index< MAP_TYPE >, localIndex >::type
Unpack( buffer_unit_type const *& buffer,
        MAP_TYPE & map,
        T_INDICES const & unpackIndices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T_FIRST, typename T_SECOND >
localIndex
Pack( buffer_unit_type *& buffer,
      std::pair< T_FIRST, T_SECOND > const & var );

//------------------------------------------------------------------------------
template< typename T_FIRST, typename T_SECOND >
localIndex
Unpack( buffer_unit_type const *& buffer,
        std::pair< T_FIRST, T_SECOND > & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
localIndex
Pack( buffer_unit_type *& buffer, InterObjectRelation< T > const & var );

template< typename... VARPACK >
localIndex
PackSize( VARPACK const &&... pack )
{
  buffer_unit_type * junk = nullptr;
  return Pack< false >( junk, pack... );
}

template< typename... VARPACK >
localIndex
PackSize( VARPACK &&... pack )
{
  buffer_unit_type * junk = nullptr;
  return Pack< false >( junk, pack... );
}

#ifdef GEOSX_USE_ARRAY_BOUNDS_CHECK
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
typename std::enable_if< !is_packable_by_index< T > && !is_map_packable_by_index< T >,
                         localIndex >::type
Pack( buffer_unit_type *& GEOSX_UNUSED_PARAM( buffer ),
      T const & GEOSX_UNUSED_PARAM( var ),
      T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to pack data type ("
               << typeid( T ).name() << ") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
typename std::enable_if< !is_packable_by_index< T > && !is_map_packable_by_index< T >,
                         localIndex >::type
Unpack( buffer_unit_type const *& GEOSX_UNUSED_PARAM( buffer ),
        T & GEOSX_UNUSED_PARAM( var ),
        T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("
               << typeid( T ).name() << ") but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename INDEX_TYPE >
localIndex
Pack( buffer_unit_type *& buffer,
      arraySlice1d< T > const & var,
      arraySlice1d< INDEX_TYPE > const & indices,
      INDEX_TYPE const length );

//------------------------------------------------------------------------------
template< typename T, typename INDEX_TYPE >
localIndex
Unpack( buffer_unit_type const *& buffer,
        arraySlice1d< T > & var,
        arraySlice1d< INDEX_TYPE > const & indices,
        INDEX_TYPE & length );

#endif /* GEOSX_USE_ARRAY_BOUNDS_CHECK */

} /* namespace bufferOps */
} /* namespace geosx */

#include "BufferOps_inline.hpp"

#endif /* GEOSX_DATAREPOSITORY_BUFFEROPS_HPP_ */
