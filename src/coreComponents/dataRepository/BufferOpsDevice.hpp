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

#ifndef GEOSX_DATAREPOSITORY_BUFFEROPSDEVICE_H_
#define GEOSX_DATAREPOSITORY_BUFFEROPSDEVICE_H_

#include "common/DataTypes.hpp"
#include "codingUtilities/Utilities.hpp"
#include "codingUtilities/static_if.hpp"
#include "codingUtilities/traits.hpp"
#include "LvArray/src/limits.hpp"
#include "BufferOps.hpp"

#include "rajaInterface/GEOS_RAJA_Interface.hpp"

#include <type_traits>

namespace geosx
{
namespace bufferOps
{
//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackDevice(
  buffer_unit_type *& buffer,
  ArrayView< T const, NDIM, USD > const & var );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T >
GEOSX_HOST_DEVICE localIndex
PackDevice( buffer_unit_type *& GEOSX_UNUSED_PARAM( buffer ),
            T const & GEOSX_UNUSED_PARAM( var ) )
{
  GEOSX_ERROR( "Trying to pack data type ("
               << LvArray::system::demangleType< T >()
               << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
PackByIndexDevice(
  buffer_unit_type *& buffer,
  ArrayView< T const, NDIM, USD > const & var,
  T_INDICES const & indices );

//------------------------------------------------------------------------------
template< bool DO_PACKING, typename T, typename T_INDICES >
localIndex
PackByIndexDevice( buffer_unit_type *& GEOSX_UNUSED_PARAM( buffer ),
                   T const & GEOSX_UNUSED_PARAM( var ),
                   T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to pack data type ("
               << LvArray::system::demangleType< T >()
               << ") on device but type is not packable by index." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackDevice(
  buffer_unit_type const *& buffer,
  ArrayView< T, NDIM, USD > const & var );

//------------------------------------------------------------------------------
template< typename T >
localIndex
UnpackDevice( buffer_unit_type const *& GEOSX_UNUSED_PARAM( buffer ),
              T & GEOSX_UNUSED_PARAM( var ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("
               << LvArray::system::demangleType< T >()
               << ") on device but type is not packable." );
  return 0;
}

//------------------------------------------------------------------------------
template< typename T, int NDIM, int USD, typename T_INDICES >
typename std::enable_if< can_memcpy< T >, localIndex >::type
UnpackByIndexDevice(
  buffer_unit_type const *& buffer,
  ArrayView< T, NDIM, USD > const & var,
  T_INDICES const & indices );

//------------------------------------------------------------------------------
template< typename T, typename T_INDICES >
localIndex
UnpackByIndexDevice( buffer_unit_type const *& GEOSX_UNUSED_PARAM( buffer ),
                     T & GEOSX_UNUSED_PARAM( var ),
                     T_INDICES const & GEOSX_UNUSED_PARAM( indices ) )
{
  GEOSX_ERROR( "Trying to unpack data type ("
               << LvArray::system::demangleType< T >()
               << ") but type is not packable by index." );
  return 0;
}

}  // namespace bufferOps
}  // namespace geosx

#endif  // GEOSX_DATAREPOSITORY_BUFFEROPSDEVICE_H_
