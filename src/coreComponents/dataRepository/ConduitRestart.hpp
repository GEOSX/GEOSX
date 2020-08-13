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
 * @file ConduitRestart.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_CONDUITRESTART_HPP_
#define GEOSX_DATAREPOSITORY_CONDUITRESTART_HPP_

// Source includes
#include "common/GeosxConfig.hpp"
#include "common/DataTypes.hpp"

// TPL includes
#include <conduit.hpp>

// System includes
#include <string>

/// @cond DO_NOT_DOCUMENT

#define CONDUIT_TYPE_INFO( T, CONDUIT_TYPE )                                 \
  template<>                                                                 \
  struct conduitTypeInfo< T >                                                \
  {                                                                          \
    using type = CONDUIT_TYPE;                                               \
    static constexpr int id = CONDUIT_TYPE##_ID;                             \
    static constexpr int sizeOfConduitType = sizeof( type );                 \
    static constexpr int numConduitValues = sizeof( T ) / sizeOfConduitType; \
    static_assert( sizeof( T ) % sizeOfConduitType == 0,                     \
                   #T " cannot be made made up of " #CONDUIT_TYPE "." );     \
  }

namespace geosx
{
namespace dataRepository
{
namespace internal
{
template< typename T >
struct conduitTypeInfo
{};

// Native integer types
CONDUIT_TYPE_INFO( char, CONDUIT_NATIVE_CHAR );
CONDUIT_TYPE_INFO( signed char, CONDUIT_NATIVE_SIGNED_CHAR );
CONDUIT_TYPE_INFO( unsigned char, CONDUIT_NATIVE_UNSIGNED_CHAR );

CONDUIT_TYPE_INFO( short, CONDUIT_NATIVE_SHORT );
CONDUIT_TYPE_INFO( int, CONDUIT_NATIVE_INT );
CONDUIT_TYPE_INFO( long, CONDUIT_NATIVE_LONG );
CONDUIT_TYPE_INFO( long long, CONDUIT_NATIVE_LONG_LONG );

CONDUIT_TYPE_INFO( unsigned short, CONDUIT_NATIVE_UNSIGNED_SHORT );
CONDUIT_TYPE_INFO( unsigned int, CONDUIT_NATIVE_UNSIGNED_INT );
CONDUIT_TYPE_INFO( unsigned long, CONDUIT_NATIVE_UNSIGNED_LONG );
CONDUIT_TYPE_INFO( unsigned long long, CONDUIT_NATIVE_UNSIGNED_LONG_LONG );

// Native floating point types
CONDUIT_TYPE_INFO( float, CONDUIT_NATIVE_FLOAT );
CONDUIT_TYPE_INFO( double, CONDUIT_NATIVE_DOUBLE );

// Tensor types
CONDUIT_TYPE_INFO( R1Tensor, CONDUIT_NATIVE_DOUBLE );

}  // namespace internal

template< typename T >
using conduitTypeInfo =
  internal::conduitTypeInfo< std::remove_const_t< std::remove_pointer_t< T > > >;

extern conduit::Node rootConduitNode;

std::string
writeRootFile( conduit::Node & root, std::string const & rootPath );

void
writeTree( std::string const & path );

void
loadTree( std::string const & path );

}  // namespace dataRepository
}  // namespace geosx

/// @endcond

#endif /* GEOSX_DATAREPOSITORY_CONDUITRESTART_HPP_ */
