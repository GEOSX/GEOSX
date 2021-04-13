/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#ifndef GEOSX_MESHUTILITIES_CORNERPOINTMESH_STRINGUTILITIES_HPP_
#define GEOSX_MESHUTILITIES_CORNERPOINTMESH_STRINGUTILITIES_HPP_

#include "common/DataTypes.hpp"
#include <sstream>
#include <vector>

// TODO: this should not be here! Move what is not CPG-specific to codingUtilities/StringUtilities

namespace geosx
{

namespace cornerPointMeshStringUtilities
{

void eclipseDataBufferToVector( string & inputBuffer, std::vector< real64 > & outputVector );
void eclipseDataBufferToVector( string & inputBuffer, std::vector< localIndex > & outputVector );
string fileToString( const string filePath );
void trim( string & str );
bool removeStringAndFollowingContentFromLine( string toBeRemoved, string & line );
void removeTab( string & v );
void removeEndOfLine( string & v );
void removeExtraSpaces( string & v );

template< typename T >
void fromStringTo( string & data, std::vector< T > & v )
{
  v.reserve( data.size());
  std::istringstream iss( data );
  T sub;
  while( iss >> sub )
  {
    v.push_back( sub );
  }
}

template< typename T >
void fromStringTo( string & data, T & v )
{
  std::istringstream iss( data );
  iss >> v;
}

} // namespace cornerPointMeshStringUtilities

} // namespace geosx

#endif
