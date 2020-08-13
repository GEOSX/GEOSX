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
 * @file StringUtilities.cpp
 */

#include "codingUtilities/StringUtilities.hpp"
//#include "codingUtilities/UnitManager.h"

#include <stdarg.h>

namespace geosx
{
namespace stringutilities
{
/**
 * String tokenizing function
 **/
string_array
Tokenize( const std::string & str, const std::string & delimiters )
{
  string_array tokens;

  if( str.length() == 0 )
  {
    tokens.emplace_back( str );
  }
  else
  {
    bool usesNonWhitespaceDelimiters = false;
    std::string::size_type i = 0;
    while( delimiters[i] && !usesNonWhitespaceDelimiters )
    {
      usesNonWhitespaceDelimiters |= !isspace( int( delimiters[i] ) );
      ++i;
    }

    if( usesNonWhitespaceDelimiters )
    {
      // do not skip multiple adjacent delimiters - indicates empty strings
      size_t lastPos = 0;

      size_t newPos = lastPos;
      while( ( newPos = str.find_first_of( delimiters, lastPos ) ) != std::string::npos )
      {
        tokens.emplace_back( str.substr( lastPos, newPos - lastPos ) );
        lastPos = newPos + 1;
      }
      tokens.emplace_back( str.substr( lastPos, str.length() - lastPos ) );
    }
    else
    {
      // whitespace delimiters
      // skip multiple adjacent delimiters
      size_t lastPos = str.find_first_not_of( delimiters, 0 );
      lastPos = ( lastPos == std::string::npos ) ? 0 : lastPos;

      size_t newPos = lastPos;
      while( ( newPos = str.find_first_of( delimiters, lastPos ) ) != std::string::npos )
      {
        tokens.emplace_back( str.substr( lastPos, newPos - lastPos ) );
        lastPos = str.find_first_not_of( delimiters, newPos );
      }
      if( lastPos != std::string::npos )
        tokens.emplace_back( str.substr( lastPos, str.length() - lastPos ) );
    }
  }
  return tokens;
}

}  // namespace stringutilities
}  // namespace geosx
