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

#include "StringUtilities.hpp"

#include "common/Logger.hpp"
#include <algorithm>
#include <fstream>

// TODO: this should not be here! Move what is not CPG-specific to codingUtilities/StringUtilities
// TODO: std::string -> string

namespace geosx
{

namespace cornerPointMeshStringUtilities
{

void eclipseDataBufferToVector( string & inputBuffer, std::vector< double > & v )
{
  string const star( "*" );
  localIndex positionOfStar( 0 );
  localIndex nRepetitions( 0 );
  double repeatedValue( 0. );

  // split the inputBuffer in chunks of information (spaces are used to split)
  std::istringstream splitBuffer( inputBuffer );
  string chunk;

  // Traverse through all chunks
  do
  {
    // Read a chunk from the buffer
    splitBuffer >> chunk;
    if( chunk.size() > 0 )
    {
      // does this word contain a star "*"
      positionOfStar = static_cast< localIndex >( chunk.find( star ) );
      if( positionOfStar > 0 )
      {
        string const nTimesValueIsRepeated( chunk.substr( 0, positionOfStar ) );
        string const valueAsString( chunk.substr( positionOfStar+1, chunk.size() ) );
        nRepetitions = std::stoi( nTimesValueIsRepeated );
        repeatedValue = std::stod( valueAsString );
        std::fill_n( back_inserter( v ), nRepetitions, repeatedValue ); // append nRepetitions of the repeatedValue as double
      }
      else
      {
        // The value is only present once (no repetition)
        repeatedValue = std::stod( chunk ); // conversion to double
        v.push_back( repeatedValue ); // append the repeatedValue once
      }
    }
    // While there is more to read
  }
  while (splitBuffer);
}

void eclipseDataBufferToVector( string & inputBuffer, std::vector< localIndex > & v )
{
  string const star( "*" );
  localIndex positionOfStar( 0 );
  localIndex nRepetitions( 0 );
  localIndex repeatedValue( 0 );

  // split the inputBuffer in chunks of information (spaces are used to split)
  std::istringstream splitBuffer( inputBuffer );

  // Traverse through all chunks
  do
  {
    // Read a word
    string word;
    splitBuffer >> word;
    if( word.size()>0 )
    {
      // does this word contain a "*"
      positionOfStar = static_cast< localIndex >( word.find( star ) );
      if( positionOfStar > 0 )
      {
        string const multiplierAsString( word.substr( 0, positionOfStar ) );
        string const valueAsString( word.substr( positionOfStar+1, word.size() ) );
        nRepetitions = std::stoi( multiplierAsString );
        repeatedValue = std::stoi( valueAsString );
        std::fill_n( back_inserter( v ), nRepetitions, repeatedValue ); // append nRepetitions of the repeatedValue as localIndex
      }
      else
      {
        // The value is only present once (no repetition)
        repeatedValue = static_cast< localIndex >( std::stoi( word ) ); // conversion to localIndex
        v.push_back( repeatedValue ); // append the repeatedValue once
      }
    }
    // While there is more to read
  }
  while (splitBuffer);
}

string fileToString( const string filePath )
{
  std::ifstream meshFile;
  string fileContent( "A" );

  //Open file
  meshFile.open( filePath );
  GEOSX_THROW_IF( !meshFile.is_open(),
                  "File could not be open",
                  InputError );

  //Transfer file content into string for easing broadcast
  std::stringstream buffer;
  buffer << meshFile.rdbuf();
  fileContent = buffer.str();

  //Close file
  meshFile.close();
  return fileContent;
}

void trim( string & str )
{
  if( str.size() > 0 )
  {
    localIndex const first = static_cast< localIndex >( str.find_first_not_of( ' ' ) );
    localIndex const last = static_cast< localIndex >(str.find_last_not_of( ' ' ));
    if( first < 0 )
    {
      str.resize( 0 );
    }
    else
    {
      str = str.substr( first, (last - first + 1));
    }

  }
}

bool removeStringAndFollowingContentFromLine( string toBeRemoved, string & line )
{
  // does the line have a "toBeRemoved" character(s)
  std::size_t pos = line.find( toBeRemoved );
  std::size_t initial_line_length = line.size();
  bool res = false;
  if( pos != string::npos )
  {
    res = true;

    std::size_t end_line_position = line.find( '\n' );
    if( end_line_position != string::npos )
    {
      // remove the character and everything afterwards
      line = line.substr( 0, pos )+line.substr( end_line_position+1, initial_line_length );
    }
    else
    {
      line = line.substr( 0, pos );
    }
  }
  return res;
}

void removeTab( string & v )
{
  std::replace( v.begin(), v.end(), '\t', ' ' );
}

void removeEndOfLine( string & v )
{
  std::replace( v.begin(), v.end(), '\r', ' ' );
  std::replace( v.begin(), v.end(), '\n', ' ' );
}

void removeExtraSpaces( string & v )
{
  v.erase( std::unique( v.begin(), v.end(),
                        []( char a, char b ) { return a == ' ' && b == ' '; } ), v.end());
}

} // end namespace CPMeshStringUtilities

} // end namespace geosx
