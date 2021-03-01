/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

#include "Path.hpp"
#include "common/Logger.hpp"

#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <vector>

namespace geosx
{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void getAbsolutePath( std::string const & path, std::string & absolutePath )
{
  char absFilePath[ PATH_MAX + 1 ];
  if( realpath( path.data(), absFilePath ) )
  {
    absolutePath = absFilePath;
  }
  else
  {
    char const * ret = getcwd( absFilePath, PATH_MAX + 1 );
    if( ret != nullptr )
    {
      GEOSX_THROW( "Could not get the absolute path for " << path << " from " << absFilePath, InputError );
    }

    GEOSX_THROW( "Could not get the absolute path for " << path, InputError );
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool isAbsolutePath( const std::string & path )
{
  return path[ 0 ] == '/';
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::istream & operator>>( std::istream & is, Path & p )
{
  is >> static_cast< std::string & >(p);
  if( !isAbsolutePath( p ) && !p.pathPrefix().empty())
  {
    getAbsolutePath( p.pathPrefix() + '/' + p, p );
  }
  return is;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void splitPath( std::string const & path, std::string & dirName, std::string & baseName )
{
  size_t pos = path.find_last_of( '/' );

  if( pos == std::string::npos )
  {
    dirName = ".";
    baseName = path;
  }
  else if( pos == 0 )
  {
    dirName = "/";
    baseName = path.substr( 1 );
    return;
  }
  else
  {
    dirName = path.substr( 0, pos );
    baseName = path.substr( pos + 1 );
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void readDirectory( std::string const & path, std::vector< std::string > & files )
{
  DIR * dirp = opendir( path.c_str() );
  struct dirent * dp;
  while( (dp = readdir( dirp )) != nullptr )
  {
    files.push_back( dp->d_name );
  }
  closedir( dirp );
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void makeDirsForPath( std::string const & path )
{
  constexpr mode_t mode = 0770;     // user and group rwx permissions

  std::string::size_type pos = 0;
  do
  {
    pos = path.find( '/', pos + 1 );
    std::string dir_name = path.substr( 0, pos );
    int const err = mkdir( dir_name.c_str(), mode );
    GEOSX_THROW_IF( err && ( errno != EEXIST ), "Failed to create a directories for " << path, std::runtime_error );
  }
  while (pos != std::string::npos);
}

} /* end namespace geosx */
