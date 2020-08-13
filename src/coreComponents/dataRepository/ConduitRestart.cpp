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
 * @file ConduitRestart.cpp
 */

// Source includes
#include "ConduitRestart.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "common/TimingMacros.hpp"
#include "common/Path.hpp"

// TPL includes
#include <conduit_relay.hpp>

namespace geosx
{
namespace dataRepository
{
conduit::Node rootConduitNode;

std::string
writeRootFile( conduit::Node & root, std::string const & rootPath )
{
  std::string rootDirName, rootFileName;
  splitPath( rootPath, rootDirName, rootFileName );

  if( MpiWrapper::Comm_rank() == 0 )
  {
    makeDirsForPath( rootPath );

    root["protocol/name"] = "hdf5";
    root["protocol/version"] = CONDUIT_VERSION;

    root["number_of_files"] = MpiWrapper::Comm_size();
    root["file_pattern"] = rootFileName + "/rank_%07d.hdf5";

    root["number_of_trees"] = 1;
    root["tree_pattern"] = "/";

    conduit::relay::io::save( root, rootPath + ".root", "hdf5" );
  }

  MpiWrapper::Barrier( MPI_COMM_GEOSX );

  std::vector< char > buffer( rootPath.size() + 64 );
  GEOSX_ERROR_IF_GE( std::snprintf( buffer.data(),
                                    buffer.size(),
                                    "%s/rank_%07d.hdf5",
                                    rootPath.data(),
                                    MpiWrapper::Comm_rank() ),
                     1024 );
  return buffer.data();
}

std::string
readRootNode( std::string const & rootPath )
{
  std::string rankFilePattern;
  if( MpiWrapper::Comm_rank() == 0 )
  {
    conduit::Node node;
    conduit::relay::io::load( rootPath + ".root", "hdf5", node );

    int const nFiles = node.fetch_child( "number_of_files" ).value();
    GEOSX_ERROR_IF_NE( nFiles, MpiWrapper::Comm_size() );

    std::string const filePattern = node.fetch_child( "file_pattern" ).as_string();

    std::string rootDirName, rootFileName;
    splitPath( rootPath, rootDirName, rootFileName );

    rankFilePattern = rootDirName + "/" + filePattern;
    GEOSX_LOG_RANK_VAR( rankFilePattern );
  }

  MpiWrapper::Broadcast( rankFilePattern, 0 );

  char buffer[1024];
  GEOSX_ERROR_IF_GE(
    std::snprintf( buffer, 1024, rankFilePattern.data(), MpiWrapper::Comm_rank() ),
    1024 );
  return buffer;
}

/* Write out a restart file. */
void
writeTree( std::string const & path )
{
  GEOSX_MARK_FUNCTION;

  conduit::Node root;
  std::string const filePathForRank = writeRootFile( root, path );
  GEOSX_LOG_RANK( "Writing out restart file at " << filePathForRank );
  conduit::relay::io::save( rootConduitNode, filePathForRank, "hdf5" );
}

void
loadTree( std::string const & path )
{
  GEOSX_MARK_FUNCTION;
  std::string const filePathForRank = readRootNode( path );
  GEOSX_LOG_RANK( "Reading in restart file at " << filePathForRank );
  conduit::relay::io::load( filePathForRank, "hdf5", rootConduitNode );
}

} /* end namespace dataRepository */
} /* end namespace geosx */
