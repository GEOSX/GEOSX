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


#include "GraphManager.hpp"

#include "mpiCommunications/SpatialPartition.hpp"
#include "GraphBase.hpp"
#include "common/TimingMacros.hpp"

namespace geosx
{

using namespace dataRepository;

GraphManager::GraphManager( std::string const & name,
                          Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::REQUIRED );
}

GraphManager::~GraphManager()
{}

Group * GraphManager::CreateChild( string const & childKey, string const & childName )
{
  GEOSX_LOG_RANK_0( "Adding Graph: " << childKey << ", " << childName );
  std::unique_ptr< GraphBase > solver = GraphBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup< GraphBase >( childName, std::move( solver ) );
}


void GraphManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from GraphBase here
  for( auto & catalogIter: GraphBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


void GraphManager::GenerateGraphs()
{
  forSubGroups< GraphBase >( [&]( GraphBase & meshGen )
  {
    meshGen.GenerateGraph();
  } );
}

} /* namespace geosx */
