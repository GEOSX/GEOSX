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

/**
 * @file MeshUtilities.cpp
 *
 */

#include "MeshUtilities.hpp"
#include "mesh/simpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "mesh/NodeManager.hpp"
#include "LvArray/src/tensorOps.hpp"

namespace geosx
{
using namespace dataRepository;

MeshUtilities::MeshUtilities()
{
  // TODO Auto-generated constructor stub

}

MeshUtilities::~MeshUtilities()
{
  // TODO Auto-generated destructor stub
}



void MeshUtilities::generateNodesets( dataRepository::Group const & geometries,
                                      NodeManager & nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const X = nodeManager.referencePosition();
  localIndex const numNodes = nodeManager.size();
  Group & sets = nodeManager.sets();

  geometries.forSubGroups< SimpleGeometricObjectBase >( [&] ( SimpleGeometricObjectBase const & object )
  {
    string const & name = object.getName();
    SortedArray< localIndex > & targetSet = sets.registerWrapper< SortedArray< localIndex > >( name ).reference();
    for( localIndex a=0; a<numNodes; ++a )
    {
      real64 nodeCoord[3] = LVARRAY_TENSOROPS_INIT_LOCAL_3( X[a] );
      if( object.isCoordInObject( nodeCoord ))
      {
        targetSet.insert( a );
      }
    }
  } );
}

} /// namespace geosx
