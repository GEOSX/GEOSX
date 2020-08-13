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
 * @file MeshUtilities.cpp
 *
 */

#include "MeshUtilities.hpp"
#include "managers/ObjectManagerBase.hpp"
#include "dataRepository/xmlWrapper.hpp"
#include "SimpleGeometricObjects/SimpleGeometricObjectBase.hpp"
#include "common/TimingMacros.hpp"
#include "mesh/NodeManager.hpp"

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

void
MeshUtilities::GenerateNodesets( dataRepository::Group const * geometries,
                                 NodeManager * const nodeManager )
{
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X =
    nodeManager->referencePosition();
  localIndex const numNodes = nodeManager->size();
  Group & sets = nodeManager->sets();

  for( int i = 0; i < geometries->GetSubGroups().size(); ++i )
  {
    SimpleGeometricObjectBase const * const object =
      geometries->GetGroup< SimpleGeometricObjectBase >( i );
    if( object != nullptr )
    {
      string name = object->getName();
      SortedArray< localIndex > & targetSet =
        sets.registerWrapper< SortedArray< localIndex > >( name )->reference();
      for( localIndex a = 0; a < numNodes; ++a )
      {
        if( object->IsCoordInObject( X[a] ) )
        {
          targetSet.insert( a );
        }
      }
    }
  }
}

}  // namespace geosx
