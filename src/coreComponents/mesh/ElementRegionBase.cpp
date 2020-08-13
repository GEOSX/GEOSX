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

#include "ElementRegionBase.hpp"

#include "common/TimingMacros.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/solid/SolidBase.hpp"

namespace geosx
{
using namespace dataRepository;

ElementRegionBase::ElementRegionBase( string const & name, Group * const parent ) :
  ObjectManagerBase( name, parent ),
  m_numericalMethod()
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  this->RegisterGroup( viewKeyStruct::elementSubRegions );

  registerWrapper( viewKeyStruct::materialListString, &m_materialList )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "List of materials present in this region" );
}

ElementRegionBase::~ElementRegionBase()
{}

}  // namespace geosx
