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
 * @file NumericalMethodsManager.cpp
 */

#include "NumericalMethodsManager.hpp"

namespace geosx
{
using namespace dataRepository;

NumericalMethodsManager::NumericalMethodsManager( string const & name,
                                                  Group * const parent ) :
  Group( name, parent ),
  m_finiteElementDiscretizationManager(
    groupKeysStruct::finiteElementDiscretizations,
    this ),
  m_finiteVolumeManager( groupKeysStruct::finiteVolumeManager, this )
{
  setInputFlags( InputFlags::OPTIONAL );

  this->RegisterGroup( groupKeysStruct::finiteElementDiscretizations,
                       &m_finiteElementDiscretizationManager );
  this->RegisterGroup( groupKeysStruct::finiteVolumeManager,
                       &m_finiteVolumeManager );
}

NumericalMethodsManager::~NumericalMethodsManager()
{}

Group *
NumericalMethodsManager::CreateChild(
  string const & GEOSX_UNUSED_PARAM( childKey ),
  string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

} /* namespace geosx */
