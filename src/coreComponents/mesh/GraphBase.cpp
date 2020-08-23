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
 * @file GraphBase.cpp
 *
 */

#include "GraphBase.hpp"


namespace geosx
{
using namespace dataRepository;

GraphBase::GraphBase( string const & name, Group * const parent ):
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );
}

GraphBase::~GraphBase()
{}

GraphBase::CatalogInterface::CatalogType & GraphBase::GetCatalog()
{
  static GraphBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

}
