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
 * @file DirichletBoundaryCondition.cpp
 */

#include "DirichletBoundaryCondition.hpp"

namespace geosx
{
using namespace dataRepository;

DirichletBoundaryCondition::DirichletBoundaryCondition( string const & name,
                                                        Group * const parent ) :
  FieldSpecificationBase( name, parent )
{
  // TODO Auto-generated constructor stub
}

DirichletBoundaryCondition::~DirichletBoundaryCondition()
{
  // TODO Auto-generated destructor stub
}

REGISTER_CATALOG_ENTRY( FieldSpecificationBase,
                        DirichletBoundaryCondition,
                        string const &,
                        Group * const )

} /* namespace geosx */
