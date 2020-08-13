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
 * @file TrilinosInterface.cpp
 */

#include "TrilinosInterface.hpp"
#include "linearAlgebra/interfaces/trilinos/TrilinosPreconditioner.hpp"

namespace geosx
{
void
TrilinosInterface::initialize( int & GEOSX_UNUSED_PARAM( argc ),
                               char **& GEOSX_UNUSED_PARAM( argv ) )
{}

void
TrilinosInterface::finalize()
{}

std::unique_ptr< PreconditionerBase< TrilinosInterface > >
TrilinosInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< TrilinosPreconditioner >( params );
}

}  // namespace geosx
