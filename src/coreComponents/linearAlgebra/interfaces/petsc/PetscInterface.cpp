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
 * @file PetscInterface.cpp
 */

#include "PetscInterface.hpp"
#include "linearAlgebra/interfaces/petsc/PetscPreconditioner.hpp"

#include <petscsys.h>

namespace geosx
{
void
PetscInterface::initialize( int & GEOSX_UNUSED_PARAM( argc ),
                            char **& GEOSX_UNUSED_PARAM( argv ) )
{
  PetscOptionsSetValue( nullptr, "-no_signal_handler", "" );
  PetscOptionsSetValue( nullptr, "-on_error_abort", "" );
  PETSC_COMM_WORLD = MPI_COMM_GEOSX;
  PetscInitializeNoArguments();
}

void
PetscInterface::finalize()
{
  PetscFinalize();
}

std::unique_ptr< PreconditionerBase< PetscInterface > >
PetscInterface::createPreconditioner( LinearSolverParameters params )
{
  return std::make_unique< PetscPreconditioner >( params );
}

}  //namespace geosx
