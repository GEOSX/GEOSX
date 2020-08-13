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
 * @file InterfaceTypes.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_

#include "common/GeosxConfig.hpp"

#ifdef GEOSX_USE_TRILINOS
  #include "linearAlgebra/interfaces/trilinos/TrilinosInterface.hpp"
#endif

#ifdef GEOSX_USE_HYPRE
  #include "linearAlgebra/interfaces/hypre/HypreInterface.hpp"
#endif

#ifdef GEOSX_USE_PETSC
  #include "linearAlgebra/interfaces/petsc/PetscInterface.hpp"
#endif

/// Macro to concatenate two strings (low level)
#define CONCAT_( A, B ) A##B

/// Macro to concatenate two strings (user level)
#define CONCAT( A, B ) CONCAT_( A, B )

namespace geosx
{
/// Alias for current interface
using LAInterface = CONCAT( GEOSX_LA_INTERFACE, Interface );

/// Alias for ParallelMatrix
using ParallelMatrix = LAInterface::ParallelMatrix;
/// Alias for ParallelVector
using ParallelVector = LAInterface::ParallelVector;
/// Alias for LinearSolver
using LinearSolver = LAInterface::LinearSolver;

/**
 * @brief High level interface to call the proper LAI initialize function.
 *
 * @param[in] argc standard argc as in any C main
 * @param[in] argv standard argv as in any C main
 */
inline void
setupLAI( int & argc, char **& argv )
{
#ifdef GEOSX_USE_TRILINOS
  TrilinosInterface::initialize( argc, argv );
#endif
#ifdef GEOSX_USE_HYPRE
  HypreInterface::initialize( argc, argv );
#endif
#ifdef GEOSX_USE_PETSC
  PetscInterface::initialize( argc, argv );
#endif
}

/**
 * @brief High level interface to call the proper LAI finalize function.
 */
inline void
finalizeLAI()
{
#ifdef GEOSX_USE_TRILINOS
  TrilinosInterface::finalize();
#endif
#ifdef GEOSX_USE_HYPRE
  HypreInterface::finalize();
#endif
#ifdef GEOSX_USE_PETSC
  PetscInterface::finalize();
#endif
}

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_INTERFACETYPES_HPP_*/
