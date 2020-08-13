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
 * @file PetscUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP
#define GEOSX_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP

#include "common/DataTypes.hpp"
#include <petscsys.h>

namespace geosx
{
/**
 * @brief Converts a non-const array from GEOSX globalIndex type to PetscInt
 * @param[in] index the input array
 * @return the converted array
 */
inline PetscInt *
toPetscInt( globalIndex * const index )
{
  return reinterpret_cast< PetscInt * >( index );
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to PetscInt
 * @param[in] index the input array
 * @return the converted array
 */
inline PetscInt const *
toPetscInt( globalIndex const * const index )
{
  return reinterpret_cast< PetscInt const * >( index );
}

}  // namespace geosx

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_PETSCUTILS_HPP
