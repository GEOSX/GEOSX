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
 * @file HypreUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_

#include "common/DataTypes.hpp"

#include <HYPRE_krylov.h>
#include <HYPRE_parcsr_ls.h>

namespace geosx
{
/**
 * @brief Convert GEOSX integer value to hypre int
 * @param index the input value
 * @return the converted value
 */
inline HYPRE_Int
toHYPRE_Int( integer const index )
{
  return LvArray::integerConversion< HYPRE_Int >( index );
}

/**
 * @brief Convert GEOSX global index value to hypre bigint
 * @param index the input value
 * @return the converted value
 */
inline HYPRE_BigInt
toHYPRE_BigInt( globalIndex const index )
{
  return LvArray::integerConversion< HYPRE_BigInt >( index );
}

/**
 * @brief Converts a non-const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt *
toHYPRE_BigInt( globalIndex * const index )
{
  return reinterpret_cast< HYPRE_BigInt * >( index );
}

/**
 * @brief Converts a const array from GEOSX globalIndex type to HYPRE_BigInt
 * @param[in] index the input array
 * @return the converted array
 */
inline HYPRE_BigInt const *
toHYPRE_BigInt( globalIndex const * const index )
{
  return reinterpret_cast< HYPRE_BigInt const * >( index );
}

/**
 * @brief Container for hypre preconditioner function pointers.
 */
struct HyprePrecFuncs
{
  /// Alias for setup function type
  using SetupFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );

  /// Alias for apply function type
  using ApplyFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );

  /// Alias for destroy function type
  using DestroyFunc = HYPRE_Int ( * )( HYPRE_Solver );

  SetupFunc setup {};          ///< pointer to setup function
  ApplyFunc apply {};          ///< pointer to apply function
  DestroyFunc destroy {};      ///< pointer to destroy function
  DestroyFunc aux_destroy {};  ///< pointer to auxillary destroy function
};

/**
 * @brief Container for hypre Krylov solver function pointers.
 */
struct HypreSolverFuncs
{
  /// Alias for set preconditioner function type
  using SetPrecondFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_PtrToParSolverFcn,
                                          HYPRE_Solver );

  /// Alias for setup function type
  using SetupFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );

  /// Alias for solve function type
  using SolveFunc = HYPRE_Int ( * )( HYPRE_Solver,
                                     HYPRE_ParCSRMatrix,
                                     HYPRE_ParVector,
                                     HYPRE_ParVector );

  /// Alias for get number of iterations function type
  using GetNumIter = HYPRE_Int ( * )( HYPRE_Solver solver,
                                      HYPRE_Int * num_iterations );

  /// Alias for get final residual norm function type
  using GetFinalNorm = HYPRE_Int ( * )( HYPRE_Solver solver, HYPRE_Real * norm );

  /// Alias for destroy function type
  using DestroyFunc = HYPRE_Int ( * )( HYPRE_Solver );

  SetPrecondFunc setPrecond {};  ///< pointer to set preconditioner function
  SetupFunc setup {};            ///< pointer to setup function
  SolveFunc solve {};            ///< pointer to solve function
  GetNumIter getNumIter {};      ///< pointer to get number of iterations function
  GetFinalNorm getFinalNorm {};  ///< pointer to get final residual norm function
  DestroyFunc destroy {};        ///< pointer to destroy function
};

/**
 * @brief Container for hypre preconditioner auxiliary data.
 */
struct HyprePrecAuxData
{
  array1d< HYPRE_Int > point_marker_array;  ///< array1d of unique tags for local degrees of freedom
};

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_INTERFACES_HYPREUTILS_HPP_*/
