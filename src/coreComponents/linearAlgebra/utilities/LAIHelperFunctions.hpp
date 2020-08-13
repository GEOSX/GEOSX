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

/*
 * LAIHelperFunctions.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_

#include "common/DataTypes.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mesh/NodeManager.hpp"
#include "mesh/ElementRegionManager.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{
/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  nodeManager       the node manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void
CreatePermutationMatrix( NodeManager const * const nodeManager,
                         localIndex const nRows,
                         localIndex const nCols,
                         int const nDofPerNode,
                         string const dofKey,
                         ParallelMatrix & permutationMatrix );

/**
 * @brief Create a permutation matrix for a given nodal variable.
 * @param[in]  elemManager       the element region manager
 * @param[in]  nRows             number of local rows in the matrix
 * @param[in]  nCols             number of local columns in the matrix
 * @param[in]  nDofPerNode       number of degrees-of-freedom per node
 * @param[in]  dofKey            DofManager key used to access dof index array
 * @param[out] permutationMatrix the target matrix
 */
void
CreatePermutationMatrix( ElementRegionManager const * const elemManager,
                         localIndex const nRows,
                         localIndex const nCols,
                         int const nDofPerNode,
                         string const DofKey,
                         ParallelMatrix & permutationMatrix );

/**
 * @brief Permute a vector.
 * @param[in] vector            the source vector
 * @param[in] permutationMatrix the permutation matrix
 * @return the permuted vector
 */
ParallelVector
PermuteVector( ParallelVector const & vector,
               ParallelMatrix const & permutationMatrix );

/**
 * @brief Permute rows and columns of a square matrix.
 * @param[in] matrix            the source matrix
 * @param[in] permutationMatrix permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix
PermuteMatrix( ParallelMatrix const & matrix,
               ParallelMatrix const & permutationMatrix );

/**
 * Permute rows and columns of a rectangular matrix
 * @param[in] matrix                 the source matrix
 * @param[in] permutationMatrixLeft  left permutation matrix
 * @param[in] permutationMatrixRight right permutation matrix
 * @return the permuted matrix
 */
ParallelMatrix
PermuteMatrix( ParallelMatrix const & matrix,
               ParallelMatrix const & permutationMatrixLeft,
               ParallelMatrix const & permutationMatrixRight );

void
PrintPermutedVector( ParallelVector const & vector,
                     ParallelMatrix const & permuationMatrix,
                     std::ostream & os );

void
PrintPermutedMatrix( ParallelMatrix const & matrix,
                     ParallelMatrix const & permutationMatrix,
                     std::ostream & os );

void
PrintPermutedMatrix( ParallelMatrix const & matrix,
                     ParallelMatrix const & permutationMatrixLeft,
                     ParallelMatrix const & permutationMatrixRight,
                     std::ostream & os );

/**
 * @brief Apply a separate component approximation (filter) to a matrix.
 * @tparam MATRIX the type of matrices
 * @param src         the source matrix
 * @param dst         the target (filtered) matrix
 * @param dofsPerNode number of degrees-of-freedom per node
 */
template< typename MATRIX >
void
SeparateComponentFilter( MATRIX const & src,
                         MATRIX & dst,
                         const localIndex dofsPerNode )
{
  GEOSX_ERROR_IF( dofsPerNode < 2, "Function requires dofsPerNode > 1" );

  const localIndex localRows = src.numLocalRows();
  const localIndex maxEntries = src.maxRowLength();
  const localIndex maxDstEntries = maxEntries / dofsPerNode;

  dst.createWithLocalSize( localRows, maxEntries, MPI_COMM_WORLD );
  dst.open();

  array1d< real64 > srcValues;
  array1d< real64 > dstValues( maxDstEntries );

  array1d< globalIndex > srcIndices;
  array1d< globalIndex > dstIndices( maxDstEntries );

  for( globalIndex row = src.ilower(); row < src.iupper(); ++row )
  {
    const globalIndex rowComponent = row % dofsPerNode;
    const localIndex rowLength = src.globalRowLength( row );
    srcIndices.resize( rowLength );
    srcValues.resize( rowLength );

    src.getRowCopy( row, srcIndices, srcValues );

    localIndex k = 0;
    for( localIndex col = 0; col < rowLength; ++col )
    {
      const globalIndex colComponent = srcIndices[col] % dofsPerNode;
      if( rowComponent == colComponent )
      {
        dstValues[k] = srcValues[col];
        dstIndices[k] = srcIndices[col];
        k++;
      }
    }
    dst.insert( row, dstIndices.data(), dstValues.data(), k );
  }
  dst.close();
}

}  // namespace LAIHelperFunctions

}  // namespace geosx

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LAIHELPERFUNCTIONS_HPP_*/
