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
 * LAIHelperFunctions.cpp
 *
 */

#include "LAIHelperFunctions.hpp"

namespace geosx
{
namespace LAIHelperFunctions
{
void
CreatePermutationMatrix( NodeManager const * const nodeManager,
                         localIndex const nRows,
                         localIndex const nCols,
                         int const nDofPerNode,
                         string const dofKey,
                         ParallelMatrix & permutationMatrix )
{
  /* Crates a permutation matrix for a given nodal variable specified by the DofKey. It consider that nDofPerNode
   * dofs are associated with each node (e.g., nDofPerNode = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the nodes.
   */

  // Create permuation matrix based on size provided.
  permutationMatrix.createWithLocalSize( nRows, nCols, 1, MPI_COMM_GEOSX );
  permutationMatrix.open();

  arrayView1d< globalIndex const > const & DofNumber =
    nodeManager->getReference< globalIndex_array >( dofKey );

  arrayView1d< globalIndex const > const & localToGlobal =
    nodeManager->localToGlobalMap();

  for( localIndex a = 0; a < nodeManager->size(); ++a )
  {
    if( DofNumber[a] >= 0 )
    {
      for( int d = 0; d < nDofPerNode; ++d )
      {
        globalIndex const rowIndex = localToGlobal[a] * nDofPerNode + d;
        globalIndex const columnIndex = DofNumber[a] + d;

        permutationMatrix.insert( rowIndex, columnIndex, 1.0 );
      }
    }
  }
  permutationMatrix.close();
  permutationMatrix.set( 1 );
}

void
CreatePermutationMatrix( ElementRegionManager const * const elemManager,
                         localIndex const nRows,
                         localIndex const nCols,
                         int const nDofPerCell,
                         string const DofKey,
                         ParallelMatrix & permutationMatrix )
{
  /* Crates a permutation matrix for a given cell centered variable specified by the DofKey. It consider that
     nDofPerNode
   * dofs are associated with each node (e.g., nDofPerNode = 3 for the displacement).
   *
   * The permutation matrix maps from the dofs ordering set by the DOFManager to the ordering based on the global
   * indexes of the cells.
   */

  // Create permuation matrix based on size provided.
  permutationMatrix.createWithLocalSize( nRows, nCols, 1, MPI_COMM_GEOSX );
  permutationMatrix.open();

  elemManager->forElementSubRegions< ElementSubRegionBase >(
    [&]( ElementSubRegionBase const & elementSubRegion ) {
      localIndex const numElems = elementSubRegion.size();
      arrayView1d< globalIndex const > const & DofNumber =
        elementSubRegion.getReference< array1d< globalIndex > >( DofKey );

      arrayView1d< globalIndex const > const & localToGlobal =
        elementSubRegion.localToGlobalMap();

      for( localIndex k = 0; k < numElems; ++k )
      {
        if( DofNumber[k] >= 0 )
        {
          for( int d = 0; d < nDofPerCell; ++d )
          {
            globalIndex const rowIndex = localToGlobal[k] * nDofPerCell + d;
            globalIndex const columnIndex = DofNumber[k] * nDofPerCell + d;

            permutationMatrix.insert( rowIndex, columnIndex, 1.0 );
          }
        }
      }
    } );
  permutationMatrix.close();
  permutationMatrix.set( 1 );
}

ParallelVector
PermuteVector( ParallelVector const & vector,
               ParallelMatrix const & permutationMatrix )
{
  ParallelVector permutedVector;

  permutedVector.createWithLocalSize( vector.localSize(), MPI_COMM_GEOSX );

  permutationMatrix.apply( vector, permutedVector );

  return permutedVector;
}

/*
 *  permutedMatrix = permutationMatrix * matrix * permutationMatrix^T;
 */
ParallelMatrix
PermuteMatrix( ParallelMatrix const & matrix,
               ParallelMatrix const & permutationMatrix )
{
  ParallelMatrix temp;
  ParallelMatrix permutedMatrix;
  // The value 24 is hardcoded and should probably be changed (It s fine for displacement).
  temp.createWithLocalSize( matrix.numLocalRows(),
                            matrix.numLocalCols(),
                            24,
                            MPI_COMM_GEOSX );

  permutedMatrix.createWithLocalSize( matrix.numLocalRows(),
                                      matrix.numLocalCols(),
                                      24,
                                      MPI_COMM_GEOSX );

  permutationMatrix.multiply( matrix, temp );
  permutationMatrix.rightMultiplyTranspose( temp, permutedMatrix );

  return permutedMatrix;
}

/*
 *  permutedMatrix = permutationMatrixLeft * matrix * permutationMatrixRight^T;
 */
ParallelMatrix
PermuteMatrix( ParallelMatrix const & matrix,
               ParallelMatrix const & permutationMatrixLeft,
               ParallelMatrix const & permutationMatrixRight )
{
  ParallelMatrix temp;
  ParallelMatrix permutedMatrix;

  temp.createWithLocalSize( matrix.numLocalRows(),
                            matrix.numLocalCols(),
                            24,
                            MPI_COMM_GEOSX );

  permutedMatrix.createWithLocalSize( matrix.numLocalRows(),
                                      matrix.numLocalCols(),
                                      24,
                                      MPI_COMM_GEOSX );

  permutationMatrixLeft.multiply( matrix, temp );
  permutationMatrixRight.rightMultiplyTranspose( temp, permutedMatrix );

  return permutedMatrix;
}

}  // namespace LAIHelperFunctions

}  // namespace geosx
