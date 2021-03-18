/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file testLinearAlgebraUtils.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_
#define GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_

#include "common/DataTypes.hpp"
#include "common/initializeEnvironment.hpp"
#include "common/MpiWrapper.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "mainInterface/ProblemManager.hpp"
#include "mesh/MeshManager.hpp"

#include <gtest/gtest.h>

namespace geosx
{
namespace testing
{

/**
 * @brief Simple scope-based manager for initialization/finalization in tests.
 */
class LinearAlgebraTestScope
{
public:

  LinearAlgebraTestScope( int argc, char * * argv )
  {
    ::testing::InitGoogleTest( &argc, argv );
    geosx::setupEnvironment( argc, argv );
    setupLAI();
  }

  ~LinearAlgebraTestScope()
  {
    finalizeLAI();
    geosx::cleanupEnvironment();
  }
};

/**
 * @name Utility functions for linear algebra unit tests.
 * @brief Functions used to construct useful matrices in the test files.
 */
///@{

/**
 * @brief Construct a square zero matrix.
 * @tparam MATRIX type of matrix
 * @param comm MPI communicator
 * @param N global size of the natrix
 * @param Z the output matrix
 */
template< typename MATRIX >
void computeZero( MPI_Comm comm,
                  globalIndex N,
                  MATRIX & Z )
{
  Z.createWithGlobalSize( N, 0, comm );
  Z.open();
  Z.close();
}

/**
 * @brief Compute the 2D Laplace operator
 * @tparam MATRIX type of matrix
 * @param comm      MPI communicator.
 * @param n         size of the nxn mesh for the square 2D Laplace operator matrix. Matrix size will be N=n^2.
 * @param laplace2D the output matrix
 *
 * This function computes the matrix corresponding to a 2D Laplace operator. These
 * matrices arise from a classical finite volume formulation on a cartesian mesh
 * (5-point stencil).
 */
template< typename MATRIX >
void compute2DLaplaceOperator( MPI_Comm comm,
                               globalIndex n,
                               MATRIX & laplace2D )
{
  // total dofs = n^2
  globalIndex N = n * n;

  // Create a matrix of global size N with 5 non-zeros per row
  laplace2D.createWithGlobalSize( N, 5, comm );

  // Allocate arrays to fill the matrix (values and columns)
  real64 values[5];
  globalIndex cols[5];

  // Open the matrix
  laplace2D.open();

  // Loop over rows to fill the matrix
  for( globalIndex i = laplace2D.ilower(); i < laplace2D.iupper(); i++ )
  {
    // Re-set the number of non-zeros for row i to 0.
    localIndex nnz = 0;

    // The left -n: position i-n
    if( i - n >= 0 )
    {
      cols[nnz] = i - n;
      values[nnz] = -1.0;
      nnz++;
    }

    // The left -1: position i-1
    if( i - 1 >= 0 )
    {
      cols[nnz] = i - 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the diagonal: position i
    cols[nnz] = i;
    values[nnz] = 4.0;
    nnz++;

    // The right +1: position i+1
    if( i + 1 < N )
    {
      cols[nnz] = i + 1;
      values[nnz] = -1.0;
      nnz++;
    }

    // The right +n: position i+n
    if( i + n < N )
    {
      cols[nnz] = i + n;
      values[nnz] = -1.0;
      nnz++;
    }

    // Set the values for row i
    laplace2D.insert( i, cols, values, nnz );
  }

  // Close the matrix (make data contiguous in memory)
  laplace2D.close();

}

/**
 * @brief Compute a 1st order FEM local stiffness matrix for a quad element.
 * @param hx element width
 * @param hy element height
 * @param E Young's modulus
 * @param nu Poisson ratio
 * @param Ke the output stiffness matrix
 */
inline void Q12d_local( real64 const & hx,
                        real64 const & hy,
                        real64 const & E,
                        real64 const & nu,
                        arraySlice2d< real64 > const & Ke )
{
  real64 fac = E / ( 1. - 2. * nu ) / ( 1. + nu );

  // Populate stiffness matrix

  // --- Fill diagonal entries
  real64 Dxx = ( fac * hx * ( 1. - 2. * nu ) ) / ( 6. * hy )
               - ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  real64 Dyy = ( fac * hy * ( 1. - 2. * nu ) ) / ( 6. * hx )
               - ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );
  for( localIndex i = 0; i < 8; i += 2 )
  {
    Ke( i, i ) = Dxx;
    Ke( i + 1, i + 1 ) = Dyy;
  }

  // --- Fill upper triangular part
  // --- --- Ke( 0, 1:7 )
  Ke( 0, 1 ) = fac / 8.;
  Ke( 0, 2 ) = ( fac * hx * ( 1. - 2. * nu ) ) / ( 12. * hy )
               + ( fac * hy * ( -1. + nu ) ) / ( 3. * hx );
  Ke( 0, 3 ) = ( fac * ( -1 + 4. * nu ) ) / 8.;
  Ke( 0, 4 ) = ( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
               + ( fac * hx * ( -1. + 2. * nu ) ) / ( 12. * hy );
  Ke( 0, 5 ) = -Ke( 0, 1 );
  Ke( 0, 6 ) = -( fac * hy * ( -1. + nu ) ) / ( 6. * hx )
               + ( fac * hx * ( -1. + 2. * nu ) ) / ( 6. * hy );
  Ke( 0, 7 ) = -( fac * ( -1. + 4. * nu ) ) / 8.;

  // --- --- Ke( 1, 2:7 )
  Ke( 1, 2 ) = Ke( 0, 7 );
  Ke( 1, 3 ) = -( fac * ( hy * hy * ( 1. - 2. * nu ) + hx * hx * ( -1. + nu ) ) ) / ( 6. * hx * hy );
  Ke( 1, 4 ) = Ke( 0, 5 );
  Ke( 1, 5 ) = ( fac * hx * ( -1. + nu ) ) / ( 6. * hy ) + ( fac * hy * ( -1. + 2. * nu ) ) / ( 12. * hx );
  Ke( 1, 6 ) = Ke( 0, 3 );
  Ke( 1, 7 ) = ( fac * hy * ( 1. - 2. * nu ) ) / ( 12. * hx )
               + ( fac * hx * ( -1. + nu ) ) / ( 3. * hy );

  // --- --- Ke( 2, 3:7 )
  Ke( 2, 3 ) = Ke( 0, 5 );
  Ke( 2, 4 ) = Ke( 0, 6 );
  Ke( 2, 5 ) = Ke( 1, 6 );
  Ke( 2, 6 ) = ( fac * hy * ( -1 + nu ) ) / ( 6. * hx ) + ( fac * hx * ( -1. + 2. * nu ) ) / ( 12. * hy );
  Ke( 2, 7 ) = Ke( 0, 1 );

  // --- --- Ke( 3, 4:7 )
  Ke( 3, 4 ) = Ke( 1, 2 );
  Ke( 3, 5 ) = Ke( 1, 7 );
  Ke( 3, 6 ) = Ke( 0, 1 );
  Ke( 3, 7 ) = Ke( 1, 5 );

  // --- --- Ke( 4, 5:7 )
  Ke( 4, 5 ) = Ke( 0, 1 );
  Ke( 4, 6 ) = Ke( 0, 2 );
  Ke( 4, 7 ) = Ke( 0, 3 );

  // --- --- Ke( 5, 6:7 )
  Ke( 5, 6 ) = Ke( 0, 7 );
  Ke( 5, 7 ) = Ke( 1, 3 );

  // --- --- Ke( 6, 7 )
  Ke( 6, 7 ) = Ke( 0, 5 );

  // --- Fill lower triangular part
  for( localIndex i = 1; i < 8; ++i )
  {
    for( localIndex j = 0; j < i; ++j )
    {
      Ke( i, j ) = Ke( j, i );
    }
  }
}

/**
 * @brief Compute the 2D elasticity (plane strain) operator
 * @tparam MATRIX type of output matrix
 * @param comm MPI communicator
 * @param domainSizeX domain size in the X-direction
 * @param domainSizeY domain size in the Y-direction
 * @param nCellsX number of cells in the X-direction
 * @param nCellsY number of cells in the Y-direction
 * @param youngModulus Young's modulus value (same for all cells)
 * @param poissonRatio Poisson's ratio value (same for all cells)
 * @param elasticity2D the output matrix
 *
 * This function computes the matrix corresponding to a 2D elasticity operator,
 * assuming plane strain conditions, based on a Q1 finite element discretization.
 * The discretized domain has dimensions domainSizeX * domainSizeY.  A regular grid
 * consting of nCellsX * nCellsY cells is used. The medium is characterized by
 * homogeneous Young's modulus and Poisson's ratio. The assembled matrix is
 * singular, meaning that Dirichlet boundary conditions have not been enforced.
 */
template< typename MATRIX >
void compute2DElasticityOperator( MPI_Comm const comm,
                                  real64 const domainSizeX,
                                  real64 const domainSizeY,
                                  globalIndex const nCellsX,
                                  globalIndex const nCellsY,
                                  real64 const youngModulus,
                                  real64 const poissonRatio,
                                  MATRIX & elasticity2D )
{
  localIndex const rank = LvArray::integerConversion< localIndex >( MpiWrapper::commRank( comm ) );
  localIndex const nproc = LvArray::integerConversion< localIndex >( MpiWrapper::commSize( comm ) );

  GEOSX_ERROR_IF( nCellsY < nproc, "Less than one cell row per processor is not supported" );

  real64 const hx = domainSizeX / nCellsX;
  real64 const hy = domainSizeY / nCellsY;

  // Compute cell partitioning
  globalIndex const nCells = nCellsX * nCellsY;
  localIndex const nLocalCells = LvArray::integerConversion< localIndex >( nCellsY / nproc * nCellsX );
  localIndex const nExtraCells = LvArray::integerConversion< localIndex >( nCells ) - nLocalCells * nproc;
  globalIndex const iCellLower = rank * nLocalCells + ( rank == 0 ? 0 : nExtraCells );
  globalIndex const iCellUpper = iCellLower + nLocalCells + ( rank == 0 ? 0 : nExtraCells );

  // Compute node partitioning
  globalIndex const iNodeLower = iCellLower / nCellsX + iCellLower;
  globalIndex const iNodeUpper = iCellUpper / nCellsX + iCellUpper + ( rank == nproc - 1 ? nCellsX + 1 : 0 );
  localIndex const nLocalNodes = iNodeUpper - iNodeLower;

  // Construct local stiffness matrix (same for all cells)
  stackArray2d< real64, 8 * 8 > Ke( 8, 8 );
  Q12d_local( hx, hy, youngModulus, poissonRatio, Ke );

  // Create a matrix of global size N with at most 18 non-zeros per row
  elasticity2D.createWithLocalSize( nLocalNodes * 2, 18, comm );

  // Open the matrix
  elasticity2D.open();

  // Loop over grid cells
  stackArray1d< globalIndex, 4 > cellNodes( 4 );
  stackArray1d< globalIndex, 8 > localDofIndex( 8 );

  for( globalIndex iCell = iCellLower; iCell < iCellUpper; ++iCell )
  {

    // Compute local DOF global indeces
    cellNodes( 0 ) = iCell / nCellsX + iCell;
    cellNodes( 1 ) = cellNodes( 0 ) + 1;
    cellNodes( 3 ) = cellNodes( 1 ) + nCellsX;
    cellNodes( 2 ) = cellNodes( 3 ) + 1;
    for( localIndex i = 0; i < 4; ++i )
    {
      localDofIndex( 2 * i ) = cellNodes( i ) * 2;
      localDofIndex( 2 * i + 1 ) = localDofIndex( 2 * i ) + 1;
    }

    // Assemble local stiffness matrix and right-hand side
    elasticity2D.insert( localDofIndex.data(), localDofIndex.data(), Ke.data(), 8, 8 );
  }

  // Close the matrix
  elasticity2D.close();
}

///@}

} // namespace testing
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_UNITTESTS_TESTLINEARALGEBRAUTILS_HPP_
