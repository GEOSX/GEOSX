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
 * @file BlockPreconditioner.cpp
 */

#include "BlockPreconditioner.hpp"
#include "codingUtilities/Utilities.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/BlockVectorWrapper.hpp"

#include <numeric>

namespace geosx
{
template< typename LAI >
BlockPreconditioner< LAI >::BlockPreconditioner(
  BlockShapeOption const shapeOption,
  SchurComplementOption const schurOption,
  BlockScalingOption const scalingOption ) :
  Base(),
  m_shapeOption( shapeOption ),
  m_schurOption( schurOption ),
  m_scalingOption( scalingOption ),
  m_matBlocks( 2, 2 ),
  m_solvers {},
  m_scaling { 1.0, 1.0 },
  m_rhs( 2 ),
  m_sol( 2 )
{}

template< typename LAI >
BlockPreconditioner< LAI >::~BlockPreconditioner() = default;

template< typename LAI >
void
BlockPreconditioner< LAI >::reinitialize( Matrix const & mat,
                                          DofManager const & dofManager )
{
  MPI_Comm const & comm = mat.getComm();

  if( m_blockDofs[1].empty() )
  {
    m_blockDofs[1] = dofManager.filterDofs( m_blockDofs[0] );
  }

  for( localIndex i = 0; i < 2; ++i )
  {
    dofManager.makeRestrictor( m_blockDofs[i], comm, false, m_restrictors[i] );
    dofManager.makeRestrictor( m_blockDofs[i], comm, true, m_prolongators[i] );
    m_rhs( i ).createWithLocalSize( m_restrictors[i].numLocalRows(), comm );
    m_sol( i ).createWithLocalSize( m_restrictors[i].numLocalRows(), comm );
  }
}

template< typename LAI >
void
BlockPreconditioner< LAI >::setupBlock(
  localIndex const blockIndex,
  std::vector< DofManager::SubComponent > blockDofs,
  std::unique_ptr< PreconditionerBase< LAI > > solver,
  real64 const scaling )
{
  GEOSX_LAI_ASSERT_GT( 2, blockIndex );
  GEOSX_LAI_ASSERT( solver );
  GEOSX_LAI_ASSERT( !blockDofs.empty() );
  GEOSX_LAI_ASSERT_GT( scaling, 0.0 );

  m_blockDofs[blockIndex] = std::move( blockDofs );
  m_solvers[blockIndex] = std::move( solver );
  m_scaling[blockIndex] = scaling;
}

template< typename LAI >
void
BlockPreconditioner< LAI >::applyBlockScaling()
{
  if( m_scalingOption != BlockScalingOption::None )
  {
    if( m_scalingOption == BlockScalingOption::FrobeniusNorm )
    {
      real64 norms[2] = { m_matBlocks( 0, 0 ).normFrobenius(),
                          m_matBlocks( 1, 1 ).normFrobenius() };
      m_scaling[0] = std::min( norms[1] / norms[0], 1.0 );
      m_scaling[1] = std::min( norms[0] / norms[1], 1.0 );
    }

    for( localIndex i = 0; i < 2; ++i )
    {
      for( localIndex j = 0; j < 2; ++j )
      {
        m_matBlocks( i, j ).scale( m_scaling[i] );
      }
    }
  }
}

template< typename LAI >
void
BlockPreconditioner< LAI >::computeSchurComplement()
{
  switch( m_schurOption )
  {
    case SchurComplementOption::None:
    {
      // nothing to do
      break;
    }
    case SchurComplementOption::FirstBlockDiagonal:
    {
      // In this case, the sparsity pattern of triple product can be denser
      // than that of (1,1)-block. Therefore, we have to add (1,1)-block to
      // the triple product result, not the other way around.
      // TODO: This is suboptimal, since we have to create multiple matrix copies.
      m_matBlocks( 0, 0 ).extractDiagonal( m_rhs( 0 ) );
      m_rhs( 0 ).reciprocal();
      Matrix mat01 = m_matBlocks( 0, 1 );  // make a copy in order to scale
      mat01.leftScale( m_rhs( 0 ) );
      Matrix mat11 = m_matBlocks( 1, 1 );  // make a copy in order to add later
      m_matBlocks( 1, 0 ).multiply( mat01, m_matBlocks( 1, 1 ) );
      m_matBlocks( 1, 1 ).addEntries( mat11 );
      break;
    }
    case SchurComplementOption::RowsumDiagonalProbing:
    {
      m_sol( 1 ).set( -1.0 );
      m_matBlocks( 0, 1 ).apply( m_sol( 1 ), m_rhs( 0 ) );
      m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
      m_matBlocks( 1, 0 ).apply( m_sol( 0 ), m_rhs( 1 ) );
      m_matBlocks( 1, 1 ).addDiagonal( m_rhs( 1 ) );
      break;
    }
    default:
    {
      GEOSX_ERROR( "BlockPreconditioner: unsupported Schur complement option" );
    }
  }
}

template< typename LAI >
void
BlockPreconditioner< LAI >::compute( Matrix const & mat,
                                     DofManager const & dofManager )
{
  // Check that user has set block solvers
  GEOSX_LAI_ASSERT( m_solvers[0] != nullptr );
  GEOSX_LAI_ASSERT( m_solvers[1] != nullptr );

  // Compare old sizes vs new matris sizes.
  // A change in size indicates a new matrix structure.
  // This is done before Base::compute() since it overwrites old sizes.
  bool const newSize = !this->ready() ||
    mat.numGlobalRows() != this->numGlobalRows() ||
    mat.numGlobalCols() != this->numGlobalRows();

  Base::compute( mat, dofManager );

  // If the matrix size/structure has changed, need to resize internal LA objects and recompute restrictors.
  if( newSize )
  {
    reinitialize( mat, dofManager );
  }

  // Extract diagonal blocks
  mat.multiplyPtAP( m_prolongators[0], m_matBlocks( 0, 0 ) );
  mat.multiplyPtAP( m_prolongators[1], m_matBlocks( 1, 1 ) );

  // Extract off-diagonal blocks only if used
  if( m_schurOption != SchurComplementOption::None &&
      m_shapeOption != BlockShapeOption::Diagonal )
  {
    mat.multiplyRAP( m_restrictors[0], m_prolongators[1], m_matBlocks( 0, 1 ) );
    mat.multiplyRAP( m_restrictors[1], m_prolongators[0], m_matBlocks( 1, 0 ) );
  }

  applyBlockScaling();
  m_solvers[0]->compute( m_matBlocks( 0, 0 ), dofManager );
  computeSchurComplement();
  m_solvers[1]->compute( m_matBlocks( 1, 1 ), dofManager );
}

template< typename LAI >
void
BlockPreconditioner< LAI >::apply( Vector const & src, Vector & dst ) const
{
  m_restrictors[0].apply( src, m_rhs( 0 ) );
  m_restrictors[1].apply( src, m_rhs( 1 ) );

  for( localIndex i = 0; i < 2; ++i )
  {
    m_rhs( i ).scale( m_scaling[i] );
  }

  // Perform a predictor step by solving (0,0) block and subtracting from 1-block rhs
  if( m_shapeOption == BlockShapeOption::LowerUpperTriangular )
  {
    m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );
    m_matBlocks( 1, 0 ).residual( m_sol( 0 ), m_rhs( 1 ), m_rhs( 1 ) );
  }

  // Solve the (1,1) block modified via Schur complement
  m_solvers[1]->apply( m_rhs( 1 ), m_sol( 1 ) );

  // Update the 0-block rhs
  if( m_shapeOption != BlockShapeOption::Diagonal )
  {
    m_matBlocks( 0, 1 ).residual( m_sol( 1 ), m_rhs( 0 ), m_rhs( 0 ) );
  }

  // Solve the (0,0) block with the current rhs
  m_solvers[0]->apply( m_rhs( 0 ), m_sol( 0 ) );

  // Combine block solutions into global solution vector
  m_prolongators[0].apply( m_sol( 0 ), dst );
  m_prolongators[1].gemv( 1.0, m_sol( 1 ), 1.0, dst );
}

template< typename LAI >
void
BlockPreconditioner< LAI >::clear()
{
  Base::clear();
  for( localIndex i = 0; i < 2; ++i )
  {
    m_restrictors[i].reset();
    m_prolongators[i].reset();
    m_solvers[i]->clear();
    m_rhs( i ).reset();
    m_sol( i ).reset();
    for( localIndex j = 0; j < 2; ++j )
    {
      m_matBlocks( i, j ).reset();
    }
  }
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class BlockPreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class BlockPreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class BlockPreconditioner< PetscInterface >;
#endif

}  // namespace geosx
