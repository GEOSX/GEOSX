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
 * @file BlockPreconditioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_
#define GEOSX_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/solvers/PreconditionerBase.hpp"
#include "linearAlgebra/utilities/BlockOperator.hpp"
#include "linearAlgebra/utilities/BlockVector.hpp"

namespace geosx
{
/**
 * @brief Type of Schur complement approximation used
 *
 * @todo Need more descriptive names for options
 */
enum class SchurComplementOption
{
  None,                  //!< No Schur complement - just block-GS/block-Jacobi preconditioner
  FirstBlockDiagonal,    //!< Approximate first block with its diagonal
  RowsumDiagonalProbing  //!< Rowsum-preserving diagonal approximation constructed with probing
};

/**
 * @brief Type of block row scaling to apply
 */
enum class BlockScalingOption
{
  None,           //!< No scaling
  FrobeniusNorm,  //!< Equilibrate Frobenius norm of the diagonal blocks
  UserProvided    //!< User-provided scaling
};

/**
 * @brief Shape of the block preconditioner
 */
enum class BlockShapeOption
{
  Diagonal,             //!< (D)^{-1}
  UpperTriangular,      //!< (DU)^{-1}
  LowerUpperTriangular  //!< (LDU)^{-1}
};

/*
 * Since formulas in Doxygen are broken with 1.8.13 and certain versions of ghostscript,
 * keeping this documentation in a separate comment block for now. Should be moved into
 * documentation of BlockPreconditioner class.
 *
 * This class implements a 2x2 block preconditioner of the form:
 * @f$
 * M^{-1} =
 * \begin{bmatrix}
 * I_{1} & -M_{1}^{-1}A_{12} \\
 *       &  I_{2}
 * \end{bmatrix}
 * \begin{bmatrix}
 * M_{1}^{-1} &            \\
 *            & M_{2}^{-1}
 * \end{bmatrix}
 * \begin{bmatrix}
 *  I_{1}            &       \\
 * -A_{21}M_{1}^{-1} & I_{2}
 * \end{bmatrix}
 * @f$
 * with
 * @f$ M_1^{-1} ~= A_{11}^{-1} @f$, @f$ M_2^{-1} ~= S_{22}^{-1} @f$, and
 * @f$ S_{22} ~= A_{22} - A_{21} A_{11}^{-1} A_{12} @f$ being the (optional) approximate Schur complement.
 *
 * The first step (predictor solve with @f$ M_{1} @f$) is optional. Enabling it results in a more powerful
 * preconditioner, but involves two applications of @f$ M_{1}^{-1} @f$ instead of just one on each solve.
 *
 * The user provides individual block solvers @f$ M_{1} @f$ and @f$ M_{2} @f$ as well as
 * a description of the block split of monolithic matrix in terms of DOF components.
 */

/**
 * @brief General 2x2 block preconditioner.
 * @tparam LAI type of linear algebra interface providing matrix/vector types
 */
template< typename LAI >
class BlockPreconditioner : public PreconditionerBase< LAI >
{
public:
  /// Alias for the base type
  using Base = PreconditionerBase< LAI >;

  /// Alias for the vector type
  using Vector = typename Base::Vector;

  /// Alias for the matrix type
  using Matrix = typename Base::Matrix;

  /**
   * @brief Constructor.
   * @param shapeOption preconditioner block shape
   * @param schurOption type of Schur complement approximation to use
   * @param scalingOption type of scaling to apply to blocks
   */
  explicit BlockPreconditioner( BlockShapeOption const shapeOption,
                                SchurComplementOption const schurOption,
                                BlockScalingOption const scalingOption );

  /**
   * @brief Destructor.
   */
  virtual ~BlockPreconditioner() override;

  /**
   * @brief Setup data for one of the two blocks.
   * @param blockIndex index of the block to set up
   * @param blockDofs choice of DoF components (from a monolithic system)
   * @param solver instance of the inner preconditioner for the block
   * @param scaling user-provided row scaling coefficient for this block
   *
   * @note While not strictly required, it is generally expected that @p blockDofs
   *       of the two blocks are non-overlapping and their union includes all
   *       DoF components in the monolithic system.
   */
  void
  setupBlock( localIndex const blockIndex,
              std::vector< DofManager::SubComponent > blockDofs,
              std::unique_ptr< PreconditionerBase< LAI > > solver,
              real64 const scaling = 1.0 );

  /**
   * @name PreconditionerBase interface methods
   */
  ///@{

  using PreconditionerBase< LAI >::compute;

  virtual void
  compute( Matrix const & mat, DofManager const & dofManager ) override;

  virtual void
  apply( Vector const & src, Vector & dst ) const override;

  virtual void
  clear() override;

  ///@}

private:
  /**
   * @brief Initialize/resize internal data structures for a new linear system.
   * @param mat the new system matrix
   * @param dofManager the new dof manager
   */
  void
  reinitialize( Matrix const & mat, DofManager const & dofManager );

  /**
   * @brief Apply block scaling to system blocks (which must be already extracted).
   */
  void
  applyBlockScaling();

  /**
   * @brief Compute and apply the Schur complement to (1,1)-block.
   */
  void
  computeSchurComplement();

  /// Shape of the block preconditioner
  BlockShapeOption m_shapeOption;

  /// Type of Schur complement to construct
  SchurComplementOption m_schurOption;

  /// Whether to scale blocks to equilibrate norms
  BlockScalingOption m_scalingOption;

  /// Description of dof components making up each of the two main blocks
  std::array< std::vector< DofManager::SubComponent >, 2 > m_blockDofs;

  /// Restriction operators for each sub-block
  std::array< Matrix, 2 > m_restrictors;

  /// Prolongation operators for each sub-block
  std::array< Matrix, 2 > m_prolongators;

  /// Matrix blocks
  BlockOperator< Vector, Matrix > m_matBlocks;

  /// Individual block preconditioners
  std::array< std::unique_ptr< PreconditionerBase< LAI > >, 2 > m_solvers;

  /// Scaling of each block
  std::array< real64, 2 > m_scaling;

  /// Internal vector of block residuals
  mutable BlockVector< Vector > m_rhs;

  /// Internal vector of block solutions
  mutable BlockVector< Vector > m_sol;
};

}  //namespace geosx

#endif  //GEOSX_LINEARALGEBRA_SOLVERS_BLOCKPRECONDITIONER_HPP_
