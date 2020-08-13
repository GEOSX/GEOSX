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
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_
#define GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
/**
 * @brief Set of parameters for a linear solver or preconditioner.
 *
 * This class holds a simple tree of linear solver options.
 * They are set to default values, but can be overwritten as needed.
 */
class LinearSolverParameters
{
public:
  integer logLevel = 0;  ///< Output level [0=none, 1=basic, 2=everything]
  string solverType =
    "direct";  ///< Solver type [direct, cg, gmres, bicgstab, preconditioner]
  string preconditionerType =
    "iluk";  ///< Preconditioner type [none, iluk, ilut, amg, mgr, block]
  integer dofsPerNode =
    1;  ///< Dofs per node (or support location) for non-scalar problems
  bool isSymmetric =
    false;  ///< Whether input matrix is symmetric (may affect choice of scheme)

  /// Krylov-method parameters
  struct Krylov
  {
    real64 relTolerance =
      1e-6;  ///< Relative convergence tolerance for iterative solvers
    integer maxIterations =
      200;  ///< Max iterations before declaring convergence failure
    integer maxRestart =
      200;                           ///< Max number of vectors in Krylov basis before restarting
    integer useAdaptiveTol = false;  ///< Use Eisenstat-Walker adaptive tolerance
    real64 weakestTol =
      1e-3;  ///< Weakest allowed tolerance when using adaptive method
  } krylov;  ///< Krylov-method parameter struct

  /// Matrix-scaling parameters
  struct Scaling
  {
    integer useRowScaling = false;  ///< Apply row scaling
    integer useRowColScaling =
      false;  ///< Apply row and column scaling (not yet implemented)
  } scaling;  ///< Matrix-scaling parameter struct

  /// Algebraic multigrid parameters
  struct AMG
  {
    integer maxLevels = 20;               ///< Maximum number of coarsening levels
    string cycleType = "V";               ///< AMG cycle type
    string smootherType = "gaussSeidel";  ///< Smoother type
    string coarseType = "direct";         ///< Coarse-level solver/smoother
    integer numSweeps = 2;                ///< Number of smoother sweeps
    string preOrPostSmoothing =
      "both";  ///< Pre and/or post smoothing [pre,post,both]
    real64 threshold =
      0.0;  ///< Threshold for "strong connections" (for classical and
            ///< smoothed-aggregation AMG)
    integer separateComponents =
      false;  ///< Apply a separate component filter before AMG construction
    string nullSpaceType =
      "constantModes";  ///< Null space type [constantModes,rigidBodyModes]
  } amg;                ///< Algebraic Multigrid (AMG) parameters

  /// Multigrid reduction parameters
  struct MGR
  {
    string strategy;  ///< Predefined MGR solution strategy (solver specific)
    integer separateComponents =
      false;                       ///< Apply a separate displacement component (SDC) filter before AMG construction
    string displacementFieldName;  ///< Displacement field name need for SDC filter
  } mgr;                           ///< Multigrid reduction (MGR) parameters

  /// Incomplete factorization parameters
  struct ILU
  {
    integer fill = 0;        ///< Fill level
    real64 threshold = 0.0;  ///< Dropping threshold
  } ilu;                     ///< Incomplete factorization parameter struct

  /// Domain decomposition parameters
  struct DD
  {
    integer overlap = 0;  ///< Ghost overlap
  } dd;                   ///< Domain decomposition parameter struct
};

} /* namespace geosx */

#endif /*GEOSX_LINEARALGEBRA_UTILITIES_LINEARSOLVERPARAMETERS_HPP_ */
