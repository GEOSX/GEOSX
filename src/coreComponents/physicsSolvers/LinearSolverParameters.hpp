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
 * @file LinearSolverParameters.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_
#define GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

namespace geosx
{

/**
 * @brief Linear solver parameters with Group capabilities
 *
 * This class is a derived version of LinearSolverParameters with
 * dataRepository::Group capabilities to allow for XML input.
 *
 * Developer note: This class exposes LAI solver settings to external users.
 * As a general philosophy, only a subset of frequently tuned parameters should
 * be exposed.  Many advanced parameters can be set by the PhysicsSolver itself
 * since it has knowledge of the underlying problem (e.g. isSymmetric = true,
 * dofsPerNode = 3, etc.).  While we want to enable power users to tune the
 * solvers, most users prefer a short list of options with good default settings
 */
class LinearSolverParametersInput : public dataRepository::Group
{
public:

  LinearSolverParametersInput() = delete;

  /// Constructor
  LinearSolverParametersInput( string const & name, Group * const parent );

  /// Copy constructor
  LinearSolverParametersInput( LinearSolverParametersInput && ) = default;

  /// Destructor
  virtual ~LinearSolverParametersInput() override = default;

  /// Catalog name
  static string catalogName() { return "LinearSolverParameters"; }

  /// Postprocessing of input
  virtual void postProcessInput() override;

  LinearSolverParameters const & get() const
  { return m_parameters; }

  LinearSolverParameters & get()
  { return m_parameters; }

  /// Keys appearing in XML
  struct viewKeyStruct
  {
    /// Solver type key
    static constexpr char const * solverTypeString() { return "solverType"; }
    /// Preconditioner type key
    static constexpr char const * preconditionerTypeString() { return "preconditionerType"; }
    /// stop if error key
    static constexpr char const * stopIfErrorString() { return "stopIfError"; }

    /// direct solver check residual key
    static constexpr char const * directCheckResidualString() { return "directCheckResidual"; }
    /// direct solver equilibrate key
    static constexpr char const * directEquilString() { return "directEquil"; }
    /// direct solver columns permutation key
    static constexpr char const * directColPermString() { return "directColPerm"; }
    /// direct solver rows permutation key
    static constexpr char const * directRowPermString() { return "directRowPerm"; }
    /// direct solver replace tiny pivot key
    static constexpr char const * directReplTinyPivotString() { return "directReplTinyPivot"; }
    /// direct solver iterative refinement key
    static constexpr char const * directIterRefString() { return "directIterRef"; }
    /// direct solver parallelism key
    static constexpr char const * directParallelString() { return "directParallel"; }

    /// Krylov max iterations key
    static constexpr char const * krylovMaxIterString() { return "krylovMaxIter"; }
    /// Krylov max iterations key
    static constexpr char const * krylovMaxRestartString() { return "krylovMaxRestart"; }
    /// Krylov tolerance key
    static constexpr char const * krylovTolString() { return "krylovTol"; }
    /// Krylov adaptive tolerance key
    static constexpr char const * krylovAdaptiveTolString() { return "krylovAdaptiveTol"; }
    /// Krylov weakest tolerance key
    static constexpr char const * krylovWeakTolString() { return "krylovWeakestTol"; }

    /// AMG number of sweeps key
    static constexpr char const * amgNumSweepsString() { return "amgNumSweeps"; }
    /// AMG smoother type key
    static constexpr char const * amgSmootherString() { return "amgSmootherType"; }
    /// AMG coarse solver key
    static constexpr char const * amgCoarseString() { return "amgCoarseSolver"; }
    /// AMG threshold key
    static constexpr char const * amgThresholdString() { return "amgThreshold"; }
    /// AMG near null space type key
    static constexpr char const * amgNullSpaceTypeString() { return "amgNullSpaceType"; }

    static constexpr char const * amgCoarseningString()         { return "amgCoarseningType";           } ///< AMG coarsening key
    static constexpr char const * amgInterpolationString()      { return "amgInterpolationType";        }   ///< AMG interpolation key
    static constexpr char const * amgNumFunctionsString()       { return "amgNumFunctions";             }   ///< AMG threshold key
    static constexpr char const * amgAggresiveNumLevelsString() { return "amgAggresiveCoarseningLevels";}             ///< AMG threshold key
    /// ILU fill key
    static constexpr char const * iluFillString() { return "iluFill"; }
    /// ILU threshold key
    static constexpr char const * iluThresholdString() { return "iluThreshold"; }

    static constexpr char const * multiscaleBasisTypeString()               { return "multiscaleBasisType"; }
    static constexpr char const * multiscaleMaxLevelsString()               { return "multiscaleMaxLevels"; }
    static constexpr char const * multiscaleMinLocalDofString()             { return "multiscaleMinLocalDof"; }
    static constexpr char const * multiscaleMinGlobalDofString()            { return "multiscaleMinGlobalDof"; }
    static constexpr char const * multiscaleNumSmootherSweepsString()       { return "multiscaleNumSmootherSweeps"; }
    static constexpr char const * multiscalePreOrPostSmoothingString()      { return "multiscalePreOrPostSmoothing"; }
    static constexpr char const * multiscaleSmootherTypeString()            { return "multiscaleSmootherType"; }
    static constexpr char const * multiscaleBoundarySets()                  { return "multiscaleBoundarySets"; }
    static constexpr char const * multiscaleDebugLevel()                    { return "multiscaleDebugLevel"; }
    static constexpr char const * multiscaleCoarseTypeString()              { return "multiscaleCoarseType"; }
    static constexpr char const * multiscaleCoarseningPartitionTypeString() { return "multiscaleCoarseningPartitionType"; }
    static constexpr char const * multiscaleCoarseningRatioString()         { return "multiscaleCoarseningRatio"; }
    static constexpr char const * multiscaleCoarseningMetisMethodString()   { return "multiscaleCoarseningMetisMethod"; }
    static constexpr char const * multiscaleCoarseningMetisUfactorString()  { return "multiscaleCoarseningMetisUfactor"; }
    static constexpr char const * multiscaleMsrsbMaxIterString()            { return "multiscaleMsrsbMaxIter"; }
    static constexpr char const * multiscaleMsrsbToleranceString()          { return "multiscaleMsrsbTolerance"; }
    static constexpr char const * multiscaleMsrsbRelaxationString()         { return "multiscaleMsrsbRelaxation"; }
  };

private:

  LinearSolverParameters m_parameters;

};

} // namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_LINEARSOLVERPARAMETERS_HPP_
