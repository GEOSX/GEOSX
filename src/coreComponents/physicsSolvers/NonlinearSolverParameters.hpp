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

#ifndef GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_
#define GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_

#include "dataRepository/Group.hpp"

namespace geosx
{
/**
 * @class Holds parameters and status for execution of nonlinear solution schemes.
 */
class NonlinearSolverParameters : public dataRepository::Group
{
public:
  /**
   * @brief Default constructor.
   */
  NonlinearSolverParameters() = delete;

  /**
   * @brief Constructor
   * @param[in] name The name of the new instantiation of this Group.
   * @param[in] parent A pointer to the parent of this Group.
   */
  NonlinearSolverParameters( std::string const & name, Group * const parent );

  /**
   * @brief Default destructor
   */
  virtual ~NonlinearSolverParameters() = default;

  /**
   * @brief Default Move Constructor
   * @param The source object of the move.
   */
  NonlinearSolverParameters( NonlinearSolverParameters && ) = default;

  /**
   * @brief The name of this object in the catalog.
   * @return A string containing the name of this object in the catalog.
   * The CatalogName is the string that will result in the creation of a new
   * NonlinearSolverParameters object when calling
   * Group::GetCatalog()::Allocate().
   */
  static string
  CatalogName()
  {
    return "NonlinearSolverParameters";
  }

  virtual void
  PostProcessInput() override;

  struct viewKeysStruct
  {
    static constexpr auto lineSearchActionString = "lineSearchAction";
    static constexpr auto lineSearchMaxCutsString = "lineSearchMaxCuts";
    static constexpr auto lineSearchCutFactorString = "lineSearchCutFactor";

    static constexpr auto newtonTolString = "newtonTol";
    static constexpr auto newtonMaxIterString = "newtonMaxIter";
    static constexpr auto newtonMinIterString = "newtonMinIter";
    static constexpr auto newtonNumIterationsString =
      "newtonNumberOfIterations";
    static constexpr auto newtonSplitOperMaxIterString =
      "newtonSplitOperMaxIter";

    static constexpr auto allowNonConvergedString = "allowNonConverged";
    static constexpr auto dtCutIterLimString = "dtCutIterLimit";
    static constexpr auto dtIncIterLimString = "dtIncIterLimit";
    static constexpr auto maxSubStepsString = "maxSubSteps";
    static constexpr auto maxTimeStepCutsString = "maxTimeStepCuts";
    static constexpr auto minNumNewtonIterationsString =
      "minNumberOfNewtonIterations";
    static constexpr auto timeStepCutFactorString = "timestepCutFactor";

  } viewKeys;

  /**
   * @brief Calculates the upper limit for the number of iterations to allow a
   * cut to the next timestep.
   * @return The scaled value of the limit (m_dtCutIterLimit * m_maxIterNewton)
   */
  integer
  dtCutIterLimit() const
  {
    return std::ceil( m_dtCutIterLimit * m_maxIterNewton );
  }

  /**
   * @brief Calculates the lower limit for the number of iterations to force an
   * increase to the next timestep.
   * @return The scaled value of the limit (m_dtIncIterLimit * m_maxIterNewton)
   */
  integer
  dtIncIterLimit() const
  {
    return std::ceil( m_dtIncIterLimit * m_maxIterNewton );
  }

  /// Flag to apply a line search.
  integer m_lineSearchAction;

  /// The maximum number of line search cuts to attempt.
  integer m_lineSearchMaxCuts;

  /// The reduction factor for each line search cut.
  real64 m_lineSearchCutFactor;

  /// The tolerance for the nonlinear convergence check.
  real64 m_newtonTol;

  /// The maximum number of nonlinear iterations that are allowed.
  integer m_maxIterNewton;

  /// The minimum number of nonlinear iterations that may be applied.
  integer m_minIterNewton;

  /// The number of nonlinear iterations that have been exectued.
  integer m_numNewtonIterations;

  /// Flag to allow for a non-converged nonlinear solution and continue with the problem.
  integer m_allowNonConverged;

  /// Fraction of the Max Newton iterations above which the solver asks for the time-step to be cut for the next dt
  real64 m_dtCutIterLimit;

  /// Fraction of the Max Newton iterations below which the solver asks for the time-step to be doubled for the next dt
  real64 m_dtIncIterLimit;

  /// Maximum number of time sub-steps allowed for the solver
  integer m_maxSubSteps;

  /// Max number of time step cuts
  integer m_maxTimeStepCuts;

  /// Factor by which the time step will be cut if a timestep cut is required.
  real64 m_timeStepCutFactor;

  /// number of times that the time-step had to be cut
  integer m_numdtAttempts;
};

} /* namespace geosx */

#endif /* GEOSX_PHYSICSSOLVERS_NONLINEARSOLVERPARAMETERS_HPP_ */
