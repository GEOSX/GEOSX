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

#include "NonlinearSolverParameters.hpp"

namespace geosx
{
using namespace dataRepository;

NonlinearSolverParameters::NonlinearSolverParameters( std::string const & name,
                                                      Group * const parent ) :
  Group( name, parent )
{
  setInputFlags( InputFlags::OPTIONAL );

  // This enables logLevel filtering
  enableLogLevelInput();

  registerWrapper( viewKeysStruct::lineSearchActionString, &m_lineSearchAction )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "How the line search is to be used. Options are: \n "
      "0 - Do not use line search.\n"
      "1 - Use line search. Allow exit from line search without achieving "
      "smaller residual than starting residual.\n"
      "2 - Use line search. If smaller residual than starting resdual is not "
      "achieved, cut time step.\n" );

  registerWrapper( viewKeysStruct::lineSearchMaxCutsString, &m_lineSearchMaxCuts )
    ->setApplyDefaultValue( 4 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Maximum number of line search cuts." );

  registerWrapper( viewKeysStruct::lineSearchCutFactorString,
                   &m_lineSearchCutFactor )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Line search cut factor. For instance, a value of 0.5 will result in the "
      "effective application of"
      " the last solution by a factor of (0.5, 0.25, 0.125, ...)" );

  registerWrapper( viewKeysStruct::newtonTolString, &m_newtonTol )
    ->setApplyDefaultValue( 1.0e-6 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "The required tolerance in order to exit the Newton iteration loop." );

  registerWrapper( viewKeysStruct::newtonMaxIterString, &m_maxIterNewton )
    ->setApplyDefaultValue( 5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Maximum number of iterations that are allowed in a Newton loop." );

  registerWrapper( viewKeysStruct::newtonMinIterString, &m_minIterNewton )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Minimum number of iterations that are required before exiting the "
      "Newton loop." );

  registerWrapper( viewKeysStruct::newtonNumIterationsString,
                   &m_numNewtonIterations )
    ->setApplyDefaultValue( 0 )
    ->setDescription( "Number of Newton's iterations." );

  registerWrapper( viewKeysStruct::allowNonConvergedString, &m_allowNonConverged )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Allow non-converged solution to be accepted. "
      "(i.e. exit from the Newton loop without achieving the desired "
      "tolerance)" );

  registerWrapper( viewKeysStruct::dtCutIterLimString, &m_dtCutIterLimit )
    ->setApplyDefaultValue( 0.7 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Fraction of the Max Newton iterations above which the solver asks for "
      "the time-step to be cut for the next dt." );

  registerWrapper( viewKeysStruct::dtIncIterLimString, &m_dtIncIterLimit )
    ->setApplyDefaultValue( 0.4 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Fraction of the Max Newton iterations below which the solver asks for "
      "the time-step to be doubled for the next dt." );

  registerWrapper( viewKeysStruct::timeStepCutFactorString, &m_timeStepCutFactor )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Factor by which the time step will be cut if a timestep cut is "
      "required." );

  registerWrapper( viewKeysStruct::maxTimeStepCutsString, &m_maxTimeStepCuts )
    ->setApplyDefaultValue( 2 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Max number of time step cuts" );

  registerWrapper( viewKeysStruct::maxSubStepsString, &m_maxSubSteps )
    ->setApplyDefaultValue( 10 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Maximum number of time sub-steps allowed for the solver" );
}

void
NonlinearSolverParameters::PostProcessInput()
{
  if( m_dtCutIterLimit <= m_dtIncIterLimit )
  {
    GEOSX_ERROR( " dtIncIterLimit should be smaller than dtCutIterLimit!!" );
  }
}

REGISTER_CATALOG_ENTRY( Group,
                        NonlinearSolverParameters,
                        std::string const &,
                        Group * const )

} /* namespace geosx */
