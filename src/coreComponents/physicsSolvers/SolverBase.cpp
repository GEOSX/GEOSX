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

#include "SolverBase.hpp"
#include "PhysicsSolverManager.hpp"

#include "common/TimingMacros.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/solvers/KrylovSolver.hpp"
#include "managers/DomainPartition.hpp"

namespace geosx
{
using namespace dataRepository;

SolverBase::SolverBase( std::string const & name, Group * const parent ) :
  ExecutableGroup( name, parent ),
  m_cflFactor(),
  m_maxStableDt { 1e99 },
  m_nextDt( 1e99 ),
  m_dofManager( name ),
  m_linearSolverParameters( groupKeyStruct::linearSolverParametersString, this ),
  m_nonlinearSolverParameters( groupKeyStruct::nonlinearSolverParametersString,
                               this )
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  // This enables logLevel filtering
  enableLogLevelInput();

  // This sets a flag to indicate that this object increments time
  this->SetTimestepBehavior( 1 );

  registerWrapper( viewKeyStruct::cflFactorString, &m_cflFactor )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Factor to apply to the `CFL condition "
      "<http://en.wikipedia.org/wiki/Courant-Friedrichs-Lewy_condition>`_"
      " when calculating the maximum allowable time step. Values should be in "
      "the interval (0,1] " );

  registerWrapper( viewKeyStruct::maxStableDtString, &m_maxStableDt )
    ->setApplyDefaultValue( 0.5 )
    ->setInputFlag( InputFlags::FALSE )
    ->setDescription( "Value of the Maximum Stable Timestep for this solver." );

  this
    ->registerWrapper( viewKeyStruct::discretizationString, &m_discretizationName )
    ->setApplyDefaultValue( "none" )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Name of discretization object (defined in the "
      ":ref:`NumericalMethodsManager`) to use for this "
      "solver. For instance, if this is a Finite Element Solver, the name of a "
      ":ref:`FiniteElement` "
      "should be specified. If this is a Finite Volume Method, the name of a "
      ":ref:`FiniteVolume` "
      "discretization should be specified." );

  registerWrapper( viewKeyStruct::targetRegionsString, &m_targetRegionNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "Allowable regions that the solver may be applied to. Note that this "
      "does not indicate that "
      "the solver will be applied to these regions, only that allocation will "
      "occur such that the "
      "solver may be applied to these regions. The decision about what regions "
      "this solver will be"
      "applied to rests in the EventManager." );

  registerWrapper( viewKeyStruct::initialDtString, &m_nextDt )
    ->setApplyDefaultValue( 1e99 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Initial time-step value required by the solver to the event manager." );

  RegisterGroup( groupKeyStruct::linearSolverParametersString,
                 &m_linearSolverParameters );
  RegisterGroup( groupKeyStruct::nonlinearSolverParametersString,
                 &m_nonlinearSolverParameters );

  m_localMatrix.setName( this->getName() + "/localMatrix" );
  m_localRhs.setName( this->getName() + "/localRhs" );
}

SolverBase::~SolverBase() = default;

Group *
SolverBase::CreateChild( string const & GEOSX_UNUSED_PARAM( childKey ),
                         string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

SolverBase::CatalogInterface::CatalogType &
SolverBase::GetCatalog()
{
  static SolverBase::CatalogInterface::CatalogType catalog;
  return catalog;
}

bool
SolverBase::CheckModelNames( array1d< string > & modelNames,
                             string const & attribute,
                             bool const allowEmpty ) const
{
  if( allowEmpty && modelNames.empty() )
  {
    return false;
  }

  // We can disable this if we want to be more strict
  if( modelNames.size() == 1 )
  {
    string const singleModelName = modelNames[0];
    modelNames.resizeDefault( m_targetRegionNames.size(), singleModelName );
  }

  GEOSX_ERROR_IF_NE_MSG( modelNames.size(),
                         m_targetRegionNames.size(),
                         getName() << ": invalid number of values in " << attribute
                                   << " attribute "
                                      "(expected one model name per target "
                                      "region, or one value for all regions)" );
  return true;
}

localIndex
SolverBase::targetRegionIndex( string const & regionName ) const
{
  auto const pos =
    std::find( m_targetRegionNames.begin(), m_targetRegionNames.end(), regionName );
  GEOSX_ERROR_IF(
    pos == m_targetRegionNames.end(),
    "Region " << regionName << " is not a target of solver " << getName() );
  return std::distance( m_targetRegionNames.begin(), pos );
}

real64
SolverBase::SolverStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                        real64 const & GEOSX_UNUSED_PARAM( dt ),
                        const integer GEOSX_UNUSED_PARAM( cycleNumber ),
                        DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  return 0;
}

void
SolverBase::Execute( real64 const time_n,
                     real64 const dt,
                     integer const cycleNumber,
                     integer const GEOSX_UNUSED_PARAM( eventCounter ),
                     real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                     Group * const domain )
{
  GEOSX_MARK_FUNCTION;
  real64 dtRemaining = dt;
  real64 nextDt = dt;

  integer const maxSubSteps = m_nonlinearSolverParameters.m_maxSubSteps;
  integer subStep = 0;

  for( ; subStep < maxSubSteps && dtRemaining > 0.0; ++subStep )
  {
    real64 const dtAccepted = SolverStep( time_n + ( dt - dtRemaining ),
                                          nextDt,
                                          cycleNumber,
                                          *domain->group_cast< DomainPartition * >() );
    /*
     * Let us check convergence history of previous solve:
     * - number of nonlinear iter.
     * - if the time-step was chopped. Then we can add some heuristics to choose next dt.
     * */
    dtRemaining -= dtAccepted;

    if( dtRemaining > 0.0 )
    {
      SetNextDt( dtAccepted, nextDt );
      nextDt = std::min( nextDt, dtRemaining );
    }

    if( getLogLevel() >= 1 && dtRemaining > 0.0 )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1,
                              getName() << ": sub-step = " << subStep
                                        << ", accepted dt = " << dtAccepted
                                        << ", remaining dt = " << dtRemaining );
    }
  }

  GEOSX_ERROR_IF( dtRemaining > 0.0,
                  "Maximum allowed number of sub-steps reached. Consider "
                  "increasing maxSubSteps." );

  // Decide what to do with the next Dt for the event running the solver.
  SetNextDt( nextDt, m_nextDt );
}

void
SolverBase::SetNextDt( real64 const & currentDt, real64 & nextDt )
{
  SetNextDtBasedOnNewtonIter( currentDt, nextDt );
}

void
SolverBase::SetNextDtBasedOnNewtonIter( real64 const & currentDt,
                                        real64 & nextDt )
{
  integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
  int const iterCutLimit = m_nonlinearSolverParameters.dtCutIterLimit();
  int const iterIncLimit = m_nonlinearSolverParameters.dtIncIterLimit();

  if( newtonIter < iterIncLimit )
  {
    // Easy convergence, let's double the time-step.
    nextDt = 2 * currentDt;
    GEOSX_LOG_LEVEL_RANK_0(
      1,
      getName() << ": Newton solver converged in less than " << iterIncLimit
                << " iterations, time-step required will be doubled." );
  }
  else if( newtonIter > iterCutLimit )
  {
    // Tough convergence let us make the time-step smaller!
    nextDt = currentDt / 2;
    GEOSX_LOG_LEVEL_RANK_0(
      1,
      getName() << ": Newton solver converged in more than " << iterCutLimit
                << " iterations, time-step required will be halved." );
  }
  else
  {
    nextDt = currentDt;
  }
}

real64
SolverBase::LinearImplicitStep( real64 const & time_n,
                                real64 const & dt,
                                integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                                DomainPartition & domain )
{
  // call setup for physics solver. Pre step allocations etc.
  // TODO: Nonlinear step does not call its own setup, need to decide on consistent behavior
  ImplicitStepSetup( time_n, dt, domain );

  // zero out matrix/rhs before assembly
  m_localMatrix.setValues< parallelDevicePolicy<> >( 0.0 );
  m_localRhs.setValues< parallelDevicePolicy<> >( 0.0 );

  // call assemble to fill the matrix and the rhs
  AssembleSystem( time_n,
                  dt,
                  domain,
                  m_dofManager,
                  m_localMatrix.toViewConstSizes(),
                  m_localRhs.toView() );

  // apply boundary conditions to system
  ApplyBoundaryConditions( time_n,
                           dt,
                           domain,
                           m_dofManager,
                           m_localMatrix.toViewConstSizes(),
                           m_localRhs.toView() );

  // Compose parallel LA matrix/rhs out of local LA matrix/rhs
  m_matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
  m_rhs.create( m_localRhs.toViewConst(), MPI_COMM_GEOSX );
  m_solution.createWithLocalSize( m_matrix.numLocalCols(), MPI_COMM_GEOSX );

  // Output the linear system matrix/rhs for debugging purposes
  DebugOutputSystem( 0.0, 0, 0, m_matrix, m_rhs );

  // Solve the linear system
  SolveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

  // Output the linear system solution for debugging purposes
  DebugOutputSolution( 0.0, 0, 0, m_solution );

  // Copy solution from parallel vector back to local
  // TODO: This step will not be needed when we teach LA vectors to wrap our pointers
  m_solution.extract( m_localSolution );

  // apply the system solution to the fields/variables
  ApplySystemSolution( m_dofManager, m_localSolution, 1.0, domain );

  // final step for completion of timestep. typically secondary variable updates and cleanup.
  ImplicitStepComplete( time_n, dt, domain );

  // return the achieved timestep
  return dt;
}

bool
SolverBase::LineSearch( real64 const & time_n,
                        real64 const & dt,
                        integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                        DomainPartition & domain,
                        DofManager const & dofManager,
                        CRSMatrixView< real64, globalIndex const > const & localMatrix,
                        arrayView1d< real64 > const & localRhs,
                        arrayView1d< real64 const > const & localSolution,
                        real64 const scaleFactor,
                        real64 & lastResidual )
{
  integer const maxNumberLineSearchCuts =
    m_nonlinearSolverParameters.m_lineSearchMaxCuts;
  real64 const lineSearchCutFactor =
    m_nonlinearSolverParameters.m_lineSearchCutFactor;

  // flag to determine if we should solve the system and apply the solution. If the line
  // search fails we just bail.
  bool lineSearchSuccess = false;

  real64 residualNorm = lastResidual;

  // scale factor is value applied to the previous solution. In this case we want to
  // subtract a portion of the previous solution.
  real64 localScaleFactor = -scaleFactor;
  real64 cumulativeScale = scaleFactor;

  // main loop for the line search.
  for( integer lineSearchIteration = 0;
       lineSearchIteration < maxNumberLineSearchCuts;
       ++lineSearchIteration )
  {
    // cut the scale factor by half. This means that the scale factors will
    // have values of -0.5, -0.25, -0.125, ...
    localScaleFactor *= lineSearchCutFactor;
    cumulativeScale += localScaleFactor;

    if( !CheckSystemSolution( domain, dofManager, localSolution, localScaleFactor ) )
    {
      GEOSX_LOG_LEVEL_RANK_0( 1,
                              "        Line search "
                                << lineSearchIteration
                                << ", solution check failed" );
      continue;
    }

    ApplySystemSolution( dofManager, localSolution, localScaleFactor, domain );

    // re-assemble system
    localMatrix.setValues< parallelDevicePolicy<> >( 0.0 );
    localRhs.setValues< parallelDevicePolicy<> >( 0.0 );
    AssembleSystem( time_n, dt, domain, dofManager, localMatrix, localRhs );

    // apply boundary conditions to system
    ApplyBoundaryConditions( time_n, dt, domain, dofManager, localMatrix, localRhs );

    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      char output[100];
      sprintf( output, "        Line search @ %0.3f:      ", cumulativeScale );
      std::cout << output;
    }

    // get residual norm
    residualNorm = CalculateResidualNorm( domain, dofManager, localRhs );

    if( getLogLevel() >= 1 && logger::internal::rank == 0 )
    {
      std::cout << std::endl;
    }

    // if the residual norm is less than the last residual, we can proceed to the
    // solution step
    if( residualNorm < lastResidual )
    {
      lineSearchSuccess = true;
      break;
    }
  }

  lastResidual = residualNorm;
  return lineSearchSuccess;
}

/**
 * @brief Eisenstat-Walker adaptive tolerance
 *
 * This method enables an inexact-Newton method is which the linear solver
 * tolerance is chosen based on the nonlinear solver convergence behavior.
 * In early Newton iterations, the search direction is usually imprecise, and
 * therefore a weak linear convergence tolerance can be chosen to minimize
 * computational cost.  As the search gets closer to the true solution, however,
 * more stringent linear tolerances are necessary to maintain quadratic convergence
 * behavior.
 *
 * The user can set the weakest tolerance allowed, with a default of 1e-3.
 * Even weaker values (e.g. 1e-2,1e-1) can be used for further speedup, but may
 * occasionally cause convergence problems.  Use this parameter with caution.  The
 * most stringent tolerance is hardcoded to 1e-8, which is sufficient for
 * most problems.
 *
 * See Eisenstat, S.C. and Walker, H.F., 1996. Choosing the forcing terms in an
 * inexact Newton method. SIAM Journal on Scientific Computing, 17(1), pp.16-32.
 *
 * @param newNewtonNorm Residual norm at current iteration
 * @param oldNewtonNorm Residual norm at previous iteration
 * @param weakestTol Weakest tolerance allowed (default 1e-3).
 * @return Adaptive tolerance recommendation
 */
real64
SolverBase::EisenstatWalker( real64 const newNewtonNorm,
                             real64 const oldNewtonNorm,
                             real64 const weakestTol )
{
  real64 const strongestTol = 1e-8;
  real64 const exponent = 2.0;
  real64 const gamma = 0.9;

  real64 normRatio = newNewtonNorm / oldNewtonNorm;
  if( normRatio > 1 ) normRatio = 1;

  real64 newKrylovTol = gamma * std::pow( normRatio, exponent );
  real64 altKrylovTol = gamma * std::pow( oldNewtonNorm, exponent );

  real64 krylovTol = std::max( newKrylovTol, altKrylovTol );
  krylovTol = std::min( krylovTol, weakestTol );
  krylovTol = std::max( krylovTol, strongestTol );

  return krylovTol;
}

real64
SolverBase::NonlinearImplicitStep( real64 const & time_n,
                                   real64 const & dt,
                                   integer const cycleNumber,
                                   DomainPartition & domain )
{
  GEOSX_MARK_FUNCTION;
  // dt may be cut during the course of this step, so we are keeping a local
  // value to track the achieved dt for this step.
  real64 stepDt = dt;

  integer const maxNewtonIter = m_nonlinearSolverParameters.m_maxIterNewton;
  integer const minNewtonIter = m_nonlinearSolverParameters.m_minIterNewton;
  real64 const newtonTol = m_nonlinearSolverParameters.m_newtonTol;

  integer const maxNumberDtCuts = m_nonlinearSolverParameters.m_maxTimeStepCuts;
  real64 const dtCutFactor = m_nonlinearSolverParameters.m_timeStepCutFactor;

  bool const allowNonConverged =
    m_nonlinearSolverParameters.m_allowNonConverged > 0;

  integer & dtAttempt = m_nonlinearSolverParameters.m_numdtAttempts;

  // a flag to denote whether we have converged
  integer isConverged = 0;

  // outer loop attempts to apply full timestep, and managed the cutting of the timestep if
  // required.
  for( dtAttempt = 0; dtAttempt < maxNumberDtCuts; ++dtAttempt )
  {
    // reset the solver state, since we are restarting the time step
    if( dtAttempt > 0 )
    {
      ResetStateToBeginningOfStep( domain );
    }

    // keep residual from previous iteration in case we need to do a line search
    real64 lastResidual = 1e99;
    integer & newtonIter = m_nonlinearSolverParameters.m_numNewtonIterations;
    real64 scaleFactor = 1.0;

    // main Newton loop
    for( newtonIter = 0; newtonIter < maxNewtonIter; ++newtonIter )
    {
      if( getLogLevel() >= 1 && logger::internal::rank == 0 )
      {
        char output[100] = { 0 };
        sprintf( output, "    Attempt: %2d, NewtonIter: %2d", dtAttempt, newtonIter );
        std::cout << output << std::endl;
      }

      // zero out matrix/rhs before assembly
      m_localMatrix.setValues< parallelDevicePolicy<> >( 0.0 );
      m_localRhs.setValues< parallelDevicePolicy<> >( 0.0 );

      // call assemble to fill the matrix and the rhs
      AssembleSystem( time_n,
                      stepDt,
                      domain,
                      m_dofManager,
                      m_localMatrix.toViewConstSizes(),
                      m_localRhs.toView() );

      // apply boundary conditions to system
      ApplyBoundaryConditions( time_n,
                               stepDt,
                               domain,
                               m_dofManager,
                               m_localMatrix.toViewConstSizes(),
                               m_localRhs.toView() );

      // TODO: maybe add scale function here?
      // Scale()

      // get residual norm
      real64 residualNorm =
        CalculateResidualNorm( domain, m_dofManager, m_localRhs.toViewConst() );

      if( getLogLevel() >= 1 && logger::internal::rank == 0 )
      {
        {
          char output[100] = { 0 };
          sprintf( output, "    ( R ) = ( %4.2e ) ; ", residualNorm );
          std::cout << output << std::endl;
        }

        if( newtonIter > 0 )
        {
          char output[100] = { 0 };
          sprintf( output,
                   "    Last LinSolve(iter,res) = ( %3d, %4.2e ) ; ",
                   m_linearSolverResult.numIterations,
                   m_linearSolverResult.residualReduction );
          std::cout << output << std::endl;
        }
      }

      // if the residual norm is less than the Newton tolerance we denote that we have
      // converged and break from the Newton loop immediately.

      if( residualNorm < newtonTol && newtonIter >= minNewtonIter )
      {
        isConverged = 1;
        break;
      }

      // do line search in case residual has increased
      if( m_nonlinearSolverParameters.m_lineSearchAction > 0 &&
          residualNorm > lastResidual )
      {
        residualNorm = lastResidual;
        bool lineSearchSuccess = LineSearch( time_n,
                                             stepDt,
                                             cycleNumber,
                                             domain,
                                             m_dofManager,
                                             m_localMatrix.toViewConstSizes(),
                                             m_localRhs.toView(),
                                             m_localSolution.toViewConst(),
                                             scaleFactor,
                                             residualNorm );

        if( !lineSearchSuccess )
        {
          if( m_nonlinearSolverParameters.m_lineSearchAction == 1 )
          {
            GEOSX_LOG_LEVEL_RANK_0( 1,
                                    "        Line search failed to produce "
                                    "reduced residual. Accepting iteration." );
          }
          else if( m_nonlinearSolverParameters.m_lineSearchAction == 2 )
          {
            // if line search failed, then break out of the main Newton loop. Timestep will be cut.
            GEOSX_LOG_LEVEL_RANK_0( 1,
                                    "        Line search failed to produce "
                                    "reduced residual. Exiting Newton Loop." );
            break;
          }
        }
      }

      // if using adaptive Krylov tolerance scheme, update tolerance.
      LinearSolverParameters::Krylov & krylovParams =
        m_linearSolverParameters.get().krylov;
      if( krylovParams.useAdaptiveTol )
      {
        krylovParams.relTolerance =
          EisenstatWalker( residualNorm, lastResidual, krylovParams.weakestTol );
      }

      // Compose parallel LA matrix/rhs out of local LA matrix/rhs
      m_matrix.create( m_localMatrix.toViewConst(), MPI_COMM_GEOSX );
      m_rhs.create( m_localRhs.toViewConst(), MPI_COMM_GEOSX );
      m_solution.createWithLocalSize( m_matrix.numLocalCols(), MPI_COMM_GEOSX );

      // Output the linear system matrix/rhs for debugging purposes
      DebugOutputSystem( time_n, cycleNumber, newtonIter, m_matrix, m_rhs );

      // Solve the linear system
      SolveSystem( m_dofManager, m_matrix, m_rhs, m_solution );

      // Output the linear system solution for debugging purposes
      DebugOutputSolution( time_n, cycleNumber, newtonIter, m_solution );

      // Copy solution from parallel vector back to local
      // TODO: This step will not be needed when we teach LA vectors to wrap our pointers
      m_solution.extract( m_localSolution );

      scaleFactor =
        ScalingForSystemSolution( domain, m_dofManager, m_localSolution );

      if( !CheckSystemSolution( domain, m_dofManager, m_localSolution, scaleFactor ) )
      {
        // TODO try chopping (similar to line search)
        GEOSX_LOG_RANK_0( "    Solution check failed. Newton loop terminated." );
        break;
      }

      // apply the system solution to the fields/variables
      ApplySystemSolution( m_dofManager, m_localSolution, scaleFactor, domain );

      lastResidual = residualNorm;
    }

    if( isConverged )
    {
      break;  // out of outer loop
    }
    else
    {
      // cut timestep, go back to beginning of step and restart the Newton loop
      stepDt *= dtCutFactor;
      GEOSX_LOG_LEVEL_RANK_0( 1, "New dt = " << stepDt );
    }
  }

  if( !isConverged )
  {
    GEOSX_LOG_RANK_0( "Convergence not achieved." );

    if( allowNonConverged )
    {
      GEOSX_LOG_RANK_0( "The accepted solution may be inaccurate." );
    }
    else
    {
      GEOSX_ERROR( "Nonconverged solutions not allowed. Terminating..." );
    }
  }

  // return the achieved timestep
  return stepDt;
}

real64
SolverBase::ExplicitStep( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                          real64 const & GEOSX_UNUSED_PARAM( dt ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "SolverBase::ExplicitStep called!. Should be overridden." );
  return 0;
}

void
SolverBase::ImplicitStepSetup( real64 const & GEOSX_UNUSED_PARAM( time_n ),
                               real64 const & GEOSX_UNUSED_PARAM( dt ),
                               DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "SolverBase::ImplicitStepSetup called!. Should be overridden." );
}

void
SolverBase::SetupDofs( DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
                       DofManager & GEOSX_UNUSED_PARAM( dofManager ) ) const
{
  GEOSX_ERROR( "SolverBase::SetupDofs called!. Should be overridden." );
}

void
SolverBase::SetupSystem( DomainPartition & domain,
                         DofManager & dofManager,
                         CRSMatrix< real64, globalIndex > & localMatrix,
                         array1d< real64 > & localRhs,
                         array1d< real64 > & localSolution,
                         bool const setSparsity )
{
  GEOSX_MARK_FUNCTION;

  dofManager.setMesh( domain, 0, 0 );

  SetupDofs( domain, dofManager );
  dofManager.reorderByRank();

  localIndex const numLocalRows = dofManager.numLocalDofs();

  SparsityPattern< globalIndex > pattern;
  if( setSparsity )
  {
    dofManager.setSparsityPattern( pattern );
    localMatrix.assimilate< parallelDevicePolicy<> >( std::move( pattern ) );
  }

  localRhs.resize( numLocalRows );
  localSolution.resize( numLocalRows );

  localMatrix.setName( this->getName() + "/localMatrix" );
  localRhs.setName( this->getName() + "/localRhs" );
  localSolution.setName( this->getName() + "/localSolution" );
}

void
SolverBase::AssembleSystem(
  real64 const GEOSX_UNUSED_PARAM( time ),
  real64 const GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR( "SolverBase::Assemble called!. Should be overridden." );
}

void
SolverBase::ApplyBoundaryConditions(
  real64 const GEOSX_UNUSED_PARAM( time ),
  real64 const GEOSX_UNUSED_PARAM( dt ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  CRSMatrixView< real64, globalIndex const > const & GEOSX_UNUSED_PARAM( localMatrix ),
  arrayView1d< real64 > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR(
    "SolverBase::ApplyBoundaryConditions called!. Should be overridden." );
}

namespace
{
/**
 * @brief Helper for debug output of linear algebra objects (matrices and vectors)
 * @tparam T type of LA object (must have stream insertion and .write() implemented)
 * @param obj                the object to output
 * @param cycleNumber        event cycle number
 * @param nonlinearIteration nonlinear iteration number
 * @param filePrefix          short filename prefix (e.g. "mat")
 * @param screenName           long name for screen output (e.g. "System matrix")
 * @param toScreen           whether to print on screen
 * @param toFile             whether to write to file
 */
template< typename T >
void
debugOutputLAObject( T const & obj,
                     real64 const & GEOSX_UNUSED_PARAM( time ),
                     integer const cycleNumber,
                     integer const nonlinearIteration,
                     string const & filePrefix,
                     string const & screenName,
                     bool const toScreen,
                     bool const toFile )
{
  if( toScreen )
  {
    string const frame( screenName.size() + 1, '=' );
    GEOSX_LOG_RANK_0( frame << "\n"
                            << screenName << ":\n"
                            << frame );
    GEOSX_LOG( obj );
  }

  if( toFile )
  {
    char filename[200] = { 0 };
    snprintf( filename,
              200,
              "%s_%06d_%02d.mtx",
              filePrefix.c_str(),
              cycleNumber,
              nonlinearIteration );
    obj.write( filename, LAIOutputFormat::MATRIX_MARKET );
    GEOSX_LOG_RANK_0( screenName << " written to " << filename );
  }
}

}  // namespace

void
SolverBase::DebugOutputSystem( real64 const & time,
                               integer const cycleNumber,
                               integer const nonlinearIteration,
                               ParallelMatrix const & matrix,
                               ParallelVector const & rhs ) const
{
  debugOutputLAObject( matrix,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_mat",
                       "System matrix",
                       getLogLevel() == 2,
                       getLogLevel() >= 3 );

  debugOutputLAObject( rhs,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_rhs",
                       "System right-hand side",
                       getLogLevel() == 2,
                       getLogLevel() >= 3 );
}

void
SolverBase::DebugOutputSolution( real64 const & time,
                                 integer const cycleNumber,
                                 integer const nonlinearIteration,
                                 ParallelVector const & solution ) const
{
  debugOutputLAObject( solution,
                       time,
                       cycleNumber,
                       nonlinearIteration,
                       getName() + "_sol",
                       "System solution",
                       getLogLevel() == 2,
                       getLogLevel() >= 3 );
}

real64
SolverBase::CalculateResidualNorm(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localRhs ) )
{
  GEOSX_ERROR(
    "SolverBase::CalculateResidualNorm called!. Should be overridden." );
  return 0;
}

void
SolverBase::SolveSystem( DofManager const & dofManager,
                         ParallelMatrix & matrix,
                         ParallelVector & rhs,
                         ParallelVector & solution )
{
  GEOSX_MARK_FUNCTION;

  //  Keep for debugging comparisons
  //  static int count = 0;
  //  if( count < 2 )
  //  {
  //  std::cout<<"************************* MATRIX *************************"<<std::endl;
  //  std::cout<<matrix<<std::endl<<std::endl;
  //  std::cout<<"************************* RHS *************************"<<std::endl;
  //  std::cout<<rhs<<std::endl<<std::endl;
  //  }

  LinearSolverParameters const & params = m_linearSolverParameters.get();

  // TODO: We probably want to keep an instance of linear solver as a member of physics solver
  //       so we can have constant access to last solve statistics, convergence history, etc.
  //       This requires unifying "LAI interface" solvers with "native" Krylov solvers somehow.

  if( params.solverType == "direct" || !m_precond )
  {
    LinearSolver solver( params );
    solver.solve( matrix, solution, rhs, &dofManager );
    m_linearSolverResult = solver.result();
  }
  else
  {
    m_precond->compute( matrix, dofManager );
    std::unique_ptr< KrylovSolver< ParallelVector > > solver =
      KrylovSolver< ParallelVector >::Create( params, matrix, *m_precond );
    solver->solve( rhs, solution );
    m_linearSolverResult = solver->result();
  }

  //  Keep for debugging comparisons
  //  if( count < 2 )
  //  {
  //  std::cout<<"************************* SOLUTION *************************"<<std::endl;
  //  std::cout<<solution<<std::endl<<std::endl;
  //  }
  //  ++count;

  GEOSX_WARNING_IF( !m_linearSolverResult.success(), "Linear solution failed" );
}

bool
SolverBase::CheckSystemSolution(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localSolution ),
  real64 const GEOSX_UNUSED_PARAM( scalingFactor ) )
{
  return true;
}

real64
SolverBase::ScalingForSystemSolution(
  DomainPartition const & GEOSX_UNUSED_PARAM( domain ),
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localSolution ) )
{
  return 1.0;
}

void
SolverBase::ApplySystemSolution(
  DofManager const & GEOSX_UNUSED_PARAM( dofManager ),
  arrayView1d< real64 const > const & GEOSX_UNUSED_PARAM( localSolution ),
  real64 const GEOSX_UNUSED_PARAM( scalingFactor ),
  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR( "SolverBase::ApplySystemSolution called!. Should be overridden." );
}

void
SolverBase::ResetStateToBeginningOfStep(
  DomainPartition & GEOSX_UNUSED_PARAM( const ) )
{
  GEOSX_ERROR(
    "SolverBase::ResetStateToBeginningOfStep called!. Should be overridden." );
}

void
SolverBase::ImplicitStepComplete( real64 const & GEOSX_UNUSED_PARAM( time ),
                                  real64 const & GEOSX_UNUSED_PARAM( dt ),
                                  DomainPartition & GEOSX_UNUSED_PARAM( domain ) )
{
  GEOSX_ERROR(
    "SolverBase::ImplicitStepComplete called!. Should be overridden." );
}

R1Tensor const
SolverBase::gravityVector() const
{
  R1Tensor rval;
  if( getParent()->group_cast< PhysicsSolverManager const * >() != nullptr )
  {
    rval = getParent()->getReference< R1Tensor >(
      PhysicsSolverManager::viewKeyStruct::gravityVectorString );
  }
  else
  {
    rval = { 0.0, 0.0, -9.81 };
  }
  return rval;
}

}  // namespace geosx
