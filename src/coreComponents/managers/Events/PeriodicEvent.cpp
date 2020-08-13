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
 * @file PeriodicEvent.cpp
 */

#include "PeriodicEvent.hpp"
#include "managers/Functions/FunctionManager.hpp"

namespace geosx
{
using namespace dataRepository;

PeriodicEvent::PeriodicEvent( const std::string & name, Group * const parent ) :
  EventBase( name, parent ),
  m_functionTarget( nullptr ),
  m_timeFrequency( -1.0 ),
  m_cycleFrequency( 1 ),
  m_targetExactTimestep( 0 ),
  m_functionName(),
  m_functionInputObject(),
  m_functionInputSetname(),
  m_functionStatOption( 0 ),
  m_eventThreshold( 0.0 )
{
  registerWrapper( viewKeyStruct::timeFrequencyString, &m_timeFrequency )
    ->setApplyDefaultValue( -1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Event application frequency (time).  Note: if this value is specified, "
      "it will override any cycle-based behavior." );

  registerWrapper( viewKeyStruct::cycleFrequencyString, &m_cycleFrequency )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Event application frequency (cycle, default)" );

  registerWrapper( viewKeyStruct::targetExactTimestepString, &m_targetExactTimestep )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If this option is set, the event will reduce its timestep requests to "
      "match the specified timeFrequency perfectly: dt_request = "
      "min(dt_request, t_last + time_frequency - time))." );

  registerWrapper( viewKeyStruct::functionNameString, &m_functionName )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Name of an optional function to evaluate when the time/cycle criteria "
      "are met."
      "If the result is greater than the specified eventThreshold, the "
      "function will continue to execute." );

  registerWrapper( viewKeyStruct::functionInputObjectString, &m_functionInputObject )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If the optional function requires an object as an input, specify its "
      "path here." );

  registerWrapper( viewKeyStruct::functionInputSetnameString,
                   &m_functionInputSetname )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If the optional function is applied to an object, specify the setname "
      "to evaluate (default = everything)." );

  registerWrapper( viewKeyStruct::functionStatOptionString, &m_functionStatOption )
    ->setApplyDefaultValue( 0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If the optional function is applied to an object, specify the statistic "
      "to compare to the eventThreshold."
      "The current options include: min, avg, and max." );

  registerWrapper( viewKeyStruct::eventThresholdString, &m_eventThreshold )
    ->setApplyDefaultValue( 0.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If the optional function is used, the event will execute if the value "
      "returned by the function exceeds this threshold." );
}

PeriodicEvent::~PeriodicEvent()
{}

void
PeriodicEvent::EstimateEventTiming( real64 const time,
                                    real64 const dt,
                                    integer const cycle,
                                    Group * domain )
{
  // Check event status
  if( cycle == 0 )
  {
    setReadyForExec();
  }
  else if( m_timeFrequency >= 0.0 )
  {
    if( dt <= 0 )
    {
      setIdle();
    }
    else
    {
      // Note: add a small value to this forecast to account for floating point errors
      real64 forecast = ( ( m_timeFrequency - ( time - m_lastTime ) ) / dt ) + 1e-10;
      setForecast( static_cast< integer >( std::min( std::max( forecast, 0.0 ), 1e9 ) ) );
    }
  }
  else
  {
    setForecast( m_cycleFrequency - ( cycle - m_lastCycle ) );
  }

  if( this->isReadyForExec() && ( !m_functionName.empty() ) )
  {
    CheckOptionalFunctionThreshold( time, dt, cycle, domain );
  }
}

void
PeriodicEvent::CheckOptionalFunctionThreshold(
  real64 const time,
  real64 const GEOSX_UNUSED_PARAM( dt ),
  integer const GEOSX_UNUSED_PARAM( cycle ),
  Group * GEOSX_UNUSED_PARAM( domain ) )
{
  // Grab the function
  FunctionManager & functionManager = FunctionManager::Instance();
  FunctionBase * function = functionManager.GetGroup< FunctionBase >( m_functionName );

  real64 result = 0.0;
  if( m_functionInputObject.empty() )
  {
    // This is a time-only function
    result = function->Evaluate( &time );
  }
  else
  {
    // Link the target object
    if( m_functionTarget == nullptr )
    {
      m_functionTarget = this->GetGroupByPath( m_functionInputObject );
    }

    // Get the set
    SortedArray< localIndex > mySet;
    if( m_functionInputSetname.empty() )
    {
      for( localIndex ii = 0; ii < m_functionTarget->size(); ++ii )
      {
        mySet.insert( ii );
      }
    }
    else
    {
      dataRepository::Group const * sets =
        m_functionTarget->GetGroup( periodicEventViewKeys.functionSetNames );
      SortedArrayView< localIndex const > const & functionSet =
        sets->getReference< SortedArray< localIndex > >( m_functionInputSetname );

      for( localIndex const index : functionSet )
      {
        mySet.insert( index );
      }
    }

    // Find the function (min, average, max)
    real64_array stats = function->EvaluateStats( m_functionTarget, time, mySet );
    result = stats[m_functionStatOption];

    // Because the function applied to an object may differ by rank, synchronize
    // (Note: this shouldn't occur very often, since it is only called if the base forecast <= 0)
#ifdef GEOSX_USE_MPI
    real64 result_global;
    MPI_Allreduce( &result, &result_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
    result = result_global;
#endif
  }

  // Forcast event
  if( result > m_eventThreshold )
  {
    setReadyForExec();
  }
  else
  {
    setIdle();
  }
}

real64
PeriodicEvent::GetEventTypeDtRequest( real64 const time )
{
  real64 requestedDt = std::numeric_limits< real64 >::max();

  if( ( m_timeFrequency > 0 ) && ( m_targetExactTimestep > 0 ) )
  {
    // Limit the timestep request to match the exact execution frequency
    real64 nextTargetTime = m_lastTime + m_timeFrequency;
    real64 tmp_t = std::nextafter( time, time + 1.0 );

    if( tmp_t < nextTargetTime )
    {
      requestedDt = std::min( requestedDt, nextTargetTime - time );
    }
    else
    {
      // Note: This should only occur on a cycle where the event is executing
      requestedDt = std::min( requestedDt, m_timeFrequency );
    }
  }

  return requestedDt;
}

void
PeriodicEvent::Cleanup( real64 const time_n,
                        integer const cycleNumber,
                        integer const GEOSX_UNUSED_PARAM( eventCounter ),
                        real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                        Group * domain )
{
  // Only call the cleanup method of the target/children if it is within its application time
  if( isActive( time_n ) )
  {
    ExecutableGroup * target = GetEventTarget();
    if( target != nullptr )
    {
      // Cleanup the target
      target->Cleanup( time_n, cycleNumber, 0, 0, domain );
    }

    // Cleanup any sub-events
    this->forSubGroups< EventBase >( [&]( EventBase & subEvent ) {
      subEvent.Cleanup( time_n, cycleNumber, 0, 0, domain );
    } );
  }
}

REGISTER_CATALOG_ENTRY( EventBase, PeriodicEvent, std::string const &, Group * const )

} /* namespace geosx */
