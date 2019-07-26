/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

 /**
  * @file PeriodicEvent.cpp
  */

#include "PeriodicEvent.hpp"
#include "managers/Functions/NewFunctionManager.hpp"

namespace geosx
{

using namespace dataRepository;


PeriodicEvent::PeriodicEvent( const std::string& name,
                              ManagedGroup * const parent ):
  EventBase(name,parent),
  m_functionTarget(nullptr),
  m_timeFrequency(-1.0),
  m_cycleFrequency(1),
  m_targetExactTimestep(0),
  m_functionName(),
  m_functionInputObject(),
  m_functionInputSetname(),
  m_functionStatOption(0),
  m_eventThreshold(0.0)
{
  RegisterViewWrapper(viewKeyStruct::timeFrequencyString, &m_timeFrequency, false )->
    setApplyDefaultValue(-1.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("event frequency (time)");

  RegisterViewWrapper(viewKeyStruct::cycleFrequencyString, &m_cycleFrequency, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("event frequency (cycle, Default)");

  RegisterViewWrapper(viewKeyStruct::targetExactTimestepString, &m_targetExactTimestep, false )->
    setApplyDefaultValue(1)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("allows timesteps to be truncated to match time frequency perfectly");

  RegisterViewWrapper(viewKeyStruct::functionNameString, &m_functionName, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Name of the symbolic math function");

  RegisterViewWrapper(viewKeyStruct::functionInputObjectString, &m_functionInputObject, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Path of the function input object (directory format)");

  RegisterViewWrapper(viewKeyStruct::functionInputSetnameString, &m_functionInputSetname, false )->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Setname of the input object (if empty, default to everything)");

  RegisterViewWrapper(viewKeyStruct::functionStatOptionString, &m_functionStatOption, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Selection of the min/avg/max for functions that target vectors");

  RegisterViewWrapper(viewKeyStruct::eventThresholdString, &m_eventThreshold, false )->
    setApplyDefaultValue(0.0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("event threshold");

}


PeriodicEvent::~PeriodicEvent()
{}






void PeriodicEvent::EstimateEventTiming(real64 const time,
                                        real64 const dt, 
                                        integer const cycle,
                                        ManagedGroup * domain)
{
  // Check event status
  if (cycle == 0)
  {
    SetForecast(0);
  } 
  else if (m_timeFrequency >= 0.0)
  {
    if (dt <= 0)
    {
      SetForecast(std::numeric_limits<integer>::max());
    } 
    else
    {
      // Note: add a small value to this forecast to account for floating point errors
      real64 forecast = ((m_timeFrequency - (time - m_lastTime)) / dt) + 1e-10;
      SetForecast(static_cast<integer>(std::min(std::max(forecast, 0.0), 1e9)));
    }
  }
  else
  {
    SetForecast(m_cycleFrequency - (cycle - m_lastCycle));
  }

  if ((this->GetForecast() <= 0) && (m_functionName.empty() == 0))
  {
    CheckOptionalFunctionThreshold(time, dt, cycle, domain);
  }
}


void PeriodicEvent::CheckOptionalFunctionThreshold(real64 const time,
                                                   real64 const dt, 
                                                   integer const cycle,
                                                   ManagedGroup * domain)
{
  // Grab the function
  NewFunctionManager * functionManager = NewFunctionManager::Instance();
  FunctionBase * function = functionManager->GetGroup<FunctionBase>(m_functionName);

  real64 result = 0.0;
  if (m_functionInputObject.empty())
  {
    // This is a time-only function
    result = function->Evaluate(&time);
  }
  else
  {
    // Link the target object
    if (m_functionTarget == nullptr)
    {
      m_functionTarget = this->GetGroupByPath(m_functionInputObject);
    }

    // Get the set
    set<localIndex> mySet;
    if (m_functionInputSetname.empty())
    {
      for(localIndex ii=0; ii<m_functionTarget->size(); ++ii)
      {
        mySet.insert(ii);
      }
    }
    else
    {
      dataRepository::ManagedGroup const * sets = m_functionTarget->GetGroup(periodicEventViewKeys.functionSetNames);
      mySet = sets->getReference< set<localIndex> >(m_functionInputSetname);
    }

    // Find the function (min, average, max)
    real64_array stats = function->EvaluateStats(m_functionTarget, time, mySet);
    result = stats[m_functionStatOption];

    // Because the function applied to an object may differ by rank, synchronize
    // (Note: this shouldn't occur very often, since it is only called if the base forecast <= 0)
#ifdef GEOSX_USE_MPI
      real64 result_global;
      MPI_Allreduce(&result, &result_global, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_GEOSX);
      result = result_global;
#endif
  }
  
  // Forcast event
  if (result > m_eventThreshold)
  {
    SetForecast(0);
  }
  else
  {
    SetForecast(std::numeric_limits<integer>::max());
  }
}


real64 PeriodicEvent::GetEventTypeDtRequest(real64 const time)
{
  real64 requestedDt = std::numeric_limits<real64>::max();

  if ((m_timeFrequency > 0) && (m_targetExactTimestep > 0))
  {
    // Limit the timestep request to match the exact execution frequency
    real64 nextTargetTime = m_lastTime + m_timeFrequency;
    real64 tmp_t = std::nextafter(time, time + 1.0);

    if (tmp_t < nextTargetTime)
    {
      requestedDt = std::min(requestedDt, nextTargetTime - time);
    }
    else
    {
      // Note: This should only occur on a cycle where the event is executing
      requestedDt = std::min(requestedDt, m_timeFrequency);
    }
  }

  return requestedDt;
}



REGISTER_CATALOG_ENTRY( EventBase, PeriodicEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
