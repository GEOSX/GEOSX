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
 * @file EventManager.cpp
 */

#include "EventManager.hpp"
#include "managers/Events/EventBase.hpp"
#include "common/TimingMacros.hpp"
#include "MPI_Communications/CommunicationTools.hpp"

#if HAVE_TRIBOLCOUPLING
#include "coupling/TribolCoupling.hpp"
#endif

namespace geosx
{

using namespace dataRepository;
using namespace cxx_utilities;


EventManager::EventManager( std::string const & name,
                            ManagedGroup * const parent ):
  ManagedGroup( name, parent),
  m_maxTime(),
  m_maxCycle(),
  m_verbosity(),
  m_time(),
  m_dt(),
  m_cycle(),
  m_currentSubEvent()
{
  setInputFlags(InputFlags::REQUIRED);
  
  RegisterViewWrapper(viewKeyStruct::maxTimeString, &m_maxTime, false )->
    setApplyDefaultValue(std::numeric_limits<real64>::max())->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum simulation time.");

  RegisterViewWrapper(viewKeyStruct::maxCycleString, &m_maxCycle, false )->
    setApplyDefaultValue(std::numeric_limits<integer>::max())->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Maximum simulation cycle.");

  RegisterViewWrapper(viewKeyStruct::verbosityString, &m_verbosity, false )->
    setApplyDefaultValue(0)->
    setInputFlag(InputFlags::OPTIONAL)->
    setDescription("Verbosity level");

  RegisterViewWrapper(viewKeyStruct::timeString, &m_time, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation time.");

  RegisterViewWrapper(viewKeyStruct::dtString, &m_dt, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation timestep.");

  RegisterViewWrapper(viewKeyStruct::cycleString, &m_cycle, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("Current simulation cycle number.");

  RegisterViewWrapper(viewKeyStruct::currentSubEventString, &m_currentSubEvent, false )->
    setRestartFlags(RestartFlags::WRITE_AND_READ)->
    setDescription("index of the current subevent.");

}


EventManager::~EventManager()
{}



ManagedGroup * EventManager::CreateChild( string const & childKey, string const & childName )
{
  GEOS_LOG_RANK_0("Adding Event: " << childKey << ", " << childName);
  std::unique_ptr<EventBase> event = EventBase::CatalogInterface::Factory( childKey, childName, this );
  return this->RegisterGroup<EventBase>( childName, std::move(event) );
}


void EventManager::ExpandObjectCatalogs()
{
  // During schema generation, register one of each type derived from EventBase here
  for (auto& catalogIter: EventBase::GetCatalog())
  {
    CreateChild( catalogIter.first, catalogIter.first );
  }
}


void EventManager::Run(dataRepository::ManagedGroup * domain)
{
  GEOSX_MARK_FUNCTION;

  integer exitFlag = 0;

  // Setup event targets, sequence indicators
  array1d<integer> eventCounters(2);
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->GetTargetReferences();
    subEvent->GetExecutionOrder(eventCounters);
  });

  // Set the progress indicators
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->SetProgressIndicator(eventCounters);
  });

  // Inform user if it appears this is a mid-loop restart
  if ((m_currentSubEvent > 0))
  {
    GEOS_LOG_RANK_0("The restart-file was written during step " << m_currentSubEvent << " of the event loop.  Resuming from that point.");
  }
#if HAVE_TRIBOLCOUPLING
  TribolCoupling::Initialize(this, domain) ;
#endif

  // Run problem
  // Note: if currentSubEvent > 0, then we are resuming from a restart file
  while(1) {
    int terminate = ((m_time < m_maxTime) && (m_cycle < m_maxCycle) && (exitFlag == 0)) ? 0 : 1 ;

    if (cycle > 0) {
#if HAVE_TRIBOLCOUPLING
       TribolCoupling::SyncTermination(&terminate) ;
#endif

       if (terminate) {
          break ;
       }

#if HAVE_TRIBOLCOUPLING
       real64 newDt ;
       TribolCoupling::SyncTimestep(&newDt) ;

       if (newDt < dt)
       {
           std::cout << "     dt: " << dt << ", coupled dt=" << newDt << std::endl;
           GEOS_ERROR( "TRIBOL coupling error" );
       }
#endif
    } else if (terminate) {
      break ;
    }
    // Determine the cycle timestep
    if (m_currentSubEvent == 0)
    {
      // The max dt request
      m_dt = m_maxTime - m_time;

      // Determine the dt requests for each event
      for ( ; m_currentSubEvent<this->numSubGroups(); ++m_currentSubEvent)
      {
        EventBase * subEvent = static_cast<EventBase *>( this->GetSubGroups()[m_currentSubEvent] );
        m_dt = std::min(subEvent->GetTimestepRequest(m_time), m_dt);
      }
      m_currentSubEvent = 0;

#ifdef GEOSX_USE_MPI
      // Find the min dt across procfesses
      GEOSX_MARK_BEGIN("EventManager::MPI calls");
      real64 dt_global;
      MPI_Allreduce(&m_dt, &dt_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_GEOSX);
      m_dt = dt_global;

      GEOSX_MARK_END("EventManager::MPI calls");
#endif
    }

    GEOS_LOG_RANK_0("Time: " << m_time << "s, dt:" << m_dt << "s, Cycle: " << m_cycle);

    // Execute 
    for ( ; m_currentSubEvent<this->numSubGroups(); ++m_currentSubEvent)
    {
      EventBase * subEvent = static_cast<EventBase *>( this->GetSubGroups()[m_currentSubEvent] );

      // Calculate the event and sub-event forecasts
      subEvent->CheckEvents(m_time, m_dt, m_cycle, domain);
      integer eventForecast = subEvent->GetForecast();

      if (m_verbosity > 0)
      {
        GEOS_LOG_RANK_0("     Event: " << m_currentSubEvent << " (" << subEvent->getName() << "), dt_request=" << subEvent->GetCurrentEventDtRequest() << ", forecast=" << eventForecast);
      }

      // Execute, signal events
      if (eventForecast == 1)
      {
        subEvent->SignalToPrepareForExecution(m_time, m_dt, m_cycle, domain);
      }

      if (eventForecast <= 0)
      {
        subEvent->Execute(m_time, m_dt, m_cycle, 0, 0, domain);
      }

      // Check the exit flag
      // Note: Currently, this is only being used by the HaltEvent
      //       If it starts being used elsewhere it may need to be synchronized
      exitFlag += subEvent->GetExitFlag();
    }

    // Increment time/cycle, reset the subevent counter
    m_time += m_dt;
    ++m_cycle;
    m_currentSubEvent = 0;
  }

  // Cleanup
  GEOS_LOG_RANK_0("Cleaning up events");
  
  this->forSubGroups<EventBase>([&]( EventBase * subEvent ) -> void
  {
    subEvent->Cleanup(m_time, m_cycle, 0, 0, domain);
  });
#if HAVE_TRIBOLCOUPLING
  TribolCoupling::Cleanup() ;
#endif
}

} /* namespace geosx */
