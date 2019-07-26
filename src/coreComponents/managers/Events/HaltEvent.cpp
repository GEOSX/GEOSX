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

#include "HaltEvent.hpp"
#include <sys/time.h>

/**
 * @file HaltEvent.cpp
 */

namespace geosx
{

using namespace dataRepository;


HaltEvent::HaltEvent( const std::string& name,
                      ManagedGroup * const parent ):
  EventBase(name,parent),
  m_externalStartTime(0.0),
  m_externalLastTime(0.0),
  m_externalDt(0.0),
  m_maxRuntime(0.0)
{
  timeval tim;
  gettimeofday(&tim, nullptr);
  m_externalStartTime = tim.tv_sec + (tim.tv_usec / 1000000.0);
  m_externalLastTime = m_externalStartTime;  

  RegisterViewWrapper(viewKeyStruct::maxRuntimeString, &m_maxRuntime, false )->
    setInputFlag(InputFlags::REQUIRED)->
    setDescription("max runtime");
}


HaltEvent::~HaltEvent()
{}


void HaltEvent::EstimateEventTiming(real64 const time,
                                     real64 const dt, 
                                     integer const cycle,
                                     ManagedGroup * domain)
{
  // Check run time
  timeval tim;
  gettimeofday(&tim, nullptr);
  real64 currentTime = tim.tv_sec + (tim.tv_usec / 1000000.0);

  // Update values
  m_externalDt = currentTime - m_externalLastTime;
  m_externalLastTime = currentTime;
  integer forecast = static_cast<integer>((m_maxRuntime - (currentTime - m_externalStartTime)) / m_externalDt);
  
  // The timing for the ranks may differ slightly, so synchronize
  // TODO: Only do the communication when you are close to the end?
#ifdef GEOSX_USE_MPI
    integer forecast_global;
    MPI_Allreduce(&forecast, &forecast_global, 1, MPI_INT, MPI_MIN, MPI_COMM_GEOSX);
    forecast = forecast_global;
#endif

  SetForecast(forecast);

  if (this->GetForecast() <= 0)
  {
    this->SetExitFlag(1);
  }
}


REGISTER_CATALOG_ENTRY( EventBase, HaltEvent, std::string const &, ManagedGroup * const )
} /* namespace geosx */
