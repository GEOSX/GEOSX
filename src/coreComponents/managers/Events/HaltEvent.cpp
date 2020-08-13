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

#include "HaltEvent.hpp"
#include <sys/time.h>

/**
 * @file HaltEvent.cpp
 */

namespace geosx
{
using namespace dataRepository;

HaltEvent::HaltEvent( const std::string & name, Group * const parent ) :
  EventBase( name, parent ),
  m_externalStartTime( 0.0 ),
  m_externalLastTime( 0.0 ),
  m_externalDt( 0.0 ),
  m_maxRuntime( 0.0 )
{
  timeval tim;
  gettimeofday( &tim, nullptr );
  m_externalStartTime = tim.tv_sec + ( tim.tv_usec / 1000000.0 );
  m_externalLastTime = m_externalStartTime;

  registerWrapper( viewKeyStruct::maxRuntimeString, &m_maxRuntime )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "The maximum allowable runtime for the job." );
}

HaltEvent::~HaltEvent()
{}

void
HaltEvent::EstimateEventTiming( real64 const GEOSX_UNUSED_PARAM( time ),
                                real64 const GEOSX_UNUSED_PARAM( dt ),
                                integer const GEOSX_UNUSED_PARAM( cycle ),
                                Group * GEOSX_UNUSED_PARAM( domain ) )
{
  // Check run time
  timeval tim;
  gettimeofday( &tim, nullptr );
  real64 currentTime = tim.tv_sec + ( tim.tv_usec / 1000000.0 );

  // Update values
  m_externalDt = currentTime - m_externalLastTime;
  m_externalLastTime = currentTime;
  integer forecast = static_cast< integer >(
    ( m_maxRuntime - ( currentTime - m_externalStartTime ) ) / m_externalDt );

  // The timing for the ranks may differ slightly, so synchronize
  // TODO: Only do the communication when you are close to the end?
#ifdef GEOSX_USE_MPI
  integer forecast_global;
  MPI_Allreduce( &forecast, &forecast_global, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD );
  forecast = forecast_global;
#endif

  setForecast( forecast );

  if( this->isReadyForExec() )
  {
    this->SetExitFlag( 1 );
  }
}

REGISTER_CATALOG_ENTRY( EventBase, HaltEvent, std::string const &, Group * const )
} /* namespace geosx */
