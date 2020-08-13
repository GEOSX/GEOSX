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
 * @file SoloEvent.cpp
 */

#include "SoloEvent.hpp"

namespace geosx
{
using namespace dataRepository;

SoloEvent::SoloEvent( const std::string & name, Group * const parent ) :
  EventBase( name, parent ),
  m_targetTime( -1.0 ),
  m_targetCycle( -1 ),
  m_targetExactTimestep( 0 )
{
  registerWrapper( viewKeyStruct::targetTimeString, &m_targetTime )
    ->setApplyDefaultValue( -1.0 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Targeted time to execute the event." );

  registerWrapper( viewKeyStruct::targetCycleString, &m_targetCycle )
    ->setApplyDefaultValue( -1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Targeted cycle to execute the event." );

  registerWrapper( viewKeyStruct::targetExactTimestepString, &m_targetExactTimestep )
    ->setApplyDefaultValue( 1 )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "If this option is set, the event will reduce its timestep requests to "
      "match the specified execution time exactly: dt_request = "
      "min(dt_request, t_target - time))." );
}

SoloEvent::~SoloEvent()
{}

void
SoloEvent::EstimateEventTiming( real64 const time,
                                real64 const dt,
                                integer const cycle,
                                Group * GEOSX_UNUSED_PARAM( domain ) )
{
  // Check event status
  if( m_lastCycle < 0 )
  {
    if( m_targetTime >= 0.0 )
    {
      if( dt <= 0 )
      {
        setIdle();
      }
      else
      {
        // Note: add a small value to this forecast to account for floating point errors
        real64 forecast = ( ( m_targetTime - time ) / dt ) + 1e-10;
        setForecast( static_cast< integer >( std::min( forecast, 1e9 ) ) );
      }
    }
    else
    {
      setForecast( m_targetCycle - cycle );
    }
  }
  else
  {
    setIdle();
  }
}

real64
SoloEvent::GetEventTypeDtRequest( real64 const time )
{
  real64 requestedDt = std::numeric_limits< real64 >::max();

  // Note: if m_lastCycle is set, then the event has already executed
  if( ( m_lastCycle < 0 ) && ( m_targetTime > 0 ) && ( m_targetExactTimestep > 0 ) )
  {
    // This extra step is necessary to prevent the event manager from
    // falling into a dt=0 loop
    real64 tmp_t = std::nextafter( time, time + 1.0 );
    if( tmp_t < m_targetTime )
    {
      requestedDt = std::min( requestedDt, m_targetTime - time );
    }
  }

  return requestedDt;
}

REGISTER_CATALOG_ENTRY( EventBase, SoloEvent, std::string const &, Group * const )
} /* namespace geosx */
