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
 * @file PeriodicEvent.hpp
 */

#ifndef GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_
#define GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_

#include "managers/Events/EventBase.hpp"

namespace geosx
{
/**
 * @class PeriodicEvent
 *
 * An event type for periodic events (using either time or cycle as a basis).
 */
class PeriodicEvent : public EventBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  PeriodicEvent( const std::string & name, Group * const parent );

  /// Destructor
  virtual ~PeriodicEvent() override;

  /**
   * @brief Catalog name interface.
   * @return This type's catalog name.
   **/
  static string
  CatalogName()
  {
    return "PeriodicEvent";
  }

  /**
   * @copydoc EventBase::EstimateEventTiming()
   * @note Estimate the expected number of cycles until an event is expected to trigger.
   *       The event frequency can be specified in terms of:
   *        - time (timeFrequency > 0, units = seconds)
   *        - or cycle (cycleFrequency >= 0, units = cycles)
   * @note In addition, there is an optional function input that will be called if the
   * the nominal forecast (based on timing) is zero.
   */
  virtual void
  EstimateEventTiming( real64 const time,
                       real64 const dt,
                       integer const cycle,
                       dataRepository::Group * domain ) override;

  /**
   * @brief Determine if an optional function f should be called, and call it if so.
   * @param time The current simulation time.
   * @param dt The current simulation time increment.
   * @param cycle The current simulation cycle.
   * @param domain The DomainPartition upcast to a Group.
   * If the event forecast is zero, and an optional function (f) is specified, then
   * this method will be called to see if the event should be triggered or ignored.
   * For example, this could be used to periodically check the condition of the mesh,
   * and trigger a cleanup if necessary.
   *
   * If functionInputObject is not specified:
   *   - The argument to the function will be the current time
   *   - The event will be executed if f(t) >= eventThreshold
   *
   * If functionInputObject is specified:
   *   - The function will be called on the object, with the arguments given by functionInputSetname
   *   - The function manager will return a set of statistics
   *   - functionStatOption selects the statistic to compare against the eventThreshold (0 = min, 1 = average, 2 = max)
   *   - The event will be executed if f(object, arguments)[stat] >= eventThreshold
   */
  void
  CheckOptionalFunctionThreshold( real64 const time,
                                  real64 const dt,
                                  integer const cycle,
                                  dataRepository::Group * domain );

  /**
   * @copydoc EventBase::GetEventTypeDtRequest()
   * Grab the next time-step.  If requested, then limit the requested
   * dt to exactly match the time frequency
   */
  virtual real64
  GetEventTypeDtRequest( real64 const time ) override;

  /**
   * @copydoc ExecutableGroup::Cleanup()
   */
  virtual void
  Cleanup( real64 const time_n,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           dataRepository::Group * domain ) override;

  /// A pointer to an optional function
  dataRepository::Group * m_functionTarget;

  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr auto timeFrequencyString = "timeFrequency";
    static constexpr auto cycleFrequencyString = "cycleFrequency";
    static constexpr auto targetExactTimestepString = "targetExactTimestep";
    static constexpr auto functionNameString = "function";
    static constexpr auto functionInputObjectString = "object";
    static constexpr auto functionInputSetnameString = "set";
    static constexpr auto functionSetNamesString = "setNames";
    static constexpr auto functionStatOptionString = "stat";
    static constexpr auto eventThresholdString = "threshold";

    dataRepository::ViewKey timeFrequency = { "timeFrequency" };
    dataRepository::ViewKey cycleFrequency = { "cycleFrequency" };
    dataRepository::ViewKey targetExactTimestep = { "targetExactTimestep" };
    dataRepository::ViewKey functionName = { "function" };
    dataRepository::ViewKey functionInputObject = { "object" };
    dataRepository::ViewKey functionInputSetname = { "set" };
    dataRepository::ViewKey functionSetNames = { "setNames" };
    dataRepository::ViewKey functionStatOption = { "stat" };
    dataRepository::ViewKey eventThreshold = { "threshold" };
  } periodicEventViewKeys;
  /// @endcond

  /// The event time frequency
  real64 m_timeFrequency;
  /// The event cycle frequency
  integer m_cycleFrequency;
  /// Whether to target the exact timestep
  integer m_targetExactTimestep;
  /// The optional function's name
  string m_functionName;
  /// The name of the optional function input object
  string m_functionInputObject;
  /// The name of the optional function input set
  string m_functionInputSetname;
  /// The optional funciton's statistic option
  integer m_functionStatOption;
  /// The event threshold
  real64 m_eventThreshold;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTS_PERIODICEVENT_HPP_ */
