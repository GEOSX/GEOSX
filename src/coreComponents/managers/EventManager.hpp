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

#ifndef GEOSX_MANAGERS_EVENTMANAGER_HPP_
#define GEOSX_MANAGERS_EVENTMANAGER_HPP_

#include "dataRepository/Group.hpp"
#include "managers/Events/EventBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const Events( "Events" );
}
}  // namespace dataRepository

/**
 * @class EventManager
 *
 * A class for managing code events.
 */
class EventManager : public dataRepository::Group
{
public:
  /**
   * @brief Constructor for the EventManager
   * @param[in] name the name of the EventManager
   * @param[in] parent group this EventManager
   */
  EventManager( std::string const & name, Group * const parent );

  /**
   * @brief Default destructor for the EventManager
   */
  virtual ~EventManager() override;

  /**
   * @brief Create a child EventManager
   * @param[in] childKey the key of the Event to be added
   * @param[in] childName the name of the child to be added
   * @return the Event
   */
  virtual Group *
  CreateChild( string const & childKey,
               string const & childName ) override;

  /**
   * @brief This method is used to expand any catalogs in the data structure
   */
  virtual void
  ExpandObjectCatalogs() override;

  /**
   * @brief The main execution loop for the code.
   * @details During each cycle, it will:
   *   - Calculate the event forecast (number of cycles until its expected execution)
   *   - Signal an event to prepare (forecast == 1)
   *   - Execute an event (forecast == 0)
   *   - Determine dt for the next cycle
   *   - Advance time, cycle, etc.
   * @param[in] domain the current DomainPartition on which the Event will be ran
   */
  void
  Run( dataRepository::Group * domain );

  /**
   * @name viewKeyStruct/groupKeyStruct
   */
  ///@{
  /// @cond DO_NOT_DOCUMENT
  struct viewKeyStruct
  {
    static constexpr auto maxTimeString = "maxTime";
    static constexpr auto maxCycleString = "maxCycle";

    static constexpr auto timeString = "time";
    static constexpr auto dtString = "dt";
    static constexpr auto cycleString = "cycle";
    static constexpr auto currentSubEventString = "currentSubEvent";

    dataRepository::ViewKey time = { "time" };
    dataRepository::ViewKey dt = { "dt" };
    dataRepository::ViewKey cycle = { "cycle" };
    dataRepository::ViewKey maxTime = { "maxTime" };
    dataRepository::ViewKey maxCycle = { "maxCycle" };
    dataRepository::ViewKey currentSubEvent = { "currentSubEvent" };
  } viewKeys;
  /// @endcond
  ///@}

  /// Alias to access the object catalog for EventBase derived types.
  using CatalogInterface =
    dataRepository::CatalogInterface< EventBase, std::string const &, Group * const >;

  /// @copydoc dataRepository::Group::GetCatalog()
  static CatalogInterface::CatalogType &
  GetCatalog();

private:
  /// Max time for a simulation
  real64 m_maxTime;

  /// Maximum number of cycles for a simulation
  integer m_maxCycle;

  /// Simulation timestamp at the beginning of the cycle
  real64 m_time;

  /// Current timestep request
  real64 m_dt;

  /// Current cycle
  integer m_cycle;

  /// Current subevent index
  integer m_currentSubEvent;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_EVENTMANAGER_HPP_ */
