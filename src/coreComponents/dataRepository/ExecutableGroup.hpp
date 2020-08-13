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
 * @file ExecutableGroup.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_EXECUTABLEGROUP_HPP_
#define GEOSX_DATAREPOSITORY_EXECUTABLEGROUP_HPP_

#include "common/DataTypes.hpp"
#include "Group.hpp"

namespace geosx
{
/**
 * @class ExecutableGroup
 *
 * Objects that are executable and/or able to request changes dt
 * should be derived from this type.
 */
class ExecutableGroup : public dataRepository::Group
{
public:
  using dataRepository::Group::Group;

  /**
   * @brief Main extension point of executable targets.
   * @param[in] time_n        current time level
   * @param[in] dt            time step to be taken
   * @param[in] cycleNumber   global cycle number
   * @param[in] eventCounter  index of event that triggered execution
   * @param[in] eventProgress fractional progress in current cycle
   * @param[in,out] domain    the physical domain up-casted to a Group.
   *
   * If the start criteria are satisfied, then the event manager
   * will call this method.
   */
  virtual void
  Execute( real64 const time_n,
           real64 const dt,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           dataRepository::Group * domain ) = 0;

  /**
   * @brief Inform the object that it expects to execute during the next timestep.
   * @param[in] time_n        current time level
   * @param[in] dt            time step to be taken
   * @param[in] cycle         global cycle number
   * @param[in,out] domain    the physical domain
   */
  virtual void
  SignalToPrepareForExecution( real64 const time_n,
                               real64 const dt,
                               integer const cycle,
                               dataRepository::Group * domain );
  /**
   * @brief Called as the code exits the main run loop.
   * @param[in] time_n        current time level
   * @param[in] cycleNumber   global cycle number
   * @param[in] eventCounter  index of event that triggered execution
   * @param[in] eventProgress fractional progress in current cycle
   * @param[in,out] domain    the physical domain
   */
  virtual void
  Cleanup( real64 const time_n,
           integer const cycleNumber,
           integer const eventCounter,
           real64 const eventProgress,
           dataRepository::Group * domain );

  /**
   * @brief Supplies the timestep request for this target to the event manager.
   * @param[in] time current time level
   * @return         desired time step size
   */
  virtual real64
  GetTimestepRequest( real64 const time )
  {
    GEOSX_UNUSED_VAR( time );
    return std::numeric_limits< integer >::max();
  }

  /**
   * @brief Set the timestep behavior for a target.
   * @param[in] behavior if positive, target does time stepping
   */
  void
  SetTimestepBehavior( integer const behavior )
  {
    m_timestepType = behavior;
  }

  /**
   * @brief Get the target's time step behavior.
   * @return @p >0 if target does time stepping, @p <=0 otherwise
   */
  integer
  GetTimestepBehavior()
  {
    return m_timestepType;
  }

private:
  integer m_timestepType = 0;
};

}  // namespace geosx

#endif /* GEOSX_DATAREPOSITORY_EXECUTABLEGROUP_HPP_ */
