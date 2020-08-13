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

/*
 * @file WellControls.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP

#include "dataRepository/Group.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellControls = "WellControls";
}
}  // namespace dataRepository

/**
 * @class WellControls
 * @brief This class describes the controls used to operate a well.
 */
class WellControls : public dataRepository::Group
{
public:
  /** Type of wells
   * Either producer or injector.
   */
  enum class Type
  {
    PRODUCER, /**< A production well */
    INJECTOR  /**< An injection well */
  };

  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control
  {
    BHP,       /**< The well operates at a specified bottom hole pressure (BHP) */
    GASRATE,   /**< The well operates at a specified gas flow rate */
    OILRATE,   /**< The well operates at a specified oil flow rate */
    WATERRATE, /**< The well operates at a specified water flow rate */
    LIQUIDRATE /**< The well operates at a specified liquid flow rate (oil + water) */
  };

  /**
   * @name Constructor / Destructor
   */
  ///@{

  /**
   * @brief Constructor for WellControls Objects.
   * @param[in] name the name of this instantiation of WellControls in the repository
   * @param[in] parent the parent group of this instantiation of WellControls
   */
  explicit WellControls( string const & name, dataRepository::Group * const parent );

  /**
   * @brief Default destructor.
   */
  ~WellControls() override;

  /**
   * @brief Deleted default constructor.
   */
  WellControls() = delete;

  /**
   * @brief Deleted copy constructor.
   */
  WellControls( WellControls const & ) = delete;

  /**
   * @brief Deleted move constructor.
   */
  WellControls( WellControls && ) = delete;

  /**
   * @brief Deleted assignment operator.
   * @return a reference to a perforation object
   */
  WellControls &
  operator=( WellControls const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  WellControls &
  operator=( WellControls && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Set the reference well elem index where the control will be enforced.
   * @param[in] refIndex reference well element index where the control will be enforced
   */
  void
  SetReferenceWellElementIndex( localIndex refIndex )
  {
    m_refWellElemIndex = refIndex;
  }

  /**
   * @brief Get the reference well element index where the control will be enforced.
   * @return a localIndex value representing the reference well element index where the control will be enforced
   */
  localIndex const &
  GetReferenceWellElementIndex() const
  {
    return m_refWellElemIndex;
  }

  /**
   * @brief Get the well type (injector or producer).
   * @return a well Type enum
   */
  Type
  GetType() const
  {
    return m_type;
  }

  /**
   * @brief Set the control type and numerical value for a well.
   * @param[in] control a Control enum with the type of control that is enforced
   * @param[in] val value for the control (depending on the control type, can be a maximum bottom hole pressure, a
   * minimum water rate...)
   */
  void
  SetControl( Control control, real64 const & val );

  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  Control
  GetControl() const
  {
    return m_currentControl;
  }

  /**
   * @brief Get the target Bottom Hole Pressure value.
   * @return a value for the target Bottom Hole Pressure
   */
  const real64 &
  GetTargetBHP() const
  {
    return m_targetBHP;
  }

  /**
   * @brief Get the target rate
   * @return the target rate
   */
  const real64 &
  GetTargetRate() const
  {
    return m_targetRate;
  }

  /**
   * @brief Const accessor for the composition of the injection rate
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > const &
  GetInjectionStream() const
  {
    return m_injectionStream;
  }

  ///@}

  /// @cond DO_NOT_DOCUMENT
  void
  Debug() const;
  /// @endcond

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the reference index (currently unused)
    static constexpr auto refWellElemIndexString = "referenceWellElementIndex";
    /// String key for the well type
    static constexpr auto typeString = "type";
    /// String key for the well control
    static constexpr auto controlString = "control";
    /// String key for the well target BHP
    static constexpr auto targetBHPString = "targetBHP";
    /// String key for the well target rate
    static constexpr auto targetRateString = "targetRate";
    /// String key for the well injection stream
    static constexpr auto injectionStreamString = "injectionStream";
    /// ViewKey for the reference index (currently unused)
    dataRepository::ViewKey referenceIndex = { refWellElemIndexString };
    /// ViewKey for the well type
    dataRepository::ViewKey type = { typeString };
    /// ViewKey for the well control
    dataRepository::ViewKey control = { controlString };
    /// ViewKey for the well target BHP
    dataRepository::ViewKey targetBHP = { targetBHPString };
    /// ViewKey for the well target rate
    dataRepository::ViewKey targetRate = { targetRateString };
    /// ViewKey for the well injection stream
    dataRepository::ViewKey injectionStream = { injectionStreamString };
  }
  /// ViewKey struct for the WellControls class
  viewKeysWellControls;

protected:
  /**
   * @brief This function provides capability to post process input values prior to
   * any other initialization operations.
   */
  virtual void
  PostProcessInput() override;

  /**
   * @brief Called by InitializePostInitialConditions() prior to initializing sub-Groups.
   * @param[in] rootGroup A group that is passed in to the initialization functions
   *                  in order to facilitate the initialization.
   */
  virtual void
  InitializePostInitialConditions_PreSubGroups(
    Group * const rootGroup ) override;

private:
  /// Well type as string
  string m_typeString;

  /// Well type (as Type enum)
  Type m_type;

  /// Reference index (currently unused)
  localIndex m_refWellElemIndex;

  /// Well controls as string
  string m_inputControlString;

  /// Well controls as a Control enum
  Control m_currentControl;

  /// Target bottom hole pressure value
  real64 m_targetBHP;

  /// Target rate value
  real64 m_targetRate;

  /// Vector with global component fractions at the injector
  array1d< real64 > m_injectionStream;
};

}  //namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
