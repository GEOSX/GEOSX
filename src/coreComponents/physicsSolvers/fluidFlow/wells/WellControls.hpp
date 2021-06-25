/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/*
 * @file WellControls.hpp
 */


#ifndef GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
#define GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "functions/TableFunction.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
static constexpr auto wellControls = "WellControls";
}
}


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
  enum class Type : integer
  {
    PRODUCER,  /**< A production well */
    INJECTOR   /**< An injection well */
  };


  /** Types of well controls
   * Used to specifiy a well's operating conditions
   */
  enum class Control : integer
  {
    BHP,  /**< The well operates at a specified bottom hole pressure (BHP) */
    PHASEVOLRATE, /**< The well operates at a specified phase volumetric flow rate */
    TOTALVOLRATE, /**< The well operates at a specified total volumetric flow rate */
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
  WellControls & operator=( WellControls const & ) = delete;

  /**
   * @brief Deleted move operator.
   * @return a reference to a perforation object
   */
  WellControls & operator=( WellControls && ) = delete;

  ///@}

  /**
   * @name Getters / Setters
   */
  ///@{

  /**
   * @brief Get the well type (injector or producer).
   * @return a well Type enum
   */
  Type getType() const { return m_type; }

  /**
   * @brief Set the control type to BHP and set a numerical value for the control.
   * @param[in] val value for the BHP control
   */
  void switchToBHPControl( real64 const & val );

  /**
   * @brief Set the control type to total rate and set a numerical value for the control.
   * @param[in] val value for the total volumetric rate
   */
  void switchToTotalRateControl( real64 const & val );

  /**
   * @brief Set the control type to phase rate and set a numerical value for the control.
   * @param[in] val value for the phase volumetric rate
   */
  void switchToPhaseRateControl( real64 const & val );

  /**
   * @brief Get the control type for the well.
   * @return the Control enum enforced at the well
   */
  Control getControl() const { return m_currentControl; }

  /**
   * @brief Getter for the reference elevation where the BHP control is enforced
   * @return the reference elevation
   */
  real64 getReferenceElevation() const { return m_refElevation; }

  /**
   * @brief Getter for the reference gravity coefficient
   * @return the reference gravity coefficient
   */
  real64 getReferenceGravityCoef() const { return m_refGravCoef; }

  /**
   * @brief Setter for the reference gravity
   */
  void setReferenceGravityCoef( real64 const & refGravCoef ) { m_refGravCoef = refGravCoef; }


  /**
   * @brief Get the target Bottom Hole Pressure value.
   * @return a value for the target Bottom Hole Pressure
   */
//  const real64 & getTargetBHP() const { return m_targetBHP; }
  real64 getTargetBHP( real64 const & currentTime )const
  {
    return m_targetBHPTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target total rate
   * @return the target total rate
   */
  real64 getTargetTotalRate( real64 const & currentTime )const
  {
    return m_targetTotalRateTable->evaluate( &currentTime );
  }

  /**
   * @brief Get the target phase rate
   * @return the target phase rate
   */
//  const real64 & getTargetPhaseRate() const { return m_targetPhaseRate; }
  real64 getTargetPhaseRate( real64 const & currentTime )const
  {
    return m_targetPhaseRateTable->evaluate( &currentTime );
  }
  /**
   * @brief Get the target phase name
   * @return the target phase name
   */
  const string & getTargetPhaseName() const { return m_targetPhaseName; }

  /**
   * @brief Const accessor for the composition of the injection rate
   * @return a global component fraction vector
   */
  arrayView1d< real64 const > getInjectionStream() const { return m_injectionStream; }

  /**
   * @brief Getter for the flag specifying whether we check rates at surface or reservoir conditions
   * @return 1 if we use surface conditions, and 0 otherwise
   */
  integer useSurfaceConditions() const { return m_useSurfaceConditions; }

  /**
   * @brief Getter for the surface pressure when m_useSurfaceConditions == 1
   * @return the surface pressure
   */
  const real64 & getSurfacePressure() const { return m_surfacePres; }

  /**
   * @brief Getter for the surface temperature when m_useSurfaceConditions == 1
   * @return the surface temperature
   */
  const real64 & getSurfaceTemperature() const { return m_surfaceTemp; }


  ///@}

  /**
   * @brief Struct to serve as a container for variable strings and keys.
   * @struct viewKeyStruct
   */
  struct viewKeyStruct
  {
    /// String key for the well reference elevation (for BHP control)
    static constexpr char const * refElevString() { return "referenceElevation"; }
    /// String key for the well type
    static constexpr char const * typeString() { return "type"; }
    /// String key for the well control
    static constexpr char const * controlString() { return "control"; }
    /// String key for the well target BHP
    static constexpr char const * targetBHPString() { return "targetBHP"; }
    /// String key for the well target rate
    static constexpr char const * targetTotalRateString() { return "targetTotalRate"; }
    /// String key for the well target phase rate
    static constexpr char const * targetPhaseRateString() { return "targetPhaseRate"; }
    /// String key for the well target phase name
    static constexpr char const * targetPhaseNameString() { return "targetPhaseName"; }
    /// String key for the well injection stream
    static constexpr char const * injectionStreamString() { return "injectionStream"; }
    /// String key for checking the rates at surface conditions
    static constexpr char const * useSurfaceConditionsString() { return "useSurfaceConditions"; }
    /// String key for the surface pressure
    static constexpr char const * surfacePressureString() { return "surfacePressure"; }
    /// String key for the surface temperature
    static constexpr char const * surfaceTemperatureString() { return "surfaceTemperature"; }
    /// string key for total rate table name
    static constexpr char const * targetTotalRateTableNameString() { return "targetTotalRateTableName"; }
    /// string key for phase rate table name
    static constexpr char const * targetPhaseRateTableNameString() { return "targetPhaseRateTableName"; }
    /// string key for BHP table name
    static constexpr char const * targetBHPTableNameString() { return "targetBHPTableName"; }
    /// ViewKey for the reference elevation
    dataRepository::ViewKey referenceElevation   = { refElevString() };
    /// ViewKey for the well type
    dataRepository::ViewKey type                 = { typeString() };
    /// ViewKey for the well control
    dataRepository::ViewKey control              = { controlString() };
    /// ViewKey for the well target BHP
    dataRepository::ViewKey targetBHP            = { targetBHPString() };
    /// ViewKey for the well target rate
    dataRepository::ViewKey targetTotalRate      = { targetTotalRateString() };
    /// ViewKey for the well target phase rate
    dataRepository::ViewKey targetPhaseRate      = { targetPhaseRateString() };
    /// ViewKey for the well target phase name
    dataRepository::ViewKey targetPhaseName      = { targetPhaseNameString() };
    /// ViewKey for the well injection stream
    dataRepository::ViewKey injectionStream      = { injectionStreamString() };
    /// ViewKey for the surface conditions flag
    dataRepository::ViewKey useSurfaceConditions = { useSurfaceConditionsString() };
    /// ViewKey for the surface pressure
    dataRepository::ViewKey surfacePressure      = { surfacePressureString() };
    /// ViewKey for the surface temperature
    dataRepository::ViewKey surfaceTemperature   = { surfaceTemperatureString() };
    /// ViewKey for the total rate table name
    dataRepository::ViewKey targetTotalRateTableName = { targetTotalRateTableNameString() };
    /// ViewKey for the phase rate table name
    dataRepository::ViewKey targetPhaseRateTableName = { targetPhaseRateTableNameString() };
    /// ViewKey for the BHP table name
    dataRepository::ViewKey targetBHPTableName       = { targetBHPTableNameString() };
  }
  /// ViewKey struct for the WellControls class
  viewKeysWellControls;

protected:

  virtual void postProcessInput() override;

  virtual void initializePostInitialConditionsPreSubGroups() override;

private:

  /// Well type (as Type enum)
  Type m_type;

  /// Reference elevation
  real64 m_refElevation;

  /// Gravity coefficient of the reference elevation
  real64 m_refGravCoef;

  /// Well controls as a Control enum
  Control m_currentControl;

  /// Target bottom hole pressure value
  real64 m_targetBHP;

  /// Target rate value
  real64 m_targetTotalRate;

  /// Target phase rate value
  real64 m_targetPhaseRate;

  /// Name of the targeted phase
  string m_targetPhaseName;

  /// Vector with global component fractions at the injector
  array1d< real64 >  m_injectionStream;

  /// Flag to decide whether rates are controlled at rates or surface conditions
  integer m_useSurfaceConditions;

  /// Surface pressure
  real64 m_surfacePres;

  /// Surface temperature
  real64 m_surfaceTemp;

  /// Total rate table name
  string m_targetTotalRateTableName;

  /// Phase rate table name
  string m_targetPhaseRateTableName;

  /// BHP table name
  string m_targetBHPTableName;

  /// Total rate table
  TableFunction * m_targetTotalRateTable;

  /// Phase rate table
  TableFunction * m_targetPhaseRateTable;

  /// BHP table
  TableFunction * m_targetBHPTable;
};

ENUM_STRINGS( WellControls::Type,
              "producer",
              "injector" );

ENUM_STRINGS( WellControls::Control,
              "BHP",
              "phaseVolRate",
              "totalVolRate" );

} //namespace geosx

#endif //GEOSX_PHYSICSSOLVERS_FLUIDFLOW_WELLS_WELLCONTROLS_HPP
