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
 * @file MultiFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{
class MultiFluidBaseUpdate
{
public:
  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex
  numElems() const
  {
    return m_phaseFraction.size( 0 );
  }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex
  numGauss() const
  {
    return m_phaseFraction.size( 1 );
  }

  /**
   * @brief Get number of fluid components.
   * @return number of components
   */
  GEOSX_HOST_DEVICE
  localIndex
  numComponents() const
  {
    return m_componentMolarWeight.size();
  }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOSX_HOST_DEVICE
  localIndex
  numPhases() const
  {
    return m_phaseFraction.size( 2 );
  }

protected:
  MultiFluidBaseUpdate(
    arrayView1d< real64 const > const & componentMolarWeight,
    bool const useMass,
    arrayView3d< real64 > const & phaseFraction,
    arrayView3d< real64 > const & dPhaseFraction_dPressure,
    arrayView3d< real64 > const & dPhaseFraction_dTemperature,
    arrayView4d< real64 > const & dPhaseFraction_dGlobalCompFraction,
    arrayView3d< real64 > const & phaseDensity,
    arrayView3d< real64 > const & dPhaseDensity_dPressure,
    arrayView3d< real64 > const & dPhaseDensity_dTemperature,
    arrayView4d< real64 > const & dPhaseDensity_dGlobalCompFraction,
    arrayView3d< real64 > const & phaseViscosity,
    arrayView3d< real64 > const & dPhaseViscosity_dPressure,
    arrayView3d< real64 > const & dPhaseViscosity_dTemperature,
    arrayView4d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
    arrayView4d< real64 > const & phaseCompFraction,
    arrayView4d< real64 > const & dPhaseCompFraction_dPressure,
    arrayView4d< real64 > const & dPhaseCompFraction_dTemperature,
    arrayView5d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
    arrayView2d< real64 > const & totalDensity,
    arrayView2d< real64 > const & dTotalDensity_dPressure,
    arrayView2d< real64 > const & dTotalDensity_dTemperature,
    arrayView3d< real64 > const & dTotalDensity_dGlobalCompFraction ) :
    m_componentMolarWeight( componentMolarWeight ),
    m_useMass( useMass ),
    m_phaseFraction( phaseFraction ),
    m_dPhaseFraction_dPressure( dPhaseFraction_dPressure ),
    m_dPhaseFraction_dTemperature( dPhaseFraction_dTemperature ),
    m_dPhaseFraction_dGlobalCompFraction( dPhaseFraction_dGlobalCompFraction ),
    m_phaseDensity( phaseDensity ),
    m_dPhaseDensity_dPressure( dPhaseDensity_dPressure ),
    m_dPhaseDensity_dTemperature( dPhaseDensity_dTemperature ),
    m_dPhaseDensity_dGlobalCompFraction( dPhaseDensity_dGlobalCompFraction ),
    m_phaseViscosity( phaseViscosity ),
    m_dPhaseViscosity_dPressure( dPhaseViscosity_dPressure ),
    m_dPhaseViscosity_dTemperature( dPhaseViscosity_dTemperature ),
    m_dPhaseViscosity_dGlobalCompFraction( dPhaseViscosity_dGlobalCompFraction ),
    m_phaseCompFraction( phaseCompFraction ),
    m_dPhaseCompFraction_dPressure( dPhaseCompFraction_dPressure ),
    m_dPhaseCompFraction_dTemperature( dPhaseCompFraction_dTemperature ),
    m_dPhaseCompFraction_dGlobalCompFraction(
      dPhaseCompFraction_dGlobalCompFraction ),
    m_totalDensity( totalDensity ),
    m_dTotalDensity_dPressure( dTotalDensity_dPressure ),
    m_dTotalDensity_dTemperature( dTotalDensity_dTemperature ),
    m_dTotalDensity_dGlobalCompFraction( dTotalDensity_dGlobalCompFraction )
  {}

  /// Default copy constructor
  MultiFluidBaseUpdate( MultiFluidBaseUpdate const & ) = default;

  /// Default move constructor
  MultiFluidBaseUpdate( MultiFluidBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  MultiFluidBaseUpdate &
  operator=( MultiFluidBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  MultiFluidBaseUpdate &
  operator=( MultiFluidBaseUpdate && ) = delete;

  arrayView1d< real64 const > m_componentMolarWeight;

  bool m_useMass;

  arrayView3d< real64 > m_phaseFraction;
  arrayView3d< real64 > m_dPhaseFraction_dPressure;
  arrayView3d< real64 > m_dPhaseFraction_dTemperature;
  arrayView4d< real64 > m_dPhaseFraction_dGlobalCompFraction;

  arrayView3d< real64 > m_phaseDensity;
  arrayView3d< real64 > m_dPhaseDensity_dPressure;
  arrayView3d< real64 > m_dPhaseDensity_dTemperature;
  arrayView4d< real64 > m_dPhaseDensity_dGlobalCompFraction;

  arrayView3d< real64 > m_phaseViscosity;
  arrayView3d< real64 > m_dPhaseViscosity_dPressure;
  arrayView3d< real64 > m_dPhaseViscosity_dTemperature;
  arrayView4d< real64 > m_dPhaseViscosity_dGlobalCompFraction;

  arrayView4d< real64 > m_phaseCompFraction;
  arrayView4d< real64 > m_dPhaseCompFraction_dPressure;
  arrayView4d< real64 > m_dPhaseCompFraction_dTemperature;
  arrayView5d< real64 > m_dPhaseCompFraction_dGlobalCompFraction;

  arrayView2d< real64 > m_totalDensity;
  arrayView2d< real64 > m_dTotalDensity_dPressure;
  arrayView2d< real64 > m_dTotalDensity_dTemperature;
  arrayView3d< real64 > m_dTotalDensity_dGlobalCompFraction;

private:
  virtual void
  Compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const > const & composition,
           arraySlice1d< real64 > const & phaseFraction,
           arraySlice1d< real64 > const & phaseDensity,
           arraySlice1d< real64 > const & phaseViscosity,
           arraySlice2d< real64 > const & phaseCompFraction,
           real64 & totalDensity ) const = 0;

  virtual void
  Compute(
    real64 const pressure,
    real64 const temperature,
    arraySlice1d< real64 const > const & composition,
    arraySlice1d< real64 > const & phaseFraction,
    arraySlice1d< real64 > const & dPhaseFraction_dPressure,
    arraySlice1d< real64 > const & dPhaseFraction_dTemperature,
    arraySlice2d< real64 > const & dPhaseFraction_dGlobalCompFraction,
    arraySlice1d< real64 > const & phaseDensity,
    arraySlice1d< real64 > const & dPhaseDensity_dPressure,
    arraySlice1d< real64 > const & dPhaseDensity_dTemperature,
    arraySlice2d< real64 > const & dPhaseDensity_dGlobalCompFraction,
    arraySlice1d< real64 > const & phaseViscosity,
    arraySlice1d< real64 > const & dPhaseViscosity_dPressure,
    arraySlice1d< real64 > const & dPhaseViscosity_dTemperature,
    arraySlice2d< real64 > const & dPhaseViscosity_dGlobalCompFraction,
    arraySlice2d< real64 > const & phaseCompFraction,
    arraySlice2d< real64 > const & dPhaseCompFraction_dPressure,
    arraySlice2d< real64 > const & dPhaseCompFraction_dTemperature,
    arraySlice3d< real64 > const & dPhaseCompFraction_dGlobalCompFraction,
    real64 & totalDensity,
    real64 & dTotalDensity_dPressure,
    real64 & dTotalDensity_dTemperature,
    arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const = 0;

  virtual void
  Update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< real64 const > const & composition ) const = 0;
};

class MultiFluidBase : public ConstitutiveBase
{
public:
  MultiFluidBase( std::string const & name, Group * const parent );

  virtual ~MultiFluidBase() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** MultiFluid-specific interface

  /**
   * @brief Maximum supported number of fluid components (species)
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr localIndex MAX_NUM_COMPONENTS = 16;

  /**
   * @brief Maximum supported number of fluid phases
   *
   * @note This puts an upper bound on memory use, allowing to optimize code better
   */
  static constexpr localIndex MAX_NUM_PHASES = 4;

  /**
   * @return number of fluid components (species) in the model
   */
  localIndex
  numFluidComponents() const
  {
    return m_componentNames.size();
  }

  /**
   * @param ic component index
   * @return name of ic-th fluid component
   */
  arrayView1d< string const > const &
  componentNames() const
  {
    return m_componentNames;
  }

  /**
   * @return number of fluid phases in the model
   */
  localIndex
  numFluidPhases() const
  {
    return m_phaseNames.size();
  }

  /**
   * @param ip phase index
   * @return name of ip-th fluid phase
   */
  arrayView1d< string const > const &
  phaseNames() const
  {
    return m_phaseNames;
  }

  /**
   * @brief Get the mass flag.
   * @return boolean value indicating whether the model is using mass-based quantities (as opposed to mole-based)
   */
  bool
  getMassFlag() const;

  /**
   * @brief Set the mass flag.
   * @param flag boolean value indicating whether the model should use mass-based quantities (as opposed to mole-based)
   *
   * @note This affects both input (compositions) and output quantities. The flag should be set prior to calling
   * any compute or state update methods.
   */
  void
  setMassFlag( bool flag );

  arrayView3d< real64 const > const &
  phaseFraction() const
  {
    return m_phaseFraction;
  }
  arrayView3d< real64 const > const &
  dPhaseFraction_dPressure() const
  {
    return m_dPhaseFraction_dPressure;
  }
  arrayView3d< real64 const > const &
  dPhaseFraction_dTemperature() const
  {
    return m_dPhaseFraction_dTemperature;
  }
  arrayView4d< real64 const > const &
  dPhaseFraction_dGlobalCompFraction() const
  {
    return m_dPhaseFraction_dGlobalCompFraction;
  }

  arrayView3d< real64 const > const &
  phaseDensity() const
  {
    return m_phaseDensity;
  }
  arrayView3d< real64 const > const &
  dPhaseDensity_dPressure() const
  {
    return m_dPhaseDensity_dPressure;
  }
  arrayView3d< real64 const > const &
  dPhaseDensity_dTemperature() const
  {
    return m_dPhaseDensity_dTemperature;
  }
  arrayView4d< real64 const > const &
  dPhaseDensity_dGlobalCompFraction() const
  {
    return m_dPhaseDensity_dGlobalCompFraction;
  }

  arrayView3d< real64 const > const &
  phaseViscosity() const
  {
    return m_phaseViscosity;
  }
  arrayView3d< real64 const > const &
  dPhaseViscosity_dPressure() const
  {
    return m_dPhaseViscosity_dPressure;
  }
  arrayView3d< real64 const > const &
  dPhaseViscosity_dTemperature() const
  {
    return m_dPhaseViscosity_dTemperature;
  }
  arrayView4d< real64 const > const &
  dPhaseViscosity_dGlobalCompFraction() const
  {
    return m_dPhaseViscosity_dGlobalCompFraction;
  }

  arrayView4d< real64 const > const &
  phaseCompFraction() const
  {
    return m_phaseCompFraction;
  }
  arrayView4d< real64 const > const &
  dPhaseCompFraction_dPressure() const
  {
    return m_dPhaseCompFraction_dPressure;
  }
  arrayView4d< real64 const > const &
  dPhaseCompFraction_dTemperature() const
  {
    return m_dPhaseCompFraction_dTemperature;
  }
  arrayView5d< real64 const > const &
  dPhaseCompFraction_dGlobalCompFraction() const
  {
    return m_dPhaseCompFraction_dGlobalCompFraction;
  }

  arrayView2d< real64 const > const &
  totalDensity() const
  {
    return m_totalDensity;
  }
  arrayView2d< real64 const > const &
  dTotalDensity_dPressure() const
  {
    return m_dTotalDensity_dPressure;
  }
  arrayView2d< real64 const > const &
  dTotalDensity_dTemperature() const
  {
    return m_dTotalDensity_dTemperature;
  }
  arrayView3d< real64 const > const &
  dTotalDensity_dGlobalCompFraction() const
  {
    return m_dTotalDensity_dGlobalCompFraction;
  }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr auto componentNamesString = "componentNames";
    static constexpr auto componentMolarWeightString = "componentMolarWeight";

    static constexpr auto phaseNamesString = "phaseNames";

    static constexpr auto phaseFractionString = "phaseFraction";  // xi_p
    static constexpr auto dPhaseFraction_dPressureString =
      "dPhaseFraction_dPressure";  // dXi_p/dP
    static constexpr auto dPhaseFraction_dTemperatureString =
      "dPhaseFraction_dTemperature";  // dXi_p/dT
    static constexpr auto dPhaseFraction_dGlobalCompFractionString =
      "dPhaseFraction_dGlobalCompFraction";  // dXi_p/dz

    static constexpr auto phaseDensityString = "phaseDensity";  // rho_p
    static constexpr auto dPhaseDensity_dPressureString =
      "dPhaseDensity_dPressure";  // dRho_p/dP
    static constexpr auto dPhaseDensity_dTemperatureString =
      "dPhaseDensity_dTemperature";  // dRho_p/dT
    static constexpr auto dPhaseDensity_dGlobalCompFractionString =
      "dPhaseDensity_dGlobalCompFraction";  // dRho_p/dz

    static constexpr auto phaseViscosityString = "phaseViscosity";  // mu_p
    static constexpr auto dPhaseViscosity_dPressureString =
      "dPhaseViscosity_dPressure";  // dMu_p/dP
    static constexpr auto dPhaseViscosity_dTemperatureString =
      "dPhaseViscosity_dTemperature";  // dMu_p/dT
    static constexpr auto dPhaseViscosity_dGlobalCompFractionString =
      "dPhaseViscosity_dGlobalCompFraction";  // dMu_p/dz

    static constexpr auto phaseCompFractionString = "phaseCompFraction";  // x_cp
    static constexpr auto dPhaseCompFraction_dPressureString =
      "dPhaseCompFraction_dPressure";  // dx_cp/dP
    static constexpr auto dPhaseCompFraction_dTemperatureString =
      "dPhaseCompFraction_dTemperature";  // dx_cp/dT
    static constexpr auto dPhaseCompFraction_dGlobalCompFractionString =
      "dPhaseCompFraction_dGlobalCompFraction";  // dx_cp/dz

    static constexpr auto totalDensityString = "totalDensity";  // rho_t
    static constexpr auto dTotalDensity_dPressureString =
      "dTotalDensity_dPressure";  // dRho_t/dP
    static constexpr auto dTotalDensity_dTemperatureString =
      "dTotalDensity_dTemperature";  // dRho_t/dT
    static constexpr auto dTotalDensity_dGlobalCompFractionString =
      "dTotalDensity_dGlobalCompFraction";  // dRho_t/dz
  } viewKeysMultiFluidBase;

protected:
  virtual void
  PostProcessInput() override;

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void
  ResizeFields( localIndex const size, localIndex const numPts );

  // flag indicating whether input/output component fractions are treated as mass fractions
  bool m_useMass;

  // general fluid composition information

  array1d< string > m_componentNames;
  array1d< real64 > m_componentMolarWeight;
  array1d< string > m_phaseNames;

  // constitutive data

  array3d< real64 > m_phaseFraction;
  array3d< real64 > m_dPhaseFraction_dPressure;
  array3d< real64 > m_dPhaseFraction_dTemperature;
  array4d< real64 > m_dPhaseFraction_dGlobalCompFraction;

  array3d< real64 > m_phaseDensity;
  array3d< real64 > m_dPhaseDensity_dPressure;
  array3d< real64 > m_dPhaseDensity_dTemperature;
  array4d< real64 > m_dPhaseDensity_dGlobalCompFraction;

  array3d< real64 > m_phaseViscosity;
  array3d< real64 > m_dPhaseViscosity_dPressure;
  array3d< real64 > m_dPhaseViscosity_dTemperature;
  array4d< real64 > m_dPhaseViscosity_dGlobalCompFraction;

  array4d< real64 > m_phaseCompFraction;
  array4d< real64 > m_dPhaseCompFraction_dPressure;
  array4d< real64 > m_dPhaseCompFraction_dTemperature;
  array5d< real64 > m_dPhaseCompFraction_dGlobalCompFraction;

  array2d< real64 > m_totalDensity;
  array2d< real64 > m_dTotalDensity_dPressure;
  array2d< real64 > m_dTotalDensity_dTemperature;
  array3d< real64 > m_dTotalDensity_dGlobalCompFraction;
};

}  //namespace constitutive

}  //namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
