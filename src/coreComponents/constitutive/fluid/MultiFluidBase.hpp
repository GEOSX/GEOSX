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
  localIndex numElems() const { return m_phaseFraction.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_phaseFraction.size( 1 ); }

  /**
   * @brief Get number of fluid components.
   * @return number of components
   */
  GEOSX_HOST_DEVICE
  localIndex numComponents() const { return m_componentMolarWeight.size(); }

  /**
   * @brief Get number of fluid phases.
   * @return number of phases
   */
  GEOSX_HOST_DEVICE
  localIndex numPhases() const { return m_phaseFraction.size( 2 ); }

protected:

  MultiFluidBaseUpdate( arrayView1d< real64 const > const & componentMolarWeight,
                        bool const useMass,
                        arrayView3d< real64 > const & phaseFraction,
                        arrayView3d< real64 > const & dPhaseFraction_dPressure,
                        arrayView3d< real64 > const & dPhaseFraction_dTemperature,
                        arrayView4d< real64 > const & dPhaseFraction_dGlobalCompFraction,
                        arrayView3d< real64 > const & phaseDensity,
                        arrayView3d< real64 > const & dPhaseDensity_dPressure,
                        arrayView3d< real64 > const & dPhaseDensity_dTemperature,
                        arrayView4d< real64 > const & dPhaseDensity_dGlobalCompFraction,
                        arrayView3d< real64 > const & phaseMassDensity,
                        arrayView3d< real64 > const & dPhaseMassDensity_dPressure,
                        arrayView3d< real64 > const & dPhaseMassDensity_dTemperature,
                        arrayView4d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
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
                        arrayView3d< real64 > const & dTotalDensity_dGlobalCompFraction,
                        arrayView3d< real64 > const & phaseEnthalpy,
                        arrayView3d< real64 > const & dPhaseEnthalpy_dPressure,
                        arrayView3d< real64 > const & dPhaseEnthalpy_dTemperature,
                        arrayView4d< real64 > const & dPhaseEnthalpy_dGlobalCompFraction,
                        arrayView3d< real64 > const & phaseInternalEnergy,
                        arrayView3d< real64 > const & dPhaseInternalEnergy_dPressure,
                        arrayView3d< real64 > const & dPhaseInternalEnergy_dTemperature,
                        arrayView4d< real64 > const & dPhaseInternalEnergy_dGlobalCompFraction )
    : m_componentMolarWeight( componentMolarWeight ),
    m_useMass( useMass ),
    m_phaseFraction( phaseFraction ),
    m_dPhaseFraction_dPressure( dPhaseFraction_dPressure ),
    m_dPhaseFraction_dTemperature( dPhaseFraction_dTemperature ),
    m_dPhaseFraction_dGlobalCompFraction( dPhaseFraction_dGlobalCompFraction ),
    m_phaseDensity( phaseDensity ),
    m_dPhaseDensity_dPressure( dPhaseDensity_dPressure ),
    m_dPhaseDensity_dTemperature( dPhaseDensity_dTemperature ),
    m_dPhaseDensity_dGlobalCompFraction( dPhaseDensity_dGlobalCompFraction ),
    m_phaseMassDensity( phaseMassDensity ),
    m_dPhaseMassDensity_dPressure( dPhaseMassDensity_dPressure ),
    m_dPhaseMassDensity_dTemperature( dPhaseMassDensity_dTemperature ),
    m_dPhaseMassDensity_dGlobalCompFraction( dPhaseMassDensity_dGlobalCompFraction ),
    m_phaseViscosity( phaseViscosity ),
    m_dPhaseViscosity_dPressure( dPhaseViscosity_dPressure ),
    m_dPhaseViscosity_dTemperature( dPhaseViscosity_dTemperature ),
    m_dPhaseViscosity_dGlobalCompFraction( dPhaseViscosity_dGlobalCompFraction ),
    m_phaseCompFraction( phaseCompFraction ),
    m_dPhaseCompFraction_dPressure( dPhaseCompFraction_dPressure ),
    m_dPhaseCompFraction_dTemperature( dPhaseCompFraction_dTemperature ),
    m_dPhaseCompFraction_dGlobalCompFraction( dPhaseCompFraction_dGlobalCompFraction ),
    m_totalDensity( totalDensity ),
    m_dTotalDensity_dPressure( dTotalDensity_dPressure ),
    m_dTotalDensity_dTemperature( dTotalDensity_dTemperature ),
    m_dTotalDensity_dGlobalCompFraction( dTotalDensity_dGlobalCompFraction ),
    m_phaseEnthalpy( phaseEnthalpy ),
    m_dPhaseEnthalpy_dPressure( dPhaseEnthalpy_dPressure ),
    m_dPhaseEnthalpy_dTemperature( dPhaseEnthalpy_dTemperature ),
    m_dPhaseEnthalpy_dGlobalCompFraction( dPhaseEnthalpy_dGlobalCompFraction ),
    m_phaseInternalEnergy( phaseInternalEnergy ),
    m_dPhaseInternalEnergy_dPressure( dPhaseInternalEnergy_dPressure ),
    m_dPhaseInternalEnergy_dTemperature( dPhaseInternalEnergy_dTemperature ),
    m_dPhaseInternalEnergy_dGlobalCompFraction( dPhaseInternalEnergy_dGlobalCompFraction )
  {}

  /// Default copy constructor
  MultiFluidBaseUpdate( MultiFluidBaseUpdate const & ) = default;

  /// Default move constructor
  MultiFluidBaseUpdate( MultiFluidBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  MultiFluidBaseUpdate & operator=( MultiFluidBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  MultiFluidBaseUpdate & operator=( MultiFluidBaseUpdate && ) = delete;

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

  arrayView3d< real64 > m_phaseMassDensity;
  arrayView3d< real64 > m_dPhaseMassDensity_dPressure;
  arrayView3d< real64 > m_dPhaseMassDensity_dTemperature;
  arrayView4d< real64 > m_dPhaseMassDensity_dGlobalCompFraction;

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

  arrayView3d< real64 > m_phaseEnthalpy;
  arrayView3d< real64 > m_dPhaseEnthalpy_dPressure;
  arrayView3d< real64 > m_dPhaseEnthalpy_dTemperature;
  arrayView4d< real64 > m_dPhaseEnthalpy_dGlobalCompFraction;

  arrayView3d< real64 > m_phaseInternalEnergy;
  arrayView3d< real64 > m_dPhaseInternalEnergy_dPressure;
  arrayView3d< real64 > m_dPhaseInternalEnergy_dTemperature;
  arrayView4d< real64 > m_dPhaseInternalEnergy_dGlobalCompFraction;

private:

  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & phaseMassDensity,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        real64 & totalDensity ) const = 0;

  virtual void compute( real64 const pressure,
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
                        arraySlice1d< real64 > const & phaseMassDensity,
                        arraySlice1d< real64 > const & dPhaseMassDensity_dPressure,
                        arraySlice1d< real64 > const & dPhaseMassDensity_dTemperature,
                        arraySlice2d< real64 > const & dPhaseMassDensity_dGlobalCompFraction,
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
                        arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseEnthalpy,
                        arraySlice1d< real64 > const & dPhaseEnthalpy_dPressure,
                        arraySlice1d< real64 > const & dPhaseEnthalpy_dTemperature,
                        arraySlice2d< real64 > const & dPhaseEnthalpy_dGlobalCompFraction,
                        arraySlice1d< real64 > const & phaseInternalEnergy,
                        arraySlice1d< real64 > const & dPhaseInternalEnergy_dPressure,
                        arraySlice1d< real64 > const & dPhaseInternalEnergy_dTemperature,
                        arraySlice2d< real64 > const & dPhaseInternalEnergy_dGlobalCompFraction ) const = 0;

  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition ) const = 0;

};

class MultiFluidBase : public ConstitutiveBase
{
public:

  MultiFluidBase( string const & name, Group * const parent );

  virtual ~MultiFluidBase() override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
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
  localIndex numFluidComponents() const { return m_componentNames.size(); }

  /**
   * @param ic component index
   * @return name of ic-th fluid component
   */
  arrayView1d< string const > componentNames() const { return m_componentNames; }

  /**
   * @return number of fluid phases in the model
   */
  localIndex numFluidPhases() const { return m_phaseNames.size(); }

  /**
   * @param ip phase index
   * @return name of ip-th fluid phase
   */
  arrayView1d< string const > phaseNames() const { return m_phaseNames; }

  /**
   * @brief Get the mass flag.
   * @return boolean value indicating whether the model is using mass-based quantities (as opposed to mole-based)
   */
  bool getMassFlag() const;

  /**
   * @brief Set the mass flag.
   * @param flag boolean value indicating whether the model should use mass-based quantities (as opposed to mole-based)
   *
   * @note This affects both input (compositions) and output quantities. The flag should be set prior to calling
   * any compute or state update methods.
   */
  void setMassFlag( bool flag );

  arrayView3d< real64 const > phaseFraction() const { return m_phaseFraction; }
  arrayView3d< real64 const > dPhaseFraction_dPressure() const { return m_dPhaseFraction_dPressure; }
  arrayView3d< real64 const > dPhaseFraction_dTemperature() const { return m_dPhaseFraction_dTemperature; }
  arrayView4d< real64 const > dPhaseFraction_dGlobalCompFraction() const { return m_dPhaseFraction_dGlobalCompFraction; }

  arrayView3d< real64 const > phaseDensity() const { return m_phaseDensity; }
  arrayView3d< real64 const > dPhaseDensity_dPressure() const { return m_dPhaseDensity_dPressure; }
  arrayView3d< real64 const > dPhaseDensity_dTemperature() const { return m_dPhaseDensity_dTemperature; }
  arrayView4d< real64 const > dPhaseDensity_dGlobalCompFraction() const { return m_dPhaseDensity_dGlobalCompFraction; }

  arrayView3d< real64 const > phaseMassDensity() const { return m_phaseMassDensity; }
  arrayView3d< real64 const > dPhaseMassDensity_dPressure() const { return m_dPhaseMassDensity_dPressure; }
  arrayView3d< real64 const > dPhaseMassDensity_dTemperature() const { return m_dPhaseMassDensity_dTemperature; }
  arrayView4d< real64 const > dPhaseMassDensity_dGlobalCompFraction() const { return m_dPhaseMassDensity_dGlobalCompFraction; }

  arrayView3d< real64 const > phaseViscosity() const { return m_phaseViscosity; }
  arrayView3d< real64 const > dPhaseViscosity_dPressure() const { return m_dPhaseViscosity_dPressure; }
  arrayView3d< real64 const > dPhaseViscosity_dTemperature() const { return m_dPhaseViscosity_dTemperature; }
  arrayView4d< real64 const > dPhaseViscosity_dGlobalCompFraction() const { return m_dPhaseViscosity_dGlobalCompFraction; }

  arrayView4d< real64 const > phaseCompFraction() const { return m_phaseCompFraction; }
  arrayView4d< real64 const > dPhaseCompFraction_dPressure() const { return m_dPhaseCompFraction_dPressure; }
  arrayView4d< real64 const > dPhaseCompFraction_dTemperature() const { return m_dPhaseCompFraction_dTemperature; }
  arrayView5d< real64 const > dPhaseCompFraction_dGlobalCompFraction() const { return m_dPhaseCompFraction_dGlobalCompFraction; }

  arrayView2d< real64 const > totalDensity() const { return m_totalDensity; }
  arrayView2d< real64 const > dTotalDensity_dPressure() const { return m_dTotalDensity_dPressure; }
  arrayView2d< real64 const > dTotalDensity_dTemperature() const { return m_dTotalDensity_dTemperature; }
  arrayView3d< real64 const > dTotalDensity_dGlobalCompFraction() const { return m_dTotalDensity_dGlobalCompFraction; }

  arrayView3d< real64 const > phaseEnthalpy() const { return m_phaseEnthalpy; }
  arrayView3d< real64 const > dPhaseEnthalpy_dPressure() const { return m_dPhaseEnthalpy_dPressure; }
  arrayView3d< real64 const > dPhaseEnthalpy_dTemperature() const { return m_dPhaseEnthalpy_dTemperature; }
  arrayView4d< real64 const > dPhaseEnthalpy_dGlobalCompFraction() const { return m_dPhaseEnthalpy_dGlobalCompFraction; }

  arrayView3d< real64 const > phaseInternalEnergy() const { return m_phaseInternalEnergy; }
  arrayView3d< real64 const > dPhaseInternalEnergy_dPressure() const { return m_dPhaseInternalEnergy_dPressure; }
  arrayView3d< real64 const > dPhaseInternalEnergy_dTemperature() const { return m_dPhaseInternalEnergy_dTemperature; }
  arrayView4d< real64 const > dPhaseInternalEnergy_dGlobalCompFraction() const { return m_dPhaseInternalEnergy_dGlobalCompFraction; }

  struct viewKeyStruct : ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * componentNamesString() { return "componentNames"; }
    static constexpr char const * componentMolarWeightString() { return "componentMolarWeight"; }

    static constexpr char const * phaseNamesString() { return "phaseNames"; }

    static constexpr char const * phaseFractionString() { return "phaseFraction"; } // xi_p
    static constexpr char const * dPhaseFraction_dPressureString() { return "dPhaseFraction_dPressure"; } // dXi_p/dP
    static constexpr char const * dPhaseFraction_dTemperatureString() { return "dPhaseFraction_dTemperature"; } // dXi_p/dT
    static constexpr char const * dPhaseFraction_dGlobalCompFractionString() { return "dPhaseFraction_dGlobalCompFraction"; } // dXi_p/dz

    static constexpr char const * phaseDensityString() { return "phaseDensity"; } // rho_p
    static constexpr char const * dPhaseDensity_dPressureString() { return "dPhaseDensity_dPressure"; } // dRho_p/dP
    static constexpr char const * dPhaseDensity_dTemperatureString() { return "dPhaseDensity_dTemperature"; } // dRho_p/dT
    static constexpr char const * dPhaseDensity_dGlobalCompFractionString() { return "dPhaseDensity_dGlobalCompFraction"; } // dRho_p/dz

    static constexpr char const * phaseMassDensityString() { return "phaseMassDensity"; } // rho_p
    static constexpr char const * dPhaseMassDensity_dPressureString() { return "dPhaseMassDensity_dPressure"; } // dRho_p/dP
    static constexpr char const * dPhaseMassDensity_dTemperatureString() { return "dPhaseMassDensity_dTemperature"; } // dRho_p/dT
    static constexpr char const * dPhaseMassDensity_dGlobalCompFractionString() { return "dPhaseMassDensity_dGlobalCompFraction"; } // dRho_p/dz

    static constexpr char const * phaseViscosityString() { return "phaseViscosity"; } // mu_p
    static constexpr char const * dPhaseViscosity_dPressureString() { return "dPhaseViscosity_dPressure"; } // dMu_p/dP
    static constexpr char const * dPhaseViscosity_dTemperatureString() { return "dPhaseViscosity_dTemperature"; } // dMu_p/dT
    static constexpr char const * dPhaseViscosity_dGlobalCompFractionString() { return "dPhaseViscosity_dGlobalCompFraction"; } // dMu_p/dz

    static constexpr char const * phaseCompFractionString() { return "phaseCompFraction"; } // x_cp
    static constexpr char const * dPhaseCompFraction_dPressureString() { return "dPhaseCompFraction_dPressure"; } // dx_cp/dP
    static constexpr char const * dPhaseCompFraction_dTemperatureString() { return "dPhaseCompFraction_dTemperature"; } // dx_cp/dT
    static constexpr char const * dPhaseCompFraction_dGlobalCompFractionString() { return "dPhaseCompFraction_dGlobalCompFraction"; } // dx_cp/dz

    static constexpr char const * totalDensityString() { return "totalDensity"; } // rho_t
    static constexpr char const * dTotalDensity_dPressureString() { return "dTotalDensity_dPressure"; } // dRho_t/dP
    static constexpr char const * dTotalDensity_dTemperatureString() { return "dTotalDensity_dTemperature"; } // dRho_t/dT
    static constexpr char const * dTotalDensity_dGlobalCompFractionString() { return "dTotalDensity_dGlobalCompFraction"; } // dRho_t/dz

    static constexpr char const * phaseEnthalpyString()                            { return "phaseEnthalpy"; }                            // h_p
    static constexpr char const * dPhaseEnthalpy_dPressureString()                 { return "dPhaseEnthalpy_dPressure"; }                 // dH_p/dP
    static constexpr char const * dPhaseEnthalpy_dTemperatureString()              { return "dPhaseEnthalpy_dTemperature"; }              // dH_p/dT
    static constexpr char const * dPhaseEnthalpy_dGlobalCompFractionString()       { return "dPhaseEnthalpy_dGlobalCompFraction"; }       // dH_p/dz

    static constexpr char const * phaseInternalEnergyString()                      { return "phaseInternalEnergy"; }                      // u_p
    static constexpr char const * dPhaseInternalEnergy_dPressureString()           { return "dPhaseInternalEnergy_dPressure"; }           // dU_p/dP
    static constexpr char const * dPhaseInternalEnergy_dTemperatureString()        { return "dPhaseInternalEnergy_dTemperature"; }        // dU_p/dT
    static constexpr char const * dPhaseInternalEnergy_dGlobalCompFractionString() { return  "dPhaseInternalEnergy_dGlobalCompFraction"; } // dU_p/dz

    static constexpr char const * useMassString() { return "useMass"; }
  };

protected:

  virtual void postProcessInput() override;

  /**
   * @brief Function called internally to resize member arrays
   * @param size primary dimension (e.g. number of cells)
   * @param numPts secondary dimension (e.g. number of gauss points per cell)
   */
  void resizeFields( localIndex const size, localIndex const numPts );

  // flag indicating whether input/output component fractions are treated as mass fractions
  int m_useMass;

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

  array3d< real64 > m_phaseMassDensity;
  array3d< real64 > m_dPhaseMassDensity_dPressure;
  array3d< real64 > m_dPhaseMassDensity_dTemperature;
  array4d< real64 > m_dPhaseMassDensity_dGlobalCompFraction;

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

  array3d< real64 > m_phaseEnthalpy;
  array3d< real64 > m_dPhaseEnthalpy_dPressure;
  array3d< real64 > m_dPhaseEnthalpy_dTemperature;
  array4d< real64 > m_dPhaseEnthalpy_dGlobalCompFraction;

  array3d< real64 > m_phaseInternalEnergy;
  array3d< real64 > m_dPhaseInternalEnergy_dPressure;
  array3d< real64 > m_dPhaseInternalEnergy_dTemperature;
  array4d< real64 > m_dPhaseInternalEnergy_dGlobalCompFraction;

};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDBASE_HPP_
