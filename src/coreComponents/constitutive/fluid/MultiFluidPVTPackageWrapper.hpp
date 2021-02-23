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
 * @file MultiFluidPVTPackageWrapper.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"

// There is something wrong in the way we are building.
// If I include, I have to add the PVT dependency so far away...
// Therefore I forward declare...
namespace pvt
{
class MultiphaseSystem;

enum class PHASE_TYPE : int;
}

#include <memory>

namespace geosx
{

namespace constitutive
{

/**
 * @brief Kernel wrapper class for MultiFluidPVTPackage.
 * @note Not thread-safe, do not use with any parallel launch policy.
 */
class MultiFluidPVTPackageWrapperUpdate final : public MultiFluidBaseUpdate
{
public:

  MultiFluidPVTPackageWrapperUpdate( pvt::MultiphaseSystem & fluid,
                                     arrayView1d< pvt::PHASE_TYPE > const & phaseTypes,
                                     arrayView1d< real64 const > const & componentMolarWeight,
                                     bool useMass,
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
    : MultiFluidBaseUpdate( componentMolarWeight,
                            useMass,
                            phaseFraction,
                            dPhaseFraction_dPressure,
                            dPhaseFraction_dTemperature,
                            dPhaseFraction_dGlobalCompFraction,
                            phaseDensity,
                            dPhaseDensity_dPressure,
                            dPhaseDensity_dTemperature,
                            dPhaseDensity_dGlobalCompFraction,
                            phaseMassDensity,
                            dPhaseMassDensity_dPressure,
                            dPhaseMassDensity_dTemperature,
                            dPhaseMassDensity_dGlobalCompFraction,
                            phaseViscosity,
                            dPhaseViscosity_dPressure,
                            dPhaseViscosity_dTemperature,
                            dPhaseViscosity_dGlobalCompFraction,
                            phaseCompFraction,
                            dPhaseCompFraction_dPressure,
                            dPhaseCompFraction_dTemperature,
                            dPhaseCompFraction_dGlobalCompFraction,
                            totalDensity,
                            dTotalDensity_dPressure,
                            dTotalDensity_dTemperature,
                            dTotalDensity_dGlobalCompFraction,
                            phaseEnthalpy,
                            dPhaseEnthalpy_dPressure,
                            dPhaseEnthalpy_dTemperature,
                            dPhaseEnthalpy_dGlobalCompFraction,
                            phaseInternalEnergy,
                            dPhaseInternalEnergy_dPressure,
                            dPhaseInternalEnergy_dTemperature,
                            dPhaseInternalEnergy_dGlobalCompFraction ),
    m_fluid( fluid ),
    m_phaseTypes( phaseTypes )
  {}

  /// Default copy constructor
  MultiFluidPVTPackageWrapperUpdate( MultiFluidPVTPackageWrapperUpdate const & ) = default;

  /// Default move constructor
  MultiFluidPVTPackageWrapperUpdate( MultiFluidPVTPackageWrapperUpdate && ) = default;

  /// Deleted copy assignment operator
  MultiFluidPVTPackageWrapperUpdate & operator=( MultiFluidPVTPackageWrapperUpdate const & ) = delete;

  /// Deleted move assignment operator
  MultiFluidPVTPackageWrapperUpdate & operator=( MultiFluidPVTPackageWrapperUpdate && ) = delete;

  virtual void compute( real64 const pressure,
                        real64 const temperature,
                        arraySlice1d< real64 const > const & composition,
                        arraySlice1d< real64 > const & phaseFraction,
                        arraySlice1d< real64 > const & phaseDensity,
                        arraySlice1d< real64 > const & phaseMassDensity,
                        arraySlice1d< real64 > const & phaseViscosity,
                        arraySlice2d< real64 > const & phaseCompFraction,
                        real64 & totalDensity ) const override;

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
                        arraySlice2d< real64 > const & dPhaseInternalEnergy_dGlobalCompFraction ) const override;

  GEOSX_FORCE_INLINE
  virtual void update( localIndex const k,
                       localIndex const q,
                       real64 const pressure,
                       real64 const temperature,
                       arraySlice1d< real64 const > const & composition ) const override
  {
    compute( pressure,
             temperature,
             composition,
             m_phaseFraction[k][q],
             m_dPhaseFraction_dPressure[k][q],
             m_dPhaseFraction_dTemperature[k][q],
             m_dPhaseFraction_dGlobalCompFraction[k][q],
             m_phaseDensity[k][q],
             m_dPhaseDensity_dPressure[k][q],
             m_dPhaseDensity_dTemperature[k][q],
             m_dPhaseDensity_dGlobalCompFraction[k][q],
             m_phaseMassDensity[k][q],
             m_dPhaseMassDensity_dPressure[k][q],
             m_dPhaseMassDensity_dTemperature[k][q],
             m_dPhaseMassDensity_dGlobalCompFraction[k][q],
             m_phaseViscosity[k][q],
             m_dPhaseViscosity_dPressure[k][q],
             m_dPhaseViscosity_dTemperature[k][q],
             m_dPhaseViscosity_dGlobalCompFraction[k][q],
             m_phaseCompFraction[k][q],
             m_dPhaseCompFraction_dPressure[k][q],
             m_dPhaseCompFraction_dTemperature[k][q],
             m_dPhaseCompFraction_dGlobalCompFraction[k][q],
             m_totalDensity[k][q],
             m_dTotalDensity_dPressure[k][q],
             m_dTotalDensity_dTemperature[k][q],
             m_dTotalDensity_dGlobalCompFraction[k][q],
             m_phaseEnthalpy[k][q],
             m_dPhaseEnthalpy_dPressure[k][q],
             m_dPhaseEnthalpy_dTemperature[k][q],
             m_dPhaseEnthalpy_dGlobalCompFraction[k][q],
             m_phaseInternalEnergy[k][q],
             m_dPhaseInternalEnergy_dPressure[k][q],
             m_dPhaseInternalEnergy_dTemperature[k][q],
             m_dPhaseInternalEnergy_dGlobalCompFraction[k][q] );
  }

private:

  pvt::MultiphaseSystem & m_fluid;

  arrayView1d< pvt::PHASE_TYPE > m_phaseTypes;

};

class MultiFluidPVTPackageWrapper : public MultiFluidBase
{
public:

  MultiFluidPVTPackageWrapper( string const & name, Group * const parent );

  virtual ~MultiFluidPVTPackageWrapper() override;

  virtual std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                Group * const parent ) const override;

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = MultiFluidPVTPackageWrapperUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( *m_fluid,
                          m_phaseTypes,
                          m_componentMolarWeight,
                          m_useMass,
                          m_phaseFraction,
                          m_dPhaseFraction_dPressure,
                          m_dPhaseFraction_dTemperature,
                          m_dPhaseFraction_dGlobalCompFraction,
                          m_phaseDensity,
                          m_dPhaseDensity_dPressure,
                          m_dPhaseDensity_dTemperature,
                          m_dPhaseDensity_dGlobalCompFraction,
                          m_phaseMassDensity,
                          m_dPhaseMassDensity_dPressure,
                          m_dPhaseMassDensity_dTemperature,
                          m_dPhaseMassDensity_dGlobalCompFraction,
                          m_phaseViscosity,
                          m_dPhaseViscosity_dPressure,
                          m_dPhaseViscosity_dTemperature,
                          m_dPhaseViscosity_dGlobalCompFraction,
                          m_phaseCompFraction,
                          m_dPhaseCompFraction_dPressure,
                          m_dPhaseCompFraction_dTemperature,
                          m_dPhaseCompFraction_dGlobalCompFraction,
                          m_totalDensity,
                          m_dTotalDensity_dPressure,
                          m_dTotalDensity_dTemperature,
                          m_dTotalDensity_dGlobalCompFraction,
                          m_phaseEnthalpy,
                          m_dPhaseEnthalpy_dPressure,
                          m_dPhaseEnthalpy_dTemperature,
                          m_dPhaseEnthalpy_dGlobalCompFraction,
                          m_phaseInternalEnergy,
                          m_dPhaseInternalEnergy_dPressure,
                          m_dPhaseInternalEnergy_dTemperature,
                          m_dPhaseInternalEnergy_dGlobalCompFraction );
  }

protected:

  virtual void postProcessInput() override;

  virtual void initializePostSubGroups() override;

  /// function that populates m_fluid ptr; to be overriden by derived classes
  virtual void createFluid() = 0;

  /// PVTPackage fluid object
  std::unique_ptr< pvt::MultiphaseSystem > m_fluid;

  /// PVTPackage phase labels
  array1d< pvt::PHASE_TYPE > m_phaseTypes;
};

} //namespace constitutive

} //namespace geosx

#endif //GEOSX_CONSTITUTIVE_FLUID_MULTIFLUIDPVTPACKAGEWRAPPER_HPP_
