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
 * @file CompositionalMultiphaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_

#include "constitutive/fluid/MultiFluidPVTPackageWrapper.hpp"

namespace PVTPackage
{
class CompositionalMultiphaseSystem;
}

namespace geosx
{
namespace constitutive
{
class CompositionalMultiphaseFluid : public MultiFluidPVTPackageWrapper
{
public:
  CompositionalMultiphaseFluid( std::string const & name, Group * const parent );

  virtual ~CompositionalMultiphaseFluid() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "CompositionalMultiphaseFluid";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  struct viewKeyStruct : MultiFluidPVTPackageWrapper::viewKeyStruct
  {
    static constexpr auto equationsOfStateString = "equationsOfState";
    static constexpr auto componentCriticalPressureString =
      "componentCriticalPressure";
    static constexpr auto componentCriticalTemperatureString =
      "componentCriticalTemperature";
    static constexpr auto componentAcentricFactorString =
      "componentAcentricFactor";
    static constexpr auto componentVolumeShiftString = "componentVolumeShift";
    static constexpr auto componentBinaryCoeffString = "componentBinaryCoeff";
  } viewKeysCompositionalMultiphaseFluid;

protected:
  virtual void
  PostProcessInput() override;

private:
  void
  createFluid() override;

  // names of equations of state to use for each phase
  string_array m_equationsOfState;

  // standard EOS component input
  array1d< real64 > m_componentCriticalPressure;
  array1d< real64 > m_componentCriticalTemperature;
  array1d< real64 > m_componentAcentricFactor;
  array1d< real64 > m_componentVolumeShift;
  array2d< real64 > m_componentBinaryCoeff;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif  //GEOSX_CONSTITUTIVE_FLUID_COMPOSITIONALMULTIPHASEFLUID_HPP_
