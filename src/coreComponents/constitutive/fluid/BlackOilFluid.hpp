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
 * @file BlackOilFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_

#include "constitutive/fluid/MultiFluidPVTPackageWrapper.hpp"

namespace PVTPackage
{
class BlackOilMultiphaseSystem;
}

namespace geosx
{
namespace constitutive
{
class BlackOilFluid : public MultiFluidPVTPackageWrapper
{
public:
  enum class FluidType
  {
    DeadOil,
    LiveOil
  };

  BlackOilFluid( std::string const & name, Group * const parent );

  virtual ~BlackOilFluid() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "BlackOilFluid";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  struct viewKeyStruct : MultiFluidPVTPackageWrapper::viewKeyStruct
  {
    static constexpr auto surfaceDensitiesString = "surfaceDensities";
    static constexpr auto tableFilesString = "tableFiles";
    static constexpr auto fluidTypeString = "fluidType";
  } viewKeysBlackOilFluid;

protected:
  virtual void
  PostProcessInput() override;

private:
  void
  createFluid() override;

  // Black-oil phase/component description
  array1d< real64 > m_surfaceDensities;

  // Black-oil table filenames
  path_array m_tableFiles;

  // Input string for type of black-oil fluid (live/dead)
  string m_fluidTypeString;

  // Type of black-oil fluid (live/dead)
  FluidType m_fluidType;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif  //GEOSX_CONSTITUTIVE_FLUID_BLACKOILFLUID_HPP_
