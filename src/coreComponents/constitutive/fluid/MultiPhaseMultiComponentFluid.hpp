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
 * @file MultiPhaseMultiComponentFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_

#include "constitutive/fluid/MultiFluidBase.hpp"

#include <memory>

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const multiPhaseMultiComponentFluid = "MultiPhaseMultiComponentFluid";
}
}  // namespace dataRepository

namespace PVTProps
{
class PVTFunction;
class FlashModel;
}  // namespace PVTProps

namespace constitutive
{
/**
 * @brief Kernel wrapper class for MultiPhaseMultiComponentFluid.
 * @note Not thread-safe, do not use with any parallel launch policy.
 */
class MultiPhaseMultiComponentFluidUpdate final : public MultiFluidBaseUpdate
{
public:
  MultiPhaseMultiComponentFluidUpdate(
    std::vector< std::shared_ptr< PVTProps::PVTFunction const > > const & phaseDensityFuns,
    std::vector< std::shared_ptr< PVTProps::PVTFunction const > > const & phaseViscosityFuns,
    std::shared_ptr< PVTProps::FlashModel const > const & flashModel,
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
    MultiFluidBaseUpdate( componentMolarWeight,
                          useMass,
                          phaseFraction,
                          dPhaseFraction_dPressure,
                          dPhaseFraction_dTemperature,
                          dPhaseFraction_dGlobalCompFraction,
                          phaseDensity,
                          dPhaseDensity_dPressure,
                          dPhaseDensity_dTemperature,
                          dPhaseDensity_dGlobalCompFraction,
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
                          dTotalDensity_dGlobalCompFraction ),
    m_phaseDensityFuns( phaseDensityFuns ),
    m_phaseViscosityFuns( phaseViscosityFuns ),
    m_flashModel( flashModel )
  {}

  /// Default copy constructor
  MultiPhaseMultiComponentFluidUpdate(
    MultiPhaseMultiComponentFluidUpdate const & ) = default;

  /// Default move constructor
  MultiPhaseMultiComponentFluidUpdate( MultiPhaseMultiComponentFluidUpdate && ) = default;

  /// Deleted copy assignment operator
  MultiPhaseMultiComponentFluidUpdate &
  operator=(
    MultiPhaseMultiComponentFluidUpdate const & ) = delete;

  /// Deleted move assignment operator
  MultiPhaseMultiComponentFluidUpdate &
  operator=(
    MultiPhaseMultiComponentFluidUpdate && ) = delete;

  virtual void
  Compute( real64 const pressure,
           real64 const temperature,
           arraySlice1d< real64 const > const & composition,
           arraySlice1d< real64 > const & phaseFraction,
           arraySlice1d< real64 > const & phaseDensity,
           arraySlice1d< real64 > const & phaseViscosity,
           arraySlice2d< real64 > const & phaseCompFraction,
           real64 & totalDensity ) const override;

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
    arraySlice1d< real64 > const & dTotalDensity_dGlobalCompFraction ) const override;

  GEOSX_FORCE_INLINE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const temperature,
          arraySlice1d< real64 const > const & composition ) const override
  {
    Compute( pressure,
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
             m_dTotalDensity_dGlobalCompFraction[k][q] );
  }

private:
  std::vector< std::shared_ptr< PVTProps::PVTFunction const > > m_phaseDensityFuns;
  std::vector< std::shared_ptr< PVTProps::PVTFunction const > > m_phaseViscosityFuns;
  std::shared_ptr< PVTProps::FlashModel const > m_flashModel;
};

class MultiPhaseMultiComponentFluid : public MultiFluidBase
{
public:
  MultiPhaseMultiComponentFluid( std::string const & name, Group * const parent );

  virtual ~MultiPhaseMultiComponentFluid() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return dataRepository::keys::multiPhaseMultiComponentFluid;
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = MultiPhaseMultiComponentFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper()
  {
    return KernelWrapper( m_phaseDensityFuns,
                          m_phaseViscosityFuns,
                          m_flashModel,
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
                          m_dTotalDensity_dGlobalCompFraction );
  }

  struct viewKeyStruct : MultiFluidBase::viewKeyStruct
  {
    static constexpr auto flashModelParaFileString = "flashModelParaFile";
    static constexpr auto phasePVTParaFilesString = "phasePVTParaFiles";
  } viewKeysMultiPhaseMultiComponentFluid;

protected:
  virtual void
  PostProcessInput() override;

  virtual void
  InitializePostSubGroups( Group * const group ) override;

private:
  void
  CreatePVTModels();

  // phase PVT parameter filenames
  path_array m_phasePVTParaFiles;

  string m_flashModelParaFile;

  // number of entries corrosponds to number of phases
  std::vector< std::shared_ptr< PVTProps::PVTFunction const > > m_phaseDensityFuns;
  std::vector< std::shared_ptr< PVTProps::PVTFunction const > > m_phaseViscosityFuns;

  std::shared_ptr< PVTProps::FlashModel const > m_flashModel;
};

}  //namespace constitutive

}  //namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_MULTIPHASEMULTICOMPONENTFLUID_HPP_
