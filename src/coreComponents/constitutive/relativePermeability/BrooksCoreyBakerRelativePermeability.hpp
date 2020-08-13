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
 * @file BrooksCoreyBakerRelativePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELPERM_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP
#define GEOSX_CONSTITUTIVE_RELPERM_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geosx
{
namespace constitutive
{
class BrooksCoreyBakerRelativePermeabilityUpdate final
  : public RelativePermeabilityBaseUpdate
{
public:
  BrooksCoreyBakerRelativePermeabilityUpdate(
    arrayView1d< real64 const > const & phaseMinVolumeFraction,
    arrayView1d< real64 const > const & waterOilRelPermExponent,
    arrayView1d< real64 const > const & waterOilRelPermMaxValue,
    arrayView1d< real64 const > const & gasOilRelPermExponent,
    arrayView1d< real64 const > const & gasOilRelPermMaxValue,
    real64 const volFracScale,
    arrayView1d< integer const > const & phaseTypes,
    arrayView1d< integer const > const & phaseOrder,
    arrayView3d< real64 > const & phaseRelPerm,
    arrayView4d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) :
    RelativePermeabilityBaseUpdate( phaseTypes,
                                    phaseOrder,
                                    phaseRelPerm,
                                    dPhaseRelPerm_dPhaseVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_waterOilRelPermExponent( waterOilRelPermExponent ),
    m_waterOilRelPermMaxValue( waterOilRelPermMaxValue ),
    m_gasOilRelPermExponent( gasOilRelPermExponent ),
    m_gasOilRelPermMaxValue( gasOilRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  /// Default copy constructor
  BrooksCoreyBakerRelativePermeabilityUpdate(
    BrooksCoreyBakerRelativePermeabilityUpdate const & ) = default;

  /// Default move constructor
  BrooksCoreyBakerRelativePermeabilityUpdate(
    BrooksCoreyBakerRelativePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  BrooksCoreyBakerRelativePermeabilityUpdate &
  operator=(
    BrooksCoreyBakerRelativePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrooksCoreyBakerRelativePermeabilityUpdate &
  operator=(
    BrooksCoreyBakerRelativePermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Compute(
    arraySlice1d< real64 const > const & phaseVolFraction,
    arraySlice1d< real64 > const & phaseRelPerm,
    arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    Compute( phaseVolFraction,
             m_phaseRelPerm[k][q],
             m_dPhaseRelPerm_dPhaseVolFrac[k][q] );
  }

private:
  /**
   * @brief Evaluate the Brooks-Corey relperm function for a given (scalar) phase saturation
   * @param[in] scaledVolFrac the scaled volume fraction for this phase
   * @param[in] dScaledVolFrac_dVolFrac the derivative of scaled volume fraction for this phase wrt to the volume
   * fraction
   * @param[out] relperm the relative permeability for this phase
   * @param[out] dRelPerm_dVolFrac the derivative of the relative permeability wrt to the volume fraction of the phase
   * @param[in] exponentInv the inverse of the exponent used in the VG model
   * @param[in] maxValue the endpoint relative permeability value
   *
   * This function evaluates the relperm function and its derivative at a given phase saturation
   * Reference: Eclipse technical description and Petrowiki
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  EvaluateBrooksCoreyFunction( real64 const & scaledVolFrac,
                               real64 const & dScaledVolFrac_dVolFrac,
                               real64 const & exponent,
                               real64 const & maxValue,
                               real64 & relPerm,
                               real64 & dRelPerm_dVolFrac );

  /**
   * @brief Interpolate the two-phase relperms to compute the three-phase relperm
   * @param[in] shiftedWaterVolFrac
   * @param[in] gasVolFrac
   * @param[out] threePhaseRelPerm
   * @param[out] dThreePhaseRelPerm_dVolFrac
   * @param[in] relPerm_wo
   * @param[in] dRelPerm_wo_dOilVolFrac
   * @param[in] relPerm_go
   * @param[in] dRelPerm_go_dOilVolFrac
   *
   * This function interpolates the two-phase relperms to compute the three-phase relperm
   * The interpolation is based on the modified Baker method, also used as default in Eclipse
   * Reference: Eclipse technical description and PetroWiki
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  InterpolateTwoPhaseRelPerms(
    real64 const & shiftedWaterVolFrac,
    real64 const & gasVolFrac,
    arraySlice1d< integer const > const & phaseOrder,
    real64 const & relPerm_wo,
    real64 const & dRelPerm_wo_dOilVolFrac,
    real64 const & relPerm_go,
    real64 const & dRelPerm_go_dOilVolFrac,
    real64 & threePhaseRelPerm,
    arraySlice1d< real64 > const & dThreePhaseRelPerm_dVolFrac );

  arrayView1d< real64 const > m_phaseMinVolumeFraction;

  arrayView1d< real64 const > m_waterOilRelPermExponent;
  arrayView1d< real64 const > m_waterOilRelPermMaxValue;

  arrayView1d< real64 const > m_gasOilRelPermExponent;
  arrayView1d< real64 const > m_gasOilRelPermMaxValue;

  real64 m_volFracScale;
};

class BrooksCoreyBakerRelativePermeability : public RelativePermeabilityBase
{
public:
  BrooksCoreyBakerRelativePermeability( std::string const & name,
                                        dataRepository::Group * const parent );

  virtual ~BrooksCoreyBakerRelativePermeability() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "BrooksCoreyBakerRelativePermeability";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyBakerRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper();

  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString =
      "phaseMinVolumeFraction";
    static constexpr auto waterOilRelPermExponentString =
      "waterOilRelPermExponent";
    static constexpr auto waterOilRelPermMaxValueString =
      "waterOilRelPermMaxValue";
    static constexpr auto gasOilRelPermExponentString = "gasOilRelPermExponent";
    static constexpr auto gasOilRelPermMaxValueString = "gasOilRelPermMaxValue";
  } viewKeysBrooksCoreyBakerRelativePermeability;

protected:
  virtual void
  PostProcessInput() override;

  array1d< real64 > m_phaseMinVolumeFraction;

  // water-oil data
  array1d< real64 > m_waterOilRelPermExponent;
  array1d< real64 > m_waterOilRelPermMaxValue;

  // gas-oil data
  array1d< real64 > m_gasOilRelPermExponent;
  array1d< real64 > m_gasOilRelPermMaxValue;

  real64 m_volFracScale;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyBakerRelativePermeabilityUpdate::Compute(
  arraySlice1d< real64 const > const & phaseVolFraction,
  arraySlice1d< real64 > const & phaseRelPerm,
  arraySlice2d< real64 > const & dPhaseRelPerm_dPhaseVolFrac ) const
{
  localIndex const NP = numPhases();

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      dPhaseRelPerm_dPhaseVolFrac[ip][jp] = 0.0;
    }
  }

  real64 const volFracScaleInv = 1.0 / m_volFracScale;
  integer const ip_water =
    m_phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
  integer const ip_oil = m_phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  integer const ip_gas = m_phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

  real64 oilRelPerm_wo = 0;               // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_wo_dOilVolFrac = 0;  // derivative w.r.t to So
  real64 oilRelPerm_go = 0;               // oil rel perm using two-phase gas-oil data
  real64 dOilRelPerm_go_dOilVolFrac = 0;  // derivative w.r.t to So

  // this function assumes that the oil phase can always be present (i.e., ip_oil > 0)

  // 1) Water and oil phase relative permeabilities using water-oil data
  if( ip_water >= 0 )
  {
    real64 const scaledWaterVolFrac =
      ( phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water] ) *
      volFracScaleInv;
    real64 const scaledOilVolFrac =
      ( phaseVolFraction[ip_oil] - m_phaseMinVolumeFraction[ip_oil] ) *
      volFracScaleInv;

    real64 const waterExponent =
      m_waterOilRelPermExponent[RelativePermeabilityBase::WaterOilPairPhaseType::WATER];
    real64 const waterMaxValue =
      m_waterOilRelPermMaxValue[RelativePermeabilityBase::WaterOilPairPhaseType::WATER];

    // water rel perm
    EvaluateBrooksCoreyFunction( scaledWaterVolFrac,
                                 volFracScaleInv,
                                 waterExponent,
                                 waterMaxValue,
                                 phaseRelPerm[ip_water],
                                 dPhaseRelPerm_dPhaseVolFrac[ip_water][ip_water] );

    real64 const oilExponent_wo =
      m_waterOilRelPermExponent[RelativePermeabilityBase::WaterOilPairPhaseType::OIL];
    real64 const oilMaxValue_wo =
      m_waterOilRelPermMaxValue[RelativePermeabilityBase::WaterOilPairPhaseType::OIL];

    // oil rel perm
    EvaluateBrooksCoreyFunction( scaledOilVolFrac,
                                 volFracScaleInv,
                                 oilExponent_wo,
                                 oilMaxValue_wo,
                                 oilRelPerm_wo,
                                 dOilRelPerm_wo_dOilVolFrac );
  }

  // 2) Gas and oil phase relative permeabilities using gas-oil data
  if( ip_gas >= 0 )
  {
    real64 const scaledGasVolFrac =
      ( phaseVolFraction[ip_gas] - m_phaseMinVolumeFraction[ip_gas] ) *
      volFracScaleInv;
    real64 const scaledOilVolFrac =
      ( phaseVolFraction[ip_oil] - m_phaseMinVolumeFraction[ip_oil] ) *
      volFracScaleInv;

    real64 const gasExponent =
      m_gasOilRelPermExponent[RelativePermeabilityBase::GasOilPairPhaseType::GAS];
    real64 const gasMaxValue =
      m_gasOilRelPermMaxValue[RelativePermeabilityBase::GasOilPairPhaseType::GAS];

    // gas rel perm
    EvaluateBrooksCoreyFunction( scaledGasVolFrac,
                                 volFracScaleInv,
                                 gasExponent,
                                 gasMaxValue,
                                 phaseRelPerm[ip_gas],
                                 dPhaseRelPerm_dPhaseVolFrac[ip_gas][ip_gas] );

    real64 const oilExponent_go =
      m_gasOilRelPermExponent[RelativePermeabilityBase::GasOilPairPhaseType::OIL];
    real64 const oilMaxValue_go =
      m_gasOilRelPermMaxValue[RelativePermeabilityBase::GasOilPairPhaseType::OIL];

    // oil rel perm
    EvaluateBrooksCoreyFunction( scaledOilVolFrac,
                                 volFracScaleInv,
                                 oilExponent_go,
                                 oilMaxValue_go,
                                 oilRelPerm_go,
                                 dOilRelPerm_go_dOilVolFrac );
  }

  // 3) Compute the "three-phase" oil relperm

  // if no gas, use water-oil data
  if( ip_gas < 0 )
  {
    phaseRelPerm[ip_oil] = oilRelPerm_wo;
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = dOilRelPerm_wo_dOilVolFrac;
  }
  // if no water, use gas-oil data
  else if( ip_water < 0 )
  {
    phaseRelPerm[ip_oil] = oilRelPerm_go;
    dPhaseRelPerm_dPhaseVolFrac[ip_oil][ip_oil] = dOilRelPerm_go_dOilVolFrac;
  }
  // if water and oil and gas can be present, use saturation-weighted interpolation
  else
  {
    real64 const shiftedWaterVolFrac =
      ( phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water] );

    InterpolateTwoPhaseRelPerms( shiftedWaterVolFrac,
                                 phaseVolFraction[ip_gas],
                                 m_phaseOrder,
                                 oilRelPerm_wo,
                                 dOilRelPerm_wo_dOilVolFrac,
                                 oilRelPerm_go,
                                 dOilRelPerm_go_dOilVolFrac,
                                 phaseRelPerm[ip_oil],
                                 dPhaseRelPerm_dPhaseVolFrac[ip_oil] );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyBakerRelativePermeabilityUpdate::EvaluateBrooksCoreyFunction(
  real64 const & scaledVolFrac,
  real64 const & dScaledVolFrac_dVolFrac,
  real64 const & exponent,
  real64 const & maxValue,
  real64 & relPerm,
  real64 & dRelPerm_dVolFrac )
{
  relPerm = 0.0;
  dRelPerm_dVolFrac = 0.0;

  if( scaledVolFrac > 0.0 && scaledVolFrac < 1.0 )
  {
    // intermediate value
    real64 const v = maxValue * std::pow( scaledVolFrac, exponent - 1.0 );

    relPerm = v * scaledVolFrac;
    dRelPerm_dVolFrac = v * exponent * dScaledVolFrac_dVolFrac;
  }
  else
  {
    relPerm = ( scaledVolFrac <= 0.0 ) ? 0.0 : maxValue;
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyBakerRelativePermeabilityUpdate::InterpolateTwoPhaseRelPerms(
  real64 const & shiftedWaterVolFrac,
  real64 const & gasVolFrac,
  arraySlice1d< integer const > const & phaseOrder,
  real64 const & relPerm_wo,
  real64 const & dRelPerm_wo_dOilVolFrac,
  real64 const & relPerm_go,
  real64 const & dRelPerm_go_dOilVolFrac,
  real64 & threePhaseRelPerm,
  arraySlice1d< real64 > const & dThreePhaseRelPerm_dVolFrac )
{
  integer const ip_water = phaseOrder[RelativePermeabilityBase::PhaseType::WATER];
  integer const ip_oil = phaseOrder[RelativePermeabilityBase::PhaseType::OIL];
  integer const ip_gas = phaseOrder[RelativePermeabilityBase::PhaseType::GAS];

  // if water phase is immobile, then use the two-phase gas-oil data only
  if( shiftedWaterVolFrac < NumericTraits< real64 >::eps )
  {
    threePhaseRelPerm = relPerm_go;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_go_dOilVolFrac;
  }
  // if gas phase is immobile, then use the two-phase water-oil data only
  else if( gasVolFrac < NumericTraits< real64 >::eps )
  {
    threePhaseRelPerm = relPerm_wo;
    dThreePhaseRelPerm_dVolFrac[ip_oil] = dRelPerm_wo_dOilVolFrac;
  }
  // if both the water phase and the gas phase are mobile,
  // then use a saturation-weighted interpolation of the two-phase oil rel perms
  else
  {
    real64 const sumRelPerm =
      ( shiftedWaterVolFrac * relPerm_wo + gasVolFrac * relPerm_go );
    real64 const dSumRelPerm_dWaterVolFrac = relPerm_wo;
    real64 const dSumRelPerm_dOilVolFrac =
      shiftedWaterVolFrac * dRelPerm_wo_dOilVolFrac +
      gasVolFrac * dRelPerm_go_dOilVolFrac;
    real64 const dSumRelPerm_dGasVolFrac = relPerm_go;

    real64 const sumVolFrac = shiftedWaterVolFrac + gasVolFrac;
    real64 const sumVolFracInv =
      1 / sumVolFrac;  // div by 0 handled by the if statement above
    real64 const dSumVolFracInv_dWaterVolFrac = -sumVolFracInv * sumVolFracInv;
    real64 const dSumVolFracInv_dGasVolFrac = dSumVolFracInv_dWaterVolFrac;

    threePhaseRelPerm = sumRelPerm * sumVolFracInv;  // three-phase oil rel perm
    dThreePhaseRelPerm_dVolFrac[ip_water] =
      dSumRelPerm_dWaterVolFrac * sumVolFracInv  // derivative w.r.t. Sw
      + sumRelPerm * dSumVolFracInv_dWaterVolFrac;
    dThreePhaseRelPerm_dVolFrac[ip_oil] =
      dSumRelPerm_dOilVolFrac * sumVolFracInv;  // derivative w.r.t. So
    dThreePhaseRelPerm_dVolFrac[ip_gas] =
      dSumRelPerm_dGasVolFrac * sumVolFracInv  // derivative w.r.t. Sg
      + sumRelPerm * dSumVolFracInv_dGasVolFrac;
  }
}

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_RELPERM_BROOKSCOREYBAKERRELATIVEPERMEABILITY_HPP
