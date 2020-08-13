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
 * @file VanGenuchtenCapillaryPressure.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_VANGENUCHTENCAPILLARYPRESSURE_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_VANGENUCHTENCAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

namespace geosx
{
namespace constitutive
{
class VanGenuchtenCapillaryPressureUpdate final
  : public CapillaryPressureBaseUpdate
{
public:
  VanGenuchtenCapillaryPressureUpdate(
    arrayView1d< real64 const > const & phaseMinVolumeFraction,
    arrayView1d< real64 const > const & phaseCapPressureExponentInv,
    arrayView1d< real64 const > const & phaseCapPressureMultiplier,
    real64 const capPressureEpsilon,
    real64 const volFracScale,
    arrayView1d< integer const > const & phaseTypes,
    arrayView1d< integer const > const & phaseOrder,
    arrayView3d< real64 > const & phaseCapPressure,
    arrayView4d< real64 > const & dPhaseCapPressure_dPhaseVolFrac ) :
    CapillaryPressureBaseUpdate( phaseTypes,
                                 phaseOrder,
                                 phaseCapPressure,
                                 dPhaseCapPressure_dPhaseVolFrac ),
    m_phaseMinVolumeFraction( phaseMinVolumeFraction ),
    m_phaseCapPressureExponentInv( phaseCapPressureExponentInv ),
    m_phaseCapPressureMultiplier( phaseCapPressureMultiplier ),
    m_capPressureEpsilon( capPressureEpsilon ),
    m_volFracScale( volFracScale )
  {}

  /// Default copy constructor
  VanGenuchtenCapillaryPressureUpdate(
    VanGenuchtenCapillaryPressureUpdate const & ) = default;

  /// Default move constructor
  VanGenuchtenCapillaryPressureUpdate( VanGenuchtenCapillaryPressureUpdate && ) = default;

  /// Deleted copy assignment operator
  VanGenuchtenCapillaryPressureUpdate &
  operator=(
    VanGenuchtenCapillaryPressureUpdate const & ) = delete;

  /// Deleted move assignment operator
  VanGenuchtenCapillaryPressureUpdate &
  operator=(
    VanGenuchtenCapillaryPressureUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Compute(
    arraySlice1d< real64 const > const & phaseVolFraction,
    arraySlice1d< real64 > const & phaseCapPres,
    arraySlice2d< real64 > const & dPhaseCapPres_dPhaseVolFrac ) const override;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          arraySlice1d< real64 const > const & phaseVolFraction ) const override
  {
    Compute( phaseVolFraction,
             m_phaseCapPressure[k][q],
             m_dPhaseCapPressure_dPhaseVolFrac[k][q] );
  }

private:
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  EvaluateVanGenuchtenFunction(
    real64 const scaledWettingVolFrac,
    real64 const dScaledWettingPhaseVolFrac_dVolFrac,
    real64 const exponentInv,
    real64 const multiplier,
    real64 const eps,
    real64 & phaseCapPressure,
    real64 & dPhaseCapPressure_dVolFrac );

  arrayView1d< real64 const > m_phaseMinVolumeFraction;
  arrayView1d< real64 const > m_phaseCapPressureExponentInv;
  arrayView1d< real64 const > m_phaseCapPressureMultiplier;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};

class VanGenuchtenCapillaryPressure : public CapillaryPressureBase
{
public:
  VanGenuchtenCapillaryPressure( std::string const & name,
                                 dataRepository::Group * const parent );

  virtual ~VanGenuchtenCapillaryPressure() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "VanGenuchtenCapillaryPressure";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = VanGenuchtenCapillaryPressureUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper();

  struct viewKeyStruct : CapillaryPressureBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString =
      "phaseMinVolumeFraction";
    static constexpr auto phaseCapPressureExponentInvString =
      "phaseCapPressureExponentInv";
    static constexpr auto phaseCapPressureMultiplierString =
      "phaseCapPressureMultiplier";
    static constexpr auto capPressureEpsilonString = "capPressureEpsilon";
  } viewKeysVanGenuchtenCapillaryPressure;

protected:
  virtual void
  PostProcessInput() override;

  array1d< real64 > m_phaseMinVolumeFraction;
  array1d< real64 > m_phaseCapPressureExponentInv;
  array1d< real64 > m_phaseCapPressureMultiplier;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
VanGenuchtenCapillaryPressureUpdate::Compute(
  arraySlice1d< real64 const > const & phaseVolFraction,
  arraySlice1d< real64 > const & phaseCapPres,
  arraySlice2d< real64 > const & dPhaseCapPres_dPhaseVolFrac ) const
{
  localIndex const NP = numPhases();

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    for( localIndex jp = 0; jp < NP; ++jp )
    {
      dPhaseCapPres_dPhaseVolFrac[ip][jp] = 0.0;
    }
  }

  // the VanGenuchten model does not support volFracScaled = 0 and = 1
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon and S > 1 - epsilon, replace the original unbounded VG curve with a bounded power-law
  // extension
  real64 const eps = m_capPressureEpsilon;
  real64 const volFracScaleInv = 1.0 / m_volFracScale;

  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = m_phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  if( ip_water >= 0 )
  {
    real64 const volFracScaled =
      ( phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water] ) *
      volFracScaleInv;
    real64 const exponentInv =
      m_phaseCapPressureExponentInv[ip_water];  // div by 0 taken care of by initialization
                                                // check
    real64 const multiplier = m_phaseCapPressureMultiplier[ip_water];

    real64 const scaledWettingVolFrac = volFracScaled;
    real64 const dScaledWettingPhaseVolFrac_dVolFrac = volFracScaleInv;

    EvaluateVanGenuchtenFunction( scaledWettingVolFrac,
                                  dScaledWettingPhaseVolFrac_dVolFrac,
                                  exponentInv,
                                  multiplier,
                                  eps,
                                  phaseCapPres[ip_water],
                                  dPhaseCapPres_dPhaseVolFrac[ip_water][ip_water] );
  }

  // then compute the oil-gas capillary pressure as a function of gas-phase vol fraction
  integer const ip_gas = m_phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  if( ip_gas >= 0 )
  {
    real64 const volFracScaled =
      ( phaseVolFraction[ip_gas] - m_phaseMinVolumeFraction[ip_gas] ) *
      volFracScaleInv;
    real64 const exponentInv =
      m_phaseCapPressureExponentInv[ip_gas];  // div by 0 taken care of by initialization
                                              // check
    real64 const multiplier =
      -m_phaseCapPressureMultiplier[ip_gas];  // for gas capillary pressure, take the opposite
                                              // of the VG function

    real64 const scaledWettingVolFrac = 1 - volFracScaled;
    real64 const dScaledWettingPhaseVolFrac_dVolFrac = -volFracScaleInv;

    EvaluateVanGenuchtenFunction( scaledWettingVolFrac,
                                  dScaledWettingPhaseVolFrac_dVolFrac,
                                  exponentInv,
                                  multiplier,
                                  eps,
                                  phaseCapPres[ip_gas],
                                  dPhaseCapPres_dPhaseVolFrac[ip_gas][ip_gas] );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
VanGenuchtenCapillaryPressureUpdate::EvaluateVanGenuchtenFunction(
  real64 const scaledWettingVolFrac,
  real64 const dScaledWettingPhaseVolFrac_dVolFrac,
  real64 const exponentInv,
  real64 const multiplier,
  real64 const eps,
  real64 & phaseCapPressure,
  real64 & dPhaseCapPressure_dVolFrac )
{
  real64 const exponent =
    1.0 / exponentInv;  // div by 0 taken care of by initialization check

  phaseCapPressure = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;

  if( scaledWettingVolFrac >= eps && scaledWettingVolFrac < 1.0 - eps )
  {
    // intermediate value
    real64 const a = 1 / pow( scaledWettingVolFrac, exponent + 1 );
    real64 const b = multiplier *
      pow( a * scaledWettingVolFrac - 1, 0.5 * ( 1 - exponentInv ) - 1 );

    phaseCapPressure = b *
      ( a * scaledWettingVolFrac - 1 );  // multiplier * ( S_w^(-1/m) - 1 )^( (1-m)/2 )
    dPhaseCapPressure_dVolFrac = -dScaledWettingPhaseVolFrac_dVolFrac * 0.5 *
      exponent * ( 1 - exponentInv ) * a * b;
  }
  else  // enforce a constant and bounded capillary pressure
  {
    phaseCapPressure = ( scaledWettingVolFrac <
                         eps )  // div by 0 taken care of by initialization check
      ? multiplier * pow( 1 / pow( eps, exponent ) - 1, 0.5 * ( 1 - exponentInv ) )
      : multiplier * pow( 1 / pow( 1 - eps, exponent ) - 1, 0.5 * ( 1 - exponentInv ) );
  }
}

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_VANGENUCHTENCAPILLARYPRESSURE_HPP
