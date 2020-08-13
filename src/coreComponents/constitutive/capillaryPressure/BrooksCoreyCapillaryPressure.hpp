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
 * @file BrooksCoreyCapillaryPressure.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
#define GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP

#include "constitutive/capillaryPressure/CapillaryPressureBase.hpp"

namespace geosx
{
namespace constitutive
{
class BrooksCoreyCapillaryPressureUpdate final : public CapillaryPressureBaseUpdate
{
public:
  BrooksCoreyCapillaryPressureUpdate(
    arrayView1d< real64 const > const & phaseMinVolumeFraction,
    arrayView1d< real64 const > const & phaseCapPressureExponentInv,
    arrayView1d< real64 const > const & phaseEntryPressure,
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
    m_phaseEntryPressure( phaseEntryPressure ),
    m_capPressureEpsilon( capPressureEpsilon ),
    m_volFracScale( volFracScale )
  {}

  /// Default copy constructor
  BrooksCoreyCapillaryPressureUpdate( BrooksCoreyCapillaryPressureUpdate const & ) =
    default;

  /// Default move constructor
  BrooksCoreyCapillaryPressureUpdate( BrooksCoreyCapillaryPressureUpdate && ) = default;

  /// Deleted copy assignment operator
  BrooksCoreyCapillaryPressureUpdate &
  operator=(
    BrooksCoreyCapillaryPressureUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrooksCoreyCapillaryPressureUpdate &
  operator=(
    BrooksCoreyCapillaryPressureUpdate && ) = delete;

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
  EvaluateBrooksCoreyFunction(
    real64 const scaledWettingVolFrac,
    real64 const dScaledWettingPhaseVolFrac_dVolFrac,
    real64 const exponentInv,
    real64 const entryPressure,
    real64 const eps,
    real64 & phaseCapPressure,
    real64 & dPhaseCapPressure_dVolFrac );

  arrayView1d< real64 const > m_phaseMinVolumeFraction;
  arrayView1d< real64 const > m_phaseCapPressureExponentInv;
  arrayView1d< real64 const > m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};

class BrooksCoreyCapillaryPressure : public CapillaryPressureBase
{
public:
  BrooksCoreyCapillaryPressure( std::string const & name,
                                dataRepository::Group * const parent );

  virtual ~BrooksCoreyCapillaryPressure() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "BrooksCoreyCapillaryPressure";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyCapillaryPressureUpdate;

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
    static constexpr auto phaseEntryPressureString = "phaseEntryPressure";
    static constexpr auto capPressureEpsilonString = "capPressureEpsilon";
  } viewKeysBrooksCoreyCapillaryPressure;

protected:
  virtual void
  PostProcessInput() override;

  array1d< real64 > m_phaseMinVolumeFraction;
  array1d< real64 > m_phaseCapPressureExponentInv;
  array1d< real64 > m_phaseEntryPressure;

  real64 m_capPressureEpsilon;
  real64 m_volFracScale;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyCapillaryPressureUpdate::Compute(
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

  real64 const volFracScaleInv = 1.0 / m_volFracScale;

  // the Brooks-Corey model does not support volFracScaled = 0,
  // hence we need an epsilon value to avoid a division by zero
  // TODO: for S < epsilon, replace the original unbounded BC curve with a bounded power-law extension
  real64 const eps = m_capPressureEpsilon;

  // compute first water-oil capillary pressure as a function of water-phase vol fraction
  integer const ip_water = m_phaseOrder[CapillaryPressureBase::PhaseType::WATER];
  if( ip_water >= 0 )
  {
    real64 const volFracScaled =
      ( phaseVolFraction[ip_water] - m_phaseMinVolumeFraction[ip_water] ) *
      volFracScaleInv;
    real64 const exponentInv = m_phaseCapPressureExponentInv[ip_water];
    real64 const entryPressure = m_phaseEntryPressure[ip_water];

    real64 const wettingVolFracScaled = volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac = volFracScaleInv;

    EvaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 exponentInv,
                                 entryPressure,
                                 eps,
                                 phaseCapPres[ip_water],
                                 dPhaseCapPres_dPhaseVolFrac[ip_water][ip_water] );
  }

  // compute first gas-oil capillary pressure as a function of gas-phase vol fraction
  integer const ip_gas = m_phaseOrder[CapillaryPressureBase::PhaseType::GAS];
  if( ip_gas >= 0 )
  {
    real64 const volFracScaled =
      ( phaseVolFraction[ip_gas] - m_phaseMinVolumeFraction[ip_gas] ) *
      volFracScaleInv;
    real64 const exponentInv = m_phaseCapPressureExponentInv[ip_gas];
    real64 const entryPressure =
      -m_phaseEntryPressure[ip_gas];  // for gas capillary pressure, take the opposite of the
                                      // VG function

    real64 const wettingVolFracScaled = 1 - volFracScaled;
    real64 const dWettingVolFracScaled_dVolFrac = -volFracScaleInv;

    EvaluateBrooksCoreyFunction( wettingVolFracScaled,
                                 dWettingVolFracScaled_dVolFrac,
                                 exponentInv,
                                 entryPressure,
                                 eps,
                                 phaseCapPres[ip_gas],
                                 dPhaseCapPres_dPhaseVolFrac[ip_gas][ip_gas] );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyCapillaryPressureUpdate::EvaluateBrooksCoreyFunction(
  real64 const scaledWettingVolFrac,
  real64 const dScaledWettingPhaseVolFrac_dVolFrac,
  real64 const exponentInv,
  real64 const entryPressure,
  real64 const eps,
  real64 & phaseCapPressure,
  real64 & dPhaseCapPressure_dVolFrac )
{
  real64 const exponent =
    1.0 / exponentInv;  // div by 0 taken care of by initialization check

  phaseCapPressure = 0.0;
  dPhaseCapPressure_dVolFrac = 0.0;

  if( scaledWettingVolFrac >= eps && scaledWettingVolFrac < 1.0 )
  {
    // intermediate value
    real64 const val = entryPressure / pow( scaledWettingVolFrac, exponent + 1 );

    phaseCapPressure =
      val * scaledWettingVolFrac;  // entryPressure * (S_w)^( - 1 / exponentInv )
    dPhaseCapPressure_dVolFrac =
      -dScaledWettingPhaseVolFrac_dVolFrac * val * exponent;
  }
  else  // enforce a constant and bounded capillary pressure
  {
    phaseCapPressure = ( scaledWettingVolFrac < eps ) ? entryPressure /
        pow( eps, exponent )  // div by 0 taken care of by initialization check
                                                      : entryPressure;
  }
}

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_CAPILLARYPRESSURE_BROOKSCOREYCAPILLARYPRESSURE_HPP
