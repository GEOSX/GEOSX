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
 * @file BrooksCoreyRelativePermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
#define GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP

#include "constitutive/relativePermeability/RelativePermeabilityBase.hpp"

namespace geosx
{
namespace constitutive
{
class BrooksCoreyRelativePermeabilityUpdate final
  : public RelativePermeabilityBaseUpdate
{
public:
  BrooksCoreyRelativePermeabilityUpdate(
    arrayView1d< real64 const > const & phaseMinVolumeFraction,
    arrayView1d< real64 const > const & phaseRelPermExponent,
    arrayView1d< real64 const > const & phaseRelPermMaxValue,
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
    m_phaseRelPermExponent( phaseRelPermExponent ),
    m_phaseRelPermMaxValue( phaseRelPermMaxValue ),
    m_volFracScale( volFracScale )
  {}

  /// Default copy constructor
  BrooksCoreyRelativePermeabilityUpdate(
    BrooksCoreyRelativePermeabilityUpdate const & ) = default;

  /// Default move constructor
  BrooksCoreyRelativePermeabilityUpdate(
    BrooksCoreyRelativePermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  BrooksCoreyRelativePermeabilityUpdate &
  operator=(
    BrooksCoreyRelativePermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  BrooksCoreyRelativePermeabilityUpdate &
  operator=(
    BrooksCoreyRelativePermeabilityUpdate && ) = delete;

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
  arrayView1d< real64 const > m_phaseMinVolumeFraction;
  arrayView1d< real64 const > m_phaseRelPermExponent;
  arrayView1d< real64 const > m_phaseRelPermMaxValue;
  real64 m_volFracScale;
};

class BrooksCoreyRelativePermeability : public RelativePermeabilityBase
{
public:
  BrooksCoreyRelativePermeability( std::string const & name,
                                   dataRepository::Group * const parent );

  virtual ~BrooksCoreyRelativePermeability() override;

  void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  //START_SPHINX_INCLUDE_00
  static std::string
  CatalogName()
  {
    return "BrooksCoreyRelativePermeability";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = BrooksCoreyRelativePermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper();

  //START_SPHINX_INCLUDE_01
  struct viewKeyStruct : RelativePermeabilityBase::viewKeyStruct
  {
    static constexpr auto phaseMinVolumeFractionString =
      "phaseMinVolumeFraction";
    static constexpr auto phaseRelPermExponentString = "phaseRelPermExponent";
    static constexpr auto phaseRelPermMaxValueString = "phaseRelPermMaxValue";

    using ViewKey = dataRepository::ViewKey;

    ViewKey phaseMinVolumeFraction = { phaseMinVolumeFractionString };
    ViewKey phaseRelPermExponent = { phaseRelPermExponentString };
    ViewKey phaseRelPermMaxValue = { phaseRelPermMaxValueString };

  } vieKeysBrooksCoreyRelativePermeability;

protected:
  virtual void
  PostProcessInput() override;

  //START_SPHINX_INCLUDE_02
  array1d< real64 > m_phaseMinVolumeFraction;
  array1d< real64 > m_phaseRelPermExponent;
  array1d< real64 > m_phaseRelPermMaxValue;

  real64 m_volFracScale;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
BrooksCoreyRelativePermeabilityUpdate::Compute(
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
  real64 const satScaleInv = 1.0 / m_volFracScale;

  for( localIndex ip = 0; ip < NP; ++ip )
  {
    real64 const satScaled =
      ( phaseVolFraction[ip] - m_phaseMinVolumeFraction[ip] ) * satScaleInv;
    real64 const exponent = m_phaseRelPermExponent[ip];
    real64 const scale = m_phaseRelPermMaxValue[ip];

    if( satScaled > 0.0 && satScaled < 1.0 )
    {
      // intermediate value
      real64 const v = scale * pow( satScaled, exponent - 1.0 );

      phaseRelPerm[ip] = v * satScaled;
      dPhaseRelPerm_dPhaseVolFrac[ip][ip] = v * exponent * satScaleInv;
    }
    else
    {
      phaseRelPerm[ip] = ( satScaled <= 0.0 ) ? 0.0 : scale;
    }
  }
}

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_BROOKSCOREYRELATIVEPERMEABILITY_HPP
