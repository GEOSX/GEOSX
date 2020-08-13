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
 * @file ProppantSlurryFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_

#include "constitutive/fluid/SlurryFluidBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @brief Update class for the proppant slurry model suitable for lambda capture.
 */
class ProppantSlurryFluidUpdate final : public SlurryFluidBaseUpdate
{
public:
  /**
   * @brief
   * @param compressibility
   * @param referenceProppantDensity
   * @param referencePressure
   * @param referenceDensity
   * @param referenceViscosity
   * @param maxProppantConcentration
   * @param defaultDensity
   * @param defaultCompressibility
   * @param defaultViscosity
   * @param nIndices
   * @param Ks
   * @param isNewtonianFluid
   * @param density
   * @param dDens_dPres
   * @param dDens_dProppantConc
   * @param dDens_dCompConc
   * @param componentDensity
   * @param dCompDens_dPres
   * @param dCompDens_dCompConc
   * @param fluidDensity
   * @param dFluidDens_dPres
   * @param dFluidDens_dCompConc
   * @param fluidViscosity
   * @param dFluidVisc_dPres
   * @param dFluidVisc_dCompConc
   * @param viscosity
   * @param dVisc_dPres
   * @param dVisc_dProppantConc
   * @param dVisc_dCompConc
   */
  ProppantSlurryFluidUpdate( real64 const compressibility,
                             real64 const referenceProppantDensity,
                             real64 const referencePressure,
                             real64 const referenceDensity,
                             real64 const referenceViscosity,
                             real64 const maxProppantConcentration,
                             arrayView1d< real64 const > const & defaultDensity,
                             arrayView1d< real64 const > const & defaultCompressibility,
                             arrayView1d< real64 const > const & defaultViscosity,
                             arrayView1d< real64 const > const & nIndices,
                             arrayView1d< real64 const > const & Ks,
                             bool const isNewtonianFluid,
                             arrayView2d< real64 > const & density,
                             arrayView2d< real64 > const & dDens_dPres,
                             arrayView2d< real64 > const & dDens_dProppantConc,
                             arrayView3d< real64 > const & dDens_dCompConc,
                             arrayView3d< real64 > const & componentDensity,
                             arrayView3d< real64 > const & dCompDens_dPres,
                             arrayView4d< real64 > const & dCompDens_dCompConc,
                             arrayView2d< real64 > const & fluidDensity,
                             arrayView2d< real64 > const & dFluidDens_dPres,
                             arrayView3d< real64 > const & dFluidDens_dCompConc,
                             arrayView2d< real64 > const & fluidViscosity,
                             arrayView2d< real64 > const & dFluidVisc_dPres,
                             arrayView3d< real64 > const & dFluidVisc_dCompConc,
                             arrayView2d< real64 > const & viscosity,
                             arrayView2d< real64 > const & dVisc_dPres,
                             arrayView2d< real64 > const & dVisc_dProppantConc,
                             arrayView3d< real64 > const & dVisc_dCompConc ) :
    SlurryFluidBaseUpdate( defaultDensity,
                           defaultCompressibility,
                           defaultViscosity,
                           nIndices,
                           Ks,
                           isNewtonianFluid,
                           density,
                           dDens_dPres,
                           dDens_dProppantConc,
                           dDens_dCompConc,
                           componentDensity,
                           dCompDens_dPres,
                           dCompDens_dCompConc,
                           fluidDensity,
                           dFluidDens_dPres,
                           dFluidDens_dCompConc,
                           fluidViscosity,
                           dFluidVisc_dPres,
                           dFluidVisc_dCompConc,
                           viscosity,
                           dVisc_dPres,
                           dVisc_dProppantConc,
                           dVisc_dCompConc ),
    m_compressibility( compressibility ),
    m_referenceProppantDensity( referenceProppantDensity ),
    m_referencePressure( referencePressure ),
    m_referenceDensity( referenceDensity ),
    m_referenceViscosity( referenceViscosity ),
    m_maxProppantConcentration( maxProppantConcentration )
  {}

  /// Default copy constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate const & ) = default;

  /// Default move constructor
  ProppantSlurryFluidUpdate( ProppantSlurryFluidUpdate && ) = default;

  /// Deleted copy assignment operator
  ProppantSlurryFluidUpdate &
  operator=( ProppantSlurryFluidUpdate const & ) = delete;

  /// Deleted move assignment operator
  ProppantSlurryFluidUpdate &
  operator=( ProppantSlurryFluidUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          real64 const pressure,
          real64 const proppantConcentration,
          arraySlice1d< real64 const > const & componentConcentration,
          real64 const shearRate,
          integer const isProppantBoundary ) const override
  {
    UpdateFluidProperty( k, q, pressure, componentConcentration, shearRate );

    Compute( proppantConcentration,
             m_fluidDensity[k][q],
             m_dFluidDens_dPres[k][q],
             m_dFluidDens_dCompConc[k][q],
             m_fluidViscosity[k][q],
             m_dFluidVisc_dPres[k][q],
             m_dFluidVisc_dCompConc[k][q],
             isProppantBoundary,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_dDens_dProppantConc[k][q],
             m_dDens_dCompConc[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q],
             m_dVisc_dProppantConc[k][q],
             m_dVisc_dCompConc[k][q] );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  UpdateFluidProperty(
    localIndex const k,
    localIndex const q,
    real64 const pressure,
    arraySlice1d< real64 const > const & componentConcentration,
    real64 const shearRate ) const override
  {
    GEOSX_UNUSED_VAR( shearRate )
    ComputeFluidDensity( pressure,
                         componentConcentration,
                         m_componentDensity[k][q],
                         m_dCompDens_dPres[k][q],
                         m_dCompDens_dCompConc[k][q],
                         m_fluidDensity[k][q],
                         m_dFluidDens_dPres[k][q],
                         m_dFluidDens_dCompConc[k][q] );

    ComputeFluidViscosity( m_componentDensity[k][q],
                           m_dCompDens_dPres[k][q],
                           m_dCompDens_dCompConc[k][q],
                           m_fluidDensity[k][q],
                           m_dFluidDens_dPres[k][q],
                           m_dFluidDens_dCompConc[k][q],
                           m_fluidViscosity[k][q],
                           m_dFluidVisc_dPres[k][q],
                           m_dFluidVisc_dCompConc[k][q] );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  UpdateComponentDensity(
    localIndex const k,
    localIndex const q,
    real64 const pressure,
    arraySlice1d< real64 const > const & componentConcentration ) const override
  {
    ComputeComponentDensity( pressure,
                             componentConcentration,
                             m_componentDensity[k][q],
                             m_dCompDens_dPres[k][q],
                             m_dCompDens_dCompConc[k][q] );
  }

private:
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  ComputeFluidDensity(
    real64 const & pressure,
    arraySlice1d< real64 const > const & componentConcentration,
    arraySlice1d< real64 > const & componentDensity,
    arraySlice1d< real64 > const & dComponentDensity_dPressure,
    arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
    real64 & fluidDensity,
    real64 & dFluidDensity_dPressure,
    arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  ComputeComponentDensity(
    real64 const pressure,
    arraySlice1d< real64 const > const & componentConcentration,
    arraySlice1d< real64 > const & componentDensity,
    arraySlice1d< real64 > const & dComponentDensity_dPressure,
    arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  ComputeFluidViscosity(
    arraySlice1d< real64 const > const & componentDensity,
    arraySlice1d< real64 const > const & dComponentDensity_dPressure,
    arraySlice2d< real64 const > const & dComponentDensity_dComponentConcentration,
    real64 const fluidDensity,
    real64 const dFluidDensity_dPressure,
    arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
    real64 & fluidViscosity,
    real64 & dFluidViscosity_dPressure,
    arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  Compute(
    real64 const & proppantConcentration,
    real64 const & fluidDensity,
    real64 const & dFluidDensity_dPressure,
    arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
    real64 const & fluidViscosity,
    real64 const & dFluidViscosity_dPressure,
    arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
    integer const & isProppantBoundary,
    real64 & density,
    real64 & dDensity_dPressure,
    real64 & dDensity_dProppantConcentration,
    arraySlice1d< real64 > const & dDensity_dComponentConcentration,
    real64 & viscosity,
    real64 & dViscosity_dPressure,
    real64 & dViscosity_dProppantConcentration,
    arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const;

  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;
};

class ProppantSlurryFluid : public SlurryFluidBase
{
public:
  ProppantSlurryFluid( std::string const & name, Group * const parent );

  virtual ~ProppantSlurryFluid() override;

  // *** ConstitutiveBase interface

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "ProppantSlurryFluid";
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) override;

  using KernelWrapper = ProppantSlurryFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper();

  // *** Data repository keys

  struct viewKeyStruct : public SlurryFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString = "compressibility";
    static constexpr auto referencePressureString = "referencePressure";
    static constexpr auto referenceDensityString = "referenceDensity";
    static constexpr auto referenceProppantDensityString =
      "referenceProppantDensity";
    static constexpr auto maxProppantConcentrationString =
      "maxProppantConcentration";
    static constexpr auto referenceViscosityString = "referenceViscosity";

  } viewKeysProppantSlurryFluid;

protected:
  virtual void
  PostProcessInput() override;

private:
  real64 m_compressibility;

  real64 m_referenceProppantDensity;

  real64 m_referencePressure;

  real64 m_referenceDensity;

  real64 m_referenceViscosity;

  real64 m_maxProppantConcentration;
};

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::ComputeFluidDensity(
  real64 const & pressure,
  arraySlice1d< real64 const > const & componentConcentration,
  arraySlice1d< real64 > const & componentDensity,
  arraySlice1d< real64 > const & dComponentDensity_dPressure,
  arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration,
  real64 & fluidDensity,
  real64 & dFluidDensity_dPressure,
  arraySlice1d< real64 > const & dFluidDensity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 const baseFluidDensity = m_referenceDensity *
    exp( m_compressibility * ( pressure - m_referencePressure ) );
  real64 const dBaseFluidDensity_dPressure = m_compressibility * baseFluidDensity;

  fluidDensity = baseFluidDensity;
  dFluidDensity_dPressure = dBaseFluidDensity_dPressure;
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensity_dComponentConcentration[c] = 0.0;
  }

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultDensity[c] *
      exp( m_defaultCompressibility[c] * ( pressure - m_referencePressure ) );

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensity_dPressure[c] =
      m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensity_dComponentConcentration[c][i] = 0.0;
    }
    dComponentDensity_dComponentConcentration[c][c] = density;

    fluidDensity +=
      componentDensity[c] - componentConcentration[c] * baseFluidDensity;
    dFluidDensity_dPressure += dComponentDensity_dPressure[c] -
      componentConcentration[c] * dBaseFluidDensity_dPressure;
    dFluidDensity_dComponentConcentration[c] +=
      dComponentDensity_dComponentConcentration[c][c] - baseFluidDensity;
  }

  // TODO: why? - Sergey
  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidDensity_dComponentConcentration[c] = 0.0;
  }
}

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::ComputeComponentDensity(
  real64 const pressure,
  arraySlice1d< real64 const > const & componentConcentration,
  arraySlice1d< real64 > const & componentDensity,
  arraySlice1d< real64 > const & dComponentDensity_dPressure,
  arraySlice2d< real64 > const & dComponentDensity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  for( localIndex c = 0; c < NC; ++c )
  {
    real64 const density = m_defaultDensity[c] *
      exp( m_defaultCompressibility[c] * ( pressure - m_referencePressure ) );

    componentDensity[c] = componentConcentration[c] * density;
    dComponentDensity_dPressure[c] =
      m_defaultCompressibility[c] * componentDensity[c];

    for( localIndex i = 0; i < NC; ++i )
    {
      dComponentDensity_dComponentConcentration[c][i] = 0.0;
    }
    dComponentDensity_dComponentConcentration[c][c] = density;
  }
}

GEOSX_HOST_DEVICE
void
ProppantSlurryFluidUpdate::ComputeFluidViscosity(
  arraySlice1d< real64 const > const & componentDensity,
  arraySlice1d< real64 const > const & dComponentDensity_dPressure,
  arraySlice2d< real64 const > const & GEOSX_UNUSED_PARAM(
    dComponentDensity_dComponentConcentration ),
  real64 const fluidDensity,
  real64 const dFluidDensity_dPressure,
  arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM(
    dFluidDensity_dComponentConcentration ),
  real64 & fluidViscosity,
  real64 & dFluidViscosity_dPressure,
  arraySlice1d< real64 > const & dFluidViscosity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  fluidViscosity = m_referenceViscosity;
  dFluidViscosity_dPressure = 0.0;

  for( localIndex c = 0; c < NC; ++c )
  {
    dFluidViscosity_dComponentConcentration[c] = 0.0;
  }

  for( localIndex c1 = 0; c1 < NC; ++c1 )
  {
    fluidViscosity += componentDensity[c1] / fluidDensity *
      ( m_defaultViscosity[c1] - m_referenceViscosity );

    dFluidViscosity_dPressure += ( dComponentDensity_dPressure[c1] / fluidDensity -
                                   componentDensity[c1] / fluidDensity /
                                     fluidDensity * dFluidDensity_dPressure ) *
      ( m_defaultViscosity[c1] - m_referenceViscosity );
  }
}

void
ProppantSlurryFluidUpdate::Compute(
  real64 const & proppantConcentration,
  real64 const & fluidDensity,
  real64 const & dFluidDensity_dPressure,
  arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM(
    dFluidDensity_dComponentConcentration ),
  real64 const & fluidViscosity,
  real64 const & dFluidViscosity_dPressure,
  arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM(
    dFluidViscosity_dComponentConcentration ),
  integer const & isProppantBoundary,
  real64 & density,
  real64 & dDensity_dPressure,
  real64 & dDensity_dProppantConcentration,
  arraySlice1d< real64 > const & dDensity_dComponentConcentration,
  real64 & viscosity,
  real64 & dViscosity_dPressure,
  real64 & dViscosity_dProppantConcentration,
  arraySlice1d< real64 > const & dViscosity_dComponentConcentration ) const
{
  localIndex const NC = numFluidComponents();

  real64 effectiveConcentration = proppantConcentration;
  if( effectiveConcentration > 0.95 * m_maxProppantConcentration ||
      isProppantBoundary == 1 )
  {
    effectiveConcentration = 0.0;
  }

  density = ( 1.0 - effectiveConcentration ) * fluidDensity +
    effectiveConcentration * m_referenceProppantDensity;
  dDensity_dPressure = ( 1.0 - effectiveConcentration ) * dFluidDensity_dPressure;

  dDensity_dProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dDensity_dComponentConcentration[c] = 0.0;
  }

  real64 const coef =
    pow( 1.0 +
           1.25 * effectiveConcentration /
             ( 1.0 - effectiveConcentration / m_maxProppantConcentration ),
         2.0 );
  viscosity = fluidViscosity * coef;
  dViscosity_dPressure = dFluidViscosity_dPressure * coef;

  dViscosity_dProppantConcentration = 0.0;
  for( localIndex c = 0; c < NC; ++c )
  {
    dViscosity_dComponentConcentration[c] = 0.0;
  }
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_PROPPANTSLURRYFLUID_HPP_ */
