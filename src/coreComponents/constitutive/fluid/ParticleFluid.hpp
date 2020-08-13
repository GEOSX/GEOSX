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
 * @file ParticleFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUID_HPP_

#include "constitutive/fluid/ParticleFluidBase.hpp"

namespace geosx
{
namespace constitutive
{
enum class ParticleSettlingModel
{
  Stokes,
  Intermediate,
  Turbulence
};

/**
 * @brief Kernel wrapper for ParticleFluid.
 */
class ParticleFluidUpdate final : public ParticleFluidBaseUpdate
{
public:
  /**
   * @brief Constructor.
   * @param particleSettlingModel
   * @param proppantDensity
   * @param fluidViscosity
   * @param proppantDiameter
   * @param hinderedSettlingCoefficient
   * @param collisionAlpha
   * @param slipConcentration
   * @param collisionBeta
   * @param sphericity
   * @param packPermeabilityCoef
   * @param isCollisionalSlip
   * @param maxProppantConcentration
   * @param settlingFactor
   * @param dSettlingFactor_dPressure
   * @param dSettlingFactor_dProppantConcentration
   * @param dSettlingFactor_dComponentConcentration
   * @param collisionFactor
   * @param dCollisionFactor_dProppantConcentration
   * @param proppantPackPermeability
   */
  ParticleFluidUpdate(
    ParticleSettlingModel const particleSettlingModel,
    real64 const proppantDensity,
    real64 const proppantDiameter,
    real64 const hinderedSettlingCoefficient,
    real64 const collisionAlpha,
    real64 const slipConcentration,
    real64 const collisionBeta,
    bool const isCollisionalSlip,
    real64 const maxProppantConcentration,
    arrayView1d< real64 > const & settlingFactor,
    arrayView1d< real64 > const & dSettlingFactor_dPressure,
    arrayView1d< real64 > const & dSettlingFactor_dProppantConcentration,
    arrayView2d< real64 > const & dSettlingFactor_dComponentConcentration,
    arrayView1d< real64 > const & collisionFactor,
    arrayView1d< real64 > const & dCollisionFactor_dProppantConcentration,
    arrayView1d< real64 > const & proppantPackPermeability ) :
    ParticleFluidBaseUpdate( isCollisionalSlip,
                             maxProppantConcentration,
                             settlingFactor,
                             dSettlingFactor_dPressure,
                             dSettlingFactor_dProppantConcentration,
                             dSettlingFactor_dComponentConcentration,
                             collisionFactor,
                             dCollisionFactor_dProppantConcentration,
                             proppantPackPermeability ),
    m_particleSettlingModel( particleSettlingModel ),
    m_proppantDensity( proppantDensity ),
    m_proppantDiameter( proppantDiameter ),
    m_hinderedSettlingCoefficient( hinderedSettlingCoefficient ),
    m_collisionAlpha( collisionAlpha ),
    m_slipConcentration( slipConcentration ),
    m_collisionBeta( collisionBeta )
  {}

  /**
   * @brief Copy constructor.
   */
  ParticleFluidUpdate( ParticleFluidUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  ParticleFluidUpdate( ParticleFluidUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  ParticleFluidUpdate &
  operator=( ParticleFluidUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  ParticleFluidUpdate &
  operator=( ParticleFluidUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Update(
    localIndex const k,
    real64 const proppantConcentration,
    real64 const fluidDensity,
    real64 const dFluidDensity_dPressure,
    arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
    real64 const fluidViscosity,
    real64 const dFluidViscosity_dPressure,
    arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration ) const override
  {
    Compute( proppantConcentration,
             fluidDensity,
             dFluidDensity_dPressure,
             dFluidDensity_dComponentConcentration,
             fluidViscosity,
             dFluidViscosity_dPressure,
             dFluidViscosity_dComponentConcentration,
             m_settlingFactor[k],
             m_dSettlingFactor_dPressure[k],
             m_dSettlingFactor_dProppantConcentration[k],
             m_dSettlingFactor_dComponentConcentration[k],
             m_collisionFactor[k],
             m_dCollisionFactor_dProppantConcentration[k] );
  }

private:
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void
  Compute(
    real64 const proppantConcentration,
    real64 const fluidDensity,
    real64 const dFluidDensity_dPressure,
    arraySlice1d< real64 const > const & dFluidDensity_dComponentConcentration,
    real64 const fluidViscosity,
    real64 const dFluidViscosity_dPressure,
    arraySlice1d< real64 const > const & dFluidViscosity_dComponentConcentration,
    real64 & settlingFactor,
    real64 & dSettlingFactor_dPressure,
    real64 & dSettlingFactor_dProppantConcentration,
    arraySlice1d< real64 > const & dSettlingFactor_dComponentConcentration,
    real64 & collisionFactor,
    real64 & dCollisionFactor_dProppantConcentration ) const;

  ParticleSettlingModel m_particleSettlingModel;

  real64 m_proppantDensity;

  real64 m_proppantDiameter;

  real64 m_hinderedSettlingCoefficient;

  real64 m_collisionAlpha;

  real64 m_slipConcentration;

  real64 m_collisionBeta;
};

class ParticleFluid : public ParticleFluidBase
{
public:
  ParticleFluid( std::string const & name, Group * const parent );

  virtual ~ParticleFluid() override;

  // *** ConstitutiveBase interface

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "ParticleFluid";
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

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ParticleFluidUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper() const;

  // *** Data repository keys

  struct viewKeyStruct : public ParticleFluidBase::viewKeyStruct
  {
    static constexpr auto fluidViscosityString = "fluidViscosity";
    static constexpr auto proppantDiameterString = "proppantDiameter";
    static constexpr auto proppantDensityString = "proppantDensity";
    static constexpr auto hinderedSettlingCoefficientString =
      "hinderedSettlingCoefficient";
    static constexpr auto collisionAlphaString = "collisionAlpha";
    static constexpr auto slipConcentrationString = "slipConcentration";
    static constexpr auto collisionBetaString = "collisionBeta";
    static constexpr auto bridgingFactorString = "bridgingFactor";
    static constexpr auto sphericityString = "sphericity";
    static constexpr auto particleSettlingModelString = "particleSettlingModel";

  } viewKeysParticleFluid;

protected:
  virtual void
  PostProcessInput() override;

private:
  string m_particleSettlingModelString;

  ParticleSettlingModel m_particleSettlingModel;

  real64 m_proppantDensity;

  real64 m_fluidViscosity;

  real64 m_proppantDiameter;

  real64 m_hinderedSettlingCoefficient;

  real64 m_collisionAlpha;

  real64 m_slipConcentration;

  real64 m_collisionBeta;

  real64 m_sphericity;

  real64 m_packPermeabilityCoef;
};

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void
ParticleFluidUpdate::Compute(
  real64 const proppantConcentration,
  real64 const fluidDensity,
  real64 const GEOSX_UNUSED_PARAM( dFluidDensity_dPressure ),
  arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM(
    dFluidDensity_dComponentConcentration ),
  real64 const fluidViscosity,
  real64 const GEOSX_UNUSED_PARAM( dFluidViscosity_dPressure ),
  arraySlice1d< real64 const > const & GEOSX_UNUSED_PARAM(
    dFluidViscosity_dComponentConcentration ),
  real64 & settlingFactor,
  real64 & dSettlingFactor_dPressure,
  real64 & dSettlingFactor_dProppantConcentration,
  arraySlice1d< real64 > const & dSettlingFactor_dComponentConcentration,
  real64 & collisionFactor,
  real64 & dCollisionFactor_dProppantConcentration ) const
{
  real64 const constCoef = 9.81 * m_proppantDiameter * m_proppantDiameter / 18.0;

  real64 singleParticleSettlingVelocity = 0.0;

  switch( m_particleSettlingModel )
  {
    case ParticleSettlingModel::Stokes:
    {
      singleParticleSettlingVelocity =
        constCoef * ( m_proppantDensity - fluidDensity ) / fluidViscosity;
      break;
    }
    case ParticleSettlingModel::Intermediate:
    {
      singleParticleSettlingVelocity = 0.2 * pow( m_proppantDiameter, 1.18 ) *
        pow( 9.81 * ( m_proppantDensity - fluidDensity ) / fluidDensity, 0.72 ) *
        pow( fluidDensity / fluidViscosity, 0.45 );
      break;
    }
    case ParticleSettlingModel::Turbulence:
    {
      singleParticleSettlingVelocity = 1.74 * pow( m_proppantDiameter, 0.5 ) *
        pow( 9.81 * ( m_proppantDensity - fluidDensity ) / fluidDensity, 0.5 );
      break;
    }
    default:
    {
    }
  }

  settlingFactor = 0.0;
  dSettlingFactor_dPressure = 0.0;
  dSettlingFactor_dProppantConcentration = 0.0;

  localIndex const NC = dSettlingFactor_dComponentConcentration.size();
  for( localIndex c = 0; c < NC; ++c )
  {
    dSettlingFactor_dComponentConcentration[c] = 0.0;
  }

  collisionFactor = 0.0;
  dCollisionFactor_dProppantConcentration = 0.0;

  if( proppantConcentration >= 0.0 &&
      proppantConcentration < m_maxProppantConcentration )
  {
    // settlingFactor
    settlingFactor = singleParticleSettlingVelocity *
      exp( -m_hinderedSettlingCoefficient * proppantConcentration );

    // collisionFactor
    // Collision model (We need to check the other models)
    if( m_isCollisionalSlip )
    {
      real64 const lambda = m_collisionAlpha -
        pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta );
      collisionFactor = ( lambda - 1.0 ) / ( 1.0 - proppantConcentration );

      // TODO: why? - Sergey
#if 0
      real64 dLambda_dC = -m_collisionBeta * pow( fabs( proppantConcentration - m_slipConcentration ), m_collisionBeta - 1.0 );
      if( proppantConcentration < m_slipConcentration )
      {
        dLambda_dC = -dLambda_dC;
      }
#endif
    }
  }
}

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_PARTICLEFLUID_HPP_ */
