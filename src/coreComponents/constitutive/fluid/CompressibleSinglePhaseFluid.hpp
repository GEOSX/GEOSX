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
 * @file CompressibleSinglePhaseFluid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_

#include "constitutive/fluid/SingleFluidBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @brief Update class for the model suitable for lambda capture.
 * @tparam DENS_EAT type of density exponent approximation
 * @tparam VISC_EAT type of viscosity exponent approximation
 */
template< ExponentApproximationType DENS_EAT, ExponentApproximationType VISC_EAT >
class CompressibleSinglePhaseUpdate final : public SingleFluidBaseUpdate
{
public:
  using DensRelationType = ExponentialRelation< real64, DENS_EAT >;
  using ViscRelationType = ExponentialRelation< real64, VISC_EAT >;

  CompressibleSinglePhaseUpdate( DensRelationType const & densRelation,
                                 ViscRelationType const & viscRelation,
                                 arrayView2d< real64 > const & density,
                                 arrayView2d< real64 > const & dDens_dPres,
                                 arrayView2d< real64 > const & viscosity,
                                 arrayView2d< real64 > const & dVisc_dPres ) :
    SingleFluidBaseUpdate( density, dDens_dPres, viscosity, dVisc_dPres ),
    m_densRelation( densRelation ),
    m_viscRelation( viscRelation )
  {}

  /// Default copy constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate const & ) = default;

  /// Default move constructor
  CompressibleSinglePhaseUpdate( CompressibleSinglePhaseUpdate && ) = default;

  /// Deleted copy assignment operator
  CompressibleSinglePhaseUpdate &
  operator=(
    CompressibleSinglePhaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  CompressibleSinglePhaseUpdate &
  operator=( CompressibleSinglePhaseUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Compute( real64 const pressure,
           real64 & density,
           real64 & viscosity ) const override
  {
    m_densRelation.Compute( pressure, density );
    m_viscRelation.Compute( pressure, viscosity );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Compute( real64 const pressure,
           real64 & density,
           real64 & dDensity_dPressure,
           real64 & viscosity,
           real64 & dViscosity_dPressure ) const override
  {
    m_densRelation.Compute( pressure, density, dDensity_dPressure );
    m_viscRelation.Compute( pressure, viscosity, dViscosity_dPressure );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          real64 const pressure ) const override
  {
    Compute( pressure,
             m_density[k][q],
             m_dDens_dPres[k][q],
             m_viscosity[k][q],
             m_dVisc_dPres[k][q] );
  }

private:
  DensRelationType m_densRelation;
  ViscRelationType m_viscRelation;
};

class CompressibleSinglePhaseFluid : public SingleFluidBase
{
public:
  CompressibleSinglePhaseFluid( std::string const & name, Group * const parent );

  virtual ~CompressibleSinglePhaseFluid() override;

  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const override;

  static std::string
  CatalogName()
  {
    return "CompressibleSinglePhaseFluid";
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

  /// Type of kernel wrapper for in-kernel update (TODO: support multiple EAT, not just linear)
  using KernelWrapper =
    CompressibleSinglePhaseUpdate< ExponentApproximationType::Linear,
                                   ExponentApproximationType::Linear >;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper
  createKernelWrapper();

  struct viewKeyStruct : public SingleFluidBase::viewKeyStruct
  {
    static constexpr auto compressibilityString = "compressibility";
    static constexpr auto viscosibilityString = "viscosibility";
    static constexpr auto referencePressureString = "referencePressure";
    static constexpr auto referenceDensityString = "referenceDensity";
    static constexpr auto referenceViscosityString = "referenceViscosity";
    static constexpr auto densityModelString = "densityModel";
    static constexpr auto viscosityModelString = "viscosityModel";
  };

protected:
  virtual void
  PostProcessInput() override;

private:
  /// scalar fluid bulk modulus parameter
  real64 m_compressibility;

  /// scalar fluid viscosity exponential coefficient
  real64 m_viscosibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// reference density parameter
  real64 m_referenceDensity;

  /// reference viscosity parameter
  real64 m_referenceViscosity;

  /// input string for type of density model (linear, quadratic, exponential)
  string m_densityModelString;

  /// input string for type of viscosity model (linear, quadratic, exponential)
  string m_viscosityModelString;

  /// type of density model (linear, quadratic, exponential)
  ExponentApproximationType m_densityModelType;

  /// type of viscosity model (linear, quadratic, exponential)
  ExponentApproximationType m_viscosityModelType;
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_FLUID_COMPRESSIBLESINGLEPHASEFLUID_HPP_ */
