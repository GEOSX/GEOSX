/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */


/**
 * @file ThermoPoroElastic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_THERMOPOROELASTIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_THERMOPOROELASTIC_HPP_

//#include "SolidBase.hpp"
#include "PoroElastic.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 * This class provides constitutive update routines for a ThermoPoroElastic solid.
 * It is derived from a solid model (the template paramer UPDATE_BASE)
 * allowing us to mix-and-match solid models.
 *
 * @tparam UPDATE_BASE A standard solid model providing stress-strain behavior
 */
template< typename UPDATE_BASE >
class ThermoPoroElasticUpdates : public UPDATE_BASE
{
public:

  /**
   * @brief Constructor
   * @tparam PARAMS A variadic list of template parameters
   * @param inputThermalStressCoefficient The thermal stress coefficient
   * @param baseParams Constructor parameters passed from base class
   */
  template< typename ... PARAMS >
  ThermoPoroElasticUpdates( real64 const & inputBiotCoefficient,
                            real64 const & inputThermalStressCoefficient,
                            real64 const & inputThermalPorosityCoefficient,
                            PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_biotCoefficient( inputBiotCoefficient ),
    m_thermalStressCoefficient( inputThermalStressCoefficient ),
    m_thermalPorosityCoefficient( inputThermalPorosityCoefficient )
  {}

  using UPDATE_BASE::getElasticStiffness;
  using UPDATE_BASE::smallStrainNoStateUpdate;
  using UPDATE_BASE::smallStrainUpdate;
  using UPDATE_BASE::smallStrainNoStateUpdate_StressOnly;
  using UPDATE_BASE::smallStrainUpdate_StressOnly;
  using UPDATE_BASE::hypoUpdate;
  using UPDATE_BASE::hyperUpdate;
  using UPDATE_BASE::saveConvergedState;

  /**
   * @brief Get Biot coefficient
   * @return Biot coefficient
   */
  GEOSX_HOST_DEVICE
  real64 getBiotCoefficient() const
  {
    return m_biotCoefficient;
  }

  /**
   * @brief Get thermal stress coefficient
   * @return thermal stress coefficient
   */
  GEOSX_HOST_DEVICE
  real64 getThermalStressCoefficient() const
  {
    return m_thermalStressCoefficient;
  }

  /**
   * @brief Get thermal porosity coefficient
   * @return thermal porosity coefficient
   */
  GEOSX_HOST_DEVICE
  real64 getThermalPorosityCoefficient() const
  {
    return m_thermalPorosityCoefficient;
  }

private:
  real64 m_biotCoefficient; ///< Scalar Biot coefficient
  real64 m_thermalStressCoefficient; ///< Scalar thermal stress coefficient
  real64 m_thermalPorosityCoefficient; ///< Scalar thermal porosity coefficient

};


/**
 * @brief (Empty) ThermoPoroElastic base class deriving from PoroElasticBase
 */
class ThermoPoroElasticBase : public PoroElasticBase
{};

/**
 * @brief Class to represent a ThermoPoroElastic solid
 * @tparam BASE A standard solid model from which this class is derived
 */
template< typename BASE >
class ThermoPoroElastic : public PoroElastic< BASE >
{
public:
  /// Alias for ThermoPoroElasticUpdates
  using KernelWrapper = ThermoPoroElasticUpdates< typename BASE::KernelWrapper >;

  using PoroElastic< BASE >::m_biotCoefficient;
  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  ThermoPoroElastic( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~ThermoPoroElastic() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "ThermoPoro" ) + BASE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * Create a kernel update object and return it
   * @return Kernel update object
   */
  KernelWrapper createKernelUpdates()
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_biotCoefficient,
                                                                       m_thermalStressCoefficient,
                                                                       m_thermalPorosityCoefficient );
  }

  /// Data view keys
  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * thermalStressCoefficientString() { return "thermalStressCoefficient"; }
    static constexpr char const * thermalPorosityCoefficientString() { return "thermalPorosityCoefficient"; }
  };


protected:

  /// scalar thermal stress coefficient
  real64 m_thermalStressCoefficient;

  /// scalar thermal porosity coefficient
  real64 m_thermalPorosityCoefficient;
};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_THERMOPOROELASTIC_HPP_ */
