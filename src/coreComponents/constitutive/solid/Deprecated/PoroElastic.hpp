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
 * @file PoroElastic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_

#include "SolidBase.hpp"
#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 * This class provides constitutive update routines for a porous solid.
 * It is derived from a solid model (the template paramer UPDATE_BASE)
 * allowing us to mix-and-match solid models.
 *
 * @tparam UPDATE_BASE A standard solid model providing stress-strain behavior
 */
template< typename UPDATE_BASE >
class PoroElasticUpdates : public UPDATE_BASE
{
public:

  /**
   * @brief Constructor
   * @tparam PARAMS A variadic list of template parameters
   * @param inputBiotCoefficient The Biot coefficient
   * @param baseParams Constructor parameters passed from base class
   */
  template< typename ... PARAMS >
  PoroElasticUpdates( real64 const & inputBiotCoefficient,
                      PARAMS && ... baseParams ):
    UPDATE_BASE( std::forward< PARAMS >( baseParams )... ),
    m_biotCoefficient( inputBiotCoefficient )
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

private:
  real64 m_biotCoefficient; ///< Scalar Biot coefficient

};


/**
 * @brief (Empty) porous base class deriving from solid base
 */
class PoroElasticBase : public SolidBase
{};

/**
 * @brief Class to represent a porous solid
 * @tparam BASE A standard solid model from which this class is derived
 */
template< typename BASE >
class PoroElastic : public BASE
{
public:

  /// Alias for PoroElasticUpdates
  using KernelWrapper = PoroElasticUpdates< typename BASE::KernelWrapper >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PoroElastic( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PoroElastic() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Poro" ) + BASE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /// Post-process XML input
  virtual void postProcessInput() override;

  /**
   * @brief Deliver a clone of this object
   * @return Pointer to clone
   * @param name Object name
   * @param parent Object's parent group
   */
  std::unique_ptr< ConstitutiveBase >
  deliverClone( string const & name,
                dataRepository::Group * const parent ) const override;

  /**
   * @brief Allocate constitutive arrays
   * @param parent Object's parent group (an element region)
   * @param numConstitutivePointsPerParentIndex (number of quadrature points per element)
   */
  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  /**
   * @brief Perform pore volume updates point-wise
   * @param[in] pres Current pressure
   * @param[in] k Element index
   * @param[in] q Quadrature index
   */
  inline virtual void
  stateUpdatePointPressure( real64 const & pres,
                            localIndex const k,
                            localIndex const q ) override
  {
    m_poreVolumeRelation.compute( pres, m_poreVolumeMultiplier[k][q], m_dPVMult_dPressure[k][q] );
  }

  /**
   * @brief Perform pore volume updates in batch model
   * @param pres Pressure array view
   * @param dPres Pressure derivative array view
   */
  virtual void stateUpdateBatchPressure( arrayView1d< real64 const > const & pres,
                                         arrayView1d< real64 const > const & dPres ) override final;

  /**
   * Create a kernel update object and return it
   * @return Kernel update object
   */
  KernelWrapper createKernelUpdates()
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_biotCoefficient );
  }

  /// Data view keys
  struct viewKeyStruct : public BASE::viewKeyStruct
  {
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * biotCoefficientString() { return "BiotCoefficient"; }
  };


protected:
  /// scalar compressibility parameter
  real64 m_compressibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  /// scalar Biot's coefficient
  real64 m_biotCoefficient;

  /// Array of pore volume multipliers
  array2d< real64 > m_poreVolumeMultiplier;

  /// Array of pore volume multiplier derivatives
  array2d< real64 > m_dPVMult_dPressure;

  /// Pore volume relationship
  ExponentialRelation< real64, ExponentApproximationType::Linear > m_poreVolumeRelation;

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
