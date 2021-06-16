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
 * @file RockBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{

class RockBaseUpdates
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_newPorosity.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_newPorosity.size( 1 ); }

  RockBaseUpdates( arrayView2d< real64 > const & newPorosity,
                   arrayView2d< real64 > const & oldPorosity,
                   arrayView2d< real64 > const & dPorosity_dPressure,
                   arrayView1d< real64 > const & referencePorosity,
                   real64 const & grainBulkModulus,
                   arrayView2d< real64 > const & grainDensity ):
    m_newPorosity( newPorosity ),
    m_oldPorosity( oldPorosity ),
    m_dPorosity_dPressure( dPorosity_dPressure ),
    m_referencePorosity ( referencePorosity ),
    m_grainBulkModulus( grainBulkModulus ),
    m_grainDensity( grainDensity )
  {}

  /// Default copy constructor
  RockBaseUpdates( RockBaseUpdates const & ) = default;

  /// Default move constructor
  RockBaseUpdates( RockBaseUpdates && ) = default;

  /// Deleted copy assignment operator
  RockBaseUpdates & operator=( RockBaseUpdates const & ) = delete;

  /// Deleted move assignment operator
  RockBaseUpdates & operator=( RockBaseUpdates && ) = delete;


  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void updatePorosity( localIndex const GEOSX_UNUSED_PARAM( k ),
                               localIndex const GEOSX_UNUSED_PARAM( q ),
                               real64 const & GEOSX_UNUSED_PARAM( pressure )  ) const
  {}

protected:
  arrayView2d< real64 > m_newPorosity;

  arrayView2d< real64 > m_oldPorosity;

  arrayView2d< real64 > m_dPorosity_dPressure;

  arrayView1d< real64 > m_referencePorosity;

  real64 m_grainBulkModulus;

  arrayView2d< real64 > m_grainDensity;
};


class RockBase : public ConstitutiveBase
{
public:
  RockBase( string const & name, Group * const parent );

  virtual ~RockBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "RockBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * newPorosityString() { return "porosity"; }
    static constexpr char const * oldPorosityString() { return "oldPorosity"; }
    static constexpr char const * dPorosity_dPressureString() { return "dPorosity_dPressure"; }
    static constexpr char const * referencePorosityString() { return "referencePorosity"; }
    static constexpr char const * defaultRefererencePorosityString() { return "defaultReferencePorosity"; }
    static constexpr char const * grainBulkModulusString() { return "grainBulkModulus"; }
    static constexpr char const * grainDensityString() { return "grainDensity"; }
    static constexpr char const * defaultGrainDensityString() { return "defaultGrainDensity"; }    ///< Density key
  } viewKeys;

  /**
   * @brief Const accessor for newPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getPorosity() const { return m_newPorosity; }

  /**
   * @brief Const/non-mutable accessor for oldPorosity.
   * @return Accessor
   */
  arrayView2d< real64 const > const  getOldPorosity() const { return m_oldPorosity; }


  /**
   * @brief Non-Const/mutable accessor for oldPorosity
   * @return Accessor
   */
  arrayView2d< real64 > const getOldPorosity() { return m_oldPorosity; }


  /**
   * @brief Const/non-mutable accessor for dPorosity_dPressure
   * @return Accessor
   */
  arrayView2d< real64 const > const  dPorosity_dPressure() const { return m_dPorosity_dPressure; }

  /**
   * @brief Non-const/Mutable accessor for density.
   * @return Accessor
   */
  arrayView2d< real64 > const getDensity()
  {
    return m_grainDensity;
  }

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return m_grainDensity;
  }

protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_newPorosity;

  array2d< real64 > m_oldPorosity;

  array2d< real64 > m_dPorosity_dPressure;

  array1d< real64 > m_referencePorosity;

  real64 m_defaultReferencePorosity;

  real64 m_grainBulkModulus;

  /// The material density at a quadrature point.
  array2d< real64 > m_grainDensity;

  /// defaultGrainDensity
  real64 m_defaultGrainDensity = 0;
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_ROCKBASE_HPP_
