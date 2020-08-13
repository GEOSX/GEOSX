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
 * @file SingleFluidBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDBASE_HPP
#define GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @brief Base class for single-phase fluid model kernel wrappers.
 */
class SingleFluidBaseUpdate
{
public:
  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex
  numElems() const
  {
    return m_density.size( 0 );
  }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex
  numGauss() const
  {
    return m_density.size( 1 );
  };

protected:
  /**
   * @brief Constructor.
   * @param density     fluid density
   * @param dDens_dPres derivative of density w.r.t. pressure
   * @param viscosity   fluid viscosity
   * @param dVisc_dPres derivative of viscosity w.r.t. pressure
   */
  SingleFluidBaseUpdate( arrayView2d< real64 > const & density,
                         arrayView2d< real64 > const & dDens_dPres,
                         arrayView2d< real64 > const & viscosity,
                         arrayView2d< real64 > const & dVisc_dPres ) :
    m_density( density ),
    m_dDens_dPres( dDens_dPres ),
    m_viscosity( viscosity ),
    m_dVisc_dPres( dVisc_dPres )
  {}

  /**
   * @brief Copy constructor.
   */
  SingleFluidBaseUpdate( SingleFluidBaseUpdate const & ) = default;

  /**
   * @brief Move constructor.
   */
  SingleFluidBaseUpdate( SingleFluidBaseUpdate && ) = default;

  /**
   * @brief Deleted copy assignment operator
   * @return reference to this object
   */
  SingleFluidBaseUpdate &
  operator=( SingleFluidBaseUpdate const & ) = delete;

  /**
   * @brief Deleted move assignment operator
   * @return reference to this object
   */
  SingleFluidBaseUpdate &
  operator=( SingleFluidBaseUpdate && ) = delete;

  /// Fluid density
  arrayView2d< real64 > m_density;

  /// Derivative of density w.r.t. pressure
  arrayView2d< real64 > m_dDens_dPres;

  /// Fluid viscosity
  arrayView2d< real64 > m_viscosity;

  /// Derivative of viscosity w.r.t. pressure
  arrayView2d< real64 > m_dVisc_dPres;

private:
  /**
   * @brief Compute fluid properties at a single point.
   * @param[in]  pressure the target pressure value
   * @param[out] density fluid density
   * @param[out] viscosity fluid viscosity
   */
  GEOSX_HOST_DEVICE
  virtual void
  Compute( real64 const pressure,
           real64 & density,
           real64 & viscosity ) const = 0;

  /**
   * @brief Compute fluid properties and derivatives at a single point.
   * @param[in]  pressure the target pressure value
   * @param[out] density fluid density
   * @param[out] dDensity_dPressure fluid density derivative w.r.t. pressure
   * @param[out] viscosity fluid viscosity
   * @param[out] dViscosity_dPressure fluid viscosity derivative w.r.t. pressure
   */
  GEOSX_HOST_DEVICE
  virtual void
  Compute( real64 const pressure,
           real64 & density,
           real64 & dDensity_dPressure,
           real64 & viscosity,
           real64 & dViscosity_dPressure ) const = 0;

  /**
   * @brief Update fluid state at a single point.
   * @param[in] k        element index
   * @param[in] q        gauss point index
   * @param[in] pressure the target pressure value
   */
  GEOSX_HOST_DEVICE
  virtual void
  Update( localIndex const k,
          localIndex const q,
          real64 const pressure ) const = 0;
};

/**
 * @brief Base class for single-phase fluid models.
 */
class SingleFluidBase : public ConstitutiveBase
{
public:
  /**
   * @brief Constructor.
   * @param name name of the group
   * @param parent pointer to parent group
   */
  SingleFluidBase( std::string const & name, Group * const parent );

  /**
   * @brief Destructor.
   */
  virtual ~SingleFluidBase() override;

  // *** ConstitutiveBase interface

  virtual void
  DeliverClone(
    string const & name,
    Group * const parent,
    std::unique_ptr< ConstitutiveBase > & clone ) const override = 0;

  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex ) override;

  // *** SingleFluid-specific interface

  arrayView2d< real64 const > const &
  density() const
  {
    return m_density;
  }
  arrayView2d< real64 > const &
  density()
  {
    return m_density;
  }

  arrayView2d< real64 const > const &
  dDensity_dPressure() const
  {
    return m_dDensity_dPressure;
  }
  arrayView2d< real64 > const &
  dDensity_dPressure()
  {
    return m_dDensity_dPressure;
  }

  arrayView2d< real64 const > const &
  viscosity() const
  {
    return m_viscosity;
  }
  arrayView2d< real64 > const &
  viscosity()
  {
    return m_viscosity;
  }

  arrayView2d< real64 const > const &
  dViscosity_dPressure() const
  {
    return m_dViscosity_dPressure;
  }
  arrayView2d< real64 > const &
  dViscosity_dPressure()
  {
    return m_dViscosity_dPressure;
  }

  real64
  defaultDensity() const
  {
    return m_defaultDensity;
  }
  real64
  defaultViscosity() const
  {
    return m_defaultViscosity;
  }

  // *** Data repository keys

  struct viewKeyStruct
  {
    static constexpr auto defaultDensityString = "defaultDensity";
    static constexpr auto densityString = "density";
    static constexpr auto dDens_dPresString = "dDensity_dPressure";

    static constexpr auto defaultViscosityString = "defaultViscosity";
    static constexpr auto viscosityString = "viscosity";
    static constexpr auto dVisc_dPresString = "dViscosity_dPressure";
  } viewKeysSingleFluidBase;

protected:
  virtual void
  PostProcessInput() override;

  real64 m_defaultDensity;
  real64 m_defaultViscosity;

  array2d< real64 > m_density;
  array2d< real64 > m_dDensity_dPressure;

  array2d< real64 > m_viscosity;
  array2d< real64 > m_dViscosity_dPressure;
};

}  //namespace constitutive

}  //namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_SINGLEFLUIDBASE_HPP
