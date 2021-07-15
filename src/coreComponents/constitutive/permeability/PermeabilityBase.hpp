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
 * @file PermeabilityBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{

class PermeabilityBaseUpdate
{
public:

  /**
   * @brief Get number of elements in this wrapper.
   * @return number of elements
   */
  GEOSX_HOST_DEVICE
  localIndex numElems() const { return m_permeability.size( 0 ); }

  /**
   * @brief Get number of gauss points per element.
   * @return number of gauss points per element
   */
  GEOSX_HOST_DEVICE
  localIndex numGauss() const { return m_permeability.size( 1 ); }


  GEOSX_HOST_DEVICE
  virtual void updateFromPorosity( localIndex const k,
                                   localIndex const q,
                                   real64 const & porosity ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( porosity );
  }

  GEOSX_HOST_DEVICE
  virtual void updateFromPressureStrain( localIndex const k,
                                         localIndex const q,
                                         real64 const & pressure,
                                         real64 const & volStrain ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( pressure );
    GEOSX_UNUSED_VAR( volStrain );
  }

  GEOSX_HOST_DEVICE
  virtual void updateFromAperture( localIndex const k,
                                   localIndex const q,
                                   real64 const & aperture ) const
  {
    GEOSX_UNUSED_VAR( k );
    GEOSX_UNUSED_VAR( q );
    GEOSX_UNUSED_VAR( aperture );
  }

protected:

  PermeabilityBaseUpdate( arrayView2d< real64 > const & permeability,
                          arrayView2d< real64 > const & dPerm_dPressure )
    : m_permeability( permeability ),
    m_dPerm_dPressure( dPerm_dPressure )
  {}

  arrayView2d< real64 > m_permeability;

  arrayView2d< real64 > m_dPerm_dPressure;
};


class PermeabilityBase : public ConstitutiveBase
{
public:
  PermeabilityBase( string const & name, Group * const parent );

  virtual ~PermeabilityBase() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "PermeabilityBase"; }

  virtual string getCatalogName() const override { return catalogName(); }

  arrayView2d< real64 const > const permeability() const { return m_permeability; }

  arrayView2d< real64 const > const dPerm_dPressure() const { return m_dPerm_dPressure; }

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * permeabilityString() { return "permeability"; }
    static constexpr char const * dPerm_dPressureString() { return "dPerm_dPressure"; }
    static constexpr char const * dPerm_dApertureString() { return "dPerm_dAperture"; }
  } viewKeys;

protected:
  virtual void postProcessInput() override;

  array2d< real64 > m_permeability;

  array2d< real64 > m_dPerm_dPressure;
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_
