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

protected:

  PermeabilityBaseUpdate( arrayView3d< real64 > const & permeability )
    : m_permeability( permeability )
  {}

  /// Default copy constructor
  PermeabilityBaseUpdate( PermeabilityBaseUpdate const & ) = default;

  /// Default move constructor
  PermeabilityBaseUpdate( PermeabilityBaseUpdate && ) = default;

  /// Deleted copy assignment operator
  PermeabilityBaseUpdate & operator=( PermeabilityBaseUpdate const & ) = delete;

  /// Deleted move assignment operator
  PermeabilityBaseUpdate & operator=( PermeabilityBaseUpdate && ) = delete;

  arrayView3d< real64 > m_permeability;

  arrayView3d< real64 > m_dPerm_dPressure;
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

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * permeabilityString() { return "permeability"; }
    static constexpr char const * dPerm_dPressureString() {return "dPerm_dPressure";}
  } viewKeys;

protected:
  virtual void postProcessInput() override;

  array3d< real64 > m_permeability;

  array3d< real64 > m_dPerm_dPressure;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYBASE_HPP_