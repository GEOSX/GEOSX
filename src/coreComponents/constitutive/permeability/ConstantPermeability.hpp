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
 * @file ConstantPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_CONSTANTPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_CONSTANTPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ConstantPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ConstantPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                              arrayView3d< real64 > const & dPerm_dPressure )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure )
  {}

  /// Default copy constructor
  ConstantPermeabilityUpdate( ConstantPermeabilityUpdate const & ) = default;

  /// Default move constructor
  ConstantPermeabilityUpdate( ConstantPermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  ConstantPermeabilityUpdate & operator=( ConstantPermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  ConstantPermeabilityUpdate & operator=( ConstantPermeabilityUpdate && ) = delete;

private:

};


class ConstantPermeability : public PermeabilityBase
{
public:
  ConstantPermeability( string const & name, Group * const parent );

  virtual ~ConstantPermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ConstantPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ConstantPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * permeabilityComponentsString() { return "permeabilityComponents"; }
  } viewKeys;

protected:
  virtual void postProcessInput() override;

private:
  R1Tensor m_permeabilityComponents;

};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
