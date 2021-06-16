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
 * @file CarmanKozenyPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_CARMANKOZENYPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class CarmanKozenyPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  CarmanKozenyPermeabilityUpdate( arrayView3d< real64 > const & permeability,
                                  arrayView3d< real64 > const & dPerm_dPressure,
                                  arrayView3d< real64 > const & dPerm_dPorosity,
                                  real64 const particleDiameter,
                                  real64 const sphericity )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_particleDiameter( particleDiameter ),
    m_sphericity( sphericity )
  {}

  /// Default copy constructor
  CarmanKozenyPermeabilityUpdate( CarmanKozenyPermeabilityUpdate const & ) = default;

  /// Default move constructor
  CarmanKozenyPermeabilityUpdate( CarmanKozenyPermeabilityUpdate && ) = default;

  /// Deleted copy assignment operator
  CarmanKozenyPermeabilityUpdate & operator=( CarmanKozenyPermeabilityUpdate const & ) = delete;

  /// Deleted move assignment operator
  CarmanKozenyPermeabilityUpdate & operator=( CarmanKozenyPermeabilityUpdate && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void compute( real64 const & porosity,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void updatePorosity( localIndex const k,
                               localIndex const q,
                               real64 const & porosity ) const override
  {
    compute( porosity,
             m_permeability[k][q],
             m_dPerm_dPorosity[k][q] );
  }

private:

  /// dPermeability_dPorosity
  arrayView3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;

};


class CarmanKozenyPermeability : public PermeabilityBase
{
public:
  CarmanKozenyPermeability( string const & name, Group * const parent );

  virtual ~CarmanKozenyPermeability() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "CarmanKozenyPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = CarmanKozenyPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper()
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dPorosity,
                          m_particleDiameter,
                          m_sphericity );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {
    static constexpr char const * dPerm_dPorosityString() { return "dPerm_dPorosity"; }
    static constexpr char const * particleDiameterString() { return "particleDiameter"; }
    static constexpr char const * sphericityString() { return "sphericity"; }
  } viewKeys;

protected:
  virtual void postProcessInput() override;

private:
  /// dPermeability_dPorosity
  array3d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;
};


GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CarmanKozenyPermeabilityUpdate::compute( real64 const & porosity,
                                              arraySlice1d< real64 > const & permeability,
                                              arraySlice1d< real64 > const & dPerm_dPorosity ) const
{
  real64 const permValue = pow( m_sphericity*m_particleDiameter, 2 ) * pow( porosity, 3 )
                           / (150 * pow( (1 - porosity), 2 ) );
  real64 const dPerm_dPorValue = pow( m_sphericity*m_particleDiameter, 2 ) * pow( porosity, 3 );

  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permValue;
    dPerm_dPorosity[i] = dPerm_dPorValue;
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
