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

  CarmanKozenyPermeabilityUpdate( arrayView2d< real64 > const & permeability,
                                  arrayView2d< real64 > const & dPerm_dPressure,
                                  arrayView2d< real64 > const & dPerm_dPorosity,
                                  real64 const particleDiameter,
                                  real64 const sphericity )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dPorosity( dPerm_dPorosity ),
    m_particleDiameter( particleDiameter ),
    m_sphericity( sphericity )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const & porosity,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dPorosity ) const;

  GEOSX_HOST_DEVICE
  virtual void updateFromPorosity( localIndex const k,
                                   localIndex const q,
                                   real64 const & porosity ) const override
  {
    GEOSX_UNUSED_VAR( q ); // to do: perm will have to be define per quadrature point.

    compute( porosity,
             m_permeability[k],
             m_dPerm_dPorosity[k] );
  }

private:

  /// dPermeability_dPorosity
  arrayView2d< real64 > m_dPerm_dPorosity;

  /// Particle diameter
  real64 m_particleDiameter;

  /// Sphericity of the particles
  real64 m_sphericity;

};


class CarmanKozenyPermeability : public PermeabilityBase
{
public:
  CarmanKozenyPermeability( string const & name, Group * const parent );

  virtual ~CarmanKozenyPermeability() override = default;

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
  KernelWrapper createKernelWrapper() const
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
  array2d< real64 > m_dPerm_dPorosity;

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
  real64 const constant = pow( m_sphericity*m_particleDiameter, 2 ) / 150;

  real64 const permValue = constant * pow( porosity, 3 )/ pow( (1 - porosity), 2 );

  real64 const dPerm_dPorValue = -constant * ( (porosity - 3) *  pow( porosity, 2 ) / pow( (1-porosity), 3 )  );

  for( localIndex i=0; i < permeability.size(); i++ )
  {
    permeability[i] = permValue;
    dPerm_dPorosity[i] = dPerm_dPorValue;
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_FRACTUREPERMEABILITY_HPP_
