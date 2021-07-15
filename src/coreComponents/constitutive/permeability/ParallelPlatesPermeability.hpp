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
 * @file ParallelPlatesPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_

#include "constitutive/permeability/PermeabilityBase.hpp"


namespace geosx
{
namespace constitutive
{

class ParallelPlatesPermeabilityUpdate : public PermeabilityBaseUpdate
{
public:

  ParallelPlatesPermeabilityUpdate( arrayView2d< real64 > const & permeability,
                                    arrayView2d< real64 > const & dPerm_dPressure,
                                    arrayView2d< real64 > const & dPerm_dAperture )
    : PermeabilityBaseUpdate( permeability, dPerm_dPressure ),
    m_dPerm_dAperture( dPerm_dAperture )
  {}

  GEOSX_HOST_DEVICE
  void compute( real64 const & effectiveAperture,
                arraySlice1d< real64 > const & permeability,
                arraySlice1d< real64 > const & dPerm_dAperture ) const;

  GEOSX_HOST_DEVICE
  void updateFromAperture( localIndex const k,
                           localIndex const q,
                           real64 const & effectiveAperture ) const override
  {
    GEOSX_UNUSED_VAR( q );

    compute( effectiveAperture,
             m_permeability[k],
             m_dPerm_dAperture[k] );
  }
private:
  arrayView2d< real64 > m_dPerm_dAperture;

};


class ParallelPlatesPermeability : public PermeabilityBase
{
public:
  ParallelPlatesPermeability( string const & name, Group * const parent );

  virtual ~ParallelPlatesPermeability() override = default;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "ParallelPlatesPermeability"; }

  virtual string getCatalogName() const override { return catalogName(); }

  /// Type of kernel wrapper for in-kernel update
  using KernelWrapper = ParallelPlatesPermeabilityUpdate;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelWrapper() const
  {
    return KernelWrapper( m_permeability,
                          m_dPerm_dPressure,
                          m_dPerm_dAperture );
  }


  struct viewKeyStruct : public PermeabilityBase::viewKeyStruct
  {} viewKeys;

private:
  array2d< real64 > m_dPerm_dAperture;

};


GEOSX_HOST_DEVICE
void ParallelPlatesPermeabilityUpdate::compute( real64 const & effectiveAperture,
                                                arraySlice1d< real64 > const & permeability,
                                                arraySlice1d< real64 > const & dPerm_dAperture ) const
{
  // Technically this is not a permeability but it's convenient to definite like this.
  real64 const perm  = effectiveAperture*effectiveAperture*effectiveAperture / 12.0;
  real64 const dPerm = 3*effectiveAperture*effectiveAperture / 12.0;

  for( int dim=0; dim < 3; dim++ )
  {
    permeability[dim]     = perm;
    dPerm_dAperture[dim]  = dPerm;
  }
}



}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PARALLELPLATESPERMEABILITY_HPP_
