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
 * @file CompressibleRock.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_
#define GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_

#include "constitutive/ExponentialRelation.hpp"
#include "RockBase.hpp"

namespace geosx
{
namespace constitutive
{

class CompressibleRockUpdates : public RockBaseUpdates
{
public:

  CompressibleRockUpdates( arrayView2d< real64 > const & newPorosity,
                           arrayView2d< real64 > const & oldPorosity,
                           arrayView2d< real64 > const & dPorosity_dPressure,
                           arrayView1d< real64 > const & referencePorosity,
                           real64 const & grainBulkModulus,
                           arrayView2d< real64 > const & grainDensity,
                           real64 const & referencePressure,
                           real64 const & compressibility ):
    RockBaseUpdates( newPorosity,
                     oldPorosity,
                     dPorosity_dPressure,
                     referencePorosity,
                     grainBulkModulus,
                     grainDensity ),
    m_referencePressure( referencePressure ),
    m_compressibility( compressibility )
  {}

  /// Default copy constructor
  CompressibleRockUpdates( CompressibleRockUpdates const & ) = default;

  /// Default move constructor
  CompressibleRockUpdates( CompressibleRockUpdates && ) = default;

  /// Deleted copy assignment operator
  CompressibleRockUpdates & operator=( CompressibleRockUpdates const & ) = delete;

  /// Deleted move assignment operator
  CompressibleRockUpdates & operator=( CompressibleRockUpdates && ) = delete;

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computePorosity( real64 const & pressure,
                        real64 & porosity,
                        real64 & dPorosity_dPressure,
                        real64 const & referencePorosity ) const
  {

    porosity            =  referencePorosity * exp( m_compressibility * (pressure - m_referencePressure) );
    dPorosity_dPressure =  m_compressibility * porosity;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  virtual void updatePorosity( localIndex const k,
                               localIndex const q,
                               real64 const & pressure ) const override
  {
    computePorosity( pressure,
                     m_newPorosity[k][q],
                     m_dPorosity_dPressure[k][q],
                     m_referencePorosity[k] );
  }

private:

  real64 m_referencePressure;

  real64 m_compressibility;

};


class CompressibleRock : public RockBase
{
public:
  CompressibleRock( string const & name, Group * const parent );

  virtual ~CompressibleRock() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static string catalogName() { return "CompressibleRock"; }

  virtual string getCatalogName() const override { return catalogName(); }

  struct viewKeyStruct : public RockBase::viewKeyStruct
  {
    static constexpr char const * referencePressureString() { return "referencePressure"; }
    static constexpr char const * compressibilityString() { return "compressibility"; }
  } viewKeys;

  using KernelWrapper = CompressibleRockUpdates;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  KernelWrapper createKernelUpdates()
  {
    return KernelWrapper( m_newPorosity,
                          m_oldPorosity,
                          m_dPorosity_dPressure,
                          m_referencePorosity,
                          m_grainBulkModulus,
                          m_grainDensity,
                          m_referencePressure,
                          m_compressibility );
  }


protected:
  virtual void postProcessInput() override;

  real64 m_referencePressure;

  real64 m_compressibility;

};


}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_POROSITY_COMPRESSIBLEROCK_HPP_
