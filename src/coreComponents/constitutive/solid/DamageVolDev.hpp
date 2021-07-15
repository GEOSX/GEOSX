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
 * @file Damage.hpp
 * @brief Overrides the SSLE constitutive updates to account for a damage varible and a volumtric-deviatoric split.
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_DAMAGEVOLDEV_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_DAMAGEVOLDEV_HPP_
#include "Damage.hpp"
#include "InvariantDecompositions.hpp"
#include "SolidBase.hpp"

namespace geosx
{
namespace constitutive
{

template< typename UPDATE_BASE >
class DamageVolDevUpdates : public DamageUpdates< UPDATE_BASE >
{
public:
  template< typename ... PARAMS >
  DamageVolDevUpdates( arrayView2d< real64 > const & inputDamage,
                       arrayView2d< real64 > const & inputStrainEnergyDensity,
                       real64 const & inputLengthScale,
                       real64 const & inputCriticalFractureEnergy,
                       real64 const & inputcriticalStrainEnergy,
                       PARAMS && ... baseParams ):
    DamageUpdates< UPDATE_BASE >( inputDamage, inputStrainEnergyDensity, inputLengthScale,
                                  inputCriticalFractureEnergy, inputcriticalStrainEnergy,
                                  std::forward< PARAMS >( baseParams )... )
  {}

  using DiscretizationOps = typename DamageUpdates< UPDATE_BASE >::DiscretizationOps;

  using DamageUpdates< UPDATE_BASE >::smallStrainUpdate;
  using DamageUpdates< UPDATE_BASE >::saveConvergedState;

  using DamageUpdates< UPDATE_BASE >::getDegradationValue;
  using DamageUpdates< UPDATE_BASE >::getDegradationDerivative;
  using DamageUpdates< UPDATE_BASE >::getDegradationSecondDerivative;
  using DamageUpdates< UPDATE_BASE >::getEnergyThreshold;

  using DamageUpdates< UPDATE_BASE >::m_strainEnergyDensity;
  using DamageUpdates< UPDATE_BASE >::m_criticalStrainEnergy;
  using DamageUpdates< UPDATE_BASE >::m_criticalFractureEnergy;
  using DamageUpdates< UPDATE_BASE >::m_lengthScale;

  GEOSX_HOST_DEVICE
  virtual void smallStrainUpdate( localIndex const k,
                                  localIndex const q,
                                  real64 const ( &strainIncrement )[6],
                                  real64 ( & stress )[6],
                                  DiscretizationOps & stiffness ) const final
  {
    // perform elastic update for "undamaged" stress

    UPDATE_BASE::smallStrainUpdate( k, q, strainIncrement, stress, stiffness );  // elastic trial update

    // compute volumetric and deviatoric strain invariants

    real64 strain[6];
    UPDATE_BASE::getElasticStrain( k, q, strain );

    real64 volStrain;
    real64 devStrain;
    real64 deviator[6];

    twoInvariant::strainDecomposition( strain,
                                       volStrain,
                                       devStrain,
                                       deviator );

    // degrade shear stiffness always
    // degrade volumetric stiffness when in tension

    real64 factor = getDegradationValue( k, q );

    stiffness.m_shearModulus *= factor;

    if( volStrain > 0 )
    {
      stiffness.m_bulkModulus *= factor;
    }

    // compute stress invariants and recompose full stress tensor

    real64 stressP = stiffness.m_bulkModulus * volStrain;
    real64 stressQ = 3 * stiffness.m_shearModulus * devStrain;

    twoInvariant::stressRecomposition( stressP,
                                       stressQ,
                                       deviator,
                                       stress );

    // update strain energy density
    // TODO: refactor as a proper history variable update.  the code below doesn't allow for rewinds.

    real64 sed = 0.5 * (stressQ * devStrain) / factor;

    if( volStrain > 0 )
    {
      sed += 0.5 * (stressP * volStrain) / factor;
    }

    if( sed > m_strainEnergyDensity( k, q ) )
    {
      m_strainEnergyDensity( k, q ) = sed;
    }
  }


  GEOSX_HOST_DEVICE
  virtual real64 getStrainEnergyDensity( localIndex const k,
                                         localIndex const q ) const override final
  {
    return m_strainEnergyDensity( k, q );
  }

};


template< typename BASE >
class DamageVolDev : public Damage< BASE >
{
public:

  using KernelWrapper = DamageVolDevUpdates< typename BASE::KernelWrapper >;

  using Damage< BASE >::m_damage;
  using Damage< BASE >::m_strainEnergyDensity;
  using Damage< BASE >::m_criticalFractureEnergy;
  using Damage< BASE >::m_lengthScale;
  using Damage< BASE >::m_criticalStrainEnergy;

  DamageVolDev( string const & name, dataRepository::Group * const parent );
  virtual ~DamageVolDev() override;


  static string catalogName() { return string( "DamageVolDev" ) + BASE::m_catalogNameString; }
  virtual string getCatalogName() const override { return catalogName(); }


  KernelWrapper createKernelUpdates() const
  {
    return BASE::template createDerivedKernelUpdates< KernelWrapper >( m_damage.toView(),
                                                                       m_strainEnergyDensity.toView(),
                                                                       m_lengthScale,
                                                                       m_criticalFractureEnergy,
                                                                       m_criticalStrainEnergy );
  }

};


}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_DAMAGEVOLDEV_HPP_ */
