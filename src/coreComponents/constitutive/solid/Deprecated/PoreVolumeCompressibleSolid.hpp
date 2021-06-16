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
 * @file PoreVolumeCompressibleSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_

#include "constitutive/ConstitutiveBase.hpp"

#include "constitutive/ExponentialRelation.hpp"

namespace geosx
{
namespace constitutive
{


class PoreVolumeCompressibleSolid : public ConstitutiveBase
{
public:
  PoreVolumeCompressibleSolid( string const & name, Group * const parent );

  virtual ~PoreVolumeCompressibleSolid() override;

  std::unique_ptr< ConstitutiveBase > deliverClone( string const & name,
                                                    Group * const parent ) const override;

  virtual void allocateConstitutiveData( dataRepository::Group & parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;


  static string catalogName() { return "PoreVolumeCompressibleSolid"; }

  virtual string getCatalogName() const override { return catalogName(); }

  virtual void stateUpdatePointPressure( real64 const & pres,
                                         localIndex const k,
                                         localIndex const q ) override final;

  virtual void stateUpdateBatchPressure( arrayView1d< real64 const > const & pres,
                                         arrayView1d< real64 const > const & dPres ) override final;

  struct viewKeyStruct : public ConstitutiveBase::viewKeyStruct
  {
    static constexpr char const * compressibilityString() { return "compressibility"; }
    static constexpr char const * referencePressureString() { return "referencePressure"; }
  };

protected:
  virtual void postProcessInput() override;

private:

  /// scalar compressibility parameter
  real64 m_compressibility;

  /// reference pressure parameter
  real64 m_referencePressure;

  array2d< real64 > m_poreVolumeMultiplier;
  array2d< real64 > m_dPVMult_dPressure;

  ExponentialRelation< real64, ExponentApproximationType::Linear > m_poreVolumeRelation;
};

}/* namespace constitutive */

} /* namespace geosx */


#endif //GEOSX_CONSTITUTIVE_SOLID_POREVOLUMECOMPRESSIBLESOLID_HPP_
