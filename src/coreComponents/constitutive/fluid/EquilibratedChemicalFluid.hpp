/*
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 * Copyright (c) 2019, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory
 *
 * LLNL-CODE-746361
 *
 * All rights reserved. See COPYRIGHT for details.
 *
 * This file is part of the GEOSX Simulation Framework.
 *
 * GEOSX is a free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License (as published by the
 * Free Software Foundation) version 2.1 dated February 1999.
 *~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

/**
 * @file EquilibratedChemicalFluid.hpp
 */
#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_EQUILIBRATEDCHEMICALFLUID_HPP_
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_EQUILIBRATEDCHEMICALFLUID_HPP_

#include "constitutive/fluid/ReactiveFluidBase.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
string const equilibratedChemicalFluid = "EquilibratedChemicalFluid";
}
}

namespace constitutive
{

class EquilibratedChemicalFluid : public ReactiveFluidBase
{
public:

  EquilibratedChemicalFluid( std::string const & name, Group * const parent );

  virtual ~EquilibratedChemicalFluid() override;

  // *** ConstitutiveBase interface

  static std::string CatalogName() { return dataRepository::keys::equilibratedChemicalFluid; }

  virtual string getCatalogName() const override { return CatalogName(); }

  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  virtual void PointUpdate( real64 const & pressure, real64 const & temperature, arraySlice1d< real64 const > const & concentration, localIndex const k ) override;

  // *** Data repository keys

  struct viewKeyStruct : public ReactiveFluidBase::viewKeyStruct
  {
    static constexpr auto databaseTypeString    = "databaseType";
    static constexpr auto databaseFileString      = "databaseFile";
    static constexpr auto activityCoefModelString  = "activityCoefModel";

    dataRepository::ViewKey databaseType    = { databaseTypeString    };
    dataRepository::ViewKey databaseFile      = { databaseFileString      };
    dataRepository::ViewKey activityCoefModel  = { activityCoefModelString  };

  } viewKeysEquilibratedChemicalFluid;

protected:
  virtual void PostProcessInput() override;

private:

  void Compute( real64 const & pressure,
                real64 const & temperature,
                arraySlice1d< const real64 > const & concentration,
                arraySlice1d< real64 > const & dependentConc,
                arraySlice2d< real64 > const & dDependentConc_dConc,
                ThermoDatabase & thermoDatabase );

  void ReadDatabase();

  void ComputeLogActCoef( real64 const & pressure,
                          real64 const & temperature,
                          real64 const & ionicStrength,
                          array1d< real64 > & logActCoef1,
                          array1d< real64 > & dLogActCoef1,
                          array1d< real64 > & logActCoef2,
                          array1d< real64 > & dLogActCoef2 );

  void ResizeFields( localIndex size );


  string m_databaseTypeString;

  Path m_databaseFileName;

  string m_activityCoefModelString;

  DatabaseType m_databaseType;

  ActivityCoefModel m_activityCoefModel;

  ThermoDatabase m_thermoDatabase;

};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_EQUILIBRATEDCHEMICALFLUID_HPP_ */
