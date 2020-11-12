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
 * @file ReactiveFluidBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_REACTIVEFLUIDBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_REACTIVEFLUIDBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "ThermoDatabases/ThermoDatabaseBase.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{

namespace constitutive
{

enum class DatabaseType
{
  EQ36,
  invalidType
};

enum class ActivityCoefModel
{
  DebyeHukel,
  BDot,
  Pitzer,
  Davies,
  invalidModel
};

class ReactiveFluidBase : public ConstitutiveBase
{
public:

  ReactiveFluidBase( std::string const & name, Group * const parent );

  virtual ~ReactiveFluidBase() override;

  // *** ConstitutiveBase interface


  virtual void allocateConstitutiveData( dataRepository::Group * const parent,
                                         localIndex const numConstitutivePointsPerParentIndex ) override;

  static constexpr localIndex MAX_NUM_SPECIES = 32;

  // *** ReactiveFluidBase-specific interface

  virtual void PointUpdate( real64 const & pressure, real64 const & temperature, arraySlice1d< real64 const > const & concentration, localIndex const k ) = 0;


  localIndex numBasisSpecies() const
  {
    return m_basisSpeciesNames.size();
  }

  const string_array & basisiSpeciesNames() const
  {
    return m_basisSpeciesNames;
  }

  const array1d< bool > & IsHplus() const
  {
    return m_isHplus;
  }


  localIndex numDependentSpecies() const
  {
    return m_dependentSpeciesNames.size();
  }

  const string_array & dependentSpeciesNames() const
  {
    return m_dependentSpeciesNames;
  }

  const array2d< real64 > & StochMatrix() const
  {
    return m_stochMatrix;
  }

  // *** Data repository keys

  struct viewKeyStruct
  {

    static constexpr auto basisSpeciesNamesString    = "basisSpeciesNames";
    static constexpr auto logActH2OString    = "logActH2O";
    static constexpr auto logFO2gString    = "logFO2g";

    static constexpr auto dependentConcString      = "dependentConc";
    static constexpr auto dDependentConc_dConcString  = "dDependentConc_dConc";

    using ViewKey = dataRepository::ViewKey;

    ViewKey basisSpeciesNames = { basisSpeciesNamesString };
    ViewKey logActH2O = { logActH2OString };
    ViewKey logFO2g = { logFO2gString };

    ViewKey dependentConc = { dependentConcString };
    ViewKey dDependentConc_dConc = { dDependentConc_dConcString };

  } viewKeysReactiveFluidBase;

protected:

  virtual void PostProcessInput() override;

  string_array m_basisSpeciesNames;
  string_array m_dependentSpeciesNames;
  array1d< bool > m_isHplus;
  array2d< real64 > m_stochMatrix;

  array2d< real64 > m_dependentConc;
  array3d< real64 > m_dDependentConc_dConc;

  real64 m_logFO2g;
  real64 m_logActH2O;

};

} //namespace constitutive

} //namespace geosx

#endif //SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_REACTIVEFLUIDBASE_HPP
