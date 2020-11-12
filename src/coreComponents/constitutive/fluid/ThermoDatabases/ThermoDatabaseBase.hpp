
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
 * @file ThermoDatabaseBase.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_THERMODATABASEBASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_THERMODATABASEBASE_HPP

#include "constitutive/ConstitutiveBase.hpp"
#include "common/DataTypes.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{

namespace constitutive
{

enum class SpeciesType
{
  Aqueous,
  Gas,
  Solid,
  Liquid,
  invalidType
};


struct Species
{
  string name;
  SpeciesType type;
  real64 MW;
  real64 DHazero;
  real64 charge;
  array1d< localIndex > speciesIndices;
  array1d< real64 > stochs;
  array1d< real64 > logKs;

};

struct ActCoefParameters
{
  array1d< real64 > pressures;
  array1d< real64 > temperatures;
  array1d< real64 > DHAs;
  array1d< real64 > DHBs;
  array1d< real64 > BDots;

};

class ThermoDatabaseBase
{
public:

  ThermoDatabaseBase( const Path & fileName ):
    m_fileName( fileName )
  {}

  virtual ~ThermoDatabaseBase(){}

  using CatalogInterface = dataRepository::CatalogInterface< ThermoDatabaseBase, Path const &, string_array const & >;

  static typename CatalogInterface::CatalogType & GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  virtual string GetCatalogName() = 0;

  virtual const array1d< Species > & GetBasisSpecies() const = 0;

  virtual const array1d< Species > & GetDependentSpecies() const = 0;

  virtual const array1d< localIndex > & GetBasisSpeciesIndices() const = 0;

  virtual const ActCoefParameters & GetActCoefParameters() const = 0;

protected:
  Path m_fileName;

};

typedef std::unique_ptr< ThermoDatabaseBase > ThermoDatabase;

}

}

#endif
