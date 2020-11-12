
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
 * @file EQ36Database.hpp
 */

#ifndef SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_EQ36DATABASE_HPP
#define SRC_COMPONENTS_CORE_SRC_CONSTITUTIVE_EQ36DATABASE_HPP

#include "ThermoDatabaseBase.hpp"

namespace geosx
{

namespace constitutive
{

class EQ36Database : public ThermoDatabaseBase
{
public:

  EQ36Database( const Path & fileName,
                const string_array & basisSpeciesNames );

  ~EQ36Database() override
  {}

  static constexpr auto m_catalogName = "EQ36";

  static string CatalogName()                    { return m_catalogName; }

  virtual string GetCatalogName() override final { return CatalogName(); }

  virtual const array1d< Species > & GetBasisSpecies() const override { return m_basisSpecies; }

  virtual const array1d< Species > & GetDependentSpecies() const override { return m_dependentSpecies; }

  virtual const array1d< localIndex > & GetBasisSpeciesIndices() const override { return m_basisSpeciesIndices; }

  virtual const ActCoefParameters & GetActCoefParameters() const override { return m_actCoefParameters; }


private:

  /** read database and generate species properties and stoch matrix */

  void CreateChemicalSystem( const string_array & basisSpeciesNames );

  array1d< Species > m_basisSpecies;
  array1d< Species > m_dependentSpecies;

  array1d< localIndex > m_basisSpeciesIndices;

  ActCoefParameters m_actCoefParameters;

};

}

}

#endif
