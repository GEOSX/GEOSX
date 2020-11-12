
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
 * @file EQ36Database.cpp
 */

#include "EQ36Database.hpp"

using namespace std;

namespace geosx
{

using namespace stringutilities;

namespace constitutive
{

EQ36Database::EQ36Database( const Path & fileName,
                            const string_array & basisSpeciesNames ):
  ThermoDatabaseBase( fileName )
{

  CreateChemicalSystem( basisSpeciesNames );

}

void EQ36Database::CreateChemicalSystem( const string_array & basisSpeciesNames )
{

  std::ifstream is( m_fileName );

  constexpr std::streamsize buf_size = 256;
  char buf[buf_size];

  unordered_map< string, int > basisSpeciesMap;
  unordered_map< string, double > DAzeroSpeciesMap;

  for( localIndex ic = 0; ic < basisSpeciesNames.size(); ++ic )
  {
    basisSpeciesMap[basisSpeciesNames[ic] ] = -1;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Temperature grid" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Pressure grid" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      m_actCoefParameters.temperatures.emplace_back( std::stod( strs[i] ));
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Pressure envelope" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      m_actCoefParameters.pressures.emplace_back( std::stod( strs[i] ));
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel A_gamma" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel A_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      m_actCoefParameters.DHAs.emplace_back( std::stod( strs[i] ));
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel B_gamma" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "Debye-Huckel B_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      m_actCoefParameters.DHBs.emplace_back( std::stod( strs[i] )*1e8 );
    }
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "B-dot" );
    if( found!=std::string::npos )
      break;
  }

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "B-dot_H" );
    if( found!=std::string::npos )
      break;

    string_array strs = Tokenize( str, " " );
    for( localIndex i = 0; i < strs.size(); ++i )
    {
      m_actCoefParameters.BDots.emplace_back( std::stod( strs[i] ));
    }
  }

  GEOSX_ERROR_IF(
    m_actCoefParameters.temperatures.size() != m_actCoefParameters.pressures.size() || m_actCoefParameters.temperatures.size() != m_actCoefParameters.DHAs.size() || m_actCoefParameters.temperatures.size() != m_actCoefParameters.DHBs.size() || m_actCoefParameters.temperatures.size() != m_actCoefParameters.BDots.size(),
    "Internal error when reading database" );

  /* read basis species */

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    auto found = str.find( "basis species" );
    if( found!=std::string::npos )
      break;
  }

  string speciesName;
  real64 MW = 0;
  real64 charge = 0;
  real64 DHazero = 0;
  int H2OIndex = -1, O2gIndex = -1;

  int count = 0;

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        auto it = basisSpeciesMap.find( speciesName );
        if( it != basisSpeciesMap.end() || speciesName == "H2O" ||speciesName == "O2(g)" )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Aqueous;
          entry.MW = MW;
          entry.DHazero = DHazero;
          entry.charge = charge;
          m_basisSpecies.emplace_back( entry );

          if( speciesName == "H2O" )
            H2OIndex = count;
          else if( speciesName == "O2(g)" )
            O2gIndex = count;
          else
            basisSpeciesMap[speciesName] = count;

          count++;

        }

        is.getline( buf, buf_size );
        speciesName = buf;
        auto found2 = speciesName.find( "auxiliary basis species" );

        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;
      }
    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {
        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );
      }
    }
  }


  localIndex idx;

  localIndex numBasisSpecies = basisSpeciesNames.size();
  m_basisSpeciesIndices.resize( numBasisSpecies + 2 );

  for( localIndex ic = 0; ic < numBasisSpecies; ++ic )
  {
    idx = basisSpeciesMap[basisSpeciesNames[ic] ];
    m_basisSpeciesIndices[ic] = idx;
    basisSpeciesMap[basisSpeciesNames[ic] ] = int(ic);
  }

  m_basisSpeciesIndices[numBasisSpecies] = H2OIndex;
  basisSpeciesMap[m_basisSpecies[H2OIndex].name] = int(numBasisSpecies);

  m_basisSpeciesIndices[numBasisSpecies + 1] = O2gIndex;
  basisSpeciesMap[m_basisSpecies[O2gIndex].name ] = int(numBasisSpecies + 1);


  /* read aux basis species */

  string_array speciesNames;
  array1d< localIndex > speciesIndices;
  array1d< real64 > stochs;
  array1d< real64 > logKs;

  //bool isSpecies;

  count = 0;

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        if( speciesIndices.size() > 1 )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Aqueous;
          entry.MW = MW;
          entry.DHazero = DHazero;
          entry.charge = charge;
          entry.speciesIndices = speciesIndices;
          entry.stochs = stochs;
          entry.logKs = logKs;

          m_dependentSpecies.emplace_back( entry );

          speciesIndices.clear();
          stochs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "aqueous species" );

        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;

      }

    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );

      }
    }

    {

      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        speciesNames.clear();
        stochs.clear();

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;



          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = basisSpeciesMap.find( speciesNames[i] );
          if( it != basisSpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( localIndex i = 0; i < strs2.size(); ++i )
              logKs.emplace_back( std::stod( strs2[i] ));

          }

        }

      }
    }

  }



  /* read aux aqueous species */

  speciesIndices.clear();
  stochs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        if( speciesIndices.size() > 1 )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Aqueous;
          entry.MW = MW;
          entry.DHazero = DHazero;
          entry.charge = charge;
          entry.speciesIndices = speciesIndices;
          entry.stochs = stochs;
          entry.logKs = logKs;

          m_dependentSpecies.emplace_back( entry );

          speciesIndices.clear();
          stochs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "solids" );
        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;

      }

    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );

      }
    }

    {

      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        speciesNames.clear();
        stochs.clear();

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;


          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = basisSpeciesMap.find( speciesNames[i] );
          if( it != basisSpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( localIndex i = 0; i < strs2.size(); ++i )
              logKs.emplace_back( std::stod( strs2[i] ));

          }

        }

      }
    }

  }

  /* read solid species */

  speciesIndices.clear();
  stochs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        if( speciesIndices.size() > 1 )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Solid;
          entry.MW = MW;
          entry.DHazero = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stochs = stochs;
          entry.logKs = logKs;

          m_dependentSpecies.emplace_back( entry );

          speciesIndices.clear();
          stochs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "liquids" );
        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;

      }

    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );

      }
    }

    {

      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        speciesNames.clear();
        stochs.clear();

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;


          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = basisSpeciesMap.find( speciesNames[i] );
          if( it != basisSpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( localIndex i = 0; i < strs2.size(); ++i )
              logKs.emplace_back( std::stod( strs2[i] ));

          }

        }

      }
    }

  }


  /* read liquid species */

  speciesIndices.clear();
  stochs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        if( speciesIndices.size() > 1 )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Liquid;
          entry.MW = MW;
          entry.DHazero = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stochs = stochs;
          entry.logKs = logKs;

          m_dependentSpecies.emplace_back( entry );

          speciesIndices.clear();
          stochs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "gases" );

        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;

      }

    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );

      }
    }

    {

      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        speciesNames.clear();
        stochs.clear();

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;


          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = basisSpeciesMap.find( speciesNames[i] );
          if( it != basisSpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( localIndex i = 0; i < strs2.size(); ++i )
              logKs.emplace_back( std::stod( strs2[i] ));

          }

        }

      }
    }

  }


  /* read gases species */

  speciesIndices.clear();
  stochs.clear();
  logKs.clear();
  speciesNames.clear();

  while( is.getline( buf, buf_size ))
  {
    std::string str( buf );
    {

      auto found = str.find( "+-------" );
      if( found!=std::string::npos )
      {

        if( speciesIndices.size() > 1 )
        {

          Species entry;
          entry.name = speciesName;
          entry.type = SpeciesType::Gas;
          entry.MW = MW;
          entry.DHazero = 0;
          entry.charge = 0;
          entry.speciesIndices = speciesIndices;
          entry.stochs = stochs;
          entry.logKs = logKs;

          m_dependentSpecies.emplace_back( entry );

          speciesIndices.clear();
          stochs.clear();
          logKs.clear();
          speciesNames.clear();

          count++;

        }

        is.getline( buf, buf_size );
        std::string str2( buf );
        string_array strs2 = Tokenize( str2, " " );
        speciesName = strs2[0];

        auto found2 = str2.find( "solid solutions" );
        if( found2 != std::string::npos )
        {
          break;
        }

      }

    }

    {
      auto found = str.find( "mol.wt." );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        MW = std::stod( strs[3] ) * 0.001;

      }

    }

    {
      auto found = str.find( "DHazero" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        DHazero = std::stod( strs[3] ) * 1e-8;

      }

    }

    {

      auto found = str.find( "charge" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        charge = std::stod( strs[2] );

      }
    }

    {

      auto found = str.find( "aqueous dissociation reaction" );
      if( found != std::string::npos )
      {

        string_array strs = Tokenize( str, " " );
        localIndex num1 = localIndex( std::stoi( strs[0] ));

        speciesNames.clear();
        stochs.clear();

        while( is.getline( buf, buf_size ))
        {

          std::string str2( buf );
          auto found2 = str2.find( "* Log K" );
          if( found2 != std::string::npos )
            break;


          string_array strs2 = Tokenize( str2, " " );
          //localIndex num2 = strs2.size();

          for( localIndex i = 0; i < strs2.size(); ++i )
          {
            if( i % 2 == 0 )
              stochs.emplace_back( std::stod( strs2[i] ));
            else
              speciesNames.emplace_back( strs2[i] );
          }

        }

        GEOSX_ERROR_IF( num1 != speciesNames.size() || num1 != stochs.size() || speciesName != speciesNames[0], "Internal error when reading database" );

        bool notFound = 0;

        speciesIndices.resize( num1 );
        speciesIndices[0] = count;

        for( localIndex i = 1; i < num1; ++i )
        {
          auto it = basisSpeciesMap.find( speciesNames[i] );
          if( it != basisSpeciesMap.end())
          {
            speciesIndices[i] = it->second;
          }
          else
          {
            notFound = 1;
            break;
          }
        }

        if( notFound )
        {
          speciesIndices.clear();
        }
        else
        {
          logKs.clear();
          while( is.getline( buf, buf_size ))
          {

            std::string str2( buf );
            auto found2 = str2.find( "*" );
            if( found2 != std::string::npos )
              break;

            string_array strs2 = Tokenize( str2, " " );

            for( localIndex i = 0; i < strs2.size(); ++i )
              logKs.emplace_back( std::stod( strs2[i] ));

          }

        }

      }
    }

  }

  is.close();

}

REGISTER_CATALOG_ENTRY( ThermoDatabaseBase,
                        EQ36Database,
                        const Path &, const string_array & )

}

}
