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
 * @file EquilibratedChemicalFluid.cpp
 */

#include "EquilibratedChemicalFluid.hpp"
#include "ThermoDatabases/ThermoDatabaseBase.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

static void interpolation( const string & name, const array1d< real64 > & values1, const array1d< real64 > & values2, const real64 & value1, real64 & value2 )
{

  if( value1 < values1[0] || value1 > values1[values1.size()-1] )
    GEOSX_ERROR( "variable: " + name + " is out of bound!" );

  localIndex idx = 0;

  for( localIndex i = 1; i < values1.size(); ++i )
    if( value1 <= values1[i] )
    {
      idx = i;
      break;
    }

  value2 = values2[idx-1] + (values2[idx] - values2[idx-1]) / (values1[idx] - values1[idx-1]) * (value1 - values1[idx-1]);


}



static DatabaseType stringToDatabaseType( string const & databaseType )
{
  if( databaseType == "EQ36" )
  {
    return DatabaseType::EQ36;
  }

  GEOSX_ERROR( "Database type not supported: " << databaseType );

  return DatabaseType::invalidType;
}

static ActivityCoefModel stringToActivityCoefModel( string const & activityCoefModel )
{
  if( activityCoefModel == "BDot" )
  {
    return ActivityCoefModel::BDot;
  }

  GEOSX_ERROR( "Activity coefficient model not supported: " << activityCoefModel );

  return ActivityCoefModel::invalidModel;
}

EquilibratedChemicalFluid::EquilibratedChemicalFluid( std::string const & name, Group * const parent ):
  ReactiveFluidBase( name, parent )
{

  registerWrapper( viewKeyStruct::databaseTypeString, &m_databaseTypeString )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Thermodynamic database" );

  registerWrapper( viewKeyStruct::databaseFileString, &m_databaseFileName )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Thermodynamic database file name" );

  registerWrapper( viewKeyStruct::activityCoefModelString, &m_activityCoefModelString )->
    setInputFlag( InputFlags::REQUIRED )->
    setDescription( "Activity coefficient model" );
}

EquilibratedChemicalFluid::~EquilibratedChemicalFluid() = default;

void EquilibratedChemicalFluid::allocateConstitutiveData( dataRepository::Group * const parent,
                                                          localIndex const numConstitutivePointsPerParentIndex )
{
  ReactiveFluidBase::allocateConstitutiveData( parent, numConstitutivePointsPerParentIndex );

  m_databaseType   = stringToDatabaseType( m_databaseTypeString );
  m_activityCoefModel = stringToActivityCoefModel( m_activityCoefModelString );

  ReadDatabase();

  bool HplusNotFound = 1;
  bool H2OFound = 0;

  localIndex const NBasis = numBasisSpecies();

  m_isHplus.resize( NBasis );

  for( localIndex id = 0; id < NBasis; ++id )
  {

    if( m_basisSpeciesNames[id] == "H+" )
    {
      HplusNotFound = 0;
      m_isHplus[id] = 1;
    }
    else
    {
      m_isHplus[id] = 0;
    }

    if( m_basisSpeciesNames[id] == "H2O" )
      H2OFound = 1;

  }

  GEOSX_ERROR_IF( HplusNotFound, "ReactiveFluidBase: H+ is not specified in basisSpeciesNames" );

  GEOSX_ERROR_IF( H2OFound, "ReactiveFluidBase: H2O cannot be specified in basisSpeciesNames" );

  ResizeFields( parent->size() );

}

void EquilibratedChemicalFluid::PostProcessInput()
{
  ReactiveFluidBase::PostProcessInput();

}

void EquilibratedChemicalFluid::ResizeFields( localIndex size )
{

  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();

  m_dependentConc.resize( size, NDependent + 1 );
  m_dDependentConc_dConc.resize( size, NDependent + 1, NBasis + 1 );

}

void EquilibratedChemicalFluid::ReadDatabase()
{
  if( !m_thermoDatabase )
    m_thermoDatabase = ( ThermoDatabaseBase::CatalogInterface::Factory( m_databaseTypeString, m_databaseFileName, m_basisSpeciesNames ) );

  const array1d< Species > & dependentSpecies = m_thermoDatabase->GetDependentSpecies();

  m_dependentSpeciesNames.clear();
  for( localIndex ic = 0; ic < dependentSpecies.size(); ++ic )
  {
    m_dependentSpeciesNames.emplace_back( dependentSpecies[ic].name );
  }

  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();

  m_stochMatrix.resize( NBasis, NDependent );

  m_stochMatrix.setValues< serialPolicy >( 0 );

  for( localIndex j = 0; j < dependentSpecies.size(); ++j )
  {

    real64 nu2 = dependentSpecies[j].stochs[0];

    bool isSolid = dependentSpecies[j].type == SpeciesType::Solid ? 1 : 0;
    bool isNonBasisGas = (dependentSpecies[j].name != "O2(g)" && dependentSpecies[j].type == SpeciesType::Gas) ? 1 : 0;

    for( localIndex i = 1; i < dependentSpecies[j].speciesIndices.size(); ++i )
    {

      localIndex ic = dependentSpecies[j].speciesIndices[i];

      if( ic >= NBasis )
        continue;

      real64 nu1 = dependentSpecies[j].stochs[i];

      if( isSolid || isNonBasisGas )
        m_stochMatrix[ic][j] = 0.0;
      else
        m_stochMatrix[ic][j] = nu1 / nu2;

    }

  }

}

void EquilibratedChemicalFluid::PointUpdate( real64 const & pressure, real64 const & temperature, arraySlice1d< real64 const > const & concentration, localIndex const k )
{

  Compute( pressure, temperature, concentration, m_dependentConc[k], m_dDependentConc_dConc[k], m_thermoDatabase );

}

void EquilibratedChemicalFluid::ComputeLogActCoef( real64 const & pressure,
                                                   real64 const & temperature,
                                                   real64 const & ionicStrength,
                                                   array1d< real64 > & logActCoef1,
                                                   array1d< real64 > & dLogActCoef1,
                                                   array1d< real64 > & logActCoef2,
                                                   array1d< real64 > & dLogActCoef2 )
{
  GEOSX_UNUSED_VAR( pressure );
  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();

  const array1d< Species > & dependentSpecies = m_thermoDatabase->GetDependentSpecies();
  const array1d< Species > & basisSpecies = m_thermoDatabase->GetBasisSpecies();
  const array1d< localIndex > & basisSpeciesIndices = m_thermoDatabase->GetBasisSpeciesIndices();

  const ActCoefParameters & actCoefParameters = m_thermoDatabase->GetActCoefParameters();

  logActCoef1.resize( NBasis );
  dLogActCoef1.resize( NBasis );

  logActCoef2.resize( NDependent );
  dLogActCoef2.resize( NDependent );

  real64 DHA, DHB, Bdot, DHazero, charge;

  interpolation( "DHA", actCoefParameters.temperatures, actCoefParameters.DHAs, temperature, DHA );

  interpolation( "DHB", actCoefParameters.temperatures, actCoefParameters.DHBs, temperature, DHB );

  static real64 C=-1.0312;
  static real64 F=0.0012806;
  static real64 G=255.9;
  static real64 E=0.4445;
  static real64 H=-0.001606;
  real64 TK = temperature + 273.15;

  if( m_activityCoefModel == ActivityCoefModel::BDot )
  {

    interpolation( "Bdot", actCoefParameters.temperatures, actCoefParameters.BDots, temperature, Bdot );

    for( localIndex i = 0; i < NBasis; ++i )
    {
      localIndex ic = basisSpeciesIndices[i];

      DHazero = basisSpecies[ic].DHazero;
      charge = basisSpecies[ic].charge;

      logActCoef1[i] = Bdot * ionicStrength - DHA * charge * charge * sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength ));

      dLogActCoef1[i] = Bdot - DHA * charge * charge *
                        (0.5 / sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength )) - 0.5 * DHazero * DHB / (1.0 + DHazero * DHB * sqrt( ionicStrength )) /
                         (1.0 + DHazero * DHB * sqrt( ionicStrength )));

    }

    for( localIndex i = 0; i < dependentSpecies.size(); ++i )
    {
      DHazero = dependentSpecies[i].DHazero;
      charge = dependentSpecies[i].charge;

      if( dependentSpecies[i].name == "CO2(aq)" || dependentSpecies[i].name == "H2(aq)" || dependentSpecies[i].name == "O2(aq)" )
      {

        logActCoef2[i] = ((C + F * TK + G / TK) * ionicStrength -(E + H * TK) * (ionicStrength / (ionicStrength+1) )) / 2.303;

        dLogActCoef2[i] = ((C + F * TK + G / TK) - (E + H * TK) * (1.0 / (ionicStrength + 1) - ionicStrength / (ionicStrength + 1) / (ionicStrength + 1))) / 2.303;

      }
      else if( strstr( dependentSpecies[i].name.c_str(), "(aq)" ) || dependentSpecies[i].type == SpeciesType::Gas || dependentSpecies[i].type == SpeciesType::Solid )
      {

        logActCoef2[i] = 0;
        dLogActCoef2[i] = 0;


      }
      else
      {

        logActCoef2[i] = Bdot * ionicStrength - DHA * charge * charge * sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength ));

        dLogActCoef2[i] = Bdot - DHA * charge * charge *
                          (0.5 / sqrt( ionicStrength ) / (1.0 + DHazero * DHB * sqrt( ionicStrength )) - 0.5 * DHazero * DHB / (1.0 + DHazero * DHB * sqrt( ionicStrength )) /
                           (1.0 + DHazero * DHB * sqrt( ionicStrength )));

      }

    }
  }
  else
    GEOSX_ERROR( "wrong activity coef model" );

}

void EquilibratedChemicalFluid::Compute( real64 const & pressure,
                                         real64 const & temperature,
                                         arraySlice1d< real64 const > const & concentration,
                                         arraySlice1d< real64 > const & dependentConc,
                                         arraySlice2d< real64 > const & dDependentConc_dConc,
                                         ThermoDatabase & thermoDatabase )
{

  localIndex const NBasis = numBasisSpecies();
  localIndex const NDependent = numDependentSpecies();

  const array1d< Species > & dependentSpecies = thermoDatabase->GetDependentSpecies();

  const array1d< Species > & basisSpecies = m_thermoDatabase->GetBasisSpecies();
  const array1d< localIndex > & basisSpeciesIndices = m_thermoDatabase->GetBasisSpeciesIndices();

  const ActCoefParameters & actCoefParameters = m_thermoDatabase->GetActCoefParameters();

  static bool firstTime = 1;

  //compute d(logActCoef)/dI

  array1d< real64 > logActCoef1;
  array1d< real64 > dLogActCoef1;

  array1d< real64 > logActCoef2;
  array1d< real64 > dLogActCoef2;


  //get ionic strength

  dependentConc[NDependent] = 0;

  for( localIndex i = 0; i < NBasis; ++i )
  {
    localIndex id = basisSpeciesIndices[i];
    real64 charge = basisSpecies[id].charge;
    real64 conc = pow( 10.0, concentration[i] );

    dependentConc[NDependent] += 0.5 * charge * charge * conc;
  }

  for( localIndex i = 0; i < dependentSpecies.size(); ++i )
  {
    real64 charge = dependentSpecies[i].charge;
    real64 conc = pow( 10.0, dependentConc[i] );

    dependentConc[NDependent] += 0.5 * charge * charge * conc;

  }

  if( firstTime )
  {
    dependentConc[NDependent] = 1.0; // Initial guess, need to fix it later
    firstTime = 0;
  }

  ComputeLogActCoef( pressure, temperature, dependentConc[NDependent], logActCoef1, dLogActCoef1, logActCoef2, dLogActCoef2 );

  real64 logK;

  for( localIndex i = 0; i < NDependent + 1; ++i )
    for( localIndex j = 0; j < NBasis + 1; ++j )
      dDependentConc_dConc[i][j] = 0;

  for( localIndex i = 0; i < dependentSpecies.size(); ++i )
  {

    if( dependentSpecies[i].type == SpeciesType::Solid )
    {
      dependentConc[i] = 0.0;
      continue;
    }


    real64 nu1 = dependentSpecies[i].stochs[0];
    interpolation( "logK", actCoefParameters.temperatures, dependentSpecies[i].logKs, temperature, logK );

    dependentConc[i] = logK / nu1 - logActCoef2[i];

    for( localIndex j = 1; j < dependentSpecies[i].speciesIndices.size(); ++j )
    {

      real64 nu2 = dependentSpecies[i].stochs[j];
      localIndex ic = dependentSpecies[i].speciesIndices[j];

      if( ic == NBasis )
      {
        dependentConc[i] -= nu2 / nu1 * m_logActH2O;
      }
      else if( ic == NBasis + 1 )
        dependentConc[i] -= nu2 / nu1 * m_logFO2g;
      else
      {
        dependentConc[i] -= nu2 / nu1 * (concentration[ic] + logActCoef1[ic]);
        dDependentConc_dConc[i][ic] -= nu2 / nu1;

      }


    }

  }


}


REGISTER_CATALOG_ENTRY( ConstitutiveBase, EquilibratedChemicalFluid, std::string const &, Group * const )

} /* namespace constitutive */

} /* namespace geosx */
