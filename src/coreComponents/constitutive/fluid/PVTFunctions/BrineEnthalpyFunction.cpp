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
 * @file BrineEnthalpyFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineEnthalpyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/CO2EnthalpyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

BrineEnthalpyFunction::BrineEnthalpyFunction( string_array const & inputPara,
                                              string_array const & componentNames,
                                              real64_array const & componentMolarWeight ):
  PVTFunction( inputPara[1], componentNames, componentMolarWeight )
{

  bool notFound = 1;

  for( localIndex i = 0; i < componentNames.size(); ++i )
  {

    if( streq( componentNames[i], "CO2" ) || streq( componentNames[i], "co2" ))
    {
      m_CO2Index = i;
      notFound = 0;
      break;
    }

  }

  GEOSX_ERROR_IF( notFound, "Component CO2 is not found!" );

  notFound = 1;

  for( localIndex i = 0; i < componentNames.size(); ++i )
  {

    if( streq( componentNames[i], "Water" ) || streq( componentNames[i], "water" ))
    {
      m_waterIndex = i;
      notFound = 0;
      break;
    }

  }

  GEOSX_ERROR_IF( notFound, "Component Water/Brine is not found!" );


  makeTable( inputPara );
  
  makeTable( inputPara );

}

void BrineEnthalpyFunction::makeTable( string_array const & inputPara )
{
  real64_array pressures;
  real64_array temperatures;

  real64 PStart, PEnd, dP;
  real64 TStart, TEnd, dT;
  real64 P, T, m;

  dT = -1.0;
  dP = -1.0;
  TStart = -1.0;
  TEnd = -1.0;
  PStart = -1.0;
  PEnd = -1.0;


  GEOSX_ERROR_IF( inputPara.size() < 9, "Invalid BrineEnthalpy input!" );

  try
  {

    PStart = stod( inputPara[2] );
    PEnd = stod( inputPara[3] );
    dP = stod( inputPara[4] );

    TStart = stod( inputPara[5] );
    TEnd = stod( inputPara[6] );
    dT = stod( inputPara[7] );

    m = stod( inputPara[8] );

  }
  catch( const std::invalid_argument & e )
  {

    GEOSX_ERROR( "Invalid BrineEnthalpy argument:" + std::string( e.what()));

  }

  P = PStart;

  while( P <= PEnd )
  {

    pressures.emplace_back( P );
    P += dP;

  }

  T = TStart;

  while( T <= TEnd )
  {

    temperatures.emplace_back( T );
    T += dT;

  }

  localIndex const nP = pressures.size();
  localIndex const nT = temperatures.size();

  real64_array2d enthalpies( nP, nT );

  calculateBrineEnthalpy( pressures, temperatures, m, enthalpies );

  m_BrineEnthalpyTable = std::make_shared< XYTable >( "BrineEnthalpyTable", pressures, temperatures, enthalpies );


  real64_array2d CO2Enthalpies( nP, nT );
  real64_array2d CO2Densities( nP, nT );

  SpanWagnerCO2DensityFunction::calculateCO2Density( pressures, temperatures, CO2Densities );

  CO2EnthalpyFunction::calculateCO2Enthalpy( pressures, temperatures, CO2Densities, CO2Enthalpies );

  m_CO2EnthalpyTable = std::make_shared< XYTable >( "CO2EnthalpyTable", pressures, temperatures, CO2Enthalpies );
  
}

void BrineEnthalpyFunction::calculateBrineEnthalpy( real64_array const & pressure, real64_array const & temperature, real64 const & m, real64_array2d const & enthalpy )
{

  GEOSX_UNUSED_VAR( pressure );
  
  static real64 a[4][3] = {
    {-9633.6, -4080.0, 286.49},
    {166.58, 68.577, -4.6856},
    {-0.90963, -0.36524, 0.0249667},
    {1.7965e-3, 7.1924e-4, -4.9e-5}
  };
    
  real64 x1, x2, h1, h2, dh, T;

  x1 = 1000.0 / (1000.0 + 58.44 * m);

  x2 = 58.44 * m  / (1000.0 + 58.44 * m);  

  for(localIndex ip = 0; ip < pressure.size(); ++ip)
  {  

    for(localIndex it = 0; it < temperature.size(); ++it)
    {  

      T = temperature[it];

      dh = 0.0;

      for(localIndex i = 0; i < 4; ++i)
        for(localIndex j = 0; j < 3; ++j)
        {  

          dh += a[i][j] * pow(T, real64(i)) * pow(m, real64(j));

        }

      dh *= 4.184 / (1000.0 + 58.44 * m);

    
      h1 = 0.12453e-4 * pow(T, 3.0) - 0.45137e-2 * pow(T, 2.0) + 4.81155 * T - 29.578;

      h2 = (-0.83624e-3 * pow(T, 3.0) + 0.16792 * pow(T, 2.0) - 25.9293 * T) * 4.184 / 58.44;

      enthalpy[ip][it] = (x1 * h1 + x2 * h2 + m * dh) * 1000.0;
  
    }
  }
}
  
void BrineEnthalpyFunction::evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition, EvalVarArgs & value, bool useMass ) const
{
  localIndex const numComponents = phaseComposition.size();  
  
  EvalArgs2D P, T, enthalpy, CO2Enthalpy;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  enthalpy = m_BrineEnthalpyTable->value( P, T );

  CO2Enthalpy = m_CO2EnthalpyTable->value( P, T );

  
  //assume there are only CO2 and brine here.
  
  EvalVarArgs enth1, enth2;

  enth1.m_var = CO2Enthalpy.m_var;
  enth1.m_der[0] = CO2Enthalpy.m_der[0];
  enth1.m_der[numComponents + 1] = CO2Enthalpy.m_der[1];    

  enth2.m_var = enthalpy.m_var;
  enth2.m_der[0] = enthalpy.m_der[0];
  enth2.m_der[numComponents + 1] = enthalpy.m_der[1];    

  real64 C = phaseComposition[m_waterIndex].m_var;

  real64 const waterMW = m_componentMolarWeight[m_waterIndex];   
  real64 const CO2MW = m_componentMolarWeight[m_CO2Index];

  if( useMass )
  {
    
    EvalVarArgs X = C * waterMW / (C * waterMW + (1.0 - C) * CO2MW);
    X.m_der[m_waterIndex+1] = 1.0;

    value = (1.0 - X ) * enth1 + X * enth2;

  }
  else 
  {

    EvalVarArgs X = C;
    X.m_der[m_waterIndex+1] = 1.0;

    value = (1.0 - X ) * enth1 / CO2MW + X * enth2 / waterMW;
  
  }  
    
}

REGISTER_CATALOG_ENTRY( PVTFunction,
                        BrineEnthalpyFunction,
                        string_array const &, string_array const &, real64_array const & )

} // namespace PVTProps
} // namespace geosx
