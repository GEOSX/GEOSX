/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file BrineCO2DensityFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineCO2DensityFunction.hpp"

namespace geosx
{
using namespace stringutilities;

namespace PVTProps
{
BrineCO2DensityFunction::BrineCO2DensityFunction(
  string_array const & inputPara,
  string_array const & componentNames,
  real64_array const & componentMolarWeight ) :
  PVTFunction( inputPara[1], componentNames, componentMolarWeight )
{
  bool notFound = 1;

  for( localIndex i = 0; i < componentNames.size(); ++i )
  {
    if( streq( componentNames[i], "CO2" ) || streq( componentNames[i], "co2" ) )
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
    if( streq( componentNames[i], "Water" ) || streq( componentNames[i], "water" ) )
    {
      m_waterIndex = i;
      notFound = 0;
      break;
    }
  }

  GEOSX_ERROR_IF( notFound, "Component Water/Brine is not found!" );

  MakeTable( inputPara );
}

void
BrineCO2DensityFunction::MakeTable( string_array const & inputPara )
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

  GEOSX_ERROR_IF( inputPara.size() < 9, "Invalid BrineCO2Density input!" );

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
    GEOSX_ERROR( "Invalid BrineCO2Density argument:" + std::string( e.what() ) );
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

  real64_array2d densities( nP, nT );

  CalculateBrineDensity( pressures, temperatures, m, densities );

  m_BrineDensityTable = std::make_shared< XYTable >( "BrineDensityTable",
                                                     pressures,
                                                     temperatures,
                                                     densities );
}

void
BrineCO2DensityFunction::Evaluation(
  EvalVarArgs const & pressure,
  EvalVarArgs const & temperature,
  arraySlice1d< EvalVarArgs const > const & phaseComposition,
  EvalVarArgs & value,
  bool useMass ) const
{
  EvalArgs2D P, T, density;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  density = m_BrineDensityTable->Value( P, T );

  constexpr real64 a = 37.51;
  constexpr real64 b = -9.585e-2;
  constexpr real64 c = 8.740e-4;
  constexpr real64 d = -5.044e-7;

  real64 temp = T.m_var;

  real64 const V =
    ( a + b * temp + c * temp * temp + d * temp * temp * temp ) * 1e-6;

  real64 const CO2MW = m_componentMolarWeight[m_CO2Index];
  real64 const waterMW = m_componentMolarWeight[m_waterIndex];

  EvalVarArgs den, C, X;

  den.m_var = density.m_var;
  den.m_der[0] = density.m_der[0];

  X = phaseComposition[m_CO2Index];

  C = X * den / ( waterMW * ( 1.0 - X ) );

  if( useMass )
  {
    value = den + CO2MW * C - C * den * V;
  }
  else
  {
    value = den / waterMW + C - C * den * V / waterMW;
  }
}

void
BrineCO2DensityFunction::CalculateBrineDensity( real64_array const & pressure,
                                                real64_array const & temperature,
                                                real64 const & salinity,
                                                real64_array2d const & density )
{
  constexpr real64 c1 = -9.9595;
  constexpr real64 c2 = 7.0845;
  constexpr real64 c3 = 3.9093;

  constexpr real64 a1 = -0.004539;
  constexpr real64 a2 = -0.0001638;
  constexpr real64 a3 = 0.00002551;

  constexpr real64 AA = -3.033405;
  constexpr real64 BB = 10.128163;
  constexpr real64 CC = -8.750567;
  constexpr real64 DD = 2.663107;

  real64 P, x;

  for( localIndex i = 0; i < pressure.size(); ++i )
  {
    P = pressure[i] / 1e5;

    for( localIndex j = 0; j < temperature.size(); ++j )
    {
      x = c1 * exp( a1 * salinity ) + c2 * exp( a2 * temperature[j] ) +
        c3 * exp( a3 * P );

      density[i][j] = ( AA + BB * x + CC * x * x + DD * x * x * x ) * 1000.0;
    }
  }
}

REGISTER_CATALOG_ENTRY( PVTFunction,
                        BrineCO2DensityFunction,
                        string_array const &,
                        string_array const &,
                        real64_array const & )

}  // namespace PVTProps
}  // namespace geosx
