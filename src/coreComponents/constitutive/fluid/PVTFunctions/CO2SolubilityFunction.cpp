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
 * @file CO2SolubilityFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/CO2SolubilityFunction.hpp"

namespace geosx
{
using namespace stringutilities;

namespace PVTProps
{
constexpr real64 minForDivision = 1e-10;

constexpr real64 T_K_f = 273.15;
constexpr real64 P_Pa_f = 1e+5;
constexpr real64 P_c = 73.773 * P_Pa_f;
constexpr real64 T_c = 304.1282;
constexpr real64 Rgas = 8.314467;
constexpr real64 V_c = Rgas * T_c / P_c;

constexpr real64 acoef[] = { 8.99288497e-2,
                             -4.94783127e-1,
                             4.77922245e-2,
                             1.03808883e-2,
                             -2.82516861e-2,
                             9.49887563e-2,
                             5.20600880e-4,
                             -2.93540971e-4,
                             -1.77265112e-3,
                             -2.51101973e-5,
                             8.93353441e-5,
                             7.88998563e-5,
                             -1.66727022e-2,
                             1.398,
                             2.96e-2 };

real64
ff( real64 const & T, real64 const & P, real64 const & V_r )
{
  real64 const P_r = P * P_Pa_f / P_c;
  real64 const T_r = ( T_K_f + T ) / T_c;

  real64 const f_Z = 1.0 +
    ( acoef[0] + acoef[1] / ( T_r * T_r ) + acoef[2] / ( T_r * T_r * T_r ) ) / V_r +
    ( acoef[3] + acoef[4] / ( T_r * T_r ) + acoef[5] / ( T_r * T_r * T_r ) ) /
      ( V_r * V_r ) +
    ( acoef[6] + acoef[7] / ( T_r * T_r ) + acoef[8] / ( T_r * T_r * T_r ) ) /
      ( V_r * V_r * V_r * V_r ) +
    ( acoef[9] + acoef[10] / ( T_r * T_r ) + acoef[11] / ( T_r * T_r * T_r ) ) /
      ( V_r * V_r * V_r * V_r * V_r ) +
    acoef[12] / ( T_r * T_r * T_r ) / ( V_r * V_r ) *
      ( acoef[13] + acoef[14] / ( V_r * V_r ) ) * exp( -acoef[14] / ( V_r * V_r ) ) -
    P_r * V_r / T_r;

  return f_Z;
}

real64
PWater( real64 const & T )
{
  constexpr real64 ccoef[] = { -38.640844, 5.8948420, 59.876516, 26.654627, 10.637097 };

  real64 const P_c_w = 220.85;  // H2O critical pressure (bars)
  real64 const T_c_w = 647.29;  // H2O critical temperature (K)
  real64 const tt = ( ( T + T_K_f ) - T_c_w ) / T_c_w;
  real64 const x = ( P_c_w * ( T + T_K_f ) / T_c_w ) *
    ( 1 + ccoef[0] * pow( -tt, 1.9 ) + ccoef[1] * tt + ccoef[2] * tt * tt +
      ccoef[3] * tt * tt * tt + ccoef[4] * tt * tt * tt * tt );

  return x;
}

real64
logF( real64 const & T, real64 const & P, real64 const & V_r )
{
  real64 const P_r = P * P_Pa_f / P_c;
  real64 const T_r = ( T_K_f + T ) / T_c;

  real64 const Z = P_r * V_r / T_r;

  real64 const log_f = Z - 1 - log( Z ) +
    ( acoef[0] + acoef[1] / T_r / T_r + acoef[2] / T_r / T_r / T_r ) / V_r +
    ( acoef[3] + acoef[4] / T_r / T_r + acoef[5] / T_r / T_r / T_r ) / 2.0 / V_r /
      V_r +
    ( acoef[6] + acoef[7] / T_r / T_r + acoef[8] / T_r / T_r / T_r ) / 4.0 / V_r /
      V_r / V_r / V_r +
    ( acoef[9] + acoef[10] / T_r / T_r + acoef[11] / T_r / T_r / T_r ) / 5.0 /
      V_r / V_r / V_r / V_r / V_r +
    acoef[12] / 2.0 / T_r / T_r / T_r / acoef[14] *
      ( acoef[13] + 1.0 -
        ( acoef[13] + 1.0 + acoef[14] / V_r / V_r ) * exp( -acoef[14] / V_r / V_r ) );

  return log_f;
}

real64
Par( real64 const & T, real64 const & P, real64 const * cc )
{
  real64 x = cc[0] + cc[1] * T + cc[2] / T + cc[3] * T * T + cc[4] / ( 630.0 - T ) +
    cc[5] * P + cc[6] * P * log( T ) + cc[7] * P / T + cc[8] * P / ( 630.0 - T ) +
    cc[9] * P * P / ( 630.0 - T ) / ( 630.0 - T ) + cc[10] * T * log( P );

  return x;
}

void
CO2Solubility( real64 const & T,
               real64 const & P,
               real64 & V_r,
               real64 ( *f )( real64 const & x1, real64 const & x2, real64 const & x3 ) )
{
  constexpr real64 eps = 1e-9;
  int count = 0;
  constexpr real64 dx = 1e-10;

  real64 dre;
  real64 Vr_int = 0.05;

  V_r = 0.75 * Rgas * ( T_K_f + T ) / ( P * P_Pa_f ) * ( 1 / V_c );

  real64 v1, v0;

  for( ;; )
  {
    if( V_r < 0.0 )
    {
      V_r = Vr_int;
      Vr_int += 0.05;
    }

    v0 = ( *f )( T, P, V_r );
    v1 = ( *f )( T, P, V_r + dx );
    dre = -v0 / ( ( v1 - v0 ) / dx );

    if( fabs( dre ) < eps ) break;

    GEOSX_ERROR_IF( count > 50,
                    "CO2Solubiltiy NR convergence fails! "
                      << "dre = " << dre << ", eps = " << eps );

    count++;

    V_r += dre;
  }
}

void
CalculateCO2Solubility( real64_array const & pressure,
                        real64_array const & temperature,
                        real64 const & salinity,
                        real64_array2d const & solubiltiy )
{
  real64 T, P, V_r, m, logK, y_CO2;

  constexpr real64 mu[] = { 28.9447706,
                            -0.0354581768,
                            -4770.67077,
                            1.02782768e-5,
                            33.8126098,
                            9.04037140e-3,
                            -1.14934031e-3,
                            -0.307405726,
                            -0.0907301486,
                            9.32713393e-4,
                            0 };

  constexpr real64 lambda[] = { -0.411370585,
                                6.07632013e-4,
                                97.5347708,
                                0,
                                0,
                                0,
                                0,
                                -0.0237622469,
                                0.0170656236,
                                0,
                                1.41335834e-5 };

  constexpr real64 zeta[] = { 3.36389723e-4,
                              -1.98298980e-5,
                              0,
                              0,
                              0,
                              0,
                              0,
                              2.12220830e-3,
                              -5.24873303e-3,
                              0,
                              0 };

  m = salinity;

  for( localIndex i = 0; i < pressure.size(); ++i )
  {
    P = pressure[i] / P_Pa_f;

    for( localIndex j = 0; j < temperature.size(); ++j )
    {
      T = temperature[j];

      CO2Solubility( T, P, V_r, &ff );

      logK = Par( T + T_K_f, P, mu ) - logF( T, P, V_r ) +
        2 * Par( T + T_K_f, P, lambda ) * m + Par( T + T_K_f, P, zeta ) * m * m;

      y_CO2 = ( P - PWater( T ) ) / P;

      solubiltiy[i][j] = y_CO2 * P / exp( logK );
    }
  }
}

CO2SolubilityFunction::CO2SolubilityFunction( string_array const & inputPara,
                                              string_array const & phaseNames,
                                              string_array const & componentNames,
                                              real64_array const & componentMolarWeight ) :
  FlashModel( inputPara[1], componentNames, componentMolarWeight )
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

  notFound = 1;

  for( localIndex i = 0; i < phaseNames.size(); ++i )
  {
    if( streq( phaseNames[i], "CO2" ) || streq( phaseNames[i], "co2" ) ||
        streq( phaseNames[i], "gas" ) || streq( phaseNames[i], "Gas" ) )
    {
      m_phaseGasIndex = i;
      notFound = 0;
      break;
    }
  }

  GEOSX_ERROR_IF( notFound, "Phase co2/gas is not found!" );

  notFound = 1;

  for( localIndex i = 0; i < phaseNames.size(); ++i )
  {
    if( streq( phaseNames[i], "Water" ) || streq( phaseNames[i], "water" ) ||
        streq( phaseNames[i], "Liquid" ) || streq( phaseNames[i], "liquid" ) )
    {
      m_phaseLiquidIndex = i;
      notFound = 0;
      break;
    }
  }

  GEOSX_ERROR_IF( notFound, "Phase water/liquid is not found!" );

  MakeTable( inputPara );
}

void
CO2SolubilityFunction::MakeTable( string_array const & inputPara )
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

  GEOSX_ERROR_IF( inputPara.size() < 9, "Invalid CO2Solubility input!" );

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
    GEOSX_ERROR( "Invalid CO2Solubility argument:" + std::string( e.what() ) );
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

  real64_array2d solubilities( nP, nT );

  CalculateCO2Solubility( pressures, temperatures, m, solubilities );

  m_CO2SolubilityTable = std::make_shared< XYTable >( "CO2SolubilityTable",
                                                      pressures,
                                                      temperatures,
                                                      solubilities );
}

void
CO2SolubilityFunction::Partition(
  EvalVarArgs const & pressure,
  EvalVarArgs const & temperature,
  arraySlice1d< EvalVarArgs const > const & compFraction,
  arraySlice1d< EvalVarArgs > const & phaseFraction,
  arraySlice2d< EvalVarArgs > const & phaseCompFraction ) const
{
  EvalArgs2D P, T, solubility;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  //solubiltiy mol/kg(water)  X = Csat/W
  solubility = m_CO2SolubilityTable->Value( P, T );

  real64 const waterMW = m_componentMolarWeight[m_waterIndex];

  solubility *= waterMW;

  EvalVarArgs X, Y;

  X.m_var = solubility.m_var;
  X.m_der[0] = solubility.m_der[0];

  //Y = C/W = z/(1-z)

  if( compFraction[m_CO2Index].m_var > 1.0 - minForDivision )
  {
    Y = compFraction[m_CO2Index] / minForDivision;
  }
  else
  {
    Y = compFraction[m_CO2Index] / ( 1.0 - compFraction[m_CO2Index] );
  }

  if( Y < X )
  {
    //liquid phase only

    phaseFraction[m_phaseLiquidIndex] = 1.0;
    phaseFraction[m_phaseGasIndex] = 0.0;

    for( localIndex c = 0; c < m_componentNames.size(); ++c )
    {
      phaseCompFraction[m_phaseLiquidIndex][c] = compFraction[c];
    }
  }
  else
  {
    // two-phase
    // liquid phase fraction = (Csat + W) / (C + W) = (Csat/W + 1) / (C/W + 1)

    phaseFraction[m_phaseLiquidIndex] = ( X + 1.0 ) / ( Y + 1.0 );
    phaseFraction[m_phaseGasIndex] = 1.0 - phaseFraction[m_phaseLiquidIndex];

    //liquid phase composition  CO2 = Csat / (Csat + W) = (Csat/W) / (Csat/W + 1)

    phaseCompFraction[m_phaseLiquidIndex][m_CO2Index] = X / ( X + 1.0 );
    phaseCompFraction[m_phaseLiquidIndex][m_waterIndex] =
      1.0 - phaseCompFraction[m_phaseLiquidIndex][m_CO2Index];

    //gas phase composition  CO2 = 1.0

    phaseCompFraction[m_phaseGasIndex][m_CO2Index] = 1.0;
    phaseCompFraction[m_phaseGasIndex][m_waterIndex] = 0.0;
  }
}

REGISTER_CATALOG_ENTRY( FlashModel,
                        CO2SolubilityFunction,
                        string_array const &,
                        string_array const &,
                        string_array const &,
                        real64_array const & )

}  // namespace PVTProps

}  // namespace geosx
