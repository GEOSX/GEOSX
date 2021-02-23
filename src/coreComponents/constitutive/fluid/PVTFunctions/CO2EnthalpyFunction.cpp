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
 * @file CO2EnthalpyFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/CO2EnthalpyFunction.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

CO2EnthalpyFunction::CO2EnthalpyFunction( string_array const & inputPara,
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

  makeTable( inputPara );

} 

void CO2EnthalpyFunction::makeTable( string_array const & inputPara )
{

  real64_array pressures;
  real64_array temperatures;

  real64 PStart, PEnd, dP;
  real64 TStart, TEnd, dT;
  real64 P, T;

  dT = -1.0;
  dP = -1.0;
  TStart = -1.0;
  TEnd = -1.0;
  PStart = -1.0;
  PEnd = -1.0;

  GEOSX_ERROR_IF( inputPara.size() < 8, "Invalid CO2Enthalpy input!" );

  try
  {

    PStart = stod( inputPara[2] );
    PEnd = stod( inputPara[3] );
    dP = stod( inputPara[4] );

    TStart = stod( inputPara[5] );
    TEnd = stod( inputPara[6] );
    dT = stod( inputPara[7] );

  }
  catch( const std::invalid_argument & e )
  {

    GEOSX_ERROR( "Invalid CO2Enthalpy argument:" + std::string( e.what()));

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

  localIndex nP = pressures.size();
  localIndex nT = temperatures.size();

  real64_array2d enthalpies( nP, nT );
  real64_array2d densities( nP, nT );

  SpanWagnerCO2DensityFunction::calculateCO2Density( pressures, temperatures, densities );

  calculateCO2Enthalpy( pressures, temperatures, densities, enthalpies );

  m_CO2EnthalpyTable = std::make_shared< XYTable >( "CO2EnthalpyTable", pressures, temperatures, enthalpies );


}
  
void CO2EnthalpyFunction::calculateCO2Enthalpy( real64_array const & pressure, real64_array const & temperature, real64_array2d const & density, real64_array2d const & enthalpy )
{

    static const real64 dc = 467.6;
    static const real64 Tc = 304.128;
    static const real64 a0[8] = {8.37304456, -3.70454304, 2.500, 1.99427042, 0.62105248, 0.41195293, 1.04028922, 0.08327678};

    static const real64 theta0[8] = {0.0, 0.0, 0.0, 3.15163, 6.11190, 6.77708, 11.32384, 27.08792};

    static const real64 ni[42] = {0.38856823203161,   2.938547594274,    -5.5867188534934,   -0.76753199592477,
                              0.31729005580416,   0.54803315897767,   0.12279411220335,   2.165896154322,   
                              1.5841735109724,   -0.23132705405503,   0.058116916431436, -0.55369137205382, 
                              0.48946615909422,  -0.024275739843501,  0.062494790501678, -0.12175860225246,
                             -0.37055685270086,  -0.016775879700426 ,-0.11960736637987,  -0.045619362508778,
                              0.035612789270346, -0.0074427727132052,-0.0017395704902432,-0.021810121289527,
                              0.024332166559236, -0.037440133423463,  0.14338715756878,  -0.13491969083286,
                             -0.02315122505348,   0.012363125492901,  0.002105832197294, -0.00033958519026368,
                              0.0055993651771592,-0.00030335118055646, -213.6548868832,    26641.569149272,
                              -24027.212204557,   -283.41603423999,    212.47284400179,    -0.66642276540751,
                              0.72608632349897,   0.055068668612842};

    static const real64 di[42] = {1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 3.0, 1.0, 2.0, 4.0,
                              5.0, 5.0, 5.0, 6.0, 6.0, 6.0, 1.0, 1.0, 4.0, 4.0,
                              4.0, 7.0, 8.0, 2.0, 3.0, 3.0, 5.0, 5.0, 6.0, 7.0,
                              8.0, 10.0, 4.0, 8.0, 2.0, 2.0, 2.0, 3.0, 3.0, 0.0,  
                              0.0, 0.0};

    static const real64 ti[42] = {0.00,  0.75,  1.00,  2.00,  0.75,  2.00,  0.75,  1.50,  1.50,  2.50,
                              0.00,  1.50,  2.00,  0.00,  1.00,  2.00,  3.00,  6.00,  3.00,  6.00,
                              8.00,  6.00,  0.00,  7.00, 12.00, 16.00, 22.00, 24.00,  16.0, 24.00,
                              8.00,  2.00, 28.00, 14.00,  1.00,  0.00,  1.00,  3.00,  3.00,  0.00,
                              0.00,  0.00};
    static const real64 ci[42] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
                              1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0,
                              2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 4.0, 4.0,
                              4.0, 4.0, 5.0, 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                              0.0, 0.0};

    static const real64 alphai[42] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 25.0,25.0,25.0,15.0,20.0, 0.0,
                                  0.0, 0.0};

    static const real64 betai[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0,325.0,300.0,300.0,275.0,275.0, 0.3,
                                  0.3, 0.3};

    static const real64 gammai[42] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                  0.0, 0.0, 0.0, 0.0, 1.16, 1.19, 1.19, 1.25,1.22, 0.0,
                                  0.0, 0.0};
    static const real64 ei[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0,
                               0.0, 0.0};

    static const real64 bi[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.875,
                               0.925, 0.875};

    static const real64 Ai[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.7,
                               0.7, 0.7};

    static const real64 Bi[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.3,
                               0.3, 1.0};

    static const real64 Ci[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 10.0,
                               10.0, 12.5};

    static const real64 Di[42] =  {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 275.0,
                               275.0, 275.0};

    real64 rd, rt, den;
    real64 Tkelvin;

    real64 phi0t;
    real64 phirt, phird;

    real64 theta, delta, R, deltard;

    for(localIndex ip = 0; ip < pressure.size(); ++ip)
      for(localIndex it = 0; it < temperature.size(); ++it)
      {
    
        den = density[ip][it];
        Tkelvin = temperature[it] + 273.15;
        rd=den/dc;
        rt=Tc/Tkelvin;

        phi0t = a0[1] + a0[2] / rt;

        for(localIndex i=3; i<8; i++)
        {   
          phi0t += a0[i] * theta0[i] * (1.0 / (1.0 - exp(-rt * theta0[i])) - 1.0);

        }

        phird = 0.0;
        phirt = 0.0;

        for(localIndex i=0; i<7; i++)
        {

          phird += ni[i] * di[i] * pow(rd, di[i] - 1.0) * pow(rt, ti[i]);

          phirt += ni[i] * ti[i] * pow(rd, di[i]) * pow(rt, ti[i] - 1.0);

        }


        for(localIndex i=7; i<34; i++)
        {
          
          phird +=  ni[i]  * exp(-pow(rd, ci[i])) * pow(rd, di[i]-1.0) * pow(rt, ti[i]) * (di[i]-ci[i] * pow(rd, ci[i]));

          phirt += ni[i] * ti[i] * pow(rd, di[i]) * pow(rt, ti[i] - 1.0) * exp(-pow(rd, ci[i]));

        }


        for(localIndex i=34; i<39; i++)
        {
          
          phird += ni[i]  * pow(rd, di[i]) * pow(rt, ti[i]) * exp(-alphai[i] * pow(rd - ei[i], 2.0) - betai[i] * pow(rt - gammai[i], 2.0)) * (di[i] / rd - 2.0 * alphai[i] * (rd - ei[i]));

          phirt += ni[i] * pow(rd, di[i]) * pow(rt, ti[i]) * exp(-alphai[i] * pow(rd - ei[i], 2.0) - betai[i] * pow(rt - gammai[i], 2.0)) * (ti[i] / rt - 2.0 * betai[i] * (rt - gammai[i]));

        }

        for(localIndex i=39; i<42; i++)
        {
          
          theta = (1.0 - rt) + Ai[i] * pow(pow(rd - 1.0, 2.0), (1.0 / (2.0 * betai[i])));
          delta = pow(theta, 2.0) + Bi[i] * pow(rd - 1.0, 2.0*alphai[i]);
          R = exp(-Ci[i] * pow(rd - 1.0, 2.0) - Di[i] * pow(rt - 1.0, 2.0));
          deltard = (rd - 1.0) * (Ai[i] * theta * 2.0 / betai[i] * pow(pow(rd - 1.0, 2.0), 0.5 / betai[i] - 1.0) + 2.0 * Bi[i] * alphai[i] * pow(pow(rd - 1.0, 2.0), alphai[i] - 1.0));

          phird += ni[i] * (pow(delta, bi[i]) * (R + rd *(-2.0 *Ci[i] * (rd - 1.0) * R)) + bi[i] * pow(delta, bi[i] - 1.0) * deltard * rd * R);

          phirt += ni[i] * rd *(-2.0 * theta * bi[i] * pow(delta, bi[i] - 1.0) * R + pow(delta, bi[i]) * (-2.0 * Di[i] * R * (rt - 1.0)));

        }

        enthalpy[ip][it] = (1.0 + rt * (phi0t + phirt) + rd * phird) * Tkelvin * 188.924;

      }

}

  
void CO2EnthalpyFunction::evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition, EvalVarArgs & value, bool useMass ) const
{

  localIndex const numComponents = phaseComposition.size();

  EvalArgs2D P, T, enth;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  enth = m_CO2EnthalpyTable->value( P, T );

  if( !useMass )
  {
    real64 CO2MW = m_componentMolarWeight[m_CO2Index];

    enth /= CO2MW;
  }

  value.m_var = enth.m_var;
  value.m_der[0] = enth.m_der[0];  
  value.m_der[numComponents+1] = enth.m_der[1];

}

REGISTER_CATALOG_ENTRY( PVTFunction,
                        CO2EnthalpyFunction,
                        string_array const &, string_array const &, real64_array const & )

} // namespace PVTProps
} // namespace geosx
