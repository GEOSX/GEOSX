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
 * @file FenghourCO2ViscosityFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/FenghourCO2ViscosityFunction.hpp"
#include "constitutive/fluid/PVTFunctions/SpanWagnerCO2DensityFunction.hpp"

namespace geosx
{
namespace PVTProps
{
FenghourCO2ViscosityFunction::FenghourCO2ViscosityFunction(
  string_array const & inputPara,
  string_array const & componentNames,
  real64_array const & componentMolarWeight ) :
  PVTFunction( inputPara[1], componentNames, componentMolarWeight )
{
  MakeTable( inputPara );
}

void
FenghourCO2ViscosityFunction::MakeTable( string_array const & inputPara )
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

  GEOSX_ERROR_IF( inputPara.size() < 8, "Invalid FenghourCO2Viscosity input!" );

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
    GEOSX_ERROR( "Invalid FenghourCO2Viscosity argument:" + std::string( e.what() ) );
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

  real64_array2d viscosities( nP, nT );
  real64_array2d densities( nP, nT );

  SpanWagnerCO2DensityFunction::CalculateCO2Density( pressures,
                                                     temperatures,
                                                     densities );

  CalculateCO2Viscosity( pressures, temperatures, densities, viscosities );

  m_CO2ViscosityTable = std::make_shared< XYTable >( "FenghourCO2ViscosityTable",
                                                     pressures,
                                                     temperatures,
                                                     viscosities );
}

void
FenghourCO2ViscosityFunction::Evaluation(
  EvalVarArgs const & pressure,
  EvalVarArgs const & temperature,
  arraySlice1d< EvalVarArgs const > const & GEOSX_UNUSED_PARAM( phaseComposition ),
  EvalVarArgs & value,
  bool GEOSX_UNUSED_PARAM( useMass ) ) const
{
  EvalArgs2D P, T, viscosity;
  P.m_var = pressure.m_var;
  P.m_der[0] = 1.0;

  T.m_var = temperature.m_var;
  T.m_der[1] = 1.0;

  viscosity = m_CO2ViscosityTable->Value( P, T );

  value.m_var = viscosity.m_var;
  value.m_der[0] = viscosity.m_der[0];
}

void
FenghourCO2ViscosityFunction::FenghourCO2Viscosity( real64 const & Tcent,
                                                    real64 const & den,
                                                    real64 & vis )
{
  constexpr real64 espar = 251.196;
  constexpr real64 esparInv = 1.0 / espar;
  constexpr real64 aa[5] = { 0.235156,
                             -0.491266,
                             5.211155e-2,
                             5.347906e-2,
                             -1.537102e-2 };
  constexpr real64 d11 = 0.4071119e-2;
  constexpr real64 d21 = 0.7198037e-4;
  constexpr real64 d64 = 0.2411697e-16;
  constexpr real64 d81 = 0.2971072e-22;
  constexpr real64 d82 = -0.1627888e-22;

  // temperature in Kelvin
  const real64 Tkelvin = Tcent + 273.15;
  // evaluate vlimit from eqns 3-5
  const real64 Tred = Tkelvin * esparInv;
  const real64 x = log( Tred );
  const real64 lnGfun =
    aa[0] + x * ( aa[1] + x * ( aa[2] + x * ( aa[3] + x * aa[4] ) ) );
  const real64 GfunInv = exp( -lnGfun );
  const real64 vlimit = 1.00697 * sqrt( Tkelvin ) * GfunInv;

  const real64 d2 = den * den;
  const real64 vxcess = den *
    ( d11 +
      den *
        ( d21 + d2 * d2 * ( d64 / ( Tred * Tred * Tred ) + d2 * ( d81 + d82 / Tred ) ) ) );

  constexpr real64 vcrit = 0.0;

  vis = 1e-6 * ( vlimit + vxcess + vcrit );
}

void
FenghourCO2ViscosityFunction::CalculateCO2Viscosity(
  real64_array const & pressure,
  real64_array const & temperature,
  real64_array2d const & density,
  real64_array2d const & viscosity )
{
  for( localIndex i = 0; i < pressure.size(); ++i )
  {
    for( localIndex j = 0; j < temperature.size(); ++j )
    {
      FenghourCO2Viscosity( temperature[j], density[i][j], viscosity[i][j] );
    }
  }
}

REGISTER_CATALOG_ENTRY( PVTFunction,
                        FenghourCO2ViscosityFunction,
                        string_array const &,
                        string_array const &,
                        real64_array const & )

}  // namespace PVTProps

}  // namespace geosx
