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
 * @file BrineViscosityFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineViscosityFunction.hpp"

namespace geosx
{
using namespace stringutilities;

namespace PVTProps
{
BrineViscosityFunction::BrineViscosityFunction(
  string_array const & inputPara,
  string_array const & componentNames,
  real64_array const & componentMolarWeight ) :
  PVTFunction( inputPara[1], componentNames, componentMolarWeight )
{
  MakeCoef( inputPara );
}

void
BrineViscosityFunction::MakeCoef( string_array const & inputPara )
{
  constexpr real64 a = 0.0816;
  constexpr real64 b = 0.0122;
  constexpr real64 c = 0.000128;
  constexpr real64 d = 0.000629;
  constexpr real64 k = -0.7;

  constexpr real64 waterVisc = 8.9e-4;  //at 25C

  real64 m = -1.0;

  GEOSX_ERROR_IF( inputPara.size() < 3, "Invalid BrineViscosity input!" );

  try
  {
    m = stod( inputPara[2] );
  }
  catch( const std::invalid_argument & e )
  {
    GEOSX_ERROR( "Invalid BrineViscosity argument:" + std::string( e.what() ) );
  }

  m_coef0 = ( 1.0 + a * m + b * m * m + c * m * m * m ) * waterVisc;

  m_coef1 = d * ( 1.0 - exp( k * m ) ) * waterVisc;
}

void
BrineViscosityFunction::Evaluation(
  EvalVarArgs const & GEOSX_UNUSED_PARAM( pressure ),
  EvalVarArgs const & temperature,
  arraySlice1d< EvalVarArgs const > const & GEOSX_UNUSED_PARAM( phaseComposition ),
  EvalVarArgs & value,
  bool GEOSX_UNUSED_PARAM( useMass ) ) const
{
  value = m_coef0 + m_coef1 * temperature;
}

REGISTER_CATALOG_ENTRY( PVTFunction,
                        BrineViscosityFunction,
                        string_array const &,
                        string_array const &,
                        real64_array const & )

}  // namespace PVTProps

}  // namespace geosx
