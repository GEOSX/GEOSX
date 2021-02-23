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
 * @file BrineInternalEnergyFunction.cpp
 */

#include "constitutive/fluid/PVTFunctions/BrineInternalEnergyFunction.hpp"

namespace geosx
{

using namespace stringutilities;

namespace PVTProps
{

BrineInternalEnergyFunction::BrineInternalEnergyFunction( string_array const & inputPara,
                                                          string_array const & componentNames,
                                                          real64_array const & componentMolarWeight ):
  PVTFunction( inputPara[1], componentNames, componentMolarWeight )
{}

void BrineInternalEnergyFunction::evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature,
                                              arraySlice1d< EvalVarArgs const > const & phaseComposition, EvalVarArgs & value, bool useMass ) const
{
  GEOSX_UNUSED_VAR( phaseComposition );
  GEOSX_UNUSED_VAR( useMass );

  localIndex const numComponents = phaseComposition.size();

  value.m_var =  0.001 * pressure.m_var + 1.0 * temperature.m_var;
  // pressure derivative
  value.m_der[0] = 0.001;
  // temperature derivative
  value.m_der[numComponents+1] = 1.0;

}



REGISTER_CATALOG_ENTRY( PVTFunction,
                        BrineInternalEnergyFunction,
                        string_array const &, string_array const &, real64_array const & )

} // namespace PVTProps
} // namespace geosx
