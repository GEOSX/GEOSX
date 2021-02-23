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


#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2INTERNALENERGYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2INTERNALENERGYFUNCTION_HPP_


#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class CO2InternalEnergyFunction : public PVTFunction
{
public:


  CO2InternalEnergyFunction( string_array const & inputPara,
                             string_array const & componentNames,
                             real64_array const & componentMolarWeight );

  ~CO2InternalEnergyFunction() override {}


  static constexpr auto m_catalogName = "CO2InternalEnergy";
  static string catalogName()                    { return m_catalogName; }
  virtual string getCatalogName() const override final { return catalogName(); }


  virtual PVTFuncType functionType() const override
  {
    return PVTFuncType::INTERNALENR;
  }

  virtual void evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;

  static void calculateCO2InternalEnergy( real64_array const & pressure, real64_array const & temperature, real64_array2d const & internalEnergy );

private:


};

}

}


#endif /* GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2INTERNALENERGYFUNCTION_HPP_ */
