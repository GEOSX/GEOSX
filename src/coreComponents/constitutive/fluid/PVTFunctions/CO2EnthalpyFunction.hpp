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


#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPYFUNCTION_HPP_


#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class CO2EnthalpyFunction : public PVTFunction
{
public:


  CO2EnthalpyFunction( string_array const & inputPara,
                       string_array const & componentNames,
                       real64_array const & componentMolarWeight );

  ~CO2EnthalpyFunction() override {}


  static constexpr auto m_catalogName = "CO2Enthalpy";
  static string catalogName()                    { return m_catalogName; }
  virtual string getCatalogName() const override final { return catalogName(); }


  virtual PVTFuncType functionType() const override
  {
    return PVTFuncType::ENTHALPY;
  }

  virtual void evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;


  static void calculateCO2Enthalpy( real64_array const & pressure, real64_array const & temperature, real64_array2d const & density,  real64_array2d const & enthalpy );

  
private:

  void makeTable( string_array const & inputPara );

  
  TableFunctionPtr m_CO2EnthalpyTable;
  localIndex m_CO2Index;
  
  
};

}

}


#endif /* GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2ENTHALPYFUNCTION_HPP_ */
