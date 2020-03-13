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
 * @file FenghourCO2ViscosityFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITYFUNCTION_HPP_

#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class FenghourCO2ViscosityFunction : public PVTFunction
{
public:

  FenghourCO2ViscosityFunction( string_array const & inputPara,
                                string_array const & componentNames,
                                real64_array const & componentMolarWeight );
  ~FenghourCO2ViscosityFunction() override
  {}

  static constexpr auto m_catalogName = "FenghourCO2Viscosity";
  static string CatalogName()                    { return m_catalogName; }
  virtual string GetCatalogName() override final { return CatalogName(); }

  virtual PVTFuncType FunctionType() const override
  {
    return PVTFuncType::VISCOSITY;

  }

  virtual void Evaluation( EvalVarArgs const & pressure,
                           EvalVarArgs const & temperature,
                           arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;


private:

  void MakeTable( string_array const & inputPara );

  void CalculateCO2Viscosity( real64_array const & pressure, real64_array const & temperature, real64_array2d const & density,
                              real64_array2d const & viscosity );

  void FenghourCO2Viscosity( real64 const & Tcent, real64 const & den, real64 & vis );

  TableFunctionPtr m_CO2ViscosityTable;
};

}

}
#endif //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FENGHOURCO2VISCOSITYFUNCTION_HPP_
