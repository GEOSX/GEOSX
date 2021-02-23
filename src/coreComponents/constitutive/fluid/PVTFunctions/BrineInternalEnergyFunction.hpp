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


#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGYFUNCTION_HPP_


#include "PVTFunctionBase.hpp"

namespace geosx
{

namespace PVTProps
{

class BrineInternalEnergyFunction : public PVTFunction
{
public:


  BrineInternalEnergyFunction( string_array const & inputPara,
                               string_array const & componentNames,
                               real64_array const & componentMolarWeight );

  ~BrineInternalEnergyFunction() override {}


  static constexpr auto m_catalogName = "BrineInternalEnergy";
  static string catalogName()                    { return m_catalogName; }
  virtual string getCatalogName() const override final { return catalogName(); }


  virtual PVTFuncType functionType() const override
  {
    return PVTFuncType::INTERNALENR;
  }

  virtual void evaluation( EvalVarArgs const & pressure, EvalVarArgs const & temperature, arraySlice1d< EvalVarArgs const > const & phaseComposition,
                           EvalVarArgs & value, bool useMass = 0 ) const override;

  static void calculateBrineInternalEnergy( real64_array const & pressure, real64_array const & temperature, real64_array2d const & internalEnergy );

private:

};

}

}


#endif /* GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_BRINEINTERNALENERGYFUNCTION_HPP_ */
