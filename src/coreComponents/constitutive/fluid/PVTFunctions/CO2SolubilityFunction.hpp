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
 * @file CO2SolubilityFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYFUNCTION_HPP_

#include "FlashModelBase.hpp"

namespace geosx
{
namespace PVTProps
{
class CO2SolubilityFunction : public FlashModel
{
public:
  CO2SolubilityFunction( string_array const & inputPara,
                         string_array const & phaseNames,
                         string_array const & componentNames,
                         real64_array const & componentMolarWeight );

  ~CO2SolubilityFunction() override
  {}

  static constexpr auto m_catalogName = "CO2Solubility";
  static string
  CatalogName()
  {
    return m_catalogName;
  }
  virtual string
  GetCatalogName() override final
  {
    return CatalogName();
  }

  virtual void
  Partition(
    EvalVarArgs const & pressure,
    EvalVarArgs const & temperature,
    arraySlice1d< EvalVarArgs const > const & compFraction,
    arraySlice1d< EvalVarArgs > const & phaseFraction,
    arraySlice2d< EvalVarArgs > const & phaseCompFraction ) const override;

private:
  void
  MakeTable( const string_array & inputPara );

  TableFunctionPtr m_CO2SolubilityTable;
  localIndex m_CO2Index;
  localIndex m_waterIndex;
  localIndex m_phaseGasIndex;
  localIndex m_phaseLiquidIndex;
};

}  // namespace PVTProps

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_CO2SOLUBILITYFUNCTION_HPP_
