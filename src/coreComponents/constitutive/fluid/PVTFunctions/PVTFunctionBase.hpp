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
 * @file PVTFunctionBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_

#include "constitutive/fluid/PVTFunctions/UtilityFunctions.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{
namespace PVTProps
{
enum class PVTFuncType
{
  UNKNOWN,
  DENSITY,
  VISCOSITY
};

class PVTFunction
{
public:
  PVTFunction( string const & name,
               string_array const & componentNames,
               real64_array const & componentMolarWeight ) :
    m_functionName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~PVTFunction()
  {}

  using CatalogInterface =
    dataRepository::CatalogInterface< PVTFunction,
                                      string_array const &,
                                      string_array const &,
                                      real64_array const & >;
  static typename CatalogInterface::CatalogType &
  GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }
  virtual string
  GetCatalogName() = 0;

  string const &
  FunctionName() const
  {
    return m_functionName;
  }

  virtual PVTFuncType
  FunctionType() const = 0;

  //phase density/viscosity
  //input: P, T, phaseCompFraction
  //output: phase density/viscoty

  virtual void
  Evaluation( EvalVarArgs const & pressure,
              EvalVarArgs const & temperature,
              arraySlice1d< EvalVarArgs const > const & phaseComposition,
              EvalVarArgs & value,
              bool useMass = 0 ) const = 0;

protected:
  string m_functionName;
  string_array m_componentNames;
  real64_array m_componentMolarWeight;
};

}  // namespace PVTProps

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_PVTFUNCTIONBASE_HPP_
