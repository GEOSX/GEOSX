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
 * @file FlashModelBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_

#include "constitutive/fluid/PVTFunctions/UtilityFunctions.hpp"
#include "codingUtilities/StringUtilities.hpp"

namespace geosx
{
namespace PVTProps
{
class FlashModel
{
public:
  FlashModel( string const & name,
              string_array const & componentNames,
              real64_array const & componentMolarWeight ) :
    m_modelName( name ),
    m_componentNames( componentNames ),
    m_componentMolarWeight( componentMolarWeight )
  {}

  virtual ~FlashModel()
  {}

  using CatalogInterface =
    dataRepository::CatalogInterface< FlashModel,
                                      string_array const &,
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
  FlashModelName() const
  {
    return m_modelName;
  }

  //partition
  //input: P, T, totalCompFraction
  //output: phaseFraction, phaseCompFraction

  virtual void
  Partition(
    EvalVarArgs const & pressure,
    EvalVarArgs const & temperature,
    arraySlice1d< EvalVarArgs const > const & compFraction,
    arraySlice1d< EvalVarArgs > const & phaseFraction,
    arraySlice2d< EvalVarArgs > const & phaseCompFraction ) const = 0;

protected:
  string m_modelName;
  string_array m_componentNames;
  real64_array m_componentMolarWeight;
};

}  // namespace PVTProps

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_FLASHMODELBASE_HPP_
