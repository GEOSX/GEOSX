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
 * @file ThermoElastic.cpp
 */

#include "ThermoElastic.hpp"

#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"

namespace geosx
{

using namespace dataRepository;
namespace constitutive
{

template< typename BASE >
ThermoElastic< BASE >::ThermoElastic( string const & name, Group * const parent ):
  BASE( name, parent ),
  m_thermalStressCoefficient()
{
  this->registerWrapper( viewKeyStruct::thermalStressCoefficientString(), &m_thermalStressCoefficient ).
    setApplyDefaultValue( 1.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Thermal stress coefficient" );
}

template< typename BASE >
ThermoElastic< BASE >::~ThermoElastic()
{}

typedef ThermoElastic< ElasticIsotropic > ThermoElasticIsotropic;
typedef ThermoElastic< ElasticTransverseIsotropic > ThermoElasticTransverseIsotropic;
typedef ThermoElastic< DruckerPrager > ThermoDruckerPrager;
typedef ThermoElastic< DruckerPragerExtended > ThermoDruckerPragerExtended;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoDruckerPragerExtended, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geosx */
