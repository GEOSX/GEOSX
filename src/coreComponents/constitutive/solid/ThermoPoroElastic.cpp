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
 * @file ThermoPoroElastic.cpp
 */

#include "ThermoPoroElastic.hpp"

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
ThermoPoroElastic< BASE >::ThermoPoroElastic( string const & name, Group * const parent ):
  PoroElastic< BASE >( name, parent ),
  m_thermalStressCoefficient(),
  m_thermalPorosityCoefficient()  
{
  this->registerWrapper( viewKeyStruct::thermalStressCoefficientString(), &m_thermalStressCoefficient ).
    setApplyDefaultValue( 1e6 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Thermal stress coefficient" );

  this->registerWrapper( viewKeyStruct::thermalPorosityCoefficientString(), &m_thermalPorosityCoefficient ).
    setApplyDefaultValue( 0.0 ).
    setInputFlag( InputFlags::REQUIRED ).
    setDescription( "Thermal porosity coefficient" );
}

template< typename BASE >
ThermoPoroElastic< BASE >::~ThermoPoroElastic()
{}

typedef ThermoPoroElastic< ElasticIsotropic > ThermoPoroElasticIsotropic;
typedef ThermoPoroElastic< ElasticTransverseIsotropic > ThermoPoroElasticTransverseIsotropic;
typedef ThermoPoroElastic< DruckerPrager > ThermoPoroDruckerPrager;
typedef ThermoPoroElastic< DruckerPragerExtended > ThermoPoroDruckerPragerExtended;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, ThermoPoroDruckerPragerExtended, string const &, Group * const )

} /* namespace constitutive */
} /* namespace geosx */
