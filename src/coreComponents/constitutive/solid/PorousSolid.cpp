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
 * @file PorousSolid.cpp
 */

#include "PorousSolid.hpp"
#include "ElasticIsotropic.hpp"
#include "ElasticTransverseIsotropic.hpp"
#include "DruckerPrager.hpp"
#include "DruckerPragerExtended.hpp"
#include "Damage.hpp"
#include "DamageSpectral.hpp"
#include "DamageVolDev.hpp"

namespace geosx
{

using namespace dataRepository;

namespace constitutive
{

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::PorousSolid( string const & name, Group * const parent ):
  CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >( name, parent )
{}

template< typename SOLID_TYPE >
PorousSolid< SOLID_TYPE >::~PorousSolid() = default;

// Register all PorousSolid model types.
typedef PorousSolid< ElasticIsotropic > PorousElasticIsotropic;
typedef PorousSolid< ElasticTransverseIsotropic > PorousElasticTransverseIsotropic;
typedef PorousSolid< DruckerPrager > PorousDruckerPrager;
typedef PorousSolid< DruckerPragerExtended > PorousDruckerPragerExtended;
typedef PorousSolid< Damage< ElasticIsotropic > > PorousDamageElasticIsotropic;
typedef PorousSolid< DamageSpectral< ElasticIsotropic > > PorousDamageSpectralElasticIsotropic;
typedef PorousSolid< DamageVolDev< ElasticIsotropic > > PorousDamageVolDevElasticIsotropic;

REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousElasticTransverseIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPrager, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDruckerPragerExtended, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageSpectralElasticIsotropic, string const &, Group * const )
REGISTER_CATALOG_ENTRY( ConstitutiveBase, PorousDamageVolDevElasticIsotropic, string const &, Group * const )

}
} /* namespace geosx */
