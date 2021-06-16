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
 * @file permeabilitySelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYSELECTOR_HPP
#define GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "ConstantPermeability.hpp"
#include "CarmanKozenyPermeability.hpp"
#include "ParallelPlatesPermeability.hpp"
#include "StrainDependentPermeability.hpp"

namespace geosx
{

namespace constitutive
{

#ifndef PASSTHROUGH_HANDLE_CASE
#define PASSTHROUGH_HANDLE_CASE( MODEL ) \
  if( dynamicCast< MODEL * >( &perm ) ) \
  { \
    lambda( static_cast< MODEL & >( perm ) ); \
    return true; \
  }

template< typename LAMBDA >
void constitutiveUpdatePassThru( PermeabilityBase const & perm,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ConstantPermeability,
                               CarmanKozenyPermeability,
                               ParallelPlatesPermeability,
                               StrainDependentPermeability >::execute( perm, std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void constitutiveUpdatePassThru( PermeabilityBase & perm,
                                 LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ConstantPermeability,
                               CarmanKozenyPermeability,
                               ParallelPlatesPermeability,
                               StrainDependentPermeability >::execute( perm, std::forward< LAMBDA >( lambda ) );
}

#undef PASSTHROUGH_HANDLE_CASE
#endif

} // namespace constitutive

} // namespace geosx

#endif //GEOSX_CONSTITUTIVE_PERMEABILITY_PERMEABILITYSELECTOR_HPP
