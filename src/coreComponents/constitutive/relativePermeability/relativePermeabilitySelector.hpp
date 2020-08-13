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
 * @file relativePermeabilitySelector.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
#define GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/relativePermeability/BrooksCoreyRelativePermeability.hpp"
#include "constitutive/relativePermeability/BrooksCoreyBakerRelativePermeability.hpp"
#include "constitutive/relativePermeability/VanGenuchtenBakerRelativePermeability.hpp"

namespace geosx
{
namespace constitutive
{
#ifndef PASSTHROUGH_HANDLE_CASE
  #define PASSTHROUGH_HANDLE_CASE( MODEL )         \
    if( dynamicCast< MODEL * >( &relPerm ) )       \
    {                                              \
      lambda( static_cast< MODEL & >( relPerm ) ); \
      return true;                                 \
    }

template< typename LAMBDA >
void
constitutiveUpdatePassThru( RelativePermeabilityBase const & relPerm,
                            LAMBDA && lambda )
{
  ConstitutivePassThruHandler<
    BrooksCoreyRelativePermeability,
    BrooksCoreyBakerRelativePermeability,
    VanGenuchtenBakerRelativePermeability >::Execute( relPerm,
                                                      std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void
constitutiveUpdatePassThru( RelativePermeabilityBase & relPerm, LAMBDA && lambda )
{
  ConstitutivePassThruHandler<
    BrooksCoreyRelativePermeability,
    BrooksCoreyBakerRelativePermeability,
    VanGenuchtenBakerRelativePermeability >::Execute( relPerm,
                                                      std::forward< LAMBDA >( lambda ) );
}

  #undef PASSTHROUGH_HANDLE_CASE
#endif

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_RELATIVEPERMEABILITY_RELATIVEPERMEABILITYSELECTOR_HPP
