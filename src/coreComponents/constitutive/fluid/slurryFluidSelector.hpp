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
 * @file slurryFluidSelector.hpp
 */
#ifndef GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDSELECTOR_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDSELECTOR_HPP_

#include "constitutive/ConstitutivePassThruHandler.hpp"
#include "constitutive/fluid/ProppantSlurryFluid.hpp"

namespace geosx
{
namespace constitutive
{
template< typename LAMBDA >
void
constitutiveUpdatePassThru( SlurryFluidBase const & fluid, LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ProppantSlurryFluid >::Execute(
    fluid,
    std::forward< LAMBDA >( lambda ) );
}

template< typename LAMBDA >
void
constitutiveUpdatePassThru( SlurryFluidBase & fluid, LAMBDA && lambda )
{
  ConstitutivePassThruHandler< ProppantSlurryFluid >::Execute(
    fluid,
    std::forward< LAMBDA >( lambda ) );
}

}  // namespace constitutive

}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_SLURRYFLUIDSELECTOR_HPP_
