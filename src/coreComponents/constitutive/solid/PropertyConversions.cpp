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
 *  @file PropertyConversions.cpp
 */

#include "PropertyConversions.hpp"

namespace geosx
{
namespace constitutive
{

PoissonRatio::PoissonRatio( BulkModulus K, ShearModulus G )
{
  PoissonRatio::value = ( 3.0 * K.value - 2.0 * G.value ) / ( 6.0 * K.value + 2.0 * G.value );
}

PoissonRatio::PoissonRatio( BulkModulus K, YoungModulus E )
{
  PoissonRatio::value = ( 3.0 * K.value - E.value ) / ( 6.0 * K.value);
}

PoissonRatio::PoissonRatio( ShearModulus G, YoungModulus E )
{
  PoissonRatio::value = 0.5 * E.value / G.value - 1.0;
}


} /* namespace constitutive */

} /* namespace geosx */
