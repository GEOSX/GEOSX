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
 * @file ExecutableGroup.cpp
 */

#include "ExecutableGroup.hpp"

namespace geosx
{
void
ExecutableGroup::SignalToPrepareForExecution(
  real64 const GEOSX_UNUSED_PARAM( time_n ),
  real64 const GEOSX_UNUSED_PARAM( dt ),
  integer const GEOSX_UNUSED_PARAM( cycle ),
  dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) )
{}

void
ExecutableGroup::Cleanup( real64 const GEOSX_UNUSED_PARAM( time_n ),
                          integer const GEOSX_UNUSED_PARAM( cycleNumber ),
                          integer const GEOSX_UNUSED_PARAM( eventCounter ),
                          real64 const GEOSX_UNUSED_PARAM( eventProgress ),
                          dataRepository::Group * GEOSX_UNUSED_PARAM( domain ) )
{}

}  // namespace geosx
