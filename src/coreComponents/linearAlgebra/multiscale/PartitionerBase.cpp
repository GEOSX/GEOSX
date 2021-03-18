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
 * @file PartitionerBase.cpp
 */

#include "PartitionerBase.hpp"

#include "linearAlgebra/multiscale/MetisPartitioner.hpp"

namespace geosx
{
namespace multiscale
{

std::unique_ptr< PartitionerBase >
PartitionerBase::create( LinearSolverParameters::Multiscale::Coarsening const params )
{
  switch( params.partitionType )
  {
    case LinearSolverParameters::Multiscale::Coarsening::PartitionType::metis:
    {
      return std::make_unique< MetisPartitioner >( params );
    }
    default:
    {
      GEOSX_THROW( "Not supported yet", std::runtime_error );
    }
  }
  return std::unique_ptr< PartitionerBase >();
}

} // namespace multiscale
} // namespace geosx
