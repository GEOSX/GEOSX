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
 * @file PartitionerBase.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP

#include "linearAlgebra/multiscale/MultiscaleMeshLevel.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <memory>

namespace geosx
{
namespace multiscale
{

class PartitionerBase
{
public:

  static std::unique_ptr< PartitionerBase >
  create( LinearSolverParameters::Multiscale::Coarsening params );

  virtual ~PartitionerBase() = default;

  virtual localIndex generate( MultiscaleMeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) = 0;

};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_PARTITIONERBASE_HPP
