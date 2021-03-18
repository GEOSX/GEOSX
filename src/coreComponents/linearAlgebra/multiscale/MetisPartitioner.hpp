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
 * @file MetisPartitioner.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_METISPARTITIONER_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_METISPARTITIONER_HPP

#include "linearAlgebra/multiscale/PartitionerBase.hpp"

namespace geosx
{
namespace multiscale
{

class MetisPartitioner : public PartitionerBase
{
public:

  MetisPartitioner( LinearSolverParameters::Multiscale::Coarsening params );

  virtual localIndex generate( MultiscaleMeshLevel const & mesh,
                               arrayView1d< localIndex > const & partition ) override;

private:

  LinearSolverParameters::Multiscale::Coarsening m_params;

};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_METISPARTITIONER_HPP
