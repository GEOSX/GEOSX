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
 * @file MultiscaleStrategy.hpp
 */

#ifndef GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALESTRATEGY_HPP
#define GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALESTRATEGY_HPP

#include "common/DataTypes.hpp"
#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/common/LinearOperator.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "mesh/MeshLevel.hpp"

#include <memory>

namespace geosx
{
namespace multiscale
{

template< typename LAI >
class MultiscaleStrategy
{
public:

  /// Alias for vector type
  using Vector = typename LAI::ParallelVector;

  /// Alias for matrix type
  using Matrix = typename LAI::ParallelMatrix;

  /// Alias for operator type
  using Operator = LinearOperator< Vector >;

  static std::unique_ptr< MultiscaleStrategy< LAI > >
  createInstance( string name, LinearSolverParameters::Multiscale params );

  explicit MultiscaleStrategy( string name,
                               LinearSolverParameters::Multiscale params );

  virtual ~MultiscaleStrategy();

  virtual Operator const & prolongation() const = 0;

  virtual Operator const & restriction() const = 0;

  virtual Matrix const & matrix() const = 0;

  virtual void initializeFineLevel( MeshLevel & mesh,
                                    DofManager const & dofManager,
                                    string const & fieldName,
                                    MPI_Comm const & comm ) = 0;

  virtual void initializeCoarseLevel( MultiscaleStrategy< LAI > & fine ) = 0;

  virtual void compute( Matrix const & fine_mat ) = 0;

protected:

  string m_name;

  LinearSolverParameters::Multiscale m_params;
};

} // namespace multiscale
} // namespace geosx

#endif //GEOSX_LINEARALGEBRA_MULTISCALE_MULTISCALESTRATEGY_HPP
