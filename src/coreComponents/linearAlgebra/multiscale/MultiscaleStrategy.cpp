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
 * @file MultiscaleStrategy.cpp
 */

#include "MultiscaleStrategy.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "MsRSBStrategy.hpp"

namespace geosx
{
namespace multiscale
{

template< typename LAI >
std::unique_ptr< MultiscaleStrategy< LAI > >
MultiscaleStrategy< LAI >::createInstance( string name,
                                           LinearSolverParameters::Multiscale params )
{
  switch( params.basisType )
  {
    case LinearSolverParameters::Multiscale::BasisType::msrsb:
    {
      return std::make_unique< MsRSBStrategy< LAI > >( std::move( name ), std::move( params ) );
    }
    default:
      GEOSX_ERROR( "Unsupported interpolation type: " << params.basisType );
  }
  return std::unique_ptr< MultiscaleStrategy< LAI > >();
}

template< typename LAI >
MultiscaleStrategy< LAI >::MultiscaleStrategy( string name,
                                               LinearSolverParameters::Multiscale params )
  : m_name( std::move( name ) ),
  m_params( std::move( params ) )
{}

template< typename LAI >
MultiscaleStrategy< LAI >::~MultiscaleStrategy() = default;

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class MultiscaleStrategy< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class MultiscaleStrategy< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class MultiscaleStrategy< PetscInterface >;
#endif

} // namespace multiscale
} // namespace geosx
