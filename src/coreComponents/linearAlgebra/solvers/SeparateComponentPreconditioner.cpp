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
 * @file SeparateComponentPreconditioner.cpp
 */

#include "SeparateComponentPreconditioner.hpp"

#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

namespace geosx
{
template< typename LAI >
SeparateComponentPreconditioner< LAI >::SeparateComponentPreconditioner(
  localIndex const numComp,
  std::unique_ptr< PreconditionerBase< LAI > > precond ) :
  Base(),
  m_numComp( numComp ),
  m_precond( std::move( precond ) )
{
  GEOSX_LAI_ASSERT_GT( m_numComp, 0 );
  GEOSX_LAI_ASSERT( m_precond );
}

template< typename LAI >
SeparateComponentPreconditioner< LAI >::~SeparateComponentPreconditioner() = default;

template< typename LAI >
void
SeparateComponentPreconditioner< LAI >::compute( Matrix const & mat,
                                                 DofManager const & dofManager )
{
  Base::compute( mat, dofManager );

  // TODO: if matrix structure hasn't changed, can just copy entries into existing m_matSC
  LAIHelperFunctions::SeparateComponentFilter( mat, m_matSC, m_numComp );
  m_precond->compute( m_matSC, dofManager );
}

template< typename LAI >
void
SeparateComponentPreconditioner< LAI >::apply( Vector const & src,
                                               Vector & dst ) const
{
  m_precond->apply( src, dst );
}

template< typename LAI >
void
SeparateComponentPreconditioner< LAI >::clear()
{
  Base::clear();
  m_precond->clear();
  m_matSC.reset();
}

// -----------------------
// Explicit Instantiations
// -----------------------
#ifdef GEOSX_USE_TRILINOS
template class SeparateComponentPreconditioner< TrilinosInterface >;
#endif

#ifdef GEOSX_USE_HYPRE
template class SeparateComponentPreconditioner< HypreInterface >;
#endif

#ifdef GEOSX_USE_PETSC
template class SeparateComponentPreconditioner< PetscInterface >;
#endif

}  // namespace geosx
