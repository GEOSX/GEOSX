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
 * @file PetscPreconditioner.cpp
 */

#include "PetscPreconditioner.hpp"

#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <petscksp.h>

#include <fenv.h>

namespace geosx
{


PetscPreconditioner::PetscPreconditioner( LinearSolverParameters params )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{}
{ }

PetscPreconditioner::~PetscPreconditioner()
{
  clear();
}

void CreatePetscAMG( LinearSolverParameters const & params, PC precond )
{
  // Default options only for the moment
  GEOSX_LAI_CHECK_ERROR( PCSetType( precond, PCGAMG ) );
  GEOSX_UNUSED_VAR( params )

  // TODO: need someone familiar with PETSc to take a look at this
#if 0
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec, PCHMG ) );
  GEOSX_LAI_CHECK_ERROR( PCHMGSetInnerPCType( prec, PCGAMG ) );

  // Set maximum number of multigrid levels
  if( m_parameters.amg.maxLevels > 0 )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetLevels( prec,
                                          LvArray::integerConversion< PetscInt >( m_parameters.amg.maxLevels ),
                                          nullptr ) );
  }

  // Set the number of sweeps
  if( m_parameters.amg.numSweeps > 1 )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetNumberSmooth( prec,
                                                LvArray::integerConversion< PetscInt >( m_parameters.amg.numSweeps ) ) );
  }

  // Set type of cycle (1: V-cycle (default); 2: W-cycle)
  if( m_parameters.amg.cycleType == "V" )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetCycleType( prec, PC_MG_CYCLE_V ) );
  }
  else if( m_parameters.amg.cycleType == "W" )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetCycleType( prec, PC_MG_CYCLE_W ) );
  }

  // Set smoother to be used (for all levels)
  PetscInt numLevels;
  PetscInt l;
  KSP smoother;
  PC smootherPC;
  GEOSX_LAI_CHECK_ERROR( PCMGGetLevels( prec, &numLevels ) );

  GEOSX_LOG_RANK_VAR( numLevels );

  for( l = 0; l < numLevels; ++l )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGGetSmoother( prec, l, &smoother ) );
    GEOSX_LAI_CHECK_ERROR( KSPSetType( smoother, KSPRICHARDSON ) );
    GEOSX_LAI_CHECK_ERROR( KSPGetPC( smoother, &smootherPC ) );

    if( m_parameters.amg.smootherType == "jacobi" )
    {
      GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCJACOBI ) );
    }
    else if( m_parameters.amg.smootherType == "gaussSeidel" )
    {
      GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCSOR ) );
    }
    else if( m_parameters.amg.smootherType.substr( 0, 3 ) == "ilu" )
    {
      // Set up additive Schwartz preconditioner
      GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCASM ) );
      GEOSX_LAI_CHECK_ERROR( PCASMSetOverlap( smootherPC, 0 ) );
      GEOSX_LAI_CHECK_ERROR( PCASMSetType( smootherPC, PC_ASM_RESTRICT ) );
      // GEOSX_LAI_CHECK_ERROR( PCSetUp( smootherPC ) );

      // Get local preconditioning context
      KSP * ksp_local;
      PetscInt n_local, first_local;
      GEOSX_LAI_CHECK_ERROR( PCASMGetSubKSP( smootherPC, &n_local, &first_local, &ksp_local ) );

      // Sanity checks
      GEOSX_LAI_ASSERT_EQ( n_local, 1 );
      GEOSX_LAI_ASSERT_EQ( first_local, MpiWrapper::Comm_rank( MPI_COMM_GEOSX ) );

      // Set up local block ILU preconditioner
      PC prec_local;
      GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp_local[0], KSPPREONLY ) );
      GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp_local[0], &prec_local ) );
      GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, PCILU ) );
      if( m_parameters.amg.smootherType == "ilu1" )
      {
        GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, 1 ) );
      }
      else
      {
        GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, 0 ) );
      }
      // GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( smootherPC ) );
    }
  }

  // Set coarsest level solver
  // TODO
  // m_parameters.amg.coarseType



  // Set aggretation threshold
  // TODO
  // m_parameters.amg.aggregationThreshold

  // TODO: add user-defined null space / rigid body mode support
  // if ( m_parameters.amg.nullSpaceType )
  // {
  //   ...
  // }
#endif
}

PCType getPetscSmootherType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, PCType > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, PCILU },
    { LinearSolverParameters::PreconditionerType::icc, PCICC },
    { LinearSolverParameters::PreconditionerType::jacobi, PCJACOBI },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Petsc smoother option: " << type );
  return typeMap.at( type );
}

void CreatePetscSmoother( LinearSolverParameters const & params, PC precond )
{
  // Set up additive Schwartz outer preconditioner
  GEOSX_LAI_CHECK_ERROR( PCSetType( precond, PCASM ) );
  GEOSX_LAI_CHECK_ERROR( PCASMSetOverlap( precond, params.dd.overlap ) );
  GEOSX_LAI_CHECK_ERROR( PCASMSetType( precond, PC_ASM_RESTRICT ) );
  GEOSX_LAI_CHECK_ERROR( PCSetUp( precond ) );

  // Get local preconditioning context
  KSP * ksp_local;
  PetscInt n_local, first_local;
  GEOSX_LAI_CHECK_ERROR( PCASMGetSubKSP( precond, &n_local, &first_local, &ksp_local ) );
  GEOSX_LAI_ASSERT_EQ( n_local, 1 );

  // Set up local block ILU preconditioner
  PC prec_local;
  GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp_local[0], KSPPREONLY ) );
  GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp_local[0], &prec_local ) );
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, getPetscSmootherType( params.preconditionerType ) ) );
  GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, params.ilu.fill ) );
}

void PetscPreconditioner::compute( PetscMatrix const & mat )
{
  Base::compute( mat );

  bool const create = m_precond == nullptr;

  // Basic setup common for all preconditioners
  if( create )
  {
    GEOSX_LAI_CHECK_ERROR( PCCreate( mat.getComm(), &m_precond ) );
  }
  GEOSX_LAI_CHECK_ERROR( PCSetOperators( m_precond, mat.unwrapped(), mat.unwrapped() ) );

  // Add specifics
  if( create )
  {
    switch( m_parameters.preconditionerType )
    {
      case LinearSolverParameters::PreconditionerType::none:
      {
        GEOSX_LAI_CHECK_ERROR( PCSetType( m_precond, PCNONE ) );
        break;
      }
      case LinearSolverParameters::PreconditionerType::jacobi:
      {
        GEOSX_LAI_CHECK_ERROR( PCSetType( m_precond, PCJACOBI ) );
        break;
      }
      case LinearSolverParameters::PreconditionerType::amg:
      {
        CreatePetscAMG( m_parameters, m_precond );
        break;
      }
      case LinearSolverParameters::PreconditionerType::iluk:
      case LinearSolverParameters::PreconditionerType::icc:
      {
        CreatePetscSmoother( m_parameters, m_precond );
        break;
      }
      default:
      {
        GEOSX_ERROR( "Preconditioner type not supported in PETSc interface: " << m_parameters.preconditionerType );
      }
    }
  }

  // To be able to use Petsc preconditioner (e.g., GAMG) we need to disable floating point exceptions
  {
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

    GEOSX_LAI_CHECK_ERROR( PCSetUp( m_precond ) );
    GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( m_precond ) );
  }
}

void PetscPreconditioner::apply( PetscVector const & src,
                                 PetscVector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  GEOSX_LAI_CHECK_ERROR( PCApply( m_precond, src.unwrapped(), dst.unwrapped() ) );
}

void PetscPreconditioner::clear()
{
  PreconditionerBase::clear();
  if( m_precond != nullptr )
  {
    PCDestroy( &m_precond );
  }
}

PC const & PetscPreconditioner::unwrapped() const
{
  return m_precond;
}

}
