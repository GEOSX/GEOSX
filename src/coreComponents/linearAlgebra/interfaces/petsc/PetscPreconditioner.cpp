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

namespace geosx
{

namespace
{

void convertRigidBodyModes( LinearSolverParameters const & params,
                            array1d< PetscVector > const & nearNullKernel,
                            MatNullSpace & nullsp )
{
  if( nearNullKernel.empty() )
  {
    nullsp = nullptr;
    return;
  }
  else
  {
    localIndex const numRBM = LvArray::integerConversion< localIndex >( nearNullKernel.size() );
    array1d< Vec > nullvecs( numRBM );
    for( localIndex i = 0; i < numRBM; ++i )
    {
      GEOSX_LAI_CHECK_ERROR( VecDuplicate( nearNullKernel[i].unwrapped(), &nullvecs[i] ) );
      GEOSX_LAI_CHECK_ERROR( VecCopy( nearNullKernel[i].unwrapped(), nullvecs[i] ) );
      GEOSX_LAI_CHECK_ERROR( VecSetBlockSize( nullvecs[i], params.dofsPerNode ) );
      GEOSX_LAI_CHECK_ERROR( VecSetUp( nullvecs[i] ) );
    }
    GEOSX_LAI_CHECK_ERROR( MatNullSpaceCreate( MPI_COMM_GEOSX, PETSC_FALSE, numRBM, nullvecs.data(), &nullsp ) );
    for( localIndex i = 0; i < numRBM; ++i )
    {
      GEOSX_LAI_CHECK_ERROR( VecDestroy( &nullvecs[i] ) );
    }
  }
}

PCType getPetscSmootherType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, PCType > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, PCILU },
    { LinearSolverParameters::PreconditionerType::ic, PCICC },
    { LinearSolverParameters::PreconditionerType::jacobi, PCJACOBI },
    { LinearSolverParameters::PreconditionerType::l1jacobi, PCJACOBI },
    { LinearSolverParameters::PreconditionerType::fgs, PCSOR },
    { LinearSolverParameters::PreconditionerType::bgs, PCSOR },
    { LinearSolverParameters::PreconditionerType::sgs, PCSOR },
    { LinearSolverParameters::PreconditionerType::l1sgs, PCSOR },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Petsc smoother option: " << type );
  return typeMap.at( type );
}

PCJacobiType getPetscJacobiType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, PCJacobiType > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::jacobi, PC_JACOBI_DIAGONAL },
    { LinearSolverParameters::PreconditionerType::l1jacobi, PC_JACOBI_ROWSUM },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Petsc smoother option: " << type );
  return typeMap.at( type );
}

MatSORType getPetscSORType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, MatSORType > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::fgs, SOR_FORWARD_SWEEP },
    { LinearSolverParameters::PreconditionerType::bgs, SOR_BACKWARD_SWEEP },
    { LinearSolverParameters::PreconditionerType::sgs, SOR_SYMMETRIC_SWEEP },
    { LinearSolverParameters::PreconditionerType::l1sgs, SOR_SYMMETRIC_SWEEP },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Petsc smoother option: " << type );
  return typeMap.at( type );
}

void createPetscSmoother( LinearSolverParameters const & params, PC const precond )
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
  PCType const prec_type = getPetscSmootherType( params.preconditionerType );
  GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, prec_type ) );
  if( !strcmp( prec_type, PCJACOBI ) )
  {
    GEOSX_LAI_CHECK_ERROR( PCJacobiSetType( prec_local, getPetscJacobiType( params.preconditionerType ) ) );
  }
  if( !strcmp( prec_type, PCSOR ) )
  {
    GEOSX_LAI_CHECK_ERROR( PCSORSetSymmetric( prec_local, getPetscSORType( params.preconditionerType ) ) );
  }
  GEOSX_LAI_CHECK_ERROR( PCFactorSetLevels( prec_local, params.ifact.fill ) );
}

PCMGCycleType getPetscMGCycleType( LinearSolverParameters::AMG::CycleType const & type )
{
  static std::map< LinearSolverParameters::AMG::CycleType, PCMGCycleType > const typeMap =
  {
    { LinearSolverParameters::AMG::CycleType::V, PC_MG_CYCLE_V },
    { LinearSolverParameters::AMG::CycleType::W, PC_MG_CYCLE_W }
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Petsc MG cycle option: " << type );
  return typeMap.at( type );
}

void createPetscAMG( LinearSolverParameters const & params,
                     PC const precond,
                     MatNullSpace const nullsp )
{
  // Default options only for the moment
  GEOSX_LAI_CHECK_ERROR( PCSetType( precond, PCGAMG ) );

  if( !params.isSymmetric )
  {
    // Usually GEOSX matrix is not symmetric, but GAMG is designed for symmetric matrices
    // In case of a general matrix, we need to compute a symmetric graph (slightly heavier
    // than the default, but necessary in case of asymmetric matrices).
    GEOSX_LAI_CHECK_ERROR( PCGAMGSetSymGraph( precond, PETSC_TRUE ) );
  }

  // Add user-defined null space / rigid body mode support
  if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes && nullsp )
  {
    Mat mat;
    PCGetOperators( precond, nullptr, &mat );
    GEOSX_LAI_CHECK_ERROR( MatSetNearNullSpace( mat, nullsp ) );
  }

  // Set type of cycle
  // GEOSX_LAI_CHECK_ERROR( PCMGSetCycleType( precond, getPetscMGCycleType( params.amg.cycleType ) ) );
  getPetscMGCycleType( params.amg.cycleType );

  // Set max number of levels
  GEOSX_LAI_CHECK_ERROR( PCGAMGSetNlevels( precond, params.amg.maxLevels ) );

  // TODO: need someone familiar with PETSc to take a look at this
#if 0
  GEOSX_LAI_CHECK_ERROR( PCSetType( precond, PCHMG ) );
  GEOSX_LAI_CHECK_ERROR( PCHMGSetInnerPCType( precond, PCGAMG ) );

  // Set maximum number of multigrid levels
  if( params.amg.maxLevels > 0 )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetLevels( precond,
                                          LvArray::integerConversion< PetscInt >( params.amg.maxLevels ),
                                          nullptr ) );
  }

  // Set the number of sweeps
  if( params.amg.numSweeps > 1 )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGSetNumberSmooth( precond,
                                                LvArray::integerConversion< PetscInt >( params.amg.numSweeps ) ) );
  }

  // Set smoother to be used (for all levels)
  PetscInt numLevels;
  PetscInt l;
  KSP smoother;
  PC smootherPC;
  GEOSX_LAI_CHECK_ERROR( PCMGGetLevels( precond, &numLevels ) );
  for( l = 0; l < numLevels; ++l )
  {
    GEOSX_LAI_CHECK_ERROR( PCMGGetSmoother( precond, l, &smoother ) );
    GEOSX_LAI_CHECK_ERROR( KSPSetType( smoother, KSPRICHARDSON ) );
    GEOSX_LAI_CHECK_ERROR( KSPGetPC( smoother, &smootherPC ) );

    if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::jacobi )
    {
      GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCJACOBI ) );
    }
    else if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::fgs )
    {
      GEOSX_LAI_CHECK_ERROR( PCSetType( smootherPC, PCSOR ) );
    }
    else if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::ilu0 )
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
      GEOSX_LAI_ASSERT_EQ( first_local, MpiWrapper::commRank( MPI_COMM_GEOSX ) );

      // Set up local block ILU preconditioner
      PC prec_local;
      GEOSX_LAI_CHECK_ERROR( KSPSetType( ksp_local[0], KSPPREONLY ) );
      GEOSX_LAI_CHECK_ERROR( KSPGetPC( ksp_local[0], &prec_local ) );
      GEOSX_LAI_CHECK_ERROR( PCSetType( prec_local, PCILU ) );
      // GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( smootherPC ) );
    }
  }
#endif
}

void createPetscDirect( LinearSolverParameters const & GEOSX_UNUSED_PARAM( params ), PC const precond )
{
  GEOSX_LAI_CHECK_ERROR( PCSetType( precond, PCLU ) );
  GEOSX_LAI_CHECK_ERROR( PCFactorSetMatSolverType( precond, MATSOLVERSUPERLU_DIST ) );
}

} // namespace

PetscPreconditioner::PetscPreconditioner( LinearSolverParameters params )
  : Base{},
  m_params( std::move( params ) ),
  m_precond{},
  m_nullsp{}
{ }

PetscPreconditioner::PetscPreconditioner( LinearSolverParameters params,
                                          array1d< Vector > const & nearNullKernel )
  : Base{},
  m_params( std::move( params ) ),
  m_precond{},
  m_nullsp{}
{
  if( m_params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
  {
    convertRigidBodyModes( params, nearNullKernel, m_nullsp );
  }
}

PetscPreconditioner::~PetscPreconditioner()
{
  PetscPreconditioner::clear();
  if( m_nullsp != nullptr )
  {
    MatNullSpaceDestroy( &m_nullsp );
  }
}

PetscMatrix const & PetscPreconditioner::setupPreconditioningMatrix( PetscMatrix const & mat )
{
  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg && m_params.amg.separateComponents )
  {
    m_precondMatrix = LAIHelperFunctions::separateComponentFilter( mat, m_params.dofsPerNode );
    return m_precondMatrix;
  }
  return mat;
}

void PetscPreconditioner::setup( PetscMatrix const & mat )
{
  PetscMatrix const & precondMat = setupPreconditioningMatrix( mat );
  Base::setup( precondMat );

  // Set dofs per node (it can be done only at the matrix level ...)
  GEOSX_LAI_CHECK_ERROR( MatSetBlockSize( precondMat.unwrapped(), m_params.dofsPerNode ) );

  bool const create = m_precond == nullptr;

  // Basic setup common for all preconditioners
  if( create )
  {
    GEOSX_LAI_CHECK_ERROR( PCCreate( precondMat.getComm(), &m_precond ) );
  }
  GEOSX_LAI_CHECK_ERROR( PCSetOperators( m_precond, mat.unwrapped(), precondMat.unwrapped() ) );

  // To be able to use PETSc solvers we need to disable floating point exceptions
  LvArray::system::FloatingPointExceptionGuard guard;

  // Add specifics
  if( create )
  {
    switch( m_params.preconditionerType )
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
        createPetscAMG( m_params, m_precond, m_nullsp );
        break;
      }
      case LinearSolverParameters::PreconditionerType::iluk:
      case LinearSolverParameters::PreconditionerType::ic:
      {
        createPetscSmoother( m_params, m_precond );
        break;
      }
      case LinearSolverParameters::PreconditionerType::direct:
      {
        createPetscDirect( m_params, m_precond );
        break;
      }
      default:
      {
        GEOSX_ERROR( "Preconditioner type not supported in PETSc interface: " << m_params.preconditionerType );
      }
    }
  }

  GEOSX_LAI_CHECK_ERROR( PCSetUp( m_precond ) );
  GEOSX_LAI_CHECK_ERROR( PCSetUpOnBlocks( m_precond ) );
}

void PetscPreconditioner::apply( Vector const & src,
                                 Vector & dst ) const
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
