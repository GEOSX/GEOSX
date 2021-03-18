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
 * @file HyprePreconditioner.cpp
 */

#include "HyprePreconditioner.hpp"
#include "HypreMGR.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"
#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>

#include <_hypre_IJ_mv.h>

#include <cfenv>

namespace geosx
{

/**
 * @brief Container for hypre preconditioner null space data
 */
struct HypreNullSpace
{
  array1d< HYPRE_ParVector > vectors; ///< Hypre vectors containing the near null kernel
};

namespace
{

HYPRE_Int getHypreAMGCycleType( LinearSolverParameters::AMG::CycleType const & type )
{
  static std::map< LinearSolverParameters::AMG::CycleType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::CycleType::V, 1 },
    { LinearSolverParameters::AMG::CycleType::W, 2 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG cycle option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreAMGRelaxationType( LinearSolverParameters::AMG::SmootherType const & type )
{
  static std::map< LinearSolverParameters::AMG::SmootherType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::SmootherType::default_, -1 },
#ifdef GEOSX_USE_HYPRE_CUDA
    { LinearSolverParameters::AMG::SmootherType::jacobi, 7 },
#else
    { LinearSolverParameters::AMG::SmootherType::jacobi, 0 },
#endif
    { LinearSolverParameters::AMG::SmootherType::fgs, 3 },
    { LinearSolverParameters::AMG::SmootherType::bgs, 4 },
    { LinearSolverParameters::AMG::SmootherType::sgs, 6 },
    { LinearSolverParameters::AMG::SmootherType::l1sgs, 8 },
    { LinearSolverParameters::AMG::SmootherType::chebyshev, 16 },
    { LinearSolverParameters::AMG::SmootherType::l1jacobi, 18 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG relaxation option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreAMGCoarseType( LinearSolverParameters::AMG::CoarseType const & type )
{
  static std::map< LinearSolverParameters::AMG::CoarseType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::CoarseType::default_, -1 },
#ifdef GEOSX_USE_HYPRE_CUDA
    { LinearSolverParameters::AMG::CoarseType::jacobi, 7 },
#else
    { LinearSolverParameters::AMG::CoarseType::jacobi, 0 },
#endif
    { LinearSolverParameters::AMG::CoarseType::fgs, 3 },
    { LinearSolverParameters::AMG::CoarseType::bgs, 4 },
    { LinearSolverParameters::AMG::CoarseType::sgs, 6 },
    { LinearSolverParameters::AMG::CoarseType::l1sgs, 8 },
    { LinearSolverParameters::AMG::CoarseType::direct, 9 },
    { LinearSolverParameters::AMG::CoarseType::chebyshev, 16 },
    { LinearSolverParameters::AMG::CoarseType::l1jacobi, 18 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG relaxation option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreAMGCoarseningType( string const & type )
{
  static std::map< string, HYPRE_Int > const typeMap =
  {
    { "CLJP", 0 },
    { "Ruge-Stueben", 3 },
    { "Falgout", 6 },
    { "CLJP", 7 },
    { "PMIS", 8 },
    { "PMISD", 9 },
    { "HMIS", 10 },
    { "CGC", 21 },
    { "CGC-E", 22 }
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG coarsening option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreILUType( LinearSolverParameters::PreconditionerType const type )
{
  static std::map< LinearSolverParameters::PreconditionerType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, 0 },
    { LinearSolverParameters::PreconditionerType::ilut, 1 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre ILU option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreILUType( LinearSolverParameters::AMG::SmootherType const type )
{
  static std::map< LinearSolverParameters::AMG::SmootherType, HYPRE_Int > const typeMap =
  {
    { LinearSolverParameters::AMG::SmootherType::ilu0, 0 },
    { LinearSolverParameters::AMG::SmootherType::ilut, 1 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre ILU option: " << type );
  return typeMap.at( type );
}

void convertRigidBodyModes( arrayView1d< HypreVector > const & nearNullKernel,
                            array1d< HYPRE_ParVector > & nullSpacePointer )
{
  if( nearNullKernel.empty() )
  {
    return;
  }
  else
  {
    localIndex dim = 0;
    if( nearNullKernel.size() == 3 )
    {
      dim = 2;
    }
    else if( nearNullKernel.size() == 6 )
    {
      dim = 3;
    }
    else
    {
      GEOSX_ERROR( "Hypre preconditioner: rigid body modes can be either 3 or 6. Current number: " << nearNullKernel.size() );
    }
    localIndex const numRotations = LvArray::integerConversion< HYPRE_Int >( nearNullKernel.size() - dim );
    nullSpacePointer.resize( numRotations );
    for( localIndex k = 0; k < numRotations; ++k )
    {
      nullSpacePointer[k] = nearNullKernel[dim+k].unwrapped();
    }
  }
}

void createAMG( LinearSolverParameters const & params,
                HypreNullSpace const & nullSpace,
                HyprePrecWrapper & precond )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &precond.ptr ) );

  // Hypre's parameters to use BoomerAMG as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( precond.ptr, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( precond.ptr, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.logLevel ) ) );;
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( precond.ptr, params.dofsPerNode ) );

  // Set maximum number of multigrid levels (default 25)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxLevels( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.maxLevels ) ) );

  // Set type of cycle (1: V-cycle (default); 2: W-cycle)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( precond.ptr, getHypreAMGCycleType( params.amg.cycleType ) ) );

  if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes && !nullSpace.vectors.empty() )
  {
    // Set of options used in MFEM
    // Nodal coarsening options (nodal coarsening is required for this solver)
    // See hypre's new_ij driver and the paper for descriptions.

    // For further information, see:
    // Improving algebraic multigrid interpolation operators for linear elasticity problems
    // A. H. Baker Tz. V. Kolev U. M. Yang
    // Numerical Linear Algebra with Applications (2010) 17 (2-3), 495-517
    // doi:10.1002/nla.688

    HYPRE_Int const nodal                 = 4; // strength reduction norm: 1, 3 or 4
    HYPRE_Int const nodal_diag            = 1; // diagonal in strength matrix: 0, 1 or 2
    HYPRE_Int const relax_coarse          = 8; // smoother on the coarsest grid: 8, 99 or 29

    // Elasticity interpolation options
    HYPRE_Int const interp_vec_variant    = 2; // 1 = GM-1, 2 = GM-2, 3 = LN
    HYPRE_Int const q_max                 = 4; // max elements per row for each Q
    HYPRE_Int const smooth_interp_vectors = 1; // smooth the rigid-body modes?

    // Optionally pre-process the interpolation matrix through iterative weight
    // refinement (this is generally applicable for any system)
    HYPRE_Int const interp_refine         = 1;

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodal( precond.ptr, nodal ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodalDiag( precond.ptr, nodal_diag ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( precond.ptr, relax_coarse, 3 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecVariant( precond.ptr, interp_vec_variant ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecQMax( precond.ptr, q_max ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothInterpVectors( precond.ptr, smooth_interp_vectors ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpRefine( precond.ptr, interp_refine ) );

    // Add user-defined null space / rigid body mode support
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVectors( precond.ptr, nullSpace.vectors.size(), nullSpace.vectors.data() ) );
  }

  // Set smoother to be used (other options available, see hypre's documentation)
  // (default "gaussSeidel", i.e. local symmetric Gauss-Seidel)

  if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::ilu0 ||
      params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::ilut )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( precond.ptr, 5 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( precond.ptr, getHypreILUType( params.amg.smootherType ) ) );
  }
  else
  {
    HYPRE_Int const relaxType = getHypreAMGRelaxationType( params.amg.smootherType );
    if( relaxType >= 0 )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( precond.ptr, relaxType ) );
    }

    if( params.amg.smootherType == LinearSolverParameters::AMG::SmootherType::chebyshev )
    {
      // Set order for Chebyshev smoother valid options 1, 2 (default), 3, 4)
      if( ( 0 < params.amg.numSweeps ) && ( params.amg.numSweeps < 5 ) )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ) ) );
      }
    }
  }

  // Coarsening options: Only PMIS is supported on GPU
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCoarsenType( precond.ptr, getHypreAMGCoarseningType( params.amg.coarseningType ) ) );

  // Interpolation options: Use options 3, 6, 14 or 15.
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpType( precond.ptr, params.amg.interpolationType ) );

  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( precond.ptr, params.amg.numFunctions ) );

  if( params.amg.aggresiveNumLevels )
  {
    HYPRE_BoomerAMGSetAggNumLevels( precond.ptr, params.amg.aggresiveNumLevels ); // agg_num_levels = 1
  }

  HYPRE_BoomerAMGSetAggInterpType( precond.ptr, 5 ); // agg_interp_type = 5,7

  // Set coarsest level solver
  HYPRE_Int const coarseType = getHypreAMGCoarseType( params.amg.coarseType );
  if( coarseType >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( precond.ptr, coarseType, 3 ) );
  }

  // TODO Why does this cause a crash?
#if !defined(GEOSX_USE_HYPRE_CUDA)
  // (by default for coarsest grid size above 5,000 superlu_dist is used)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetDSLUThreshold( precond.ptr, 5000 ) );
#endif

  // Set the number of sweeps
  switch( params.amg.preOrPostSmoothing )
  {
    case LinearSolverParameters::AMG::PreOrPost::both:
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ) ) );
      break;
    }
    case LinearSolverParameters::AMG::PreOrPost::pre:
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ), 1 ) );
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, 0, 2 ) );
      break;
    }
    case LinearSolverParameters::AMG::PreOrPost::post:
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, 0, 1 ) );
      GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.amg.numSweeps ), 2 ) );
      break;
    }
    default:
    {
      GEOSX_ERROR( "Unsupported Hypre/AMG smoothing direction: " << params.amg.preOrPostSmoothing );
    }
  }

  // Set strength of connection
  if( params.amg.threshold > 0.0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( precond.ptr, params.amg.threshold ) );
  }

  precond.setup = HYPRE_BoomerAMGSetup;
  precond.solve = HYPRE_BoomerAMGSolve;
  precond.destroy = HYPRE_BoomerAMGDestroy;
}

void createILU( LinearSolverParameters const & params,
                HyprePrecWrapper & precond )
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &precond.ptr ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( precond.ptr, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( precond.ptr, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( precond.ptr, getHypreILUType( params.preconditionerType ) ) );

  if( params.ifact.fill >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( precond.ptr, LvArray::integerConversion< HYPRE_Int >( params.ifact.fill ) ) );
  }
  if( params.ifact.threshold >= 0 && params.preconditionerType == LinearSolverParameters::PreconditionerType::ilut )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetDropThreshold( precond.ptr, params.ifact.threshold ) );
  }

  precond.setup = HYPRE_ILUSetup;
  precond.solve = HYPRE_ILUSolve;
  precond.destroy = HYPRE_ILUDestroy;
}

} // namespace

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params )
  : Base{},
  m_params( std::move( params ) ),
  m_nullSpace( std::make_unique< HypreNullSpace >() )
{}

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params,
                                          arrayView1d< HypreVector > const & nearNullKernel )
  : HyprePreconditioner( std::move( params ) )
{
  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg &&
      m_params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
  {
    convertRigidBodyModes( nearNullKernel, m_nullSpace->vectors );
  }
}

HyprePreconditioner::~HyprePreconditioner()
{
  HyprePreconditioner::clear();
}

void HyprePreconditioner::create( DofManager const * const dofManager )
{
  switch( m_params.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_precond->setup = (HyprePrecWrapper::SetupFunc) hypre_ParKrylovIdentitySetup;
      m_precond->solve = (HyprePrecWrapper::SetupFunc) hypre_ParKrylovIdentity;
      break;
    }
    case LinearSolverParameters::PreconditionerType::jacobi:
    {
      m_precond->setup = HYPRE_ParCSRDiagScaleSetup;
      m_precond->solve = HYPRE_ParCSRDiagScale;
      break;
    }
    case LinearSolverParameters::PreconditionerType::fgs:
    {
      m_precond->solve = hypre::ForwardGaussSeidelSolve;
      break;
    }
    case LinearSolverParameters::PreconditionerType::bgs:
    {
      m_precond->solve = hypre::BackwardGaussSeidelSolve;
      break;
    }
    case LinearSolverParameters::PreconditionerType::sgs:
    {
      m_precond->solve = hypre::SymmetricGaussSeidelSolve;
      break;
    }
    case LinearSolverParameters::PreconditionerType::amg:
    {
      createAMG( m_params, *m_nullSpace, *m_precond );
      break;
    }
    case LinearSolverParameters::PreconditionerType::mgr:
    {
      m_mgrData = std::make_unique< HypreMGRData >();
      hypre::mgr::createMGR( m_params, dofManager, *m_precond, *m_mgrData );
      break;
    }
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
    {
      createILU( m_params, *m_precond );
      break;
    }
    case LinearSolverParameters::PreconditionerType::direct:
    {
      m_precond->setup = hypre::SuperLUDistSetup;
      m_precond->solve = hypre::SuperLUDistSolve;
      m_precond->destroy = hypre::SuperLUDistDestroy;
      break;
    }
    default:
    {
      GEOSX_ERROR( "Preconditioner type not supported in hypre interface: " << m_params.preconditionerType );
    }
  }
}

HypreMatrix const & HyprePreconditioner::setupPreconditioningMatrix( HypreMatrix const & mat )
{
  Stopwatch timer( m_componentFilterTime );
  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::mgr && m_params.mgr.separateComponents )
  {
    GEOSX_LAI_ASSERT_MSG( mat.dofManager() != nullptr, "MGR preconditioner requires a DofManager instance" );
    HypreMatrix Pu;
    HypreMatrix Auu;
    mat.dofManager()->makeRestrictor( { { m_params.mgr.displacementFieldName, { 3, true } } }, mat.getComm(), true, Pu );
    mat.multiplyPtAP( Pu, Auu );
    m_precondMatrix = LAIHelperFunctions::separateComponentFilter( Auu, m_params.dofsPerNode );
  }
  else if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg && m_params.amg.separateComponents )
  {
    m_precondMatrix = LAIHelperFunctions::separateComponentFilter( mat, m_params.dofsPerNode );
    return m_precondMatrix;
  }
  return mat;
}

void HyprePreconditioner::setup( Matrix const & mat )
{
  if( !m_precond )
  {
    m_precond = std::make_unique< HyprePrecWrapper >();
    create( mat.dofManager() );
  }

  HypreMatrix const & precondMat = setupPreconditioningMatrix( mat );
  create( mat.dofManager() );
  Base::setup( precondMat );

  // To be able to use Hypre preconditioner (e.g., BoomerAMG) we need to disable floating point exceptions
  {
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

    // Perform setup of the MGR mechanics F-solver with SDC matrix, if used
    if( m_mgrData && m_mgrData->mechSolver.ptr && m_mgrData->mechSolver.setup )
    {
      GEOSX_LAI_CHECK_ERROR( m_mgrData->mechSolver.setup( m_mgrData->mechSolver.ptr, m_precondMatrix.unwrapped(), nullptr, nullptr ) );
    }

    // Perform setup of the main solver, if needed
    if( m_precond->setup )
    {
      GEOSX_LAI_CHECK_ERROR( m_precond->setup( m_precond->ptr, precondMat.unwrapped(), nullptr, nullptr ) );
    }
  }
}

void HyprePreconditioner::apply( Vector const & src,
                                 Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.localSize(), numLocalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.localSize(), numLocalRows() );

  // Needed to avoid accumulation inside HYPRE solver phase
  dst.zero();

  GEOSX_LAI_CHECK_ERROR( m_precond->solve( m_precond->ptr, matrix().unwrapped(), src.unwrapped(), dst.unwrapped() ) );
}

void HyprePreconditioner::clear()
{
  Base::clear();
  if( m_precond && m_precond->ptr && m_precond->destroy )
  {
    GEOSX_LAI_CHECK_ERROR( m_precond->destroy( m_precond->ptr ) );
  }
  if( m_mgrData && m_mgrData->coarseSolver.ptr && m_mgrData->coarseSolver.destroy )
  {
    GEOSX_LAI_CHECK_ERROR( m_mgrData->coarseSolver.destroy( m_mgrData->coarseSolver.ptr ) );
  }
  if( m_mgrData && m_mgrData->mechSolver.ptr && m_mgrData->mechSolver.destroy )
  {
    GEOSX_LAI_CHECK_ERROR( m_mgrData->mechSolver.destroy( m_mgrData->mechSolver.ptr ) );
  }
  m_precond.reset();
  m_mgrData.reset();
}

HyprePrecWrapper const & HyprePreconditioner::unwrapped() const
{
  GEOSX_LAI_ASSERT( m_precond );
  return *m_precond;
}

}
