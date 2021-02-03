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
#include "HypreMGRStrategies.hpp"

#include "linearAlgebra/DofManager.hpp"
#include "linearAlgebra/interfaces/hypre/HypreUtils.hpp"
#include "linearAlgebra/utilities/LinearSolverParameters.hpp"

#include <_hypre_utilities.h>
#include <_hypre_parcsr_ls.h>
#include <_hypre_IJ_mv.h>
#include <krylov.h>

#include <fenv.h>

namespace geosx
{

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params,
                                          DofManager const * const dofManager )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{},
  m_functions( std::make_unique< HyprePrecFuncs >() ),
  m_auxData{},
  m_dofManager( dofManager ),
  m_ready{},
  m_nearNullKernel{}
{}

HyprePreconditioner::HyprePreconditioner( LinearSolverParameters params,
                                          array1d< HypreVector > const & nearNullKernel,
                                          DofManager const * const dofManager )
  : Base{},
  m_parameters( std::move( params ) ),
  m_precond{},
  m_functions( std::make_unique< HyprePrecFuncs >() ),
  m_auxData{},
  m_dofManager( dofManager ),
  m_ready{},
  m_nearNullKernel( &nearNullKernel )
{}

HyprePreconditioner::~HyprePreconditioner()
{
  clear();
}

namespace
{

HYPRE_Int getHypreAMGCycleType( string const & type )
{
  static std::map< string, HYPRE_Int > const typeMap =
  {
    { "V", 1 },
    { "W", 2 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG cycle option: " << type );
  return typeMap.at( type );
}

HYPRE_Int getHypreAMGRelaxationType( string const & type )
{
  static std::map< string, HYPRE_Int > const typeMap =
  {
    { "jacobi", 0 },
    { "hybridForwardGaussSeidel", 3 },
    { "hybridBackwardGaussSeidel", 4 },
    { "hybridSymmetricGaussSeidel", 6 },
    { "gaussSeidel", 6 },
    { "L1hybridSymmetricGaussSeidel", 8 },
    { "chebyshev", 16 },
    { "L1jacobi", 18 },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Hypre AMG relaxation option: " << type );
  return typeMap.at( type );
}

void ConvertRigidBodyModes( array1d< HypreVector > const & nearNullKernel,
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
    localIndex const numRotations = toHYPRE_Int( nearNullKernel.size() - dim );
    nullSpacePointer.resize( numRotations );
    void * object;
    for( localIndex k = 0; k < numRotations; ++k )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_IJVectorGetObject( nearNullKernel[dim+k].unwrappedIJ(), &object ) );
      nullSpacePointer[k] = (HYPRE_ParVector) object;
    }
  }
}

}

void HyprePreconditioner::createHyprePreconditioner( DofManager const * const dofManager )
{
  switch( m_parameters.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_functions->setup = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentitySetup;
      m_functions->apply = (HYPRE_PtrToParSolverFcn) hypre_ParKrylovIdentity;
      break;
    }
    case LinearSolverParameters::PreconditionerType::jacobi:
    {
      m_functions->setup = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScaleSetup;
      m_functions->apply = (HYPRE_PtrToParSolverFcn) HYPRE_ParCSRDiagScale;
      break;
    }
    case LinearSolverParameters::PreconditionerType::amg:
    {
      createAMG();
      break;
    }
    case LinearSolverParameters::PreconditionerType::mgr:
    {
      createMGR( dofManager );
      break;
    }
    case LinearSolverParameters::PreconditionerType::iluk:
    {
      createILU();
      break;
    }
    case LinearSolverParameters::PreconditionerType::ilut:
    {
      createILUT();
      break;
    }
    default:
    {
      GEOSX_ERROR( "Preconditioner type not supported in hypre interface: " << m_parameters.preconditionerType );
    }
  }
}

void HyprePreconditioner::createAMG()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &m_precond ) );

  m_auxData = std::unique_ptr< HyprePrecAuxData >( new HyprePrecAuxData() );
  if( m_nearNullKernel != nullptr && m_auxData->nullSpacePointer.empty() && m_parameters.amg.nullSpaceType == "rigidBodyModes" )
  {
    ConvertRigidBodyModes( *m_nearNullKernel, m_auxData->nullSpacePointer );
  }

  // Hypre's parameters to use BoomerAMG as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( m_precond, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( m_precond, toHYPRE_Int( m_parameters.logLevel ) ) );;
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumFunctions( m_precond, m_parameters.dofsPerNode ) );

  // Set maximum number of multigrid levels (default 25)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxLevels( m_precond, toHYPRE_Int( m_parameters.amg.maxLevels ) ) );

  // Set type of cycle (1: V-cycle (default); 2: W-cycle)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleType( m_precond, getHypreAMGCycleType( m_parameters.amg.cycleType ) ) );

  if( m_parameters.amg.nullSpaceType == "rigidBodyModes" && m_auxData->nullSpacePointer.size() > 0 )
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

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodal( m_precond, nodal ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNodalDiag( m_precond, nodal_diag ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleRelaxType( m_precond, relax_coarse, 3 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecVariant( m_precond, interp_vec_variant ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVecQMax( m_precond, q_max ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothInterpVectors( m_precond, smooth_interp_vectors ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpRefine( m_precond, interp_refine ) );

    // Add user-defined null space / rigid body mode support
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetInterpVectors( m_precond, m_auxData->nullSpacePointer.size(), m_auxData->nullSpacePointer.data() ) );
  }

  // Set smoother to be used (other options available, see hypre's documentation)
  // (default "gaussSeidel", i.e. local symmetric Gauss-Seidel)
  if( m_parameters.amg.smootherType.substr( 0, 3 ) == "ilu" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetSmoothType( m_precond, 5 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( m_precond, 0 ) );
    if( m_parameters.amg.smootherType == "ilu1" )
    {
      GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, 1 ) );
    }
  }
  else
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxType( m_precond, getHypreAMGRelaxationType( m_parameters.amg.smootherType ) ) );

    if( m_parameters.amg.smootherType == "chebyshev" )
    {
      // Set order for Chebyshev smoother valid options 1, 2 (default), 3, 4)
      if( ( 0 < m_parameters.amg.numSweeps ) && ( m_parameters.amg.numSweeps < 5 ) )
      {
        GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetChebyOrder( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ) ) );
      }
    }
  }

  // Set coarsest level solver
  // (by default for coarsest grid size above 5,000 superlu_dist is used)
  GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetDSLUThreshold( m_precond, 5000 ) );
  if( m_parameters.amg.coarseType == "direct" )
  {
    GEOSX_LAI_CHECK_ERROR( hypre_BoomerAMGSetCycleRelaxType( m_precond, 9, 3 ) );
  }

  // Set the number of sweeps
  if( m_parameters.amg.preOrPostSmoothing == "both" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ) ) );
  }
  else if( m_parameters.amg.preOrPostSmoothing == "pre" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ), 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, 0, 2 ) );
  }
  else if( m_parameters.amg.preOrPostSmoothing == "post" )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, 0, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetCycleNumSweeps( m_precond, toHYPRE_Int( m_parameters.amg.numSweeps ), 2 ) );
  }

  // Set strength of connection
  if( m_parameters.amg.threshold > 0.0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetStrongThreshold( m_precond, m_parameters.amg.threshold ) );
  }

  m_functions->setup = HYPRE_BoomerAMGSetup;
  m_functions->apply = HYPRE_BoomerAMGSolve;
  m_functions->destroy = HYPRE_BoomerAMGDestroy;
}

void HyprePreconditioner::createILU()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &m_precond ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( m_precond, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( m_precond, 0 ) );

  if( m_parameters.ilu.fill >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, toHYPRE_Int( m_parameters.ilu.fill ) ) );
  }

  m_functions->setup = HYPRE_ILUSetup;
  m_functions->apply = HYPRE_ILUSolve;
  m_functions->destroy = HYPRE_ILUDestroy;
}

void HyprePreconditioner::createMGR( DofManager const * const dofManager )
{
  GEOSX_ERROR_IF( dofManager == nullptr, "MGR preconditioner requires a DofManager instance" );

  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRCreate( &m_precond ) );

  // Hypre's parameters to use MGR as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTol( m_precond, 0.0 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPrintLevel( m_precond, toHYPRE_Int( m_parameters.logLevel ) ) );

  array1d< localIndex > numComponentsPerField = dofManager->numComponentsPerField();
  array1d< localIndex > numLocalDofsPerField = dofManager->numLocalDofsPerField();

  if( m_auxData == nullptr )
  {
    m_auxData = std::unique_ptr< HyprePrecAuxData >( new HyprePrecAuxData() );
  }
  m_auxData->point_marker_array = computeLocalDofComponentLabels( numComponentsPerField,
                                                                  numLocalDofsPerField );

  if( m_parameters.logLevel >= 1 )
  {
    GEOSX_LOG_RANK_0( numComponentsPerField );
  }
  if( m_parameters.logLevel >= 2 )
  {
    GEOSX_LOG_RANK_VAR( numLocalDofsPerField );
  }
  if( m_parameters.logLevel >= 3 )
  {
    GEOSX_LOG_RANK_VAR( computeLocalDofComponentLabels( numComponentsPerField,
                                                        numLocalDofsPerField ) );
  }

  HYPRE_Int mgr_bsize;
  HYPRE_Int mgr_nlevels;
  std::vector< HYPRE_Int > mgr_num_cindexes;
  std::vector< std::vector< HYPRE_Int > > lv_cindexes;
  std::vector< HYPRE_Int * > mgr_cindexes;

  std::vector< HYPRE_Int > mgr_coarse_grid_method;
  std::vector< HYPRE_Int > mgr_level_interp_type;
  std::vector< HYPRE_Int > mgr_level_frelax_method;

  if( ( m_parameters.mgr.strategy == "Poroelastic" ) | ( m_parameters.mgr.strategy == "Hydrofracture" ) )
  {
    // Note: at the moment we assume single-phase flow poroelasticity
    //
    // dofLabel: 0 = displacement, x-component
    // dofLabel: 1 = displacement, y-component
    // dofLabel: 2 = displacement, z-component
    // dofLabel: 3 = pressure
    //
    // Ingredients
    //
    // 1. F-points displacement (0,1,2), C-points pressure (3)
    // 2. F-points smoother: AMG, single V-cycle, separate displacemente components
    // 3. C-points coarse-grid/Schur complement solver: boomer AMG
    // 4. Global smoother: none

    mgr_nlevels = 1;
    mgr_bsize = 4;

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = 1; // Eliminate displacement components

    lv_cindexes.resize( mgr_nlevels );
    lv_cindexes[0].push_back( 3 );

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );


    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[0] = 2; //diagonal scaling (Jacobi)

    mgr_coarse_grid_method.resize( mgr_nlevels );
    if( m_parameters.mgr.strategy == "Poroelastic" )
    {
      mgr_coarse_grid_method[0] = 1; //diagonal sparsification
    }
    else
    {
      mgr_coarse_grid_method[0] = 0; //Galerkin coarse grid computation using RAP
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( aux_precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( aux_precond, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( aux_precond, 1 ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxMethod( m_precond, 2 ) ); // AMG V-cycle
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( m_precond, 0 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( m_precond, mgr_coarse_grid_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_BoomerAMGDestroy;
  }
  else if( m_parameters.mgr.strategy == "SinglePhaseHybridFVM" )
  {
    // Labels description stored in point_marker_array
    // dofLabel: 0 = cell-centered pressure
    // dofLabel: 1 = face-centered pressure
    //
    // Ingredients
    //
    // 1. F-points cell-centered pressure, C-points face-centered pressure
    // 2. F-points smoother: l1-Jacobi
    // 3. C-points coarse-grid/Schur complement solver: boomer AMG
    // 4. Global smoother: none

    mgr_nlevels = 1;
    mgr_bsize = 2;

    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[0] = 2; // Diagonal scaling (Jacobi)

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = 1; // Eliminate cell-centered pressure

    lv_cindexes.resize( mgr_nlevels );
    lv_cindexes[0].push_back( 1 );

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( aux_precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( aux_precond, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( aux_precond, 1 ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetFRelaxMethod( m_precond, 18 ) ); // l1-Jacobi
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_BoomerAMGDestroy;

  }
  else if( m_parameters.mgr.strategy == "SinglePhaseReservoirFVM" )
  {
    // Labels description stored in point_marker_array
    // dofLabel: 0 = reservoir pressure (numResLabels = 1)
    // dofLabel: 1 = well pressure
    // dofLabel: 2 = well rate (numWellLabels = 2)
    //
    // Ingredients
    //
    // 1. F-points cell-centered pressure, C-points well vars
    // 2. F-points smoother: boomer AMG
    // 3. C-points coarse-grid/Schur complement solver: direct solver
    // 4. Global smoother: none

    mgr_nlevels = 1;
    mgr_bsize = 3;

    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[1] = 2;

    mgr_coarse_grid_method.resize( mgr_nlevels );
    mgr_coarse_grid_method[0] = 0;

    mgr_level_frelax_method.resize( mgr_nlevels );
    mgr_level_frelax_method[0] = 2; // AMG V-cycle

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = 2; // Eliminate cell-centered pressure

    lv_cindexes.resize( mgr_nlevels );
    lv_cindexes[0].push_back( 1 );
    lv_cindexes[0].push_back( 2 );

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( m_precond, mgr_level_frelax_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( m_precond, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( m_precond, 15 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( m_precond, mgr_coarse_grid_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_MGRDirectSolverDestroy;
  }
  else if( m_parameters.mgr.strategy == "SinglePhaseReservoirHybridFVM" )
  {
    // Labels description stored in point_marker_array
    // dofLabel: 0 = cell-centered reservoir pressure
    // dofLabel: 1 = face-centered reservoir pressure
    // dofLabel: 2 = well pressure
    // dofLabel: 3 = well rate
    //
    // Ingredients
    //
    // 2-level MGR reduction strategy which seems to work well for 2 components
    // 1st level: eliminate the cell-centered reservoir pressure
    // 2nd level: eliminate the face-centered reservoir pressure (Lagrange multiplier)
    // The coarse grid is the well vars
    // The coarse grid solved with a direct solver

    mgr_nlevels = 2;
    mgr_bsize = 4;

    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[0] = 2;
    mgr_level_interp_type[1] = 2;

    mgr_coarse_grid_method.resize( mgr_nlevels );
    mgr_coarse_grid_method[0] = 1;
    mgr_coarse_grid_method[1] = 0;

    mgr_level_frelax_method.resize( mgr_nlevels );
    mgr_level_frelax_method[0] = 0;
    mgr_level_frelax_method[1] = 2; // AMG V-cycle

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = mgr_bsize - 1; // eliminate the cell-centered pressure in the reservoir block
    mgr_num_cindexes[1] = mgr_bsize - 2; // eliminate the face-centered pressure in the reservoir block

    lv_cindexes.resize( mgr_nlevels );
    for( int cid=0; cid < mgr_bsize; cid++ )
    {
      // All points except the reservoir pressure
      if( cid != 0 )
      {
        lv_cindexes[0].push_back( cid );
      }
    }
    for( auto & cid : lv_cindexes[0] )
    {
      // Keep the well vars
      if( cid != 1 )
      {
        lv_cindexes[1].push_back( cid );
      }
    }

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( aux_precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( aux_precond, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( aux_precond, 1 ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( m_precond, mgr_level_frelax_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( m_precond, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( m_precond, 15 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( m_precond, mgr_coarse_grid_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_MGRDirectSolverDestroy;
  }
  else if( m_parameters.mgr.strategy == "CompositionalMultiphaseFlow" )
  {
    // Labels description stored in point_marker_array
    //             0 = pressure
    //             1 = density
    //           ... = densities
    // numLabels - 1 = density
    //
    // 2-level MGR reduction strategy which seems to work well for 2 components
    // 1st level: eliminate the reservoir density associated with the volume constraint
    // 2nd level: eliminate the pressure
    // The coarse grid solved with ILU(0)
    //
    // TODO:
    // - Experiment with block Jacobi for F-relaxation/interpolation of the reservoir densities
    // - Explore ways to reduce onto the pressure variable and use AMG for coarse-grid solve
    HYPRE_Int numLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );

    mgr_bsize = numLabels;
    mgr_nlevels = 2;

    /* options for solvers at each level */
    HYPRE_Int mgr_gsmooth_type = 16; // ILU(0)
    HYPRE_Int mgr_num_gsmooth_sweeps = 1;

    mgr_level_frelax_method.resize( mgr_nlevels );
    mgr_level_frelax_method[0] = 0; // Jacobi
    mgr_level_frelax_method[1] = 2; // AMG V-cycle

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = mgr_bsize - 1; // eliminate the last density reservoir block
    mgr_num_cindexes[1] = mgr_bsize - 2; // eliminate pressure

    lv_cindexes.resize( mgr_nlevels );
    for( int cid=0; cid < mgr_bsize; cid++ )
    {
      // All points except the last density
      // which corresponds to the volume constraint equation
      if( cid < numLabels - 1 )
      {
        lv_cindexes[0].push_back( cid );
      }
    }
    for( auto & cid : lv_cindexes[0] )
    {
      // eliminate pressure
      if( cid != 0 )
      {
        lv_cindexes[1].push_back( cid );
      }
    }

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }


    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetType( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( aux_precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( aux_precond, 0.0 ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( m_precond, mgr_level_frelax_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( m_precond, mgr_gsmooth_type ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, mgr_num_gsmooth_sweeps ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_ILUSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_ILUSetup,
                                aux_precond )
      );
    m_functions->aux_destroy = HYPRE_ILUDestroy;
  }
  else if( m_parameters.mgr.strategy == "CompositionalMultiphaseReservoir" )
  {
    // Labels description stored in point_marker_array
    //                0 = reservoir pressure
    //                1 = reservoir density
    //              ... = ... (reservoir densities)
    // numResLabels - 1 = reservoir density
    //     numResLabels = well pressure
    // numResLabels + 1 = well density
    //              ... = ... (well densities)
    // numResLabels + numWellLabels - 2 = well density
    // numResLabels + numWellLabels - 1 = well rate
    //
    // 3-level MGR reduction strategy which seems to work well for 2 components
    // 1st level: eliminate the reservoir density associated with the volume constraint
    // 2nd level: eliminate the rest of the reservoir densities
    // 3rd level: eliminate the pressure
    // The coarse grid is the well block and solved with a direct solver
    //
    // TODO:
    // - Use block Jacobi for F-relaxation/interpolation of the reservoir densities (2nd level)

    HYPRE_Int numResLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[0] );
    HYPRE_Int numWellLabels = LvArray::integerConversion< HYPRE_Int >( numComponentsPerField[1] );

    mgr_bsize = numResLabels + numWellLabels;
    mgr_nlevels = 3;

    /* options for solvers at each level */
    HYPRE_Int mgr_gsmooth_type = 16;
    HYPRE_Int mgr_num_gsmooth_sweeps = 0;

    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[0] = 2;
    mgr_level_interp_type[1] = 2;
    mgr_level_interp_type[2] = 2;

    mgr_coarse_grid_method.resize( mgr_nlevels );
    mgr_coarse_grid_method[0] = 1;
    mgr_coarse_grid_method[1] = 1;
    mgr_coarse_grid_method[2] = 0;

    mgr_level_frelax_method.resize( mgr_nlevels );
    mgr_level_frelax_method[0] = 0; // Jacobi
    mgr_level_frelax_method[1] = 0; // Jacobi
    mgr_level_frelax_method[2] = 2; // AMG V-cycle

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = mgr_bsize - 1; // eliminate the last density in the reservoir block
    mgr_num_cindexes[1] = mgr_bsize - numResLabels + 1; // eliminate all densities reservoir block
    mgr_num_cindexes[2] = mgr_bsize - numResLabels; // eliminate reservoir block

    lv_cindexes.resize( mgr_nlevels );
    for( int cid=0; cid < mgr_bsize; cid++ )
    {
      // All points except the last reservoir density
      // which corresponds to the volume constraint equation
      if( cid != numResLabels - 1 )
      {
        lv_cindexes[0].push_back( cid );
      }
    }
    for( auto & cid : lv_cindexes[0] )
    {
      // eliminate the rest of the reservoir densities
      if( cid == 0 || cid >= numResLabels )
      {
        lv_cindexes[1].push_back( cid );
      }
    }
    for( auto & cid : lv_cindexes[1] )
    {
      // eliminate the reservoir pressure
      if( cid != 0 )
      {
        lv_cindexes[2].push_back( cid );
      }
    }

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRDirectSolverCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByPointMarkerArray( m_precond, mgr_bsize, mgr_nlevels,
                                                                  mgr_num_cindexes.data(),
                                                                  mgr_cindexes.data(),
                                                                  m_auxData->point_marker_array.data() ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( m_precond, mgr_level_frelax_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetTruncateCoarseGridThreshold( m_precond, 1e-14 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetPMaxElmts( m_precond, 15 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( m_precond, mgr_coarse_grid_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( m_precond, mgr_gsmooth_type ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, mgr_num_gsmooth_sweeps ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_MGRDirectSolverSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_MGRDirectSolverDestroy;
  }
  else if( m_parameters.mgr.strategy == "LagrangianContactMechanics" )
  {
    // Contact mechanis with face-centered lagrangian multipliers
    //
    // dofLabel: 0 = displacement, x-component
    // dofLabel: 1 = displacement, y-component
    // dofLabel: 2 = displacement, z-component
    // dofLabel: 3 = pressure
    //
    // Ingredients
    //
    // 1. F-points displacement (0,1,2), C-points pressure (3)
    // 2. F-points smoother: AMG, single V-cycle, separate displacemente components
    // 3. C-points coarse-grid/Schur complement solver: boomer AMG
    // 4. Global smoother: none

    mgr_nlevels = 1;
    mgr_bsize = 2;

    HYPRE_Int mgr_relax_type = 0;
    HYPRE_Int mgr_gsmooth_type = 0;
    HYPRE_Int mgr_num_gsmooth_sweeps = 1;
    HYPRE_Int mgr_num_relax_sweeps = 0; // skip F-relaxation

    std::vector< HYPRE_Int > mgr_level_restrict_type( mgr_bsize );
    mgr_level_restrict_type[0] = 0;
    HYPRE_Int mgr_num_restrict_sweeps = 0;
    HYPRE_Int mgr_num_interp_sweeps = 0;

    mgr_num_cindexes.resize( mgr_nlevels );
    mgr_num_cindexes[0] = 1; // Eliminate lagrangian components

    std::vector< HYPRE_BigInt > mgr_idx_array( mgr_bsize );
    mgr_idx_array[0] = toHYPRE_BigInt( dofManager->rankOffset() );
    mgr_idx_array[1] = toHYPRE_BigInt( dofManager->rankOffset() + dofManager->numLocalDofs( m_parameters.mgr.displacementFieldName ) );

    lv_cindexes.resize( mgr_nlevels );
    lv_cindexes[0].push_back( 0 );

    mgr_cindexes.resize( mgr_nlevels );
    for( HYPRE_Int iLevel = 0; iLevel < mgr_nlevels; ++iLevel )
    {
      mgr_cindexes[iLevel] = lv_cindexes[iLevel].data();
    }

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCpointsByContiguousBlock( m_precond, mgr_bsize, mgr_nlevels,
                                                                 mgr_idx_array.data(),
                                                                 mgr_num_cindexes.data(),
                                                                 mgr_cindexes.data() ) );


    mgr_level_interp_type.resize( mgr_nlevels );
    mgr_level_interp_type[0] = 2; // diagonal scaling (Jacobi)

    mgr_level_frelax_method.resize( mgr_nlevels );
    mgr_level_frelax_method[0] = 0; // Jacobi

    mgr_coarse_grid_method.resize( mgr_nlevels );
    mgr_coarse_grid_method[0] = 0; // all

    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGCreate( &aux_precond ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetPrintLevel( aux_precond, 0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetMaxIter( aux_precond, 1 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetTol( aux_precond, 0.0 ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_BoomerAMGSetRelaxOrder( aux_precond, 1 ) );

    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelRestrictType( m_precond, mgr_level_restrict_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRestrictSweeps( m_precond, mgr_num_restrict_sweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumInterpSweeps( m_precond, mgr_num_interp_sweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNumRelaxSweeps( m_precond, mgr_num_relax_sweeps ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetNonCpointsToFpoints( m_precond, 1 ));
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetRelaxType( m_precond, mgr_relax_type ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelFRelaxMethod( m_precond, mgr_level_frelax_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetLevelInterpType( m_precond, mgr_level_interp_type.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetCoarseGridMethod( m_precond, mgr_coarse_grid_method.data() ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetGlobalsmoothType( m_precond, mgr_gsmooth_type ) );
    GEOSX_LAI_CHECK_ERROR( HYPRE_MGRSetMaxGlobalsmoothIters( m_precond, mgr_num_gsmooth_sweeps ) );
    GEOSX_LAI_CHECK_ERROR(
      HYPRE_MGRSetCoarseSolver( m_precond,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSolve,
                                (HYPRE_PtrToParSolverFcn)HYPRE_BoomerAMGSetup,
                                aux_precond )
      );

    m_functions->aux_destroy = HYPRE_BoomerAMGDestroy;
  }
  else
  {
    GEOSX_ERROR( "Unsupported MGR strategy: " << m_parameters.mgr.strategy );
  }
  m_functions->setup = HYPRE_MGRSetup;
  m_functions->apply = HYPRE_MGRSolve;
  m_functions->destroy = HYPRE_MGRDestroy;
}

void HyprePreconditioner::createILUT()
{
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUCreate( &m_precond ) );

  // Hypre's parameters to use ParCSR ILU as a preconditioner
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetMaxIter( m_precond, 1 ) );
  GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetTol( m_precond, 0.0 ) );

  if( m_parameters.ilu.fill >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetLevelOfFill( m_precond, toHYPRE_Int( m_parameters.ilu.fill ) ) );
  }
  if( m_parameters.ilu.threshold >= 0 )
  {
    GEOSX_LAI_CHECK_ERROR( HYPRE_ILUSetDropThreshold( m_precond,
                                                      m_parameters.ilu.threshold ) );
  }

  m_functions->setup = HYPRE_ILUSetup;
  m_functions->apply = HYPRE_ILUSolve;
  m_functions->destroy = HYPRE_ILUDestroy;
}

void HyprePreconditioner::create()
{
  if( !m_ready )
  {
    // Basic setup for functions
    if( !m_functions )
    {
      m_functions = std::make_unique< HyprePrecFuncs >();
    }

    // Basic setup common for all preconditioners
    if( m_precond == nullptr )
    {
      createHyprePreconditioner( m_dofManager );
    }
    m_ready = true;
  }
}

void HyprePreconditioner::compute( Matrix const & mat )
{
  create();

  PreconditionerBase::compute( mat );

  // To be able to use Hypre preconditioner (e.g., BoomerAMG) we need to disable floating point exceptions
  {
    LvArray::system::FloatingPointExceptionGuard guard( FE_ALL_EXCEPT );

    GEOSX_LAI_CHECK_ERROR( m_functions->setup( m_precond, mat.unwrapped(), nullptr, nullptr ) );
  }
}

void HyprePreconditioner::apply( Vector const & src,
                                 Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  // Needed to avoid accumulation inside HYPRE solver phase
  dst.zero();

  GEOSX_LAI_CHECK_ERROR( m_functions->apply( m_precond, this->matrix().unwrapped(), src.unwrapped(), dst.unwrapped() ) );
}

void HyprePreconditioner::clear()
{
  PreconditionerBase::clear();
  if( m_precond != nullptr && m_functions && m_functions->destroy != nullptr )
  {
    GEOSX_LAI_CHECK_ERROR( m_functions->destroy( m_precond ) );
    m_precond = nullptr;
  }
  if( aux_precond != nullptr && m_functions && m_functions->aux_destroy != nullptr )
  {
    GEOSX_LAI_CHECK_ERROR( m_functions->aux_destroy( aux_precond ) );
    aux_precond = nullptr;
  }
  m_functions.reset();
  m_ready = false;
}

HYPRE_Solver const & HyprePreconditioner::unwrapped() const
{
  return m_precond;
}

HyprePrecFuncs const & HyprePreconditioner::unwrappedFuncs() const
{
  GEOSX_LAI_ASSERT( m_functions );
  return *m_functions;
}

}
