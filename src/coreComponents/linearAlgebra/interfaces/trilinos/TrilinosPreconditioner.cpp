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
 * @file TrilinosPreconditioner.cpp
 */

#include "TrilinosPreconditioner.hpp"

#include "linearAlgebra/utilities/LAIHelperFunctions.hpp"

#include <Epetra_Operator.h>
#include <Epetra_FEVector.h>
#include <Epetra_FECrsMatrix.h>
#include <ml_MultiLevelPreconditioner.h>
#include <Ifpack.h>

namespace geosx
{

TrilinosPreconditioner::TrilinosPreconditioner( LinearSolverParameters params )
  : Base{},
  m_params( std::move( params ) ),
  m_precond{},
  m_nullSpacePointer{}
{}

void convertRigidBodyModes( array1d< EpetraVector > const & nearNullKernel,
                            array2d< real64 > & nullSpacePointer )
{
  if( nearNullKernel.empty() )
  {
    return;
  }
  else
  {
    // ML expects the nullspace modes in a single double vector.
    localIndex const size = nearNullKernel[0].localSize();
    localIndex const numRBM = LvArray::integerConversion< integer >( nearNullKernel.size() );
    nullSpacePointer.resize( numRBM, size );
    for( localIndex k = 0; k < numRBM; ++k )
    {
      for( localIndex j = 0; j < size; ++j )
      {
        nullSpacePointer( k, j ) = nearNullKernel[k].get( nearNullKernel[k].getGlobalRowID( j ) );
      }
    }
  }
}

TrilinosPreconditioner::TrilinosPreconditioner( LinearSolverParameters params,
                                                array1d< Vector > const & nearNullKernel )
  : Base{},
  m_params( std::move( params ) ),
  m_precond{},
  m_nullSpacePointer{}
{
  if( m_params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes )
  {
    convertRigidBodyModes( nearNullKernel, m_nullSpacePointer );
  }
}

TrilinosPreconditioner::~TrilinosPreconditioner() = default;


namespace
{

string getMLCycleType( LinearSolverParameters::AMG::CycleType const & value )
{
  static std::map< LinearSolverParameters::AMG::CycleType, string > const optionMap =
  {
    { LinearSolverParameters::AMG::CycleType::V, "MGV" },
    { LinearSolverParameters::AMG::CycleType::W, "MGW" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported Trilinos/ML cycle option: " << value );
  return optionMap.at( value );
}

string getMLSmootherType( LinearSolverParameters::AMG::SmootherType const & value )
{
  static std::map< LinearSolverParameters::AMG::SmootherType, string > const optionMap =
  {
    { LinearSolverParameters::AMG::SmootherType::default_, "Chebyshev" },
    { LinearSolverParameters::AMG::SmootherType::jacobi, "Jacobi" },
    { LinearSolverParameters::AMG::SmootherType::l1jacobi, "Jacobi" },
    { LinearSolverParameters::AMG::SmootherType::fgs, "Gauss-Seidel" },
    { LinearSolverParameters::AMG::SmootherType::bgs, "Gauss-Seidel" },
    { LinearSolverParameters::AMG::SmootherType::sgs, "symmetric Gauss-Seidel" },
    { LinearSolverParameters::AMG::SmootherType::l1sgs, "symmetric Gauss-Seidel" },
    { LinearSolverParameters::AMG::SmootherType::chebyshev, "Chebyshev" },
    { LinearSolverParameters::AMG::SmootherType::ic0, "IC" },
    { LinearSolverParameters::AMG::SmootherType::ilu0, "ILU" },
    { LinearSolverParameters::AMG::SmootherType::ilut, "ILUT" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported Trilinos/ML smoother option: " << value );
  return optionMap.at( value );
}

string getMLCoarseType( LinearSolverParameters::AMG::CoarseType const & value )
{
  static std::map< LinearSolverParameters::AMG::CoarseType, string > const optionMap =
  {
    { LinearSolverParameters::AMG::CoarseType::default_, "Amesos-KLU" },
    { LinearSolverParameters::AMG::CoarseType::jacobi, "Jacobi" },
    { LinearSolverParameters::AMG::CoarseType::l1jacobi, "Jacobi" },
    { LinearSolverParameters::AMG::CoarseType::fgs, "Gauss-Seidel" },
    { LinearSolverParameters::AMG::CoarseType::bgs, "Gauss-Seidel" },
    { LinearSolverParameters::AMG::CoarseType::sgs, "symmetric Gauss-Seidel" },
    { LinearSolverParameters::AMG::CoarseType::l1sgs, "symmetric Gauss-Seidel" },
    { LinearSolverParameters::AMG::CoarseType::chebyshev, "Chebyshev" },
    { LinearSolverParameters::AMG::CoarseType::direct, "Amesos-KLU"},
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported Trilinos/ML coarse solver option: " << value );
  return optionMap.at( value );
}

string getMLPreOrPostSmoothingType( LinearSolverParameters::AMG::PreOrPost const & value )
{
  static std::map< LinearSolverParameters::AMG::PreOrPost, string > const optionMap =
  {
    { LinearSolverParameters::AMG::PreOrPost::pre, "pre" },
    { LinearSolverParameters::AMG::PreOrPost::post, "post" },
    { LinearSolverParameters::AMG::PreOrPost::both, "both" }
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0, "Unsupported Trilinos/ML smoothing direction option: " << value );
  return optionMap.at( value );
}

std::unique_ptr< Epetra_Operator >
createMLOperator( LinearSolverParameters const & params,
                  Epetra_RowMatrix const & matrix,
                  array2d< real64 > const & nullSpacePointer )
{
  Teuchos::ParameterList list;
  GEOSX_LAI_CHECK_ERROR( ML_Epetra::SetDefaults( params.isSymmetric ? "SA" : "NSSA", list ) );

  list.set( "ML output", params.logLevel );
  list.set( "max levels", params.amg.maxLevels );
  list.set( "aggregation: type", "Uncoupled" );
  list.set( "aggregation: threshold", params.amg.threshold );
  list.set( "PDE equations", params.dofsPerNode );
  list.set( "smoother: sweeps", params.amg.numSweeps );
  list.set( "prec type", getMLCycleType( params.amg.cycleType ) );
  list.set( "smoother: type", getMLSmootherType( params.amg.smootherType ) );
  list.set( "smoother: pre or post", getMLPreOrPostSmoothingType( params.amg.preOrPostSmoothing ) );
  list.set( "coarse: type", getMLCoarseType( params.amg.coarseType ) );

  if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::constantModes )
  {
    list.set( "null space: type", "default vectors" );
  }
  else if( params.amg.nullSpaceType == LinearSolverParameters::AMG::NullSpaceType::rigidBodyModes && !nullSpacePointer.empty() )
  {
    // Add user-defined null space / rigid body mode support
    list.set( "null space: type", "pre-computed" );
    list.set( "null space: vectors", nullSpacePointer.data() );
    list.set( "null space: dimension", LvArray::integerConversion< integer >( nullSpacePointer.size( 0 ) ) );
  }

  std::unique_ptr< Epetra_Operator > precond =
    std::make_unique< ML_Epetra::MultiLevelPreconditioner >( matrix, list );

  return precond;
}

Ifpack::EPrecType getIfpackPrecondType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, Ifpack::EPrecType > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::iluk, Ifpack::ILU },
    { LinearSolverParameters::PreconditionerType::ilut, Ifpack::ILUT },
    { LinearSolverParameters::PreconditionerType::ic, Ifpack::IC },
    { LinearSolverParameters::PreconditionerType::ict, Ifpack::ICT },
    { LinearSolverParameters::PreconditionerType::chebyshev, Ifpack::CHEBYSHEV },
    { LinearSolverParameters::PreconditionerType::jacobi, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::fgs, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::bgs, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::sgs, Ifpack::POINT_RELAXATION },
    { LinearSolverParameters::PreconditionerType::direct, Ifpack::AMESOS }
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return typeMap.at( type );
}
string getIfpackRelaxationType( LinearSolverParameters::PreconditionerType const & type )
{
  static std::map< LinearSolverParameters::PreconditionerType, string > const typeMap =
  {
    { LinearSolverParameters::PreconditionerType::jacobi, "Jacobi" },
    { LinearSolverParameters::PreconditionerType::fgs, "Gauss-Seidel" },
    { LinearSolverParameters::PreconditionerType::sgs, "symmetric Gauss-Seidel" },
  };

  GEOSX_LAI_ASSERT_MSG( typeMap.count( type ) > 0, "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return typeMap.at( type );
}

std::unique_ptr< Epetra_Operator >
createIfpackOperator( LinearSolverParameters const & params,
                      Epetra_RowMatrix const & matrix )
{
  Ifpack::EPrecType const type = getIfpackPrecondType( params.preconditionerType );
  std::unique_ptr< Ifpack_Preconditioner > precond( Ifpack::Create( type,
                                                                    const_cast< Epetra_RowMatrix * >( &matrix ),
                                                                    params.dd.overlap ) );
  Teuchos::ParameterList list;
  switch( type )
  {
    case Ifpack::POINT_RELAXATION:
    {
      list.set( "relaxation: type", getIfpackRelaxationType( params.preconditionerType ) );
      break;
    }
    case Ifpack::ILU:
    case Ifpack::IC:
    {
      list.set( "fact: level-of-fill", params.ifact.fill );
      break;
    }
    case Ifpack::ILUT:
    {
      list.set( "fact: ilut level-of-fill", std::max( real64( params.ifact.fill ), 1.0 ) );
      break;
    }
    case Ifpack::ICT:
    {
      list.set( "fact: ict level-of-fill", std::max( real64( params.ifact.fill ), 1.0 ) );
      break;
    }
    default:
    {
      break;
    }
  }

  GEOSX_LAI_CHECK_ERROR( precond->SetParameters( list ) );
  GEOSX_LAI_CHECK_ERROR( precond->Compute() );

  return precond;
}

} // namespace

EpetraMatrix const & TrilinosPreconditioner::setupPreconditioningMatrix( EpetraMatrix const & mat )
{
  if( m_params.preconditionerType == LinearSolverParameters::PreconditionerType::amg && m_params.amg.separateComponents )
  {
    m_precondMatrix = LAIHelperFunctions::separateComponentFilter( mat, m_params.dofsPerNode );
    return m_precondMatrix;
  }
  return mat;
}

void TrilinosPreconditioner::setup( Matrix const & mat )
{
  EpetraMatrix const & precondMat = setupPreconditioningMatrix( mat );
  Base::setup( precondMat );

  // To be able to use Trilinos solvers we need to disable floating point exceptions
  // In particular, we want to avoid a ML FPE crashes on Mac
  LvArray::system::FloatingPointExceptionGuard guard;

  switch( m_params.preconditionerType )
  {
    case LinearSolverParameters::PreconditionerType::amg:
    {
      m_precond = createMLOperator( m_params, precondMat.unwrapped(), m_nullSpacePointer );
      break;
    }
    case LinearSolverParameters::PreconditionerType::jacobi:
    case LinearSolverParameters::PreconditionerType::fgs:
    case LinearSolverParameters::PreconditionerType::bgs:
    case LinearSolverParameters::PreconditionerType::sgs:
    case LinearSolverParameters::PreconditionerType::chebyshev:
    case LinearSolverParameters::PreconditionerType::iluk:
    case LinearSolverParameters::PreconditionerType::ilut:
    case LinearSolverParameters::PreconditionerType::ic:
    case LinearSolverParameters::PreconditionerType::ict:
    case LinearSolverParameters::PreconditionerType::direct:
    {
      m_precond = createIfpackOperator( m_params, precondMat.unwrapped() );
      break;
    }
    case LinearSolverParameters::PreconditionerType::none:
    {
      m_precond.reset();
      break;
    }
    default:
    {
      GEOSX_ERROR( "Preconditioner type not supported in Trilinos interface: " << m_params.preconditionerType );
    }
  }
}

void TrilinosPreconditioner::apply( Vector const & src,
                                    Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  if( !m_precond )
  {
    dst.copy( src );
  }
  else
  {
    GEOSX_LAI_ASSERT( m_precond );
    m_precond->ApplyInverse( src.unwrapped(), dst.unwrapped() );
  }
}

void TrilinosPreconditioner::clear()
{
  PreconditionerBase::clear();
  m_precond.reset();
}

Epetra_Operator const & TrilinosPreconditioner::unwrapped() const
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

Epetra_Operator & TrilinosPreconditioner::unwrapped()
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

}
