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
 * @file PreconditionerAMG.cpp
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
TrilinosPreconditioner::TrilinosPreconditioner( LinearSolverParameters params ) :
  Base {},
  m_parameters( std::move( params ) ),
  m_precond {}
{}

TrilinosPreconditioner::~TrilinosPreconditioner() = default;

namespace
{
string const &
getMLCycleType( string const & value )
{
  static std::map< string, string > const optionMap = {
    { "V", "MGV" },
    { "W", "MGW" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0,
                        "Unsupported Trilinos/ML cycle option: " << value );
  return optionMap.at( value );
}

string const &
getMLSmootherType( string const & value )
{
  static std::map< string, string > const optionMap = {
    { "jacobi", "Jacobi" },
    { "blockJacobi", "block Jacobi" },
    { "gaussSeidel", "Gauss-Seidel" },
    { "blockGaussSeidel", "block Gauss-Seidel" },
    { "chebyshev", "Chebyshev" },
    { "icc", "IC" },
    { "ilu", "ILU" },
    { "ilut", "ILUT" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0,
                        "Unsupported Trilinos/ML smoother option: " << value );
  return optionMap.at( value );
}

string const &
getMLCoarseType( string const & value )
{
  static std::map< string, string > const optionMap = {
    { "jacobi", "Jacobi" },
    { "gaussSeidel", "Gauss-Seidel" },
    { "blockGaussSeidel", "block Gauss-Seidel" },
    { "chebyshev", "Chebyshev" },
    { "direct", "Amesos-KLU" },
  };

  GEOSX_LAI_ASSERT_MSG( optionMap.count( value ) > 0,
                        "Unsupported Trilinos/ML coarse solver option: " << value );
  return optionMap.at( value );
}

std::unique_ptr< Epetra_Operator >
CreateMLOperator( LinearSolverParameters const & params,
                  Epetra_RowMatrix const & matrix )
{
  GEOSX_LAI_ASSERT_EQ( params.preconditionerType, "amg" );

  Teuchos::ParameterList list;
  GEOSX_LAI_CHECK_ERROR(
    ML_Epetra::SetDefaults( params.isSymmetric ? "SA" : "NSSA", list ) );

  list.set( "ML output", params.logLevel );
  list.set( "max levels", params.amg.maxLevels );
  list.set( "aggregation: type", "Uncoupled" );
  list.set( "aggregation: threshold", params.amg.threshold );
  list.set( "PDE equations", params.dofsPerNode );
  list.set( "smoother: sweeps", params.amg.numSweeps );
  list.set( "prec type", getMLCycleType( params.amg.cycleType ) );
  list.set( "smoother: type", getMLSmootherType( params.amg.smootherType ) );
  list.set( "smoother: pre or post", params.amg.preOrPostSmoothing );
  list.set( "coarse: type", getMLCoarseType( params.amg.coarseType ) );

  std::unique_ptr< Epetra_Operator > precond =
    std::make_unique< ML_Epetra::MultiLevelPreconditioner >( matrix, list );

  return precond;
}

Ifpack::EPrecType
getIfpackPrecondType( string const & type )
{
  static std::map< string, Ifpack::EPrecType > const typeMap = {
    { "iluk", Ifpack::ILU },
    { "ilut", Ifpack::ILUT },
    { "icc", Ifpack::IC },
    { "ict", Ifpack::ICT },
    { "jacobi", Ifpack::POINT_RELAXATION },
    { "gs", Ifpack::POINT_RELAXATION },
    { "sgs", Ifpack::POINT_RELAXATION },
  };

  GEOSX_LAI_ASSERT_MSG(
    typeMap.count( type ) > 0,
    "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return typeMap.at( type );
}
string
getIfpackRelaxationType( string const & type )
{
  static std::map< string, string > const typeMap = {
    { "jacobi", "Jacobi" },
    { "gs", "Gauss-Seidel" },
    { "sgs", "symmetric Gauss-Seidel" },
  };

  GEOSX_LAI_ASSERT_MSG(
    typeMap.count( type ) > 0,
    "Unsupported Trilinos/Ifpack preconditioner option: " << type );
  return typeMap.at( type );
}

std::unique_ptr< Epetra_Operator >
CreateIfpackOperator(
  LinearSolverParameters const & params,
  Epetra_RowMatrix const & matrix )
{
  Ifpack::EPrecType const type = getIfpackPrecondType( params.preconditionerType );
  std::unique_ptr< Ifpack_Preconditioner > precond(
    Ifpack::Create( type,
                    const_cast< Epetra_RowMatrix * >( &matrix ),
                    params.dd.overlap ) );
  Teuchos::ParameterList list;
  if( type == Ifpack::POINT_RELAXATION )
  {
    list.set( "relaxation: type",
              getIfpackRelaxationType( params.preconditionerType ) );
  }
  else
  {
    list.set( "fact: level-of-fill", params.ilu.fill );
    list.set( "fact: ilut level-of-fill", std::max( real64( params.ilu.fill ), 1.0 ) );
    list.set( "fact: ict level-of-fill", std::max( real64( params.ilu.fill ), 1.0 ) );
  }

  GEOSX_LAI_CHECK_ERROR( precond->SetParameters( list ) );
  GEOSX_LAI_CHECK_ERROR( precond->Compute() );

  return precond;
}

}  // namespace

void
TrilinosPreconditioner::compute( Matrix const & mat )
{
  Base::compute( mat );

  if( m_parameters.preconditionerType == "amg" )
  {
    m_precond = CreateMLOperator( m_parameters, mat.unwrapped() );
  }
  else if( m_parameters.preconditionerType == "mgr" )
  {
    GEOSX_ERROR(
      "MGR preconditioner available only through the hypre interface" );
  }
  else if( m_parameters.preconditionerType == "iluk" ||
           m_parameters.preconditionerType == "ilut" ||
           m_parameters.preconditionerType == "icc" ||
           m_parameters.preconditionerType == "ict" )
  {
    m_precond = CreateIfpackOperator( m_parameters, mat.unwrapped() );
  }
  else if( m_parameters.preconditionerType == "none" )
  {
    m_precond.reset();
  }
  else
  {
    GEOSX_ERROR(
      "Unsupported preconditioner type: " << m_parameters.preconditionerType );
  }
}

void
TrilinosPreconditioner::apply( Vector const & src, Vector & dst ) const
{
  GEOSX_LAI_ASSERT( ready() );
  GEOSX_LAI_ASSERT( src.ready() );
  GEOSX_LAI_ASSERT( dst.ready() );
  GEOSX_LAI_ASSERT_EQ( src.globalSize(), this->numGlobalCols() );
  GEOSX_LAI_ASSERT_EQ( dst.globalSize(), this->numGlobalRows() );

  if( m_parameters.preconditionerType == "none" )
  {
    dst.copy( src );
  }
  else
  {
    GEOSX_LAI_ASSERT( m_precond );
    m_precond->ApplyInverse( src.unwrapped(), dst.unwrapped() );
  }
}

void
TrilinosPreconditioner::clear()
{
  PreconditionerBase::clear();
  m_precond.reset();
}

Epetra_Operator const &
TrilinosPreconditioner::unwrapped() const
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

Epetra_Operator &
TrilinosPreconditioner::unwrapped()
{
  GEOSX_LAI_ASSERT( ready() );
  return *m_precond;
}

}  // namespace geosx
