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
 * @file CompositeFunction.cpp
 */

#include "FunctionManager.hpp"
#include "CompositeFunction.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
std::string const functionNames = "functionNames";
std::string const variableNames = "variableNames";
std::string const expression = "expression";
}  // namespace keys
}  // namespace dataRepository

using namespace dataRepository;

CompositeFunction::CompositeFunction( const std::string & name, Group * const parent ) :
  FunctionBase( name, parent ),
#ifdef GEOSX_USE_MATHPRESSO
  parserContext(),
  parserExpression(),
#endif
  m_numSubFunctions(),
  m_subFunctions()
{
  registerWrapper( keys::functionNames, &m_functionNames )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "List of source functions. The order must match the variableNames "
      "argument." );

  registerWrapper( keys::variableNames, &m_variableNames )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "List of variables in expression" );

  registerWrapper( keys::expression, &m_expression )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Composite math expression" );
}

CompositeFunction::~CompositeFunction()
{}

void
CompositeFunction::InitializeFunction()
{
#ifdef GEOSX_USE_MATHPRESSO
  // Register variables
  for( localIndex ii = 0; ii < m_variableNames.size(); ++ii )
  {
    parserContext.addVariable( m_variableNames[ii].c_str(),
                               static_cast< int >( ii * sizeof( double ) ) );
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.),
  // compile
  parserContext.addBuiltIns();
  mathpresso::Error err = parserExpression.compile( parserContext,
                                                    m_expression.c_str(),
                                                    mathpresso::kNoOptions );
  GEOSX_ERROR_IF( err != mathpresso::kErrorOk, "JIT Compiler Error" );

  // Grab pointers to sub functions
  FunctionManager & functionManager = FunctionManager::Instance();
  m_numSubFunctions =
    LvArray::integerConversion< localIndex >( m_functionNames.size() );
  for( localIndex ii = 0; ii < m_numSubFunctions; ++ii )
  {
    m_subFunctions.emplace_back(
      functionManager.GetGroup< FunctionBase >( m_functionNames[ii] ) );
  }
#else
  GEOSX_ERROR( "GEOSX was not configured with mathpresso!" );
#endif
}

void
CompositeFunction::Evaluate( dataRepository::Group const * const group,
                             real64 const time,
                             SortedArrayView< localIndex const > const & set,
                             real64_array & result ) const
{
#ifdef GEOSX_USE_MATHPRESSO
  // Evaluate each of the subFunctions independently and place the results into
  // a temporary field
  array1d< real64_array > subFunctionResults;
  for( localIndex ii = 0; ii < m_numSubFunctions; ++ii )
  {
    real64_array tmp( result.size() );
    m_subFunctions[ii]->Evaluate( group, time, set, tmp );
    subFunctionResults.emplace_back( std::move( tmp ) );
  }

  // Evaluate the symbolic math
  forAll< serialPolicy >( set.size(), [&, set]( localIndex const i ) {
    localIndex const ii = set[i];
    real64 functionResults[m_maxNumSubFunctions];
    for( localIndex jj = 0; jj < m_numSubFunctions; ++jj )
    {
      functionResults[jj] = subFunctionResults[jj][ii];
    }
    result[ii] =
      parserExpression.evaluate( reinterpret_cast< void * >( functionResults ) );
  } );
#else
  GEOSX_ERROR( "GEOSX was not configured with mathpresso!" );
#endif
}

real64
CompositeFunction::Evaluate( real64 const * const input ) const
{
#ifdef GEOSX_USE_MATHPRESSO
  real64 functionResults[m_maxNumSubFunctions];

  for( localIndex ii = 0; ii < m_numSubFunctions; ++ii )
  {
    functionResults[ii] = m_subFunctions[ii]->Evaluate( input );
  }

  return parserExpression.evaluate( reinterpret_cast< void * >( functionResults ) );
#else
  GEOSX_ERROR( "GEOSX was not configured with mathpresso!" );
  return 0;
#endif
}

REGISTER_CATALOG_ENTRY( FunctionBase,
                        CompositeFunction,
                        std::string const &,
                        Group * const )

}  // namespace geosx
