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
 * @file SymbolicFunction.cpp
 */

#include "SymbolicFunction.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
std::string const variableNames = "variableNames";
std::string const expression = "expression";
}  // namespace keys
}  // namespace dataRepository

using namespace dataRepository;

SymbolicFunction::SymbolicFunction( const std::string & name, Group * const parent ) :
  FunctionBase( name, parent )
#ifdef GEOSX_USE_MATHPRESSO
  ,
  parserContext(),
  parserExpression()
#endif
{
  registerWrapper< string_array >( keys::variableNames )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription(
      "List of variables in expression.  The order must match the evaluate "
      "argument" );

  registerWrapper< string >( keys::expression )
    ->setInputFlag( InputFlags::REQUIRED )
    ->setDescription( "Symbolic math expression" );
}

SymbolicFunction::~SymbolicFunction()
{}

void
SymbolicFunction::InitializeFunction()
{
#ifdef GEOSX_USE_MATHPRESSO
  // Register variables
  string_array & variables = getReference< string_array >( keys::variableNames );
  for( localIndex ii = 0; ii < variables.size(); ++ii )
  {
    parserContext.addVariable( variables[ii].c_str(),
                               static_cast< int >( ii * sizeof( double ) ) );
  }

  // Add built in constants/functions (PI, E, sin, cos, ceil, exp, etc.),
  // compile
  parserContext.addBuiltIns();
  std::string const & expression = getReference< std::string >( keys::expression );
  mathpresso::Error err = parserExpression.compile( parserContext,
                                                    expression.c_str(),
                                                    mathpresso::kNoOptions );
  GEOSX_ERROR_IF( err != mathpresso::kErrorOk, "JIT Compiler Error" );
#else
  GEOSX_ERROR( "GEOSX was not built with mathpresso!" );
#endif
}

REGISTER_CATALOG_ENTRY( FunctionBase,
                        SymbolicFunction,
                        std::string const &,
                        Group * const )

}  // namespace geosx
