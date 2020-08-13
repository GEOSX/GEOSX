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
 * @file SymbolicFunction.hpp
 */

#ifndef GEOSX_MANAGERS_FUNCTIONS_SYMBOLICFUNCTION_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_SYMBOLICFUNCTION_HPP_

#include "FunctionBase.hpp"

#ifdef GEOSX_USE_MATHPRESSO
  #include <mathpresso/mathpresso.h>
#endif

namespace geosx
{
/**
 * @class SymbolicFunction
 *
 * An interface for an arbitrary symbolic function
 */
class SymbolicFunction : public FunctionBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  SymbolicFunction( const std::string & name, dataRepository::Group * const parent );

  /**
   * @brief The destructor
   */
  virtual ~SymbolicFunction() override;

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string
  CatalogName()
  {
    return "SymbolicFunction";
  }

  /**
   * @brief Initialize the table function
   */
  virtual void
  InitializeFunction() override;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  inline void
  Evaluate( dataRepository::Group const * const group,
            real64 const time,
            SortedArrayView< localIndex const > const & set,
            real64_array & result ) const override final
  {
    FunctionBase::EvaluateT< SymbolicFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  inline real64
  Evaluate( real64 const * const input ) const override final
  {
#ifdef GEOSX_USE_MATHPRESSO
    return parserExpression.evaluate(
      reinterpret_cast< void * >( const_cast< real64 * >( input ) ) );
#else
    GEOSX_ERROR( "GEOSX was not built with mathpresso!" );
    return 0;
#endif
  }

private:
  // Symbolic math driver objects
#ifdef GEOSX_USE_MATHPRESSO
  mathpresso::Context parserContext;
  mathpresso::Expression parserExpression;
#endif
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_SYMBOLICFUNCTION_HPP_ */
