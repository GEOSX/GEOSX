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
 * @file FunctionBase.hpp
 */

#ifndef GEOSX_MANAGERS_FUNCTIONS_FUNCTIONBASE_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_FUNCTIONBASE_HPP_

#include "dataRepository/Group.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"

namespace geosx
{
namespace dataRepository
{
namespace keys
{
/**
 * @brief The key for inputVarNames
 * @return the key
 */
string const inputVarNames( "inputVarNames" );
}  // namespace keys
}  // namespace dataRepository

/**
 * @class FunctionBase
 *
 * An object for interfacing with arbitrary N-dimensional functions.
 */
class FunctionBase : public dataRepository::Group
{
public:
  /// @copydoc geosx::dataRepository::Group::Group( std::string const & name, Group * const parent )
  FunctionBase( const std::string & name, dataRepository::Group * const parent );

  /**
   * @brief destructor
   */
  virtual ~FunctionBase() override;

  /**
   * @brief Static Factory Catalog Functions
   * @return the catalog name
   */
  static string
  CatalogName()
  {
    return "FunctionBase";
  }

  /**
   * @brief Function initialization
   */
  virtual void
  InitializeFunction()
  {}

  /**
   * @brief Test to see if the function is a 1D function of time
   * @return integer value:
   *         - 0 is the function does not have time as parameter
   *         - 1 is the function has time as one of the parameters
   *         - 2 is the function has time as only parameter
   */
  integer
  isFunctionOfTime() const;

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void
  Evaluate( dataRepository::Group const * const group,
            real64 const time,
            SortedArrayView< localIndex const > const & set,
            real64_array & result ) const = 0;

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function evaluation
   */
  virtual real64
  Evaluate( real64 const * const input ) const = 0;

  /// Alias for the catalog interface
  using CatalogInterface =
    dataRepository::CatalogInterface< FunctionBase, std::string const &, Group * const >;

  /**
   * @brief return the catalog entry for the function
   * @return the catalog entry
   */
  static CatalogInterface::CatalogType &
  GetCatalog()
  {
    static CatalogInterface::CatalogType catalog;
    return catalog;
  }

  /**
   * @brief This generates statistics by applying a function to an object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @return An array holding the min, average, max values of the results
   */
  real64_array
  EvaluateStats( dataRepository::Group const * const group,
                 real64 const time,
                 SortedArray< localIndex > const & set ) const;

protected:
  /// names for the input variables
  string_array m_inputVarNames;

  /**
   * @brief Method to apply an function with an arbitrary type of output
   * @tparam LEAF the return type
   * @param[in] group a pointer to the object holding the function arguments
   * @param[in] time current time
   * @param[in] set the subset of nodes to apply the function to
   * @param[out] result the results
   */
  template< typename LEAF >
  void
  EvaluateT( dataRepository::Group const * const group,
             real64 const time,
             SortedArrayView< localIndex const > const & set,
             real64_array & result ) const;

  virtual void
  PostProcessInput() override
  {
    InitializeFunction();
  }
};

template< typename LEAF >
void
FunctionBase::EvaluateT( dataRepository::Group const * const group,
                         real64 const time,
                         SortedArrayView< localIndex const > const & set,
                         real64_array & result ) const
{
  real64 const * input_ptrs[4];

  arrayView1d< string const > const & inputVarNames =
    this->getReference< string_array >( dataRepository::keys::inputVarNames );

  localIndex const numVars =
    LvArray::integerConversion< localIndex >( inputVarNames.size() );
  GEOSX_ERROR_IF( numVars > 4, "Number of variables is: " << numVars );

  localIndex varSize[4];
  int timeVar[4] = { 1, 1, 1, 1 };
  for( auto varIndex = 0; varIndex < numVars; ++varIndex )
  {
    string const & varName = inputVarNames[varIndex];

    if( varName == "time" )
    {
      input_ptrs[varIndex] = &time;
      varSize[varIndex] = 1;
      timeVar[varIndex] = 0;
    }
    else
    {
      dataRepository::WrapperBase const & wrapper =
        *( group->getWrapperBase( varName ) );
      input_ptrs[varIndex] =
        reinterpret_cast< double const * >( wrapper.voidPointer() );
      varSize[varIndex] = wrapper.elementByteSize() / sizeof( double );
    }
  }

  integer count = 0;
  forAll< serialPolicy >( set.size(), [&, set]( localIndex const i ) {
    localIndex const index = set[i];
    double input[4];
    int c = 0;
    for( int a = 0; a < numVars; ++a )
    {
      for( int b = 0; b < varSize[a]; ++b )
      {
        input[c] = input_ptrs[a][( index * varSize[a] + b ) * timeVar[a]];
        ++c;
      }
    }

    // TODO: Check this line to make sure it is correct
    result[count] = static_cast< LEAF const * >( this )->Evaluate( input );
    ++count;
  } );
}
} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_FUNCTIONBASE_HPP_ */
