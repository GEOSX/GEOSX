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
 * @file FunctionBase.cpp
 */

#include "FunctionBase.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
using namespace dataRepository;

FunctionBase::FunctionBase( const std::string & name, Group * const parent ) :
  Group( name, parent ),
  m_inputVarNames()
{
  setInputFlags( InputFlags::OPTIONAL_NONUNIQUE );

  registerWrapper( keys::inputVarNames, &m_inputVarNames )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setSizedFromParent( 0 )
    ->setDescription( "Name of fields are input to function." );
}

FunctionBase::~FunctionBase()
{}

integer
FunctionBase::isFunctionOfTime() const
{
  integer rval = 0;
  arrayView1d< string const > const & inputVarNames =
    this->getReference< string_array >( dataRepository::keys::inputVarNames );
  localIndex numVars =
    LvArray::integerConversion< localIndex >( inputVarNames.size() );

  if( numVars == 1 )
  {
    if( inputVarNames[0] == "time" )
    {
      rval = 2;
    }
  }
  else
  {
    for( auto varIndex = 0; varIndex < numVars; ++varIndex )
    {
      if( inputVarNames[varIndex] == "time" )
      {
        rval = 1;
        break;
      }
    }
  }
  return rval;
}

real64_array
FunctionBase::EvaluateStats( dataRepository::Group const * const group,
                             real64 const time,
                             SortedArray< localIndex > const & set ) const
{
  localIndex N = set.size();
  real64_array sub( N );
  Evaluate( group, time, set.toViewConst(), sub );

  real64_array result( 3 );
  result[0] = 1e10;   // min
  result[1] = 0.0;    // avg
  result[2] = -1e10;  // max
  for( localIndex ii = 0; ii < N; ii++ )
  {
    result[0] = std::min( result[0], sub[ii] );
    result[1] += sub[ii];
    result[2] = std::max( result[0], sub[ii] );
  }
  result[1] /= N;

  return result;
}

}  // namespace geosx
