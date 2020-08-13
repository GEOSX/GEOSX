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
 * @file UtilityFunctions.cpp
 */

#include "constitutive/fluid/PVTFunctions/UtilityFunctions.hpp"

namespace geosx
{
namespace PVTProps
{
EvalArgs2D
XYTable::Value( EvalArgs2D const & x, EvalArgs2D const & y ) const
{
  localIndex_array xIndices( 2 ), yIndices( 2 );

  //find i index

  if( x <= m_x[0] )
  {
    xIndices[0] = 0;
  }
  else if( x >= m_x[m_x.size() - 2] )
  {
    xIndices[0] = m_x.size() - 2;
  }
  else
  {
    for( localIndex i = 1; i < m_x.size(); ++i )
    {
      if( x <= m_x[i] )
      {
        xIndices[0] = i - 1;
        break;
      }
    }
  }

  xIndices[1] = xIndices[0] + 1;

  EvalArgs2D xWeight =
    ( x - m_x[xIndices[0]] ) / ( m_x[xIndices[1]] - m_x[xIndices[0]] );

  //find j index

  if( y <= m_y[0] )
  {
    yIndices[0] = 0;
  }
  else if( y >= m_y[m_y.size() - 2] )
  {
    yIndices[0] = m_y.size() - 2;
  }
  else
  {
    for( localIndex i = 1; i < m_y.size(); ++i )
    {
      if( y <= m_y[i] )
      {
        yIndices[0] = i - 1;
        break;
      }
    }
  }

  yIndices[1] = yIndices[0] + 1;

  EvalArgs2D yWeight =
    ( y - m_y[yIndices[0]] ) / ( m_y[yIndices[1]] - m_y[yIndices[0]] );

  EvalArgs2D out =
    m_value[xIndices[0]][yIndices[0]] * ( 1.0 - xWeight ) * ( 1.0 - yWeight ) +
    m_value[xIndices[0]][yIndices[1]] * ( 1.0 - xWeight ) * yWeight +
    m_value[xIndices[1]][yIndices[1]] * xWeight * yWeight +
    m_value[xIndices[1]][yIndices[0]] * xWeight * ( 1.0 - yWeight );

  return out;
}

template< class T >
T
XTable::GetValue( T const & x ) const
{
  localIndex idx = 0;
  //find i index

  if( x <= m_x[0] )
  {
    idx = 0;
  }
  else if( x >= m_x[m_x.size() - 1] )
  {
    idx = m_x.size() - 2;
  }
  else
  {
    for( localIndex i = 1; i < m_x.size(); ++i )
    {
      if( x <= m_x[i] )
      {
        idx = i - 1;
        break;
      }
    }
  }

  T weight = ( x - m_x[idx] ) / ( m_x[idx + 1] - m_x[idx] );

  T out = m_value[idx] * ( 1.0 - weight ) + m_value[idx + 1] * weight;

  return out;
}

}  // namespace PVTProps
}  // namespace geosx
