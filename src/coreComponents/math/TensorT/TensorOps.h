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

#ifndef TENSOR_OPS_H_
#define TENSOR_OPS_H_

#include <cmath>

typedef double realT;

inline int
e_ijk( const int i, const int j, const int k )
{
  int r_value = -10000;

  // check to see if any of the indices are the same
  if( ( i == j ) || ( i == k ) || ( j == k ) )
    r_value = 0;
  else if( i == 1 )
  {
    if( j == 2 )
      r_value = 1;
    else
      r_value = -1;
  }
  else if( j == 1 )
  {
    if( i == 2 )
      r_value = -1;
    else
      r_value = 1;
  }
  else if( k == 1 )
  {
    if( i == 2 )
      r_value = 1;
    else
      r_value = -1;
  }
  return r_value;
}

inline int
d_ij( const int i, const int j )
{
  if( i == j )
    return 1;
  else
    return 0;
}

inline int
is_equal( const realT & a, const realT & b, const int num_digits )
{
  int rval = 0;
  const realT mean = ( a + b ) * 0.5;
  const realT diff = fabs( a - b );
  if( diff < mean * pow( 10.0, -num_digits ) ) rval = 1;

  return rval;
}

#endif /* _CONSTANTS_H_ */
