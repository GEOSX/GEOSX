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

#ifndef _KINEMATICS_H_
#define _KINEMATICS_H_

#include "common/DataTypes.hpp"
#include "common/GeosxMacros.hpp"
#include "LvArray/src/tensorOps.hpp"

//*****************************************************************************
//***** DECLARATIONS **********************************************************
//*****************************************************************************

namespace geosx
{
template< int N, typename ARRAY_2D >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE void CalculateGradients(
  real64 ( &gradient0 )[3][3],
  real64 ( &gradient1 )[3][3],
  real64 const ( &var0 )[N][3],
  real64 const ( &var1 )[N][3],
  ARRAY_2D const & dNdX )
{
  LvArray::tensorOps::AiBj< 3, 3 >( gradient0, var0[0], dNdX[0] );
  LvArray::tensorOps::AiBj< 3, 3 >( gradient1, var1[0], dNdX[0] );

  for( int a = 1; a < N; ++a )
  {
    LvArray::tensorOps::plusAiBj< 3, 3 >( gradient0, var0[a], dNdX[a] );
    LvArray::tensorOps::plusAiBj< 3, 3 >( gradient1, var1[a], dNdX[a] );
  }
}

GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void HughesWinget( real64 ( &Rot )[3][3], real64 ( &Dadt )[6], real64 const ( &G )[3][3] )
{
  //Dadt = 0.5*(G + GT);
  Dadt[0] = G[0][0];
  Dadt[1] = G[1][1];
  Dadt[2] = G[2][2];

  Dadt[3] = 0.5 * ( G[2][1] + G[1][2] );
  Dadt[4] = 0.5 * ( G[2][0] + G[0][2] );
  Dadt[5] = 0.5 * ( G[0][1] + G[1][0] );

  //Omega = 0.5*(G - GT);
  real64 const w12 = 0.5 * ( G[0][1] - G[1][0] );
  real64 const w13 = 0.5 * ( G[0][2] - G[2][0] );
  real64 const w23 = 0.5 * ( G[1][2] - G[2][1] );

  real64 const w12w12div4 = 0.25 * w12 * w12;
  real64 const w13w13div4 = 0.25 * w13 * w13;
  real64 const w23w23div4 = 0.25 * w23 * w23;
  real64 const w12w13div2 = 0.5 * ( w12 * w13 );
  real64 const w12w23div2 = 0.5 * ( w12 * w23 );
  real64 const w13w23div2 = 0.5 * ( w13 * w23 );
  real64 const invDetIplusOmega =
    1.0 / ( 1 + ( w12w12div4 + w13w13div4 + w23w23div4 ) );

  Rot[0][0] = ( 1.0 + ( -w12w12div4 - w13w13div4 + w23w23div4 ) ) * invDetIplusOmega;
  Rot[0][1] = ( w12 - w13w23div2 ) * invDetIplusOmega;
  Rot[0][2] = ( w13 + w12w23div2 ) * invDetIplusOmega;

  Rot[1][0] = ( -w12 - w13w23div2 ) * invDetIplusOmega;
  Rot[1][1] = ( 1.0 + ( -w12w12div4 + w13w13div4 - w23w23div4 ) ) * invDetIplusOmega;
  Rot[1][2] = ( w23 - w12w13div2 ) * invDetIplusOmega;

  Rot[2][0] = ( -w13 + w12w23div2 ) * invDetIplusOmega;
  Rot[2][1] = ( -w23 - w12w13div2 ) * invDetIplusOmega;
  Rot[2][2] = ( 1.0 + ( w12w12div4 - w13w13div4 - w23w23div4 ) ) * invDetIplusOmega;
}
}  // namespace geosx

#endif
