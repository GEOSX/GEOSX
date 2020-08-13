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

/*
 * File: Interpolation.hpp
 */

#ifndef INTERPOLATION_HPP_
#define INTERPOLATION_HPP_

namespace geosx
{
namespace interpolation
{
static real64
ParabolicInterpolationThreePoints( real64 const lambdac,
                                   real64 const lambdam,
                                   real64 const ff0,
                                   real64 const ffT,
                                   real64 const ffm )
{
  // Apply three-point safeguarded parabolic model for a line search.
  //
  // C. T. Kelley, April 1, 2003
  //
  // This code comes with no guarantee or warranty of any kind.
  //
  // function lambdap = parab3p(lambdac, lambdam, ff0, ffT, ffm)
  //
  // input:
  //       lambdac = current steplength
  //       lambdam = previous steplength
  //       ff0 = value of \| F(x_c) \|^2
  //       ffc = value of \| F(x_c + \lambdac d) \|^2
  //       ffm = value of \| F(x_c + \lambdam d) \|^2
  //
  // output:
  //       lambdap = new value of lambda given parabolic model
  //
  // internal parameters:
  //       sigma0 = .1, sigma1 = .5, safeguarding bounds for the linesearch

  // Set internal parameters.
  real64 const sigma0 = 0.1;
  real64 const sigma1 = 0.5;

  // Compute coefficients of interpolation polynomial.
  // p(lambda) = ff0 + (c1 lambda + c2 lambda^2)/d1
  // d1 = (lambdac - lambdam)*lambdac*lambdam < 0
  //      so, if c2 > 0 we have negative curvature and default to
  //      lambdap = sigma1 * lambda.
  real64 const c2 = lambdam * ( ffT - ff0 ) - lambdac * ( ffm - ff0 );
  if( c2 >= 0.0 )
  {
    return ( sigma1 * lambdac );
  }
  real64 const c1 =
    lambdac * lambdac * ( ffm - ff0 ) - lambdam * lambdam * ( ffT - ff0 );
  real64 lambdap = -c1 * 0.5;
  if( lambdap > sigma0 * lambdac * c2 )
  {
    lambdap = sigma0 * lambdac * c2;
  }
  if( lambdap < sigma1 * lambdac * c2 )
  {
    lambdap = sigma1 * lambdac * c2;
  }
  lambdap /= c2;
  return lambdap;
}

}  // namespace interpolation

}  // namespace geosx

#endif  // INTERPOLATION_HPP
