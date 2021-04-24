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
 * @file SolidModelHelperOrthotropic.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSORTHOTROPIC_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSORTHOTROPIC_HPP_

#include "SolidModelDiscretizationOps.hpp"

namespace geosx
{
namespace constitutive
{

struct SolidModelDiscretizationOpsOrthotropic : public SolidModelDiscretizationOps
{
  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void upperBTDB( BASIS_GRADIENT const & gradN,
                  real64 const & detJxW,
                  real64 ( &elementStiffness )[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagBTDB( BASIS_GRADIENT const & gradN,
                 real64 const & detJxW,
                 real64 ( &diagElementStiffness )[NUM_SUPPORT_POINTS*3] );

  template< int NUM_SUPPORT_POINTS,
            typename BASIS_GRADIENT >
  GEOSX_HOST_DEVICE
  void diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                       real64 const & detJxW,
                       real64 ( &diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] );

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void scaleParams( real64 const scale )
  {
    m_c11 *= scale;
    m_c12 *= scale;
    m_c13 *= scale;
    m_c22 *= scale;
    m_c23 *= scale;
    m_c33 *= scale;
    m_c44 *= scale;
    m_c55 *= scale;
    m_c66 *= scale;
  }

  real64 m_c11;
  real64 m_c12;
  real64 m_c13;
  real64 m_c22;
  real64 m_c23;
  real64 m_c33;
  real64 m_c44;
  real64 m_c55;
  real64 m_c66;

};

#if __GNUC__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wshadow"
#endif

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsOrthotropic::
  upperBTDB( BASIS_GRADIENT const & gradN,
             real64 const & detJxW,
             real64 (& elementStiffness)[NUM_SUPPORT_POINTS *3][NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c12 = this->m_c12 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c22 = this->m_c22 * detJxW;
  real64 const c23 = this->m_c23 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c55 = this->m_c55 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelDiscretizationOps::upperBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                elementStiffness,
                                                                [ c11,
                                                                  c12,
                                                                  c13,
                                                                  c22,
                                                                  c23,
                                                                  c33,
                                                                  c44,
                                                                  c55,
                                                                  c66 ] GEOSX_HOST_DEVICE
                                                                  ( int const a,
                                                                  int const b,
                                                                  real64 const (&gradNa_gradNb)[3][3],
                                                                  real64 (& elementStiffness)[NUM_SUPPORT_POINTS*3][NUM_SUPPORT_POINTS*3] )
  {
    elementStiffness[a*3+0][b*3+0] += c11 * gradNa_gradNb[0][0] + c66 * gradNa_gradNb[1][1] + c55 * gradNa_gradNb[2][2];
    elementStiffness[a*3+0][b*3+1] += c12 * gradNa_gradNb[0][1] + c66 * gradNa_gradNb[1][0];
    elementStiffness[a*3+0][b*3+2] += c13 * gradNa_gradNb[0][2] + c55 * gradNa_gradNb[2][0];
    elementStiffness[a*3+1][b*3+0] += c66 * gradNa_gradNb[0][1] + c12 * gradNa_gradNb[1][0];
    elementStiffness[a*3+1][b*3+1] += c66 * gradNa_gradNb[0][0] + c22 * gradNa_gradNb[1][1] + c44 * gradNa_gradNb[2][2];
    elementStiffness[a*3+1][b*3+2] += c23 * gradNa_gradNb[1][2] + c44 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+0] += c55 * gradNa_gradNb[0][2] + c13 * gradNa_gradNb[2][0];
    elementStiffness[a*3+2][b*3+1] += c44 * gradNa_gradNb[1][2] + c23 * gradNa_gradNb[2][1];
    elementStiffness[a*3+2][b*3+2] += c55 * gradNa_gradNb[0][0] + c44 * gradNa_gradNb[1][1] + c33 * gradNa_gradNb[2][2];
  } );
}


template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsOrthotropic::diagBTDB( BASIS_GRADIENT const & gradN,
                                                       real64 const & detJxW,
                                                       real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS *3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c22 = this->m_c22 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c55 = this->m_c55 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;
  SolidModelDiscretizationOps::diagBTDB< NUM_SUPPORT_POINTS,
                                         3 >( gradN,
                                              diagElementStiffness,
                                              [ c11,
                                                c22,
                                                c33,
                                                c44,
                                                c55,
                                                c66 ] GEOSX_HOST_DEVICE
                                                ( const int a,
                                                const real64 (& gradN_gradN)[3],
                                                real64 (& diagElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagElementStiffness[a*3+0] += c11 * gradN_gradN[0] + c66 * gradN_gradN[1] + c55 * gradN_gradN[2];
    diagElementStiffness[a*3+1] += c66 * gradN_gradN[0] + c22 * gradN_gradN[1] + c44 * gradN_gradN[2];
    diagElementStiffness[a*3+2] += c55 * gradN_gradN[0] + c44 * gradN_gradN[1] + c33 * gradN_gradN[2];
  } );
}

template< int NUM_SUPPORT_POINTS,
          typename BASIS_GRADIENT >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SolidModelDiscretizationOpsOrthotropic::diagRowSumBTDB( BASIS_GRADIENT const & gradN,
                                                             real64 const & detJxW,
                                                             real64 ( & diagSumElementStiffness )[NUM_SUPPORT_POINTS*3] )
{
  real64 const c11 = this->m_c11 * detJxW;
  real64 const c12 = this->m_c12 * detJxW;
  real64 const c13 = this->m_c13 * detJxW;
  real64 const c22 = this->m_c22 * detJxW;
  real64 const c23 = this->m_c23 * detJxW;
  real64 const c33 = this->m_c33 * detJxW;
  real64 const c44 = this->m_c44 * detJxW;
  real64 const c55 = this->m_c55 * detJxW;
  real64 const c66 = this->m_c66 * detJxW;

  SolidModelDiscretizationOps::diagRowSumBTDB< NUM_SUPPORT_POINTS >( gradN,
                                                                     diagSumElementStiffness,
                                                                     [ c11,
                                                                       c12,
                                                                       c13,
                                                                       c22,
                                                                       c23,
                                                                       c33,
                                                                       c44,
                                                                       c55,
                                                                       c66 ] GEOSX_HOST_DEVICE
                                                                       ( int const a,
                                                                       real64 const (&gradNa_gradNb)[3][3],
                                                                       real64 (& diagSumElementStiffness)[NUM_SUPPORT_POINTS*3] )
  {
    diagSumElementStiffness[a*3+0] = diagSumElementStiffness[a*3+0] +
                                     c11 * gradNa_gradNb[0][0] +
                                     c66 * gradNa_gradNb[1][1] +
                                     c55 * gradNa_gradNb[2][2] +
                                     c12 * gradNa_gradNb[0][1] +
                                     c66 * gradNa_gradNb[1][0] +
                                     c13 * gradNa_gradNb[0][2] +
                                     c55 * gradNa_gradNb[2][0];

    diagSumElementStiffness[a*3+1] = diagSumElementStiffness[a*3+1] +
                                     c66 * gradNa_gradNb[0][1] +
                                     c12 * gradNa_gradNb[1][0] +
                                     c66 * gradNa_gradNb[0][0] +
                                     c22 * gradNa_gradNb[1][1] +
                                     c44 * gradNa_gradNb[2][2] +
                                     c23 * gradNa_gradNb[1][2] +
                                     c44 * gradNa_gradNb[2][1];

    diagSumElementStiffness[a*3+2] = diagSumElementStiffness[a*3+2] +
                                     c55 * gradNa_gradNb[0][2] +
                                     c13 * gradNa_gradNb[2][0] +
                                     c44 * gradNa_gradNb[1][2] +
                                     c23 * gradNa_gradNb[2][1] +
                                     c55 * gradNa_gradNb[0][0] +
                                     c44 * gradNa_gradNb[1][1] +
                                     c33 * gradNa_gradNb[2][2];
  } );
}

#if __GNUC__
#pragma GCC diagnostic pop
#endif

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_SOLIDMODELDISCRETIZATIONOPSORTHOTROPIC_HPP_ */
