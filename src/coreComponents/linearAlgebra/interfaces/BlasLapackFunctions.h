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
 * @file BlasLapackFunctions.h
 */

#ifndef GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_
#define GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_

/// This macro provide a flexible interface for Fortran naming convetion for compiled objects
#ifdef FORTRAN_MANGLE_NO_UNDERSCORE
  #define FORTRAN_MANGLE( name ) name
#else
  #define FORTRAN_MANGLE( name ) name##_
#endif

#ifdef __cplusplus
extern "C" {

  #define GEOSX_dasum FORTRAN_MANGLE( dasum )
double
GEOSX_dasum( int const * N, double const * DX, int const * INCX );

  #define GEOSX_daxpy FORTRAN_MANGLE( daxpy )
void
GEOSX_daxpy( int const * N,
             double const * DA,
             double const * DX,
             int const * INCX,
             double * DY,
             int const * INCY );

  #define GEOSX_dcopy FORTRAN_MANGLE( dcopy )
void
GEOSX_dcopy( int const * N,
             double const * DX,
             int const * INCX,
             double * DY,
             int const * INCY );

  #define GEOSX_ddot FORTRAN_MANGLE( ddot )
double
GEOSX_ddot( int const * N,
            double const * DX,
            int const * INCX,
            double const * DY,
            int const * INCY );

  #define GEOSX_idamax FORTRAN_MANGLE( idamax )
int
GEOSX_idamax( int const * N, double const * DX, int const * INCX );

  #define GEOSX_dgemm FORTRAN_MANGLE( dgemm )
void
GEOSX_dgemm( char const * TRANSA,
             char const * TRANSB,
             int const * M,
             int const * N,
             int const * K,
             double const * ALPHA,
             double const * A,
             int const * LDA,
             double const * B,
             int const * LDB,
             double const * BETA,
             double * C,
             int const * LDC );

  #define GEOSX_dgetrf FORTRAN_MANGLE( dgetrf )
void
GEOSX_dgetrf( int const * M,
              int const * N,
              double * A,
              int const * LDA,
              int * IPIV,
              int * INFO );

  #define GEOSX_dgetri FORTRAN_MANGLE( dgetri )
void
GEOSX_dgetri( int const * N,
              double * A,
              int const * LDA,
              int const * IPIV,
              double * WORK,
              int const * LWORK,
              int * INFO );

  #define GEOSX_dlange FORTRAN_MANGLE( dlange )
double
GEOSX_dlange( char const * NORM,
              int const * M,
              int const * N,
              double const * A,
              int const * LDA,
              double * WORK );

  #define GEOSX_dlarnv FORTRAN_MANGLE( dlarnv )
void
GEOSX_dlarnv( int const * IDIST, int * ISEED, int const * N, double * X );

  #define GEOSX_dnrm2 FORTRAN_MANGLE( dnrm2 )
double
GEOSX_dnrm2( int const * N, double const * X, int const * INCX );

  #define GEOSX_dscal FORTRAN_MANGLE( dscal )
void
GEOSX_dscal( int const * N, double const * DA, double * DX, int const * INCX );

  #define GEOSX_dgesvd FORTRAN_MANGLE( dgesvd )
void
GEOSX_dgesvd( char const * JOBU,
              char const * JOBVT,
              int const * M,
              int const * N,
              double * A,
              int const * LDA,
              double * S,
              double * U,
              int const * LDU,
              double * VT,
              int const * LDVT,
              double * WKOPT,
              int const * LWORK,
              int * INFO );
}
#endif

#endif  //GEOSX_LINEARALGEBRA_INTERFACES_BLASLAPACKFUNCTIONS_HPP_
