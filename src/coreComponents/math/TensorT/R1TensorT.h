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
 * @brief This file contains the definition of the R1TensorT class
 * @file R1TensorT.h
 */

#ifndef R1_TENSOR_T_H_
#define R1_TENSOR_T_H_

#include "TensorBaseT.h"
#include "LvArray/src/ArraySlice.hpp"

#include <cstdlib>

/**
 * @brief R1TensorT is a rank-1 tensor object type
 * @tparam T_dim length of tensor index
 *
 * R1TensorT derives from TensorBaseT, and defines basic operations that can be
 * done on a rank-1 tensor, as well as operations that result in a rank-1
 * tensor.
 */
template< int T_dim >
class R1TensorT : public TensorBaseT< T_dim >
{
public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  R1TensorT() = default;
  R1TensorT( R1TensorT const & ) = default;
  ~R1TensorT() = default;

  /**
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const realT data ) :
    TensorBaseT< T_dim >( data )
  {}

  /**
   * @param[in] data naked array used for initialization of t_data
   */
  explicit R1TensorT( realT * const data ) :
    TensorBaseT< T_dim >( data )
  {}

  /**
   * @param[in] data use for initialization of t_data
   */
  explicit R1TensorT( const int data ) :
    TensorBaseT< T_dim >( realT( data ) )
  {}

  template< int USD >
  R1TensorT( LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src ) :
    TensorBaseT< T_dim >()
  {
    *this = src;
  }

  //**** CONSTRUCTORS AND DESTRUCTORS
  // *******************************************

  /**
   * Explicit constructors - will throw compile-time errors if not called with
   * the correct dimension
   */
  R1TensorT( realT x, realT y );  //2D only

  GEOSX_HOST_DEVICE
  R1TensorT( realT x, realT y, realT z );  //3D only

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  R1TensorT &
  operator=( const int & rhs );

  /// assignment to all data to a realT
  GEOSX_HOST_DEVICE
  R1TensorT &
  operator=( const realT & rhs );

  /// assignment to another R1TensorT
  R1TensorT &
  operator=( const R1TensorT & rhs ) = default;

  template< int USD >
  GEOSX_HOST_DEVICE CONSTEXPR_WITH_NDEBUG inline R1TensorT &
  operator=(
    LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), T_dim );

    for( int i = 0; i < T_dim; ++i )
    {
      this->t_data[i] = src[i];
    }

    return *this;
  }

  using TensorBaseT< T_dim >::operator+=;

  template< int USD >
  GEOSX_HOST_DEVICE CONSTEXPR_WITH_NDEBUG inline R1TensorT &
  operator+=(
    LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), T_dim );

    for( int i = 0; i < T_dim; ++i )
    {
      this->t_data[i] += src[i];
    }

    return *this;
  }

  template< int USD >
  GEOSX_HOST_DEVICE constexpr inline R1TensorT &
  operator+=(
    LvArray::ArraySlice< realT, 1, USD, std::ptrdiff_t > const & src )
  {
    return ( *this ) +=
      reinterpret_cast< LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & >(
             src );
  }

  using TensorBaseT< T_dim >::operator-=;

  template< int USD >
  GEOSX_HOST_DEVICE CONSTEXPR_WITH_NDEBUG inline R1TensorT &
  operator-=(
    LvArray::ArraySlice< realT const, 1, USD, std::ptrdiff_t > const & src )
  {
    GEOSX_ASSERT_EQ( src.size(), T_dim );

    for( int i = 0; i < T_dim; ++i )
    {
      this->t_data[i] -= src[i];
    }

    return *this;
  }

  //***** ACCESS OPERATORS ****************************************************
  /// const access to data
  GEOSX_HOST_DEVICE inline realT
  operator()( const int i ) const
  {
    return this->t_data[i];
  }

  /// non-const access to data
  GEOSX_HOST_DEVICE inline realT &
  operator()( const int i )
  {
    return this->t_data[i];
  }

  /// const access to data
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE realT
  operator[]( const int i ) const
  {
    return this->t_data[i];
  }

  /// non-const access to data
  GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE realT &
  operator[]( const int i )
  {
    return this->t_data[i];
  }

  realT
  ProductOfSquares() const;

  /// Hadamard product between two Rank1 tensors
  void
  AiBi( const R1TensorT< T_dim > & A, const R1TensorT< T_dim > & B );

  /// cross product of 2 rank1 tensors
  GEOSX_HOST_DEVICE
  void
  Cross( const R1TensorT< T_dim > & a, const R1TensorT< T_dim > & b );

  //****** TENSOR OPERATIONS **************************************************
  /// take the L2 norm of the tensor
  GEOSX_HOST_DEVICE
  realT
  L2_Norm( void ) const;

  /// get the unit vector
  R1TensorT< T_dim >
  UnitVector( void ) const
  {
    realT n = this->L2_Norm();
    return ( n > 0.0 ) ? ( *this / n ) : *this;
  }

  /// Normalize the vector
  GEOSX_HOST_DEVICE
  realT
  Normalize( void )
  {
    realT n = this->L2_Norm();
    if( n > 0.0 ) *this /= n;
    return n;
  }

  /// sum the components of the tensor
  inline realT
  Sum( void ) const;

  //***** OUTPUT **************************************************************
  /// output
  void
  print( std::ostream & os ) const;

  // define cross product
  friend inline GEOSX_HOST_DEVICE R1TensorT< T_dim >
  Cross( const R1TensorT< T_dim > & a,
         const R1TensorT< T_dim > & b )
  {
    R1TensorT< T_dim > c;
    c.Cross( a, b );
    return c;
  }

private:
};

//*****************************************************************************
//***** END DECLARATION *******************************************************
//*****************************************************************************

/// Explicit 2D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R1TensorT< 2 >::R1TensorT( realT x, realT y ) :
  TensorBaseT< 2 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
}

/// Explicit 3D constructor
///
/// Template specialisation - if templated on another dimension constructor will
// throw a compile time error.
template<>
inline R1TensorT< 3 >::R1TensorT( realT x, realT y, realT z ) :
  TensorBaseT< 3 >()
{
  this->t_data[0] = x;
  this->t_data[1] = y;
  this->t_data[2] = z;
}

template< int T_dim >
void
R1TensorT< T_dim >::print( std::ostream & os ) const
{
  for( int i = 0; i < T_dim; ++i ) os << ( *this )( i ) << '\t';
}

//*****************************************************************************
//***** R1TensorT Member Function Definition **********************************
//*****************************************************************************

//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template< int T_dim >
inline R1TensorT< T_dim > &
R1TensorT< T_dim >::operator=( const int & rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

/**
 * @param[in] rhs value to set each member of t_data to
 * @return reference to *this
 */
template< int T_dim >
GEOSX_HOST_DEVICE inline R1TensorT< T_dim > &
R1TensorT< T_dim >::operator=(
  const realT & rhs )
{
  TensorBaseT< T_dim >::operator=( rhs );
  return *this;
}

//***** TENSOR OPERATORS ******************************************************

/**
 * @return L2 norm of tensor
 */
template< int T_dim >
GEOSX_HOST_DEVICE inline realT
R1TensorT< T_dim >::L2_Norm( void ) const
{
  realT norm = 0.0;

  for( int i = 1; i <= T_dim; ++i )
    norm += this->t_data[i - 1] * this->t_data[i - 1];
  norm = sqrt( norm );

  return norm;
}

/**
 * @return sum of the tensor components
 */
template< int T_dim >
inline realT
R1TensorT< T_dim >::Sum( void ) const
{
  realT sum = 0.0;

  for( int i = 1; i <= T_dim; ++i ) sum += this->t_data[i - 1];

  return sum;
}

//***** MULTIPLICATION OPERATORS **********************************************

/**
 * @param[in] A rank2 tensor
 * @param[in] B rank1 tensor
 * @return none
 */
template< int T_dim >
inline void
R1TensorT< T_dim >::AiBi( const R1TensorT< T_dim > & A,
                          const R1TensorT< T_dim > & B )
{
  for( int i = 0; i < T_dim; i++ ) this->t_data[i] = A.t_data[i] * B.t_data[i];
}

/**
 * @param[in] a rank1 tensor
 * @param[in] b rank1 tensor
 * @return none
 *
 * this function takes the cross product of two rank1 tensors and places the
 * result into this->tdata
 */
template< int T_dim >
GEOSX_HOST_DEVICE inline void
R1TensorT< T_dim >::Cross( const R1TensorT< T_dim > & a,
                           const R1TensorT< T_dim > & b )
{
  if( T_dim == 3 )
  {
    this->t_data[0] = a.t_data[1] * b.t_data[2] - a.t_data[2] * b.t_data[1];
    this->t_data[1] = -( a.t_data[0] * b.t_data[2] - a.t_data[2] * b.t_data[0] );
    this->t_data[2] = a.t_data[0] * b.t_data[1] - a.t_data[1] * b.t_data[0];
  }
  //  else
  //    std::cout << "R1TensorT not implemented for nsdof>3";
}

#endif
