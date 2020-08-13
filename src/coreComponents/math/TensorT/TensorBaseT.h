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
 * @brief This file contains the definition of the TensorBaseT class
 * @file TensorBaseT.h
 */

#ifndef TENSOR_BASE_T_H_
#define TENSOR_BASE_T_H_
#include <iostream>
#include <cstring>
#include <string>
#include <sstream>
#include <exception>
#include <limits>
#include "TensorOps.h"
#include "common/Logger.hpp"
#include "common/GeosxMacros.hpp"

/**
 * @brief TensorBaseT is the base class for the tensor library.
 * @tparam T_length length
 *
 * TensorBaseT defines basic operations on the data for use by the derived
 * class.
 */
template< int T_length >
class TensorBaseT
{
  friend std ::istream &
  operator>>( std::istream & in, TensorBaseT< T_length > & t )
  {
    realT * tp = t.Data();
    for( int ii = 0; ii < T_length; ++ii )
    {
      while( in.peek() == ',' || in.peek() == ' ' )
      {
        in.ignore();
      }

      in >> tp[ii];
    }
    return in;
  }

  friend std ::ostream &
  operator<<( std::ostream & out,
              const TensorBaseT< T_length > & t )
  {
    const realT * tp = t.Data();
    for( int ii = 0; ii < T_length; ++ii )
    {
      if( ii > 0 ) out << " ";
      out << tp[ii];
    }
    return out;
  }

public:
  //**** CONSTRUCTORS AND DESTRUCTORS ******************************************

  /// default constructor
  GEOSX_HOST_DEVICE
  TensorBaseT( void );

  /// constructor initialized by single value
  GEOSX_HOST_DEVICE
  explicit TensorBaseT( const realT data );

  /// constructor initialized by raw data
  GEOSX_HOST_DEVICE
  explicit TensorBaseT( const realT data[T_length] );

  /// constructor initialized by another TensorBaseT object
  TensorBaseT( const TensorBaseT< T_length > & rhs ) = default;

  /// non-virtual destructor. This means that this class in NOT intended to be
  /// used as a polymorphically.
  ~TensorBaseT() = default;

  //***** ASSIGNMENT OPERATORS *************************************************
  /// assignment of all data to an integer
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator=( const int & rhs );

  /// assignment to all data to a realT
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator=( const realT & rhs );

  /// assignment to another TensorBaseT
  TensorBaseT &
  operator=( const TensorBaseT & rhs ) = default;

  /// add a realT to data
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator+=( const realT & rhs );

  /// subtract a realT from data
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator-=( const realT & rhs );

  /// multiply each entry in t_data by a realT
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator*=( const realT & rhs );

  /// divide each entry in t_data by a realT
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator/=( const realT & rhs );

  /// add another tensor
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator+=( const TensorBaseT & rhs );

  /// subtract a tensor
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator-=( const TensorBaseT & rhs );

  /// multiply by a tensor (data component by component)
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator*=( const TensorBaseT & rhs );

  /// divide by a tensor (data component by component)
  GEOSX_HOST_DEVICE
  TensorBaseT &
  operator/=( const TensorBaseT & rhs );

  bool
  operator<( const TensorBaseT & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] >= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool
  operator<=( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] > rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool
  operator>( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] <= rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool
  operator>=( const TensorBaseT< T_length > & rhs ) const
  {
    bool rval = true;
    for( int i = 0; i < T_length; ++i )
    {
      if( t_data[i] < rhs.t_data[i] )
      {
        rval = false;
      }
    }
    return rval;
  }

  bool
  operator==( const TensorBaseT< T_length > & rhs ) const
  {
    for( int i = 0; i < T_length; ++i )
    {
      if( ( t_data[i] > rhs.t_data[i] ) || ( t_data[i] < rhs.t_data[i] ) )
      {
        return false;
      }
    }
    return true;
  }

  /// function to add the product of a scalar and tensor
  GEOSX_HOST_DEVICE
  inline void
  plus_cA( const realT & c, const TensorBaseT< T_length > & A )
  {
    for( int i = 0; i < T_length; ++i ) t_data[i] = t_data[i] + c * A.t_data[i];
  }

  /// function to take the product of a scalar and tensor
  GEOSX_HOST_DEVICE
  inline void
  cA( const realT & c, const TensorBaseT< T_length > & A )
  {
    for( int i = 0; i < T_length; ++i ) t_data[i] = c * A.t_data[i];
  }

  /// function to take the quotient of a tensor by a scalar
  GEOSX_HOST_DEVICE
  inline void
  Adivc( const realT & c, const TensorBaseT< T_length > & A )
  {
    for( int i = 0; i < T_length; ++i ) t_data[i] = A.t_data[i] / c;
  }

  //***** OUTPUT **************************************************************

  /// function to cast data array to float
  GEOSX_HOST_DEVICE
  void
  CastDataToFloat( float rval[T_length] ) const
  {
    for( int a = 0; a < T_length; ++a )
    {
      rval[a] = static_cast< float >( t_data[a] );
    }
  }

  inline void
  StrVal( const std::string & str )
  {
    std::istringstream iss( str, std::istringstream::in );
    for( int i = 0; i < T_length; i++ )
      GEOSX_ERROR_IF( !( iss >> t_data[i] ), "Error" );
  }

  /*
   /// ouput function
   virtual void print( ostream& os ) const = 0;

   /// stream function
   friend ostream &operator<<( ostream &os, const TensorBaseT< T_length >& A )
   {
    A.print( os );
    return os;
   }
 */
  //***** DATA MEMBERS ********************************************************
protected:
  /// Tensor data array
  realT t_data[T_length];

  //***** MEMBER ACCESS *******************************************************
public:
  /**
   * @return gives a non-const realT* which points to t_data
   * @brief returns a non-const realT* which points to t_data
   */
  GEOSX_HOST_DEVICE inline constexpr realT *
  Data( void )
  {
    return t_data;
  }

  /**
   * @return gives a const realT* which points to t_data
   * @brief gives a const realT* which points to t_data
   */
  GEOSX_HOST_DEVICE inline constexpr const realT *
  Data( void ) const
  {
    return t_data;
  }

  /**
   * @return gives a non-const realT* which points to t_data
   * @brief returns a non-const realT* which points to t_data
   */
  realT *
  begin( void )
  {
    return t_data;
  }

  /**
   * @return gives a const realT* which points to t_data
   * @brief gives a const realT* which points to t_data
   */
  GEOSX_HOST_DEVICE
  const realT *
  begin( void ) const
  {
    return t_data;
  }

  /**
   * @return gives a non-const realT* which points to the past-the-end element
   * of t_data
   * @brief returns a non-const realT* which points to the past-the-end element
   * of t_data
   */
  GEOSX_HOST_DEVICE
  realT *
  end( void )
  {
    return t_data + T_length;
  }

  /**
   * @return gives a const realT* which points to the past-the-end element of
   * t_data
   * @brief gives a const realT* which points to the past-the-end element of
   * t_data
   */
  GEOSX_HOST_DEVICE
  const realT *
  end( void ) const
  {
    return t_data + T_length;
  }

  /**
   * @return the number of data entries (Length) for the tensor
   * @brief gives the number of data entries (Length) for the tensor
   */
  GEOSX_HOST_DEVICE constexpr static int
  Length( void )
  {
    return T_length;
  }

  /**
   * @return the maximum single value in the t_data
   * @brief gives the maximum single value in the t_data
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  realT
  MaxVal( void ) const
  {
    realT rval = 0;
    for( int i = 0; i < T_length; ++i )
      if( fabs( t_data[i] ) > rval ) rval = fabs( t_data[i] );
    return rval;
  }

  /**
   * @return the minimum single value in the t_data
   * @brief gives the minimum single value in the t_data
   */
  realT
  MinVal( void ) const
  {
    realT rval = std::numeric_limits< realT >::max();
    for( int i = 0; i < T_length; ++i )
      if( fabs( t_data[i] ) < rval ) rval = fabs( t_data[i] );
    return rval;
  }

  void
  SetMax( const TensorBaseT< T_length > & newval )
  {
    for( int i = 0; i < T_length; ++i )
      t_data[i] = ( t_data[i] < newval.t_data[i] ) ? newval.t_data[i] : t_data[i];
  }

  void
  SetMin( const TensorBaseT< T_length > & newval )
  {
    for( int i = 0; i < T_length; ++i )
      t_data[i] = ( t_data[i] > newval.t_data[i] ) ? newval.t_data[i] : t_data[i];
  }

  friend inline GEOSX_HOST_DEVICE realT
  Dot( const TensorBaseT< T_length > & A,
       const TensorBaseT< T_length > & B )
  {
    realT rval = 0;
    for( int i = 0; i < T_length; ++i )
    {
      rval += A.t_data[i] * B.t_data[i];
    }
    return rval;
  }

private:
  //  TensorBaseT(TensorBaseT<T_length>&);
};
//*****************************************************************************
//***** END DECLARATION *******************************************************
//*****************************************************************************

//*****************************************************************************
//***** TensorBaseT Member Function Definition ********************************
//*****************************************************************************

//**** CONSTRUCTORS AND DESTRUCTORS *******************************************

/**
 * @return none
 */
template< int T_length >
GEOSX_HOST_DEVICE
TensorBaseT< T_length >::TensorBaseT( void )  //:
{
  *this = 0.0;
}

/**
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template< int T_length >
TensorBaseT< T_length >::TensorBaseT( const realT data )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] = data;
}

/**
 * @param[in] data naked array used for initialization of t_data
 * @return none
 */
template< int T_length >
TensorBaseT< T_length >::TensorBaseT( const realT data[T_length] )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] = data[i];
}

//***** ASSIGNMENT OPERATORS **************************************************

/**
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template< int T_length >
GEOSX_FORCE_INLINE TensorBaseT< T_length > &
TensorBaseT< T_length >::operator=(
  const int & rhs )
{
  operator=( static_cast< realT >( rhs ) );
  return *this;
}

/**
 * @param[in] rhs value to set each member of t_data to
 * @return none
 */
template< int T_length >
GEOSX_FORCE_INLINE GEOSX_HOST_DEVICE TensorBaseT< T_length > &
TensorBaseT< T_length >::operator=( const realT & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] = rhs;
  return *this;
}

/**
 * @param[in] rhs value to add to t_data
 * @return none
 */
template< int T_length >
inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator+=( const realT & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] += rhs;
  return *this;
}

template<>
inline TensorBaseT< 3 > &
TensorBaseT< 3 >::operator+=( const realT & rhs )
{
  t_data[0] += rhs;
  t_data[1] += rhs;
  t_data[2] += rhs;
  return *this;
}

/**
 * @param[in] rhs value to subtract from t_data
 * @return none
 */
template< int T_length >
inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator-=( const realT & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] -= rhs;
  return *this;
}

/**
 * @param[in] rhs value to multiply t_data with
 * @return none
 */
template< int T_length >
GEOSX_HOST_DEVICE GEOSX_FORCE_INLINE TensorBaseT< T_length > &
TensorBaseT< T_length >::operator*=( const realT & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] *= rhs;
  return *this;
}

template<>
inline TensorBaseT< 3 > &
TensorBaseT< 3 >::operator*=( const realT & rhs )
{
  t_data[0] *= rhs;
  t_data[1] *= rhs;
  t_data[2] *= rhs;
  return *this;
}

/**
 * @param[in] rhs value to divide t_data with
 * @return none
 */
template< int T_length >
inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator/=( const realT & rhs )
{
  const realT irhs = 1 / rhs;
  operator*=( irhs );
  return *this;
}

/**
 * @param[in] rhs tensor to add
 * @return none
 */
template< int T_length >
GEOSX_FORCE_INLINE TensorBaseT< T_length > &
TensorBaseT< T_length >::operator+=(
  const TensorBaseT< T_length > & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] += rhs.t_data[i];
  return *this;
}

template<>
GEOSX_HOST_DEVICE inline TensorBaseT< 3 > &
TensorBaseT< 3 >::operator+=(
  const TensorBaseT< 3 > & rhs )
{
  t_data[0] += rhs.t_data[0];
  t_data[1] += rhs.t_data[1];
  t_data[2] += rhs.t_data[2];
  return *this;
}

/**
 * @param[in] rhs tensor to subract
 * @return none
 */
template< int T_length >
GEOSX_HOST_DEVICE inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator-=(
  const TensorBaseT< T_length > & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] -= rhs.t_data[i];
  return *this;
}

/**
 * @param[in] rhs tensor to multiply by
 * @return none
 */
template< int T_length >
inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator*=(
  const TensorBaseT< T_length > & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] *= rhs.t_data[i];
  return *this;
}

/**
 * @param[in] rhs tensor to divide by
 * @return none
 */
template< int T_length >
inline TensorBaseT< T_length > &
TensorBaseT< T_length >::operator/=(
  const TensorBaseT< T_length > & rhs )
{
  for( int i = 0; i < T_length; ++i ) t_data[i] /= rhs.t_data[i];
  return *this;
}

/*
   template<int N>
   TensorBaseT<N>&& operator+(TensorBaseT<N> &&src1, const TensorBaseT<N> &src2)
   {
   for (int i=0; i<N; ++i)
    src1[i] += src2[i];

   return src1;
   }*/

#endif
