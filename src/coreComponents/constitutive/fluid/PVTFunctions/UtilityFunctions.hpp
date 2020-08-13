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
 * @file TableFunction.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_UTILITYFUNCTION_HPP_
#define GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_UTILITYFUNCTION_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"

namespace geosx
{
namespace PVTProps
{
constexpr localIndex MAX_VAR_DIM = 10;

template< typename T, int Dim >
class EvalArgs
{
public:
  EvalArgs()
  {
    m_var = 0.0;
    for( int i = 0; i < Dim; ++i )
    {
      m_der[i] = 0;
    }
  }

  ~EvalArgs() = default;

  EvalArgs( const EvalArgs & arg )
  {
    m_var = arg.m_var;
    for( int i = 0; i < Dim; ++i )
    {
      m_der[i] = arg.m_der[i];
    }
  }

  EvalArgs( const T & var )
  {
    m_var = var;
    for( int i = 0; i < Dim; ++i )
    {
      m_der[i] = 0;
    }
  }

  EvalArgs &
  operator+=( const EvalArgs & arg )
  {
    this->m_var += arg.m_var;
    for( localIndex i = 0; i < Dim; ++i )
    {
      this->m_der[i] += arg.m_der[i];
    }

    return *this;
  }

  EvalArgs &
  operator-=( const EvalArgs & arg )
  {
    this->m_var -= arg.m_var;
    for( localIndex i = 0; i < Dim; ++i )
    {
      this->m_der[i] -= arg.m_der[i];
    }

    return *this;
  }

  EvalArgs &
  operator*=( const EvalArgs & arg )
  {
    const T & u = this->m_var;
    const T & v = arg.m_var;
    for( localIndex i = 0; i < Dim; ++i )
    {
      const T & uDer = this->m_der[i];
      const T & vDer = arg.m_der[i];

      this->m_der[i] = ( v * uDer + u * vDer );
    }

    this->m_var *= v;
    return *this;
  }

  EvalArgs &
  operator/=( const EvalArgs & arg )
  {
    const T & u = this->m_var;
    const T & v = arg.m_var;
    for( localIndex i = 0; i < Dim; ++i )
    {
      const T & uDer = this->m_der[i];
      const T & vDer = arg.m_der[i];

      this->m_der[i] = ( v * uDer - u * vDer ) / ( v * v );
    }

    this->m_var /= v;
    return *this;
  }

  EvalArgs &
  Exponent( const EvalArgs & arg )
  {
    const T & v = arg.m_var;
    this->m_var = exp( v );
    for( localIndex i = 0; i < Dim; ++i )
    {
      const T & vDer = arg.m_der[i];
      this->m_der[i] = this->m_var * vDer;
    }

    return *this;
  }

  EvalArgs
  operator+( const EvalArgs & arg ) const
  {
    EvalArgs result( *this );
    result += arg;
    return result;
  }

  EvalArgs
  operator-( const EvalArgs & arg ) const
  {
    EvalArgs result( *this );
    result -= arg;
    return result;
  }

  EvalArgs
  operator-() const
  {
    EvalArgs result;
    result.m_var = -this->m_var;
    for( localIndex i = 0; i < Dim; ++i )
    {
      result.m_der[i] = -this->m_der[i];
    }

    return result;
  }

  EvalArgs
  operator*( const EvalArgs & arg ) const
  {
    EvalArgs result( *this );
    result *= arg;
    return result;
  }

  EvalArgs
  operator/( const EvalArgs & arg ) const
  {
    EvalArgs result( *this );
    result /= arg;
    return result;
  }

  EvalArgs &
  operator=( const EvalArgs & arg )
  {
    this->m_var = arg.m_var;
    for( int i = 0; i < Dim; ++i )
    {
      m_der[i] = arg.m_der[i];
    }

    return *this;
  }

  bool
  operator==( const EvalArgs & arg ) const
  {
    if( this->m_var != arg._m_var ) return false;

    for( localIndex i = 0; i < Dim; ++i )
    {
      if( this->m_der[i] != arg.m_der[i] ) return false;
    }

    return true;
  }

  bool
  operator!=( const EvalArgs & arg ) const
  {
    return !( *this == arg );
  }

  bool
  operator>( const EvalArgs & arg ) const
  {
    return this->m_var > arg.m_var;
  }

  bool
  operator<( const EvalArgs & arg ) const
  {
    return this->m_var < arg.m_var;
  }

  bool
  operator>=( const EvalArgs & arg ) const
  {
    return this->m_var >= arg.m_var;
  }

  bool
  operator<=( const EvalArgs & arg ) const
  {
    return this->m_var <= arg.m_var;
  }

  EvalArgs &
  operator+=( const T & arg )
  {
    this->m_var += arg;
    return *this;
  }

  EvalArgs &
  operator-=( const T & arg )
  {
    this->m_var -= arg;
    return *this;
  }

  EvalArgs &
  operator*=( const T & arg )
  {
    for( localIndex i = 0; i < Dim; ++i )
    {
      this->m_der[i] *= arg;
    }

    this->m_var *= arg;
    return *this;
  }

  EvalArgs &
  operator/=( const T & arg )
  {
    for( localIndex i = 0; i < Dim; ++i )
    {
      this->m_der[i] /= arg;
    }

    this->m_var /= arg;
    return *this;
  }

  EvalArgs
  operator+( const T & arg ) const
  {
    EvalArgs result( *this );
    result += arg;
    return result;
  }

  EvalArgs
  operator-( const T & arg ) const
  {
    EvalArgs result( *this );
    result -= arg;
    return result;
  }

  EvalArgs
  operator*( const T & arg ) const
  {
    EvalArgs result( *this );
    result *= arg;
    return result;
  }

  EvalArgs
  operator/( const T & arg ) const
  {
    EvalArgs result( *this );
    result /= arg;
    return result;
  }

  EvalArgs &
  operator=( const T & arg )
  {
    m_var = arg;
    for( int i = 0; i < Dim; ++i )
    {
      m_der[i] = 0;
    }
    return *this;
  }

  bool
  operator==( const T & arg ) const
  {
    return this->m_var == arg;
  }

  bool
  operator!=( const T & arg ) const
  {
    return !( *this == arg );
  }

  bool
  operator>( const T & arg ) const
  {
    return this->m_var > arg;
  }

  bool
  operator<( const T & arg ) const
  {
    return this->m_var < arg;
  }

  bool
  operator>=( const T & arg ) const
  {
    return this->m_var >= arg;
  }
  bool
  operator<=( const T & arg ) const
  {
    return this->m_var <= arg;
  }

  T m_var;
  T m_der[Dim];
};

template< class T, int Dim >
inline EvalArgs< T, Dim >
operator+( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  EvalArgs< T, Dim > result( arg2 );
  result += arg1;
  return result;
}

template< class T, int Dim >
inline EvalArgs< T, Dim >
operator-( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  EvalArgs< T, Dim > result( arg2 );
  result.m_var = arg1 - result.m_var;

  for( localIndex i = 0; i < Dim; ++i )
  {
    result.m_der[i] = -result.m_der[i];
  }

  return result;
}

template< class T, int Dim >
inline EvalArgs< T, Dim >
operator*( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  EvalArgs< T, Dim > result( arg2 );
  result *= arg1;
  return result;
}

template< class T, int Dim >
inline EvalArgs< T, Dim >
operator/( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  EvalArgs< T, Dim > result( arg2 );

  T coef = -( arg1 / result.m_var / result.m_var );
  result.m_var = arg1 / result.m_var;

  for( localIndex i = 0; i < Dim; ++i )
  {
    result.m_der[i] = coef * result.m_der[i];
  }

  return result;
}

template< class T, int Dim >
inline bool
operator==( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return ( arg2 == arg1 );
}

template< class T, int Dim >
inline bool
operator!=( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return ( arg2 != arg1 );
}

template< class T, int Dim >
inline bool
operator>( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return arg1 > arg2.m_var;
}

template< class T, int Dim >
inline bool
operator<( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return arg1 < arg2.m_var;
}

template< class T, int Dim >
inline bool
operator>=( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return arg1 >= arg2.m_var;
}

template< class T, int Dim >
inline bool
operator<=( const T & arg1, const EvalArgs< T, Dim > & arg2 )
{
  return arg1 <= arg2.m_var;
}

typedef EvalArgs< real64, 1 > EvalArgs1D;
typedef EvalArgs< real64, 2 > EvalArgs2D;
typedef EvalArgs< real64, 3 > EvalArgs3D;
typedef EvalArgs< real64, MAX_VAR_DIM > EvalVarArgs;

class TableFunctionBase
{
public:
  virtual const string &
  TableName() const = 0;
  virtual ~TableFunctionBase()
  {}

  virtual EvalArgs1D
  Value( const EvalArgs1D & x ) const = 0;

  virtual EvalArgs2D
  Value( const EvalArgs2D & x ) const = 0;

  virtual EvalArgs2D
  Value( const EvalArgs2D & x, const EvalArgs2D & y ) const = 0;

  virtual void
  Print() const = 0;
};

typedef std::shared_ptr< TableFunctionBase > TableFunctionPtr;

class XYTable : public TableFunctionBase
{
public:
  XYTable( std::string const & tableName,
           real64_array const & x,
           real64_array const & y,
           real64_array2d const & value ) :
    m_tableName( tableName ),
    m_x( x ),
    m_y( y ),
    m_value( value )
  {}

  ~XYTable()
  {}

  real64_array &
  XArray()
  {
    return m_x;
  }

  real64_array &
  YArray()
  {
    return m_y;
  }

  real64_array2d &
  ValueArray()
  {
    return m_value;
  }

  virtual string const &
  TableName() const
  {
    return m_tableName;
  }

  virtual EvalArgs1D
  Value( EvalArgs1D const & ) const
  {
    return 0;
  }

  virtual EvalArgs2D
  Value( EvalArgs2D const & ) const
  {
    return 0;
  }

  virtual EvalArgs2D
  Value( EvalArgs2D const & x, EvalArgs2D const & y ) const;

  virtual void
  Print() const
  {}

private:
  std::string m_tableName;
  real64_array m_x;
  real64_array m_y;
  real64_array2d m_value;
};

class XTable : public TableFunctionBase
{
public:
  XTable( string const & tableName, real64_array const & x, real64_array const & value ) :
    m_tableName( tableName ),
    m_x( x ),
    m_value( value )
  {}
  ~XTable()
  {}

  real64_array &
  XArray()
  {
    return m_x;
  }

  real64_array &
  ValueArray()
  {
    return m_value;
  }

  virtual string const &
  TableName() const
  {
    return m_tableName;
  }

  virtual EvalArgs1D
  Value( EvalArgs1D const & x ) const
  {
    return GetValue< EvalArgs1D >( x );
  }

  virtual EvalArgs2D
  Value( EvalArgs2D const & x ) const
  {
    return GetValue< EvalArgs2D >( x );
  }

  virtual EvalArgs2D
  Value( EvalArgs2D const &, EvalArgs2D const & ) const
  {
    return 0;
  }

private:
  template< class T >
  T
  GetValue( T const & x ) const;

  virtual void
  Print() const
  {}

  string m_tableName;
  real64_array m_x;
  real64_array m_value;
};

}  // namespace PVTProps
}  // namespace geosx

#endif  //GEOSX_CONSTITUTIVE_FLUID_PVTFUNCTIONS_UTILITYFUNCTION_HPP_
