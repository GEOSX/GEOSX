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
 * @file ExponentialRelations.hpp
 */

#ifndef GEOSX_CONSITUTIVE_EXPONENTIALRELATION_HPP_
#define GEOSX_CONSITUTIVE_EXPONENTIALRELATION_HPP_

#include <common/DataTypes.hpp>

#include <cmath>

namespace geosx
{
namespace constitutive
{
/**
 * @enum Enumeration describing available approximations of the exponent
 */
enum class ExponentApproximationType
{
  Full,
  Linear,
  Quadratic
};

namespace detail
{
// bring into scope for ADL
using std::exp;
using std::log;
using std::sqrt;

/**
 * @brief Specializations of this struct provide compute methods for different exponent approximations
 * @tparam T the data type
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full)
 * @tparam INVERSE whether to compute the inverse of the exponent
 */
template< typename T, ExponentApproximationType EAT, bool INVERSE = false >
struct ExponentialCompute
{};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Full, false >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * exp( alpha * ( x - x0 ) );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & x,
           T & y,
           T & dy_dx )
  {
    y = y0 * exp( alpha * ( x - x0 ) );
    dy_dx = alpha * y;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Full, true >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + log( y / y0 ) / alpha;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & y,
           T & x,
           T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * log( y / y0 );
    dx_dy = alpha_inv / y;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Quadratic, false >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    const T z = T( 1.0 ) + alpha * ( x - x0 );
    y = y0 / T( 2.0 ) * ( T( 1.0 ) + z * z );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & x,
           T & y,
           T & dy_dx )
  {
    const T z = T( 1.0 ) + alpha * ( x - x0 );
    y = y0 / T( 2.0 ) * ( T( 1.0 ) + z * z );
    dy_dx = alpha * y0 * z;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Quadratic, true >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ) );
    x = x0 + ( z - 1 ) / alpha;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & y,
           T & x,
           T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    const T z = sqrt( T( 2.0 ) * y / y0 - T( 1.0 ) );
    x = x0 + alpha_inv * ( z - 1 );
    dx_dy = alpha_inv / ( y0 * z );
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Linear, false >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & x, T & y )
  {
    y = y0 * ( T( 1.0 ) + alpha * ( x - x0 ) );
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & x,
           T & y,
           T & dy_dx )
  {
    y = y0 * ( T( 1.0 ) + alpha * ( x - x0 ) );
    dy_dx = alpha * y0;
  }
};

template< typename T >
struct ExponentialCompute< T, ExponentApproximationType::Linear, true >
{
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0, const T & y0, const T & alpha, const T & y, T & x )
  {
    x = x0 + ( y / y0 - 1 ) / alpha;
  }

  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  static void
  Compute( const T & x0,
           const T & y0,
           const T & alpha,
           const T & y,
           T & x,
           T & dx_dy )
  {
    const T alpha_inv = T( 1.0 ) / alpha;
    x = x0 + alpha_inv * ( y / y0 - 1 );
    dx_dy = alpha_inv / y0;
  }
};

}  // namespace detail

/**
 * @class ExponentialRelation
 *
 * Describes a simple exponential relationship with three parameters of the form: \f$ y = y0 * exp(a * (x - x0)) \f$.
 * Optional linear and quadratic Taylor series approximation of the exponent.
 * Implements point and batch direct (y(x)) and inverse (x(y)) compute functions and provides derivatives.
 *
 * @tparam T scalar real-valued type used in computation
 * @tparam EAT the type/order of exponent approximation (linear, quadratic, full)
 */
template< typename T, ExponentApproximationType EAT >
class ExponentialRelation
{
public:
  static constexpr ExponentApproximationType type = EAT;

  // *** constructors ***

  /**
   * @brief Default constructor. Sets \f$ x0 = 0, y0 = 1, alpha = 1 \f$
   */
  GEOSX_HOST_DEVICE
  ExponentialRelation() :
    ExponentialRelation( T( 0 ), T( 1 ), T( 1 ) )
  {}

  /**
   * @brief Constructor with arguments
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  GEOSX_HOST_DEVICE
  ExponentialRelation( T x0, T y0, T alpha )
  {
    SetCoefficients( x0, y0, alpha );
  }

  // *** setters ***

  /**
   * @brief Sets the relation coefficients
   * @param x0 variable shift
   * @param y0 scaling coefficient
   * @param alpha exponential coefficient
   */
  GEOSX_HOST_DEVICE
  void
  SetCoefficients( T x0, T y0, T alpha )
  {
    m_x0 = x0;
    m_y0 = y0;
    m_alpha = alpha;
  }

  // *** no-derivative computes ***

  /**
   * @brief Compute value
   * @param x
   * @param y
   */
  GEOSX_HOST_DEVICE
  void
  Compute( const T & x, T & y ) const
  {
    detail::ExponentialCompute< T, EAT >::Compute( m_x0, m_y0, m_alpha, x, y );
  }

  /**
   * @brief Compute inverse value
   * @param y
   * @param x
   */
  GEOSX_HOST_DEVICE
  void
  Inverse( const T & y, T & x ) const
  {
    detail::ExponentialCompute< T, EAT, true >::Compute( m_x0, m_y0, m_alpha, y, x );
  }

  // *** derivative computes ***

  /**
   * @brief Compute value and derivative
   * @param x
   * @param y
   * @param dy_dx
   */
  GEOSX_HOST_DEVICE
  void
  Compute( const T & x, T & y, T & dy_dx ) const
  {
    detail::ExponentialCompute< T, EAT >::Compute( m_x0, m_y0, m_alpha, x, y, dy_dx );
  }

  /**
   * @brief Compute inverse value and derivative
   * @param y
   * @param x
   * @param dx_dy
   */
  GEOSX_HOST_DEVICE
  void
  Inverse( const T & y, T & x, T & dx_dy ) const
  {
    detail::ExponentialCompute< T, EAT, true >::Compute( m_x0, m_y0, m_alpha, y, x, dx_dy );
  }

private:
  /// variable shift
  T m_x0;

  /// scaling coefficient
  T m_y0;

  /// exponential coefficient
  T m_alpha;
};

inline ExponentApproximationType
stringToExponentType( string const & model )
{
  static std::map< string, ExponentApproximationType > const approxTypes = {
    { "linear", ExponentApproximationType::Linear },
    { "quadratic", ExponentApproximationType::Quadratic },
    { "exponential", ExponentApproximationType::Full },
  };
  auto const it = approxTypes.find( model );
  GEOSX_ERROR_IF( it == approxTypes.end(), "Model type not supported: " << model );
  return it->second;
}

template< ExponentApproximationType EAT >
struct ExponentApproximationTypeWrapper
{
  static constexpr ExponentApproximationType value = EAT;
};

/**
 * @brief Function that implements the switchyard for ExponentApproximationType
 * @tparam LAMBDA the type of the user-provided generic lambda
 * @param type type/order of the approximation (linear, quadratic, full)
 * @param lambda user-provided generic lambda to be called with an instance of EAT wrapper type
 */
template< typename T, typename LAMBDA >
void
ExponentApproximationTypeSwitchBlock( ExponentApproximationType const type,
                                      T const & x0,
                                      T const & y0,
                                      T const & alpha,
                                      LAMBDA lambda )
{
  switch( type )
  {
    case ExponentApproximationType::Full:
    {
      return lambda(
        ExponentialRelation< T, ExponentApproximationType::Full >( x0, y0, alpha ) );
    }
    case ExponentApproximationType::Quadratic:
    {
      return lambda(
        ExponentialRelation< T, ExponentApproximationType::Quadratic >( x0, y0, alpha ) );
    }
    case ExponentApproximationType::Linear:
    {
      return lambda(
        ExponentialRelation< T, ExponentApproximationType::Linear >( x0, y0, alpha ) );
    }
    default:
    {
      GEOSX_ERROR(
        "ExponentApproximationTypeSwitchBlock() ExponentApproximationType is "
        "invalid!" );
    }
  }
}

/**
 * @brief Same as above, but with non-default constructor arguments (relation coefficients)
 * @tparam T scalar real-valued type
 * @tparam LAMBDA type of the user-provided generic lambda
 * @param type runtime choice of ExponentApproximationType
 * @param x0 variable shift
 * @param y0 scaling coefficient
 * @param alpha exponential coefficient
 * @param lambda user-provided generic lambda that will be passed an instance of ExponentialRelation<T,?>
 */
template< typename T, typename LAMBDA >
void
makeExponentialRelation( ExponentApproximationType type,
                         T const & x0,
                         T const & y0,
                         T const & alpha,
                         LAMBDA && lambda )
{
  ExponentApproximationTypeSwitchBlock( type,
                                        x0,
                                        y0,
                                        alpha,
                                        std::forward< LAMBDA >( lambda ) );
}

}  // namespace constitutive

}  // namespace geosx

#endif  // GEOSX_CONSITUTIVE_EXPONENTIALRELATION_HPP_
