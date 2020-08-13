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
 * @file DefaultValue.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_DEFAULTVALUE_HPP_
#define GEOSX_DATAREPOSITORY_DEFAULTVALUE_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "codingUtilities/traits.hpp"

namespace geosx
{
namespace dataRepository
{
namespace internal
{
/**
 * @struct is_defaultable
 * @tparam T type to check
 * @brief trait to determine if type @p T should have a default value
 */
template< typename T >
struct is_defaultable
{
  /// attribute to set what type is able to contain a default value
  static constexpr bool value = std::is_arithmetic< T >::value ||
    std::is_same< T, string >::value || std::is_same< T, Path >::value ||
    traits::is_tensorT< T > || std::is_enum< T >::value;
};

/**
 * @struct Helper
 * @tparam T type to check
 * @tparam ENABLE template parameter for use with SFINAE
 *
 * default implementation of struct to return if a type \p T has a default value.
 */
template< typename T, typename ENABLE = void >
struct Helper
{
  /// attribute to indicate whether type \p T has a default value
  static constexpr bool has_default_value = false;

  /// alias for default value type (void be default)
  using value_type = void;
};

/**
 * @struct Helper
 * @tparam T type to check
 *
 * Specialization of Helper struct to return if a type @p T has a default
 * value. This specialization specifically tests the type itself. Contains
 * a member to hold a default value.
 */
template< typename T >
struct Helper< T, std::enable_if_t< is_defaultable< T >::value > >
{
  /// attribute to indicate whether type @p T has a default value
  static constexpr bool has_default_value = true;

  /// alias for the type T
  using value_type = T;

  /// a member to hold a default value for the type @p T
  value_type value = value_type();
};

/**
 * @struct Helper
 * @tparam T type to check
 *
 * Specialization of Helper struct to return if a type @p T has a default
 * value. This specialization specifically tests the type has an alias
 * named "value_type" as is the case for stl containers and GEOSX
 * containers.
 */
template< typename T >
struct Helper< T,
               std::enable_if_t< traits::HasAlias_value_type< T > &&
                                 is_defaultable< typename T::value_type >::value &&
                                 !traits::is_string< T > > >
{
  /// attribute to indicate whether type @p T has a default value
  static constexpr bool has_default_value = true;

  /// alias for the type T
  using value_type = typename T::value_type;

  /// a member to hold a default value for the type @p T
  value_type value = value_type();
};

template< typename T >
std::enable_if_t< !Helper< T >::has_default_value, std::ostream & >
operator<<(
  std::ostream & stream,
  Helper< T > const & GEOSX_UNUSED_PARAM( value ) )
{
  return stream;
}

template< typename T >
std::enable_if_t< Helper< T >::has_default_value, std::ostream & >
operator<<(
  std::ostream & stream,
  Helper< T > const & value )
{
  return stream << value.value;
}

}  // namespace internal

/**
 * @tparam T the type to check
 * @brief A templated alias to hold default values.
 */
template< typename T >
using DefaultValue = internal::Helper< T >;

}  // namespace dataRepository
}  // namespace geosx

#endif /* GEOSX_DATAREPOSITORY_DEFAULTVALUE_HPP_ */
