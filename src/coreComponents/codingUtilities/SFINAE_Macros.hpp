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
 * @file sfinae.hpp
 */

#ifndef SRC_CODINGUTILITIES_SFINAE_HPP_
#define SRC_CODINGUTILITIES_SFINAE_HPP_

#include "LvArray/src/Macros.hpp"
#include "LvArray/src/typeManipulation.hpp"

/**
 * @brief Macro that expands to a static constexpr bool templated on two types that is only true when
 *        the first type has a method @p NAME which takes arguments @p __VA_ARGS__ which may contain the second
 *        template type. The name of the boolean variable is HasMemberFunction_ ## @p NAME.
 * @param NAME The name of the method to look for.
 * @param T The name of the template type in the argument list.
 * @param __VA_ARGS__ The argument list to call the method with.
 * @note The class type is available through the name CLASS.
 * @note This doesn't check that the templated class method instantiation is valid, only that a
 *       matching method was found.
 */
#define HAS_MEMBER_FUNCTION_TEMPLATE_NO_RTYPE( NAME, T, ... ) \
  IS_VALID_EXPRESSION_2( HasMemberFunction_##NAME,            \
                         CLASS,                               \
                         T,                                   \
                         std::declval< CLASS >().NAME( __VA_ARGS__ ) )

/**
 * @brief Macro that expands to a static constexpr bool templated on a type that is only true when
 *        the type has a method @p NAME which takes arguments @p __VA_ARGS__ and whose return value is convertible
 *        to @p RTYPE. The name of the boolean variable is HasMemberFunction_ ## @p NAME.
 * @param NAME The name of the method to look for.
 * @param RTYPE The type to attempt to convert the return value to.
 * @param __VA_ARGS__ The argument list to call the method with.
 * @note The class type is available through the name CLASS.
 */
#define HAS_MEMBER_FUNCTION( NAME, RTYPE, ... )                                   \
  IS_VALID_EXPRESSION(                                                            \
    HasMemberFunction_##NAME,                                                     \
    CLASS,                                                                        \
    std::is_convertible< decltype( std::declval< CLASS >().NAME( __VA_ARGS__ ) ), \
                         RTYPE >::value )

/**
 * @brief Macro that expands to a static constexpr bool templated on a type that is only true when
 *        the type has a an alias @p NAME. The name of the boolean variable is HasAlias_ ## @p NAME.
 */
#define HAS_ALIAS( NAME )               \
  IS_VALID_EXPRESSION( HasAlias_##NAME, \
                       CLASS,           \
                       !std::is_enum< typename CLASS::NAME >::value )

#endif /* SRC_CODINGUTILITIES_SFINAE_HPP_ */
