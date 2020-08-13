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
 * @file ConstitutivePassThruHandler.hpp
 */
#ifndef GEOSX_CONSTITUTIVEPASSTHRUHANDLER_HPP
#define GEOSX_CONSTITUTIVEPASSTHRUHANDLER_HPP

#include "codingUtilities/traits.hpp"
#include "common/DataTypes.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @brief Struct to help in dynamic dispatch based on type of contitutive model.
 * @tparam TYPES pack of derived constitutive types to handle
 */
template< typename... TYPES >
struct ConstitutivePassThruHandler
{};

/**
 * @brief Specialization for an empty type pack (always errors)
 */
template<>
struct ConstitutivePassThruHandler<>
{
  template< typename BASE, typename LAMBDA >
  static void
  Execute( BASE & relation, LAMBDA lambda )
  {
    GEOSX_UNUSED_VAR( relation )
    GEOSX_UNUSED_VAR( lambda )
    GEOSX_ERROR(
      "The constitutive model "
      << relation.getName() << " was not dispatched."
      << "The model type does not match the list of supported types." );
  }
};

/**
 * @brief Specialization for a non-empty type pack
 * @tparam TYPE first type to handle
 * @tparam TYPES rest of the type list
 */
template< typename TYPE, typename... TYPES >
struct ConstitutivePassThruHandler< TYPE, TYPES... >
{
  template< typename BASE, typename LAMBDA >
  static void
  Execute( BASE & relation, LAMBDA && lambda )
  {
    using Derived = add_const_if_t< TYPE, std::is_const< BASE >::value >;

    if( dynamicCast< Derived * >( &relation ) )
    {
      lambda( static_cast< Derived & >( relation ) );
    }
    else
    {
      ConstitutivePassThruHandler< TYPES... >::Execute(
        relation,
        std::forward< LAMBDA >( lambda ) );
    }
  }
};

}  //namespace constitutive

}  //namespace geosx

#endif  //GEOSX_CONSTITUTIVEPASSTHRUHANDLER_HPP
