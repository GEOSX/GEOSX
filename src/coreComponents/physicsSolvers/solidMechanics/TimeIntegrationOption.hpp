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

#pragma once

/// System includes
#include <ostream>
#include <istream>

/**
 * @file TimeIntegrationOption.hpp
 */

namespace geosx
{
/**
 * @enum TimeIntegrationOption
 *
 * The options for time integration
 */
enum class TimeIntegrationOption : int
{
  QuasiStatic,      //!< QuasiStatic
  ImplicitDynamic,  //!< ImplicitDynamic
  ExplicitDynamic   //!< ExplicitDynamic
};

/**
 * @brief Function to get a TimeIntegrationOption enum from an int
 * @param val int that represents the TimeIntegrationOption
 * @return The TimeIntegrationOption that corresponds to the input
 */
inline TimeIntegrationOption
toTimeIntegrationOption( std::string const & val )
{
  if( val == "TimeIntegrationOption::QuasiStatic" || val == "QuasiStatic" )
  {
    return TimeIntegrationOption::QuasiStatic;
  }
  if( val == "TimeIntegrationOption::ImplicitDynamic" || val == "ImplicitDynamic" )
  {
    return TimeIntegrationOption::ImplicitDynamic;
  }
  if( val == "TimeIntegrationOption::ExplicitDynamic" || val == "ExplicitDynamic" )
  {
    return TimeIntegrationOption::ExplicitDynamic;
  }
  else
  {
    GEOSX_ERROR( "Could not parse " << val << " into a TimeIntegrationOption." );
    return TimeIntegrationOption::QuasiStatic;
  }
}

inline std::istream &
operator>>( std::istream & is, TimeIntegrationOption & option )
{
  std::string value;
  is >> value;
  option = toTimeIntegrationOption( value );
  return is;
}

inline std::ostream &
operator<<( std::ostream & os,
            TimeIntegrationOption const & option )
{
  switch( option )
  {
    case TimeIntegrationOption::QuasiStatic:
    {
      os << "TimeIntegrationOption::QuasiStatic";
      return os;
    }
    case TimeIntegrationOption::ImplicitDynamic:
    {
      os << "TimeIntegrationOption::ImplicitDynamic";
      return os;
    }
    case TimeIntegrationOption::ExplicitDynamic:
    {
      os << "TimeIntegrationOption::ExplicitDynamic";
      return os;
    }
    default:
    {
      GEOSX_ERROR( "TimeIntegrationOption " << static_cast< int >( option )
                                            << " not recognized." );
      return os;
    }
  }
}

}  // namespace geosx
