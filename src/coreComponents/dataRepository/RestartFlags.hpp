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

#ifndef GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_
#define GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_

/**
 * @file RestartFlags.hpp
 */

namespace geosx
{
namespace dataRepository
{
/**
 * @enum RestartFlags
 *
 * A scoped enum for the restart options.
 */
enum class RestartFlags : int
{
  NO_WRITE,       ///< Doe no write into restart
  WRITE,          ///< Write into restart
  WRITE_AND_READ  ///< Write and read from restart
};

/**
 * @enum PlotLevel
 *
 * A scoped enum for the Plot options.
 */
enum class PlotLevel : int
{
  LEVEL_0,  ///< Write to plot always
  LEVEL_1,  ///< Write to plot when plotLevel>=1 is specified in input
  LEVEL_2,  ///< Write to plot when plotLevel>=2 is specified in input
  LEVEL_3,  ///< Write to plot when plotLevel>=3 is specified in input
  NOPLOT    ///< Do not ever write to plot file
};

/**
 * @brief Function to get a PlotLevel enum from an int
 * @param val int that represents the PlotLevel
 * @return The PlotLevel that corresponds to the input
 */
inline PlotLevel
toPlotLevel( int const val )
{
  switch( val )
  {
    case static_cast< int >( PlotLevel::LEVEL_0 ):
    {
      return PlotLevel::LEVEL_0;
    }
    case static_cast< int >( PlotLevel::LEVEL_1 ):
    {
      return PlotLevel::LEVEL_1;
    }
    case static_cast< int >( PlotLevel::LEVEL_2 ):
    {
      return PlotLevel::LEVEL_2;
    }
    case static_cast< int >( PlotLevel::LEVEL_3 ):
    {
      return PlotLevel::LEVEL_3;
    }
    case static_cast< int >( PlotLevel::NOPLOT ):
    {
      return PlotLevel::NOPLOT;
    }
    default:
    {
      GEOSX_ERROR( "Could not parse " << val << " into a PlotLevel." );
      return PlotLevel::NOPLOT;
    }
  }
}

/**
 * @brief Reads a PlotLevel enum from a stream.
 * @param is The stream to read from.
 * @param plotLevel The PlotLevel to write to.
 * @return The stream.
 */
inline std::istream &
operator>>( std::istream & is, PlotLevel & plotLevel )
{
  int value;
  is >> value;
  plotLevel = toPlotLevel( value );
  return is;
}

/**
 * @brief Writes a plot level to a stream.
 * @param os The stream to read from.
 * @param plotLevel the PlotLevel to write.
 * @return The stream.
 */
inline std::ostream &
operator<<( std::ostream & os, PlotLevel const & plotLevel )
{
  return os << static_cast< int >( plotLevel );
}

}  // namespace dataRepository
}  // namespace geosx

#endif /* GEOSX_DATAREPOSITORY_RESTARTFLAGS_HPP_ */
