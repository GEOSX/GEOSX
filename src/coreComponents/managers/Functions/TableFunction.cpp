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
 * @file TableFunction.cpp
 */

#include "TableFunction.hpp"
#include "common/DataTypes.hpp"
#include <algorithm>

namespace geosx
{
namespace dataRepository
{
namespace keys
{
std::string const tableCoordinates = "coordinates";
std::string const tableValues = "values";
std::string const tableInterpolation = "interpolation";
std::string const coordinateFiles = "coordinateFiles";
std::string const voxelFile = "voxelFile";
std::string const valueType = "valueType";
}  // namespace keys
}  // namespace dataRepository

using namespace dataRepository;

TableFunction::TableFunction( const std::string & name, Group * const parent ) :
  FunctionBase( name, parent ),
  m_tableCoordinates1D(),
  m_coordinateFiles(),
  m_voxelFile(),
  m_interpolationMethodString(),
  m_interpolationMethod( InterpolationType::Linear ),
  m_coordinates(),
  m_values(),
  m_dimensions( 0 ),
  m_size(),
  m_indexIncrement(),
  m_corners(),
  m_numCorners( 0 )
{
  registerWrapper( keys::tableCoordinates, &m_tableCoordinates1D )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Coordinates inputs for 1D tables" );

  registerWrapper( keys::tableValues, &m_values )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Values for 1D tables" );

  registerWrapper( keys::coordinateFiles, &m_coordinateFiles )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "List of coordinate file names for ND Table" );

  registerWrapper( keys::voxelFile, &m_voxelFile )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription( "Voxel file name for ND Table" );

  registerWrapper( keys::tableInterpolation, &m_interpolationMethodString )
    ->setInputFlag( InputFlags::OPTIONAL )
    ->setDescription(
      "Interpolation method (options = linear, nearest, upper, lower)" )
    ->setApplyDefaultValue( "linear" );
}

TableFunction::~TableFunction()
{}

template< typename T >
void
TableFunction::parse_file( array1d< T > & target,
                           string const & filename,
                           char delimiter )
{
  std::ifstream inputStream( filename.c_str() );
  std::string lineString;
  T value;

  GEOSX_ERROR_IF( !inputStream, "Could not read input file: " << filename );

  // Read the file
  // TODO: Update this to handle large parallel jobs
  while( !inputStream.eof() )
  {
    std::getline( inputStream, lineString );
    std::istringstream ss( lineString );

    while( ss.peek() == delimiter || ss.peek() == ' ' )
    {
      ss.ignore();
    }
    while( ss >> value )
    {
      target.emplace_back( value );
      while( ss.peek() == delimiter || ss.peek() == ' ' )
      {
        ss.ignore();
      }
    }
  }

  inputStream.close();
}

void
TableFunction::InitializeFunction()
{
  // Read in data
  if( m_coordinateFiles.empty() )
  {
    // 1D Table
    m_dimensions = 1;
    m_coordinates.emplace_back( m_tableCoordinates1D );
    m_size.emplace_back( m_tableCoordinates1D.size() );

    // Check to make sure that the table dimensions match
    GEOSX_ERROR_IF(
      m_size[0] != m_values.size(),
      "1D Table function coordinates and values must have the same length." );
  }
  else
  {
    // ND Table
    m_dimensions =
      LvArray::integerConversion< localIndex >( m_coordinateFiles.size() );
    m_coordinates.resize( m_dimensions );

    parse_file( m_values, m_voxelFile, ',' );
    for( localIndex ii = 0; ii < m_dimensions; ++ii )
    {
      parse_file( m_coordinates[ii], m_coordinateFiles[ii], ',' );
      m_size.emplace_back( m_coordinates[ii].size() );
    }
  }

  // Parse the interpolation method string
  if( m_interpolationMethodString == "linear" )
  {
    m_interpolationMethod = InterpolationType::Linear;
  }
  else if( m_interpolationMethodString == "nearest" )
  {
    m_interpolationMethod = InterpolationType::Nearest;
  }
  else if( m_interpolationMethodString == "upper" )
  {
    m_interpolationMethod = InterpolationType::Upper;
  }
  else if( m_interpolationMethodString == "lower" )
  {
    m_interpolationMethod = InterpolationType::Lower;
  }
  else
  {
    GEOSX_ERROR(
      "Unrecognized interpolation type: " << m_interpolationMethodString );
  }

  reInitializeFunction();
}

void
TableFunction::reInitializeFunction()
{
  // Setup index increment (assume data is in Fortran array order)
  localIndex increment = 1;
  m_indexIncrement.resize( m_dimensions );
  for( localIndex ii = 0; ii < m_dimensions; ++ii )
  {
    m_size[ii] = m_coordinates[ii].size();
    m_indexIncrement[ii] = increment;
    increment *= m_size[ii];
  }

  // Error checking
  GEOSX_ERROR_IF( increment != m_values.size(), "Table dimensions do not match!" );

  // Build a quick map to help with linear interpolation
  m_numCorners = static_cast< localIndex >( pow( 2, m_dimensions ) );
  for( localIndex ii = 0; ii < m_numCorners; ++ii )
  {
    for( localIndex jj = 0; jj < m_dimensions; ++jj )
    {
      m_corners[jj][ii] = int( ii / pow( 2, jj ) ) % 2;
    }
  }
}

real64
TableFunction::Evaluate( real64 const * const input ) const
{
  real64 result = 0.0;

  // Linear interpolation
  if( m_interpolationMethod == InterpolationType::Linear )
  {
    localIndex bounds[m_maxDimensions][2];
    real64 weights[m_maxDimensions][2];

    // Determine position, weights
    for( localIndex ii = 0; ii < m_dimensions; ++ii )
    {
      if( input[ii] <= m_coordinates[ii][0] )
      {
        // Coordinate is to the left of this axis
        bounds[ii][0] = 0;
        bounds[ii][1] = 0;
        weights[ii][0] = 0;
        weights[ii][1] = 1;
      }
      else if( input[ii] >= m_coordinates[ii][m_size[ii] - 1] )
      {
        // Coordinate is to the right of this axis
        bounds[ii][0] = m_size[ii] - 1;
        bounds[ii][1] = bounds[ii][0];
        weights[ii][0] = 1;
        weights[ii][1] = 0;
      }
      else
      {
        // Find the coordinate index
        ///TODO make this fast
        // Note: lower_bound uses a binary search...  If we assume coordinates are
        // evenly spaced, we can speed things up considerably
        auto lower = std::lower_bound( m_coordinates[ii].begin(),
                                       m_coordinates[ii].end(),
                                       input[ii] );
        bounds[ii][1] = LvArray::integerConversion< localIndex >(
          std::distance( m_coordinates[ii].begin(), lower ) );
        bounds[ii][0] = bounds[ii][1] - 1;

        real64 dx =
          m_coordinates[ii][bounds[ii][1]] - m_coordinates[ii][bounds[ii][0]];
        weights[ii][0] =
          1.0 - ( input[ii] - m_coordinates[ii][bounds[ii][0]] ) / dx;
        weights[ii][1] = 1.0 - weights[ii][0];
      }
    }

    // Calculate the result
    for( localIndex ii = 0; ii < m_numCorners; ++ii )
    {
      // Find array index
      localIndex tableIndex = 0;
      for( localIndex jj = 0; jj < m_dimensions; ++jj )
      {
        tableIndex += bounds[jj][m_corners[jj][ii]] * m_indexIncrement[jj];
      }

      // Determine weighted value
      real64 cornerValue = m_values[tableIndex];
      for( localIndex jj = 0; jj < m_dimensions; ++jj )
      {
        cornerValue *= weights[jj][m_corners[jj][ii]];
      }
      result += cornerValue;
    }
  }
  // Nearest, Upper, Lower interpolation methods
  else
  {
    // Determine the index to the nearest table entry
    localIndex tableIndex = 0;
    for( localIndex ii = 0; ii < m_dimensions; ++ii )
    {
      // Determine the index along each table axis
      localIndex subIndex = 0;

      if( input[ii] <= m_coordinates[ii][0] )
      {
        // Coordinate is to the left of the table axis
        subIndex = 0;
      }
      else if( input[ii] >= m_coordinates[ii][m_size[ii] - 1] )
      {
        // Coordinate is to the right of the table axis
        subIndex = m_size[ii] - 1;
      }
      else
      {
        // Coordinate is within the table axis
        // Note: std::distance will return the index of the upper table vertex
        auto lower = std::lower_bound( m_coordinates[ii].begin(),
                                       m_coordinates[ii].end(),
                                       input[ii] );
        subIndex = LvArray::integerConversion< localIndex >(
          std::distance( m_coordinates[ii].begin(), lower ) );

        // Interpolation types:
        //   - Nearest returns the value of the closest table vertex
        //   - Upper returns the value of the next table vertex
        //   - Lower returns the value of the previous table vertex
        if( m_interpolationMethod == InterpolationType::Nearest )
        {
          if( ( input[ii] - m_coordinates[ii][subIndex - 1] ) <=
              ( m_coordinates[ii][subIndex] - input[ii] ) )
          {
            --subIndex;
          }
        }
        else if( m_interpolationMethod == InterpolationType::Lower )
        {
          if( subIndex > 0 )
          {
            --subIndex;
          }
        }
      }

      // Increment the global table index
      tableIndex += subIndex * m_indexIncrement[ii];
    }

    // Retrieve the nearest value
    result = m_values[tableIndex];
  }

  return result;
}

REGISTER_CATALOG_ENTRY( FunctionBase, TableFunction, std::string const &, Group * const )

}  // namespace geosx
