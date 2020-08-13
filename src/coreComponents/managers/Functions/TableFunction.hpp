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

#ifndef GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_
#define GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_

#include "FunctionBase.hpp"

namespace geosx
{
/**
 * @class TableFunction
 *
 * An interface for a dense table-based function
 */
class TableFunction : public FunctionBase
{
public:
  /**
   * @brief The constructor
   * @param[in] name the name of this object manager
   * @param[in] parent the parent Group
   */
  TableFunction( const std::string & name, dataRepository::Group * const parent );

  /**
   * @brief The destructor
   */
  virtual ~TableFunction() override;

  /**
   * @brief The catalog name interface
   * @return name of the TableFunction in the FunctionBase catalog
   */
  static string
  CatalogName()
  {
    return "TableFunction";
  }

  /**
   * @brief Parse a table file.
   *
   * @tparam T The type for table or axis values.
   * @param[in] target The place to store values.
   * @param[in] filename The name of the file to read.
   * @param[in] delimiter The delimiter used for file entries.
   */
  template< typename T >
  void
  parse_file( array1d< T > & target, string const & filename, char delimiter );

  /**
   * @brief Initialize the table function
   */
  virtual void
  InitializeFunction() override;

  /**
   * @brief Build the maps used to evaluate the table function
   */
  void
  reInitializeFunction();

  /**
   * @brief Method to evaluate a function on a target object
   * @param group a pointer to the object holding the function arguments
   * @param time current time
   * @param set the subset of nodes to apply the function to
   * @param result an array to hold the results of the function
   */
  virtual void
  Evaluate( dataRepository::Group const * const group,
            real64 const time,
            SortedArrayView< localIndex const > const & set,
            real64_array & result ) const override final
  {
    FunctionBase::EvaluateT< TableFunction >( group, time, set, result );
  }

  /**
   * @brief Method to evaluate a function
   * @param input a scalar input
   * @return the function result
   */
  virtual real64
  Evaluate( real64 const * const input ) const override final;

  /**
   * @brief Get the table axes definitions
   * @return a reference to an array of arrays that define each table axis
   */
  array1d< real64_array > const &
  getCoordinates() const
  {
    return m_coordinates;
  }

  /**
   * @copydoc getCoordinates() const
   */
  array1d< real64_array > &
  getCoordinates()
  {
    return m_coordinates;
  }

  /**
   * @brief Get the table values
   * @return a reference to the 1d array of table values.  For ND arrays, values are stored in Fortran order.
   */
  array1d< real64 > const &
  getValues() const
  {
    return m_values;
  }

  /**
   * @copydoc getValues() const
   */
  array1d< real64 > &
  getValues()
  {
    return m_values;
  }

  /// Enumerator of available interpolation types
  enum class InterpolationType
  {
    Linear,
    Nearest,
    Upper,
    Lower
  };

private:
  /// Coordinates for 1D table
  real64_array m_tableCoordinates1D;

  /// List of table coordinate file names
  string_array m_coordinateFiles;

  /// Table voxel file names
  string m_voxelFile;

  /// Table interpolation method input string
  string m_interpolationMethodString;

  /// Table interpolation method
  InterpolationType m_interpolationMethod;

  /// An array of table axes
  array1d< real64_array > m_coordinates;

  /// Table values (in fortran order)
  real64_array m_values;

  /// Maximum number of table dimensions
  static localIndex constexpr m_maxDimensions = 4;

  /// Number of active table dimensions
  localIndex m_dimensions;

  /// Size of the table
  localIndex_array m_size;

  /// Array used to locate values within ND tables
  localIndex_array m_indexIncrement;

  /**
   * @brief The corners of the box that surround the value in N dimensions
   * m_corners should be of size m_maxDimensions x (2^m_maxDimensions)
   */
  localIndex m_corners[m_maxDimensions][16];

  /// The number of active table corners
  localIndex m_numCorners;
};

} /* namespace geosx */

#endif /* GEOSX_MANAGERS_FUNCTIONS_TABLEFUNCTION_HPP_ */
