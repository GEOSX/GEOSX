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
 * @file xmlWrapper.hpp
 */

#ifndef GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_
#define GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_

// Source includes
#include "common/DataTypes.hpp"
#include "dataRepository/DefaultValue.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "LvArray/src/output.hpp"
#include "LvArray/src/input.hpp"

// TPL includes
#include <pugixml.hpp>

// System includes
#include <algorithm>
#include <sstream>

namespace geosx
{
/**
 * @class xmlWrapper
 *
 * This class wraps/provides facilities to process entries from an xml file into the appropriate
 * data types. xmlWrapper provides some aliases that will wrap the underlying xml package being
 * used to extract data from the xml file, and a set of static member functions that facilitate
 * the parsing of string data from the xml into the variables that will hold those values in the
 * code.
 */
class xmlWrapper
{
public:
  /// Alias for the type of xml document.
  using xmlDocument = pugi::xml_document;

  /// Alias for the type of the result from an xml parse attempt.
  using xmlResult = pugi::xml_parse_result;

  /// Alias for the type of an xml node.
  using xmlNode = pugi::xml_node;

  /// Alias for the type of an xml attribute.
  using xmlAttribute = pugi::xml_attribute;

  /// Alias for the type variant of an xml node.
  using xmlTypes = pugi::xml_node_type;

  /**
   * @brief constexpr variable to hold name for inserting the file path into the xml file.
   *
   * This is used because we would like the option to hold the file path in the xml structure.
   */
  static constexpr auto filePathString = "filePath";

  /**
   * @name Constructors/destructor/assignment.
   */
  ///@{

  /// Deleted default constructor.
  xmlWrapper() = delete;

  /// Deleted default default destructor.
  ~xmlWrapper() = delete;

  /// Deleted copy constructor.
  xmlWrapper( xmlWrapper const & ) = delete;

  /// Deleted move constructor.
  xmlWrapper( xmlWrapper && ) = delete;

  ///@}

  /**
   * @brief Function to add xml nodes from included files.
   * @param targetNode the node for which to look for included children specifications
   *
   * This function looks for a "Included" node under the targetNode, loops over all subnodes under the "Included"
   * node, and then parses the file specified in those subnodes taking all the nodes in the file and adding them to
   * the targetNode.
   */
  static void
  addIncludedXML( xmlNode & targetNode );

  /**
   * @name String to variable parsing.
   *
   * These functions take in @p value and parse that string based on the type of
   * @p target. The function implementation should provide sufficient error checking
   * in the case that @p value is formatted incorrectly for the type specified in @p target.
   */
  ///@{

  /**
   * @brief Parse a string and fill a variable with the value(s) in the string.
   * @tparam T the type of variable fill with string value
   * @param[out] target the object to read values into
   * @param[in]  value  the string that contains the data to be parsed into target
   * @return void.
   */
  template< typename T >
  static std::enable_if_t< traits::CanStreamInto< std::istringstream, T > >
  StringToInputVariable( T & target, string const & value )
  {
    std::istringstream ss( value );
    ss >> target;
  }

  /**
   * @brief Parse a string and fill a R1Tensor with the value(s) in the string.
   * @param[out] target the object to read values into
   * @param[in]  value  the string that contains the data to be parsed into target
   */
  static void
  StringToInputVariable( R1Tensor & target, string const & value );

  /**
   * @brief Parse a string and fill an Array with the value(s) in the string.
   * @tparam T    data type of the array
   * @tparam NDIM number of dimensions of the array
   * @tparam PERMUTATION the permutation of the array
   * @param[out] array the array to read values into
   * @param[in]  value the string that contains the data to be parsed into target
   * @return void.
   */
  template< typename T, int NDIM, typename PERMUTATION >
  static std::enable_if_t< traits::CanStreamInto< std::istringstream, T > >
  StringToInputVariable( Array< T, NDIM, PERMUTATION > & array, string const & value )
  {
    LvArray::input::stringToArray( array, value );
  }

  ///@}

  /// Defines a static constexpr bool canParseVariable that is true iff the template parameter T
  /// is a valid argument to StringToInputVariable.
  IS_VALID_EXPRESSION( canParseVariable,
                       T,
                       StringToInputVariable( std::declval< T & >(), std::string() ) );

  /**
   * @name Attribute extraction from XML nodes.
   */
  ///@{

  /**
   * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
   * @tparam T             the type of variable fill with xml attribute.
   * @tparam T_DEF         type of the default value for @p rval
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T, typename T_DEF = T >
  static std::enable_if_t< canParseVariable< T >, bool >
  ReadAttributeAsType(
    T & rval,
    string const & name,
    xmlNode const & targetNode,
    T_DEF const & defVal )
  {
    pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );
    if( !xmlatt.empty() )
    {
      // parse the string/attribute into a value
      StringToInputVariable( rval, xmlatt.value() );
    }
    else
    {
      // set the value to the default value
      equate( rval, defVal );
    }
    return true;
  }

  /**
   * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] required   whether or not the value is required
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static std::enable_if_t< canParseVariable< T >, bool >
  ReadAttributeAsType(
    T & rval,
    string const & name,
    xmlNode const & targetNode,
    bool const required )
  {
    pugi::xml_attribute xmlatt = targetNode.attribute( name.c_str() );

    bool const success = !( xmlatt.empty() && required );

    if( success )
    {
      // parse the string/attribute into a value
      StringToInputVariable( rval, xmlatt.value() );
    }
    return success;
  }

  /**
   * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static std::enable_if_t< !dataRepository::DefaultValue< T >::has_default_value, bool >
  ReadAttributeAsType( T & rval,
                       string const & name,
                       xmlNode const & targetNode,
                       dataRepository::DefaultValue< T > const & )
  {
    return ReadAttributeAsType( rval, name, targetNode, false );
  }

  /**
   * @brief Extract attribute in an xml tree, and translate its value into a typed variable.
   * @tparam T             the type of variable fill with xml attribute.
   * @param[out] rval      the variable to fill with value
   * @param[in] name       the name of the xml attribute to process
   * @param[in] targetNode the xml node that should contain the attribute
   * @param[in] defVal     default value of @p rval (or entries of @p rval, if it is an array)
   * @return boolean value indicating whether the value was successfully read from XML.
   */
  template< typename T >
  static
    typename std::enable_if_t< dataRepository::DefaultValue< T >::has_default_value, bool >
    ReadAttributeAsType( T & rval,
                         string const & name,
                         xmlNode const & targetNode,
                         dataRepository::DefaultValue< T > const & defVal )
  {
    return ReadAttributeAsType( rval, name, targetNode, defVal.value );
  }

  /**
   * @brief Stub that for unreadable types that errors out.
   * @return false.
   */
  template< typename T, typename U >
  static std::enable_if_t< !canParseVariable< T >, bool >
  ReadAttributeAsType( T &, string const &, xmlNode const &, U const & )
  {
    GEOSX_ERROR( "Cannot parse the given type "
                 << LvArray::system::demangleType< T >() );
    return false;
  }

  ///@}

private:
  /**
   * @brief Set @p lhs equal to @p rhs.
   * @tparam T The type of @p lhs and @p rhs.
   * @param lhs The value to set to @p rhs.
   * @param rhs The value to set @p lhs to.
   */
  template< typename T >
  static void
  equate( T & lhs, T const & rhs )
  {
    lhs = rhs;
  }

  /**
   * @brief Set the entries of @p lhs equal to @p rhs.
   * @tparam T The type of the values in @p lhs and @p rhs.
   * @tparam NDIM The dimension of @p lhs.
   * @tparam PERM The permutation of @p rhs.
   * @param lhs The array of value to set to @p rhs.
   * @param rhs The value to set @p lhs to.
   */
  template< typename T, int NDIM, typename PERM >
  static void
  equate( Array< T, NDIM, PERM > const & lhs, T const & rhs )
  {
    lhs.template setValues< serialPolicy >( rhs );
  }
};

} /* namespace geosx */

#endif /*GEOSX_DATAREPOSITORY_XMLWRAPPER_HPP_ */
