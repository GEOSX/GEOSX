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

#ifndef GEOSX_COMMON_PATH_HPP
#define GEOSX_COMMON_PATH_HPP

// System includes
#include <string>
#include <vector>

namespace geosx
{
/*!
 * @brief Class describing a file Path
 * @details Purpose of this class is to
 * be used as a type to specify file path within the XML input files,
 * through the operator >>
 */
class Path : public std::string
{
public:
  /// Default constructor.
  Path() :
    std::string()
  {}

  /**
   * @brief Copy constructor, creates a copy of @p src.
   * @param src the Path to copy.
   */
  Path( Path const & src ) :
    std::string( src )
  {}

  /// Destructor.
  ~Path()
  {}

  /**
   * @brief Copy Constructor
   * @param rhs Reference to the Path that will be copied.
   * @return *this
   */
  Path &
  operator=( Path const & rhs )
  {
    std::string::operator=( rhs );
    return *this;
  }

  using std::string::string;

  /*!
   * @brief Get the path prefix of the file
   * @details The path prefix is usually a folder path
   * in which the XML file is located
   * @return the path prefix
   */
  static std::string &
  pathPrefix()
  {
    static std::string m_pathPrefix;
    return m_pathPrefix;
  }
};

/*!
 * @brief Gets the absolute path of a file
 * @param[in] path the relative path to the file
 * @param[out] absolutePath the absolute path to the file
 */
void
getAbsolutePath( std::string const & path, std::string & absolutePath );

/*!
 * @brief Tells wether the path is absolute of not
 * @param[in] path the input path
 * @retval true if the path is absolute
 * @retval false if the path is relative
 */
bool
isAbsolutePath( const std::string & path );

/*!
 * @brief Operator use with the class Path while parsing the XML file
 * @param[in,out] is the input stream
 * @param[in,out] p the path that will be set to an absolute path relative to the xml file
 * @return the input stream
 */
std::istream &
operator>>( std::istream & is, Path & p );

/*!
 * @brief Split the path in two parts : directory name and file name
 * @param[in] path the input path
 * @param[out] dirName name of the directory
 * @param[out] baseName the name of the file
 */
void
splitPath( std::string const & path,
           std::string & dirName,
           std::string & baseName );

/*!
 * @brief List all the files of one directory
 * @details Taken from http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
 * @param[in] path path to the directory
 * @param[out] files vector containing allt the file path
 */
void
readDirectory( std::string const & path, std::vector< std::string > & files );

/*!
 * @brief Make directories for @p path.
 * @param path The path to create.
 * @details Everything in @p path is intended to be a directory. If a directory in the path
 *   already exists nothing is done. if a directory doesn't exist it is created.
 */
void
makeDirsForPath( std::string const & path );

} /* end namespace geosx */

#endif
