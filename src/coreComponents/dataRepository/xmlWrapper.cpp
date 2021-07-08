/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file xmlWrapper.cpp
 */

#include "xmlWrapper.hpp"

namespace geosx
{
using namespace dataRepository;

template< typename T, int SIZE >
void xmlWrapper::stringToInputVariable( Tensor< T, SIZE > & target, string const & inputValue )
{
  std::istringstream ss( inputValue );

  real64 value;
  int count = 0;
  while( ss.peek() == '{' || ss.peek() == ' ' )
  {
    ss.ignore();
  }
  while( !((ss>>value).fail()) )
  {
    target[count++] = value;
    while( ss.peek() == ',' || ss.peek() == ' ' || ss.peek() == '}' )
    {
      ss.ignore();
    }
  }
  GEOSX_ERROR_IF( count!=SIZE, "incorrect number of components specified for Tensor" );
}


template void xmlWrapper::stringToInputVariable( Tensor< real64, 3 > & target, string const & inputValue );
template void xmlWrapper::stringToInputVariable( Tensor< real64, 6 > & target, string const & inputValue );

void xmlWrapper::addIncludedXML( xmlNode & targetNode )
{

  xmlNode rootNode = targetNode.root();
  string path = rootNode.child( filePathString ).attribute( filePathString ).value();

  xmlNode includedNode = targetNode.child( "Included" );

  for( xmlWrapper::xmlNode childNode=includedNode.first_child(); childNode; childNode=childNode.next_sibling())
  {
    // Get the child tag and name
    string childName = childNode.name();
    GEOSX_ERROR_IF( childName!="File", "Child nodes of \"Included\" should be named \"File\" " );

    string filePathName = childNode.attribute( "name" ).value();
    if( filePathName[0] != '/' )
    {
      filePathName = path + filePathName;
    }
    xmlDocument includedXmlDocument;
    xmlResult result;
    result = includedXmlDocument.load_file( filePathName.c_str());
    GEOSX_ERROR_IF( !result, "Attempt to include file ("<<filePathName.c_str()<<") failed\n" );

    // To validate correctly, included files should contain the root Problem node
    xmlNode includedRootNode = includedXmlDocument.child( "Problem" );
    for( xmlNode importNode=includedRootNode.first_child(); importNode; importNode=importNode.next_sibling())
    {
      targetNode.append_copy( importNode );
    }
  }
}

void xmlWrapper::buildMultipleInputXML( string_array const & inputFileList, string & inputFileName )
{
  if( inputFileList.size() == 1 )
  {
    inputFileName = inputFileList[0];
  }
  else
  {
    // Write the composite xml file
    inputFileName = "composite_input.xml";
    xmlWrapper::xmlDocument compositeTree;
    xmlWrapper::xmlNode compositeRoot = compositeTree.append_child( "Problem" );
    xmlWrapper::xmlNode includedRoot = compositeRoot.append_child( "Included" );

    for( auto & fileName: inputFileList )
    {
      xmlWrapper::xmlNode fileNode = includedRoot.append_child( "File" );
      fileNode.append_attribute( "name" ) = fileName.c_str();
    }

    compositeTree.save_file( inputFileName.c_str() );
  }
}


} /* namespace geosx */
