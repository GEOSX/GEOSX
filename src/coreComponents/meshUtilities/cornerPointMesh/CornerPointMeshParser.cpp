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
 * @file CPMeshParser.cpp
 */

#include "CornerPointMeshParser.hpp"

#include "meshUtilities/cornerPointMesh/utilities/StringUtilities.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include <iostream>

namespace geosx
{

namespace cornerPointMesh
{

CornerPointMeshParser::CornerPointMeshParser( string const & name )
  : m_meshName( name )
{}


void CornerPointMeshParser::readNumberOfCells( Path const & filePath,
                                               localIndex & nX,
                                               localIndex & nY,
                                               localIndex & nZ )
{
  // Note: This is definitely not what we want to do (reading the full file into a string)
  //       For now, to get a quick start, I ported this from PAMELA. I will improve that in a second step
  string const fileContent = cornerPointMeshStringUtilities::fileToString( filePath );
  std::istringstream meshFile;
  meshFile.str( fileContent );

  nX = 0;
  nY = 0;
  nZ = 0;

  // first read the dimensions of the mesh
  std::string line, buffer;
  while( getline( meshFile, line ) )
  {
    cornerPointMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", line );
    cornerPointMeshStringUtilities::removeExtraSpaces( line );
    cornerPointMeshStringUtilities::removeEndOfLine( line );
    cornerPointMeshStringUtilities::removeTab( line );
    cornerPointMeshStringUtilities::trim( line );
    if( line == "SPECGRID" || line == "DIMENS" )
    {
      std::vector< localIndex > bufInt;
      buffer = extractDataBelowKeyword( meshFile );
      cornerPointMeshStringUtilities::fromStringTo( buffer, bufInt );

      nX = bufInt[0];
      nY = bufInt[1];
      nZ = bufInt[2];
    }
  }
  GEOSX_THROW_IF( nX <= 0 || nY <= 0 || nZ <= 0,
                  "Neither SPECGRID nor DIMENS was found",
                  InputError );

}

void CornerPointMeshParser::readMesh( Path const & filePath,
                                      CornerPointMeshDimensions const & dims )
{
  // Note: This is definitely not what we want to do (reading the full file into a string)
  //       For now, to get a quick start, I ported this from PAMELA. I will improve that in a second step
  string const fileContent = cornerPointMeshStringUtilities::fileToString( filePath );
  std::istringstream meshFile;
  meshFile.str( fileContent );

  bool foundCOORD = false;
  bool foundZCORN = false;
  bool foundACTNUM = false;
  bool foundPERMX = false;
  bool foundPERMY = false;
  bool foundPERMZ = false;
  bool foundREGIONS = false;

  std::string line;
  while( getline( meshFile, line ) )
  {
    cornerPointMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", line );
    cornerPointMeshStringUtilities::removeExtraSpaces( line );
    cornerPointMeshStringUtilities::removeEndOfLine( line );
    cornerPointMeshStringUtilities::removeTab( line );
    cornerPointMeshStringUtilities::trim( line );

    // TODO: shorten, there must be a way to factorize all these functions
    //       at least, we shoud be able to merge readLocalACTNUM and readLocalPROP

    // at this point, SPECGRID has already been read, we skip it
    if( line == "COORD" )
    {
      foundCOORD = true;
      readLocalCOORD( meshFile, dims );
    }
    else if( line == "ZCORN" )
    {
      foundZCORN = true;
      readLocalZCORN( meshFile, dims );
    }
    else if( line == "ACTNUM" )
    {
      foundACTNUM = true;
      readLocalACTNUM( meshFile, dims );
    }
    else if( line == "PERMX" )
    {
      foundPERMX = true;
      readLocalPROP( meshFile, dims, m_permx );
    }
    else if( line == "PERMY" )
    {
      foundPERMY = true;
      readLocalPROP( meshFile, dims, m_permy );
    }
    else if( line == "PERMZ" )
    {
      foundPERMZ = true;
      readLocalPROP( meshFile, dims, m_permz );
    }
    else if( line == "PORO" )
    {
      readLocalPROP( meshFile, dims, m_poro );
    }
    else if( line == "REGIONS" )
    {
      foundREGIONS = true;
      readLocalPROP( meshFile, dims, m_regionId );
    }
  }

  if( !foundACTNUM )
  {
    fillLocalPROP( dims, m_actnum, static_cast< localIndex >( 1 ) ); // default actnum value = 1
  }
  if( !foundREGIONS )
  {
    fillLocalPROP( dims, m_regionId, static_cast< localIndex >( 0 ) ); // default region id = 0
  }

  // TODO: error out if one of them is missing
  if( foundPERMX && foundPERMY && foundPERMZ )
  {
    m_perm.resizeDimension< 0, 1 >( m_permx.size(), 3 );
    for( localIndex i = 0; i < m_permx.size(); ++i )
    {
      m_perm( i, 0 ) = m_permx( i );
      m_perm( i, 1 ) = m_permy( i );
      m_perm( i, 2 ) = m_permz( i );
    }
  }

  GEOSX_THROW_IF( !foundCOORD || !foundZCORN,
                  "At least one of the following keywords was not found: COORD, ZCORN",
                  InputError );


}

void CornerPointMeshParser::readLocalCOORD( std::istringstream & meshFile,
                                            CornerPointMeshDimensions const & dims )
{
  localIndex const nX = dims.nX();
  localIndex const nY = dims.nY();

  localIndex const maxSize = 6*(nX+1)*(nY+1); // TODO: change
  std::vector< real64 > tmp;
  tmp.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  cornerPointMeshStringUtilities::eclipseDataBufferToVector( buffer, tmp ); // TODO: change to save only what is necessary

  localIndex const nXLocal = dims.nXLocal();
  localIndex const nYLocal = dims.nYLocal();

  localIndex const iMinLocal = dims.iMinLocal();
  localIndex const jMinLocal = dims.jMinLocal();

  m_coord.resize( 6*(nXLocal+1)*(nYLocal+1) );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex j = 0; j < nYLocal+1; ++j )
  {
    for( localIndex i = 0; i < nXLocal+1; ++i )
    {
      localIndex const iPillarLocal = j*(nXLocal+1)+i;
      localIndex const iPillar = (j+jMinLocal)*(nX+1)+(i+iMinLocal);
      for( localIndex pos = 0; pos < 6; ++pos )
      {
        m_coord( 6*iPillarLocal+pos ) = tmp[ 6*iPillar+pos ]; // TODO: change to avoid doing that
      }
    }
  }
}

void CornerPointMeshParser::readLocalZCORN( std::istringstream & meshFile,
                                            CornerPointMeshDimensions const & dims )
{
  localIndex const nX = dims.nX();
  localIndex const nY = dims.nY();
  localIndex const nZ = dims.nZ();

  localIndex const maxSize = 8*nX*nY*nZ; // TODO: change
  std::vector< real64 > tmp;
  tmp.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  cornerPointMeshStringUtilities::eclipseDataBufferToVector( buffer, tmp ); // TODO: change to save only what is necessary

  localIndex const nXLocal = dims.nXLocal();
  localIndex const nYLocal = dims.nYLocal();
  localIndex const nZLocal = dims.nZLocal();

  localIndex const iMinLocal = dims.iMinLocal();
  localIndex const jMinLocal = dims.jMinLocal();
  localIndex const kMinLocal = 0;

  m_zcorn.resize( 8*nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const iXmLocal = k*8*nXLocal*nYLocal + j*4*nXLocal + 2*i;
        localIndex const iXpLocal = iXmLocal + 2*nXLocal;
        localIndex const iXm = (k+kMinLocal)*8*nX*nY + (j+jMinLocal)*4*nX + 2*(i+iMinLocal);
        localIndex const iXp = iXm + 2*nX;

        m_zcorn( iXmLocal ) = tmp[ iXm ];
        m_zcorn( iXmLocal + 1 ) = tmp[ iXm + 1 ];
        m_zcorn( iXpLocal ) = tmp[ iXp ];
        m_zcorn( iXpLocal + 1 ) = tmp[ iXp + 1 ];
        m_zcorn( iXmLocal + 4*nXLocal*nYLocal ) = tmp[ iXm + 4*nX*nY ];
        m_zcorn( iXmLocal + 4*nXLocal*nYLocal + 1 ) = tmp[ iXm + 4*nX*nY + 1 ];
        m_zcorn( iXpLocal + 4*nXLocal*nYLocal ) = tmp[ iXp + 4*nX*nY ];
        m_zcorn( iXpLocal + 4*nXLocal*nYLocal + 1 ) = tmp[ iXp + 4*nX*nY + 1 ];
      }
    }
  }
}

void CornerPointMeshParser::readLocalACTNUM( std::istringstream & meshFile,
                                             CornerPointMeshDimensions const & dims )
{
  localIndex const nX = dims.nX();
  localIndex const nY = dims.nY();
  localIndex const nZ = dims.nZ();

  localIndex const maxSize = nX*nY*nZ; // TODO: change
  std::vector< real64 > tmp;
  tmp.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  cornerPointMeshStringUtilities::eclipseDataBufferToVector( buffer, tmp ); // TODO: change to save only what is necessary
  std::replace( tmp.begin(), tmp.end(), 2, 0 );
  std::replace( tmp.begin(), tmp.end(), 3, 0 );

  localIndex const nXLocal = dims.nXLocal();
  localIndex const nYLocal = dims.nYLocal();
  localIndex const nZLocal = dims.nZLocal();

  localIndex const iMinLocal = dims.iMinLocal();
  localIndex const jMinLocal = dims.jMinLocal();
  localIndex const kMinLocal = 0;

  m_actnum.resize( nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs
  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const eiLocal = k*nXLocal*nYLocal + j*nXLocal + i;
        localIndex const ei = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + i+iMinLocal;
        m_actnum( eiLocal ) = tmp[ ei ];
      }
    }
  }
}

template< typename T >
void CornerPointMeshParser::readLocalPROP( std::istringstream & meshFile,
                                           CornerPointMeshDimensions const & dims,
                                           array1d< T > & prop )
{
  localIndex const nX = dims.nX();
  localIndex const nY = dims.nY();
  localIndex const nZ = dims.nZ();

  localIndex const maxSize = nX*nY*nZ; // TODO: change
  std::vector< real64 > tmp;
  tmp.reserve( maxSize );

  string buffer = extractDataBelowKeyword( meshFile ); // TODO: change to save only what is necessary
  cornerPointMeshStringUtilities::eclipseDataBufferToVector( buffer, tmp ); // TODO: change to save only what is necessary

  localIndex const nXLocal = dims.nXLocal();
  localIndex const nYLocal = dims.nYLocal();
  localIndex const nZLocal = dims.nZLocal();

  localIndex const iMinLocal = dims.iMinLocal();
  localIndex const jMinLocal = dims.jMinLocal();
  localIndex const kMinLocal = 0;

  prop.resize( nXLocal*nYLocal*nZLocal );

  // Note: Below I mimic what I will ultimately do: not storing the full file in a string,
  //       but instead going through the file and saving only what the MPI rank needs

  for( localIndex k = 0; k < nZLocal; ++k )
  {
    for( localIndex j = 0; j < nYLocal; ++j )
    {
      for( localIndex i = 0; i < nXLocal; ++i )
      {
        localIndex const eiLocal = k*nXLocal*nYLocal + j*nXLocal + i;
        localIndex const ei = (k+kMinLocal)*nX*nY + (j+jMinLocal)*nX + i+iMinLocal;
        prop( eiLocal ) = tmp[ ei ];
      }
    }
  }
}

template
void CornerPointMeshParser::readLocalPROP< real64 >( std::istringstream & meshFile,
                                                     CornerPointMeshDimensions const & dims,
                                                     array1d< real64 > & prop );
template
void CornerPointMeshParser::readLocalPROP< localIndex >( std::istringstream & meshFile,
                                                         CornerPointMeshDimensions const & dims,
                                                         array1d< localIndex > & prop );


template< typename T >
void CornerPointMeshParser::fillLocalPROP( CornerPointMeshDimensions const & dims,
                                           array1d< T > & prop,
                                           T defaultValue )
{
  localIndex const nXLocal = dims.nXLocal();
  localIndex const nYLocal = dims.nYLocal();
  localIndex const nZLocal = dims.nZLocal();

  prop.resize( nXLocal*nYLocal*nZLocal );
  prop.template setValues< serialPolicy >( defaultValue );
}

template
void CornerPointMeshParser::fillLocalPROP< real64 >( CornerPointMeshDimensions const & dims,
                                                     array1d< real64 > & prop,
                                                     real64 defaultValue );
template
void CornerPointMeshParser::fillLocalPROP< localIndex >( CornerPointMeshDimensions const & dims,
                                                         array1d< localIndex > & prop,
                                                         localIndex defaultValue );

std::string CornerPointMeshParser::extractDataBelowKeyword( std::istringstream & stringBlock )
{
  char keywordEnd = '/';
  std::string chunk;
  std::streampos pos;
  std::vector< std::string > res;
  getline( stringBlock, chunk, keywordEnd );
  stringBlock.clear();

  cornerPointMeshStringUtilities::removeStringAndFollowingContentFromLine( "--", chunk );
  cornerPointMeshStringUtilities::removeTab( chunk );
  cornerPointMeshStringUtilities::removeEndOfLine( chunk );
  cornerPointMeshStringUtilities::trim( chunk );
  res.push_back( chunk );
  stringBlock.ignore( 10, '\n' );
  getline( stringBlock, chunk );
  cornerPointMeshStringUtilities::removeTab( chunk );
  cornerPointMeshStringUtilities::removeEndOfLine( chunk );
  return res[0];
}

REGISTER_CATALOG_ENTRY( CornerPointMeshParser, CornerPointMeshParser, string const & )

} // namespace cornerPointMesh

} // namespace geosx
