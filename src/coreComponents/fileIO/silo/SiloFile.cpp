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
 * @file SiloFile.cpp
 */


#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"

#include "SiloFile.hpp"

#include "codingUtilities/Utilities.hpp"
#include "common/DataTypes.hpp"
#include "constitutive/ConstitutiveManager.hpp"
#include "constitutive/fluid/SingleFluidBase.hpp"
#include "constitutive/fluid/MultiFluidBase.hpp"
#include "constitutive/contact/ContactRelationBase.hpp"
#include "constitutive/NullModel.hpp"
#include "mesh/DomainPartition.hpp"
#include "mesh/MeshBody.hpp"
#include "common/MpiWrapper.hpp"

#include <iostream>

#include <sstream>
#include <algorithm>
#include <iterator>
#include <sys/stat.h>



#if !defined(GEOSX_USE_MPI)
int MPI_Comm_size( MPI_Comm, int * size ) { *size=1; return 0; }
int MPI_Comm_rank( MPI_Comm, int * rank ) { *rank=1; return 0; }

int MPI_Ssend( const void *, int, MPI_Datatype, int, int,
               MPI_Comm )
{
  return 0;
}

int MPI_Recv( void * buf, int, MPI_Datatype, int, int,
              MPI_Comm, MPI_Status * )
{
  *reinterpret_cast< int * >(buf) = 0;
  return 0;
}
#endif
#include "pmpio.h"

/// forward declaration of NodeManagerT for use as a template argument
class NodeManager;

namespace geosx
{

/**
 *
 * @return
 */
namespace SiloFileUtilities
{

template<> int DB_TYPE< int >()
{
  return DB_INT;
}
template<> int DB_TYPE< unsigned int >()
{
  return DB_INT;
}
template<> int DB_TYPE< float >()
{
  return DB_FLOAT;
}
template<> int DB_TYPE< real64 >()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE< R1Tensor >()
{
  return DB_DOUBLE;
}
template<> int DB_TYPE< unsigned long >()
{
  return DB_LONG;
}
template<> int DB_TYPE< long >()
{
  return DB_LONG;
}
template<> int DB_TYPE< long long >()
{
  return DB_LONG_LONG;
}
template<> int DB_TYPE< string >()
{
  return DB_CHAR;
}

template<> int GetNumberOfVariablesInField< int >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< unsigned int >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< unsigned long >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< long >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< float >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< real64 >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< long long unsigned int >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< long long int >()
{
  return 1;
}
template<> int GetNumberOfVariablesInField< R1Tensor >()
{
  return R1Tensor::SIZE;
}
template<> int GetNumberOfVariablesInField< string >()
{
  return 1;
}

template< typename TYPE >
void SetVariableNames( string const & fieldName, string_array & varnamestring, char const * varnames[] )
{
  varnamestring.resize( GetNumberOfVariablesInField< TYPE >());
  varnamestring[0] = fieldName;
  varnames[0] = varnamestring[0].c_str();
}
template void SetVariableNames< int >( string const & fieldName,
                                       string_array & varnamestring,
                                       char const * varnames[] );
template void SetVariableNames< unsigned long >( string const & fieldName,
                                                 string_array & varnamestring,
                                                 char const * varnames[] );
template void SetVariableNames< real64 >( string const & fieldName,
                                          string_array & varnamestring,
                                          char const * varnames[] );
template void SetVariableNames< long long unsigned int >( string const & fieldName,
                                                          string_array & varnamestring,
                                                          char const * varnames[] );



template<>
void SetVariableNames< R1Tensor >( string const & fieldName,
                                   string_array & varnamestring,
                                   char const * varnames[] )
{
  varnamestring.resize( GetNumberOfVariablesInField< R1Tensor >());
  varnamestring[0] = fieldName + "_1";
  varnamestring[1] = fieldName + "_2";
  varnamestring[2] = fieldName + "_3";
  varnames[0] = const_cast< char * >( varnamestring[0].c_str() );
  varnames[1] = const_cast< char * >( varnamestring[1].c_str() );
  varnames[2] = const_cast< char * >( varnamestring[2].c_str() );
}


template<> int GetTensorRank< int >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< unsigned int >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< long >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< unsigned long >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< real32 >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< real64 >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< R1Tensor >()
{
  return DB_VARTYPE_VECTOR;
}
template<> int GetTensorRank< long long unsigned int >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< long long int >()
{
  return DB_VARTYPE_SCALAR;
}
template<> int GetTensorRank< string >()
{
  return DB_VARTYPE_SCALAR;
}



}



using namespace constitutive;
using namespace dataRepository;

// *********************************************************************************************************************
/// Default Constructor
SiloFile::SiloFile():
  m_dbFilePtr( nullptr ),
  m_dbBaseFilePtr( nullptr ),
  m_numGroups( 1 ),
  m_baton( nullptr ),
  m_driver( DB_HDF5 ),
  m_plotFileRoot( "plot" ),
  m_restartFileRoot( "restart" ),
  m_siloDirectory( "siloFiles" ),
  m_fileName(),
  m_baseFileName(),
  m_emptyMeshes(),
//  m_emptyMaterials(),
  m_emptyVariables(),
  m_writeEdgeMesh( 0 ),
  m_writeFaceMesh( 0 ),
  m_writeCellElementMesh( 1 ),
  m_writeFaceElementMesh( 1 ),
  m_plotLevel( dataRepository::PlotLevel::LEVEL_1 ),
  m_ghostFlags( true )
{
  DBSetAllowEmptyObjects( 1 );
}

// *********************************************************************************************************************
/// Destructor
SiloFile::~SiloFile()
{}

// *********************************************************************************************************************

void SiloFile::makeSiloDirectories()
{

  int rank=0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif

  if( rank==0 )
  {
    makeDirsForPath( joinPath( m_siloDirectory, m_siloDataSubDirectory ) );
  }
}

/**
 *
 */
void SiloFile::initialize( int const MPI_PARAM( numGroups ) )
{
  makeSiloDirectories();

#ifdef GEOSX_USE_MPI
  // Ensure all procs agree on numGroups, driver and file_ext
  m_numGroups = numGroups;
#else
  m_numGroups = 1;
#endif
  MpiWrapper::bcast( &m_numGroups, 1, 0, MPI_COMM_GEOSX );
//  MPI_Bcast( const_cast<int*>(&m_driver), 1, MPI_INT, 0, MPI_COMM_GEOSX);
  // Initialize PMPIO, pass a pointer to the driver type as the user data.
  m_baton = PMPIO_Init( m_numGroups,
                        PMPIO_WRITE,
                        MPI_COMM_GEOSX,
                        1,
                        PMPIO_DefaultCreate,
                        PMPIO_DefaultOpen,
                        PMPIO_DefaultClose,
                        const_cast< int * >(&m_driver));
}

// *********************************************************************************************************************
/**
 *
 */
void SiloFile::finish()
{
  PMPIO_Finish( m_baton );
}

int SiloFile::groupRank( int const i ) const
{
  return PMPIO_GroupRank( m_baton, i );
}

// *********************************************************************************************************************
/**
 *
 * @param domainNumber
 * @param cycleNum
 * @param fileName
 * @param nsName
 */
void SiloFile::waitForBatonWrite( int const domainNumber,
                                  int const cycleNum,
                                  integer const eventCounter,
                                  bool const isRestart )
{

  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  int const groupRank = PMPIO_GroupRank( m_baton, rank );
  char fileName[200] = { 0 };
  char baseFileName[200] = { 0 };
  char dirName[200] = { 0 };



  if( isRestart )
  {
    // The integrated test repo does not use the eventProgress indicator, so skip it for now
    sprintf( baseFileName, "%s_%06d", m_restartFileRoot.c_str(), cycleNum );
    sprintf( fileName, "%s%s%s_%06d.%03d",
             m_siloDataSubDirectory.c_str(), "/", m_restartFileRoot.c_str(), cycleNum, groupRank );
  }
  else
  {
    sprintf( baseFileName, "%s_%06d%02d",
             m_plotFileRoot.c_str(),
             cycleNum,
             eventCounter );

    sprintf( fileName,
             "%s_%06d%02d.%03d",
             m_plotFileRoot.c_str(),
             cycleNum,
             eventCounter,
             groupRank );
  }
  sprintf( dirName, "domain_%05d", domainNumber );

  string const dataFilePathAndName = joinPath( m_siloDirectory, m_siloDataSubDirectory, fileName );
  m_dbFilePtr = static_cast< DBfile * >( PMPIO_WaitForBaton( m_baton, dataFilePathAndName.c_str(), dirName ) );

  m_fileName = fileName;
  m_baseFileName = baseFileName;

  if( rank==0 )
  {
    m_dbBaseFilePtr = DBCreate( (m_siloDirectory + "/"+ m_baseFileName).c_str(), DB_CLOBBER, DB_LOCAL, nullptr, DB_HDF5 );
//    m_dbBaseFilePtr = DBOpen( m_baseFileName.c_str(), DB_HDF5, DB_APPEND );
  }
}


void SiloFile::waitForBaton( int const domainNumber, string const & restartFileName )
{

  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  int const groupRank = PMPIO_GroupRank( m_baton, rank );
  char fileName[200] = { 0 };
  char baseFileName[200] = { 0 };
  char dirName[200] = { 0 };


  sprintf( baseFileName, "%s", restartFileName.c_str());
  if( groupRank == 0 )
    sprintf( fileName, "%s", restartFileName.c_str());
  else
  {
    if( m_siloDirectory.empty())
    {
      sprintf( fileName, "%s.%03d", restartFileName.c_str(), groupRank );
    }
    else
    {
      sprintf( fileName, "%s%s%s.%03d", m_siloDirectory.c_str(), "/", restartFileName.c_str(), groupRank );
    }

  }

  sprintf( dirName, "domain_%05d", domainNumber );

  m_dbFilePtr = (DBfile *) PMPIO_WaitForBaton( m_baton, fileName, dirName );

  m_fileName = fileName;
  m_baseFileName = baseFileName;
}
/**
 *
 */
void SiloFile::handOffBaton()
{
  PMPIO_HandOffBaton( m_baton, m_dbFilePtr );

  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  if( rank==0 )
  {
    DBClose( m_dbBaseFilePtr );
  }
}

/**
 *
 * @param meshName
 * @param nnodes
 * @param coords
 * @param globalNodeNum
 * @param numRegions
 * @param shapecnt
 * @param meshConnectivity
 * @param globalElementNum
 * @param ghostFlag
 * @param shapetype
 * @param shapesize
 * @param cycleNumber
 * @param problemTime
 */
void SiloFile::writeMeshObject( string const & meshName,
                                const localIndex nnodes,
                                real64 * coords[3],
                                globalIndex const * const globalNodeNum,
                                char const * const ghostNodeFlag,
                                char const * const ghostZoneFlag,
                                int const numShapes,
                                int const * shapecnt,
                                const localIndex * const * const meshConnectivity,
                                globalIndex const * const * const globalElementNum,
                                int const * const shapetype,
                                int const * const shapesize,
                                int const cycleNumber,
                                real64 const problemTime )
{

  const DBdatatype datatype = DB_DOUBLE;
  int const one = 1;

  DBoptlist * optlist = DBMakeOptlist( 5 );
  if( globalNodeNum!=nullptr )
  {
    if( std::is_same< globalIndex, int >::value || std::is_same< globalIndex, long long >::value )
    {
      DBAddOption( optlist, DBOPT_NODENUM, const_cast< globalIndex * >(globalNodeNum));
      if( std::is_same< globalIndex, long long >::value )
      {
        DBAddOption( optlist, DBOPT_LLONGNZNUM, const_cast< int * >(&one) );
      }
    }
  }
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));
  if( m_ghostFlags )
  {
    DBAddOption( optlist, DBOPT_GHOST_NODE_LABELS, const_cast< char * >(ghostNodeFlag) );
  }

  int numTotZones = 0;
  int lnodelist = 0;
  for( int i = 0; i < numShapes; ++i )
  {
    numTotZones += shapecnt[i];
    //  if shapesize <= 0, that signals that we are using arbitrary polygons.
    if( shapesize[i] > 0 )
      lnodelist += shapecnt[i] * shapesize[i];
    else
      lnodelist += -shapesize[i];
  }


  if( numTotZones == 0 )
  {
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );
    string emptyObject = pwd;
    emptyObject += "/" + meshName;
    m_emptyMeshes.emplace_back( emptyObject );
  }
  else
  {

    string zonelistName;
    zonelistName = meshName + "_zonelist";

    DBPutUcdmesh( m_dbFilePtr, meshName.c_str(), 3, nullptr, (float * *) coords, nnodes, numTotZones,
                  zonelistName.c_str(), nullptr, datatype, optlist );

    DBClearOptlist( optlist );

    array1d< integer > nodelist( lnodelist );
    globalIndex_array globalZoneNumber( lnodelist );

    int count = 0;
    int elemCount = 0;
    for( int i = 0; i < numShapes; ++i )
    {
      int n;
      if( shapesize[i] > 0 )
        n = shapecnt[i] * shapesize[i];
      else
        n = -shapesize[i];
      for( int j = 0; j < n; ++j )
      {
        nodelist[count++] = meshConnectivity[i][j];
      }
    }

    if( globalElementNum != nullptr )
    {
      for( int i = 0; i < numShapes; ++i )
      {
        if( std::is_same< globalIndex, int >::value || std::is_same< globalIndex, long long >::value )
        {
          for( int j = 0; j < shapecnt[i]; ++j )
          {
            globalZoneNumber[elemCount++] = globalElementNum[i][j];
          }
          // write zonelist
          DBAddOption( optlist, DBOPT_ZONENUM, const_cast< globalIndex * >(globalZoneNumber.data()));
          if( std::is_same< globalIndex, long long >::value )
          {
            DBAddOption( optlist, DBOPT_LLONGNZNUM, const_cast< int * >(&one));
          }
        }
      }
    }

    integer_array shapesize2( numShapes );
    for( int i = 0; i < numShapes; ++i )
    {
      if( shapesize[i] < 0 )
      {
        shapesize2[i] = 0;
      }
      else
        shapesize2[i] = shapesize[i];
    }

    int hi_offset = 0;


    if( m_ghostFlags )
    {
      DBAddOption( optlist, DBOPT_GHOST_ZONE_LABELS, const_cast< char * >( ghostZoneFlag ) );
    }

    DBPutZonelist2( m_dbFilePtr, zonelistName.c_str(), numTotZones, 3, nodelist.data(), lnodelist, 0, 0,
                    hi_offset, const_cast< int * >(shapetype), const_cast< int * >(shapesize2.data()),
                    const_cast< int * >(shapecnt), numShapes,
                    optlist );

    DBClearOptlist( optlist );
  }

  // write multimesh object
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  if( rank == 0 )
  {
    DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
    DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

    writeMultiXXXX( DB_UCDMESH, DBPutMultimesh, 0, meshName, cycleNumber, "/", optlist );
  }

  DBFreeOptlist( optlist );
}


void SiloFile::writeBeamMesh( string const & meshName,
                              const localIndex nnodes,
                              real64 * coords[3],
                              const localIndex_array & node1,
                              const localIndex_array & node2,
                              int const cycleNumber,
                              real64 const problemTime )
{
  // Connectivity.
  integer_array nodelist;
  nodelist.reserve( 2*node1.size());
  for( localIndex i = 0; i < node1.size(); ++i )
  {
    nodelist.emplace_back( static_cast< int >(node1[i]));
    nodelist.emplace_back( static_cast< int >(node2[i]));
  }

  writeBeamMesh( meshName, nnodes, coords, nodelist, cycleNumber, problemTime );
}

void SiloFile::writeBeamMesh( string const & meshName,
                              const localIndex nnodes,
                              real64 * coords[3],
                              integer_array & nodelist,
                              int const cycleNumber,
                              real64 const problemTime )
{
  const DBdatatype datatype = DB_DOUBLE;

  DBoptlist * optlist = DBMakeOptlist( 4 );
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

  int lnodelist = nodelist.size();
  int nshapetypes = 1;
  int nzones = lnodelist/2;
  int shapecnt[] = {nzones};
  int shapesize[] = {2};
  int shapetype[] = {DB_ZONETYPE_BEAM};

  // Write out connectivity information.
  int const origin = 0;
  int const lo_offset = 0;
  int const hi_offset = 0;
  int const ndims = 3;
  const char * coordnames[3]  = { "xcoords", "ycoords", "zcoords" };

  // Write out connectivity information.
  DBPutZonelist2( m_dbFilePtr, "zonelist", nzones, ndims, lnodelist > 0 ? &nodelist[0] : nullptr,
                  lnodelist, origin, lo_offset, hi_offset,
                  shapetype, shapesize, shapecnt, nshapetypes, optlist );

  // Write an unstructured mesh.
  DBPutUcdmesh( m_dbFilePtr, meshName.c_str(), ndims, const_cast< char * * >(coordnames), &coords[0], nnodes, nzones,
                "zonelist", nullptr, datatype, nullptr );
  DBClearOptlist( optlist );

  //----write multimesh object
  {
    int rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  #endif
    if( rank == 0 )
    {
      DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
      DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));
      writeMultiXXXX( DB_UCDMESH, DBPutMultimesh, 0, meshName.c_str(), cycleNumber, "/", optlist );
    }
  }

  //---free the option list
  DBFreeOptlist( optlist );
}



void SiloFile::writePointMesh( string const & meshName,
                               const localIndex numPoints,
                               real64 * coords[3],
                               int const cycleNumber,
                               real64 const problemTime )
{
  const DBdatatype datatype = DB_DOUBLE;
  DBoptlist * optlist = DBMakeOptlist( 2 );
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));
  DBPutPointmesh ( m_dbFilePtr, meshName.c_str(), 3, coords, numPoints, datatype, optlist );


  //----write multimesh object
  {
    int rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  #endif
    if( rank == 0 )
    {
      writeMultiXXXX( DB_POINTMESH, DBPutMultimesh, 0, meshName.c_str(), cycleNumber, "/", optlist );
    }
  }

}

void SiloFile::writeMaterialMapsFullStorage( ElementRegionBase const & elemRegion,
                                             string const & meshName,
                                             string_array const & regionMaterialList,
                                             int const cycleNumber,
                                             real64 const problemTime )
{
  string const name = meshName + "_materials";

  int const nmat = regionMaterialList.size();

  if( nmat > 0 )
  {
    array1d< int > matnos( nmat );
    std::vector< string > materialNameStrings( nmat );
    array1d< char const * > materialNames( nmat+1 );
    materialNames.back() = nullptr;

    for( int matIndex=0; matIndex<nmat; ++matIndex )
    {
      matnos[matIndex] = matIndex;
      materialNameStrings[matIndex] = regionMaterialList[matIndex];
      materialNames[matIndex] = materialNameStrings[matIndex].c_str();
    }

    int ndims = 1;
    int dims = 0;
    int mixlen=0;

    elemRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
    {
      if( nmat > 1 )
      {
        mixlen += subRegion.size() * nmat;
      }
      dims += subRegion.size();
    } );

    array1d< integer > matlist( dims );
    array1d< integer > mix_zone( mixlen );
    array1d< integer > mix_mat( mixlen );
    array1d< integer > mix_next( mixlen );
    array1d< real64 > mix_vf( mixlen );

    int elemCount = 0;
    int mixCount = 0;

    if( nmat > 0 )
    {
      elemRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
      {
        if( nmat == 1 )
        {
          for( localIndex k = 0; k < subRegion.size(); ++k )
          {
            matlist[elemCount++] = 0;
          }
        }
        else
        {
          for( localIndex k = 0; k < subRegion.size(); ++k )
          {
            matlist[elemCount++] = -(mixCount+1);
            for( localIndex a=0; a<nmat; ++a )
            {
              mix_zone[mixCount] = k;
              mix_mat[mixCount] = a;
              mix_vf[mixCount] = 1.0/nmat;
              if( a == nmat-1 )
              {
                mix_next[mixCount] = 0;
              }
              else
              {
                mix_next[mixCount] = mixCount+2;
              }
              ++mixCount;
            }
          }
        }
      } );
    }

//    if( nmat == 0 )
//    {
//      char pwd[256];
//      DBGetDir(m_dbFilePtr, pwd);
//      string emptyObject = pwd;
//      emptyObject += "/" + fieldName;
//      m_emptyVariables.emplace_back(emptyObject);
//    }
//    else
    {
      DBoptlist * optlist = DBMakeOptlist( 3 );
      DBAddOption( optlist, DBOPT_MATNAMES, materialNames.data());
      DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
      DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

      DBPutMaterial( m_dbFilePtr,
                     name.c_str(),
                     meshName.c_str(),
                     nmat,
                     matnos.data(),
                     matlist.data(),
                     &dims,
                     ndims,
                     mix_next.data(),
                     mix_mat.data(),
                     mix_zone.data(),
                     mix_vf.data(),
                     mixlen,
                     DB_DOUBLE,
                     optlist );

      DBFreeOptlist( optlist );
    }
    // write multimesh object
    int rank = 0;
  #ifdef GEOSX_USE_MPI
    MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
  #endif
    if( rank == 0 )
    {

      int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );

      array1d< string > vBlockNames( size );
      std::vector< char * > BlockNames( size );
      char tempBuffer[1024];
      char currentDirectory[256];

      DBGetDir( m_dbBaseFilePtr, currentDirectory );
      DBSetDir( m_dbBaseFilePtr, "/" );

      for( int i = 0; i < size; ++i )
      {
        int groupRank = PMPIO_GroupRank( m_baton, i );

        /* this mesh block is another file */
        sprintf( tempBuffer,
                 "%s%s%s.%03d:/domain_%05d/%s",
                 m_siloDataSubDirectory.c_str(),
                 "/",
                 m_baseFileName.c_str(),
                 groupRank,
                 i,
                 name.c_str() );

        vBlockNames[i] = tempBuffer;
        BlockNames[i] = const_cast< char * >( vBlockNames[i].c_str() );
      }

      {
        DBoptlist * optlist = DBMakeOptlist( 5 );
        DBAddOption( optlist, DBOPT_MATNAMES, materialNames.data());
        DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
        DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));
        DBAddOption( optlist, DBOPT_NMATNOS, const_cast< int * >(&nmat) );
        DBAddOption( optlist, DBOPT_MATNOS, matnos.data() );

        DBPutMultimat( m_dbBaseFilePtr, name.c_str(), size, BlockNames.data(),
                       const_cast< DBoptlist * >(optlist));
        DBFreeOptlist( optlist );

      }

      DBSetDir( m_dbBaseFilePtr, currentDirectory );

    }


    string subDirectory = meshName + "_MaterialFields";
    string rootDirectory = "/" + subDirectory;

    {
      string shortsubdir( subDirectory );
      string::size_type pos = subDirectory.find_last_of( "//" );

      if( pos != shortsubdir.npos )
      {
        shortsubdir.erase( 0, pos+1 );
      }


      makeSubDirectory( shortsubdir, rootDirectory );
      DBSetDir( m_dbFilePtr, shortsubdir.c_str());

    }

    std::set< std::pair< string, WrapperBase const * > > fieldNames;
    for( localIndex matI=0; matI<nmat; ++matI )
    {
      Group const & constitutiveModel = elemRegion.getSubRegion( 0 ).getConstitutiveModel( regionMaterialList[matI] );

      for( auto const & wrapperIter : constitutiveModel.wrappers() )
      {
        auto const & wrapper = wrapperIter.second;

        if( wrapper->getPlotLevel() < m_plotLevel )
        {
          std::type_info const & typeID = wrapper->getTypeId();

          if( typeID == typeid( array2d< real64 > ) ||
              typeID == typeid( array3d< real64 > ) ||
              typeID == typeid( array4d< real64 > ) )
          {
            fieldNames.insert( std::make_pair( wrapper->getName(), wrapper ) );
          }
        }
      }
    }

    for( auto fieldName : fieldNames )
    {
      if( fieldName.second->getTypeId() == typeid( array2d< real64 >))
      {
        writeMaterialDataField2d< real64, real64 >( meshName,
                                                    fieldName.first,
                                                    elemRegion,
                                                    DB_ZONECENT,
                                                    cycleNumber,
                                                    problemTime,
                                                    rootDirectory,
                                                    regionMaterialList );
      }
      if( fieldName.second->getTypeId() == typeid( array3d< real64 >))
      {
        writeMaterialDataField3d< real64, real64 >( meshName,
                                                    fieldName.first,
                                                    elemRegion,
                                                    DB_ZONECENT,
                                                    cycleNumber,
                                                    problemTime,
                                                    rootDirectory,
                                                    regionMaterialList );
      }
      if( fieldName.second->getTypeId() == typeid( array4d< real64 >))
      {
        writeMaterialDataField4d< real64, real64 >( meshName,
                                                    fieldName.first,
                                                    elemRegion,
                                                    DB_ZONECENT,
                                                    cycleNumber,
                                                    problemTime,
                                                    rootDirectory,
                                                    regionMaterialList );
      }

    }


    if( rank == 0 )
    {
      DBSetDir( m_dbBaseFilePtr, subDirectory.c_str() );
      DBtoc * const siloTOC = DBGetToc ( m_dbBaseFilePtr );
      bool stressFound = false;
      for( int ivar=0; ivar<siloTOC->nmultivar; ++ivar )
      {
        string varName = siloTOC->multivar_names[ivar];
        if( varName == "stress_11" )
        {
          stressFound = true;
        }
      }
      if( stressFound )
      {
        writeStressVarDefinition( subDirectory );
      }
      DBSetDir( m_dbBaseFilePtr, ".." );
    }

    DBSetDir( m_dbFilePtr, ".." );
  }
}

void SiloFile::writeMaterialVarDefinition( string const & subDir,
                                           string const & matDir,
                                           localIndex const matIndex,
                                           string const & fieldName )
{
  string const expressionName = matDir + "/" + fieldName;
  const char * expObjName = expressionName.c_str();
  const char * const names[1] = { expObjName };
  int const types[1] = { DB_VARTYPE_SCALAR };
  string const definition = "value_for_material(<" + subDir + "/" + fieldName + ">, " + std::to_string( matIndex ) + ")";
  const char * const defns[1] = { definition.c_str() };
  DBPutDefvars( m_dbBaseFilePtr,
                expObjName,
                1,
                names,
                types,
                defns,
                nullptr );
}

void SiloFile::clearEmptiesFromMultiObjects( int const cycleNum )
{

  int const size = MpiWrapper::commSize( MPI_COMM_GEOSX );
  int const rank = MpiWrapper::commRank( MPI_COMM_GEOSX );

  string sendbufferVars;
  string sendbufferMesh;
//  string sendbufferMaterials;

  if( rank != 0 )
  {
    for( string const & emptyObject : m_emptyVariables )
    {
      sendbufferVars += emptyObject + ' ';
    }

    for( string const & emptyObject : m_emptyMeshes )
    {
      sendbufferMesh += emptyObject + ' ';
    }

//    for( string_array::const_iterator emptyObject=m_emptyMaterials.begin() ;
//         emptyObject!=m_emptyMaterials.end() ; ++emptyObject )
//    {
//      sendbufferMesh += *emptyObject + ' ';
//    }

  }

  int sizeOfSendBufferVars = sendbufferVars.size();
  int sizeOfSendBufferMesh = sendbufferMesh.size();

  integer_array rcounts( size );
  integer_array displs( size );
  MpiWrapper::gather( &sizeOfSendBufferVars, 1, rcounts.data(), 1, 0, MPI_COMM_GEOSX );

  int sizeOfReceiveBuffer = 0;
  displs[0] = 0;
  for( int i=1; i<size; ++i )
  {
    displs[i] = displs[i-1]+rcounts[i-1];
    sizeOfReceiveBuffer += rcounts[i];
  }
  string receiveBufferVars( sizeOfReceiveBuffer, '\0' );

  MpiWrapper::gatherv ( &sendbufferVars[0],
                        sizeOfSendBufferVars,
                        &receiveBufferVars[0],
                        rcounts.data(),
                        displs.data(),
                        0,
                        MPI_COMM_GEOSX );


  MpiWrapper::gather( &sizeOfSendBufferMesh, 1, rcounts.data(), 1, 0, MPI_COMM_GEOSX );

  int sizeOfReceiveBufferMesh = 0;
  displs[0] = 0;
  for( int i=1; i<size; ++i )
  {
    displs[i] = displs[i-1]+rcounts[i-1];
    sizeOfReceiveBufferMesh += rcounts[i];
  }
  string receiveBufferMesh( sizeOfReceiveBufferMesh, '\0' );

  MpiWrapper::gatherv ( &sendbufferMesh[0],
                        sizeOfSendBufferMesh,
                        &receiveBufferMesh[0],
                        rcounts.data(),
                        displs.data(),
                        0,
                        MPI_COMM_GEOSX );



  if( rank== 0 )
  {
    std::istringstream iss( receiveBufferVars );
    copy( std::istream_iterator< string >( iss ),
          std::istream_iterator< string >(),
          std::back_inserter( m_emptyVariables ));

    std::istringstream issm( receiveBufferMesh );
    copy( std::istream_iterator< string >( issm ),
          std::istream_iterator< string >(),
          std::back_inserter( m_emptyMeshes ));
  }

  if( rank == 0 )
  {
    string baseFilePathAndName = m_siloDirectory + "/" + m_baseFileName;
    DBfile *siloFile = DBOpen( baseFilePathAndName.c_str(), DB_UNKNOWN, DB_APPEND );
    string empty( "EMPTY" );

    for( string const & emptyObject : m_emptyVariables )
    {
      size_t pathBegin = emptyObject.find_first_of( '/', 1 );
      size_t pathEnd = emptyObject.find_last_of( '/' );
      string domainString( emptyObject, 1, pathBegin );
      string pathToMultiObj( emptyObject, pathBegin, pathEnd - pathBegin );
      string varName( emptyObject, pathEnd + 1 );

      DBSetDir( siloFile, pathToMultiObj.c_str());

      DBmultivar * multiVar = DBGetMultivar( siloFile, varName.c_str());

      if( multiVar != nullptr )
      {
        array1d< const char * > newvarnames( multiVar->nvars );

        for( int i = 0; i < multiVar->nvars; ++i )
        {

          string path( multiVar->varnames[i] );
          if( path.find( domainString ) != string::npos )
          {
            newvarnames( i ) = empty.c_str();
          }
          else
          {
            newvarnames( i ) = multiVar->varnames[i];
          }
        }

        DBoptlist *optlist = DBMakeOptlist( 5 );
        DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNum));
        //      DBAddOption( optlist, DBOPT_DTIME,
        // const_cast<real64*>(&problemTime) );
        DBAddOption( optlist, DBOPT_REGION_PNAMES, multiVar->region_pnames );
        DBAddOption( optlist, DBOPT_TENSOR_RANK, &multiVar->tensor_rank );
        DBAddOption( optlist, DBOPT_MMESH_NAME, multiVar->mmesh_name );

        DBPutMultivar( siloFile, varName.c_str(), multiVar->nvars,
                       const_cast< char * * >(newvarnames.data()), multiVar->vartypes, optlist );
        DBFreeOptlist( optlist );
        DBFreeMultivar( multiVar );
      }
    }


    for( string const & emptyObject : m_emptyMeshes )
    {
      size_t pathBegin = emptyObject.find_first_of( '/', 1 );
      size_t pathEnd = emptyObject.find_last_of( '/' );
      string domainString( emptyObject, 1, pathBegin );
      string pathToMultiObj( emptyObject, pathBegin, pathEnd - pathBegin );
      string varName( emptyObject, pathEnd + 1 );

      if( !(pathToMultiObj.compare( "" )) )
      {
        pathToMultiObj = "/";
      }

      DBSetDir( siloFile, pathToMultiObj.c_str());

      DBmultimesh * multiMesh = DBGetMultimesh( siloFile, varName.c_str());

      if( multiMesh != nullptr )
      {
        array1d< const char * > newmeshnames( multiMesh->nblocks );

        for( int i = 0; i < multiMesh->nblocks; ++i )
        {

          string path( multiMesh->meshnames[i] );
          if( path.find( domainString ) != string::npos )
          {
            newmeshnames( i ) = empty.c_str();
          }
          else
          {
            newmeshnames( i ) = multiMesh->meshnames[i];
          }
        }

        DBoptlist *optlist = DBMakeOptlist( 2 );
        DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNum));
//        DBAddOption( optlist, DBOPT_DTIME, const_cast<real64*>(&problemTime) );

        DBPutMultimesh( siloFile, varName.c_str(), multiMesh->nblocks,
                        const_cast< char * * >(newmeshnames.data()), multiMesh->meshtypes, optlist );
        DBFreeOptlist( optlist );
        DBFreeMultimesh( multiMesh );
      }
    }
    DBClose( siloFile );
  }
  m_emptyVariables.clear();
  m_emptyMeshes.clear();

}

void SiloFile::writeGroupSilo( Group const & group,
                               string const & siloDirName,
                               string const & meshname,
                               int const centering,
                               int const cycleNum,
                               real64 const problemTime,
                               bool const isRestart,
                               const localIndex_array & mask )
{

  string subDirectory = siloDirName;
  string rootDirectory = "/" + siloDirName;

  {
    string shortsubdir( siloDirName );
    string::size_type pos = siloDirName.find_last_of( "//" );

    if( pos != shortsubdir.npos )
    {
      shortsubdir.erase( 0, pos+1 );
    }

    makeSubDirectory( shortsubdir, rootDirectory );
    DBSetDir( m_dbFilePtr, shortsubdir.c_str());
  }

  writeWrappersToSilo< real64 >( meshname,
                                 group.wrappers(),
                                 centering,
                                 cycleNum,
                                 problemTime,
                                 isRestart,
                                 rootDirectory,
                                 mask );



  DBSetDir( m_dbFilePtr, ".." );

}

void SiloFile::writeElementRegionSilo( ElementRegionBase const & elemRegion,
                                       string const & siloDirName,
                                       string const & meshName,
                                       int const cycleNum,
                                       real64 const problemTime,
                                       bool const isRestart )
{
  // TODO: This is a hack.
  conduit::Node conduitNode;
  dataRepository::Group fakeGroup( elemRegion.getName(), conduitNode );

  localIndex numElems = 0;
  std::vector< std::map< string, WrapperBase const * > > viewPointers;

  viewPointers.resize( elemRegion.numSubRegions() );
  elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
  {
    numElems += subRegion.size();

    for( auto const & wrapperIter : subRegion.wrappers() )
    {
      WrapperBase const * const wrapper = wrapperIter.second;

      if( wrapper->getPlotLevel() < m_plotLevel )
      {
        // the field name is the key to the map
        string const & fieldName = wrapper->getName();

        viewPointers[esr][fieldName] = wrapper;

        std::type_info const & typeID = wrapper->getTypeId();

        rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeID ),
                                        false,
                                        [&]( auto array, auto GEOSX_UNUSED_PARAM( Type ) )
        {
          typedef decltype( array ) arrayType;
          Wrapper< arrayType > const & sourceWrapper = dynamicCast< Wrapper< arrayType > const & >( *wrapper );
          traits::ViewTypeConst< arrayType > const sourceArray = sourceWrapper.reference();

          Wrapper< arrayType > & newWrapper = fakeGroup.registerWrapper< arrayType >( fieldName );
          newWrapper.setPlotLevel( PlotLevel::LEVEL_0 );
          arrayType & newarray = newWrapper.reference();
          newarray.resize( arrayType::NDIM, sourceArray.dims() );
        } );
      }
    }
  } );
  fakeGroup.resize( numElems );


  for( auto & wrapperIter : fakeGroup.wrappers() )
  {
    WrapperBase * const wrapper = wrapperIter.second;
    string const & fieldName = wrapper->getName();
    std::type_info const & typeID = wrapper->getTypeId();

    rtTypes::applyArrayTypeLambda2( rtTypes::typeID( typeID ),
                                    false,
                                    [&]( auto array, auto GEOSX_UNUSED_PARAM( scalar ) )
    {
      using arrayType =  decltype( array );
      using WrapperType = Wrapper< arrayType >;
      WrapperType & wrapperT = dynamicCast< Wrapper< arrayType > & >( *wrapper );
      arrayType & targetArray = wrapperT.reference();

      localIndex counter = 0;
      elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
        [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
      {
        // check if the field actually exists / plotted on the current subregion
        if( viewPointers[esr].count( fieldName ) > 0 )
        {
          WrapperType const & sourceWrapper = dynamic_cast< WrapperType const & >( *(viewPointers[esr].at( fieldName )) );
          auto const sourceArray = sourceWrapper.reference().toViewConst();

          localIndex const offset = counter * targetArray.strides()[ 0 ];
          GEOSX_ERROR_IF_GT( sourceArray.size(), targetArray.size() - offset );

          for( localIndex i = 0; i < sourceArray.size(); ++i )
          {
            targetArray.data()[ i + offset ] = sourceArray.data()[ i ];
          }

          counter += sourceArray.size( 0 );
        }
        else
        {
          counter += subRegion.size();
        }
      } );
    } );
  }

  writeGroupSilo( fakeGroup,
                  siloDirName,
                  meshName,
                  DB_ZONECENT,
                  cycleNum,
                  problemTime,
                  isRestart,
                  localIndex_array() );

  fakeGroup.getConduitNode().parent()->remove( fakeGroup.getName() );
}


void SiloFile::writeDomainPartition( DomainPartition const & domain,
                                     int const cycleNum,
                                     real64 const problemTime,
                                     bool const isRestart )
{

  MeshLevel const & mesh = domain.getMeshBody( 0 ).getMeshLevel( 0 );
  writeMeshLevel( mesh, cycleNum, problemTime, isRestart );

  if( isRestart )
  {
//    siloFile.DBWriteWrapper("m_globalDomainNumber",m_globalDomainNumber);
  }

}


void SiloFile::writeElementMesh( ElementRegionBase const & elementRegion,
                                 NodeManager const & nodeManager,
                                 string const & meshName,
                                 const localIndex numNodes,
                                 real64 * coords[3],
                                 globalIndex const * const globalNodeNum,
                                 char const * const ghostNodeFlag,
                                 int const cycleNumber,
                                 real64 const problemTime,
                                 bool & writeArbitraryPolygon )
{
  localIndex numElementShapes = 0;
  localIndex numElements = 0;

  elementRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
  {
    ++numElementShapes;
    numElements += subRegion.size();
  } );

  //if( numElements>0 )
  {
    array1d< localIndex * > meshConnectivity( numElementShapes );
    array1d< globalIndex * > globalElementNumbers( numElementShapes );
    array1d< integer > shapecnt( numElementShapes );
    array1d< integer > shapetype( numElementShapes );
    array1d< integer > shapesize( numElementShapes );
    array1d< char > ghostZoneFlag;


    std::vector< FixedOneToManyRelation > elementToNodeMap( numElementShapes );

    int count = 0;

    elementRegion.forElementSubRegions< CellElementSubRegion,
                                        FaceElementSubRegion,
                                        WellElementSubRegion >( [&]( auto const & elementSubRegion )
    {
      typename TYPEOFREF( elementSubRegion ) ::NodeMapType const & elemsToNodes = elementSubRegion.nodeList();

      // TODO HACK. this isn't correct for variable relations.
      elementToNodeMap[count].resize( elementSubRegion.size(), elementSubRegion.numNodesPerElement( 0 ) );

      arrayView1d< integer const > const & elemGhostRank = elementSubRegion.ghostRank();


      string const & elementType = elementSubRegion.getElementTypeString();
      std::vector< int > const & nodeOrdering = elementSubRegion.getVTKNodeOrdering();
      for( localIndex k = 0; k < elementSubRegion.size(); ++k )
      {
        integer numNodesPerElement = LvArray::integerConversion< int >( elementSubRegion.numNodesPerElement( k ));
        for( localIndex a = 0; a < numNodesPerElement; ++a )
        {
          elementToNodeMap[count]( k, a ) = elemsToNodes[k][nodeOrdering[a]];
        }

        if( elemGhostRank[k] >= 0 )
        {
          ghostZoneFlag.emplace_back( 1 );
        }
        else
        {
          ghostZoneFlag.emplace_back( 0 );
        }
      }

      meshConnectivity[count] = elementToNodeMap[count].data();

      //        globalElementNumbers[count] = elementRegion.localToGlobalMap().data();
      shapecnt[count] = static_cast< int >(elementSubRegion.size());


      if( !elementType.compare( 0, 4, "C3D8" ) )
      {
        shapetype[count] = DB_ZONETYPE_HEX;
      }
      else if( !elementType.compare( 0, 4, "C3D4" ) )
      {
        shapetype[count] = DB_ZONETYPE_TET;
      }
      else if( !elementType.compare( 0, 4, "C3D6" ) )
      {
        shapetype[count] = DB_ZONETYPE_PRISM;
        writeArbitraryPolygon = true;
      }
      else if( !elementType.compare( 0, 4, "C3D5" ) )
      {
        shapetype[count] = DB_ZONETYPE_PYRAMID;
        writeArbitraryPolygon = true;
      }
      else if( !elementType.compare( 0, 4, "BEAM" ) )
      {
        shapetype[count] = DB_ZONETYPE_BEAM;
      }
      shapesize[count] = LvArray::integerConversion< int >( elementSubRegion.numNodesPerElement( 0 ) );
      ++count;
    } );

    string_array
      regionSolidMaterialList = elementRegion.getConstitutiveNames< constitutive::RockBase >();
//    string_array const
//    regionSolidMaterialList2 = elementRegion.getConstitutiveNames< constitutive::PoreVolumeCompressibleSolid >();
//
//    for( string const & entry : regionSolidMaterialList2 )
//    {
//      regionSolidMaterialList.emplace_back( entry );
//    }
    localIndex const numSolids = regionSolidMaterialList.size();

    string_array regionFluidMaterialList = elementRegion.getConstitutiveNames< constitutive::SingleFluidBase >();
    string_array regionMultiPhaseFluidList = elementRegion.getConstitutiveNames< constitutive::MultiFluidBase >();

    for( string const & matName : regionMultiPhaseFluidList )
    {
      regionFluidMaterialList.emplace_back( matName );
    }
    localIndex const numFluids = regionFluidMaterialList.size();

    string_array
      fractureContactMaterialList = elementRegion.getConstitutiveNames< constitutive::ContactRelationBase >();

    localIndex const numContacts = fractureContactMaterialList.size();

    string_array
      nullModelMaterialList = elementRegion.getConstitutiveNames< constitutive::NullModel >();

    localIndex const numNullModels = nullModelMaterialList.size();

    if( numSolids + numFluids + numContacts + numNullModels > 0 )
    {
      writeMeshObject( meshName,
                       numNodes,
                       coords,
                       globalNodeNum,
                       ghostNodeFlag,
                       ghostZoneFlag.data(),
                       LvArray::integerConversion< int >( numElementShapes ),
                       shapecnt.data(),
                       meshConnectivity.data(),
                       nullptr /*globalElementNumbers.data()*/,
                       shapetype.data(),
                       shapesize.data(),
                       cycleNumber,
                       problemTime );

      writeGroupSilo( nodeManager,
                      meshName + "_NodalFields",
                      meshName,
                      DB_NODECENT,
                      cycleNumber,
                      problemTime,
                      false,
                      localIndex_array() );

      writeElementRegionSilo( elementRegion,
                              meshName + "_ElementFields",
                              meshName,
                              cycleNumber,
                              problemTime,
                              false );
    }

    if( numSolids > 0 )
    {
      string const solidMeshName = meshName + "_Solid";
      writeMeshObject( solidMeshName,
                       numNodes,
                       coords,
                       globalNodeNum,
                       ghostNodeFlag,
                       ghostZoneFlag.data(),
                       LvArray::integerConversion< int >( numElementShapes ),
                       shapecnt.data(),
                       meshConnectivity.data(),
                       nullptr /*globalElementNumbers.data()*/,
                       shapetype.data(),
                       shapesize.data(),
                       cycleNumber,
                       problemTime );

      writeMaterialMapsFullStorage( elementRegion,
                                    solidMeshName,
                                    regionSolidMaterialList,
                                    cycleNumber,
                                    problemTime );
    }

    if( numFluids > 0 )
    {
      string const fluidMeshName = meshName + "_Fluid";
      writeMeshObject( fluidMeshName,
                       numNodes,
                       coords,
                       globalNodeNum,
                       ghostNodeFlag,
                       ghostZoneFlag.data(),
                       LvArray::integerConversion< int >( numElementShapes ),
                       shapecnt.data(),
                       meshConnectivity.data(),
                       nullptr /*globalElementNumbers.data()*/,
                       shapetype.data(),
                       shapesize.data(),
                       cycleNumber,
                       problemTime );

      writeMaterialMapsFullStorage( elementRegion,
                                    fluidMeshName,
                                    regionFluidMaterialList,
                                    cycleNumber,
                                    problemTime );
    }

    if( numContacts > 0 )
    {
      string const contactMeshName = meshName + "_Contact";
      writeMeshObject( contactMeshName,
                       numNodes,
                       coords,
                       globalNodeNum,
                       ghostNodeFlag,
                       ghostZoneFlag.data(),
                       LvArray::integerConversion< int >( numElementShapes ),
                       shapecnt.data(),
                       meshConnectivity.data(),
                       nullptr /*globalElementNumbers.data()*/,
                       shapetype.data(),
                       shapesize.data(),
                       cycleNumber,
                       problemTime );

      writeMaterialMapsFullStorage( elementRegion,
                                    contactMeshName,
                                    fractureContactMaterialList,
                                    cycleNumber,
                                    problemTime );
    }
  }
}

void SiloFile::writeMeshLevel( MeshLevel const & meshLevel,
                               int const cycleNum,
                               real64 const problemTime,
                               bool const isRestart )
{

  NodeManager const & nodeManager = meshLevel.getNodeManager();
  localIndex const numNodes = nodeManager.size();


  string const ghostNodeName = "ghostNodeFlag";
  string const ghostZoneName = "ghostZoneFlag";

  arrayView1d< integer const > const & nodeGhostRank = nodeManager.ghostRank();
  array1d< char > ghostNodeFlag( nodeGhostRank.size() );
  array1d< char > ghostZoneFlag;

  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & referencePosition = nodeManager.referencePosition();

  bool writeArbitraryPolygon( false );
  string const meshName( "MeshLevel" );

  //set the nodal coordinate data structure
  real64 * coords[3];
  array1d< real64 > xcoords( numNodes );
  array1d< real64 > ycoords( numNodes );
  array1d< real64 > zcoords( numNodes );
  for( localIndex a = 0; a < numNodes; ++a )
  {
    xcoords[a] = referencePosition( a, 0 );
    ycoords[a] = referencePosition( a, 1 );
    zcoords[a] = referencePosition( a, 2 );

    if( nodeGhostRank[a] >=0 )
    {
      ghostNodeFlag[a] = 1;
    }
  }

  if( nodeManager.hasWrapper( keys::TotalDisplacement ) )
  {
    arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const & totalDisplacement = nodeManager.totalDisplacement();
    for( localIndex a = 0; a < numNodes; ++a )
    {
      xcoords[a] += totalDisplacement( a, 0 );
      ycoords[a] += totalDisplacement( a, 1 );
      zcoords[a] += totalDisplacement( a, 2 );
    }
  }

  coords[0] = xcoords.data();
  coords[1] = ycoords.data();
  coords[2] = zcoords.data();



  ElementRegionManager const & elementManager = meshLevel.getElemManager();
  elementManager.forElementRegions( [&]( ElementRegionBase const & elemRegion )
  {
    string const & regionName = elemRegion.getName();

    writeElementMesh( elemRegion,
                      nodeManager,
                      regionName,
                      numNodes,
                      coords,
                      nodeManager.localToGlobalMap().data(),
                      ghostNodeFlag.data(),
                      cycleNum,
                      problemTime,
                      writeArbitraryPolygon );
  } );


  if( m_writeFaceMesh )
  {
    FaceManager const & faceManager = meshLevel.getFaceManager();
    localIndex const numFaces = faceManager.size();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    // face mesh
    const string facemeshName( "face_mesh" );

    if( writeArbitraryPolygon )
    {
      const int numFaceTypes = 1;
      int dbZoneType = DB_ZONETYPE_POLYGON;
      // See a discussion of silo's arbitrary polygon implementation at
      // https://visitbugs.ornl.gov/projects/7/wiki/Arbitrary_Polygons_and_Polyhedra_in_Silo
      // It is not documented in silo manual.
      array1d< localIndex * > faceConnectivity( numFaceTypes );
      array1d< globalIndex const * > globalFaceNumbers( numFaceTypes );
      std::vector< int > fshapecnt( numFaceTypes );
      std::vector< int > fshapetype( numFaceTypes );
      std::vector< int > fshapesize( numFaceTypes );

      array1d< array1d< localIndex > > faceToNodeMapCopy( numFaceTypes );
      {
        for( localIndex k = 0; k < numFaces; ++k )
        {
          faceToNodeMapCopy[0].emplace_back( faceToNodeMap.sizeOfArray( k ));
          for( localIndex const a : faceToNodeMap[ k ] )
          {
            faceToNodeMapCopy[0].emplace_back( a );
          }
        }

        faceConnectivity[0] = faceToNodeMapCopy[0].data();

        globalFaceNumbers[0] = faceManager.localToGlobalMap().data();
        fshapecnt[0] = numFaces;
        fshapetype[0] = dbZoneType;
        fshapesize[0] = 0;
      }
      int lnodelist = faceToNodeMapCopy[0].size();

      writePolygonMeshObject( facemeshName, numNodes, coords,
                              nodeManager.localToGlobalMap().data(), numFaceTypes,
                              fshapecnt.data(), faceConnectivity.data(), globalFaceNumbers.data(),
                              nullptr, fshapetype.data(), fshapesize.data(), cycleNum, problemTime, lnodelist );


    }
    else  //The old way
    {
      const int numFaceTypes = 1;
      int numNodesPerFace = faceToNodeMap.sizeOfArray( 0 ); // TODO assumes all faces have same number of nodes
      int dbZoneType = DB_ZONETYPE_POLYGON;
      if( numNodesPerFace == 3 )
      {
        dbZoneType = DB_ZONETYPE_TRIANGLE;
      }
      else if( numNodesPerFace == 4 )
      {
        dbZoneType = DB_ZONETYPE_QUAD;
      }
      else if( numNodesPerFace == 2 )
      {
        dbZoneType = DB_ZONETYPE_BEAM;
      }

      array1d< localIndex * > faceConnectivity( numFaceTypes );
      array1d< globalIndex const * > globalFaceNumbers( numFaceTypes );
      std::vector< int > fshapecnt( numFaceTypes );
      std::vector< int > fshapetype( numFaceTypes );
      std::vector< int > fshapesize( numFaceTypes );

      array1d< array2d< localIndex > > faceToNodeMapCopy( numFaceTypes );


      for( int faceType = 0; faceType < numFaceTypes; ++faceType )
      {
        faceToNodeMapCopy[faceType].resize( numFaces, numNodesPerFace );

        for( localIndex k = 0; k < numFaces; ++k )
        {
          for( int a = 0; a < numNodesPerFace; ++a )
          {
            faceToNodeMapCopy[faceType][k][a] = faceToNodeMap( k, a );
          }
        }

        faceConnectivity[faceType] = faceToNodeMapCopy[faceType].data();

        globalFaceNumbers[faceType] = faceManager.localToGlobalMap().data();
        fshapecnt[faceType] = numFaces;
        fshapetype[faceType] = dbZoneType;
        fshapesize[faceType] = numNodesPerFace;
      }

      writeMeshObject( facemeshName,
                       numNodes,
                       coords,
                       nodeManager.localToGlobalMap().data(),
                       nullptr,
                       nullptr,
                       numFaceTypes,
                       fshapecnt.data(),
                       faceConnectivity.data(),
                       globalFaceNumbers.data(),
                       fshapetype.data(),
                       fshapesize.data(),
                       cycleNum,
                       problemTime );
    }

    writeGroupSilo( faceManager,
                    "FaceFields",
                    facemeshName,
                    DB_ZONECENT,
                    cycleNum,
                    problemTime,
                    isRestart,
                    localIndex_array());

  }


  if( m_writeEdgeMesh )
  {
    // write edges
    FaceManager const & faceManager = meshLevel.getFaceManager();
    ArrayOfArraysView< localIndex const > const & faceToNodeMap = faceManager.nodeList().toViewConst();

    EdgeManager const & edgeManager = meshLevel.getEdgeManager();
    localIndex const numEdges = edgeManager.size();

    const string edgeMeshName( "edge_mesh" );

    const int numEdgeTypes = 1;
    const int numNodesPerEdge = 2;
    int dbZoneType = DB_ZONETYPE_BEAM;

    array1d< localIndex * > edgeConnectivity( numEdgeTypes );
    array1d< globalIndex const * > globalEdgeNumbers( numEdgeTypes );
    std::vector< int > eshapecnt( numEdgeTypes );
    std::vector< int > eshapetype( numEdgeTypes );
    std::vector< int > eshapesize( numEdgeTypes );

    array1d< array2d< localIndex > > edgeToNodeMap( numEdgeTypes );


    for( int edgeType = 0; edgeType < numEdgeTypes; ++edgeType )
    {
      edgeToNodeMap[edgeType].resize( numEdges, numNodesPerEdge );

      for( localIndex k = 0; k < numEdges; ++k )
      {
        for( int a = 0; a < numNodesPerEdge; ++a )
        {
          if( faceToNodeMap.sizeOfArray( 0 ) == 2 && a > 0 )
          {
            edgeToNodeMap[edgeType][k][a] = edgeManager.nodeList()[k][0];
          }
          else
          {
            edgeToNodeMap[edgeType][k][a] = edgeManager.nodeList()[k][a];
          }
        }
      }

      edgeConnectivity[edgeType] = edgeToNodeMap[edgeType].data();

      globalEdgeNumbers[edgeType] = edgeManager.localToGlobalMap().data();
      eshapecnt[edgeType] = numEdges;
      eshapetype[edgeType] = dbZoneType;
      eshapesize[edgeType] = numNodesPerEdge;
    }

    writeMeshObject( edgeMeshName,
                     numNodes,
                     coords,
                     nodeManager.localToGlobalMap().data(),
                     nullptr,
                     nullptr,
                     numEdgeTypes,
                     eshapecnt.data(),
                     edgeConnectivity.data(),
                     globalEdgeNumbers.data(),
                     eshapetype.data(),
                     eshapesize.data(),
                     cycleNum,
                     problemTime );

    writeGroupSilo( edgeManager,
                    "EdgeFields",
                    edgeMeshName,
                    DB_ZONECENT,
                    cycleNum,
                    problemTime,
                    isRestart,
                    localIndex_array());
  }

}

// Arbitrary polygon. Have to deal with this separately
void SiloFile::writePolygonMeshObject( const string & meshName,
                                       const localIndex nnodes,
                                       real64 * coords[3],
                                       const globalIndex *,
                                       const int numRegions,
                                       const int * shapecnt,
                                       const localIndex * const * const meshConnectivity,
                                       const globalIndex * const * const globalElementNum,
                                       const int * const * const,
                                       const int * const shapetype,
                                       const int * const shapesize,
                                       const int cycleNumber,
                                       const real64 problemTime,
                                       const int lnodelist )
{

  const DBdatatype datatype = DB_DOUBLE;


//  DBfacelist* facelist;
//  string facelistName;
//  facelistName = meshName + "_facelist";

  DBoptlist * optlist = DBMakeOptlist( 4 );
//  DBAddOption(optlist, DBOPT_NODENUM, const_cast<globalIndex*> (globalNodeNum));
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

  int numTotZones = shapecnt[0];
  if( numTotZones == 0 )
  {
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );
    string emptyObject = pwd;
    emptyObject += "/" + meshName;
    m_emptyMeshes.emplace_back( emptyObject );
  }
  else
  {
    string zonelistName;
    zonelistName = meshName + "_zonelist";


    DBPutUcdmesh( m_dbFilePtr, meshName.c_str(), 3, nullptr, (float * *) coords, nnodes, numTotZones,
                  zonelistName.c_str(), nullptr, datatype, optlist );

    DBClearOptlist( optlist );

    std::vector< int > nodelist( lnodelist );
    std::vector< globalIndex > globalZoneNumber( lnodelist );

    int elemCount = 0;
    for( int j = 0; j < lnodelist; ++j )
    {
      nodelist[j] = meshConnectivity[0][j];
    }

    {
      if( globalElementNum != nullptr )
      {
        for( int j = 0; j < shapecnt[0]; ++j )
        {
          globalZoneNumber[elemCount++] = globalElementNum[0][j];
        }
        // write zonelist
        //      DBAddOption(optlist, DBOPT_ZONENUM, const_cast<globalIndex*> (globalZoneNumber.data()));
//        if (type_name<globalIndex>::name() == type_name<long long>::name())
//          DBAddOption(optlist, DBOPT_LLONGNZNUM, const_cast<int*> (&one));
      }
    }



    int hi_offset = 0;

    DBPutZonelist2( m_dbFilePtr,
                    zonelistName.c_str(),
                    numTotZones,
                    3,
                    nodelist.data(),
                    lnodelist,
                    0,
                    0,
                    hi_offset,
                    const_cast< int * >(shapetype),
                    const_cast< int * >(shapesize),
                    const_cast< int * >(shapecnt),
                    numRegions,
                    optlist );

    DBClearOptlist( optlist );

  }

  // write multimesh object
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
#endif
  if( rank == 0 )
  {
    DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
    DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

    writeMultiXXXX( DB_UCDMESH, DBPutMultimesh, 0, meshName, cycleNumber, "/", optlist );
  }

  DBFreeOptlist( optlist );
}


template< typename TYPE, int NDIM, int USD >
int getTensorRankOfArray( ArrayView< TYPE const, NDIM, USD > const & field )
{
  int rval = 1;
  if( NDIM==2 )
  {
    localIndex const lastDim = field.size( 1 );
    if( lastDim == 3 )
    {
      rval = DB_VARTYPE_VECTOR;
    }
    else if( lastDim == 6 )
    {
      rval = DB_VARTYPE_SYMTENSOR;
    }
    else if( lastDim == 9 )
    {
      rval = DB_VARTYPE_TENSOR;
    }
  }
  return rval;
}

template< typename OUTPUTTYPE >
void SiloFile::writeWrappersToSilo( string const & meshname,
                                    const dataRepository::Group::wrapperMap & wrappers,
                                    int const centering,
                                    int const cycleNum,
                                    real64 const problemTime,
                                    bool const GEOSX_UNUSED_PARAM( isRestart ),
                                    string const & multiRoot,
                                    const localIndex_array & GEOSX_UNUSED_PARAM( mask ) )
{

  // iterate over all entries in the member map
  for( auto const & wrapperIter : wrappers )
  {
    auto const & wrapper = wrapperIter.second;

    if( wrapper->getPlotLevel() <= m_plotLevel )
    {
      // the field name is the key to the map
      string const & fieldName = wrapper->getName();

      std::type_info const & typeID = wrapper->getTypeId();

      // TODO This is wrong. problem with uniqueness
      if( typeID==typeid(array1d< real64 >) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< array1d< real64 > > const & >( *wrapper );
        this->writeDataField< real64 >( meshname.c_str(), fieldName,
                                        wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(array2d< real64 >) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< array2d< real64 > > const & >( *wrapper );

        arrayView2d< real64 const > const & array = wrapperT.reference();
        this->writeDataField< real64 >( meshname.c_str(),
                                        fieldName,
                                        array,
//                                      getTensorRankOfArray(array),
                                        centering,
                                        cycleNum,
                                        problemTime,
                                        multiRoot );

        writeVectorVarDefinition( fieldName, multiRoot );
      }
      if( typeID==typeid(array2d< real64, RAJA::PERM_JI >) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< array2d< real64, RAJA::PERM_JI > > const & >( *wrapper );

        arrayView2d< real64 const, LvArray::typeManipulation::getStrideOneDimension( RAJA::PERM_JI {} ) > const &
        array = wrapperT.reference();
        this->writeDataField< real64 >( meshname.c_str(),
                                        fieldName,
                                        array,
//                                      getTensorRankOfArray(array),
                                        centering,
                                        cycleNum,
                                        problemTime,
                                        multiRoot );
        writeVectorVarDefinition( fieldName, multiRoot );
      }
      if( typeID==typeid(array3d< real64 >) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< array3d< real64 > > const & >( *wrapper );
        this->writeDataField< real64 >( meshname.c_str(), fieldName,
                                        wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(r1_array) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< r1_array > const & >( *wrapper );
        this->writeDataField< real64 >( meshname.c_str(), fieldName,
                                        wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(integer_array) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< integer_array > const & >( *wrapper );
        this->writeDataField< integer >( meshname.c_str(), fieldName,
                                         wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(localIndex_array) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< localIndex_array > const & >( *wrapper );
        this->writeDataField< localIndex >( meshname.c_str(), fieldName,
                                            wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
      if( typeID==typeid(globalIndex_array) )
      {
        auto const & wrapperT = dynamic_cast< dataRepository::Wrapper< globalIndex_array > const & >( *wrapper );
        this->writeDataField< globalIndex >( meshname.c_str(), fieldName,
                                             wrapperT.reference(), centering, cycleNum, problemTime, multiRoot );
      }
    }
  }
}


template< typename CBF >
void SiloFile::writeMultiXXXX( const DBObjectType type,
                               CBF DBPutMultiCB,
                               int const centering,
                               string const name,
                               const int,
                               string const & multiRoot,
                               const DBoptlist * optlist )
{
  (void)centering;

  int size = 1;
#ifdef GEOSX_USE_MPI
  MPI_Comm_size( MPI_COMM_GEOSX, &size );
#endif

  string_array vBlockNames( size );
  array1d< char * > BlockNames( size );
  array1d< int > blockTypes( size );
  char tempBuffer[1024];
  char currentDirectory[256];

  DBGetDir( m_dbBaseFilePtr, currentDirectory );
  DBSetDir( m_dbBaseFilePtr, multiRoot.c_str());


  string multiRootString( multiRoot );
  if( !(multiRootString.compare( "/" )) )
  {
    multiRootString.clear();
  }

  for( int i = 0; i < size; ++i )
  {

    sprintf( tempBuffer,
             "%s%s%s.%03d:/domain_%05d%s/%s",
             m_siloDataSubDirectory.c_str(),
             "/",
             m_baseFileName.c_str(),
             groupRank( i ),
             i,
             multiRootString.c_str(),
             name.c_str());

    vBlockNames[i] = tempBuffer;
    BlockNames[i] = const_cast< char * >( vBlockNames[i].c_str() );
    blockTypes[i] = type;
  }

  string multiName = name;
  DBPutMultiCB( m_dbBaseFilePtr, multiName.c_str(), size, BlockNames.data(), blockTypes.data(),
                const_cast< DBoptlist * >(optlist));

  DBSetDir( m_dbBaseFilePtr, currentDirectory );

}


/**
 *
 * @param meshName
 * @param fieldName
 * @param field
 * @param centering
 * @param cycleNumber
 * @param problemTime
 */
template< typename OUTTYPE, typename TYPE >
void SiloFile::writeDataField( string const & meshName,
                               string const & fieldName,
                               arrayView1d< TYPE const > const & field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot )
{
  int const nvars = SiloFileUtilities::GetNumberOfVariablesInField< TYPE >();
  int nels = LvArray::integerConversion< int >( field.size());

  int const meshType = getMeshType( meshName );


  DBoptlist *optlist = DBMakeOptlist( 5 );
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

  array1d< char const * > varnames( nvars );
  array1d< void * > vars( nvars );


  string_array varnamestring( nvars );
  array1d< array1d< OUTTYPE > > castedField( nvars );


  for( int i = 0; i < nvars; ++i )
  {
    if( std::is_same< OUTTYPE, TYPE >::value )
    {
      vars[i] = const_cast< void * >(static_cast< void const * >(field.data()+i));
    }
    else
    {
      castedField[i].resize( nels );
      vars[i] = static_cast< void * >( (castedField[i]).data() );
      for( int k = 0; k < nels; ++k )
      {
        castedField[i][k] = SiloFileUtilities::CastField< OUTTYPE >( field[k], i );
      }
    }
  }

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if( nels == 0 )
  {
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );
    string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.emplace_back( emptyObject );
  }
  else
  {

    SiloFileUtilities::SetVariableNames< TYPE >( fieldName, varnamestring, varnames.data() );


    int err = -2;
//    if( meshType == DB_UCDMESH )
//    {
    err = DBPutUcdvar( m_dbFilePtr,
                       fieldName.c_str(),
                       meshName.c_str(),
                       nvars,
                       varnames.data(),
                       reinterpret_cast< float * * >(vars.data()),
                       nels,
                       nullptr,
                       0,
                       SiloFileUtilities::DB_TYPE< OUTTYPE >(),
                       centering,
                       optlist );
//    }
//    else if( meshType == DB_POINTMESH )
//    {
//      err = DBPutPointvar( m_dbFilePtr,
//                           fieldName.c_str(),
//                           meshName.c_str(),
//                           nvars,
//                           reinterpret_cast<float**>(vars.data()),
//                           nels,
//                           SiloFileUtilities::DB_TYPE<OUTTYPE>(),
//                           optlist);
//    }
    if( err < 0 )
    {
      if( err < -1 )
      {
        GEOSX_ERROR( "unhandled case in SiloFile::WriteDataField A\n" );
      }
      else
      {
        GEOSX_ERROR( "unhandled failure in adding variable during SiloFile::WriteDataField\n" );
      }
    }
  }

  // write multimesh object
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  if( rank == 0 )
  {
    int tensorRank = SiloFileUtilities::GetTensorRank< TYPE >();
    DBAddOption( optlist, DBOPT_TENSOR_RANK, const_cast< int * >(&tensorRank));
    DBAddOption( optlist, DBOPT_MMESH_NAME, const_cast< char * >(meshName.c_str()));

    DBObjectType vartype = DB_INVALID_OBJECT;

    if( meshType == DB_UCDMESH )
    {
      vartype = DB_UCDVAR;
    }
    else if( meshType == DB_POINTMESH )
    {
      vartype = DB_POINTVAR;
    }
    else if( meshType == DB_QUADCURV )
    {
      vartype = DB_QUADVAR;
    }
    else
    {
      vartype = DB_UCDVAR;
    }


    writeMultiXXXX( vartype, DBPutMultivar, centering, fieldName.c_str(), cycleNumber, multiRoot,
                    optlist );
  }

  DBFreeOptlist( optlist );

}

template< typename OUTTYPE, typename TYPE, int USD >
void SiloFile::writeDataField( string const & meshName,
                               string const & fieldName,
                               arrayView2d< TYPE const, USD > const & field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot )
{
  int const primaryDimIndex = 0;
  int const secondaryDimIndex = 1;

  localIndex const npts = field.size( primaryDimIndex );
  localIndex const nvar = field.size( secondaryDimIndex );

  array1d< TYPE > data( npts );
  localIndex indices[2];

  for( localIndex ivar = 0; ivar < nvar; ++ivar )
  {
    indices[secondaryDimIndex] = ivar;
    // make a copy of the ivar'th slice of data
    for( localIndex ip = 0; ip < npts; ++ip )
    {
      indices[primaryDimIndex] = ip;
      data[ip] = field( indices[0], indices[1] );
    }
    writeDataField< OUTTYPE >( meshName,
                               fieldName + "_" + std::to_string( ivar ),
                               data.toViewConst(),
                               centering,
                               cycleNumber,
                               problemTime,
                               multiRoot );
  }
}

template< typename OUTTYPE, typename TYPE, int USD >
void SiloFile::writeDataField( string const & meshName,
                               string const & fieldName,
                               arrayView3d< TYPE const, USD > const & field,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot )
{
  int const primaryDimIndex = 0;
  int const secondaryDimIndex1 = 1;
  int const secondaryDimIndex2 = 2;

  localIndex const npts  = field.size( primaryDimIndex );
  localIndex const nvar1 = field.size( secondaryDimIndex1 );
  localIndex const nvar2 = field.size( secondaryDimIndex2 );

  array1d< TYPE > data( npts );
  localIndex indices[3];

  for( localIndex ivar = 0; ivar < nvar1; ++ivar )
  {
    indices[secondaryDimIndex1] = ivar;
    for( localIndex jvar = 0; jvar < nvar2; ++jvar )
    {
      indices[secondaryDimIndex2] = jvar;
      // make a copy of the ivar/jvar'th slice of data
      for( localIndex ip = 0; ip < npts; ++ip )
      {
        indices[primaryDimIndex] = ip;
        data[ip] = field( indices[0], indices[1], indices[2] );
      }
      writeDataField< OUTTYPE >( meshName,
                                 fieldName + "_" + std::to_string( ivar ) + "_" + std::to_string( jvar ),
                                 data.toViewConst(),
                                 centering,
                                 cycleNumber,
                                 problemTime,
                                 multiRoot );
    }
  }
}


template< typename OUTTYPE, typename TYPE, int USD >
void createSiloVarArrays( ArrayView< TYPE const, 1, USD > const & field,
                          array1d< array1d< OUTTYPE > > & castedField,
                          array1d< void * > & vars )
{
  castedField[0].resize( field.size());
  vars[0] = static_cast< void * >( (castedField[0]).data() );
  for( localIndex k=0; k<field.size(); ++k )
  {
    castedField[0][k] = static_cast< OUTTYPE >(field( k ));
  }
}

template< typename OUTTYPE, typename TYPE, int USD >
void createSiloVarArrays( ArrayView< TYPE const, 2, USD > const & field,
                          array1d< array1d< OUTTYPE > > & castedField,
                          array1d< void * > & vars )
{
  for( localIndex i=0; i<field.size( 1 ); ++i )
  {
    castedField[i].resize( field.size( 0 ));
    vars[i] = static_cast< void * >( (castedField[i]).data() );
    for( localIndex k=0; k<field.size( 0 ); ++k )
    {
      castedField[i][k] = static_cast< OUTTYPE >(field( k, i ));
    }
  }
}

template< typename OUTTYPE, typename TYPE, int USD >
void createSiloVarArrays( ArrayView< TYPE const, 3, USD > const & field,
                          array1d< array1d< OUTTYPE > > & castedField,
                          array1d< void * > & vars )
{
  for( localIndex i=0; i<field.size( 1 ); ++i )
  {
    for( localIndex j=0; j<field.size( 2 ); ++j )
    {
      castedField[i*(field.size( 2 ))+j].resize( field.size( 0 ));
      vars[i*(field.size( 2 ))+j] = static_cast< void * >( (castedField[i*(field.size( 2 ))+j]).data() );
      for( localIndex k=0; k<field.size( 0 ); ++k )
      {
        castedField[i*(field.size( 2 ))+j][k] = static_cast< OUTTYPE >(field( k, i, j ));
      }
    }
  }
}


template< typename OUTTYPE, typename TYPE, int NDIM, int USD >
void SiloFile::writeDataField( string const & meshName,
                               string const & fieldName,
                               ArrayView< TYPE const, NDIM, USD > const & field,
                               int const,
                               int const centering,
                               int const cycleNumber,
                               real64 const problemTime,
                               string const & multiRoot )
{
  int nvars = 1;
  for( int i=1; i<NDIM; ++i )
  {
    nvars *= field.size( i );
  }
  int const nels = LvArray::integerConversion< int >( field.size( 0 ));

  int const meshType = getMeshType( meshName );

  DBoptlist *optlist = DBMakeOptlist( 5 );
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));

  array1d< char const * > varnames( nvars );
  array1d< void * > vars( nvars );

  string_array varnamestring( nvars );
  array1d< array1d< OUTTYPE > > castedField( nvars );

  createSiloVarArrays( field, castedField, vars );

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if( nels == 0 )
  {
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );
    string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.emplace_back( emptyObject );
  }
  else
  {


    varnamestring.resize( nvars );
    for( int i=0; i<nvars; ++i )
    {
      varnamestring[i] = fieldName + "_" + std::to_string( i );
      varnames[i] = const_cast< char * >( varnamestring[i].c_str() );
    }

    SiloFileUtilities::SetVariableNames< TYPE >( fieldName, varnamestring, varnames.data() );


    int err = -2;
//    if( meshType == DB_UCDMESH )
//    {
    err = DBPutUcdvar( m_dbFilePtr,
                       fieldName.c_str(),
                       meshName.c_str(),
                       nvars,
                       varnames.data(),
                       reinterpret_cast< float * * >(vars.data()),
                       nels,
                       nullptr,
                       0,
                       SiloFileUtilities::DB_TYPE< OUTTYPE >(),
                       centering,
                       optlist );
//    }
//    else if( meshType == DB_POINTMESH )
//    {
//      err = DBPutPointvar( m_dbFilePtr,
//                           fieldName.c_str(),
//                           meshName.c_str(),
//                           nvars,
//                           reinterpret_cast<float**>(vars.data()),
//                           nels,
//                           SiloFileUtilities::DB_TYPE<OUTTYPE>(),
//                           optlist);
//    }
    if( err < 0 )
    {
      if( err < -1 )
      {
        GEOSX_ERROR( "unhandled case in SiloFile::WriteDataField A\n" );
      }
      else
      {
        GEOSX_ERROR( "unhandled failure in adding variable during SiloFile::WriteDataField\n" );
      }
    }
  }

  // write multimesh object
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  if( rank == 0 )
  {
//    int tensorRank = siloTensorRank;
//    DBAddOption(optlist, DBOPT_TENSOR_RANK, const_cast<int*> (&tensorRank));
    DBAddOption( optlist, DBOPT_MMESH_NAME, const_cast< char * >(meshName.c_str()));

    DBObjectType vartype = DB_INVALID_OBJECT;

    if( meshType == DB_UCDMESH )
    {
      vartype = DB_UCDVAR;
    }
    else if( meshType == DB_POINTMESH )
    {
      vartype = DB_POINTVAR;
    }
    else if( meshType == DB_QUADCURV )
    {
      vartype = DB_QUADVAR;
    }
    else
    {
      vartype = DB_UCDVAR;
    }


    writeMultiXXXX( vartype, DBPutMultivar, centering, fieldName.c_str(), cycleNumber, multiRoot,
                    optlist );
  }

  DBFreeOptlist( optlist );

}


template< typename OUTTYPE, typename TYPE >
void SiloFile::writeMaterialDataField2d( string const & meshName,
                                         string const & fieldName,
                                         ElementRegionBase const & elemRegion,
                                         int const centering,
                                         int const cycleNumber,
                                         real64 const problemTime,
                                         string const & multiRoot,
                                         string_array const & materialNames )
{
  array1d< array1d< array2d< TYPE > > > fieldData( elemRegion.numSubRegions());
  array1d< array1d< arrayView2d< TYPE const > > > field( elemRegion.numSubRegions());


  elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
  {
    localIndex const nmat = materialNames.size();
    field[esr].resize( nmat );
    fieldData[esr].resize( nmat );
    for( int matIndex=0; matIndex<nmat; ++matIndex )
    {
      Group const & constitutiveModel = subRegion.getConstitutiveModel( materialNames[matIndex] );

      if( constitutiveModel.hasWrapper( fieldName ) )
      {
        arrayView2d< TYPE const > const & fieldView = constitutiveModel.getReference< array2d< TYPE > >( fieldName );

        fieldData[esr][matIndex].resize( fieldView.size( 0 ), 1 );
        field[esr][matIndex] = fieldData[esr][matIndex].toViewConst();
        for( localIndex k = 0; k<fieldView.size( 0 ); ++k )
        {
          fieldData[esr][matIndex]( k, 0 ) = 0;
          for( localIndex q=0; q<fieldView.size( 1 ); ++q )
          {
            fieldData[esr][matIndex]( k, 0 ) += fieldView( k, q );
          }
          fieldData[esr][matIndex]( k, 0 ) *= 1.0 / fieldView.size( 1 );
        }
      }
    }
  } );

  writeMaterialDataField< OUTTYPE, TYPE >( meshName,
                                           fieldName,
                                           field,
                                           elemRegion,
                                           centering,
                                           cycleNumber,
                                           problemTime,
                                           multiRoot,
                                           materialNames );
}

template< typename OUTTYPE, typename TYPE >
void SiloFile::writeMaterialDataField( string const & meshName,
                                       string const & fieldName,
                                       array1d< array1d< arrayView2d< TYPE const > > > const & field,
                                       ElementRegionBase const & elemRegion,
                                       int const centering,
                                       int const cycleNumber,
                                       real64 const problemTime,
                                       string const & multiRoot,
                                       string_array const & materialNames )
{
  int const nvars = SiloFileUtilities::GetNumberOfVariablesInField< TYPE >();
  int const meshType = getMeshType( meshName );

//  double missingValue = 0.0;
  DBoptlist *optlist = DBMakeOptlist( 5 );
  DBAddOption( optlist, DBOPT_CYCLE, const_cast< int * >(&cycleNumber));
  DBAddOption( optlist, DBOPT_DTIME, const_cast< real64 * >(&problemTime));
//  DBAddOption(optlist, DBOPT_MISSING_VALUE, &missingValue);

  char * regionpnames[ 100 ];

  localIndex numElemsInRegion = 0;
  elemRegion.forElementSubRegions< ElementSubRegionBase >( [&] ( ElementSubRegionBase const & subRegion )
  {
    numElemsInRegion += subRegion.size();
  } );

  // if the number of elements is zero, then record the path to the var. This
  // will be used later to delete the entry
  // from the multivar.
  if( numElemsInRegion == 0 )
  {
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );
    string emptyObject = pwd;
    emptyObject += "/" + fieldName;
    m_emptyVariables.emplace_back( emptyObject );
  }
  else
  {
    string_array varnamestring( nvars );
    array1d< char const * > varnames( nvars );

    SiloFileUtilities::SetVariableNames< TYPE >( fieldName, varnamestring, varnames.data() );

    int nels = 0;
    localIndex mixlen = 0;
    localIndex const numMatInRegion = materialNames.size();

    elemRegion.forElementSubRegions< ElementSubRegionBase >( [&]( ElementSubRegionBase const & subRegion )
    {
      nels += subRegion.size();
      for( localIndex matIndex=0; matIndex<numMatInRegion; ++matIndex )
      {
        mixlen += subRegion.size();
      }
    } );

    array1d< void * > vars( nvars );
    array1d< array1d< OUTTYPE > > varsData( nvars );

    array1d< void * > mixvars( nvars );
    array1d< array1d< OUTTYPE > > mixvarsData( nvars );

    for( int a=0; a<nvars; ++a )
    {
      varsData[a].resize( nels );
      mixvarsData[a].resize( mixlen );

      vars[a] = static_cast< void * >(varsData[a].data());
      mixvars[a] = static_cast< void * >(mixvarsData[a].data());
    }

    array1d< localIndex > mixlen2( nvars );
    array1d< localIndex > nels2( nvars );

    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      if( numMatInRegion == 1 )
      {
        for( int i = 0; i < nvars; ++i )
        {
          for( localIndex k = 0; k < subRegion.size(); ++k )
          {
            localIndex const numQ = field[esr][0].size( 1 );
            varsData[i][nels2[i]] = 0;
            mixvarsData[i][mixlen2[i]] = 0;
            for( localIndex q=0; q < numQ; ++q )
            {
              varsData[i][nels2[i]] += SiloFileUtilities::CastField< OUTTYPE >( field[esr][0][k][q], i );
              mixvarsData[i][mixlen2[i]] += SiloFileUtilities::CastField< OUTTYPE >( field[esr][0][k][q], i );
            }
            varsData[i][nels2[i]] /= numQ;
            mixvarsData[i][mixlen2[i]] /= numQ;

            ++nels2[i];
            ++mixlen2[i];
          }
        }
      }
      else if( numMatInRegion > 1 )
      {
        for( int i = 0; i < nvars; ++i )
        {
          for( localIndex k = 0; k < subRegion.size(); ++k )
          {
            for( localIndex a=0; a<numMatInRegion; ++a )
            {
              if( field[esr][a].size() > 0 )
              {
                localIndex const numQ = field[esr][0].size( 1 );
                mixvarsData[i][mixlen2[i]] = 0;
                for( localIndex q=0; q < numQ; ++q )
                {
                  mixvarsData[i][mixlen2[i]] = SiloFileUtilities::CastField< OUTTYPE >( field[esr][a][k][q], i );
                }
                mixvarsData[i][mixlen2[i]] /= numQ;
                ++mixlen2[i];
              }
              else
              {
                mixvarsData[i][mixlen2[i]++] = 0.0;
              }
            }

            if( field[esr][0].size() > 0 )
            {
              localIndex const numQ = field[esr][0].size( 1 );
              varsData[i][nels2[i]] = 0;
              for( localIndex q=0; q < numQ; ++q )
              {
                varsData[i][nels2[i]] += SiloFileUtilities::CastField< OUTTYPE >( field[esr][0][k][q], i );
              }
              varsData[i][nels2[i]] /= numQ;
              ++nels2[i];
            }
            else
            {
              varsData[i][nels2[i]++] = 0.0;
            }
          }
        }
      }
    } );

    for( localIndex a=0; a<materialNames.size(); ++a )
    {
      regionpnames[a] = const_cast< char * >(materialNames[a].c_str());
    }
    regionpnames[materialNames.size()] = nullptr;
    DBAddOption( optlist, DBOPT_REGION_PNAMES, &regionpnames );

    int err = -2;
//    if( meshType == DB_UCDMESH )
//    {
    err = DBPutUcdvar( m_dbFilePtr,
                       fieldName.c_str(),
                       meshName.c_str(),
                       nvars,
                       varnames.data(),
                       vars.data(),
                       nels,
                       mixvars.data(),
                       mixlen,
                       SiloFileUtilities::DB_TYPE< OUTTYPE >(),
                       centering,
                       optlist );
//    }
//    else if( meshType == DB_POINTMESH )
//    {
//      err = DBPutPointvar( m_dbFilePtr,
//                           fieldName.c_str(),
//                           meshName.c_str(),
//                           nvars,
//                           vars.data(),
//                           nels,
//                           SiloFileUtilities::DB_TYPE<OUTTYPE>(),
//                           optlist);
//    }
    if( err < 0 )
    {
      if( err < -1 )
      {
        GEOSX_ERROR( "unhandled case in SiloFile::WriteDataField A\n" );
      }
      else
      {
        GEOSX_ERROR( "unhandled failure in adding variable during SiloFile::WriteDataField\n" );
      }
    }
  }

  // write multimesh object
  int rank = 0;
#ifdef GEOSX_USE_MPI
  MPI_Comm_rank( MPI_COMM_GEOSX, &rank );
#endif
  if( rank == 0 )
  {
    int tensorRank = SiloFileUtilities::GetTensorRank< TYPE >();
    DBAddOption( optlist, DBOPT_TENSOR_RANK, const_cast< int * >(&tensorRank));
    DBAddOption( optlist, DBOPT_MMESH_NAME, const_cast< char * >(meshName.c_str()));

    DBObjectType vartype = DB_INVALID_OBJECT;

    if( meshType == DB_UCDMESH )
    {
      vartype = DB_UCDVAR;
    }
    else if( meshType == DB_POINTMESH )
    {
      vartype = DB_POINTVAR;
    }
    else if( meshType == DB_QUADCURV )
    {
      vartype = DB_QUADVAR;
    }
    else
    {
      vartype = DB_UCDVAR;
//      GEOSX_ERROR("unhandled case in SiloFile::WriteDataField B\n");
    }

    writeMultiXXXX( vartype, DBPutMultivar, centering, fieldName.c_str(), cycleNumber, multiRoot,
                    optlist );
  }

  DBFreeOptlist( optlist );
}


template< typename OUTTYPE, typename TYPE >
void SiloFile::writeMaterialDataField3d( string const & meshName,
                                         string const & fieldName,
                                         ElementRegionBase const & elemRegion,
                                         int const centering,
                                         int const cycleNumber,
                                         real64 const problemTime,
                                         string const & multiRoot,
                                         string_array const & materialNames )
{
  localIndex const numSubRegions = elemRegion.numSubRegions();
  localIndex const numMat = materialNames.size();

  array1d< array1d< array2d< TYPE > > >    fieldCopy( numSubRegions );
  array1d< array1d< arrayView2d< TYPE const > > > fieldView( numSubRegions );

  // resize the container and find the maximum size along 3rd dimension across sub-region data
  localIndex nvar = 0;
  elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
  {
    fieldCopy[esr].resize( numMat );
    fieldView[esr].resize( numMat );
    for( localIndex matIndex = 0; matIndex<numMat; ++matIndex )
    {
      arrayView3d< TYPE const > const &
      fieldData = subRegion.getConstitutiveModel( materialNames[matIndex] ).getReference< array3d< TYPE > >( fieldName );
      if( fieldData.size() > 0 )
      {
        fieldCopy[esr][matIndex].resize( fieldData.size( 0 ), fieldData.size( 1 ) );
        nvar = std::max( nvar, fieldData.size( 2 ) );
      }
      fieldView[esr][matIndex] = fieldCopy[esr][matIndex];
    }
  } );

  // loop over variables and copy into the container
  for( localIndex ivar = 0; ivar < nvar; ++ivar )
  {
    elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
      [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
    {
      for( localIndex matIndex = 0; matIndex < numMat; ++matIndex )
      {
        arrayView3d< TYPE const > const &
        fieldData = subRegion.getConstitutiveModel( materialNames[matIndex] ).getReference< array3d< TYPE > >( fieldName );

        if( fieldData.size() > 0 && ivar < fieldData.size( 2 ))
        {
          arrayView2d< TYPE > & fieldCopyData = fieldCopy[esr][matIndex];
          for( localIndex ei = 0; ei < fieldData.size( 0 ); ++ei )
          {
            for( localIndex q = 0; q < fieldData.size( 1 ); ++q )
            {
              fieldCopyData[ei][q] = fieldData[ei][q][ivar];
            }
          }
        }
      }
    } );


    string component = std::to_string( ivar );
    if( fieldName=="stress" )
    {
      if( ivar==0 )
      {
        component = "11";
      }
      else if( ivar==1 )
      {
        component = "22";
      }
      else if( ivar==2 )
      {
        component = "33";
      }
      else if( ivar==3 )
      {
        component = "23";
      }
      else if( ivar==4 )
      {
        component = "13";
      }
      else if( ivar==5 )
      {
        component = "12";
      }
    }
    string componentFieldName = fieldName + "_" + component;
    writeMaterialDataField< real64 >( meshName,
                                      componentFieldName,
                                      fieldView,
                                      elemRegion,
                                      centering,
                                      cycleNumber,
                                      problemTime,
                                      multiRoot,
                                      materialNames );
  }
}


template< typename OUTTYPE, typename TYPE >
void SiloFile::writeMaterialDataField4d( string const & meshName,
                                         string const & fieldName,
                                         ElementRegionBase const & elemRegion,
                                         int const centering,
                                         int const cycleNumber,
                                         real64 const problemTime,
                                         string const & multiRoot,
                                         string_array const & materialNames )
{
  localIndex const numSubRegions = elemRegion.numSubRegions();
  localIndex const numMat = materialNames.size();

  array1d< array1d< array2d< TYPE > > >    fieldCopy( numSubRegions );
  array1d< array1d< arrayView2d< TYPE const > > > fieldView( numSubRegions );

  // resize the container and find the maximum size along 3rd/4th dimension across sub-region data
  localIndex nvar1 = 0;
  localIndex nvar2 = 0;

  elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
    [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
  {
    fieldCopy[esr].resize( numMat );
    fieldView[esr].resize( numMat );
    for( localIndex matIndex = 0; matIndex<numMat; ++matIndex )
    {
      arrayView4d< TYPE const > const &
      fieldData = subRegion.getConstitutiveModel( materialNames[matIndex] ).getReference< array4d< TYPE > >( fieldName );
      if( fieldData.size() > 0 )
      {
        fieldCopy[esr][matIndex].resize( fieldData.size( 0 ), fieldData.size( 1 ) );
        nvar1 = std::max( nvar1, fieldData.size( 2 ) );
        nvar2 = std::max( nvar2, fieldData.size( 3 ) );
      }
      fieldView[esr][matIndex] = fieldCopy[esr][matIndex];
    }
  } );

  // loop over variables and copy into the container
  for( localIndex ivar = 0; ivar < nvar1; ++ivar )
  {
    for( localIndex jvar = 0; jvar < nvar2; ++jvar )
    {
      elemRegion.forElementSubRegionsIndex< ElementSubRegionBase >(
        [&]( localIndex const esr, ElementSubRegionBase const & subRegion )
      {
        for( localIndex matIndex = 0; matIndex < numMat; ++matIndex )
        {
          arrayView4d< TYPE const > const &
          fieldData = subRegion.getConstitutiveModel( materialNames[matIndex] ).getReference< array4d< TYPE > >( fieldName );

          if( fieldData.size() > 0 && ivar < fieldData.size( 2 ) && jvar < fieldData.size( 3 ))
          {
            arrayView2d< TYPE > & fieldCopyData = fieldCopy[esr][matIndex];
            for( localIndex ei = 0; ei < fieldData.size( 0 ); ++ei )
            {
              for( localIndex q = 0; q < fieldData.size( 1 ); ++q )
              {
                fieldCopyData[ei][q] = fieldData[ei][q][ivar][jvar];
              }
            }
          }
        }
      } );
      writeMaterialDataField< real64 >( meshName,
                                        fieldName + "_" + std::to_string( ivar ) + "_" + std::to_string( jvar ),
                                        fieldView,
                                        elemRegion,
                                        centering,
                                        cycleNumber,
                                        problemTime,
                                        multiRoot,
                                        materialNames );
    }
  }
}

void SiloFile::writeStressVarDefinition( string const & MatDir )
{
  {
    string const expressionName = "/" + MatDir + "/stress";
    const char * expObjName = expressionName.c_str();
    const char * const names[1] = { expObjName };
    int const types[1] = { DB_VARTYPE_TENSOR };
    string const definition = "{{<"+ MatDir + "/stress_11>,<"+ MatDir + "/stress_12>,<"+ MatDir + "/stress_13>},"
                                                                                                  "{<"+ MatDir + "/stress_12>,<"+ MatDir + "/stress_22>,<"+ MatDir + "/stress_23>},"
                                                                                                                                                                     "{<"+ MatDir + "/stress_13>,<"+
                              MatDir + "/stress_23>,<"+ MatDir + "/stress_33>}}";
    const char * const defns[1] = { definition.c_str() };
    DBPutDefvars( m_dbBaseFilePtr,
                  expObjName,
                  1,
                  names,
                  types,
                  defns,
                  nullptr );
  }

  {
    string const expressionName = "/" + MatDir + "/meanStress";
    const char * expObjName = expressionName.c_str();
    const char * const names[1] = { expObjName };
    int const types[1] = { DB_VARTYPE_SCALAR };
    string const definition = "(<"+ MatDir + "/stress_11>+<"+ MatDir + "/stress_22>+<"+ MatDir + "/stress_33>)/3";
    const char * const defns[1] = { definition.c_str() };
    DBPutDefvars( m_dbBaseFilePtr,
                  expObjName,
                  1,
                  names,
                  types,
                  defns,
                  nullptr );
  }

  {
    string const expressionName = "/" + MatDir + "/principalStress";
    const char * expObjName = expressionName.c_str();
    const char * const names[1] = { expObjName };
    int const types[1] = { DB_VARTYPE_VECTOR };
    string const definition = "eigenvalue(<"+MatDir+"/stress>)";
    const char * const defns[1] = { definition.c_str() };
    DBPutDefvars( m_dbBaseFilePtr,
                  expObjName,
                  1,
                  names,
                  types,
                  defns,
                  nullptr );
  }

  {
    string const expressionName = "/" + MatDir + "/principalStressDirections";
    const char * expObjName = expressionName.c_str();
    const char * const names[1] = { expObjName };
    int const types[1] = { DB_VARTYPE_TENSOR };
    string const definition = "eigenvector(<"+MatDir+"/stress>)";
    const char * const defns[1] = { definition.c_str() };
    DBPutDefvars( m_dbBaseFilePtr,
                  expObjName,
                  1,
                  names,
                  types,
                  defns,
                  nullptr );
  }

  for( int i=0; i<3; ++i )
  {
    string const expressionName = "/" + MatDir + "/principalStressVector" + std::to_string( i );
    const char * expObjName = expressionName.c_str();
    const char * const names[1] = { expObjName };
    int const types[1] = { DB_VARTYPE_VECTOR };
    string const definition = "transpose(<"+MatDir+"/principalStressDirections>)[" + std::to_string( i ) + "] * <"+ MatDir + "/principalStress>["+
                              std::to_string( i ) + "]";
    const char * const defns[1] = { definition.c_str() };
    DBPutDefvars( m_dbBaseFilePtr,
                  expObjName,
                  1,
                  names,
                  types,
                  defns,
                  nullptr );
  }
}


void SiloFile::writeVectorVarDefinition( string const & fieldName,
                                         string const & subDirectory )
{
  if( MpiWrapper::commRank( MPI_COMM_GEOSX ) == 0 )
  {
    DBSetDir( m_dbBaseFilePtr, subDirectory.c_str() );
    DBtoc * const siloTOC = DBGetToc ( m_dbBaseFilePtr );
    int numComponents = 0;
    for( int ivar=0; ivar<siloTOC->nmultivar; ++ivar )
    {
      string const varName = siloTOC->multivar_names[ivar];
      if( ( varName == fieldName + "_0" ) ||
          ( varName == fieldName + "_1" ) ||
          ( varName == fieldName + "_2" ) )
      {
        ++numComponents;
      }
    }

    if( numComponents==3 )
    {
      string const expressionName = subDirectory + "/" + fieldName;
      const char * expObjName = expressionName.c_str();
      const char * const names[1] = { expObjName };
      int const types[1] = { DB_VARTYPE_VECTOR };
      string const definition = "{<"+ subDirectory + "/" + fieldName + "_0>,<"+ subDirectory + "/" + fieldName + "_1>,<"+ subDirectory + "/" + fieldName +
                                "_2>}";
      const char * const defns[1] = { definition.c_str() };
      DBPutDefvars( m_dbBaseFilePtr,
                    expObjName,
                    1,
                    names,
                    types,
                    defns,
                    nullptr );
    }

    DBSetDir( m_dbBaseFilePtr, ".." );
  }


}


int SiloFile::getMeshType( string const & meshName ) const
{
  int meshType = -1;
  {
    // in order to get mesh type, we might have to go up a few levels in the
    // silo directory structure
    // before we can find the mesh.
    char pwd[256];
    DBGetDir( m_dbFilePtr, pwd );

    for( int i=0; i<3; ++i )
    {
      meshType = DBInqMeshtype( m_dbFilePtr, meshName.c_str());
      if( meshType != -1 && meshType != 610 )
      {
        break;
      }
      else
      {
        DBSetDir( m_dbFilePtr, ".." );
      }
    }
    DBSetDir( m_dbFilePtr, pwd );
  }
  return meshType;
}

}
#pragma GCC diagnostic pop
