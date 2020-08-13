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
 * @file SiloFile.hpp
 */

#ifndef GEOSX_FILEIO_SILO_SILOFILE_HPP_
#define GEOSX_FILEIO_SILO_SILOFILE_HPP_

#include "common/DataTypes.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "mpiCommunications/MpiWrapper.hpp"
#include "mesh/ElementRegionManager.hpp"
#include "mesh/CellElementSubRegion.hpp"
#include "mesh/FaceElementSubRegion.hpp"
#include "mesh/InterObjectRelation.hpp"

#include "silo.h"

/// _PMPIO_baton_t struct forward declaration
struct _PMPIO_baton_t;

/// Type alias for _PMPIO_baton_t struct
typedef _PMPIO_baton_t PMPIO_baton_t;

namespace geosx
{
class DomainPartition;
class MeshLevel;

namespace constitutive
{
class ConstitutiveManager;
}

// *********************************************************************************************************************
// *********************************************************************************************************************
/**
 * This class serves as a wrapper to isolate the code from the specifics of SILO
 * output/input. Its members contain all the necessary information for reading/writing a group of SILO files.
 */
class SiloFile
{
public:
  /// Default Constructor
  SiloFile();

  /// Destructor
  virtual ~SiloFile();

  /**
   * @brief function to setup directories where silo files will be written
   */
  void
  MakeSiloDirectories();

  /**
   * @brief Initializes silo for input/output
   * @param numGroups number of individual Silo files to generate
   */
  void
  Initialize( int const numGroups = 1 );

  /**
   * @brief finishes/closes up the silo interface
   */
  void
  Finish();

  /**
   *  @brief obtain the group number of the calling processor, indexed from zero.
   *  @param i rank of calling processor in the MPI communicator
   *  @return group number of the calling processor, indexed from zero
   */
  int
  groupRank( int const i ) const;

  /**
   * @brief Wait for the Baton when writing using PMPIO
   * @param domainNumber domain partition number
   * @param cycleNum  cycle number of simulation
   * @param eventCounter Counter to indicate the event number during the current timestep.
   * @param isRestart whether or not we are writing a restart file
   *
   * This function requests the write baton from silo PMPIO. The involves determining
   * the file names, and opening the file for write.
   */
  void
  WaitForBatonWrite( int const domainNumber,
                     int const cycleNum,
                     integer const eventCounter,
                     bool const isRestart );

  /**
   * @brief Wait for the Baton when reading using PMPIO
   * @param domainNumber domain partition number
   * @param restartFileName base of the restart file to open
   */
  void
  WaitForBaton( int const domainNumber, string const & restartFileName );

  /**
   * @brief Hand off the Baton when done writing to file
   */
  void
  HandOffBaton();

  /**
   * @brief Make a subdirectory within the silo file.
   * @param subdir the new directory name
   * @param rootdir the root directory path
   */
  void
  MakeSubDirectory( string const & subdir, string const & rootdir )
  {
    int const rank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );

    // char dirname[100];
    if( rank == 0 )
    {
      //      DBGetDir(m_dbBaseFilePtr, dirname );
      DBMkDir( m_dbBaseFilePtr, rootdir.c_str() );
    }

    //    DBGetDir (m_dbFilePtr, dirname );
    DBMkDir( m_dbFilePtr, subdir.c_str() );
  }

  /**
   * @brief Write out a single silo mesh object
   * @param meshName name of the mesh in the silo db
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param globalNodeNum array to the global node numbers. This might be redundant as there is a field for this.
   * @param ghostNodeName Character array for whether or not a node is a ghost. ( 0 = local, 1 = ghost )
   * @param ghostZoneName Character array for whether or not a zone is a ghost. ( 0 = local, 1 = ghost )
   * @param numShapes number of element zone type (i.e. number of zone types with different topology)
   * @param shapecnt pointer to array that contains the number of zones per shape type
   * @param meshConnectivity pointer to array that contains the zone to element map for each  zone type
   * @param globalElementNum pointer to array of global zone numbers for each shape type
   * @param shapetype pointer to array containing the shape types
   * @param shapesize pointer to array containing the number of nodes in each zone in the shape types
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   *
   * This function takes in the the required data to call a silo::DBPutUcdMesh() and a
   * DBPutZonelist2, and calls those functions to create a silo mesh object. In addition
   * the MultiVar is written in the root file.
   */
  void
  WriteMeshObject( string const & meshName,
                   const localIndex nnodes,
                   real64 * coords[3],
                   const globalIndex * globalNodeNum,
                   char const * const ghostNodeName,
                   char const * const ghostZoneName,
                   int const numShapes,
                   int const * shapecnt,
                   const localIndex * const * const meshConnectivity,
                   const globalIndex * const * const globalElementNum,
                   int const * const shapetype,
                   int const * const shapesize,
                   int const cycleNumber,
                   real64 const problemTime );

  /**
   * @todo Verify: documentation missing / incomplete. The TPL version of doxygen on Travis cannot parse
   * unnamed parameters, @p dummy parameter introduced to remove warning
   *
   * @param meshName name of the mesh in the silo db
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param dummy1 unused parameter
   * @param numRegions
   * @param shapecnt pointer to array that contains the number of zones per shape type
   * @param meshConnectivity pointer to array that contains the zone to element map for each  zone type
   * @param globalElementNum pointer to array of global zone numbers for each shape type
   * @param dummy2 unused parameter
   * @param shapetype pointer to array containing the shape types
   * @param shapesize pointer to array containing the number of nodes in each zone in the shape types
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param lnodelist
   */
  void
  WritePolygonMeshObject( const std::string & meshName,
                          const localIndex nnodes,
                          realT * coords[3],
                          const globalIndex * dummy1,
                          const int numRegions,
                          const int * shapecnt,
                          const localIndex * const * const meshConnectivity,
                          const globalIndex * const * const globalElementNum,
                          const int * const * const dummy2,
                          const int * const shapetype,
                          const int * const shapesize,
                          const int cycleNumber,
                          const realT problemTime,
                          const int lnodelist );

  /**
 * @brief write a domain parititon out to silo file
 * @param domain the domain partition to write
 * @param cycleNum the current cycle number
 * @param problemTime the current problem time
 * @param isRestart whether or not we want to write restart only data
 */
  void
  WriteDomainPartition( DomainPartition const & domain,
                        int const cycleNum,
                        real64 const problemTime,
                        bool const isRestart );

  /**
   * @todo Verify: documentation missing / incomplete.
   *
   * @param elementRegion the element region that holds the data to be written to the silo file
   * @param nodeManager the NodeManager containing the nodes of the domain to be output
   * @param meshName name of the mesh to write
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param globalNodeNum array to the global node numbers. This might be redundant as there is a field for this.
   * @param ghostNodeFlag
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param writeArbitraryPolygon
   */
  void
  WriteElementMesh( ElementRegionBase const & elementRegion,
                    NodeManager const * const nodeManager,
                    string const & meshName,
                    const localIndex nnodes,
                    real64 * coords[3],
                    globalIndex const * const globalNodeNum,
                    char const * const ghostNodeFlag,
                    int const cycleNumber,
                    real64 const problemTime,
                    bool & writeArbitraryPolygon );

  /**
   * @brief write a mesh level out to the silo file
   * @param meshLevel the meshLevel to write out
   * @param cycleNum the current cycle number
   * @param problemTime the current problem time
   * @param isRestart whether or not we want to write restart only data
   */
  void
  WriteMeshLevel( MeshLevel const * const meshLevel,
                  int const cycleNum,
                  real64 const problemTime,
                  bool const isRestart );

  /**
   * @param meshName name of the mesh to write
   * @param numPoints number of points to write
   * @param coords array[3] of pointers to x, y, and z.
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void
  WritePointMesh( string const & meshName,
                  const localIndex numPoints,
                  real64 * coords[3],
                  int const cycleNumber,
                  real64 const problemTime );

  /**
   * @brief write a beam mesh to silo file
   * @param meshName name of the mesh to write
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param node1 array of the first node for each element
   * @param node2 array of the second node for each element
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void
  WriteBeamMesh( string const & meshName,
                 const localIndex nnodes,
                 real64 * coords[3],
                 const localIndex_array & node1,
                 const localIndex_array & node2,
                 int const cycleNumber,
                 real64 const problemTime );

  /**
   *
   * @param meshName name of the mesh to write
   * @param nnodes number of nodes
   * @param coords array[3] of pointers to x, y, and z.
   * @param nodelist nodal connectivity array
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void
  WriteBeamMesh( string const & meshName,
                 const localIndex nnodes,
                 real64 * coords[3],
                 integer_array & nodelist,
                 int const cycleNumber,
                 real64 const problemTime );

  /**
   * @param elementRegion the element region that holds the data to be written to the silo file
   * @param meshName name of the mesh to write
   * @param regionMaterialList  region material list
   * @param cycleNumber current cycle number
   * @param problemTime current problem time
   */
  void
  WriteMaterialMapsFullStorage( ElementRegionBase const & elementRegion,
                                string const & meshName,
                                string_array const & regionMaterialList,
                                int const cycleNumber,
                                real64 const problemTime );
  /**
   *
   * @param group the group that holds the data to be written to the silo file
   * @param siloDirName the name of the silo directory to put this data into (i.e. nodalFields, elemFields)
   * @param meshname the name of the mesh attach this write to
   * @param centering the centering of the data (e.g. DB_ZONECENT)
   * @param cycleNum current cycle number
   * @param problemTime current problem time
   * @param isRestart write restart only data
   * @param mask indices to write out to the silo file
   */
  void
  WriteGroupSilo( dataRepository::Group const * group,
                  string const & siloDirName,
                  string const & meshname,
                  int const centering,
                  int const cycleNum,
                  real64 const problemTime,
                  bool const isRestart,
                  const localIndex_array & mask );

  /**
   * @param elemRegion the element region that holds the data to be written to the silo file
   * @param siloDirName the name of the silo directory to put this data into (i.e. nodalFields, elemFields)
   * @param meshName the name of the mesh attach this write to
   * @param cycleNum current cycle number
   * @param problemTime current problem time
   * @param isRestart write restart only data
   */
  void
  WriteElementRegionSilo( ElementRegionBase const & elemRegion,
                          string const & siloDirName,
                          string const & meshName,
                          int const cycleNum,
                          real64 const problemTime,
                          bool const isRestart );
  /**
   * Writes the contents of a group of Wrapper objects
   * @tparam the output variable type
   * @param meshname the name of the mesh attach this write to
   * @param wrappers a group of wrappers
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNum the current cycle number
   * @param problemTime the current problem time
   * @param isRestart write restart only data
   * @param multiRoot location to write the multivar entries
   * @param mask indices to write out to the silo file
   */
  template< typename OUTPUTTYPE >
  void
  WriteWrappersToSilo( string const & meshname,
                       const dataRepository::Group::wrapperMap & wrappers,
                       int const centering,
                       int const cycleNum,
                       real64 const problemTime,
                       bool const isRestart,
                       string const & multiRoot,
                       const localIndex_array & mask );

  /**
   *
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   */
  template< typename OUTTYPE, typename TYPE >
  void
  WriteDataField( string const & meshName,
                  string const & fieldName,
                  arrayView1d< TYPE const > const & field,
                  int const centering,
                  int const cycleNumber,
                  real64 const problemTime,
                  string const & multiRoot );

  /**
   *
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   */
  template< typename OUTTYPE, typename TYPE, int USD >
  void
  WriteDataField( string const & meshName,
                  string const & fieldName,
                  arrayView2d< TYPE const, USD > const & field,
                  int const centering,
                  int const cycleNumber,
                  real64 const problemTime,
                  string const & multiRoot );

  /**
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   */
  template< typename OUTTYPE, typename TYPE, int USD >
  void
  WriteDataField( string const & meshName,
                  string const & fieldName,
                  arrayView3d< TYPE const, USD > const & field,
                  int const centering,
                  int const cycleNumber,
                  real64 const problemTime,
                  string const & multiRoot );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param siloTensorRank <B>****** UNUSED IN THE IMPLEMENTATION ****** </B>
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   */
  template< typename OUTTYPE, typename TYPE, int NDIM, int USD >
  void
  WriteDataField( string const & meshName,
                  string const & fieldName,
                  ArrayView< TYPE const, NDIM, USD > const & field,
                  int const siloTensorRank,
                  int const centering,
                  int const cycleNumber,
                  real64 const problemTime,
                  string const & multiRoot );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param field field data
   * @param elemRegion the element region that holds the data to be written to the silo file
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   * @param materialNames material names
   */
  template< typename OUTTYPE, typename TYPE >
  void
  WriteMaterialDataField( string const & meshName,
                          string const & fieldName,
                          array1d< array1d< arrayView2d< TYPE const > > > const & field,
                          ElementRegionBase const & elemRegion,
                          int const centering,
                          int const cycleNumber,
                          real64 const problemTime,
                          string const & multiRoot,
                          string_array const & materialNames );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param elemRegion the element region that holds the data to be written to the silo file
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   * @param materialNames material names
   */
  template< typename OUTTYPE, typename TYPE >
  void
  WriteMaterialDataField2d( string const & meshName,
                            string const & fieldName,
                            ElementRegionBase const & elemRegion,
                            int const centering,
                            int const cycleNumber,
                            real64 const problemTime,
                            string const & multiRoot,
                            string_array const & materialNames );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param elemRegion the element region that holds the data to be written to the silo file
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   * @param materialNames material names
   */
  template< typename OUTTYPE, typename TYPE >
  void
  WriteMaterialDataField3d( string const & meshName,
                            string const & fieldName,
                            ElementRegionBase const & elemRegion,
                            int const centering,
                            int const cycleNumber,
                            real64 const problemTime,
                            string const & multiRoot,
                            string_array const & materialNames );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param meshName the name of the mesh attach this write to
   * @param fieldName name of the field to write
   * @param elemRegion the element region that holds the data to be written to the silo file
   * @param centering the silo centering to use for this operation (DB_NODECENT, DB_ZONECENT)
   * @param cycleNumber the current cycle number
   * @param problemTime the current problem time
   * @param multiRoot location to write the multivar entries
   * @param materialNames material names
   */
  template< typename OUTTYPE, typename TYPE >
  void
  WriteMaterialDataField4d( string const & meshName,
                            string const & fieldName,
                            ElementRegionBase const & elemRegion,
                            int const centering,
                            int const cycleNumber,
                            real64 const problemTime,
                            string const & multiRoot,
                            string_array const & materialNames );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param subDir
   * @param matDir
   * @param matIndex
   * @param fieldName
   */
  void
  WriteMaterialVarDefinition( string const & subDir,
                              string const & matDir,
                              localIndex const matIndex,
                              string const & fieldName );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param MatDir
   */
  void
  WriteStressVarDefinition( string const & MatDir );

  /**
   * @todo Verify: documentation missing / incomplete
   * @param fieldName vector field name
   * @param subDirectory
   */
  void
  WriteVectorVarDefinition( string const & fieldName,
                            string const & subDirectory );

  /**
   * find the silo mesh type that we are attempting to reference
   * @param meshName the name of the mesh object we are attaching data to
   * @return integer that results from a call to DBInqMeshtype()
   */
  int
  GetMeshType( string const & meshName ) const;

  /**
   * write out a multi mesh or multivar object assocaited with a mesh or var object.
   * @param type the object type (DB_UCDVAR, ...)
   * @param DBPutMultiCB a callback function for the placement of the multi object. ( DBPutMultivar, DBPutMultimesh...)
   * @param centering the centering of the data (DB_NODECENT, DB_ZONECENT, ...)
   * @param name the name of the multivar
   * @param cycleNumber the current cycle number
   * @param multiRoot the location in the silo file to put hte multiXXXX
   * @param optlist the option list assocaited with the multiXXXX
   */
  template< typename CBF >
  void
  WriteMultiXXXX( const DBObjectType type,
                  CBF DBPutMultiCB,
                  int const centering,
                  string const name,
                  int const cycleNumber,
                  string const & multiRoot,
                  const DBoptlist * optlist = nullptr );

  /**
   * function to clear any empty multi-objects
   * @param cycleNum
   *
   * When we write our multimesh and multivar objects, we do it assuming that there is a non-empty multimesh
   * or multivar object on each domain. This is incorrect, so we must modify the rootfile to remove empty references
   * to the multivar or multimesh objects and replace their path with "EMPTY"
   */
  void
  ClearEmptiesFromMultiObjects( int const cycleNum );

  /**
   * @brief Sets the number of individual Silo files to generate
   * @param numGroups number of individual Silo files to generate
   */
  void
  setNumGroups( int const numGroups )
  {
    m_numGroups = numGroups;
  }

  /**
   * @brief Sets the plot level option
   * @param plotLevel the plot level desired value
   */
  void
  setPlotLevel( int const plotLevel )
  {
    m_plotLevel = dataRepository::toPlotLevel( plotLevel );
  }

  /**
   * @brief Sets the edge mesh output option
   * @param val if 1, the edge mesh is written to a Silo file
   */
  void
  setWriteEdgeMesh( int const val )
  {
    m_writeEdgeMesh = val;
  }

  /**
   * @brief Sets the face mesh output option
   * @param val if 1, the face mesh is written to a Silo file
   */
  void
  setWriteFaceMesh( int const val )
  {
    m_writeFaceMesh = val;
  }

  /**
   * @brief Sets the cell element mesh output option
   * @param val if 1 the cell element mesh is written to a Silo file
   */
  void
  setWriteCellElementMesh( int const val )
  {
    m_writeCellElementMesh = val;
  }

  /**
   * @brief Sets the face element mesh output option
   * @param val if 1 the face element mesh is written to a Silo file
   */
  void
  setWriteFaceElementMesh( int const val )
  {
    m_writeFaceElementMesh = val;
  }

  /**
   * @brief Sets root of the filename that will be read/written
   * @param fileRoot root of the filename
   */
  void
  setPlotFileRoot( string const & fileRoot )
  {
    m_plotFileRoot = fileRoot;
  }

private:
  /// pointer to the DBfile that this class is working on
  DBfile * m_dbFilePtr;

  /// pointer to the Master DBfile
  DBfile * m_dbBaseFilePtr;

  /// total number of "parallel" files to write out
  int m_numGroups;

  /// the pmpio baton. A processor needs this to write to the file.
  PMPIO_baton_t * m_baton;

  /// which database to use. DB_PDB or DB_HDF5
  int const m_driver;

  /// root of the filename that we will be reading/writing
  string m_plotFileRoot;

  string m_restartFileRoot;

  string const m_siloDirectory = "siloFiles";

  string const m_siloDataSubDirectory = "data";

  string m_fileName;

  string m_baseFileName;

  std::vector< std::string > m_emptyMeshes;
  //  string_array m_emptyMaterials;
  std::vector< std::string > m_emptyVariables;

  integer m_writeEdgeMesh;
  integer m_writeFaceMesh;
  integer m_writeCellElementMesh;
  integer m_writeFaceElementMesh;

  dataRepository::PlotLevel m_plotLevel;

  bool m_ghostFlags;
};

/**
 * Namespace to hold some utilities needed by the functions in SiloFile.
 */
namespace SiloFileUtilities
{
/**
 * @tparam OUTTYPE the type of data to write out (int,float,real64)
 * @return the integer identifier associated with the enum DBdatatype defined in
 * the silo.h file.
 *
 * This templated function is a "specialization only" definition. There is no
 * general definition, only specializations for predetermined data types.
 */
template< typename OUTTYPE >
int
DB_TYPE();

/**
 * @tparam TYPE the data type in question
 * @return the number of "variables" in a data field. For instance, 1 for a int,
 * float or real64.
 * 1 for a scalar
 * 3 for a R1Tensor...etc.
 */
template< typename TYPE >
int
GetNumberOfVariablesInField();

/**
 * @tparam TYPE the data type in question
 * @return the silo DB_VARTYPE specifier that describes the rank of the variable.
 */
template< typename TYPE >
int
GetTensorRank();

/**
 * @tparam OUTTYPE the type to cast to
 * @tparam TYPE the type to cast from
 * @param field the value to cast
 * @param i the component of the varaible to cast, assuming there is dimensionaliy to the variable
 * @return the casted value
 */
template< typename OUTTYPE, typename TYPE >
OUTTYPE
CastField( const TYPE & field, int const i = 0 );  // avoids compiler
                                                   // warning

template< typename OUTTYPE, typename TYPE >
OUTTYPE
CastField( const TYPE & field, int const i )
{
  return field.Data()[i];
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline int
CastField< int, int >( const int & field, int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return field;
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline long int
CastField< long int, long int >( const long int & field,
                                 int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return field;
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline int
CastField< int, long int >( const long int & field, int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return LvArray::integerConversion< int >( field );
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline long long int
CastField< long long int, long long int >(
  const long long int & field,
  int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return field;
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline int
CastField< int, long long int >( const long long int & field,
                                 int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return LvArray::integerConversion< int >( field );
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline real64
CastField< real64, real64 >( const real64 & field, int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return field;
}

/**
 * @todo Verify: the TPL version of doxygen on Travis cannot parse unnamed parameters, @p dummy
 *       parameter introduced to remove warning
 * @param field the value to cast
 * @param dummy unused parameter
 * @return the casted value
 */
template<>
inline float
CastField< float, real64 >( const real64 & field, int const dummy )
{
  GEOSX_UNUSED_VAR( dummy );
  return static_cast< float >( field );
}

/**
 * @tparam the type of the field
 * @param fieldName the name of the field
 * @param varnamestring an array of strings that hold the generated names
 * @param varnames a pointer to each of the strings that hold the generated names
 *
 * This function sets variable names for a field. In some cases, this will just fill varnamestring with
 * the input variable, and in other cases such as in the case of a tensor, the names will have a suffix
 * with the component of the tensor.
 */
template< typename TYPE >
void
SetVariableNames( string const & fieldName,
                  string_array & varnamestring,
                  char const * varnames[] );

}  // namespace SiloFileUtilities

}  // namespace geosx
#endif /* GEOSX_FILEIO_SILO_SILOFILE_HPP_ */
