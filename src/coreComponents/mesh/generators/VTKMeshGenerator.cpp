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
 * @file VTKMeshGenerator.cpp
 */

#include "VTKMeshGenerator.hpp"

#include "Elements/Element.hpp"
#include "MeshDataWriters/Variable.hpp"
#include "mesh/DomainPartition.hpp"

#include "mesh/mpiCommunications/PartitionBase.hpp"
#include "mesh/mpiCommunications/SpatialPartition.hpp"
#include "common/MpiWrapper.hpp"
#include "Mesh/MeshFactory.hpp"

#include "MeshDataWriters/MeshParts.hpp"

#include "mesh/MeshBody.hpp"

#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLPUnstructuredGridReader.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkRedistributeDataSetFilter.h>
#include <vtkGenerateGlobalIds.h>

#include <math.h>
namespace geosx
{
using namespace dataRepository;



VTKMeshGenerator::VTKMeshGenerator( string const & name, Group * const parent ):
  MeshGeneratorBase( name, parent )
{

  registerWrapper( viewKeyStruct::filePathString(), &m_filePath ).
    setInputFlag( InputFlags::REQUIRED ).
    setRestartFlags( RestartFlags::NO_WRITE ).
    setDescription( "path to the mesh file" );
}

void VTKMeshGenerator::postProcessInput()
{
}

Group * VTKMeshGenerator::createChild( string const & GEOSX_UNUSED_PARAM( childKey ), string const & GEOSX_UNUSED_PARAM( childName ) )
{
  return nullptr;
}

/**
 * @brief Check if two boxes are intersected
 * @param[in] box1 the first box [xMin xMax yMin yMax zMin zMax]
 * @param[in] box2 the second box [xMin xMax yMin yMax zMin zMax]
 * @return true if the boxes are intersecting, false otherwise
 */
bool checkIntersection(double box1[6], double box2[6])
{
    //GEOSX_LOG_RANK(box2[0] << " and " box1[0] and)
    return box1[0] < box2[1] && box1[1]> box2[0] && box1[2] < box2[3] && box1[3]> box2[2] && box1[4] < box2[5] && box1[5]> box2[4];
}

void VTKMeshGenerator::generateMesh( DomainPartition & domain )
{
  vtkSmartPointer<vtkMultiProcessController> controller = getVTKController();
  vtkMultiProcessController::SetGlobalController(controller);

  vtkSmartPointer<vtkUnstructuredGrid> loadedMesh = loadVTKMesh();
  MpiWrapper::barrier();

  redistributeMesh(loadedMesh);
  MpiWrapper::barrier();

  GEOSX_LOG_RANK_0(" Successfully redistributed the mesh on the " << controller->GetNumberOfProcesses() << " MPI processes");

  Group & meshBodies = domain.getMeshBodies();
  MeshBody & meshBody = meshBodies.registerGroup< MeshBody >( this->getName() );

  MeshLevel & meshLevel0 = meshBody.registerGroup< MeshLevel >( string( "Level0" ));
  NodeManager & nodeManager = meshLevel0.getNodeManager();

  double globalLength = writeMeshNodes( nodeManager );
  meshBody.setGlobalLengthScale( globalLength);

  computePotentialNeighborLists( domain, globalLength );

  localIndex numHex =0;
  localIndex numTet = 0;
  localIndex numWedge = 0;
  localIndex numPyr = 0;

  // region id link to their number of cells
  std::map<int,localIndex> regionsHex;
  std::map<int,localIndex> regionsTetra;
  std::map<int,localIndex> regionsWedges;
  std::map<int,localIndex> regionsPyramids;

  // surface and greionids
  std::set<int> regions;
  std::set<int> surfaces; //local to this rank
  std::vector<int> allSurfaces; // global to the whole mesh

  countCellsAndFaces( numHex, numTet, numWedge, numPyr, regionsHex, regionsTetra, regionsWedges, regionsPyramids, regions, surfaces, allSurfaces);

  std::vector< vtkDataArray * > arraysToBeImported = findArrayToBeImported();

  writeCellBlocks( domain, numHex, numTet, numWedge, numPyr,regionsHex, regionsTetra, regionsWedges, regionsPyramids, regions, arraysToBeImported );

  writeSurfaces( nodeManager, allSurfaces);

  GEOSX_LOG_RANK_0("Mesh was loaded successfully");
}

vtkSmartPointer<vtkMultiProcessController> VTKMeshGenerator::getVTKController()
{
  #ifdef GEOSX_USE_MPI
    vtkNew<vtkMPIController> controller;
    vtkMPICommunicatorOpaqueComm vtkGeosxComm(&MPI_COMM_GEOSX);
    vtkNew<vtkMPICommunicator> communicator;
    communicator->InitializeExternal( &vtkGeosxComm);
    controller->SetCommunicator(communicator);
  #else
    vtkNew<vtkDummyController> controller;
  #endif
  return controller;
}

vtkSmartPointer<vtkUnstructuredGrid> VTKMeshGenerator::loadVTKMesh()
{
  string_array filePathTokenized = stringutilities::tokenize( m_filePath, "."); //TODO maybe code a method in Path to get the file extension?
  string extension = filePathTokenized[filePathTokenized.size() - 1];
  vtkSmartPointer< vtkUnstructuredGrid > loadedMesh = vtkSmartPointer< vtkUnstructuredGrid >::New();

  auto read = [&](auto vtkUgReader) { // auto can be multiple types in the same function
    vtkUgReader->SetFileName( m_filePath.c_str() );
    vtkUgReader->Update();
    loadedMesh = vtkUgReader->GetOutput();
  };

  if( extension == "pvtu" )
  {
    read(vtkSmartPointer<vtkXMLPUnstructuredGridReader>::New());
  }
  else
  {
    if( MpiWrapper::commRank() == 0 )
    {
    if( extension == "vtk")
      {
        read(vtkSmartPointer<vtkUnstructuredGridReader>::New());
      }
      else if( extension == "vtu")
      {
        read(vtkSmartPointer<vtkXMLUnstructuredGridReader>::New());
      }
      else
      {
        GEOSX_ERROR( extension << " is not a recognized extension for using the VTK reader with GEOSX. Please use .vtk or .vtu" );
      }
    } 
  }
  GEOSX_LOG_RANK( loadedMesh->GetNumberOfCells());

  return loadedMesh;
}

void VTKMeshGenerator::redistributeMesh( vtkUnstructuredGrid * loadedMesh)
{
  // Redistribute data all over the available ranks
  vtkNew<vtkRedistributeDataSetFilter> rdsf;
  rdsf->SetInputDataObject( loadedMesh );
  rdsf->SetNumberOfPartitions( MpiWrapper::commSize() );
  rdsf->GenerateGlobalCellIdsOn();
  rdsf->Update();
  MpiWrapper::barrier();

  // Generate global IDs for vertices and cells
  vtkNew<vtkGenerateGlobalIds> generator;
  generator->SetInputDataObject(vtkUnstructuredGrid::SafeDownCast( rdsf->GetOutputDataObject(0)));
  generator->Update();
  m_vtkMesh = vtkUnstructuredGrid::SafeDownCast( generator->GetOutputDataObject(0));
  GEOSX_LOG_RANK( m_vtkMesh->GetNumberOfCells());
}

double VTKMeshGenerator::writeMeshNodes(NodeManager & nodeManager) const
{
  nodeManager.resize( m_vtkMesh->GetNumberOfPoints() );
  arrayView1d< globalIndex > const & nodeLocalToGlobal = nodeManager.localToGlobalMap();

  // Writing the points
  GEOSX_ERROR_IF( m_vtkMesh->GetNumberOfPoints() == 0, "Mesh is empty, aborting");
  GEOSX_LOG_RANK_0("Writing " << m_vtkMesh->GetNumberOfPoints() << " vertices");

  arrayView2d< real64, nodes::REFERENCE_POSITION_USD > const & X = nodeManager.referencePosition();
  Group & nodeSets = nodeManager.sets();
  SortedArray< localIndex > & allNodes  = nodeSets.registerWrapper< SortedArray< localIndex > >( string( "all" ) ).reference();
  real64 xMax[3] = { std::numeric_limits< real64 >::min() };
  real64 xMin[3] = { std::numeric_limits< real64 >::max() };

  vtkPointData * pointData = m_vtkMesh->GetPointData();
  int globalPointIdIndex = -1;
  vtkAbstractArray * globalPointIdArray = pointData->GetAbstractArray("GlobalPointIds", globalPointIdIndex);
  GEOSX_ERROR_IF(globalPointIdIndex == -1, "Array \"vtkGlobalPointIds\" not found" );
  GEOSX_ERROR_IF(globalPointIdArray->GetDataType() != VTK_ID_TYPE, "Array name \"vtkGlobalPointIds\" do not have the right underlying type");
  vtkIdTypeArray * globalPointIdDataArray =  nullptr;
  if( globalPointIdIndex >= 0)
  {
    globalPointIdDataArray = vtkIdTypeArray::SafeDownCast( globalPointIdArray );
  }
  for( vtkIdType v = 0; v < m_vtkMesh->GetNumberOfPoints(); v++)
  {
    double * point = m_vtkMesh->GetPoint( v );
    nodeLocalToGlobal[v] = globalPointIdDataArray->GetValue( v );
    for( integer i = 0; i < 3; i++)
    {
      X(v,i) = point[i];
      if( X( v, i ) > xMax[i] )
      {
        xMax[i] = X( v, i );
      }
      if( X( v, i ) < xMin[i] )
      {
        xMin[i] = X( v, i );
      }
    }
    allNodes.insert( v );
  }
  LvArray::tensorOps::subtract< 3 >( xMax, xMin );
  return LvArray::tensorOps::l2Norm< 3 >( xMax );
}

void VTKMeshGenerator::computePotentialNeighborLists( DomainPartition & domain, double globalLength)
{
  double * rankBounds = m_vtkMesh->GetBounds();
  double haloSize = globalLength/1000.;
  /// compute the halo 
  for(unsigned int i = 0; i < 6;i++)
  {
    rankBounds[i] += i%2?-haloSize:haloSize;
  }

  std::vector<double> allBounds(6 * MpiWrapper::commSize());
  MpiWrapper::allgather<double, double>(rankBounds, 6, allBounds.data(),6,MPI_COMM_GEOSX);

  //Checking the intersection
  std::set< int > rankNeighbor;
  for(int i = 0; i < MpiWrapper::commSize(); i++)
  {
    if( i == MpiWrapper::commRank())
    {
      continue;
    }
    std::vector< double >::const_iterator first = allBounds.begin() + 6*i;
    std::vector< double >::const_iterator last = allBounds.begin() + 6*(i+1);
    std::vector< double > curBound(first,last);
    if(checkIntersection(rankBounds, curBound.data()))
    {
      domain.getMetisNeighborList().insert(i);
    }
  }
}

vtkIntArray * VTKMeshGenerator::getAttributeDataArray()
{
  vtkCellData * cellData = m_vtkMesh->GetCellData();
  int attributeIndex = -1;
  vtkAbstractArray * attributeArray = cellData->GetAbstractArray("attribute", attributeIndex);
  GEOSX_ERROR_IF(attributeIndex >=0 && attributeArray->GetDataType() != VTK_INT, "Array name \"attribute\" is reserved to store the regions " 
                  << "used by GEOSX. The underlying type has to be an INT");
  vtkIntArray * attributeDataArray =  nullptr;
  if( attributeIndex >= 0)
  {
    attributeDataArray = vtkIntArray::SafeDownCast( attributeArray );
  }
  return attributeDataArray;
}
  vtkIdTypeArray * VTKMeshGenerator::getCellGlobalIdDataArray()
  {
    int globalCellIdsIndex = -1;
    vtkCellData * cellData = m_vtkMesh->GetCellData();
    vtkAbstractArray * globalCellIdsArray = cellData->GetAbstractArray("GlobalCellIds", globalCellIdsIndex);
    GEOSX_ERROR_IF(globalCellIdsIndex >=0 && globalCellIdsArray->GetDataType() != VTK_ID_TYPE, "Array name \"GlobalCellIds\" is not present or do not have " 
                  << " the underlying type id type");
    vtkIdTypeArray * globalCellIdsDataArray =  nullptr;
    if( globalCellIdsIndex >= 0)
    {
      globalCellIdsDataArray = vtkIdTypeArray::SafeDownCast( globalCellIdsArray );
    }
    return globalCellIdsDataArray;
  }
void VTKMeshGenerator::countCellsAndFaces( localIndex & numHex, localIndex & numTet, localIndex & numWedge, localIndex & numPyr,
                                           std::map<int,localIndex> & regionsHex, std::map<int,localIndex> & regionsTetra,
                                           std::map<int,localIndex> & regionsWedges, std::map<int,localIndex> & regionsPyramids,
                                           std::set< int > & regions, std::set< int > & surfaces, std::vector< int > & allSurfaces)
{
  vtkIntArray * attributeDataArray =  getAttributeDataArray();

  auto countCell = [&](localIndex c, localIndex & numElem, std::map<int,localIndex> & regionsElem ) { // auto can be multiple types in the same function
    numElem++;
    if (attributeDataArray != nullptr)
    {
      int region = attributeDataArray->GetValue(c);
      GEOSX_LOG_RANK_0_IF( region < 0,"Attribute value " << region << " not supported, please use value >=0 to describe regions");
      regions.insert(region);
      std::map< int, localIndex>::iterator it = regionsElem.find(region);
      if( it == regionsElem.end())
      {
        regionsElem[region] = 1;
      }
      else
      {
        it->second++;
      }
    }
  };


  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    if( m_vtkMesh->GetCellType(c) == VTK_HEXAHEDRON )
    {
      countCell(c, numHex, regionsHex);
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_TETRA )
    {
      countCell(c, numTet, regionsTetra);
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_WEDGE )
    {
      countCell(c, numWedge, regionsWedges);
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_PYRAMID )
    {
      countCell(c, numPyr, regionsPyramids);
    }
    else if( m_vtkMesh->GetCellType(c) == VTK_TRIANGLE || m_vtkMesh->GetCellType(c) == VTK_QUAD )
    {
      if (attributeDataArray != nullptr)
      {
        surfaces.insert(attributeDataArray->GetValue(c));
      }
    }
  }


  // Communicate all the surfaces attributes to the ranks
  std::vector< int > const surfaceVector(surfaces.cbegin(), surfaces.cend());
  int numSurfaceId = surfaceVector.size();
  int totalNbSurfaceId = 0;
  MpiWrapper::allReduce( &numSurfaceId, &totalNbSurfaceId, 1, MPI_SUM,MPI_COMM_GEOSX );
  allSurfaces.resize(totalNbSurfaceId);
  MpiWrapper::allgather<int, int>(surfaceVector.data(),  surfaceVector.size(), allSurfaces.data(),surfaceVector.size(),MPI_COMM_GEOSX);

  std::sort(allSurfaces.begin(), allSurfaces.end());
  auto last = std::unique( allSurfaces.begin(), allSurfaces.end());
  allSurfaces.erase(last, allSurfaces.end());
}

std::vector< vtkDataArray * > VTKMeshGenerator::findArrayToBeImported()
{
  std::vector< vtkDataArray * > arrayToBeImported;
  vtkCellData * cellData = m_vtkMesh->GetCellData();
  for( int i = 0; i < cellData->GetNumberOfArrays(); i++)
  {
    vtkAbstractArray * curArray = cellData->GetAbstractArray(i);
    // We support only float and double type. Good assumption?
    int dataType = curArray->GetDataType();
    if( !curArray->IsNumeric() )
    {
      GEOSX_LOG_RANK_0( curArray->GetName() << " is not a numeric array");
      continue;
    }
    if( dataType != VTK_FLOAT && dataType != VTK_DOUBLE )
    {
      GEOSX_LOG_RANK_0( "Underlying data Type " << curArray->GetDataTypeAsString() << " for array "
                        << curArray->GetName() << " is not supported by GEOSX and will be ignored "
                        << "(Only double and float are supported)");
      continue;
    }
    GEOSX_LOG_RANK_0("Importing array " << curArray->GetName() );
    arrayToBeImported.push_back( vtkDataArray::SafeDownCast(curArray) );
  }
  return arrayToBeImported;
}

void VTKMeshGenerator::writeCellBlocks( DomainPartition& domain, localIndex numHex, localIndex numTet, localIndex numWedge, localIndex numPyr,
                                        std::map<int,localIndex> & regionsHex, std::map<int,localIndex> & regionsTetra,
                                        std::map<int,localIndex> & regionsWedges, std::map<int,localIndex> & regionsPyramids,
                                        std::set< int > & regions, std::vector< vtkDataArray * > & arraysToBeImported)
{
  CellBlockManager & cellBlockManager = domain.getGroup< CellBlockManager >( keys::cellManager );

  if( !regions.empty() )
  {
    for (int region : regions)
    {
      std::map< int, localIndex>::iterator it;
      it = regionsHex.find(region);
      if( it != regionsHex.end())
      {
        writeCellBlock("hexahedron", it->second, it->first, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported);
      }
      it = regionsTetra.find(region);
      if( it != regionsTetra.end())
      {
        writeCellBlock("tetrahedron", it->second, it->first, VTK_TETRA, cellBlockManager, arraysToBeImported);
      }
      it = regionsWedges.find(region);
      if( it != regionsWedges.end())
      {
        writeCellBlock("wedges", it->second, it->first, VTK_WEDGE, cellBlockManager, arraysToBeImported);
      }
      it = regionsPyramids.find(region);
      if( it != regionsPyramids.end())
        writeCellBlock("pyramids", it->second, it->first, VTK_PYRAMID, cellBlockManager, arraysToBeImported);
    }
  }
  else
  {
      writeCellBlock("hexahedron", numHex, -1, VTK_HEXAHEDRON, cellBlockManager, arraysToBeImported);
      writeCellBlock("tetrahedron", numTet, -1, VTK_TETRA, cellBlockManager, arraysToBeImported);
      writeCellBlock("wedges", numWedge, -1, VTK_WEDGE, cellBlockManager, arraysToBeImported);
      writeCellBlock("pyramids", numPyr, -1, VTK_PYRAMID, cellBlockManager, arraysToBeImported);
  }
}


void VTKMeshGenerator::writeSurfaces( NodeManager & nodeManager, std::vector<int> const & allSurfaces )
{
  Group & nodeSets = nodeManager.sets();
  vtkIntArray * attributeDataArray = getAttributeDataArray();
  for( int surface : allSurfaces)
  {
    SortedArray< localIndex > & curNodeSet  = nodeSets.registerWrapper< SortedArray< localIndex > >( std::to_string( surface) ).reference();
    for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = m_vtkMesh->GetCell(c);
      if( m_vtkMesh->GetCellType(c) == VTK_TRIANGLE ||  m_vtkMesh->GetCellType(c) == VTK_QUAD)
      {
        if( attributeDataArray->GetValue(c) == surface )
        {
          for( localIndex v = 0; v < currentCell->GetNumberOfPoints(); v++)
          {
            curNodeSet.insert(currentCell->GetPointId(v) );
          }
        } 
      }
    }
  }
}

localIndex VTKMeshGenerator::getNumberOfPoints( int cellType )
{
  if( cellType == VTK_HEXAHEDRON )
  {
    return 8;
  }
  else if( cellType == VTK_TETRA )
  {
    return 4;
  }
  else if( cellType == VTK_WEDGE )
  {
    return 6;
  }
  else if( cellType == VTK_PYRAMID )
  {
    return 5;
  }
  else
  {
    GEOSX_ERROR( cellType << " is not a recognized cell type to be used with the VTKMeshGenerator");
    return -1;
  }

}

void VTKMeshGenerator::writeHexahedronVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal  )
{
  localIndex cellCount = 0;
  vtkIntArray * attributeDataArray =  getAttributeDataArray();
  vtkIdTypeArray * globalCellIdsDataArray =getCellGlobalIdDataArray();
  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = m_vtkMesh->GetCell(c);
    if( currentCell->GetCellType() == VTK_HEXAHEDRON )
    {
      if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(4);
		  cellToVertex[cellCount][5] = currentCell->GetPointId(5);
		  cellToVertex[cellCount][6] = currentCell->GetPointId(7);
		  cellToVertex[cellCount][7] = currentCell->GetPointId(6);
      localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
      cellCount++;
    }
  }
}

void VTKMeshGenerator::writeWedgeVertices( CellBlock::NodeMapType & cellToVertex, int region_id,  arrayView1d< globalIndex > const & localToGlobal )
{
  localIndex cellCount = 0;
  vtkIntArray * attributeDataArray =  getAttributeDataArray();
  vtkIdTypeArray * globalCellIdsDataArray = getCellGlobalIdDataArray();
  for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
  {
    vtkCell * currentCell = m_vtkMesh->GetCell(c);
    if( currentCell->GetCellType() == VTK_WEDGE )
    {
      if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
		  cellToVertex[cellCount][0] = currentCell->GetPointId(0);
		  cellToVertex[cellCount][1] = currentCell->GetPointId(3);
		  cellToVertex[cellCount][2] = currentCell->GetPointId(1);
		  cellToVertex[cellCount][3] = currentCell->GetPointId(4);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(2);
		  cellToVertex[cellCount][4] = currentCell->GetPointId(5);
      localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
      cellCount++;
    }
  }
}
void VTKMeshGenerator::writeCellBlock( string const & name, localIndex numCells, int region_id, int cellType,
                                       CellBlockManager & cellBlockManager,
                                       std::vector< vtkDataArray * > const & arraysTobeImported )
{

  if( numCells > 0)
  {
    assert(region_id >= -1 );
    string const cellBlockName = region_id != -1  ? std::to_string(region_id) + "_" + name : name;

    GEOSX_LOG_RANK_0( "Writing " << numCells << " cells to CellBlock " << cellBlockName);
    CellBlock & cellBlock = cellBlockManager.getGroup( keys::cellBlocks ).registerGroup< CellBlock >( cellBlockName );
    localIndex numPointsPerCell = getNumberOfPoints( cellType );
    cellBlock.setElementType( "C3D" + std::to_string(numPointsPerCell));
    CellBlock::NodeMapType & cellToVertex = cellBlock.nodeList();
    cellBlock.resize( numCells );
    arrayView1d< globalIndex > const & localToGlobal = cellBlock.localToGlobalMap();

    vtkIntArray * attributeDataArray =  getAttributeDataArray();


    /// Writing properties
    for( vtkDataArray * array : arraysTobeImported )
    {
      int dimension = array->GetNumberOfComponents();
      if( dimension == 1 )
      {
        real64_array & property = cellBlock.addProperty< real64_array >( array->GetName() );
        localIndex cellCount = 0;
        for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = m_vtkMesh->GetCell(c);
          if( currentCell->GetCellType() == cellType )
          {
            if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
            property[cellCount++] = array->GetTuple1(c);
          }
        }
      }
      else
      {
        array2d< real64 > & property = cellBlock.addProperty< array2d< real64 > >( array->GetName() );
        property.resizeDimension< 1 >( dimension );
        localIndex cellCount = 0;
        for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
        {
          vtkCell * currentCell = m_vtkMesh->GetCell(c);
          if( currentCell->GetCellType() == cellType )
          {
            if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id ) continue;
            for(int i = 0; i < dimension; i++ )
            {
              property(cellCount,i) = array->GetComponent(c,i);
            }
            cellCount++;
          }
        }
      }
    }

    /// Writing connectivity and Local to Global
    cellToVertex.resize( numCells, numPointsPerCell );
    if( cellType == VTK_HEXAHEDRON ) // Special case for hexahedron because of the ordering
    {
      writeHexahedronVertices( cellToVertex, region_id, localToGlobal );
      return;
    }
    if( cellType == VTK_WEDGE ) // Special case for hexahedron because of the ordering
    {
      writeWedgeVertices( cellToVertex, region_id, localToGlobal );
      return;
    }
    vtkIdTypeArray * globalCellIdsDataArray = getCellGlobalIdDataArray();
    localIndex cellCount = 0;
    for( vtkIdType c = 0; c < m_vtkMesh->GetNumberOfCells(); c++)
    {
      vtkCell * currentCell = m_vtkMesh->GetCell(c);
      if( currentCell->GetCellType() == cellType )
      {
        if( attributeDataArray != nullptr && attributeDataArray->GetValue(c) != region_id) continue;
        for( localIndex v = 0; v < numPointsPerCell; v++)
        {
          cellToVertex[cellCount][v] = currentCell->GetPointId(v);
        }
        localToGlobal[cellCount] = globalCellIdsDataArray->GetValue(c);
        cellCount++;
      }
    }
  }
}
REGISTER_CATALOG_ENTRY( MeshGeneratorBase, VTKMeshGenerator, string const &, Group * const )
}
