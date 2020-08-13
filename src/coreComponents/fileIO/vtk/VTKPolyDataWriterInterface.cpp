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

// Source includes
#include "VTKPolyDataWriterInterface.hpp"
#include "dataRepository/Group.hpp"

// TPL includes
#include <vtkUnstructuredGrid.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkExtentTranslator.h>

// System includes
#include <unordered_set>
#include <sys/stat.h>

namespace geosx
{
namespace vtk
{
/*!
 * @brief Get the DataSet file path
 * @details the DataSet file path is the path to the .vtu mesh file
 * @param[in] er The ElementRegion
 * @param[in] time the time-step
 * @param[in] rank the rank to be written
 */
string
GetDataSetFilePath( ElementRegionBase const & er, double time, int rank )
{
  return std::to_string( time ) + "/" +
    stringutilities::PadValue( rank,
                               std::to_string( MpiWrapper::Comm_size() ).size() ) +
    "_" + er.getName() + ".vtu";
}

const std::map< string, int > geosx2VTKCellTypes = { { "C3D4", VTK_TETRA },
                                                     { "C3D8", VTK_HEXAHEDRON },
                                                     { "C3D6", VTK_WEDGE },
                                                     { "C3D5", VTK_PYRAMID } };

/*!
 * @brief Gets the VTK cell identifier
 * @param[in] elementType the type of the element (using the abaqus nomenclature)
 * @return the VTK cell identifier
 */
int
ToVTKCellType( const string & elementType )
{
  int vtkIdentifier = VTK_EMPTY_CELL;
  try
  {
    vtkIdentifier = geosx2VTKCellTypes.at( elementType );
  }
  catch( const std::out_of_range & outOfRange )
  {
    GEOSX_ERROR( "Element type " << elementType
                                 << " not recognized for VTK output " );
  }
  return vtkIdentifier;
}

VTKPolyDataWriterInterface::VTKPolyDataWriterInterface( string const & outputName ) :
  m_outputFolder( outputName ),
  m_pvd( outputName + ".pvd" ),
  m_previousCycle( -1 )
{
  int const mpiRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0733;
    int errorCode = mkdir( outputName.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << outputName
                               << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create main directory " << outputName
                                                      << " for the VTK output" );
      }
    }
  }
  MpiWrapper::Barrier();
}

vtkSmartPointer< vtkPoints >
VTKPolyDataWriterInterface::GetVTKPoints(
  NodeManager const & nodeManager ) const
{
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( nodeManager.size() );
  auto connectivity = nodeManager.referencePosition();
  for( localIndex v = 0; v < nodeManager.size(); v++ )
  {
    points->SetPoint( v, connectivity[v][0], connectivity[v][1], connectivity[v][2] );
  }
  return points;
}

std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::GetWell( WellElementSubRegion const & esr,
                                     NodeManager const & nodeManager ) const
{
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( esr.size() + 1 );
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.size() );
  localIndex numberOfNodesPerElement = esr.numNodesPerElement();
  GEOSX_ERROR_IF_NE( numberOfNodesPerElement, 2 );
  std::vector< vtkIdType > connectivity( numberOfNodesPerElement );
  for( localIndex edge = 0; edge < esr.size(); edge++ )
  {
    localIndex firstPoint = esr.nodeList()[edge][0];
    auto point = nodeManager.referencePosition()[firstPoint];
    points->SetPoint( edge, point[0], point[1], point[2] );
    connectivity[0] = edge;
    connectivity[1] = edge +
      1;  // We can do that because of the pattern in which the wells are stored
    cellsArray->InsertNextCell( numberOfNodesPerElement, connectivity.data() );
  }
  if( esr.size() > 0 )
  {
    localIndex lastPoint = esr.nodeList()[esr.size() - 1][1];
    auto point = nodeManager.referencePosition()[lastPoint];
    points->SetPoint( esr.size(), point[0], point[1], point[2] );
  }
  return std::make_pair( points, cellsArray );
}
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::GetSurface( FaceElementSubRegion const & esr,
                                        NodeManager const & nodeManager ) const
{
  // Get unique node set composing the surface
  auto & nodeListPerElement = esr.nodeList();
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( esr.size() );
  std::unordered_map< localIndex, localIndex > geosx2VTKIndexing;
  geosx2VTKIndexing.reserve( esr.size() * esr.numNodesPerElement() );
  localIndex nodeIndexInVTK = 0;
  std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
  std::vector< int > vtkOrdering = esr.getVTKNodeOrdering();
  for( localIndex ei = 0; ei < esr.size(); ei++ )
  {
    auto & elem = nodeListPerElement[ei];
    for( localIndex i = 0; i < elem.size(); i++ )
    {
      auto const & VTKIndexPos = geosx2VTKIndexing.find( elem[vtkOrdering[i]] );
      if( VTKIndexPos == geosx2VTKIndexing.end() )
      {
        connectivity[i] = geosx2VTKIndexing[elem[vtkOrdering[i]]] =
          nodeIndexInVTK++;
      }
      else
      {
        connectivity[i] = VTKIndexPos->second;
      }
    }
    cellsArray->InsertNextCell( elem.size(), connectivity.data() );
  }

  vtkSmartPointer< vtkPoints > points = vtkPoints::New();
  points->SetNumberOfPoints( geosx2VTKIndexing.size() );
  for( auto nodeIndex : geosx2VTKIndexing )
  {
    auto point = nodeManager.referencePosition()[nodeIndex.first];
    points->SetPoint( nodeIndex.second, point[0], point[1], point[2] );
  }
  return std::make_pair( points, cellsArray );
}
std::pair< vtkSmartPointer< vtkPoints >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::GetEmbeddedSurface(
  EmbeddedSurfaceSubRegion const & esr,
  ElementRegionManager const & elemManager,
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  vtkSmartPointer< vtkPoints > points = vtkPoints::New();

  array1d< R1Tensor > intersectionPoints;
  array1d< localIndex > connectivityList;
  array1d< int > offSet, typesList;
  // Get "nodes" relative to the fracture subregion
  esr.getIntersectionPoints( nodeManager,
                             edgeManager,
                             elemManager,
                             intersectionPoints,
                             connectivityList,
                             offSet );

  points->SetNumberOfPoints( intersectionPoints.size() );
  for( localIndex pointIndex = 0; pointIndex < intersectionPoints.size();
       pointIndex++ )
  {
    points->SetPoint( pointIndex,
                      intersectionPoints[pointIndex][0],
                      intersectionPoints[pointIndex][1],
                      intersectionPoints[pointIndex][2] );
  }

  cellsArray->SetNumberOfCells( esr.size() );
  for( localIndex cellIndex = 0; cellIndex < esr.size(); cellIndex++ )
  {
    std::vector< vtkIdType > connectivity( offSet[cellIndex + 1] - offSet[cellIndex] );
    for( localIndex nodeCellIndex = 0;
         LvArray::integerConversion< size_t >( nodeCellIndex ) < connectivity.size();
         nodeCellIndex++ )
    {
      connectivity[nodeCellIndex] =
        connectivityList[offSet[cellIndex] + nodeCellIndex];
    }
    cellsArray->InsertNextCell( connectivity.size(), connectivity.data() );
  }
  return std::make_pair( points, cellsArray );
}

std::pair< std::vector< int >, vtkSmartPointer< vtkCellArray > >
VTKPolyDataWriterInterface::GetVTKCells( CellElementRegion const & er ) const
{
  vtkSmartPointer< vtkCellArray > cellsArray = vtkCellArray::New();
  cellsArray->SetNumberOfCells( er.getNumberOfElements< CellElementRegion >() );
  std::vector< int > cellType;
  cellType.reserve( er.getNumberOfElements< CellElementRegion >() );
  er.forElementSubRegions< CellElementSubRegion >(
    [&]( CellElementSubRegion const & esr ) {
      std::vector< vtkIdType > connectivity( esr.numNodesPerElement() );
      std::vector< int > vtkOrdering = esr.getVTKNodeOrdering();
      int vtkCellType = ToVTKCellType( esr.GetElementTypeString() );
      for( localIndex c = 0; c < esr.size(); c++ )
      {
        for( std::size_t i = 0; i < connectivity.size(); i++ )
        {
          connectivity[i] = esr.nodeList( c, vtkOrdering[i] );
        }
        cellType.push_back( vtkCellType );
        cellsArray->InsertNextCell( esr.numNodesPerElement(), connectivity.data() );
      }
    } );
  return std::make_pair( cellType, cellsArray );
}

void
VTKPolyDataWriterInterface::WriteField( WrapperBase const & wrapperBase,
                                        vtkSmartPointer< VTKGEOSXData > data,
                                        localIndex size,
                                        localIndex & count ) const
{
  std::type_info const & typeID = wrapperBase.get_typeid();
  if( typeID == typeid( r1_array ) )
  {
    // We need a special case for the R1 array
    // Because those array are not stored the same way than the classical array2d which is
    // preferable to store a vector on each element.
    data->SetNumberOfComponents( 3 );
  }
  rtTypes::ApplyArrayTypeLambda2(
    rtTypes::typeID( typeID ),
    true,
    [&]( auto array, auto GEOSX_UNUSED_PARAM( Type ) ) -> void {
      typedef decltype( array ) arrayType;
      Wrapper< arrayType > const & wrapperT = Wrapper< arrayType >::cast( wrapperBase );
      typename arrayType::ViewTypeConst const & sourceArray = wrapperT.reference();
      if( typeID != typeid( r1_array ) )
      {
        integer nbOfComponents = 1;
        for( localIndex i = 1; i < arrayType::NDIM; i++ )
        {
          nbOfComponents = nbOfComponents * sourceArray.size( i );
        }
        data->SetNumberOfComponents( nbOfComponents );
      }
      for( localIndex i = 0; i < size; i++ )
      {
        LvArray::forValuesInSlice( sourceArray[i], [&]( auto const & value ) {
          data->CustomInsertValue( count++, value );
        } );
      }
    } );
}

void
VTKPolyDataWriterInterface::WriteNodeFields(
  vtkSmartPointer< vtkPointData > const pointdata,
  NodeManager const & nodeManager ) const
{
  for( auto const & wrapperIter : nodeManager.wrappers() )
  {
    auto const & wrapper = *wrapperIter.second;
    if( wrapper.getPlotLevel() <= m_plotLevel )
    {
      vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
      data->SetNumberOfValues( nodeManager.size() );
      data->SetName( wrapper.getName().c_str() );
      ;
      localIndex count = 0;
      WriteField( wrapper, data, nodeManager.size(), count );
      pointdata->AddArray( data );
    }
  }
}

template< class SUBREGION >
void
VTKPolyDataWriterInterface::WriteElementFields(
  vtkSmartPointer< vtkCellData > const celldata,
  ElementRegionBase const & er ) const
{
  std::unordered_set< string > allFields;
  er.forElementSubRegions< SUBREGION >( [&]( auto const & esr ) {
    for( auto const & wrapperIter : esr.wrappers() )
    {
      auto const * const wrapper = wrapperIter.second;
      if( wrapper->getPlotLevel() <= m_plotLevel )
      {
        allFields.insert( wrapperIter.first );
      }
    }
  } );

  for( auto const & field : allFields )
  {
    vtkSmartPointer< VTKGEOSXData > data = VTKGEOSXData::New();
    data->SetNumberOfValues( er.getNumberOfElements< SUBREGION >() );
    data->SetName( field.c_str() );

    localIndex count = 0;
    er.forElementSubRegions< SUBREGION >( [&]( auto const & esr ) {
      auto const & wrapper = *esr.getWrapperBase( field );
      WriteField( wrapper, data, esr.size(), count );
    } );
    celldata->AddArray( data );
  }
}
void
VTKPolyDataWriterInterface::WriteCellElementRegions(
  real64 time,
  ElementRegionManager const & elemManager,
  NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< CellElementRegion >(
    [&]( CellElementRegion const & er ) -> void {
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKPoints = GetVTKPoints( nodeManager );
      ug->SetPoints( VTKPoints );
      auto VTKCells = GetVTKCells( er );
      ug->SetCells( VTKCells.first.data(), VTKCells.second );
      WriteElementFields< CellElementSubRegion >( ug->GetCellData(), er );
      WriteNodeFields( ug->GetPointData(), nodeManager );
      WriteUnstructuredGrid( ug, time, er.getName() );
    } );
}

void
VTKPolyDataWriterInterface::WriteWellElementRegions(
  real64 time,
  ElementRegionManager const & elemManager,
  NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< WellElementRegion >(
    [&]( WellElementRegion const & er ) -> void {
      auto esr = er.GetSubRegion( 0 )->group_cast< WellElementSubRegion const * >();
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKWell = GetWell( *esr, nodeManager );
      ug->SetPoints( VTKWell.first );
      ug->SetCells( VTK_LINE, VTKWell.second );
      WriteElementFields< WellElementSubRegion >( ug->GetCellData(), er );
      WriteUnstructuredGrid( ug, time, er.getName() );
    } );
}

void
VTKPolyDataWriterInterface::WriteFaceElementRegions(
  real64 time,
  ElementRegionManager const & elemManager,
  NodeManager const & nodeManager ) const
{
  elemManager.forElementRegions< FaceElementRegion >(
    [&]( FaceElementRegion const & er ) -> void {
      auto esr = er.GetSubRegion( 0 )->group_cast< FaceElementSubRegion const * >();
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();
      auto VTKSurface = GetSurface( *esr, nodeManager );
      ug->SetPoints( VTKSurface.first );
      if( esr->numNodesPerElement() == 8 )
      {
        ug->SetCells( VTK_HEXAHEDRON, VTKSurface.second );
      }
      else if( esr->numNodesPerElement() == 6 )
      {
        ug->SetCells( VTK_WEDGE, VTKSurface.second );
      }
      else
      {
        GEOSX_ERROR( "Elements with "
                     << esr->numNodesPerElement() << " nodes can't be output "
                     << "in the FaceElementRegion " << er.getName() );
      }
      WriteElementFields< FaceElementSubRegion >( ug->GetCellData(), er );
      WriteUnstructuredGrid( ug, time, er.getName() );
    } );
}

void
VTKPolyDataWriterInterface::WriteEmbeddedSurfaceElementRegions(
  real64 time,
  ElementRegionManager const & elemManager,
  NodeManager const & nodeManager,
  EdgeManager const & edgeManager ) const
{
  elemManager.forElementRegions< EmbeddedSurfaceRegion >(
    [&]( EmbeddedSurfaceRegion const & er ) -> void {
      auto esr =
        er.GetSubRegion( 0 )->group_cast< EmbeddedSurfaceSubRegion const * >();
      vtkSmartPointer< vtkUnstructuredGrid > ug = vtkUnstructuredGrid::New();

      auto VTKEmbeddedSurface =
        GetEmbeddedSurface( *esr, elemManager, nodeManager, edgeManager );
      ug->SetPoints( VTKEmbeddedSurface.first );
      ug->SetCells( VTK_POLYGON, VTKEmbeddedSurface.second );

      WriteElementFields< EmbeddedSurfaceSubRegion >( ug->GetCellData(), er );
      WriteUnstructuredGrid( ug, time, er.getName() );
    } );
}

void
VTKPolyDataWriterInterface::WriteVTMFile( real64 time,
                                          ElementRegionManager const & elemManager,
                                          VTKVTMWriter const & vtmWriter ) const
{
  int const mpiRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  int const mpiSize = MpiWrapper::Comm_size( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    auto writeSubBlocks = [&]( ElementRegionBase const & er ) {
      vtmWriter.AddSubBlock( er.getCatalogName(), er.getName() );
      for( int i = 0; i < mpiSize; i++ )
      {
        string const dataSetFile = GetDataSetFilePath( er, time, i );
        vtmWriter.AddDataToSubBlock( er.getCatalogName(),
                                     er.getName(),
                                     dataSetFile,
                                     i );
      }
    };
    // Cells
    vtmWriter.AddBlock( CellElementRegion::CatalogName() );
    elemManager.forElementRegions< CellElementRegion >( writeSubBlocks );

    // Wells
    vtmWriter.AddBlock( WellElementRegion::CatalogName() );
    elemManager.forElementRegions< WellElementRegion >( writeSubBlocks );

    // Surfaces
    vtmWriter.AddBlock( FaceElementRegion::CatalogName() );
    elemManager.forElementRegions< FaceElementRegion >( writeSubBlocks );

    // Embedded Surfaces
    vtmWriter.AddBlock( EmbeddedSurfaceRegion::CatalogName() );
    elemManager.forElementRegions< EmbeddedSurfaceRegion >( writeSubBlocks );

    vtmWriter.Save();
  }
}

void
VTKPolyDataWriterInterface::CreateTimeStepSubFolder( real64 time ) const
{
  int const mpiRank = MpiWrapper::Comm_rank( MPI_COMM_GEOSX );
  if( mpiRank == 0 )
  {
    mode_t mode = 0773;
    string const timeStepSubFolder = GetTimeStepSubFolder( time );
    int errorCode = mkdir( timeStepSubFolder.c_str(), mode );
    if( errorCode == -1 )
    {
      if( errno == EEXIST )
      {
        GEOSX_WARNING( "Path " << timeStepSubFolder
                               << " already exists from a previous simulation" );
      }
      else
      {
        GEOSX_ERROR( "Fail to create the timestep directory "
                     << timeStepSubFolder << " for the VTK output" );
      }
    }
  }
  MpiWrapper::Barrier();
}

void
VTKPolyDataWriterInterface::WriteUnstructuredGrid(
  vtkSmartPointer< vtkUnstructuredGrid > ug,
  double time,
  string const & name ) const
{
  string timeStepSubFolder =
    VTKPolyDataWriterInterface::GetTimeStepSubFolder( time );
  vtkSmartPointer< vtkXMLUnstructuredGridWriter > vtuWriter =
    vtkXMLUnstructuredGridWriter::New();
  vtuWriter->SetInputData( ug );
  string vtuFilePath = timeStepSubFolder + "/" +
    stringutilities::PadValue( MpiWrapper::Comm_rank(),
                               std::to_string( MpiWrapper::Comm_size() ).size() ) +
    "_" + name + ".vtu";
  vtuWriter->SetFileName( vtuFilePath.c_str() );
  if( m_outputMode == VTKOutputMode::BINARY )
  {
    vtuWriter->SetDataModeToBinary();
  }
  else if( m_outputMode == VTKOutputMode::ASCII )
  {
    vtuWriter->SetDataModeToAscii();
  }
  vtuWriter->Write();
}

string
VTKPolyDataWriterInterface::GetTimeStepSubFolder( real64 time ) const
{
  return m_outputFolder + "/" + std::to_string( time );
}

void
VTKPolyDataWriterInterface::Write( real64 time,
                                   integer cycle,
                                   DomainPartition const & domain )
{
  CreateTimeStepSubFolder( time );
  ElementRegionManager const & elemManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getElemManager();
  NodeManager const & nodeManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getNodeManager();
  EdgeManager const & edgeManager =
    *domain.getMeshBody( 0 )->getMeshLevel( 0 )->getEdgeManager();
  WriteCellElementRegions( time, elemManager, nodeManager );
  WriteWellElementRegions( time, elemManager, nodeManager );
  WriteFaceElementRegions( time, elemManager, nodeManager );
  WriteEmbeddedSurfaceElementRegions( time, elemManager, nodeManager, edgeManager );
  string vtmPath = GetTimeStepSubFolder( time ) + ".vtm";
  VTKVTMWriter vtmWriter( vtmPath );
  WriteVTMFile( time, elemManager, vtmWriter );
  if( cycle != m_previousCycle )
  {
    m_pvd.AddData( time, vtmPath );
    m_pvd.Save();
  }
  m_previousCycle = cycle;
}
}  // namespace vtk
}  // namespace geosx
