# Specify all headers
set( mesh_headers
     AggregateElementSubRegion.hpp
     BufferOps.hpp
     CellBlockManager.hpp
     CellBlock.hpp
     CellElementRegion.hpp
     CellElementSubRegion.hpp
     DomainPartition.hpp
     EdgeManager.hpp
     ElementSubRegionBase.hpp
     ElementRegionBase.hpp
     ElementRegionManager.hpp
     EmbeddedSurfaceSubRegion.hpp
     ExtrinsicMeshData.hpp
     FaceElementSubRegion.hpp
     FaceManager.hpp
     InterObjectRelation.hpp
     MeshBody.hpp
     MeshLevel.hpp
     MeshManager.hpp
     ObjectManagerBase.hpp
     Perforation.hpp
     PerforationData.hpp
     NodeManager.hpp
     EmbeddedSurfaceNodeManager.hpp
     SurfaceElementRegion.hpp
     SurfaceElementSubRegion.hpp
     ToElementRelation.hpp
     WellElementSubRegion.hpp
     WellElementRegion.hpp
     generators/MeshGeneratorBase.hpp
     generators/InternalMeshGenerator.hpp
     generators/InternalWellGenerator.hpp
     generators/InternalWellboreGenerator.hpp
     mpiCommunications/CommID.hpp
     mpiCommunications/CommunicationTools.hpp
     mpiCommunications/MPI_iCommData.hpp
     mpiCommunications/NeighborCommunicator.hpp
     mpiCommunications/PartitionBase.hpp
     mpiCommunications/SpatialPartition.hpp
     mpiCommunications/NeighborData.hpp
     simpleGeometricObjects/GeometricObjectManager.hpp
     simpleGeometricObjects/SimpleGeometricObjectBase.hpp
     simpleGeometricObjects/Box.hpp
     simpleGeometricObjects/ThickPlane.hpp
     simpleGeometricObjects/BoundedPlane.hpp
     utilities/ComputationalGeometry.hpp
     utilities/MeshMapUtilities.hpp
     utilities/MeshUtilities.hpp
     utilities/StructuredGridUtilities.hpp
  )

# Specify all sources
set(mesh_sources
    AggregateElementSubRegion.cpp
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellElementRegion.cpp
    CellElementSubRegion.cpp
    DomainPartition.cpp
    EdgeManager.cpp
    ElementSubRegionBase.cpp
    ElementRegionBase.cpp
    ElementRegionManager.cpp
    EmbeddedSurfaceSubRegion.cpp
    FaceElementSubRegion.cpp
    FaceManager.cpp
    MeshBody.cpp
    MeshLevel.cpp
    MeshManager.cpp
    NodeManager.cpp
    EmbeddedSurfaceNodeManager.cpp
    ObjectManagerBase.cpp
    PerforationData.cpp
    Perforation.cpp
    SurfaceElementRegion.cpp
    SurfaceElementSubRegion.cpp
    ToElementRelation.cpp
    WellElementRegion.cpp
    WellElementSubRegion.cpp
    generators/InternalMeshGenerator.cpp
    generators/InternalWellGenerator.cpp
    generators/InternalWellboreGenerator.cpp
    generators/MeshGeneratorBase.cpp
    mpiCommunications/CommID.cpp
    mpiCommunications/CommunicationTools.cpp
    mpiCommunications/MPI_iCommData.cpp    
    mpiCommunications/NeighborCommunicator.cpp
    mpiCommunications/PartitionBase.cpp
    mpiCommunications/SpatialPartition.cpp
    simpleGeometricObjects/GeometricObjectManager.cpp
    simpleGeometricObjects/SimpleGeometricObjectBase.cpp
    simpleGeometricObjects/Box.cpp
    simpleGeometricObjects/Cylinder.cpp
    simpleGeometricObjects/ThickPlane.cpp
    simpleGeometricObjects/BoundedPlane.cpp
    utilities/ComputationalGeometry.cpp
    utilities/MeshUtilities.cpp
   )

set(dependencyList schema dataRepository constitutive metis)

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_PAMELA )
    message(STATUS "Adding PAMELAMeshGenerator sources and headers")
    list( APPEND mesh_headers generators/PAMELAMeshGenerator.hpp )
    list( APPEND mesh_sources generators/PAMELAMeshGenerator.cpp)
    add_subdirectory(PAMELA)
    list( APPEND dependencyList PAMELA )
endif()

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX mesh 
                       EXCLUDES PAMELA )
