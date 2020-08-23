# Specify all headers
set(mesh_headers
    AggregateElementSubRegion.hpp
    BufferOps.hpp
    CellBlockManager.hpp
    CellBlock.hpp
    CellElementRegion.hpp
    CellElementSubRegion.hpp
    EdgeManager.hpp
    ElementSubRegionBase.hpp
    ElementRegionBase.hpp
    ElementRegionManager.hpp
    FaceElementSubRegion.hpp
    FaceElementRegion.hpp
    EmbeddedSurfaceSubRegion.hpp
    EmbeddedSurfaceRegion.hpp
    WellElementSubRegion.hpp
    WellElementRegion.hpp
    ExtrinsicMeshData.hpp
    FaceManager.hpp
    InterObjectRelation.hpp
    GraphBase.hpp
    GraphFromText.hpp
    GraphManager.hpp
    Vertice.hpp
    Edge.hpp
    MeshBody.hpp
    MeshLevel.hpp
    NodeManager.hpp
    ToElementRelation.hpp
  )

# Specify all sources
set(mesh_sources
    AggregateElementSubRegion.cpp
    BufferOps.cpp
    CellBlockManager.cpp
    CellBlock.cpp
    CellElementRegion.cpp
    CellElementSubRegion.cpp
    EdgeManager.cpp
    ElementSubRegionBase.cpp
    ElementRegionBase.cpp
    ElementRegionManager.cpp
    FaceElementRegion.cpp
    FaceElementSubRegion.cpp
    EmbeddedSurfaceRegion.cpp
    EmbeddedSurfaceSubRegion.cpp
    WellElementRegion.cpp
    WellElementSubRegion.cpp
    FaceManager.cpp
    GraphBase.cpp
    GraphFromText.cpp
    GraphManager.cpp
    Vertice.cpp
    Edge.cpp
    MeshBody.cpp
    MeshLevel.cpp
    NodeManager.cpp
    ToElementRelation.cpp
   )

if( BUILD_OBJ_LIBS)
  set(dependencyList managers metis)
else()
  set(dependencyList common metis)
endif()


if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  mesh
                 SOURCES               ${mesh_sources}
                 HEADERS               ${mesh_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( mesh PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX mesh )
