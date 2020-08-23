#
# Specify all headers
#
set( finiteVolume_headers
     BoundaryStencil.hpp
     CellElementStencilMPFA.hpp
     CellElementStencilTPFA.hpp
     FaceElementStencil.hpp
     FiniteVolumeManager.hpp
     FluxApproximationBase.hpp
     TwoPointFluxApproximation.hpp
     TwoPointFluxApproximationWithGraph.hpp
     FluxStencil.hpp
     HybridFVMInnerProduct.hpp
   )
#
# Specify all sources
#
set( finiteVolume_sources
     BoundaryStencil.cpp
     CellElementStencilMPFA.cpp
     CellElementStencilTPFA.cpp
     FaceElementStencil.cpp
     FiniteVolumeManager.cpp
     FluxApproximationBase.cpp
     TwoPointFluxApproximation.cpp
     TwoPointFluxApproximationWithGraph.cpp
   )

if ( BUILD_OBJ_LIBS )
  set( dependencyList dataRepository codingUtilities managers)
else()
  set( dependencyList common)
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  finiteVolume
                 SOURCES               ${finiteVolume_sources}
                 HEADERS               ${finiteVolume_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )
               
target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX finiteVolume )
