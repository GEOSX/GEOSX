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
     FluxStencil.hpp
     HybridMimeticDiscretization.hpp
     MimeticInnerProductDispatch.hpp     
     mimeticInnerProducts/MimeticInnerProductBase.hpp
     mimeticInnerProducts/QuasiTPFAInnerProduct.hpp
     mimeticInnerProducts/QuasiRTInnerProduct.hpp
     mimeticInnerProducts/TPFAInnerProduct.hpp
     mimeticInnerProducts/SimpleInnerProduct.hpp
     mimeticInnerProducts/BdVLMInnerProduct.hpp
     mimeticInnerProducts/MimeticInnerProductHelpers.hpp
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
     HybridMimeticDiscretization.cpp
   )

set( dependencyList mesh codingUtilities fieldSpecification )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  finiteVolume
                 SOURCES               ${finiteVolume_sources}
                 HEADERS               ${finiteVolume_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )
               
target_include_directories( finiteVolume PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX finiteVolume )
