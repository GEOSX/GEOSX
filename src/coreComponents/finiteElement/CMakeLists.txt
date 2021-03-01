#
# Specify all headers
#
set( finiteElement_headers
     FiniteElementDiscretization.hpp
     FiniteElementDiscretizationManager.hpp
     FiniteElementDispatch.hpp
     Kinematics.h
     kernelInterface/KernelBase.hpp
     kernelInterface/ImplicitKernelBase.hpp
     FiniteElementDiscretizationManager.hpp
     FiniteElementDispatch.hpp
     elementFormulations/FiniteElementBase.hpp
     elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp
     elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp
     elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp
     elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp
     elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp
     elementFormulations/LagrangeBasis1.hpp
     elementFormulations/LagrangeBasis2.hpp
     ${CMAKE_BINARY_DIR}/include/kernelJITCompileCommands.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     FiniteElementDiscretizationManager.cpp
     kernelInterface/kernelJIT.cpp
   )

if( BUILD_OBJ_LIBS)
  set( dependencyList dataRepository codingUtilities)
else()
  set( dependencyList common)
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  finiteElement
                 SOURCES               ${finiteElement_sources}
                 HEADERS               ${finiteElement_headers}
                 DEPENDS_ON            ${dependencyList} jitti
                 OBJECT                ${buildAsObj}
               )

target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

set( linkDirectories "" )
set( linkLibraries "" )

if( ENABLE_CUDA )
    list( APPEND linkDirectories ${CUDA_TOOLKIT_ROOT_DIR}/lib64 )
    list( APPEND linkLibraries cudart_static cudadevrt )
endif()

add_custom_command( OUTPUT ${CMAKE_BINARY_DIR}/include/kernelJITCompileCommands.hpp
                    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                    COMMAND python ${CMAKE_CURRENT_LIST_DIR}/../LvArray/src/jitti/generateCompileCommandsHeader.py
                                   ${CMAKE_BINARY_DIR}/compile_commands.json
                                   --cpp ${CMAKE_CURRENT_LIST_DIR}/kernelInterface/kernelJIT.cpp
                                   --hpp ${CMAKE_BINARY_DIR}/include/kernelJITCompileCommands.hpp
                                   --include ${CMAKE_BINARY_DIR}/include
                                   --linker ${CMAKE_CXX_COMPILER}
                                   --linkDirectories ${linkDirectories}
                                   --linkLibraries ${linkLibraries} )

geosx_add_code_checks( PREFIX finiteElement )

add_subdirectory( unitTests )

