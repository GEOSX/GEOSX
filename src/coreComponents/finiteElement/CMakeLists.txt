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
     elementFormulations/LagrangeBasis3GL.hpp
     elementFormulations/LagrangeBasis5GL.hpp
     elementFormulations/Q3_Hexahedron_Lagrange_GaussLobatto.hpp
     elementFormulations/Q5_Hexahedron_Lagrange_GaussLobatto.hpp
   )
#
# Specify all sources
#
set( finiteElement_sources
     FiniteElementDiscretization.cpp
     FiniteElementDiscretizationManager.cpp
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
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( finiteElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX finiteElement )

add_subdirectory( unitTests )

