#
# Specify all headers
#
set(fieldSpecification_headers
    FieldSpecificationBase.hpp
    FieldSpecificationManager.hpp
    DirichletBoundaryCondition.hpp
    SourceFluxBoundaryCondition.hpp
    TractionBoundaryCondition.hpp
   )


#
# Specify all sources
#

set(fieldSpecification_sources
    FieldSpecificationBase.cpp
    FieldSpecificationManager.cpp
    DirichletBoundaryCondition.cpp
    SourceFluxBoundaryCondition.cpp
    TractionBoundaryCondition.cpp
   )



set( dependencyList functions RAJA linearAlgebra conduit )

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if ( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

if ( ENABLE_MKL )
  set( dependencyList ${dependencyList} mkl )
endif()

blt_add_library( NAME                  fieldSpecification
                 SOURCES               ${fieldSpecification_sources}
                 HEADERS               ${fieldSpecification_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( fieldSpecification PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX fieldSpecification )
