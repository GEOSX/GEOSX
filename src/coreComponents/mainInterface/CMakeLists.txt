#
# Specify all headers
#
set(mainInterface_headers
    ProblemManager.hpp
    initialization.hpp
    GeosxState.hpp
   )


#
# Specify all sources
#

set(mainInterface_sources
    ProblemManager.cpp
    initialization.cpp
    GeosxState.cpp
   )

set( dependencyList physicsSolvers discretizationMethods fieldSpecification linearAlgebra events fileIO optionparser )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  mainInterface
                 SOURCES               ${mainInterface_sources}
                 HEADERS               ${mainInterface_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( mainInterface PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX mainInterface )
