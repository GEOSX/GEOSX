#
# Specify all headers
#
set(functions_headers
    FunctionBase.hpp
    TableFunction.hpp
    FunctionManager.hpp
   )


#
# Specify all sources
#

set(functions_sources
    FunctionBase.cpp
    TableFunction.cpp
    FunctionManager.cpp
   )

if( ENABLE_MATHPRESSO )
  list( APPEND functions_headers
        SymbolicFunction.hpp
        CompositeFunction.hpp )
  list( APPEND functions_sources
        SymbolicFunction.cpp
        CompositeFunction.cpp )
endif()


set( dependencyList dataRepository )

if( ENABLE_MATHPRESSO )
   set( dependencyList ${dependencyList} mathpresso )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  functions
                 SOURCES               ${functions_sources}
                 HEADERS               ${functions_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( functions PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX functions )
