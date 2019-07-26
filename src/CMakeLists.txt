message( "Entering /src/CMakeLists.txt")

cmake_minimum_required(VERSION 3.9)

################################
# GEOSX
################################
project(geosx LANGUAGES C CXX)
set(BLT_CXX_STD "c++14" CACHE STRING "Version of C++ standard")
set( ENABLE_WARNINGS_AS_ERRORS "ON" CACHE PATH "")
set(CMAKE_ENABLE_EXPORTS ON)

include(cmake/thirdparty/setupThirdPartyDirectory.cmake)



################################
# Include blt
################################
include(cmake/blt/SetupBLT.cmake)

################################
# Include standard build system logic and options
################################
include(cmake/CMakeBasics.cmake)
set(CMAKE_DEBUG_POSTFIX "" CACHE STRING "" FORCE)


################################
# Add Thirdparty Builtin Libs
################################
add_subdirectory(thirdparty)


################################
# Add components
################################
#add_subdirectory(components)

if(ENABLE_TRIBOL)
   add_definitions(-DHAVE_TRIBOLCOUPLING)
endif()

set( externalComponentsLinkList "" ) 
add_subdirectory(externalComponents)
add_subdirectory(coreComponents)


################################
# Add docs
################################
if (ENABLE_DOCS)
  add_subdirectory(docs)
endif()


#add_code_check_targets(uncrustify.cfg)


################################
# install targets from build
################################s
install(DIRECTORY ${PROJECT_BINARY_DIR}/bin DESTINATION . OPTIONAL)
install(DIRECTORY ${PROJECT_BINARY_DIR}/include DESTINATION . OPTIONAL)


################################
# Create header of configuration options
################################
include(cmake/GeosxConfig.cmake)


################################
# Add main
################################
set( extraComponentsLinkList "")
if( ENABLE_PYTHON )
  if( PYTHON_DIR )
    set( extraComponentsLinkList ${extraComponentsLinkList} python_interp )
#    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=1)
    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -L${PYTHON_DIR}/lib)
  else()
    MESSAGE(FATAL_ERROR "ENABLE_PYTHON is true, but PYTHON_DIR is not defined.")
  endif()
else()
#  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DGEOSX_USE_PYTHON=0)
endif()

if( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

blt_add_executable(NAME geosx
                   SOURCES main/main.cpp
                   DEPENDS_ON geosx_core sidre
                              ${extraComponentsLinkList}
                              ${externalComponentsLinkList}
                  )

target_include_directories( geosx PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)



add_custom_target( geosx_update_rst_tables
                   ALL
                   COMMAND bash ../scripts/updateRstTables.bash > /dev/null 2>&1
                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                   DEPENDS geosx
                 )

if( ${SPHINX_FOUND} )
  blt_add_sphinx_target( geosx_docs )
endif()

message( "Exiting /src/CMakeLists.txt")
