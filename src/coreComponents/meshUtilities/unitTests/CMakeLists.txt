#
# Specify list of tests
#

if(ENABLE_PAMELA)
set( gtest_geosx_tests
     ${gtest_geosx_tests}
    testMeshImport.cpp
   )

set( GMSH_FILE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/toy_model.msh )
set( ECLIPSE_FILE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/toy_model.GRDECL )
configure_file( ${CMAKE_CURRENT_SOURCE_DIR}/meshFileNames.hpp.in ${CMAKE_BINARY_DIR}/include/tests/meshFileNames.hpp )
endif(ENABLE_PAMELA)

set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if ( ENABLE_MPI )
  set ( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_OPENMP )
    set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
            SOURCES ${test}
            OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
            DEPENDS_ON ${dependencyList}
            )

    blt_add_test( NAME ${test_name}
            COMMAND ${test_name} ${CMAKE_CURRENT_LIST_DIR}
            )

endforeach()
