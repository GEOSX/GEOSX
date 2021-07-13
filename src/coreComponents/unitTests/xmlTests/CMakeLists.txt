#
# Specify list of tests
#

set(gtest_geosx_tests
   testXML.cpp
   )

#
# Add gtest C++ based tests
#
set( dependencyList gtest optionparser)

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

foreach(test ${gtest_geosx_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} "-lz"
                        )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name} -i ${CMAKE_CURRENT_LIST_DIR}/basic_input.xml
                  )
endforeach()
