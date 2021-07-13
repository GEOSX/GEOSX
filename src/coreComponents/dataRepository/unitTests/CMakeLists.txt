#
# Specify list of tests
#

set( dataRepository_tests
     testDefaultValue.cpp
     testWrapper.cpp
     testXmlWrapper.cpp
     testBufferOps.cpp
   )

set( dependencyList gtest dataRepository )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${dataRepository_tests})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList}
                        )

    blt_add_test( NAME ${test_name}
                  COMMAND ${test_name} )
endforeach()
