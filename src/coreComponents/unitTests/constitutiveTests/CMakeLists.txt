#
# Specify list of tests
#

set( gtest_geosx_tests
     testDamage.cpp
     testRelPerm.cpp
     testCapillaryPressure.cpp
   )

set( gtest_triaxial_xmls
     testTriaxial_elasticIsotropic.xml
     testTriaxial_druckerPragerExtended.xml
   )

set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core )
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testMultiFluid.cpp
          testCO2BrinePVTModels.cpp )

    set( dependencyList ${dependencyList} PVTPackage )
endif()

if( ENABLE_CUDA )
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
                  COMMAND ${test_name}
                )

endforeach()


#
# Add triaxial xml based tests
#

blt_add_executable( NAME testTriaxial
                    SOURCES testTriaxial.cpp
                    OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                    DEPENDS_ON ${dependencyList} "-lz"
                  )

foreach(test ${gtest_triaxial_xmls})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_test( NAME ${test_name} 
                  COMMAND testTriaxial -i ${CMAKE_CURRENT_LIST_DIR}/${test}
                  )
endforeach()
