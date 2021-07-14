#
# Specify list of tests
#

set( gtest_geosx_tests
     testRelPerm.cpp
     testCapillaryPressure.cpp
   )

set( gtest_triaxial_xmls
     testTriaxial_druckerPragerExtended.xml
     testTriaxial_elasticIsotropic.xml
     testTriaxial_elasticIsotropicPressureDependent.xml
     testTriaxial_modifiedCamClay.xml
     testTriaxial_delftEggUseLinear.xml
     testTriaxial_delftEggLoadPathDry.xml
     testTriaxial_delftEggLoadPathWet.xml
     testTriaxial_modifiedCamClayVolumetric.xml
     testTriaxial_delftEggVolumetric.xml
     testTriaxial_delftEggLoadPathDryUseLinear.xml
     testTriaxial_delftEggLoadPathWetUseLinear.xml
   )

set( dependencyList gtest geosx_core conduit )

if( ENABLE_MPI )
    set( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testMultiFluid.cpp
          testCO2BrinePVTModels.cpp )

    set( dependencyList ${dependencyList} PVTPackage )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
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
