#
# Specify list of tests
#

set( gtest_geosx_tests
     testDelftEgg.cpp
     testDruckerPrager.cpp
     testLinearElasticIsotropic.cpp
     testLinearElasticAnisotropic.cpp
     testRelPerm.cpp
     testCapillaryPressure.cpp 
   )

set( dependencyList gtest )

if ( GEOSX_BUILD_SHARED_LIBS )
  set (dependencyList ${dependencyList} geosx_core)
else()
  set (dependencyList ${dependencyList} ${geosx_core_libs} )
endif()

if( ENABLE_MPI )
    set( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_PVTPackage )
    set( gtest_geosx_tests ${gtest_geosx_tests} testMultiFluid.cpp )
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
