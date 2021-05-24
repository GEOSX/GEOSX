#
# Specify list of tests
#

set( gtest_geosx_tests
     testSinglePhaseBaseKernels.cpp
     testSinglePhaseFVMKernels.cpp     
     testSinglePhaseHybridFVMKernels.cpp
   )

set( dependencyList geosx_core gtest )

if ( ENABLE_MPI )
  set ( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()


if( ENABLE_PVTPackage )
    list( APPEND gtest_geosx_tests
          testCompMultiphaseFlow.cpp
          testCompMultiphaseFlowHybrid.cpp
          testCompMultiphaseFlowUpwind.cpp
          testCompMultiphaseFlowHelpers.cpp )

    set( dependencyList ${dependencyList} PVTPackage )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${gtest_geosx_tests})
  get_filename_component( test_name ${test} NAME_WE )

  blt_add_executable( NAME ${test_name}
                      SOURCES ${test}
                      OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                      DEPENDS_ON ${dependencyList} )

  blt_add_test( NAME ${test_name}
                COMMAND ${test_name} )
endforeach()

# For some reason, BLT is not setting CUDA language for these source files
if ( ENABLE_CUDA )
  set_source_files_properties( ${gtest_geosx_tests} PROPERTIES LANGUAGE CUDA )
endif()
