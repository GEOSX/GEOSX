#
# Specify list of tests
#

set(testSources
    testFiniteElementBase.cpp
    testH1_QuadrilateralFace_Lagrange1_GaussLegendre2.cpp
    testH1_Hexahedron_Lagrange1_GaussLegendre2.cpp
    testH1_Tetrahedron_Lagrange1_Gauss1.cpp
    testH1_Wedge_Lagrange1_Gauss6.cpp
    testH1_Pyramid_Lagrange1_Gauss5.cpp
    testH1_TriangleFace_Lagrange1_Gauss1.cpp
   )

set( dependencyList gtest finiteElement)

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

#
# Add gtest C++ based tests
#
foreach(test ${testSources})
    get_filename_component( test_name ${test} NAME_WE )
    blt_add_executable( NAME ${test_name}
                        SOURCES ${test}
                        OUTPUT_DIR ${TEST_OUTPUT_DIRECTORY}
                        DEPENDS_ON ${dependencyList} )

    blt_add_test( NAME ${test_name} 
                  COMMAND ${test_name}
                  )

    if( ( NOT ${ENABLE_CUDA} ) AND ( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" ) )
        set_target_properties(${test_name} PROPERTIES COMPILE_FLAGS "-Wno-everything")
    endif()

endforeach()

if ( ENABLE_CUDA )
#blt_add_target_compile_flags( TO testShapeFunctionKernel
#                              FLAGS "--maxrregcount=64") 
endif()
