message("\nProcessing geosxOptions.cmake")


if( BLT_CXX_STD STREQUAL c++11)
    MESSAGE(FATAL_ERROR "c++11 is enabled. GEOSX requires c++14")
endif(BLT_CXX_STD STREQUAL c++11)

if(NOT BLT_CXX_STD STREQUAL c++14)
    MESSAGE(FATAL_ERROR "c++14 is NOT enabled. GEOSX requires c++14")
endif(NOT BLT_CXX_STD STREQUAL c++14)

if( (CMAKE_CXX_STANDARD EQUAL 11) OR (CMAKE_CXX_STANDARD EQUAL 14) )
  add_definitions("-DUSE_CXX11")
endif()


message( "ENABLE_CONTAINERARRAY_RETURN_PTR = ${ENABLE_CONTAINERARRAY_RETURN_PTR}" )
if( ENABLE_CONTAINERARRAY_RETURN_PTR )
  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DCONTAINERARRAY_RETURN_PTR=1)  
else()
  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -DCONTAINERARRAY_RETURN_PTR=0)  
endif()


message("CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")




#if( (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")  AND (CMAKE_CXX_STANDARD EQUAL 14) )
#    blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -std=c++14)
#endif()


	blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT "${OpenMP_CXX_FLAGS}")
  blt_append_custom_compiler_flag( FLAGS_VAR CMAKE_CXX_FLAGS
                                   GNU "-Wall -pedantic-errors -Wno-abi -Wextra  -Wshadow -Wfloat-equal	-Wcast-align	-Wpointer-arith	-Wwrite-strings	-Wcast-qual	-Wswitch-default  -Wno-vla  -Wno-switch-default  -Wno-unused-parameter  -Wno-unused-variable  -Wno-unused-function" 
                                   CLANG "-Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded -Wno-missing-prototypes -Wno-covered-switch-default -Wno-double-promotion -Wno-documentation -Wno-switch-enum -Wno-sign-conversion -Wno-unused-parameter -Wno-unused-variable -Wno-reserved-id-macro" 
                                  )

  blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -rdynamic)
    set(CMAKE_EXE_LINKER_FLAGS "-rdynamic")
    
    

message("CMAKE_CXX_FLAGS = ${CMAKE_CXX_FLAGS}")

message("Leaving geosxOptions.cmake\n")