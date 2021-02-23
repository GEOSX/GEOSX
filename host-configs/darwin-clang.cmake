site_name(HOST_NAME)
set(CONFIG_NAME "${HOST_NAME}-darwin-x86_64-clang@apple-mp" CACHE PATH "") 
message("CONFIG_NAME = ${CONFIG_NAME}")

set(CMAKE_C_COMPILER "/usr/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/clang++" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/local/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/local/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "/usr/local/bin/mpirun" CACHE PATH "")

set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE)

set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)

set(CUDA_ENABLED "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

option(ENABLE_CALIPER "Enables CALIPER" OFF)

set(BLAS_LIBRARIES /usr/lib/libblas.dylib CACHE PATH "" FORCE)
set(LAPACK_LIBRARIES /usr/lib/liblapack.dylib CACHE PATH "" FORCE )

set(ENABLE_DOXYGEN OFF CACHE BOOL "" FORCE)

#set( DOXYGEN_EXECUTABLE /usr/local/bin/doxygen CACHE PATH "" FORCE )
#set( SPHINX_EXECUTABLE /usr/local/bin/sphinx-build CACHE PATH "" FORCE )

set(ENABLE_GEOSX_PTP ON CACHE BOOL "" FORCE)

set(GEOSX_TPL_DIR "/Users/cusini1/geosx/thirdPartyLibs/install-darwin-clang-release" CACHE PATH "" FORCE)
include(${CMAKE_CURRENT_LIST_DIR}/tpls.cmake)
