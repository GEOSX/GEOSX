include(${CMAKE_CURRENT_LIST_DIR}/../../host-configs/LLNL/lassen-clang@upstream.cmake)

set(CONFIG_NAME "lassen-clang@upstream-NoMPI" CACHE PATH "" FORCE)
set(GEOSX_TPL_ROOT_DIR /usr/gapps/GEOSX/thirdPartyLibs CACHE PATH "" FORCE)
set(GEOSX_TPL_DIR ${GEOSX_TPL_ROOT_DIR}/2020-05-01/install-${CONFIG_NAME}-release CACHE PATH "" FORCE)

set(DOXYGEN_EXECUTABLE /usr/bin/doxygen CACHE PATH "")

set(DOXYGEN_EXECUTABLE /usr/bin/doxygen CACHE PATH "")

set(ENABLE_MPI OFF CACHE BOOL "" FORCE)
unset(MPI_ROOT CACHE )
unset(MPI_C_COMPILER     CACHE )
unset(MPI_CXX_COMPILER    CACHE )
unset(MPI_Fortran_COMPILER  CACHE )
unset(MPIEXEC                 CACHE )
unset(MPIEXEC_NUMPROC_FLAG   CACHE )
set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC OFF CACHE BOOL "")

set(ENABLE_PAMELA OFF CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_GEOSX_PTP OFF CACHE BOOL "" FORCE)
set(ENABLE_PARMETIS OFF CACHE BOOL "" FORCE)
set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE)
set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER} CACHE STRING "" FORCE)
