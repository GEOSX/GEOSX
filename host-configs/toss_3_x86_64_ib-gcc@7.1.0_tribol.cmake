

set(CONFIG_NAME "tribol-toss_3_x86_64_ib-gcc@7.1.0" CACHE PATH "") 

#set(TPL_DIR "/usr/gapps/GEOS/geosx/2018_05_17_08_46_08" CACHE PATH "" )
#include("${TPL_DIR}/${CONFIG_NAME}.cmake")

#set(ATK_DIR "/usr/gapps/GEOS/geosx/axom/toss_3_x86_64_ib-gcc@7.1.0-release" CACHE PATH "")

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-7.1.0/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-7.1.0/bin/g++" CACHE PATH "")

set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(ENABLE_MPI ON CACHE BOOL "" FORCE)

set(MPI_HOME             "/usr/tce/packages/mvapich2/mvapich2-2.2-gcc-7.1.0" CACHE PATH "")
set(MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "")
set(MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "")
set(MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "")

set(MPIEXEC              "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")


#set( GEOSX_TPL_ROOT_DIR "/usr/gapps/GEOS/geosx/thirdPartyLibs/" CACHE PATH "" )
set( GEOSX_TPL_ROOT_DIR "../../thirdPartyLibs/" CACHE PATH "" )
get_filename_component(ABS_TPL_ROOT_DIR ${GEOSX_TPL_ROOT_DIR} ABSOLUTE)
get_filename_component( INSTALL_DIR_NAME "${CMAKE_INSTALL_PREFIX}" NAME)
set(GEOSX_TPL_DIR ${ABS_TPL_ROOT_DIR}/${INSTALL_DIR_NAME} CACHE PATH "")
set(SPHINX_EXECUTABLE "/usr/bin/sphinx-build" CACHE PATH "" FORCE)
set(DOXYGEN_EXECUTABLE "/usr/bin/doxygen" CACHE PATH "" FORCE )

set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )


#######################################
# RAJA/CHAI SETUP
#######################################
#set(RAJA_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/raja/" CACHE PATH "" FORCE )
#set(CHAI_DIR "/usr/gapps/GEOS/geosx/cab/gcc-4.9.3/chai/" CACHE PATH "" FORCE )

option(ENABLE_TRIBOL "Enables TRIBOL" ON)
option( RAJA_ENABLE_TBB "" OFF)
option(ENABLE_CHAI "Enables CHAI" OFF)

option( ENABLE_CALIPER "Enables CALIPER" OFF )

set(CUDA_ENABLED      "OFF"       CACHE PATH "" FORCE)
set(CHAI_BUILD_TYPE   "cpu-no-rm" CACHE PATH "" FORCE)
set(CHAI_ARGS         ""          CACHE PATH "" FORCE)

set(ENABLE_OPENMP     "OFF"        CACHE PATH "" FORCE)
set(RAJA_ENABLE_OPENMP "OFF"        CACHE PATH "" FORCE)
 
set(NUM_PROC "36" CACHE PATH "")
