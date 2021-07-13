####################################
# 3rd Party Dependencies
####################################

macro(find_and_register)
    set(singleValueArgs NAME HEADER)
    set(multiValueArgs INCLUDE_DIRECTORIES 
                       LIBRARY_DIRECTORIES
                       LIBRARIES 
                       EXTRA_LIBRARIES 
                       DEPENDS )

    ## parse the arguments
    cmake_parse_arguments(arg
                          "${options}" "${singleValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT DEFINED arg_NAME)
        message(FATAL_ERROR "The find_and_register required parameter NAME specifies the name of the library to register.")
    endif()

    if(NOT DEFINED arg_INCLUDE_DIRECTORIES)
        message(FATAL_ERROR "The find_and_register required parameter INCLUDE_DIRECTORIES specifies the directories to search for the given header.")
    endif()

    if(NOT DEFINED arg_LIBRARY_DIRECTORIES)
        message(FATAL_ERROR "The find_and_register required parameter LIBRARY_DIRECTORIES specifies the directories to search for the given libraries.")
    endif()

    if(NOT DEFINED arg_HEADER)
        message(FATAL_ERROR "The find_and_register required parameter HEADER specifies the header to search for.")
    endif()

    if(NOT DEFINED arg_LIBRARIES)
        message(FATAL_ERROR "The find_and_register required parameter LIBRARIES specifies the libraries to search for.")
    endif()

    find_path(${arg_NAME}_INCLUDE_DIR ${arg_HEADER}
              PATHS ${arg_INCLUDE_DIRECTORIES}
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    if(${arg_NAME}_INCLUDE_DIR STREQUAL ${arg_NAME}_INCLUDE_DIR-NOTFOUND)
        message(FATAL_ERROR "Could not find '${arg_HEADER}' in '${arg_INCLUDE_DIRECTORIES}'")
    endif()

    blt_find_libraries(FOUND_LIBS ${arg_NAME}_LIBRARIES
                       NAMES ${arg_LIBRARIES}
                       PATHS ${arg_LIBRARY_DIRECTORIES}
                       REQUIRED ON)

    blt_import_library(NAME ${arg_NAME}
                         INCLUDES ${${arg_NAME}_INCLUDE_DIR}
                         LIBRARIES ${${arg_NAME}_LIBRARIES} ${arg_EXTRA_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON
                         DEPENDS_ON ${arg_DEPENDS})

endmacro(find_and_register)


set(thirdPartyLibs "")

################################
# BLAS/LAPACK
################################

include(cmake/thirdparty/FindMathLibraries.cmake)

blt_import_library(NAME blas
                     TREAT_INCLUDES_AS_SYSTEM ON
                     LIBRARIES ${BLAS_LIBRARIES})

blt_import_library(NAME lapack
                     DEPENDS_ON blas
                     TREAT_INCLUDES_AS_SYSTEM ON
                     LIBRARIES ${LAPACK_LIBRARIES})

################################
# Intel MKL
################################
if(ENABLE_MKL)
    message(STATUS "Using Intel MKL")

    blt_import_library(NAME mkl
                         INCLUDES ${MKL_INCLUDE_DIRS}
                         LIBRARIES ${MKL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS mkl)
    set(thirdPartyLibs ${thirdPartyLibs} mkl)

################################
# IBM ESSL
################################
elseif(ENABLE_ESSL)
    message(STATUS "Using up IBM ESSL")

    blt_import_library(NAME essl
                         INCLUDES ${ESSL_INCLUDE_DIRS}
                         LIBRARIES ${ESSL_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(TRILINOS_DEPENDS essl)
    set(thirdPartyLibs ${thirdPartyLibs} essl)
else()
    set(TRILINOS_DEPENDS blas lapack)
    set(thirdPartyLibs ${thirdPartyLibs} blas lapack)
endif()

################################
# HDF5
################################
if(DEFINED HDF5_DIR)
    message(STATUS "HDF5_DIR = ${HDF5_DIR}")

    set(HDF5_ROOT ${HDF5_DIR})
    set(HDF5_USE_STATIC_LIBRARIES FALSE)
    set(HDF5_NO_FIND_PACKAGE_CONFIG_FILE ON)
    include(FindHDF5)

    blt_import_library(NAME hdf5
                         INCLUDES ${HDF5_INCLUDE_DIRS}
                         LIBRARIES ${HDF5_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    set(ENABLE_HDF5 ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} hdf5)
else()
    message(FATAL_ERROR "GEOSX requires hdf5, set HDF5_DIR to the hdf5 installation directory.")
endif()

################################
# Conduit
################################
if(DEFINED CONDUIT_DIR)
    message(STATUS "CONDUIT_DIR = ${CONDUIT_DIR}")

    find_package(Conduit REQUIRED
                 PATHS ${CONDUIT_DIR}/lib/cmake
                 NO_DEFAULT_PATH)

    set( CONDUIT_TARGETS conduit conduit_relay conduit_blueprint )
    foreach( targetName ${CONDUIT_TARGETS} )
        get_target_property( includeDirs 
                             ${targetName}
                             INTERFACE_INCLUDE_DIRECTORIES)
                             
        set_property(TARGET ${targetName} 
                     APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                     ${includeDirs})
    endforeach()
    set(thirdPartyLibs ${thirdPartyLibs} conduit::conduit )
else()
    message(FATAL_ERROR "GEOSX requires conduit, set CONDUIT_DIR to the conduit installation directory.")
endif()

################################
# SILO
################################
if(DEFINED SILO_DIR)
    message(STATUS "SILO_DIR = ${SILO_DIR}")

    find_and_register(NAME silo
                      INCLUDE_DIRECTORIES ${SILO_DIR}/include
                      LIBRARY_DIRECTORIES ${SILO_DIR}/lib
                      HEADER silo.h
                      LIBRARIES siloh5
                      DEPENDS hdf5)

    set(ENABLE_SILO ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} silo)
else()
    message(FATAL_ERROR "GEOSX requires Silo, set SILO_DIR to the Silo installation directory.")
endif()

################################
# PUGIXML
################################
if(DEFINED PUGIXML_DIR)
    message(STATUS "PUGIXML_DIR = ${PUGIXML_DIR}")

    find_package(pugixml REQUIRED
                 PATHS ${PUGIXML_DIR}
                 NO_DEFAULT_PATH)

    set(ENABLE_PUGIXML ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} pugixml)
else()
    message(FATAL_ERROR "GEOSX requires pugixml, set PUGIXML_DIR to the pugixml installation directory.")
endif()

################################
# RAJA
################################
if(DEFINED RAJA_DIR)
    message(STATUS "RAJA_DIR = ${RAJA_DIR}")

    find_package(RAJA REQUIRED
                 PATHS ${RAJA_DIR}
                 NO_DEFAULT_PATH)

    get_target_property(RAJA_INCLUDE_DIRS RAJA INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(RAJA
                          PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${RAJA_INCLUDE_DIRS}")

    set(ENABLE_RAJA ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} RAJA )
else()
    message(FATAL_ERROR "GEOSX requires RAJA, set RAJA_DIR to the RAJA installation directory.")
endif()

################################
# Umpire
################################
if(DEFINED UMPIRE_DIR)
    message(STATUS "UMPIRE_DIR = ${UMPIRE_DIR}")

    find_package(umpire REQUIRED
                 PATHS ${UMPIRE_DIR}
                 NO_DEFAULT_PATH)

    set(ENABLE_UMPIRE ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} umpire)
else()
    message(FATAL_ERROR "GEOSX requires Umpire, set UMPIRE_DIR to the Umpire installation directory.")
endif()

################################
# CHAI
################################
if(DEFINED CHAI_DIR)
    message(STATUS "CHAI_DIR = ${CHAI_DIR}")

    find_package(chai REQUIRED
                 PATHS ${CHAI_DIR}
                 NO_DEFAULT_PATH)

    get_target_property(CHAI_INCLUDE_DIRS chai INTERFACE_INCLUDE_DIRECTORIES)
    set_target_properties(chai
                          PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${CHAI_INCLUDE_DIRS}")

    set(ENABLE_CHAI ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} chai)
else()
    message(FATAL_ERROR "GEOSX requires CHAI, set CHAI_DIR to the CHAI installation directory.")
endif()

################################
# Adiak
################################
if(DEFINED ADIAK_DIR)
    message(STATUS "ADIAK_DIR = ${ADIAK_DIR}")

    find_package(adiak REQUIRED
                 PATHS ${ADIAK_DIR}
                 NO_DEFAULT_PATH)

    set_property(TARGET adiak
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 ${adiak_INCLUDE_DIR} )
    set_property(TARGET adiak
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 ${adiak_INCLUDE_DIR} )

    set(ENABLE_ADIAK ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} adiak)
else()
    if(ENABLE_ADIAK)
        message(WARNING "ENABLE_ADIAK is ON but ADIAK_DIR isn't defined.")
    endif()

    set(ENABLE_ADIAK OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Adiak.")
endif()

################################
# Caliper
################################
if(DEFINED CALIPER_DIR)
    message(STATUS "CALIPER_DIR = ${CALIPER_DIR}")

    find_package(caliper REQUIRED
                 PATHS ${CALIPER_DIR}
                 NO_DEFAULT_PATH)

    set_property(TARGET caliper
                 APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                 ${caliper_INCLUDE_PATH} )

    set_property(TARGET caliper
                 APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                 ${caliper_INCLUDE_PATH} )

    set(ENABLE_CALIPER ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} caliper)
else()
    if(ENABLE_CALIPER)
        message(WARNING "ENABLE_CALIPER is ON but CALIPER_DIR isn't defined.")
    endif()

    set(ENABLE_CALIPER OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Caliper.")
endif()

################################
# MATHPRESSO
################################
if(DEFINED MATHPRESSO_DIR)
    message(STATUS "MATHPRESSO_DIR = ${MATHPRESSO_DIR}")

    find_and_register(NAME mathpresso
                      INCLUDE_DIRECTORIES ${MATHPRESSO_DIR}/include
                      LIBRARY_DIRECTORIES ${MATHPRESSO_DIR}/lib
                      HEADER mathpresso/mathpresso.h
                      LIBRARIES mathpresso)

    set(ENABLE_MATHPRESSO ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} mathpresso)
else()
    if(ENABLE_MATHPRESSO)
        message(WARNING "ENABLE_MATHPRESSO is ON but MATHPRESSO_DIR isn't defined.")
    endif()

    set(ENABLE_MATHPRESSO OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using mathpresso.")
endif()

################################
# METIS
################################
if(DEFINED METIS_DIR)
    message(STATUS "METIS_DIR = ${METIS_DIR}")

    find_and_register(NAME metis
                      INCLUDE_DIRECTORIES ${METIS_DIR}/include
                      LIBRARY_DIRECTORIES ${METIS_DIR}/lib
                      HEADER metis.h
                      LIBRARIES metis)

    set(ENABLE_METIS ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} metis)
else()
    if(ENABLE_METIS)
        message(WARNING "ENABLE_METIS is ON but METIS_DIR isn't defined.")
    endif()

    set(ENABLE_METIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using METIS.")
endif()

################################
# PARMETIS
################################
if(DEFINED PARMETIS_DIR)
    message(STATUS "PARMETIS_DIR = ${PARMETIS_DIR}")

    find_and_register(NAME parmetis
                      INCLUDE_DIRECTORIES ${PARMETIS_DIR}/include
                      LIBRARY_DIRECTORIES ${PARMETIS_DIR}/lib
                      HEADER parmetis.h
                      LIBRARIES parmetis
                      DEPENDS metis)

    set(ENABLE_PARMETIS ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} parmetis)
else()
    if(ENABLE_PARMETIS)
        message(WARNING "ENABLE_PARMETIS is ON but PARMETIS_DIR isn't defined.")
    endif()

    set(ENABLE_PARMETIS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using ParMETIS.")
endif()

################################
# SUPERLU_DIST
################################
if(DEFINED SUPERLU_DIST_DIR)
    message(STATUS "SUPERLU_DIST_DIR = ${SUPERLU_DIST_DIR}")

    find_and_register(NAME superlu_dist
                      INCLUDE_DIRECTORIES ${SUPERLU_DIST_DIR}/include
                      LIBRARY_DIRECTORIES ${SUPERLU_DIST_DIR}/lib PATHS ${SUPERLU_DIST_DIR}/lib64
                      HEADER superlu_defs.h
                      LIBRARIES superlu_dist
                      DEPENDS parmetis blas lapack)

    set(ENABLE_SUPERLU_DIST ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} superlu_dist)
else()
    if(ENABLE_SUPERLU_DIST)
        message(WARNING "ENABLE_SUPERLU_DIST is ON but SUPERLU_DIST_DIR isn't defined.")
    endif()

    set(ENABLE_SUPERLU_DIST OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using superlu_dist.")
endif()

################################
# SUITESPARSE
################################
if(DEFINED SUITESPARSE_DIR)
    message(STATUS "SUITESPARSE_DIR = ${SUITESPARSE_DIR}")

    find_and_register(NAME suitesparse
                      INCLUDE_DIRECTORIES ${SUITESPARSE_DIR}/include
                      LIBRARY_DIRECTORIES ${SUITESPARSE_DIR}/lib ${SUITESPARSE_DIR}/lib64
                      HEADER umfpack.h
                      LIBRARIES umfpack
                      DEPENDS blas lapack)

    set(ENABLE_SUITESPARSE ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} suitesparse)
else()
    if(ENABLE_SUITESPARSE)
        message(WARNING "ENABLE_SUITESPARSE is ON but SUITESPARSE_DIR isn't defined.")
    endif()

    set(ENABLE_SUITESPARSE OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using SuiteSparse.")
endif()

################################
# HYPRE
################################
if(DEFINED HYPRE_DIR AND ENABLE_HYPRE)
    message(STATUS "HYPRE_DIR = ${HYPRE_DIR}")

    if( ENABLE_HYPRE_CUDA )
        set( EXTRA_LIBS ${CUDA_cusparse_LIBRARY} ${CUDA_curand_LIBRARY} )
    endif()

    find_and_register(NAME hypre
                      INCLUDE_DIRECTORIES ${HYPRE_DIR}/include
                      LIBRARY_DIRECTORIES ${HYPRE_DIR}/lib 
                      HEADER HYPRE.h
                      LIBRARIES HYPRE
                      EXTRA_LIBRARIES ${EXTRA_LIBS}
                      DEPENDS blas lapack superlu_dist)

    # if( ENABLE_CUDA AND ( NOT ENABLE_HYPRE_CUDA ) )
    #   set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
    #   if( GEOSX_LA_INTERFACE STREQUAL "Hypre")
    #     message( FATAL_ERROR "Hypre LAI selected, but ENABLE_HYPRE_CUDA not ON while ENABLE_CUDA is ON.")
    #   endif()
    # else()
    #   set(ENABLE_HYPRE ON CACHE BOOL "")
    # endif()
    
    set(ENABLE_HYPRE ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} hypre)
else()
    if(ENABLE_HYPRE)
        message(WARNING "ENABLE_HYPRE is ON but HYPRE_DIR isn't defined.")
    endif()

    set(ENABLE_HYPRE OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using HYPRE.")
endif()

################################
# TRILINOS
################################
if(DEFINED TRILINOS_DIR AND ENABLE_TRILINOS)
    message(STATUS "TRILINOS_DIR = ${TRILINOS_DIR}")
  
    include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)
  
    list(REMOVE_ITEM Trilinos_LIBRARIES "gtest")
    list(REMOVE_DUPLICATES Trilinos_LIBRARIES)

    blt_import_library(NAME trilinos
                         DEPENDS_ON ${TRILINOS_DEPENDS}
                         INCLUDES ${Trilinos_INCLUDE_DIRS} 
                         LIBRARIES ${Trilinos_LIBRARIES}
                         TREAT_INCLUDES_AS_SYSTEM ON)

    # This conditional is due to the lack of mixedInt support on hypre GPU.
    # This can be removed when support is added into hypre.
    if( NOT ENABLE_HYPRE_CUDA )
        set(ENABLE_TRILINOS ON CACHE BOOL "")
    endif()
    set(thirdPartyLibs ${thirdPartyLibs} trilinos)
else()
    if(ENABLE_TRILINOS)
        message(WARNING "ENABLE_TRILINOS is ON but TRILINOS_DIR isn't defined.")
    endif()

    set(ENABLE_TRILINOS OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using Trilinos")
endif()

###############################
# PETSC
###############################
if(DEFINED PETSC_DIR AND ENABLE_PETSC)
    message(STATUS "PETSC_DIR = ${PETSC_DIR}")

    find_and_register(NAME petsc
                      INCLUDE_DIRECTORIES ${PETSC_DIR}/include
                      LIBRARY_DIRECTORIES ${PETSC_DIR}/lib
                      HEADER petscvec.h
                      LIBRARIES petsc
                      DEPENDS metis superlu_dist blas lapack)

    set(ENABLE_PETSC ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} petsc)
else()
    if(ENABLE_PETSC)
        message(WARNING "ENABLE_PETSC is ON but PETSC_DIR isn't defined.")
    endif()

    set(ENABLE_PETSC OFF CACHE BOOL "" FORCE)
    message(STATUS "Not using PETSc")
endif()

################################
# VTK
################################
if(DEFINED VTK_DIR)
    message(STATUS "VTK_DIR = ${VTK_DIR}")
    find_package(VTK REQUIRED
                 PATHS ${VTK_DIR}
                 NO_DEFAULT_PATH)

    set( VTK_TARGETS VTK::WrappingTools VTK::WrapHierarchy VTK::WrapPython VTK::WrapPythonInit VTK::ParseJava VTK::WrapJava VTK::loguru VTK::kwiml VTK::vtksys VTK::utf8 VTK::CommonCore VTK::CommonMath VTK::CommonTransforms VTK::CommonMisc VTK::CommonSystem VTK::CommonDataModel VTK::CommonExecutionModel VTK::doubleconversion VTK::lz4 VTK::lzma VTK::zlib VTK::IOCore VTK::expat VTK::IOXMLParser VTK::IOXML )
    foreach( targetName ${VTK_TARGETS} )

        get_target_property( includeDirs ${targetName}  INTERFACE_INCLUDE_DIRECTORIES)
    
        set_property(TARGET ${targetName} 
                     APPEND PROPERTY INTERFACE_SYSTEM_INCLUDE_DIRECTORIES
                     ${includeDirs})
    endforeach()
    
    set(ENABLE_VTK ON CACHE BOOL "")
    set(thirdPartyLibs ${thirdPartyLibs} vtk)
else()
    if(ENABLE_VTK)
        message(WARNING "ENABLE_VTK is ON but VTK_DIR isn't defined.")
    endif()

    set(ENABLE_VTK OFF CACHE BOOL "")
    message(STATUS "Not using VTK")
endif()

################################
# uncrustify
################################
if(UNCRUSTIFY_FOUND)
    message(STATUS "UNCRUSTIFY_EXECUTABLE = ${UNCRUSTIFY_EXECUTABLE}")

    # targets for verifying formatting
    if( NOT TARGET uncrustify_check )
        add_custom_target(uncrustify_check)
        add_dependencies(check uncrustify_check)
    endif()

    # targets for modifying formatting
    if( NOT TARGET uncrustify_style )
        add_custom_target(uncrustify_style)
        add_dependencies(style uncrustify_style)
    endif()
else()
    message(STATUS "Not using uncrustify.")
endif()

################################
# doxygen
################################
if(DOXYGEN_FOUND)
    message(STATUS "DOXYGEN_EXECUTABLE = ${DOXYGEN_EXECUTABLE}")
else()
    message(STATUS "Not using doxygen.")
endif()

################################
# Python
################################
if(ENABLE_PYGEOSX)
    message(STATUS "Python3_EXECUTABLE=${Python3_EXECUTABLE}")
    find_package(Python3 REQUIRED
                 COMPONENTS Development NumPy)

    message(STATUS "Python3_INCLUDE_DIRS = ${Python3_INCLUDE_DIRS}")
    message(STATUS "Python3_LIBRARY_DIRS = ${Python3_LIBRARY_DIRS}")
    message(STATUS "Python3_NumPy_INCLUDE_DIRS = ${Python3_NumPy_INCLUDE_DIRS}")

    if(DEFINED ENABLE_PYLVARRAY AND NOT ENABLE_PYLVARRAY)
        message(FATAL_ERROR "Cannot build pygeosx without pylvarray")
    else()
        set(ENABLE_PYLVARRAY ON CACHE BOOL "")
    endif()

    set(thirdPartyLibs ${thirdPartyLibs} Python3::Python Python3::NumPy)
else()
    message(STATUS "Not building pygeosx.")
endif()

################################
# LAI
################################
string(TOUPPER "${GEOSX_LA_INTERFACE}" upper_LAI)
if(NOT ENABLE_${upper_LAI})
  message(FATAL_ERROR "${GEOSX_LA_INTERFACE} LA interface is selected, but ENABLE_${upper_LAI} is OFF")
endif()
option(GEOSX_LA_INTERFACE_${upper_LAI} "${upper_LAI} LA interface is selected" ON)

message(STATUS "thirdPartyLibs = ${thirdPartyLibs}")
