message(STATUS "\nProcessing SetupGeosxThirdParty.cmake")


####################################
# 3rd Party Dependencies
####################################

set( thirdPartyLibs "")


################################
# Conduit
################################
if (NOT EXISTS ${CONDUIT_DIR})
    message(STATUS "Using conduit from thirdPartyLibs")
    set(CONDUIT_DIR ${GEOSX_TPL_DIR}/conduit CACHE PATH "")
endif()

if (EXISTS ${CONDUIT_DIR})
  message(STATUS "CONDUIT_DIR = ${CONDUIT_DIR}" )
  include( cmake/thirdparty/FindConduit.cmake )
  blt_register_library( NAME conduit
                        INCLUDES ${CONDUIT_INCLUDE_DIRS} 
                        LIBRARIES  conduit
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME conduit_blueprint
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES conduit_blueprint
                        TREAT_INCLUDES_AS_SYSTEM ON )

  blt_register_library( NAME conduit_relay
                        INCLUDES ${CONDUIT_INCLUDE_DIRS}
                        LIBRARIES conduit_relay
                        TREAT_INCLUDES_AS_SYSTEM ON )
                        
  set( CONDUIT_FOUND ON CACHE BOOL "" )
  set( thirdPartyLibs ${thirdPartyLibs} conduit conduit_blueprint conduit_relay )
else()
  set( CONDUIT_FOUND OFF CACHE BOOL "" )
  message(STATUS "Not using conduit" )
endif()

################################
# AXOM
################################
include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindATK.cmake)


################################
# HDF5
################################
if( EXISTS ${HDF5_DIR})
    message(STATUS "Using system HDF5 found at ${HDF5_DIR}")
else()
    message(STATUS "Using HDF5 from thirdPartyLibs")
    set(HDF5_DIR ${GEOSX_TPL_DIR}/hdf5)
endif()

if (HDF5_DIR)
  include(cmake/thirdparty/FindHDF5.cmake)
  blt_register_library(NAME hdf5
                       INCLUDES ${HDF5_INCLUDE_DIRS}
                       LIBRARIES ${HDF5_LIBRARIES} dl
                       TREAT_INCLUDES_AS_SYSTEM ON )
                       
  set( thirdPartyLibs ${thirdPartyLibs} hdf5 )

endif()

################################
# SILO
################################
if( EXISTS ${SILO_DIR})
    message(STATUS "Using system SILO found at ${SILO_DIR}")
else()
    message(STATUS "Using SILO from thirdPartyLibs")
    set(SILO_DIR ${GEOSX_TPL_DIR}/silo)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindSilo.cmake)
if (NOT SILO_FOUND)
    message(FATAL_ERROR "SILO not found in ${SILO_DIR}. Maybe you need to build it")
endif()
blt_register_library( NAME silo
                      INCLUDES ${SILO_INCLUDE_DIRS}
                      LIBRARIES ${SILO_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON
                      DEPENDS_ON hdf5 )

set( thirdPartyLibs ${thirdPartyLibs} silo )

################################
# RAJA
################################
if( EXISTS ${RAJA_DIR})
    message(STATUS "Using system RAJA found at ${RAJA_DIR}")
else()
    message(STATUS "Using RAJA from thirdPartyLibs")
    set(RAJA_DIR ${GEOSX_TPL_DIR}/raja)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindRAJA.cmake)
if (NOT RAJA_FOUND)
    message(FATAL_ERROR "RAJA not found in ${RAJA_DIR}. Maybe you need to build it")
endif()    
blt_register_library( NAME raja
                      INCLUDES ${RAJA_INCLUDE_DIRS}
                      LIBRARIES ${RAJA_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )

set( thirdPartyLibs ${thirdPartyLibs} raja )  

################################
# CHAI
################################
include(${CMAKE_SOURCE_DIR}/cmake/thirdparty/FindCHAI.cmake)

################################
# FPARSER
################################
if( USE_FPARSER )

message(STATUS "Using FPARSER from thirdPartyLibs")
set(FPARSER_INSTALL_DIR ${GEOSX_TPL_DIR}/fparser)

find_path( FPARSER_INCLUDE_DIRS fparser.h
           PATHS  ${FPARSER_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( FPARSER_LIBRARY NAMES fparser
              PATHS ${FPARSER_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FPARSER  DEFAULT_MSG
                                  FPARSER_INCLUDE_DIRS
                                  FPARSER_LIBRARY )
if (NOT FPARSER_FOUND)
    message(STATUS "FPARSER not found in ${FPARSER_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME fparser
                      INCLUDES ${FPARSER_INCUDE_DIRS} 
                      LIBRARIES ${FPARSER_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )

set( thirdPartyLibs ${thirdPartyLibs} fparser )  

endif()

################################
# CALIPER
################################
if( ENABLE_CALIPER )
    message( STATUS "setting up caliper" )

    if( EXISTS ${CALIPER_DIR} )
        message( STATUS "Found system caliper" )
        message( STATUS "Using system CALIPER found at ${CALIPER_DIR}")
        set(CALIPER_FOUND TRUE)    
    else()
        message(STATUS "Using CALIPER from thirdPartyLibs")
        set(CALIPER_DIR ${GEOSX_TPL_DIR}/caliper)
        set(CALIPER_FOUND TRUE)
    endif()

    find_path( CALIPER_INCLUDE_DIRS caliper/Caliper.h
            PATHS  ${CALIPER_DIR}/include
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)

    set( caliper_lib_list caliper-mpi caliper )

    message(STATUS "looking for libs in ${CALIPER_DIR}")
    blt_find_libraries( FOUND_LIBS CALIPER_LIBRARIES
                        NAMES ${caliper_lib_list}
                        PATHS ${CALIPER_DIR}/lib ${CALIPER_DIR}/lib64
                    )

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(CALIPER  DEFAULT_MSG
                                    CALIPER_INCLUDE_DIRS
                                    CALIPER_LIBRARIES )

    if (NOT CALIPER_FOUND)
        message(FATAL_ERROR "CALIPER not found in ${CALIPER_DIR}. Maybe you need to build it")
    else() 
        message(STATUS "CALIPER_INCLUDE_DIRS = ${CALIPER_INCLUDE_DIRS}")
        message(STATUS "CALIPER_LIBRARIES = ${CALIPER_LIBRARIES}")    
    endif()
    blt_register_library( NAME caliper
                        INCLUDES ${CALIPER_INCLUDE_DIRS}
                        LIBRARIES ${CALIPER_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} caliper )  
endif()

################################
# ASMJIT / MATHPRESSO
################################
if (NOT EXISTS ${MATHPRESSO_DIR})
    message(STATUS "Using mathpresso from thirdPartyLibs")
    set(MATHPRESSO_DIR ${GEOSX_TPL_DIR}/mathpresso CACHE PATH "")
endif()

if (EXISTS ${MATHPRESSO_DIR})
    message( STATUS "setting up asmjit" )
    set(ASMJIT_LOCAL_DIR ${PROJECT_BINARY_DIR}/thirdparty/asmjit/src/asmjit)
    set(ASMJIT_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/thirdparty/asmjit)

    message( STATUS "setting up MathPresso" )
    find_path( MATHPRESSO_INCLUDE_DIRS mathpresso/mathpresso.h
            PATHS  ${MATHPRESSO_DIR}/include
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)

    find_library( MATHPRESSO_LIBRARY NAMES mathpresso
                PATHS ${MATHPRESSO_DIR}/lib
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)


    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(MATHPRESSO  DEFAULT_MSG
                                    MATHPRESSO_INCLUDE_DIRS
                                    MATHPRESSO_LIBRARY )
                                    
    if (NOT MATHPRESSO_FOUND)
        message(FATAL_ERROR "MATHPRESSO not found in ${MATHPRESSO_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME mathpresso
                        INCLUDES ${MATHPRESSO_INCLUDE_DIRS}
                        LIBRARIES ${MATHPRESSO_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )
    
    set(MATHPRESSO_FOUND ON CACHE BOOL "")
    set( thirdPartyLibs ${thirdPartyLibs} mathpresso )
else()
    set(MATHPRESSO_FOUND OFF CACHE BOOL "")
    message(STATUS "Not using mathpresso")
endif()

################################
# PUGIXML
################################
message( STATUS "setting up pugixml" )
set(PUGIXML_DIR ${GEOSX_TPL_DIR}/pugixml)


find_path( PUGIXML_INCLUDE_DIRS NAMES pugixml.hpp
           PATHS  ${PUGIXML_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)


if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    find_library( PUGIXML_LIBRARY NAMES pugixml
                  PATHS ${PUGIXML_DIR}/lib
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)

else()
    find_library( PUGIXML_LIBRARY NAMES pugixml
                  PATHS ${PUGIXML_DIR}/lib64 ${PUGIXML_DIR}/lib
                  NO_DEFAULT_PATH
                  NO_CMAKE_ENVIRONMENT_PATH
                  NO_CMAKE_PATH
                  NO_SYSTEM_ENVIRONMENT_PATH
                  NO_CMAKE_SYSTEM_PATH)
endif()

find_package_handle_standard_args(PUGIXML  DEFAULT_MSG
                                  PUGIXML_INCLUDE_DIRS
                                  PUGIXML_LIBRARY )


if (NOT PUGIXML_FOUND)
    message(FATAL_ERROR "PUGIXML not found in ${PUGIXML_DIR}. Maybe you need to build it")
endif()

blt_register_library( NAME pugixml
                      INCLUDES ${PUGIXML_INCLUDE_DIRS}
                      LIBRARIES ${PUGIXML_LIBRARY}
                      TREAT_INCLUDES_AS_SYSTEM ON )

set( thirdPartyLibs ${thirdPartyLibs} pugixml )  



################################
# BLAS/LAPACK
################################

include( cmake/thirdparty/FindMathLibraries.cmake )

blt_register_library( NAME blas
                      TREAT_INCLUDES_AS_SYSTEM ON
                      LIBRARIES ${BLAS_LIBRARIES}
                      )

blt_register_library( NAME lapack
                      DEPENDS_ON blas
                      TREAT_INCLUDES_AS_SYSTEM ON
                      LIBRARIES ${LAPACK_LIBRARIES}
                      )

################################
# Intel MKL
################################
if (ENABLE_MKL)
    message( STATUS "setting up Intel MKL" )

    blt_register_library( NAME mkl
                          INCLUDES ${MKL_INCLUDE_DIRS}
                          LIBRARIES ${MKL_LIBRARIES}
                          TREAT_INCLUDES_AS_SYSTEM ON )
    
    set( TRILINOS_DEPENDS mkl )
    set( thirdPartyLibs ${thirdPartyLibs} mkl )

################################
# IBM ESSL
################################
elseif (ENABLE_ESSL)
    message( STATUS "setting up IBM ESSL" )

    blt_register_library( NAME essl
                          INCLUDES ${ESSL_INCLUDE_DIRS}
                          LIBRARIES ${ESSL_LIBRARIES}
                          TREAT_INCLUDES_AS_SYSTEM ON )
    
    set( TRILINOS_DEPENDS essl )
    set( thirdPartyLibs ${thirdPartyLibs} essl )
else()
    set( TRILINOS_DEPENDS blas lapack )
    set( thirdPartyLibs ${thirdPartyLibs} blas lapack )
endif()

################################
# TRILINOS
################################
if( ENABLE_TRILINOS )

  if(EXISTS ${TRILINOS_DIR})
  
  else()
      message( STATUS "setting up TRILINOS" )
      set(TRILINOS_DIR ${GEOSX_TPL_DIR}/trilinos)
  endif()
  
  include(${TRILINOS_DIR}/lib/cmake/Trilinos/TrilinosConfig.cmake)
  
  list(REMOVE_ITEM Trilinos_LIBRARIES "gtest")
  list(REMOVE_DUPLICATES Trilinos_LIBRARIES)
  message(STATUS "Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
  message(STATUS "Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
  
  blt_register_library( NAME trilinos
                        DEPENDS_ON ${TRILINOS_DEPENDS}
                        INCLUDES ${Trilinos_INCLUDE_DIRS} 
                        LIBRARIES ${Trilinos_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )
  set( thirdPartyLibs ${thirdPartyLibs} trilinos )

endif()

################################
# METIS
################################
if( ENABLE_METIS )
    message( STATUS "setting up METIS" )

    set(METIS_DIR ${GEOSX_TPL_DIR}/metis)

    find_path( METIS_INCLUDE_DIRS metis.h
       PATHS  ${METIS_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( METIS_LIBRARY NAMES metis
          PATHS ${METIS_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message(STATUS "METIS_INCLUDE_DIRS = ${METIS_INCLUDE_DIRS}" )
    message(STATUS "METIS_LIBRARY = ${METIS_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(METIS DEFAULT_MSG METIS_INCLUDE_DIRS METIS_LIBRARY )
    if (NOT METIS_FOUND)
        message(FATAL_ERROR "METIS not found in ${METIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME metis
                            INCLUDES ${METIS_INCLUDE_DIRS} 
                LIBRARIES ${METIS_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} metis )
endif()

################################
# PARMETIS
################################
if( ENABLE_PARMETIS )
    message( STATUS "setting up PARMETIS" )

    set(PARMETIS_DIR ${GEOSX_TPL_DIR}/parmetis)

    find_path( PARMETIS_INCLUDE_DIRS parmetis.h
           PATHS  ${PARMETIS_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( PARMETIS_LIBRARY NAMES parmetis
              PATHS ${PARMETIS_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message(STATUS "PARMETIS_INCLUDE_DIRS = ${PARMETIS_INCLUDE_DIRS}" )
    message(STATUS "PARMETIS_LIBRARY = ${PARMETIS_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(PARMETIS  DEFAULT_MSG
                                    PARMETIS_INCLUDE_DIRS
                                    PARMETIS_LIBRARY )
    if (NOT PARMETIS_FOUND)
        message(FATAL_ERROR "PARMETIS not found in ${PARMETIS_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME parmetis
                          INCLUDES ${PARMETIS_INCLUDE_DIRS} 
                          LIBRARIES ${PARMETIS_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} parmetis )
endif()

################################
# SUPERLU_DIST
################################
if( ENABLE_SUPERLU_DIST)
    message( STATUS "setting up SUPERLU_DIST" )

    set(SUPERLU_DIST_DIR ${GEOSX_TPL_DIR}/superlu_dist)

    find_path( SUPERLU_DIST_INCLUDE_DIRS superlu_defs.h
        PATHS  ${SUPERLU_DIST_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

   find_library( SUPERLU_DIST_LIBRARY NAMES superlu_dist
          PATHS ${SUPERLU_DIST_DIR}/lib PATHS ${SUPERLU_DIST_DIR}/lib64
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)

    message(STATUS "SUPERLU_DIST_INCLUDE_DIRS = ${SUPERLU_DIST_INCLUDE_DIRS}" )
    message(STATUS "SUPERLU_DIST_LIBRARY = ${SUPERLU_DIST_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(SUPERLU_DIST  DEFAULT_MSG
                                SUPERLU_DIST_INCLUDE_DIRS
                    SUPERLU_DIST_LIBRARY )
    if (NOT SUPERLU_DIST_FOUND)
        message(FATAL_ERROR "SUPERLU_DIST not found in ${SUPERLU_DIST_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME superlu_dist
                          DEPENDS_ON lapack blas
                          INCLUDES ${SUPERLU_DIST_INCLUDE_DIRS} 
                          LIBRARIES ${SUPERLU_DIST_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} superlu_dist )
endif()

################################
# HYPRE
################################
if( ENABLE_HYPRE )
    message( STATUS "setting up HYPRE" )

    set(HYPRE_DIR ${GEOSX_TPL_DIR}/hypre)

    find_path( HYPRE_INCLUDE_DIRS HYPRE.h
            PATHS  ${HYPRE_DIR}/include
            NO_DEFAULT_PATH
            NO_CMAKE_ENVIRONMENT_PATH
            NO_CMAKE_PATH
            NO_SYSTEM_ENVIRONMENT_PATH
            NO_CMAKE_SYSTEM_PATH)

    find_library( HYPRE_LIBRARY NAMES HYPRE
                PATHS ${HYPRE_DIR}/lib
                NO_DEFAULT_PATH
                NO_CMAKE_ENVIRONMENT_PATH
                NO_CMAKE_PATH
                NO_SYSTEM_ENVIRONMENT_PATH
                NO_CMAKE_SYSTEM_PATH)

    message(STATUS "HYPRE_INCLUDE_DIRS = ${HYPRE_INCLUDE_DIRS}" )
    message(STATUS "HYPRE_LIBRARY = ${HYPRE_LIBRARY}" )
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(HYPRE  DEFAULT_MSG
                                    HYPRE_INCLUDE_DIRS
                                    HYPRE_LIBRARY )
    if (NOT HYPRE_FOUND)
        message(FATAL_ERROR "HYPRE not found in ${HYPRE_DIR}. Maybe you need to build it")
    endif()

    blt_register_library( NAME hypre
                          DEPENDS_ON superlu_dist blas lapack
                          INCLUDES ${HYPRE_INCLUDE_DIRS}
                          LIBRARIES ${HYPRE_LIBRARY}
                          TREAT_INCLUDES_AS_SYSTEM ON )

    set( thirdPartyLibs ${thirdPartyLibs} hypre )
endif()


################################
# unsrustify
################################
message(STATUS "UNCRUSTIFY_FOUND = ${UNCRUSTIFY_FOUND}")
if(UNCRUSTIFY_FOUND)
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
    
endif()

################################
# TRIBOL
################################
if (ENABLE_TRIBOL)
  set(CHAI_DIR ${GEOSX_TPL_DIR}/chai)
  set(TRIBOL_DIR ${GEOSX_TPL_DIR}/tribol)
  set(VISTA_DIR ${GEOSX_TPL_DIR}/vista)
  set(RXB_DIR ${GEOSX_TPL_DIR}/rxb)
  set(WORLDS_CORE_DIR ${GEOSX_TPL_DIR}/worlds_core)
  set(LLNL_GLOBALID_DIR ${GEOSX_TPL_DIR}/LLNL_GlobalID)

  # ale3d chai version differs from github version
  include (${CHAI_DIR}/lib/cmake/chai-targets.cmake)
  set (CHAI_INCLUDE_DIRS ${CHAI_DIR}/include)
  blt_register_library( NAME chai
                        INCLUDES ${CHAI_INCLUDE_DIRS}
                        LIBRARIES ${CHAI_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

  include(cmake/thirdparty/FindLLNL_GlobalID.cmake)
  blt_register_library( NAME llnl_globalid
                        INCLUDES ${LLNL_GLOBALID_INCLUDE_DIRS} 
                        TREAT_INCLUDES_AS_SYSTEM ON )

  # Register extra axom libraries needed for tribol
  blt_register_library( NAME axom_utils
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  axom_utils
                        TREAT_INCLUDES_AS_SYSTEM ON)
                        
  blt_register_library( NAME mint
                        INCLUDES ${ATK_INCLUDE_DIRS} 
                        LIBRARIES  mint
                        TREAT_INCLUDES_AS_SYSTEM ON)
                        
  set( thirdPartyLibs ${thirdPartyLibs} mint axom_utils)
  include(cmake/thirdparty/FindTribol.cmake)
  blt_register_library( NAME tribol
                        DEPENDS_ON slic mint
                        INCLUDES ${TRIBOL_INCLUDE_DIRS} 
                        LIBRARIES  ${TRIBOL_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

  include(cmake/thirdparty/FindVista.cmake)
  blt_register_library( NAME vista
                        DEPENDS_ON llnl_globalid chai
                        INCLUDES ${VISTA_INCLUDE_DIRS} 
                        LIBRARIES  ${VISTA_LIBRARY}
                        TREAT_INCLUDES_AS_SYSTEM ON )

  include(cmake/thirdparty/FindRXB.cmake)
  blt_register_library( NAME rxb
                        DEPENDS_ON llnl_globalid mpi
                        INCLUDES ${RXB_INCLUDE_DIRS} 
                        LIBRARIES  ${RXB_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )

  include(cmake/thirdparty/FindWorldsCore.cmake)
  blt_register_library( NAME worlds_core
                        DEPENDS_ON llnl_globalid rxb vista mpi
                        INCLUDES ${WORLDS_CORE_INCLUDE_DIRS} 
                        LIBRARIES  ${WORLDS_CORE_LIBRARIES}
                        TREAT_INCLUDES_AS_SYSTEM ON )

else()
  message( "Not using tribol" )
endif()

message(STATUS "Leaving SetupGeosxThirdParty.cmake\n")

