###############################################################################
#
# Setup Vista
# This file defines:
#  VISTA_FOUND - If Vista was found
#  VISTA_INCLUDE_DIRS - The Vista include directories
#  VISTA_LIBRARY - The Vista library

# first Check for VISTA_DIR

if(NOT VISTA_DIR)
    MESSAGE(FATAL_ERROR "Could not find Vista. Vista support needs explicit VISTA_DIR")
endif()

#find includes
find_path( VISTA_INCLUDE_DIRS vista/Vista.h
           PATHS  ${VISTA_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( VISTA_LIBRARY NAMES vista libvista
              PATHS ${VISTA_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set VISTA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(VISTA  DEFAULT_MSG
                                  VISTA_INCLUDE_DIRS
                                  VISTA_LIBRARY )
