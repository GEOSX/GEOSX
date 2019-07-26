###############################################################################
#
# Setup TRIBOL
# This file defines:
#  TRIBOL_FOUND - If TRIBOL was found
#  TRIBOL_INCLUDE_DIRS - The TRIBOL include directories
#  TRIBOL_LIBRARY - The TRIBOL library

# first Check for TRIBOL_DIR

if(NOT TRIBOL_DIR)
    MESSAGE(FATAL_ERROR "Could not find TRIBOL. TRIBOL support needs explicit TRIBOL_DIR")
endif()

#find includes
find_path( TRIBOL_INCLUDE_DIRS tribol/tribol.hpp
           PATHS  ${TRIBOL_DIR}/include/
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

find_library( TRIBOL_LIBRARY NAMES tribol
              PATHS ${TRIBOL_DIR}/lib
              NO_DEFAULT_PATH
              NO_CMAKE_ENVIRONMENT_PATH
              NO_CMAKE_PATH
              NO_SYSTEM_ENVIRONMENT_PATH
              NO_CMAKE_SYSTEM_PATH)


include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set TRIBOL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(TRIBOL  DEFAULT_MSG
                                  TRIBOL_INCLUDE_DIRS
                                  TRIBOL_LIBRARY )
