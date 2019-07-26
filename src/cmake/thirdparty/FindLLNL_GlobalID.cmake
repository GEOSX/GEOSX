###############################################################################
#
# Setup LLNL_GlobalID
# This file defines:
#  LLNL_GLOBALID_FOUND - If LLNL_GlobalID was found
#  LLNL_GLOBALID_INCLUDE_DIRS - The LLNL_GlobalID include directories

# first Check for LLNL_GLOBALID_DIR

if(NOT LLNL_GLOBALID_DIR)
    MESSAGE(FATAL_ERROR "Could not find LLNL_GlobalID. LLNL_GlobalID support needs explicit LLNL_GLOBALID_DIR")
endif()

#find includes
find_path( LLNL_GLOBALID_INCLUDE_DIRS LLNL_GlobalID.h
           PATHS  ${LLNL_GLOBALID_DIR}/include/ ${LLNL_GLOBALID_DIR}
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set LLNL_GLOBALID_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(LLNL_GLOBALID  DEFAULT_MSG
                                  LLNL_GLOBALID_INCLUDE_DIRS )
