###############################################################################
#
# Setup RXB
# This file defines:
#  RXB_FOUND - If RXB was found
#  RXB_INCLUDE_DIRS - The RXB include directories
#  RXB_LIBRARIES - The RXB library

# first Check for RXB_DIR

if(NOT RXB_DIR)
    MESSAGE(FATAL_ERROR "Could not find RXB. RXB support needs explicit RXB_DIR")
endif()

#find includes
find_path( RXB_INCLUDE_DIRS rxb/RCB.h
           PATHS  ${RXB_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

set(_path ${RXB_DIR}/lib/cmake/rxb-targets.cmake)

if(NOT EXISTS ${_path})
     MESSAGE(STATUS "Could not find RXB cmake include file (${_path})")
else()
     include(${_path})
     set(_targets bisectree rcb rib)
     set(RXB_LIBRARIES)
     foreach(_target ${_targets})
         if(TARGET ${_target})
            list(APPEND RXB_LIBRARIES ${_target})
         endif()
     endforeach()
     set(RXB_INCLUDE_DIRS ${RXB_DIR}/include)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set RXB_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(RXB  DEFAULT_MSG
                                  RXB_INCLUDE_DIRS
                                  RXB_LIBRARIES )
