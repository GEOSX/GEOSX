###############################################################################
#
# Setup WORLDS_CORE
# This file defines:
#  WORLDS_CORE_FOUND - If WORLDS_CORE was found
#  WORLDS_CORE_INCLUDE_DIRS - The WORLDS_CORE include directories
#  WORLDS_CORE_LIBRARIES - The WORLDS_CORE library

# first Check for WORLDS_CORE_DIR

if(NOT WORLDS_CORE_DIR)
    MESSAGE(FATAL_ERROR "Could not find WORLDS_CORE. WORLDS_CORE support needs explicit WORLDS_CORE_DIR")
endif()

#find includes
find_path( WORLDS_CORE_INCLUDE_DIRS worlds/WorldBase.h
           PATHS  ${WORLDS_CORE_DIR}/include
           NO_DEFAULT_PATH
           NO_CMAKE_ENVIRONMENT_PATH
           NO_CMAKE_PATH
           NO_SYSTEM_ENVIRONMENT_PATH
           NO_CMAKE_SYSTEM_PATH)

set(_path ${WORLDS_CORE_DIR}/lib/cmake/worlds_core-targets.cmake)

if(NOT EXISTS ${_path})
     MESSAGE(STATUS "Could not find WORLDS_CORE cmake include file (${_path})")
else()
     include(${_path})
     set(_targets worlds_core slide_world)
     set(WORLDS_CORE_LIBRARIES)
     foreach(_target ${_targets})
         if(TARGET ${_target})
            list(APPEND WORLDS_CORE_LIBRARIES ${_target})
         endif()
     endforeach()
     set(WORLDS_CORE_INCLUDE_DIRS ${WORLDS_CORE_DIR}/include)
endif()

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set WORLDS_CORE_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(WORLDS_CORE  DEFAULT_MSG
                                  WORLDS_CORE_INCLUDE_DIRS
                                  WORLDS_CORE_LIBRARIES )
