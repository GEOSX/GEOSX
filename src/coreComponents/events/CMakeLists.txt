#
# Specify all headers
#
set(events_headers
    EventManager.hpp
    EventBase.hpp
    PeriodicEvent.hpp
    HaltEvent.hpp
    SoloEvent.hpp
    tasks/TasksManager.hpp
    tasks/TaskBase.hpp
   )


#
# Specify all sources
#

set(events_sources
    EventManager.cpp
    EventBase.cpp
    PeriodicEvent.cpp
    HaltEvent.cpp
    SoloEvent.cpp
    tasks/TaskBase.cpp
    tasks/TasksManager.cpp
   )

set( dependencyList fileIO )

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif( )


blt_add_library( NAME                  events
                 SOURCES               ${events_sources}
                 HEADERS               ${events_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( events PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX events )
