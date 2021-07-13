#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    BufferOpsDevice.hpp
    ExecutableGroup.hpp
    Group.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    MappedVector.hpp
    ObjectCatalog.hpp
    ReferenceWrapper.hpp
    RestartFlags.hpp
    ConduitRestart.hpp
    wrapperHelpers.hpp
    Wrapper.hpp
    WrapperBase.hpp
    xmlWrapper.hpp
    HistoryDataSpec.hpp
    )

#
# Specify all sources
#
set(dataRepository_sources
    BufferOpsDevice.cpp
    ExecutableGroup.cpp
    Group.cpp
    ConduitRestart.cpp
    Wrapper.cpp
    WrapperBase.cpp
    xmlWrapper.cpp
    )

set (dependencyList codingUtilities )

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

if( ENABLE_PYGEOSX )
  set( dependencyList ${dependencyList} pylvarray )
endif()

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
