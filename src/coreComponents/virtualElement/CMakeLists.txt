#
# Specify all headers
#
set( virtualElement_headers
  ConformingVirtualElementOrder1.hpp
  ConformingVirtualElementOrder1_impl.hpp
  VirtualElementBase.hpp
  )
#
# Specify all sources
#
set( virtualElement_sources
   )

set( dependencyList mesh finiteElement )

blt_add_library( NAME                  virtualElement
#                 SOURCES               ${virtualElement_sources}
                 HEADERS               ${virtualElement_headers}
                 DEPENDS_ON            ${dependencyList}
#                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

#target_include_directories( virtualElement PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks( PREFIX virtualElement )
