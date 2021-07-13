

#
# Specify all headers
#
set( fileIO_headers
     Outputs/OutputManager.hpp
     Outputs/OutputBase.hpp
     Outputs/SiloOutput.hpp
     Outputs/RestartOutput.hpp
     Outputs/TimeHistoryOutput.hpp
     Outputs/BlueprintOutput.hpp
     Outputs/PythonOutput.hpp
     silo/SiloFile.hpp
     timeHistory/TimeHistHDF.hpp
     timeHistory/TimeHistoryCollection.hpp
     timeHistory/PackCollection.hpp
     timeHistory/HistoryIO.hpp
     )

#
# Specify all sources
#
set( fileIO_sources
     Outputs/OutputManager.cpp
     Outputs/OutputBase.cpp
     Outputs/SiloOutput.cpp
     Outputs/RestartOutput.cpp
     Outputs/TimeHistoryOutput.cpp
     Outputs/BlueprintOutput.cpp
     Outputs/PythonOutput.cpp
     silo/SiloFile.cpp
     timeHistory/PackCollection.cpp
     timeHistory/TimeHistHDF.cpp
      )

set( dependencyList mesh constitutive silo hdf5 )


if( ENABLE_MPI )
  add_subdirectory(coupling/hdf5_interface)
  set( dependencyList ${dependencyList} mpi hdf5_interface)
  list( APPEND fileIO_headers coupling/ChomboCoupler.hpp Outputs/ChomboIO.hpp)
  list( APPEND fileIO_sources coupling/ChomboCoupler.cpp Outputs/ChomboIO.cpp)
endif() 

if( ENABLE_VTK )
  list( APPEND fileIO_headers
               vtk/VTKPVDWriter.hpp
               vtk/VTKVTMWriter.hpp
               vtk/VTKPolyDataWriterInterface.hpp
               Outputs/VTKOutput.hpp )

  list( APPEND fileIO_sources
               vtk/VTKPVDWriter.cpp
               vtk/VTKVTMWriter.cpp
               vtk/VTKPolyDataWriterInterface.cpp
               Outputs/VTKOutput.cpp )


  list( APPEND dependencyList
               VTK::CommonCore
               VTK::CommonSystem
               VTK::CommonMisc
               VTK::IOCore
               VTK::IOXML
               VTK::expat
               VTK::doubleconversion
               VTK::lz4
               VTK::lzma
               VTK::zlib
               VTK::loguru )
endif()

if ( ENABLE_CUDA )
  list( APPEND dependencyList  cuda )
endif()

blt_add_library( NAME                  fileIO
                 SOURCES               ${fileIO_sources}
                 HEADERS               ${fileIO_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )


target_include_directories( fileIO PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX fileIO )
