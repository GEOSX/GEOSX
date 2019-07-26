#
# Specify all headers
#
set(managers_headers
    DomainPartition.hpp
    EventManager.hpp
    Events/EventBase.hpp
    Events/PeriodicEvent.hpp
    Events/HaltEvent.hpp
    Events/SoloEvent.hpp
    Outputs/OutputManager.hpp
    Outputs/OutputBase.hpp
    Outputs/SiloOutput.hpp
    Outputs/RestartOutput.hpp
    Outputs/BlueprintOutput.hpp
    Outputs/ChomboIO.hpp
    Functions/FunctionBase.hpp
    Functions/SymbolicFunction.hpp
    Functions/TableFunction.hpp
    Functions/CompositeFunction.hpp
    Functions/NewFunctionManager.hpp
    ObjectManagerBase.hpp
    ProblemManager.hpp
    NumericalMethodsManager.hpp
    FieldSpecification/FieldSpecificationBase.hpp
    FieldSpecification/FieldSpecificationManager.hpp
    FieldSpecification/DirichletBoundaryCondition.hpp
    FieldSpecification/SourceFluxBoundaryCondition.hpp)    


#
# Specify all sources
#
set(managers_sources
  	DomainPartition.cpp
    EventManager.cpp
    Events/EventBase.cpp
    Events/PeriodicEvent.cpp
    Events/HaltEvent.cpp
    Events/SoloEvent.cpp
    Outputs/OutputManager.cpp
    Outputs/OutputBase.cpp
    Outputs/SiloOutput.cpp
    Outputs/RestartOutput.cpp
    Outputs/BlueprintOutput.cpp
    Outputs/ChomboIO.cpp
    Functions/FunctionBase.cpp
    Functions/SymbolicFunction.cpp
    Functions/TableFunction.cpp
    Functions/CompositeFunction.cpp
    Functions/NewFunctionManager.cpp
    ObjectManagerBase.cpp
    ProblemManager.cpp
    NumericalMethodsManager.cpp
    FieldSpecification/FieldSpecificationBase.cpp
    FieldSpecification/FieldSpecificationManager.cpp
    FieldSpecification/DirichletBoundaryCondition.cpp
    FieldSpecification/SourceFluxBoundaryCondition.cpp)

set( dependencyList dataRepository fileIO optionparser RAJA trilinos )
if( ENABLE_MATHPRESSO )
   set( dependencyList ${dependencyList} mathpresso )
endif()

if( ENABLE_TRIBOL )
   set( dependencyList ${dependencyList} coupling )
endif()

blt_add_library( NAME                  managers
                 SOURCES               ${managers_sources}
                 HEADERS               ${managers_headers}
                 DEPENDS_ON            ${dependencyList}
               )

target_include_directories( managers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

add_subdirectory( unitTests )

geosx_add_code_checks( PREFIX managers )
