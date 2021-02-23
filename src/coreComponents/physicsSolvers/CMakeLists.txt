add_subdirectory (GEOSX_PTP)

#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     fluidFlow/CompositionalMultiphaseFlow.hpp
     fluidFlow/CompositionalMultiphaseFlowKernels.hpp
     fluidFlow/IsothermalCompositionalMultiphaseFlowKernels.hpp
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/ProppantTransport.hpp
     fluidFlow/ProppantTransportKernels.hpp
     fluidFlow/SinglePhaseBase.hpp
     fluidFlow/SinglePhaseFVM.hpp
     fluidFlow/SinglePhaseHybridFVM.hpp     
     fluidFlow/SinglePhaseBaseKernels.hpp
     fluidFlow/SinglePhaseFVMKernels.hpp          
     fluidFlow/SinglePhaseHybridFVMKernels.hpp     
     fluidFlow/SinglePhaseProppantBase.hpp
     fluidFlow/SinglePhaseProppantBaseKernels.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp     
     fluidFlow/wells/CompositionalMultiphaseWell.hpp
     fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
     fluidFlow/wells/WellControls.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/LagrangianContactSolver.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/PoroelasticSolver.hpp
     multiphysics/PoroelasticSolverEmbeddedFractures.hpp
     multiphysics/ReservoirSolverBase.hpp
     multiphysics/SinglePhaseReservoir.hpp
     multiphysics/CompositionalMultiphaseReservoir.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/PhaseFieldDamageFEM.hpp
     solidMechanics/SolidMechanicsEmbeddedFractures.hpp
     solidMechanics/SolidMechanicsLagrangianFEM.hpp
     solidMechanics/SolidMechanicsLagrangianSSLE.hpp
     solidMechanics/SolidMechanicsLagrangianFEMKernels.hpp
     solidMechanics/SolidMechanicsEFEMKernels.hpp
     solidMechanics/SolidMechanicsEFEMKernelsHelper.hpp
     solidMechanics/SolidMechanicsPoroElasticKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainQuasiStaticKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainImplicitNewmarkKernel.hpp
     solidMechanics/SolidMechanicsSmallStrainExplicitNewmarkKernel.hpp
     surfaceGeneration/SurfaceGenerator.hpp
     surfaceGeneration/EmbeddedSurfaceGenerator.hpp
     )

#
# Specify solver sources
#
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
     fluidFlow/CompositionalMultiphaseFlow.cpp
     fluidFlow/CompositionalMultiphaseFlowKernels.cpp
     fluidFlow/IsothermalCompositionalMultiphaseFlowKernels.cpp     
     fluidFlow/FlowSolverBase.cpp
     fluidFlow/ProppantTransport.cpp
     fluidFlow/ProppantTransportKernels.cpp
     fluidFlow/SinglePhaseBase.cpp
     fluidFlow/SinglePhaseBaseKernels.cpp
     fluidFlow/SinglePhaseFVM.cpp
     fluidFlow/SinglePhaseFVMKernels.cpp
     fluidFlow/SinglePhaseHybridFVM.cpp     
     fluidFlow/SinglePhaseProppantBase.cpp     
     fluidFlow/wells/WellSolverBase.cpp
     fluidFlow/wells/SinglePhaseWell.cpp          
     fluidFlow/wells/CompositionalMultiphaseWell.cpp     
     fluidFlow/wells/WellControls.cpp     
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/LagrangianContactSolver.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/PoroelasticSolver.cpp
     multiphysics/PoroelasticSolverEmbeddedFractures.cpp
     multiphysics/ReservoirSolverBase.cpp
     multiphysics/SinglePhaseReservoir.cpp
     multiphysics/CompositionalMultiphaseReservoir.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp
     solidMechanics/SolidMechanicsEmbeddedFractures.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     )

if( BUILD_OBJ_LIBS)
  set( dependencyList common constitutive dataRepository linearAlgebra )
else()
  set( dependencyList common linearAlgebra )
endif()


set( externalComponentDeps "" )
if( ENABLE_GEOSX_PTP )
  message(STATUS "Adding GEOSX_PTP to physics solver dependencies")
  list( APPEND externalComponentDeps GEOSX_PTP )
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            ${dependencyList} ${externalComponentDeps}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)
endif()


geosx_add_code_checks( PREFIX physicsSolvers )

add_subdirectory( fluidFlow/unitTests )
add_subdirectory( fluidFlow/wells/unitTests )
