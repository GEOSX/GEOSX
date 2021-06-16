add_subdirectory (GEOSX_PTP)

#
# Specify solver headers
#
set( physicsSolvers_headers
     LinearSolverParameters.hpp
     NonlinearSolverParameters.hpp
     PhysicsSolverManager.hpp
     SolverBase.hpp
     fluidFlow/HybridFVMHelperKernels.hpp          
     fluidFlow/FlowSolverBase.hpp
     fluidFlow/FlowSolverBaseKernels.hpp
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
     fluidFlow/FluxKernelsHelper.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/WellSolverBase.hpp
     fluidFlow/wells/SinglePhaseWell.hpp
     fluidFlow/wells/SinglePhaseWellKernels.hpp     
     fluidFlow/wells/WellControls.hpp
     multiphysics/FlowProppantTransportSolver.hpp
     multiphysics/HydrofractureSolver.hpp
     multiphysics/LagrangianContactSolver.hpp
     multiphysics/PhaseFieldFractureSolver.hpp
     multiphysics/PoroelasticSolver.hpp
     multiphysics/PoroelasticSolverEmbeddedFractures.hpp
     multiphysics/ReservoirSolverBase.hpp
     multiphysics/SinglePhaseReservoir.hpp
     multiphysics/SinglePhasePoroelasticKernel.hpp
     multiphysics/SinglePhasePoroelasticEFEMKernel.hpp
     multiphysics/SinglePhasePoroelasticFluxKernels.hpp
     simplePDE/LaplaceBaseH1.hpp
     simplePDE/LaplaceFEM.hpp
     simplePDE/LaplaceFEMKernels.hpp
     simplePDE/LaplaceVEM.hpp
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
     wavePropagation/AcousticWaveEquationSEM.hpp
     )

#
# Specify solver sources
#
set( physicsSolvers_sources
     LinearSolverParameters.cpp
     NonlinearSolverParameters.cpp
     PhysicsSolverManager.cpp
     SolverBase.cpp
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
     fluidFlow/wells/WellControls.cpp     
     multiphysics/FlowProppantTransportSolver.cpp
     multiphysics/HydrofractureSolver.cpp
     multiphysics/LagrangianContactSolver.cpp
     multiphysics/PhaseFieldFractureSolver.cpp
     multiphysics/PoroelasticSolver.cpp
     multiphysics/PoroelasticSolverEmbeddedFractures.cpp
     multiphysics/ReservoirSolverBase.cpp
     multiphysics/SinglePhaseReservoir.cpp
     multiphysics/SinglePhasePoroelasticFluxKernels.cpp
     simplePDE/LaplaceBaseH1.cpp
     simplePDE/LaplaceFEM.cpp
     simplePDE/LaplaceVEM.cpp
     simplePDE/PhaseFieldDamageFEM.cpp
     solidMechanics/SolidMechanicsEmbeddedFractures.cpp
     solidMechanics/SolidMechanicsLagrangianFEM.cpp
     solidMechanics/SolidMechanicsLagrangianSSLE.cpp
     surfaceGeneration/SurfaceGenerator.cpp
     surfaceGeneration/EmbeddedSurfaceGenerator.cpp
     wavePropagation/AcousticWaveEquationSEM.cpp
     )

set( dependencyList constitutive mesh linearAlgebra discretizationMethods events )


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

if( ENABLE_PVTPackage )
    list( APPEND physicsSolvers_headers
         fluidFlow/CompositionalMultiphaseBase.hpp
         fluidFlow/CompositionalMultiphaseFVM.hpp
         fluidFlow/CompositionalMultiphaseHybridFVM.hpp
         fluidFlow/CompositionalMultiphaseBaseKernels.hpp
         fluidFlow/CompositionalMultiphaseFVMKernels.hpp
         fluidFlow/CompositionalMultiphaseHybridFVMKernels.hpp
         fluidFlow/wells/CompositionalMultiphaseWell.hpp
         fluidFlow/wells/CompositionalMultiphaseWellKernels.hpp
         multiphysics/MultiphasePoroelasticSolver.hpp
         multiphysics/CompositionalMultiphaseReservoir.hpp
         multiphysics/MultiphasePoroelasticKernel.hpp
         )

    list( APPEND physicsSolvers_sources
         fluidFlow/CompositionalMultiphaseBase.cpp
         fluidFlow/CompositionalMultiphaseFVM.cpp
         fluidFlow/CompositionalMultiphaseHybridFVM.cpp
         fluidFlow/CompositionalMultiphaseBaseKernels.cpp
         fluidFlow/CompositionalMultiphaseFVMKernels.cpp
         fluidFlow/CompositionalMultiphaseHybridFVMKernels.cpp
         multiphysics/CompositionalMultiphaseReservoir.cpp
         multiphysics/MultiphasePoroelasticSolver.cpp
         fluidFlow/wells/CompositionalMultiphaseWell.cpp
         )
endif()

blt_add_library( NAME                  physicsSolvers
                 SOURCES               ${physicsSolvers_sources}
                 HEADERS               ${physicsSolvers_headers}
                 DEPENDS_ON            ${dependencyList} ${externalComponentDeps}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )
if( externalComponentDeps )
  target_include_directories( physicsSolvers PUBLIC ${CMAKE_SOURCE_DIR}/externalComponents)
endif()


geosx_add_code_checks( PREFIX physicsSolvers )
