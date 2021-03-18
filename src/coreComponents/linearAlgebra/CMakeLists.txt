#
# Specify all headers
#
set( linearAlgebra_headers
     common/common.hpp
     common/LinearOperator.hpp
     common/LinearSolverBase.hpp
     common/PreconditionerBase.hpp
     common/traits.hpp
     interfaces/dense/BlasLapackFunctions.h
     interfaces/dense/BlasLapackLA.hpp
     interfaces/InterfaceTypes.hpp
     interfaces/MatrixBase.hpp
     interfaces/VectorBase.hpp
     multiscale/MultiscaleMeshData.hpp
     multiscale/MultiscaleMeshUtils.hpp
     multiscale/MetisPartitioner.hpp
     multiscale/MultiscaleMeshLevel.hpp
     multiscale/MultiscaleMeshObjectManager.hpp
     multiscale/MultiscaleStrategy.hpp
     multiscale/MsRSBStrategy.hpp
     multiscale/PartitionerBase.hpp
     multiscale/PreconditionerMultiscale.hpp
     solvers/BicgstabSolver.hpp
     solvers/BlockPreconditioner.hpp
     solvers/CgSolver.hpp
     solvers/GmresSolver.hpp
     solvers/KrylovSolver.hpp
     solvers/KrylovUtils.hpp
     solvers/PreconditionerIdentity.hpp
     solvers/PreconditionerJacobi.hpp
     solvers/PreconditionerBlockJacobi.hpp
     solvers/SeparateComponentPreconditioner.hpp
     utilities/Arnoldi.hpp
     utilities/BlockOperatorView.hpp
     utilities/BlockOperatorWrapper.hpp
     utilities/BlockOperator.hpp
     utilities/BlockVectorView.hpp
     utilities/BlockVectorWrapper.hpp
     utilities/BlockVector.hpp
     utilities/ComponentMask.hpp
     utilities/LAIHelperFunctions.hpp
     utilities/LinearSolverParameters.hpp
     utilities/LinearSolverResult.hpp
     utilities/InverseNormalOperator.hpp
     utilities/NormalOperator.hpp
     utilities/TransposeOperator.hpp
     DofManager.hpp
     DofManagerHelpers.hpp )

#
# Specify all sources
#
set( linearAlgebra_sources
     interfaces/dense/BlasLapackLA.cpp
     multiscale/MultiscaleMeshUtils.cpp
     multiscale/MetisPartitioner.cpp
     multiscale/MultiscaleMeshLevel.cpp
     multiscale/MultiscaleMeshObjectManager.cpp
     multiscale/MultiscaleStrategy.cpp
     multiscale/MsRSBStrategy.cpp
     multiscale/PartitionerBase.cpp
     multiscale/PreconditionerMultiscale.cpp
     solvers/BicgstabSolver.cpp
     solvers/BlockPreconditioner.cpp
     solvers/CgSolver.cpp
     solvers/GmresSolver.cpp
     solvers/KrylovSolver.cpp
     solvers/SeparateComponentPreconditioner.cpp
     DofManager.cpp )

set( dependencyList mesh blas lapack RAJA )


list( APPEND linearAlgebra_headers
      interfaces/direct/SuiteSparse.hpp
      interfaces/direct/SuperLUDist.hpp )

list( APPEND linearAlgebra_sources
      interfaces/direct/SuiteSparse.cpp
      interfaces/direct/SuperLUDist.cpp )

list( APPEND dependencyList suitesparse superlu_dist )

if( ENABLE_TRILINOS )

    list( APPEND linearAlgebra_headers
          interfaces/trilinos/EpetraExport.hpp
          interfaces/trilinos/EpetraMatrix.hpp
          interfaces/trilinos/EpetraVector.hpp
          interfaces/trilinos/EpetraUtils.hpp
          interfaces/trilinos/TrilinosPreconditioner.hpp
          interfaces/trilinos/TrilinosSolver.hpp
          interfaces/trilinos/TrilinosInterface.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/trilinos/EpetraExport.cpp
          interfaces/trilinos/EpetraMatrix.cpp
          interfaces/trilinos/EpetraVector.cpp
          interfaces/trilinos/TrilinosPreconditioner.cpp
          interfaces/trilinos/TrilinosSolver.cpp
          interfaces/trilinos/TrilinosInterface.cpp )

    list( APPEND dependencyList trilinos )

endif()

if( ENABLE_HYPRE )

    list( APPEND linearAlgebra_headers
          interfaces/hypre/HypreExport.hpp
          interfaces/hypre/HypreMGR.hpp
          interfaces/hypre/HypreMatrix.hpp
          interfaces/hypre/HypreVector.hpp
          interfaces/hypre/HyprePreconditioner.hpp
          interfaces/hypre/HypreSolver.hpp
          interfaces/hypre/HypreInterface.hpp
          interfaces/hypre/HypreUtils.hpp 
          interfaces/hypre/mgrStrategies/SinglePhasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/HybridSinglePhasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/MultiphasePoromechanics.hpp
          interfaces/hypre/mgrStrategies/SinglePhaseReservoirFVM.hpp	  
          interfaces/hypre/mgrStrategies/SinglePhaseHybridFVM.hpp
          interfaces/hypre/mgrStrategies/SinglePhaseReservoirHybridFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseHybridFVM.hpp
          interfaces/hypre/mgrStrategies/CompositionalMultiphaseReservoirHybridFVM.hpp
          interfaces/hypre/mgrStrategies/Hydrofracture.hpp
          interfaces/hypre/mgrStrategies/LagrangianContactMechanics.hpp)

    list( APPEND linearAlgebra_sources
          interfaces/hypre/HypreExport.cpp
          interfaces/hypre/HypreMGR.cpp
          interfaces/hypre/HypreMatrix.cpp
          interfaces/hypre/HypreVector.cpp
          interfaces/hypre/HyprePreconditioner.cpp
          interfaces/hypre/HypreSolver.cpp
          interfaces/hypre/HypreInterface.cpp )

    list( APPEND dependencyList hypre )

endif()

if( ENABLE_PETSC )

    list( APPEND linearAlgebra_headers
          interfaces/petsc/PetscExport.hpp
          interfaces/petsc/PetscVector.hpp
          interfaces/petsc/PetscMatrix.hpp
          interfaces/petsc/PetscPreconditioner.hpp
          interfaces/petsc/PetscSolver.hpp
          interfaces/petsc/PetscInterface.hpp
          interfaces/petsc/PetscUtils.hpp )

    list( APPEND linearAlgebra_sources
          interfaces/petsc/PetscExport.cpp
          interfaces/petsc/PetscVector.cpp
          interfaces/petsc/PetscMatrix.cpp
          interfaces/petsc/PetscPreconditioner.cpp
          interfaces/petsc/PetscSolver.cpp
          interfaces/petsc/PetscInterface.cpp )

    list( APPEND dependencyList petsc )

endif()

if( ENABLE_OPENMP )
    list( APPEND dependencyList openmp )
endif()

if ( ENABLE_CUDA )
    list( APPEND dependencyList cuda )
endif()

if ( ENABLE_MKL )
    list( APPEND dependencyList mkl )
endif()

if( ENABLE_MPI )
    list( APPEND dependencyList mpi )
endif()

# We depend on METIS for multiscale coarsening
if( ENABLE_METIS )
    list( APPEND dependencyList metis )
endif()

blt_add_library( NAME                  linearAlgebra
                 SOURCES               ${linearAlgebra_sources}
                 HEADERS               ${linearAlgebra_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
               )

target_include_directories( linearAlgebra PUBLIC ${CMAKE_CURRENT_LIST_DIR} )

add_subdirectory(unitTests)

geosx_add_code_checks( PREFIX linearAlgebra )
