/**
 * @file GeosxConfig.hpp
 *
 * GEOSX build configuration file.
 * Contains a CMake-generated list of macros that define a particular build configuration.
 */

#ifndef GEOSX_COMMON_CONFIG_HPP
#define GEOSX_COMMON_CONFIG_HPP

/// GEOSX major version number
#define GEOSX_VERSION_MAJOR 0
/// GEOSX minor version number
#define GEOSX_VERSION_MINOR 0
/// GEOSX patch version number
#define GEOSX_VERSION_PATCH 31
/// GEOSX full version number string
#define GEOSX_VERSION_FULL  "0.0.31"

/// Enables floating point execptions
#define GEOSX_USE_FPE

/// Enables bounds check in LvArray classes (CMake option ARRAY_BOUNDS_CHECK)
#define GEOSX_USE_ARRAY_BOUNDS_CHECK

/// Enables use of Caliper (CMake option ENABLE_CALIPER)
#define GEOSX_USE_CALIPER

/// Enables use of CHAI (CMake option ENABLE_CHAI)
#define GEOSX_USE_CHAI

/// Enables use of Mathpresso library (CMake option ENABLE_MATHPRESSO)
#define GEOSX_USE_MATHPRESSO

/// Enables use of MPI (CMake option ENABLE_MPI)
#define GEOSX_USE_MPI

/// Enables use of OpenMP (CMake option ENABLE_OPENMP)
#define GEOSX_USE_OPENMP

/// Enables use of CUDA (CMake option ENABLE_CUDA)
#define GEOSX_USE_CUDA

/// Enables use of Python (CMake option ENABLE_PYTHON)
#define GEOSX_USE_PYTHON

/// Enables use of GEOSX PTP module (CMake option ENABLE_GEOSX_PTP)
#define USE_GEOSX_PTP

/// Enables use of RAJA (CMake option ENABLE_RAJA)
#define GEOSX_USE_RAJA

/// Enables use of sys/time.h based timers (CMake option ENABLE_TIMERS)
#define GEOSX_USE_TIMERS

/// Enables use of additional debugging interface for TotalView (Cmake option ENABLE_TOTALVIEW_OUTPUT)
#define GEOSX_USE_TOTALVIEW_OUTPUT

/// Enables use of Intel MKL (CMake option ENABLE_MKL)
#define GEOSX_USE_MKL

/// Enables use of Trilinos library (CMake option ENABLE_TRILINOS)
#define GEOSX_USE_TRILINOS

/// Enables use of Hypre library (CMake option ENABLE_HYPRE)
#define GEOSX_USE_HYPRE

/// Enables use of PETSc library (CMake option ENABLE_PETSC)
#define GEOSX_USE_PETSC

/// Choice of global linear algebra interface (CMake option GEOSX_LA_INTERFACE)
#define GEOSX_LA_INTERFACE Hypre
/// Macro defined when Trilinos interface is selected
/* #undef GEOSX_LA_INTERFACE_TRILINOS */
/// Macro defined when Hypre interface is selected
#define GEOSX_LA_INTERFACE_HYPRE
/// Macro defined when PETSc interface is selected
/* #undef GEOSX_LA_INTERFACE_PETSC */

/// Platform-dependent mangling of fortran function names (CMake option FORTRAN_MANGLE_NO_UNDERSCORE)
#define FORTRAN_MANGLE_NO_UNDERSCORE

/// USE OF SEPARATION COEFFICIENT IN FRACTURE FLOW
#define GEOSX_USE_SEPARATION_COEFFICIENT

/// CMake option CMAKE_BUILD_TYPE
#define GEOSX_CMAKE_BUILD_TYPE "Release"

#endif  /* GEOSX_CONFIG_HPP */

