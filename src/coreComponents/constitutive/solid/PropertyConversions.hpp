/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file PropertyConversions.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP
#define GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP

namespace geosx
{

namespace constitutive
{

/// @namespace Namespace to collect common property conversion functions (elastic, poroelastic, etc.)
namespace conversions
{

/// @namespace Bulk modulus and shear modulus as input
namespace BulkModAndShearMod
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Young's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungMod( real64 const & K, real64 const & G )
{
  return 9 * K * G / ( 3 * K + G );
}

/**
 * @brief Compute Poisson's ratio
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return Poisson's ratio
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & K, real64 const & G )
{
  return ( 3 * K - 2 * G ) / ( 2 * ( 3 * K + G ) );
}

/**
 * @brief Compute First Lamé parameter
 * @param[in] K Bulk modulus
 * @param[in] G Shear modulus
 * @return First Lamé parameter
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toFirstLame( real64 const & K, real64 const & G )
{
  return K - 2 * G / 3;
}

} /* namespace BulkModeAndShearMod */

/// @namespace Young's modulus and Poisson ratio as input
namespace YoungModAndPoissonRatio
{

/**
 * @brief Compute bulk modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & E, real64 const & nu )
{
  return E / (3 * ( 1 - 2*nu ) );
}

/**
 * @brief Compute bulk modulus
 * @param[in] E Young's modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & E, real64 const & nu )
{
  return E / (2 * ( 1 + nu ) );
}

} /* namespace YoungModAndPoissonRatio*/

/// @namespace Shear modulus and Poisson's ratio as input
namespace ShearModAndPoissonRatio
{

/**
 * @brief Compute bulk modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & G, real64 const & nu )
{
  return 2 * G * ( 1 + nu ) / ( 3 * ( 1 - 2 * nu ) );
}

/**
 * @brief Compute Young's modulus
 * @param[in] G Shear modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungMod( real64 const & G, real64 const & nu )
{
  return 2 * G * ( 1 + nu );
}

} /* namespace ShearModAndPoissonRatio*/

/// @namespace Bulk modulus and Poisson's ratio as input
namespace BulkModAndPoissonRatio
{

/**
 * @brief Compute Young's modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Young's modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toYoungMod( real64 const & K, real64 const & nu )
{
  return 3 * K * ( 1 - 2 * nu );
}

/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] nu Poisson's ratio
 * @return Shear modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & K, real64 const & nu )
{
  return 3 * K * ( 1 - 2 * nu) / ( 2 * ( 1 + nu ) );
}

} /* namespace BulkModAndPoissonRatio */

/// @namespace Bulk modulus and Young's modulus
namespace BulkModAndYoungMod
{

/**
 * @brief Compute Shear modulus
 * @param[in] K Bulk modulus
 * @param[in] E Young's ratio
 * @return Shear modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toShearMod( real64 const & K, real64 const & E )
{
  return 3 * K * E / ( 9 * K - E );
}

/**
 * @brief Compute Poisson ratio
 * @param[in] K Bulk modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & K, real64 const & E )
{
  return ( 3 * K - E ) / ( 6 * K);
}

} /* namespace BulkModAndYoungMod */

/// @namespace Shear modulus and Young's modulus
namespace ShearModAndYoungMod
{
/**
 * @brief Compute Poisson ratio
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Poisson's ratio
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toPoissonRatio( real64 const & G, real64 const & E )
{
  return 0.5 * E / G - 1.0;
}

/**
 * @brief Compute Bulk modulus
 * @param[in] G Shear modulus
 * @param[in] E Young's modulus
 * @return Bulk modulus
 */
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
real64 toBulkMod( real64 const & G, real64 const & E )
{
  return E * G / ( 3 * ( 3 * G - E ) );
}

} /* namespace ShearModAndYoungMod*/

} /* namespace conversions */

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_PROPERTYCONVERSIONS_HPP */
