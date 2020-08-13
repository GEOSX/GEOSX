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
 * @file ExtrinsicMeshData.hpp
 */

#ifndef GEOSX_EXTRINSIC_MESH_DATA_HPP_
#define GEOSX_EXTRINSIC_MESH_DATA_HPP_

#include "codingUtilities/traits.hpp"
#include "dataRepository/RestartFlags.hpp"
#include "common/DataTypes.hpp"

/**
 * @brief Generates a traits struct.
 * @param NAME Name of the traits struct.
 * @param KEY The string literal that will be used as the key to register and
 *   lookup the data in the repository.
 * @param TYPE The type of data that will be registered.
 * @param DEFAULT The default value for the data.
 * @param PLOTLEVEL The default plot level for the wrapper that will contain the data.
 * @param RESTARTFLAG The default restart flag for the wrapper that contains the data.
 * @param DESCRIPTION A string literal that contains a description of the data for
 *   use in sphinx documentation.
 */
#define EXTRINSIC_MESH_DATA_TRAIT( NAME,                                             \
                                   KEY,                                              \
                                   TYPE,                                             \
                                   DEFAULT,                                          \
                                   PLOTLEVEL,                                        \
                                   RESTARTFLAG,                                      \
                                   DESCRIPTION )                                     \
  /** @struct NAME                                                                */ \
  /** @brief Trait struct for NAME data                                           */ \
  struct NAME                                                                        \
  {                                                                                  \
    /** The key for registration with the data repository.                      */   \
    static constexpr auto key = KEY;                                                 \
    /** The actual type to be registered.                                       */   \
    using type = TYPE;                                                               \
    /** The template type T for registration of a container<T>.                 */   \
    using dataType = internal::typeHelper_t< TYPE >;                                 \
    /** The dataRepository::DefaultValue for NAME.                              */   \
    static constexpr dataType defaultValue = DEFAULT;                                \
    /** The default dataRepository::PlotLevel for NAME.                         */   \
    static constexpr auto plotLevel = dataRepository::PlotLevel::PLOTLEVEL;          \
    /** The default dataRepository::RestartFlags for NAME.                      */   \
    static constexpr auto restartFlag =                                              \
      dataRepository::RestartFlags::RESTARTFLAG;                                     \
    /** Description of the NAME data for use in sphinx documentation            */   \
    static constexpr auto description = DESCRIPTION;                                 \
  }

namespace geosx
{
/**
 * A scope for extrinsic mesh data traits.
 */
namespace extrinsicMeshData
{
namespace internal
{
template< typename TYPE,
          bool HAS_TYPE =
            std::enable_if_t< traits::HasAlias_value_type< TYPE >, std::true_type >::value >
struct typeHelper
{
  using type = typename TYPE::value_type;
};

template< typename TYPE >
struct typeHelper< TYPE, false >
{
  using type =
    TYPE;  //typename std::enable_if< std::is_fundamental<TYPE>::value, TYPE>::type;
};

template< typename T >
using typeHelper_t = typename typeHelper< T >::type;
}  // namespace internal

EXTRINSIC_MESH_DATA_TRAIT(
  ParentIndex,
  "parentIndex",
  array1d< localIndex >,
  -1,
  LEVEL_2,
  WRITE_AND_READ,
  "Index of parent within the mesh object it is registered on." );

EXTRINSIC_MESH_DATA_TRAIT(
  ChildIndex,
  "childIndex",
  array1d< localIndex >,
  -1,
  LEVEL_2,
  WRITE_AND_READ,
  "Index of child within the mesh object it is registered on." );

EXTRINSIC_MESH_DATA_TRAIT(
  DegreeFromCrack,
  "degreeFromCrack",
  array1d< integer >,
  -1,
  LEVEL_1,
  WRITE_AND_READ,
  "Distance to the crack in terms of topological distance. "
  "(i.e. how many nodes are along the path to the closest "
  "node that is on the crack surface." );

EXTRINSIC_MESH_DATA_TRAIT(
  DegreeFromCrackTip,
  "degreeFromCrackTip",
  array1d< integer >,
  100000,
  LEVEL_1,
  WRITE_AND_READ,
  "Distance to the crack tip in terms of topological distance. "
  "(i.e. how many nodes are along the path to the closest "
  "node that is on the crack surface." );

EXTRINSIC_MESH_DATA_TRAIT( SIFNode,
                           "SIFNode",
                           array1d< real64 >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Calculated Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( RuptureTime,
                           "ruptureTime",
                           array1d< real64 >,
                           1.0e9,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Time that the object was ruptured/split." );

EXTRINSIC_MESH_DATA_TRAIT(
  RuptureRate,
  "ruptureRate",
  array1d< real64 >,
  1.0e99,
  LEVEL_0,
  WRITE_AND_READ,
  "Rate of rupture in terms of number of objects split per time." );

EXTRINSIC_MESH_DATA_TRAIT(
  SIF_I,
  "SIF_I",
  array1d< real64 >,
  -1,
  LEVEL_1,
  WRITE_AND_READ,
  "Calculated mode 1 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT(
  SIF_II,
  "SIF_II",
  array1d< real64 >,
  -1,
  LEVEL_1,
  WRITE_AND_READ,
  "Calculated mode 2 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT(
  SIF_III,
  "SIF_III",
  array1d< real64 >,
  -1,
  LEVEL_1,
  WRITE_AND_READ,
  "Calculated mode 3 Stress Intensity Factor on the node." );

EXTRINSIC_MESH_DATA_TRAIT( RuptureState,
                           "ruptureState",
                           array1d< integer >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Rupture state of the face: \n 0=not ready for "
                           "rupture \n 1=ready for rupture \n 2=ruptured." );

EXTRINSIC_MESH_DATA_TRAIT( SIFonFace,
                           "SIFonFace",
                           array1d< real64 >,
                           1,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "Calculated Stress Intensity Factor on the face." );

/// The template type T for registration of a container<T>.

/** @struct K_IC
 * @brief Trait struct for K_IC data
 */
struct K_IC
{
  /// The key for registration with the data repository.
  static constexpr auto key = "K_IC";

  /// The template type T for registration of a container<T>.
  using dataType = R1Tensor;

  /// The actual type to be registered.
  using type = array1d< dataType >;

  /// The dataRepository::DefaultValue for K_IC.
  static constexpr real64 defaultValue = 1e99;

  /// The default dataRepository::PlotLevel for K_IC.
  static constexpr auto plotLevel = dataRepository::PlotLevel::LEVEL_0;

  /// The default dataRepository::RestartFlags for K_IC.
  static constexpr auto restartFlag =
    dataRepository::RestartFlags::WRITE_AND_READ;

  /// Description of the K_IC data for use in sphinx documentation
  static constexpr auto description =
    "Critical Stress Intensity Factor :math:`K_{IC}` in the plane of the face.";
};

EXTRINSIC_MESH_DATA_TRAIT( K_IC_00,
                           "K_IC_00",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_01,
                           "K_IC_01",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_02,
                           "K_IC_02",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 0-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_10,
                           "K_IC_10",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_11,
                           "K_IC_11",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_12,
                           "K_IC_12",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 1-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_20,
                           "K_IC_20",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 0-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_21,
                           "K_IC_21",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 1-direction." );

EXTRINSIC_MESH_DATA_TRAIT( K_IC_22,
                           "K_IC_22",
                           array1d< real64 >,
                           -1,
                           NOPLOT,
                           WRITE_AND_READ,
                           ":math:`K_{IC}` on 2-plane, in 2-direction." );

EXTRINSIC_MESH_DATA_TRAIT( PrimaryCandidateFace,
                           "primaryCandidateFace",
                           array1d< localIndex >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "??" );

EXTRINSIC_MESH_DATA_TRAIT( IsFaceSeparable,
                           "isFaceSeparable",
                           array1d< integer >,
                           0,
                           LEVEL_0,
                           WRITE_AND_READ,
                           "A flag to mark if the face is separable." );

}  // namespace extrinsicMeshData
}  // namespace geosx

#endif /* GEOSX_MESH_MESHFIELDS_HPP_ */
