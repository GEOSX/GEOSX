/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2020 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2020 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2020 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All rights reserved
 *
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file InternalWellboreGenerator.hpp
 */

#ifndef GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP
#define GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP

#include "codingUtilities/EnumStrings.hpp"
#include "dataRepository/Group.hpp"
#include "InternalMeshGenerator.hpp"

namespace geosx
{

/**
 * @class InternalWellboreGenerator
 * @brief The InternalWellboreGenerator class is a class generating internal wellbore mesh.
 */
class InternalWellboreGenerator : public InternalMeshGenerator
{
public:

  /**
   * @brief Main constructor for InternalWellboreGenerator.
   * @param[in] name of the InternalWellboreGenerator
   * @param[in] parent point to the parent Group of the InternalWellboreGenerator
   */
  InternalWellboreGenerator( string const & name, Group * const parent );

  ~InternalWellboreGenerator() override = default;

  /**
   * @brief Return the name of the InternalWellboreGenerator in object Catalog.
   * @return string that contains the key name to InternalWellboreGenerator in the Catalog
   */
  static string catalogName() { return "InternalWellbore"; }

protected:

  void reduceNumNodesForPeriodicBoundary( SpatialPartition & partition,
                                          integer ( &numNodes )[3] ) override final;

  void setNodeGlobalIndicesOnPeriodicBoundary( SpatialPartition & partition,
                                               int ( & index )[3] ) override final;

  void setConnectivityForPeriodicBoundaries( int ( & globalIJK )[3],
                                             integer const ( &numNodesInDir )[3],
                                             int const ( &firstElemIndexInPartition )[3],
                                             localIndex ( &nodeOfBox )[8] ) override final;

  void coordinateTransformation( array2d< real64, nodes::REFERENCE_POSITION_PERM > & X,
                                 std::map< string, SortedArray< localIndex > > & nodeSets ) override final;

  virtual inline bool isCartesian() const override final
  {
    return false;
  }

  ///@cond DO_NOT_DOCUMENT
  struct viewKeyStruct : public InternalMeshGenerator::viewKeyStruct
  {
    constexpr static char const * radiusString() { return "radius"; }
    constexpr static char const * thetaString() { return "theta"; }
    constexpr static char const * rOutString() { return "rOut"; }
    constexpr static char const * rElemsString() { return "nr"; }
    constexpr static char const * tElemsString() { return "nt"; }
    constexpr static char const * rBiasString() { return "rBias"; }
    constexpr static char const * trajectoryString() { return "trajectory"; }
    constexpr static char const * cartesianOuterBoundaryString() { return "useCartesianOuterBoundary"; }
    constexpr static char const * autoSpaceRadialElemsString() { return "autoSpaceRadialElems"; }
  };
  /// @endcond

  void postProcessInput() override final;

private:

  /// Trajectory defined by coordinates of the wellbore centers
  array2d< real64 > m_trajectory;

  /// Parameter defining whether or not the outer boundary is cartesian
  int m_cartesianOuterBoundary;

  /// Parameter defining whether or not the full wellbore geometry is considered: max theta = 360
  bool m_isFullAnnulus;

  /// Parameter for automatic increasing element size in the radial direction
  array1d< int > m_autoSpaceRadialElems;

  /// Radial coordinates defining the wellbore geometry
  array1d< real64 > & m_radialCoords;

};

} /* namespace geosx */

#endif /* GEOSX_MESHUTILITIES_INTERNALWELLBOREGENERATOR_HPP */
