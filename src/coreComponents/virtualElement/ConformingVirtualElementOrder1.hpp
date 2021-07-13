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
 * @file ConformingVirtualElementOrder1.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
#define GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_

#include "VirtualElementBase.hpp"

namespace geosx
{
namespace virtualElement
{
template< localIndex MAXCELLNODES, localIndex MAXFACENODES >
class ConformingVirtualElementOrder1 final : public VirtualElementBase
{
public:
  using InputNodeCoords = arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD >;
  using InputCellToNodeMap = arrayView2d< localIndex const, cells::NODE_MAP_USD >;
  using InputCellToFaceMap = arrayView2d< localIndex const >;
  using InputFaceToNodeMap = ArrayOfArraysView< localIndex const >;
  using InputFaceToEdgeMap = ArrayOfArraysView< localIndex const >;
  using InputEdgeToNodeMap = arrayView2d< localIndex const >;

  static constexpr localIndex maxSupportPoints = MAXCELLNODES;
  static constexpr localIndex numQuadraturePoints = 1;

  struct BasisData
  {
    static constexpr localIndex maxSupportPoints = MAXCELLNODES;
    static constexpr localIndex numQuadraturePoints = 1;
    localIndex numSupportPoints;
    real64 quadratureWeight;
    real64 basisFunctionsIntegralMean[MAXCELLNODES];
    real64 stabilizationMatrix[MAXCELLNODES][MAXCELLNODES];
    real64 basisDerivativesIntegralMean[MAXCELLNODES][3];
  };

private:
  GEOSX_HOST_DEVICE
  static void
    computeFaceIntegrals( InputNodeCoords const & nodesCoords,
                          localIndex const (&faceToNodes)[MAXFACENODES],
                          localIndex const (&faceToEdges)[MAXFACENODES],
                          localIndex const & numFaceVertices,
                          real64 const & faceArea,
                          real64 const (&faceCenter)[3],
                          real64 const (&faceNormal)[3],
                          InputEdgeToNodeMap const & edgeToNodes,
                          real64 const & invCellDiameter,
                          real64 const (&cellCenter)[3],
                          real64 ( &basisIntegrals )[MAXFACENODES],
                          real64 ( &threeDMonomialIntegrals )[3] );

public:

  BasisData m_basisData;

  ConformingVirtualElementOrder1() = default;

  virtual ~ConformingVirtualElementOrder1() = default;

  GEOSX_HOST_DEVICE
  static void
  computeProjectors( localIndex const & cellIndex,
                     InputNodeCoords const & nodesCoords,
                     InputCellToNodeMap const & cellToNodeMap,
                     InputCellToFaceMap const & elementToFaceMap,
                     InputFaceToNodeMap const & faceToNodeMap,
                     InputFaceToEdgeMap const & faceToEdgeMap,
                     InputEdgeToNodeMap const & edgeToNodeMap,
                     arrayView2d< real64 const > const faceCenters,
                     arrayView2d< real64 const > const faceNormals,
                     arrayView1d< real64 const > const faceAreas,
                     real64 const (&cellCenter)[3],
                     real64 const & cellVolume,
                     BasisData & basisData
                     );

  GEOSX_HOST_DEVICE
  void processLocalGeometry( localIndex const & cellIndex,
                             InputNodeCoords const & nodesCoords,
                             InputCellToNodeMap const & cellToNodeMap,
                             InputCellToFaceMap const & elementToFaceMap,
                             InputFaceToNodeMap const & faceToNodeMap,
                             InputFaceToEdgeMap const & faceToEdgeMap,
                             InputEdgeToNodeMap const & edgeToNodeMap,
                             arrayView2d< real64 const > const faceCenters,
                             arrayView2d< real64 const > const faceNormals,
                             arrayView1d< real64 const > const faceAreas,
                             real64 const (&cellCenter)[3],
                             real64 const & cellVolume
                             )
  {
    computeProjectors( cellIndex,
                       nodesCoords,
                       cellToNodeMap,
                       elementToFaceMap,
                       faceToNodeMap,
                       faceToEdgeMap,
                       edgeToNodeMap,
                       faceCenters,
                       faceNormals,
                       faceAreas,
                       cellCenter,
                       cellVolume,
                       m_basisData );
  }

  GEOSX_HOST_DEVICE
  static constexpr localIndex getMaxSupportPoints()
  {
    return maxSupportPoints;
  }

  GEOSX_HOST_DEVICE
  localIndex getNumQuadraturePoints() const override
  {
    return numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  localIndex getNumSupportPoints() const override
  {
    return m_basisData.numSupportPoints;
  }

  GEOSX_HOST_DEVICE
  void calcN( localIndex const q,
              real64 ( & N )[maxSupportPoints] ) const
  {
    GEOSX_UNUSED_VAR( q );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      N[i] = m_basisData.basisFunctionsIntegralMean[i];
    }
  }

  GEOSX_HOST_DEVICE
  real64 calcGradN( localIndex const q,
                    real64 ( & gradN )[maxSupportPoints][3] ) const
  {
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        gradN[i][j] = m_basisData.basisDerivativesIntegralMean[i][j];
      }
    }
    return transformedQuadratureWeight( q );
  }

  GEOSX_HOST_DEVICE
  real64 transformedQuadratureWeight( localIndex const GEOSX_UNUSED_PARAM( q ) ) const
  {
    return m_basisData.quadratureWeight;
  }

  GEOSX_HOST_DEVICE
  real64 calcStabilizationValue( localIndex const iBasisFunction,
                                 localIndex const jBasisFunction ) const
  {
    return m_basisData.stabilizationMatrix[iBasisFunction][jBasisFunction];
  }

  GEOSX_HOST_DEVICE
  static localIndex getMaxSupportPoints( BasisData const & basisData )
  {
    return basisData.maxSupportPoints;
  }

  GEOSX_HOST_DEVICE
  static localIndex getNumQuadraturePoints( BasisData const & basisData )
  {
    return basisData.numQuadraturePoints;
  }

  GEOSX_HOST_DEVICE
  static localIndex getNumSupportPoints( BasisData const & basisData )
  {
    return basisData.numSupportPoints;
  }

  GEOSX_HOST_DEVICE
  static real64 calcStabilizationValue( localIndex const iBasisFunction,
                                        localIndex const jBasisFunction,
                                        BasisData const & basisData )
  { return basisData.stabilizationMatrix[iBasisFunction][jBasisFunction]; }

  GEOSX_HOST_DEVICE
  static real64 calcGradN( localIndex const q,
                           BasisData const & basisData,
                           real64 ( & gradN )[maxSupportPoints][3] )
  {
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      for( localIndex j = 0; j < 3; ++j )
      {
        gradN[i][j] = basisData.basisDerivativesIntegralMean[i][j];
      }
    }
    return transformedQuadratureWeight( q, basisData );
  }

  GEOSX_HOST_DEVICE
  static void calcN( localIndex const q,
                     BasisData const & basisData,
                     real64 ( & N )[maxSupportPoints] )
  {
    GEOSX_UNUSED_VAR( q );
    for( localIndex i = 0; i < maxSupportPoints; ++i )
    {
      N[i] = basisData.basisFunctionsIntegralMean[i];
    }
  }

  GEOSX_HOST_DEVICE
  static real64 transformedQuadratureWeight( localIndex const GEOSX_UNUSED_PARAM( q ),
                                             BasisData const & basisData )
  {
    return basisData.quadratureWeight;
  }


public:

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE >
  GEOSX_HOST_DEVICE
  static real64 computeDiameter( POINT_COORDS_TYPE points,
                                 localIndex const & numPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[numPoint][i] - points[numOthPoint][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
        {
          diameter = candidateDiameter;
        }
      }
      diameter = LvArray::math::sqrt< real64 >( diameter );
    }
    return diameter;
  }

  template< localIndex DIMENSION, typename POINT_COORDS_TYPE, typename POINT_SELECTION_TYPE >
  GEOSX_HOST_DEVICE
  static real64 computeDiameter( POINT_COORDS_TYPE points,
                                 POINT_SELECTION_TYPE selectedPoints,
                                 localIndex const & numSelectedPoints )
  {
    real64 diameter = 0;
    for( localIndex numPoint = 0; numPoint < numSelectedPoints; ++numPoint )
    {
      for( localIndex numOthPoint = 0; numOthPoint < numPoint; ++numOthPoint )
      {
        real64 candidateDiameter = 0.0;
        for( localIndex i = 0; i < DIMENSION; ++i )
        {
          real64 coordDiff = points[selectedPoints[numPoint]][i] -
                             points[selectedPoints[numOthPoint]][i];
          candidateDiameter += coordDiff * coordDiff;
        }
        if( diameter < candidateDiameter )
        {
          diameter = candidateDiameter;
        }
      }
      diameter = LvArray::math::sqrt< real64 >( diameter );
    }
    return diameter;
  }

};
}
}

#include "ConformingVirtualElementOrder1_impl.hpp"

#endif // GEOSX_VIRTUALELEMENT_CONFORMINGVIRTUALELEMENTORDER1_HPP_
