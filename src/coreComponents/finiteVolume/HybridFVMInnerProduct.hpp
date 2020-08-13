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
 * @file HybridFVMInnerProduct.hpp
 */

#ifndef GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
#define GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP

#include "common/DataTypes.hpp"
#include "rajaInterface/GEOS_RAJA_Interface.hpp"
#include "linearAlgebra/interfaces/InterfaceTypes.hpp"
#include "meshUtilities/ComputationalGeometry.hpp"

namespace geosx
{
namespace HybridFVMInnerProduct
{
/******************************** Helpers ********************************/

/**
 * @struct HybridFVMInnerProductHelper
 * @brief Helper struct handling inner product for hybrid finite volume schemes.
 */
struct HybridFVMInnerProductHelper
{
  /**
   * @brief Create a full tensor from an array.
   * @param[in] values the input array
   * @param[out] result the full tensor
   */
  GEOSX_HOST_DEVICE
  static void
  MakeFullTensor( real64 const ( &values )[3], real64 ( &result )[3][3] )
  {
    LvArray::tensorOps::fill< 3, 3 >( result, 0.0 );
    result[0][0] = values[0];
    result[1][1] = values[1];
    result[2][2] = values[2];
  }

  /**
   * @brief Orthonormalize a set of three vectors
   * @tparam NF vector space dimensionality
   * @param[in,out] q0 first vector
   * @param[in,out] q1 second vector
   * @param[in,out] q2 third vector
   * @param[out] cellToFaceMat a copy of in/out vectors stacked into a matrix
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE static void Orthonormalize( real64 ( &q0 )[NF],
                                                real64 ( &q1 )[NF],
                                                real64 ( &q2 )[NF],
                                                real64 ( &cellToFaceMat )[NF][3] )
  {
    // modified Gram-Schmidt algorithm

    // q0
    LvArray::tensorOps::scale< NF >( q0, 1.0 / LvArray::tensorOps::l2Norm< NF >( q0 ) );

    // q1
    real64 const q0Dotq1 = LvArray::tensorOps::AiBi< NF >( q0, q1 );
    LvArray::tensorOps::scaledAdd< NF >( q1, q0, -q0Dotq1 );
    LvArray::tensorOps::scale< NF >( q1, 1.0 / LvArray::tensorOps::l2Norm< NF >( q1 ) );

    // q2
    real64 const q0Dotq2 = LvArray::tensorOps::AiBi< NF >( q0, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q0, -q0Dotq2 );
    real64 const q1Dotq2 = LvArray::tensorOps::AiBi< NF >( q1, q2 );
    LvArray::tensorOps::scaledAdd< NF >( q2, q1, -q1Dotq2 );
    LvArray::tensorOps::scale< NF >( q2, 1.0 / LvArray::tensorOps::l2Norm< NF >( q2 ) );

    for( localIndex i = 0; i < NF; ++i )
    {
      cellToFaceMat[i][0] = q0[i];
      cellToFaceMat[i][1] = q1[i];
      cellToFaceMat[i][2] = q2[i];
    }
  }
};

/******************************** TPFA Kernel ********************************/

/**
 * @struct TPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the TPFA scheme.
 */
struct TPFACellInnerProductKernel
{
  /**
   * @brief In a given element, recompute the transmissibility matrix in a cell using TPFA.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE static void
  Compute(
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
    ArrayOfArraysView< localIndex const > const & faceToNodes,
    arraySlice1d< localIndex const > const & elemToFaces,
    arraySlice1d< real64 const > const & elemCenter,
    real64 const ( &elemPerm )[3],
    real64 const & lengthTolerance,
    arraySlice2d< real64 > const & transMatrix )
  {
    real64 const areaTolerance = lengthTolerance * lengthTolerance;
    real64 const weightTolerance = 1e-30 * lengthTolerance;

    real64 permTensor[3][3] = { { 0 } };

    // we are ready to compute the transmissibility matrix
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      for( localIndex jfaceLoc = 0; jfaceLoc < NF; ++jfaceLoc )
      {
        // for now, TPFA trans
        if( ifaceLoc == jfaceLoc )
        {
          real64 faceCenter[3], faceNormal[3], faceConormal[3], cellToFaceVec[3];

          // 1) compute the face geometry data: center, normal, vector from cell center to face center
          real64 const faceArea = computationalGeometry::Centroid_3DPolygon(
            faceToNodes[elemToFaces[ifaceLoc]],
            nodePosition,
            faceCenter,
            faceNormal,
            areaTolerance );

          LvArray::tensorOps::copy< 3 >( cellToFaceVec, faceCenter );
          LvArray::tensorOps::subtract< 3 >( cellToFaceVec, elemCenter );

          if( LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceNormal ) < 0.0 )
          {
            LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
          }

          real64 const c2fDistance =
            LvArray::tensorOps::normalize< 3 >( cellToFaceVec );

          // 2) assemble full coefficient tensor from principal axis/components
          HybridFVMInnerProductHelper::MakeFullTensor( elemPerm, permTensor );
          LvArray::tensorOps::hadamardProduct< 3 >( faceConormal,
                                                    elemPerm,
                                                    faceNormal );

          // 3) compute the one-sided face transmissibility
          transMatrix[ifaceLoc][jfaceLoc] =
            LvArray::tensorOps::AiBi< 3 >( cellToFaceVec, faceConormal );
          transMatrix[ifaceLoc][jfaceLoc] *= faceArea / c2fDistance;
          transMatrix[ifaceLoc][jfaceLoc] =
            LvArray::math::max( transMatrix[ifaceLoc][jfaceLoc], weightTolerance );
        }
        else
        {
          transMatrix[ifaceLoc][jfaceLoc] = 0;
        }
      }
    }
  }
};

/******************************** Quasi TPFA Kernel ********************************/

/**
 * @struct QTPFACellInnerProductKernel
 * @brief Struct handling cell inner product in the quasi TPFA scheme.
 */
struct QTPFACellInnerProductKernel
{
  /**
   * @brief In a given element, recompute the transmissibility matrix using a consistent inner product.
   * @param[in] nodePosition the position of the nodes
   * @param[in] faceToNodes the map from the face to their nodes
   * @param[in] elemToFaces the maps from the one-sided face to the corresponding face
   * @param[in] elemCenter the center of the element
   * @param[in] elemVolume the volume of the element
   * @param[in] elemPerm the permeability in the element
   * @param[in] tParam parameter used in the transmissibility matrix computations
   * @param[in] lengthTolerance the tolerance used in the trans calculations
   * @param[inout] transMatrix the output
   *
   * When tParam = 2, we obtain a scheme that reduces to TPFA
   * on orthogonal meshes, but remains consistent on non-orthogonal meshes
   */
  template< localIndex NF >
  GEOSX_HOST_DEVICE static void
  Compute(
    arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & nodePosition,
    ArrayOfArraysView< localIndex const > const & faceToNodes,
    arraySlice1d< localIndex const > const & elemToFaces,
    arraySlice1d< real64 const > const & elemCenter,
    real64 const & elemVolume,
    real64 const ( &elemPerm )[3],
    real64 const & tParam,
    real64 const & lengthTolerance,
    arraySlice2d< real64 > const & transMatrix )
  {
    real64 const areaTolerance = lengthTolerance * lengthTolerance;

    real64 cellToFaceMat[NF][3] = { { 0 } };
    real64 normalsMat[NF][3] = { { 0 } };
    real64 permMat[3][3] = { { 0 } };

    real64 work_dimByNumFaces[3][NF] = { { 0 } };
    real64 worka_numFacesByNumFaces[NF][NF] = { { 0 } };
    real64 workb_numFacesByNumFaces[NF][NF] = { { 0 } };
    real64 workc_numFacesByNumFaces[NF][NF] = { { 0 } };

    real64 q0[NF], q1[NF], q2[NF];

    // 1) fill the matrices cellToFaceMat and normalsMat row by row
    for( localIndex ifaceLoc = 0; ifaceLoc < NF; ++ifaceLoc )
    {
      real64 faceCenter[3], faceNormal[3];

      // compute the face geometry data: center, normal, vector from cell center to face center
      real64 const faceArea = computationalGeometry::Centroid_3DPolygon(
        faceToNodes[elemToFaces[ifaceLoc]],
        nodePosition,
        faceCenter,
        faceNormal,
        areaTolerance );

      q0[ifaceLoc] = faceCenter[0] - elemCenter[0];
      q1[ifaceLoc] = faceCenter[1] - elemCenter[1];
      q2[ifaceLoc] = faceCenter[2] - elemCenter[2];

      real64 const dotProduct = q0[ifaceLoc] * faceNormal[0] +
        q1[ifaceLoc] * faceNormal[1] + q2[ifaceLoc] * faceNormal[2];
      if( dotProduct < 0.0 )
      {
        LvArray::tensorOps::scale< 3 >( faceNormal, -faceArea );
      }
      else
      {
        LvArray::tensorOps::scale< 3 >( faceNormal, faceArea );
      }

      normalsMat[ifaceLoc][0] = faceNormal[0];
      normalsMat[ifaceLoc][1] = faceNormal[1];
      normalsMat[ifaceLoc][2] = faceNormal[2];
    }

    // 2) assemble full coefficient tensor from principal axis/components
    HybridFVMInnerProductHelper::MakeFullTensor( elemPerm, permMat );

    // 3) compute N K N'
    LvArray::tensorOps::AikBjk< 3, NF, 3 >( work_dimByNumFaces, permMat, normalsMat );
    LvArray::tensorOps::AikBkj< NF, NF, 3 >( transMatrix,
                                             normalsMat,
                                             work_dimByNumFaces );

    // 4) compute the orthonormalization of the matrix cellToFaceVec
    HybridFVMInnerProductHelper::Orthonormalize< NF >( q0, q1, q2, cellToFaceMat );

    // 5) compute P_Q = I - QQ'
    // note: we compute -P_Q and then at 6) ( - P_Q ) D ( - P_Q )
    LvArray::tensorOps::addIdentity< NF >( worka_numFacesByNumFaces, -1 );
    LvArray::tensorOps::plusAikAjk< NF, 3 >( worka_numFacesByNumFaces,
                                             cellToFaceMat );

    // 6) compute P_Q D P_Q where D = diag(diag(N K N'))
    // 7) compute T = ( N K N' + t U diag(diag(N K N')) U ) / elemVolume
    real64 const scale = tParam / elemVolume;
    for( localIndex i = 0; i < NF; ++i )
    {
      workb_numFacesByNumFaces[i][i] = scale * transMatrix[i][i];
    }

    LvArray::tensorOps::AikBkj< NF, NF, NF >( workc_numFacesByNumFaces,
                                              workb_numFacesByNumFaces,
                                              worka_numFacesByNumFaces );
    LvArray::tensorOps::scale< NF, NF >( transMatrix, 1 / elemVolume );
    LvArray::tensorOps::plusAikBkj< NF, NF, NF >( transMatrix,
                                                  worka_numFacesByNumFaces,
                                                  workc_numFacesByNumFaces );
  }
};

}  // namespace HybridFVMInnerProduct

}  // namespace geosx

#endif  //GEOSX_PHYSICSSOLVERS_FINITEVOLUME_HYBRIDFVMINNERPRODUCT_HPP
