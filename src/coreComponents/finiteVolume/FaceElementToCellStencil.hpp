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
 * @file FaceElementToCellStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/**
 * @struct FaceElementToCellStencil_Traits
 * Struct to predeclare the types and constexpr values of FaceElementToCellStencil so that they may be used in
 * StencilBase.
 */
struct FaceElementToCellStencil_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = array2d< localIndex >;

  /// The array view type for the stencil indices
  using IndexContainerViewType = arrayView2d< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = arrayView2d< localIndex const >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = array2d< real64 >;

  /// The array view type for the stencil weights
  using WeightContainerViewType = arrayView2d< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = arrayView2d< real64 const >;

  /// Number of points the flux is between (always 2 for TPFA)
  static constexpr localIndex NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil (this is 2 for TPFA)
  static constexpr localIndex MAX_STENCIL_SIZE = 2;
};


class FaceElementToCellStencilWrapper : public StencilWrapperBase< FaceElementToCellStencil_Traits >,
  public FaceElementToCellStencil_Traits
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::MaterialViewAccessor< VIEWTYPE >;

  FaceElementToCellStencilWrapper( IndexContainerType & elementRegionIndices,
                                   IndexContainerType & elementSubRegionIndices,
                                   IndexContainerType & elementIndices,
                                   WeightContainerType & weights,
                                   arrayView2d< real64 > const & faceNormal,
                                   arrayView2d< real64 > const & cellToFaceVec,
                                   arrayView1d< real64 > const & transMultiplier )
    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_faceNormal( faceNormal ),
    m_cellToFaceVec( cellToFaceVec ),
    m_transMultiplier( transMultiplier )
  {}

  /// Default copy constructor
  FaceElementToCellStencilWrapper( FaceElementToCellStencilWrapper const & ) = default;

  /// Default move constructor
  FaceElementToCellStencilWrapper( FaceElementToCellStencilWrapper && ) = default;

  /// Deleted copy assignment operator
  FaceElementToCellStencilWrapper & operator=( FaceElementToCellStencilWrapper const & ) = delete;

  /// Deleted move assignment operator
  FaceElementToCellStencilWrapper & operator=( FaceElementToCellStencilWrapper && ) = delete;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return MAX_STENCIL_SIZE;
  }

  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex numPointsInFlux( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return NUM_POINT_IN_FLUX;
  }

  template< typename PERMTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeTransmissibility( localIndex iconn,
                                PERMTYPE permeability,
                                PERMTYPE dPerm_dPressure,
                                real64 ( &transmissibility )[2],
                                real64 ( &dTrans_dPressure )[2] ) const;

  template< typename PERMTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeTransmissibility( localIndex iconn,
                                PERMTYPE permeability,
                                PERMTYPE dPerm_dPressure,
                                PERMTYPE dPerm_dAperture,
                                real64 ( &transmissibility )[2],
                                real64 ( &dTrans_dPressure )[2],
                                real64 ( &dTrans_dAperture )[2] ) const;

private:

  arrayView2d< real64 > m_faceNormal;
  arrayView2d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
};

/**
 * @class FaceElementToCellStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class FaceElementToCellStencil : public StencilBase< FaceElementToCellStencil_Traits, FaceElementToCellStencil >,
  public FaceElementToCellStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  FaceElementToCellStencil();

  virtual void move( LvArray::MemorySpace const space ) override final;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;



  void addVectors( real64 const & transMultiplier,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[3] );

  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = FaceElementToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  StencilWrapper createStencilWrapper()
  {
    return StencilWrapper( m_elementRegionIndices,
                           m_elementSubRegionIndices,
                           m_elementIndices,
                           m_weights,
                           m_faceNormal,
                           m_cellToFaceVec,
                           m_transMultiplier );
  }


  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex stencilSize( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return MAX_STENCIL_SIZE;
  }

private:
  array2d< real64 > m_faceNormal;
  array2d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;

};

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FaceElementToCellStencilWrapper::computeTransmissibility( localIndex iconn,
                                                               PERMTYPE permeability,
                                                               PERMTYPE dPerm_dPressure,
                                                               real64 (& transmissibility)[2],
                                                               real64 (& dTrans_dPressure )[2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  real64 halfTrans = m_weights[iconn][0];

  real64 faceConormal[3];

  LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, permeability[er0][esr0][ei0][0], m_faceNormal[iconn] );
  halfTrans *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal ) * m_transMultiplier[iconn];

  transmissibility[0] = halfTrans;
  transmissibility[1] = -halfTrans;

  dTrans_dPressure[0] = 0.0 * dPerm_dPressure[er0][esr0][ei0][0][0];
  dTrans_dPressure[1] = 0.0;
}

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void FaceElementToCellStencilWrapper::computeTransmissibility( localIndex iconn,
                                                               PERMTYPE permeability,
                                                               PERMTYPE dPerm_dPressure,
                                                               PERMTYPE dPerm_dAperture,
                                                               real64 (& transmissibility)[2],
                                                               real64 (& dTrans_dPressure )[2],
                                                               real64 (& dTrans_dAperture )[2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  real64 halfTrans = m_weights[iconn][0];

  real64 faceConormal[3];

  LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, permeability[er0][esr0][ei0][0], m_faceNormal[iconn] );
  halfTrans *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn], faceConormal ) * m_transMultiplier[iconn];

  transmissibility[0] = halfTrans;
  transmissibility[1] = -halfTrans;

  dTrans_dPressure[0] = 0.0 * dPerm_dPressure[er0][esr0][ei0][0][0];
  dTrans_dPressure[1] = 0.0;

  dTrans_dAperture[0] = 0.0 * dPerm_dAperture[er0][esr0][ei0][0][0];
  dTrans_dAperture[1] = 0.0;
}


} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_FACEELEMENTTOCELLSTENCIL_HPP_ */
