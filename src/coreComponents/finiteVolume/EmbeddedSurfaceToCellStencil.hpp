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
 * @file EmbeddedSurfaceToCellStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/**
 * @struct EmbeddedSurfaceToCellStencil_Traits
 * Struct to predeclare the types and constexpr values of EmbeddedSurfaceToCellStencil so that they may be used in
 * StencilBase.
 */
struct EmbeddedSurfaceToCellStencil_Traits
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


class EmbeddedSurfaceToCellStencilWrapper : public StencilWrapperBase< EmbeddedSurfaceToCellStencil_Traits >,
  public EmbeddedSurfaceToCellStencil_Traits
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::MaterialViewAccessor< VIEWTYPE >;

  EmbeddedSurfaceToCellStencilWrapper( IndexContainerType & elementRegionIndices,
                                       IndexContainerType & elementSubRegionIndices,
                                       IndexContainerType & elementIndices,
                                       WeightContainerType & weights )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights )
  {}

  /// Default copy constructor
  EmbeddedSurfaceToCellStencilWrapper( EmbeddedSurfaceToCellStencilWrapper const & ) = default;

  /// Default move constructor
  EmbeddedSurfaceToCellStencilWrapper( EmbeddedSurfaceToCellStencilWrapper && ) = default;

  /// Deleted copy assignment operator
  EmbeddedSurfaceToCellStencilWrapper & operator=( EmbeddedSurfaceToCellStencilWrapper const & ) = delete;

  /// Deleted move assignment operator
  EmbeddedSurfaceToCellStencilWrapper & operator=( EmbeddedSurfaceToCellStencilWrapper && ) = delete;

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

};

/**
 * @class EmbeddedSurfaceToCellStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class EmbeddedSurfaceToCellStencil : public StencilBase< EmbeddedSurfaceToCellStencil_Traits, EmbeddedSurfaceToCellStencil >,
  public EmbeddedSurfaceToCellStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  EmbeddedSurfaceToCellStencil();

  virtual void move( LvArray::MemorySpace const space ) override final;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  /**
   * @brief Add an entry to the stencil.
   * @param[in] numPts The number of points in the stencil entry
   * @param[in] cellCenterToEdgeCenter vectors pointing from the cell center to the edge center
   * @param[in] connectorIndex The index of the connector element that the stencil acts across
   */
  void add( localIndex const numPts,
            R1Tensor const * const cellCenterToEdgeCenter,
            localIndex const connectorIndex );


  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = EmbeddedSurfaceToCellStencilWrapper;

  /**
   * @brief Create an update kernel wrapper.
   * @return the wrapper
   */
  StencilWrapper createStencilWrapper()
  {
    return StencilWrapper( m_elementRegionIndices,
                           m_elementSubRegionIndices,
                           m_elementIndices,
                           m_weights );
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

};

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void EmbeddedSurfaceToCellStencilWrapper::computeTransmissibility( localIndex iconn,
                                                                   PERMTYPE permeability,
                                                                   PERMTYPE dPerm_dPressure,
                                                                   real64 (& transmissibility)[2],
                                                                   real64 (& dTrans_dPressure )[2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

//  localIndex const er1  =  m_elementRegionIndices[iconn][1];
//  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
//  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * LvArray::tensorOps::l2Norm< 3 >( permeability[er0][esr0][ei0][0] );
//  real64 const t1 = m_weights[iconn][1] * permeability[er1][esr1][ei1][0][0];

  real64 const harmonicWeight   = t0; // *t1 / (t0+t1);

  real64 const value =  harmonicWeight;

  transmissibility[0] = value;
  transmissibility[1] = -value;

  dTrans_dPressure[0] = 0.0 * dPerm_dPressure[er0][esr0][ei0][0][0];
  dTrans_dPressure[1] = 0.0;
}

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void EmbeddedSurfaceToCellStencilWrapper::computeTransmissibility( localIndex iconn,
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

  //  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  //  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  //  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * LvArray::tensorOps::l2Norm< 3 >( permeability[er0][esr0][ei0][0] );
  //  real64 const t1 = m_weights[iconn][1] * permeability[er1][esr1][ei1][0][0];

  real64 const harmonicWeight   = t0; // *t1 / (t0+t1);

  real64 const value =  harmonicWeight;

  transmissibility[0] = value;
  transmissibility[1] = -value;

  dTrans_dPressure[0] = 0.0 * dPerm_dPressure[er0][esr0][ei0][0][0];
  dTrans_dPressure[1] = 0.0;

  dTrans_dAperture[0] = 0.0 * dPerm_dAperture[er0][esr0][ei0][0][0];
  dTrans_dAperture[1] = 0.0;
}
} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_EMBEDDEDSURFACETOCELLSTENCIL_HPP_ */
