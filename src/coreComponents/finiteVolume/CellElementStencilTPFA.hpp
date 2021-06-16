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
 * @file CellElementStencilTPFA.hpp
 */

#ifndef GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_
#define GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_

#include "StencilBase.hpp"

namespace geosx
{

/**
 * @struct CellElementStencilTPFA_Traits
 * Struct to predeclare the types and constexpr values of CellElementStencilTPFA so that they may be used in
 * StencilBase.
 */
struct CellElementStencilTPFA_Traits
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

class CellElementStencilTPFAWrapper : public StencilWrapperBase< CellElementStencilTPFA_Traits >,
  public CellElementStencilTPFA_Traits
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using PermeabilityViewAccessor = ElementRegionManager::MaterialViewAccessor< VIEWTYPE >;

  CellElementStencilTPFAWrapper( IndexContainerType & elementRegionIndices,
                                 IndexContainerType & elementSubRegionIndices,
                                 IndexContainerType & elementIndices,
                                 WeightContainerType & weights,
                                 arrayView2d< real64 > const & faceNormal,
                                 arrayView3d< real64 > const & cellToFaceVec,
                                 arrayView1d< real64 > const & transMultiplier )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_faceNormal( faceNormal ),
    m_cellToFaceVec( cellToFaceVec ),
    m_transMultiplier( transMultiplier )
  {}

  /// Default copy constructor
  CellElementStencilTPFAWrapper( CellElementStencilTPFAWrapper const & ) = default;

  /// Default move constructor
  CellElementStencilTPFAWrapper( CellElementStencilTPFAWrapper && ) = default;

  /// Deleted copy assignment operator
  CellElementStencilTPFAWrapper & operator=( CellElementStencilTPFAWrapper const & ) = delete;

  /// Deleted move assignment operator
  CellElementStencilTPFAWrapper & operator=( CellElementStencilTPFAWrapper && ) = delete;


  template< typename PERMTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeTransmissibility( localIndex iconn,
                                PERMTYPE permeability,
                                PERMTYPE dPerm_dPressure,
                                real64 ( &transmissibility )[2],
                                real64 ( &dTrans_dPressure )[2] ) const;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size( 0 ); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
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


private:

  arrayView2d< real64 > m_faceNormal;
  arrayView3d< real64 > m_cellToFaceVec;
  arrayView1d< real64 > m_transMultiplier;
};


/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA : public StencilBase< CellElementStencilTPFA_Traits, CellElementStencilTPFA >,
  public CellElementStencilTPFA_Traits
{
public:

  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::MaterialViewAccessor< VIEWTYPE >;

  /**
   * @brief Default constructor.
   */
  CellElementStencilTPFA();

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  void addVectors( real64 const & transMultiplier,
                   real64 const (&faceNormal)[3],
                   real64 const (&cellToFaceVec)[2][3] );

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

  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = CellElementStencilTPFAWrapper;

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

private:
  array2d< real64 > m_faceNormal;
  array3d< real64 > m_cellToFaceVec;
  array1d< real64 > m_transMultiplier;

};

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void CellElementStencilTPFAWrapper::computeTransmissibility( localIndex iconn,
                                                             PERMTYPE permeability,
                                                             PERMTYPE dPerm_dPressure,
                                                             real64 (& transmissibility)[2],
                                                             real64 (& dTrans_dPressure )[2] ) const
{
  GEOSX_UNUSED_VAR( dPerm_dPressure );

  real64 halfTrans[2];

  // real64 const tolerance = 1e-30 * lengthTolerance; // TODO: choice of constant based on physics?

  for( localIndex i =0; i<2; i++ )
  {
    localIndex const er  =  m_elementRegionIndices[iconn][i];
    localIndex const esr =  m_elementSubRegionIndices[iconn][i];
    localIndex const ei  =  m_elementIndices[iconn][i];

    halfTrans[i] = m_weights[iconn][i];

    // Proper computation
    real64 faceNormal[3], faceConormal[3];

    LvArray::tensorOps::copy< 3 >( faceNormal, m_faceNormal[iconn] );

    if( LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceNormal ) < 0.0 )
    {
      LvArray::tensorOps::scale< 3 >( faceNormal, -1 );
    }

    LvArray::tensorOps::hadamardProduct< 3 >( faceConormal, permeability[er][esr][ei][0], faceNormal );
    halfTrans[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );

    // correct negative weight issue arising from non-K-orthogonal grids
    if( halfTrans[i] < 0.0 )
    {
      LvArray::tensorOps::hadamardProduct< 3 >( faceConormal,
                                                permeability[er][esr][ei][0],
                                                m_cellToFaceVec[iconn][i] );
      halfTrans[i] = m_weights[iconn][i];
      halfTrans[i] *= LvArray::tensorOps::AiBi< 3 >( m_cellToFaceVec[iconn][i], faceConormal );
    }
  }

  // Do harmonic and arithmetic averaging
  real64 const product = halfTrans[0]*halfTrans[1];
  real64 const sum = halfTrans[0]+halfTrans[1];

  real64 const harmonicWeight   = sum > 0 ? product / sum : 0.0;
  real64 const arithmeticWeight = sum / 2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member if it is really necessary

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;
  for( localIndex ke = 0; ke < 2; ++ke )
  {
    transmissibility[ke] = m_transMultiplier[iconn] * value * (ke == 0 ? 1 : -1);
  }

  dTrans_dPressure[0] = 0.0;
  dTrans_dPressure[1] = 0.0;
}


} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
