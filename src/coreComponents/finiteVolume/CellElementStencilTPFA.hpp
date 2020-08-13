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

/**
 * @class CellElementStencilTPFA
 *
 * Provides management of the interior stencil points when using Two-Point flux approximation.
 */
class CellElementStencilTPFA
  : public StencilBase< CellElementStencilTPFA_Traits, CellElementStencilTPFA >,
    public CellElementStencilTPFA_Traits
{
public:
  /**
   * @brief Default constructor.
   */
  CellElementStencilTPFA();

  virtual void
  add( localIndex const numPts,
       localIndex const * const elementRegionIndices,
       localIndex const * const elementSubRegionIndices,
       localIndex const * const elementIndices,
       real64 const * const weights,
       localIndex const connectorIndex ) override final;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex
  size() const override final
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex
  stencilSize( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return MAX_STENCIL_SIZE;
  }
};

} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_CELLELEMENTSTENCILTPFA_HPP_ */
