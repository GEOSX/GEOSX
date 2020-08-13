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
 * @file BoundaryStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_BOUNDARYSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_BOUNDARYSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/**
 * @struct BoundaryStencil_Traits
 * Struct to predeclare the types and consexpr values of BoundaryStencil so that they may be used in StencilBase.
 */
struct BoundaryStencil_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = array2d< localIndex >;

  /// The array view type for the stencil indices
  using IndexContainerViewType = IndexContainerType::ViewType;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = IndexContainerType::ViewTypeConst;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = array2d< real64 >;

  /// The array view type for the stencil weights
  using WeightContainerViewType = WeightContainerType::ViewType;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = WeightContainerType::ViewTypeConst;

  /// Number of points the flux is between (always 2 for TPFA)
  static constexpr localIndex NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil (this is 2 for TPFA)
  static constexpr localIndex MAX_STENCIL_SIZE = 2;
};

/**
 * @class BoundaryStencil
 *
 * Provides management of the boundary stencil points
 * (stencils used to prescribe boundary conditions on domain boundaries, i.e. faces)
 */
class BoundaryStencil
  : public StencilBase< BoundaryStencil_Traits, BoundaryStencil >,
    public BoundaryStencil_Traits
{
public:
  /**
   * @brief Defines the order of element/face in the stencil.
   */
  struct Order
  {
    static constexpr localIndex ELEM = 0;  ///< Order of element index in stencil
    static constexpr localIndex FACE = 1;  ///< Order of face index in stencil
  };

  /// default constructor
  BoundaryStencil();

  virtual void
  add( localIndex const numPts,
       localIndex const * const elementRegionIndices,
       localIndex const * const elementSubRegionIndices,
       localIndex const * const elementIndices,
       real64 const * const weights,
       localIndex const connectorIndex ) override final;

  virtual localIndex
  size() const override final
  {
    return m_elementRegionIndices.size( 0 );
  }

  /**
   * @brief Gives the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  constexpr localIndex
  stencilSize( localIndex const index ) const
  {
    GEOSX_UNUSED_VAR( index )
    return MAX_STENCIL_SIZE;
  }
};

}  // namespace geosx

#endif  //GEOSX_FINITEVOLUME_BOUNDARYSTENCIL_HPP_
