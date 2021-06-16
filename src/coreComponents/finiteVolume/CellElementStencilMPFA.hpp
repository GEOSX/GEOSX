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
 * @file CellElementStencilMPFA.hpp
 */

#ifndef GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_
#define GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_

#include "StencilBase.hpp"

namespace geosx
{


/**
 * @struct CellElementStencilMPFA_Traits
 * Struct to predeclare the types and constexpr values of CellElementStencilMPFA so that they may be used in
 * StencilBase.
 */
struct CellElementStencilMPFA_Traits
{
  /// The array type that will be used to store the indices of the stencil contributors
  using IndexContainerType = ArrayOfArrays< localIndex >;

  /// The array view type for the stencil indices
  using IndexContainerViewType = ArrayOfArraysView< localIndex >;

  /// The array view to const type for the stencil indices
  using IndexContainerViewConstType = ArrayOfArraysView< localIndex const >;

  /// The array type that is used to store the weights of the stencil contributors
  using WeightContainerType = ArrayOfArrays< real64 >;

  /// The array view type for the stencil weights
  using WeightContainerViewType = ArrayOfArraysView< real64 >;

  /// The array view to const type for the stencil weights
  using WeightContainerViewConstType = ArrayOfArraysView< real64 const >;

  ///
  using CoefficientAccessorViewConstType = ArrayOfArraysView< real64 const >;

  /// Number of points the flux is between (always 2)
  static localIndex constexpr NUM_POINT_IN_FLUX = 2;

  /// Maximum number of points in a stencil
  static localIndex constexpr MAX_STENCIL_SIZE = 18;

};

/**
 * @class CellElementStencilMPFA
 *
 * Provides management of the interior stencil points when using a Multi-point flux approximation.
 */
class CellElementStencilMPFA : public StencilBase< CellElementStencilMPFA_Traits, CellElementStencilMPFA >,
  public CellElementStencilMPFA_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  CellElementStencilMPFA();

  virtual void reserve( localIndex const size ) override final;

  virtual void add( localIndex const numPts,
                    localIndex const * const elementRegionIndices,
                    localIndex const * const elementSubRegionIndices,
                    localIndex const * const elementIndices,
                    real64 const * const weights,
                    localIndex const connectorIndex ) override final;

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of points in a stencil entry.
   * @param[in] index of the stencil entry for which to query the size
   * @return the size of a stencil entry
   */
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }


  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  constexpr localIndex numPointsInFlux( localIndex index ) const
  {
    GEOSX_UNUSED_VAR( index );
    return NUM_POINT_IN_FLUX;
  }

};

} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_CELLELEMENTSTENCILMPFA_HPP_ */
