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
 * @file FaceElementStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_

#include "StencilBase.hpp"

namespace geosx
{
/// @cond DO_NOT_DOCUMENT
// TODO remove! This option allows for the creation of new mass inside a newly
// created FaceElement. The new mass will be equal to:
// creationMass = defaultDensity * defaultAperture * faceArea.
// If 0, then the beginning of step density is artificially set to zero...which
// may cause some newton convergence problems.
#define ALLOW_CREATION_MASS 1

// TODO remove! This option sets the pressure in a newly created FaceElement to
// be the lowest value of all attached non-new FaceElements.
#define SET_CREATION_PRESSURE 1

// TODO remove! This option sets the nodal displacements attached a newly
// created FaceElement to some scalar fraction of the aperture of the
// lowest attached non-new FaceElements.
#define SET_CREATION_DISPLACEMENT 0
/// @endcond

/**
 * @struct FaceElementStencil_Traits
 * Struct to predeclare the types and constexpr values of FaceElementStencil so that they may be used in
 * StencilBase.
 */
struct FaceElementStencil_Traits
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

  /// Number of points the flux is between (normally 2)
  static localIndex constexpr NUM_POINT_IN_FLUX = 6;

  /// Maximum number of points in a stencil
  static localIndex constexpr MAX_STENCIL_SIZE = 6;
};

/**
 * @class FaceElementStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class FaceElementStencil
  : public StencilBase< FaceElementStencil_Traits, FaceElementStencil >,
    public FaceElementStencil_Traits
{
public:
  /**
   * @brief Default constructor.
   */
  FaceElementStencil();

  virtual void
  move( LvArray::MemorySpace const space ) override final;

  virtual void
  add( localIndex const numPts,
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
  void
  add( localIndex const numPts,
       R1Tensor const * const cellCenterToEdgeCenter,
       localIndex const connectorIndex );

  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex
  size() const override final
  {
    return m_elementRegionIndices.size();
  }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  localIndex
  stencilSize( localIndex index ) const
  {
    return m_elementRegionIndices.sizeOfArray( index );
  }

  /**
   * @brief Give the array of vectors pointing from the cell center to the edge center.
   * @return The array of vectors pointing from the cell center to the edge center
   */
  ArrayOfArraysView< R1Tensor const > const &
  getCellCenterToEdgeCenters() const
  {
    return m_cellCenterToEdgeCenters.toViewConst();
  }

private:
  ArrayOfArrays< R1Tensor > m_cellCenterToEdgeCenters;
};

} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_FACEELEMENTSTENCIL_HPP_ */
