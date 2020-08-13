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
 * @file TwoPointFluxApproximation.hpp
 */

#ifndef GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
#define GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_

#include "finiteVolume/FluxApproximationBase.hpp"

namespace geosx
{
/**
 * @class TwoPointFluxApproximation
 *
 * Provides management of the interior stencil points when using a two-point flux approximation.
 */
class TwoPointFluxApproximation : public FluxApproximationBase
{
public:
  /**
   * @brief Static Factory Catalog Functions.
   * @return the catalog name
   */
  static std::string
  CatalogName()
  {
    return "TwoPointFluxApproximation";
  }

  TwoPointFluxApproximation() = delete;

  /**
   * @brief Constructor.
   * @param name the name of the TwoPointFluxApproximation in the data repository
   * @param parent the parent group of this group.
   */
  TwoPointFluxApproximation( std::string const & name,
                             dataRepository::Group * const parent );

protected:
  virtual void
  registerCellStencil( Group & stencilGroup ) const override;

  virtual void
  computeCellStencil( MeshLevel & mesh ) const override;

  virtual void
  registerFractureStencil( Group & stencilGroup ) const override;

  virtual void
  addToFractureStencil( MeshLevel & mesh,
                        string const & faceElementRegionName,
                        bool const initFlag ) const override;

  virtual void
  registerBoundaryStencil( Group & stencilGroup,
                           string const & setName ) const override;

  virtual void
  computeBoundaryStencil(
    MeshLevel & mesh,
    string const & setName,
    SortedArrayView< localIndex const > const & faceSet ) const override;
};

}  // namespace geosx

#endif  //GEOSX_FINITEVOLUME_TWOPOINTFLUXAPPROXIMATION_HPP_
