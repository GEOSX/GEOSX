/*
 * ------------------------------------------------------------------------------------------------------------
 * SPDX-License-Identifier: LGPL-2.1-only
 *
 * Copyright (c) 2018-2019 Lawrence Livermore National Security LLC
 * Copyright (c) 2018-2019 The Board of Trustees of the Leland Stanford Junior University
 * Copyright (c) 2018-2019 Total, S.A
 * Copyright (c) 2019-     GEOSX Contributors
 * All right reserved
 * See top level LICENSE, COPYRIGHT, CONTRIBUTORS, NOTICE, and ACKNOWLEDGEMENTS files for details.
 * ------------------------------------------------------------------------------------------------------------
 */

/**
 * @file VirtualElementBase.hpp
 */

#ifndef GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
#define GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_

#include "mesh/FaceManager.hpp"
#include "mesh/EdgeManager.hpp"
#include "mesh/NodeManager.hpp"

namespace geosx
{
namespace virtualElement
{
/**
 * @brief Base class for VEM implementations.
 */
class VirtualElementBase
{
public:

  /// Default constructor.
  VirtualElementBase() = default;

  /// Default destructor.
  virtual ~VirtualElementBase() = default;

  /**
   * @brief Get the shape function projected derivatives at a given quadrature point.
   * @details The output contains the integral mean of derivatives.
   * @tparam LEAF Type of the derived virtual element implementation.
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projected derivatives. Size will be @ref
   * getNumSupportPoints() x 3.
   */
  template< typename LEAF, typename LEAFBASISDATA >
  GEOSX_HOST_DEVICE
  static real64 getGradN( localIndex const q,
                          LEAFBASISDATA const & basisData,
                          real64 ( & gradN )[LEAF::maxSupportPoints][3] )
  {
    return LEAF::calcGradN( q, basisData, gradN );
  }

  /**
   * @brief Get a value of the stabilization matrix.
   * @tparam LEAF Type of the derived virtual element implementation.
   * @param iBasisFunction The row index.
   * @param jBasisFunction The column index.
   * @return The requested value.
   */
  template< typename LEAF, typename LEAFBASISDATA >
  GEOSX_HOST_DEVICE
  static real64 getStabilizationValue( localIndex const iBasisFunction,
                                       localIndex const jBasisFunction,
                                       LEAFBASISDATA const & basisData
                                       )
  {
    return LEAF::calcStabilizationValues( iBasisFunction, jBasisFunction, basisData );
  }

  /**
   * @brief Get the shape function projections at a given quadrature point.
   * @details The output contains the integral mean of functions
   * @tparam LEAF Type of the derived virtual element implementation.
   * @param q The quadrature point index.
   * @param gradN Return array of the shape function projections. Size will be @ref getNumSupportPoints().
   */
  template< typename LEAF, typename LEAFBASISDATA >
  GEOSX_HOST_DEVICE
  void getN( localIndex const q,
             LEAFBASISDATA const & basisData,
             real64 ( & N )[LEAF::maxSupportPoints] )
  {
    return LEAF::calcN( q, basisData, N );
  }

  /**
   * @brief Get the integration weight for a quadrature point.
   * @tparam LEAF Type of the derived virtual element implementation.
   * @param q Index of the quadrature point.
   * @return The weight.
   */
  template< typename LEAF, typename LEAFBASISDATA >
  GEOSX_HOST_DEVICE
  real64 getTransformedQuadratureWeight( localIndex const q, LEAFBASISDATA const & basisData ) const
  {
    return LEAF::transformedQuadratureWeight( q, basisData );
  }

  /**
   * @brief Getter for the number of quadrature points per element.
   * @tparam LEAF Type of the derived virtual element implementation.
   * @return The number of quadrature points per element.
   */
  template< typename LEAF, typename LEAFBASISDATA >
  localIndex getNumQuadraturePoints( LEAFBASISDATA const & basisData ) const
  {
    return LEAF::getNumQuadraturePoints( basisData );
  }

  /**
   * @brief Getter for the number of support points per element.
   * @tparam LEAF Type of the derived virtual element implementation.
   * @return The number of support points per element.
   */
  template< typename LEAF, typename LEAFBASISDATA >
  localIndex getNumSupportPoints( LEAFBASISDATA const & basisData ) const
  {
    return LEAF::getNumSupportPoints( basisData );
  }
};
}
}

#endif // GEOSX_VIRTUALELEMENT_VIRTUALELEMENTBASE_HPP_
