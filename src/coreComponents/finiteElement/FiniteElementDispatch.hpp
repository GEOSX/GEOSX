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
 * @file FiniteElementDispatch.hpp
 */

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_

#include "elementFormulations/H1_Hexahedron_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Pyramid_Lagrange1_Gauss5.hpp"
#include "elementFormulations/H1_QuadrilateralFace_Lagrange1_GaussLegendre2.hpp"
#include "elementFormulations/H1_Tetrahedron_Lagrange1_Gauss1.hpp"
#include "elementFormulations/H1_TriangleFace_Lagrange1_Gauss1.hpp"
#include "elementFormulations/H1_Wedge_Lagrange1_Gauss6.hpp"

namespace geosx
{
namespace finiteElement
{
struct ParentElementTypeStrings
{
  static constexpr auto Tetrahedon = "C3D4";
  static constexpr auto Pyramid = "C3D5";
  static constexpr auto Prism = "C3D6";
  static constexpr auto Hexahedron = "C3D8";
  static constexpr auto Quadralateral = "C2D4";
  static constexpr auto Triangle = "C2D3";
  static constexpr auto Polyhedral = "POLYHEDRAL";
  static constexpr auto Polytope = "POLYTOPE";
};

template< typename LAMBDA >
void
dispatch3D( FiniteElementBase const & input, LAMBDA && lambda )
{
  if( auto const * const ptr1 =
        dynamic_cast< H1_Hexahedron_Lagrange1_GaussLegendre2 const * >( &input ) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 =
             dynamic_cast< H1_Wedge_Lagrange1_Gauss6 const * >( &input ) )
  {
    lambda( *ptr2 );
  }
  else if( auto const * const ptr3 =
             dynamic_cast< H1_Tetrahedron_Lagrange1_Gauss1 const * >( &input ) )
  {
    lambda( *ptr3 );
  }
  else if( auto const * const ptr4 =
             dynamic_cast< H1_Pyramid_Lagrange1_Gauss5 const * >( &input ) )
  {
    lambda( *ptr4 );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch3D() is not implemented for input of "
                 << typeid( input ).name() );
  }
}

template< typename LAMBDA >
void
dispatch2D( FiniteElementBase const & input, LAMBDA && lambda )
{
  if( auto const * const ptr1 =
        dynamic_cast< H1_QuadrilateralFace_Lagrange1_GaussLegendre2 const * >( &input ) )
  {
    lambda( *ptr1 );
  }
  else if( auto const * const ptr2 =
             dynamic_cast< H1_TriangleFace_Lagrange1_Gauss1 const * >( &input ) )
  {
    lambda( *ptr2 );
  }
  else
  {
    GEOSX_ERROR( "finiteElement::dispatch2D() is not implemented for input of: "
                 << typeid( input ).name() );
  }
}

}  // namespace finiteElement
}  // namespace geosx

#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISPATCH_HPP_ */
