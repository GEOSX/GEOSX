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
 * @file FiniteElementSpace.hpp
 */

#ifndef GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_
#define GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_

#include "common/TimingMacros.hpp"
#include "dataRepository/Group.hpp"
#include "dataRepository/Wrapper.hpp"
#include "LvArray/src/tensorOps.hpp"
#include "FiniteElementDispatch.hpp"


namespace geosx
{

// TODO remove when these quantities are placed inside the FiniteElementBase
// class.
namespace dataRepository
{
namespace keys
{
string const dNdX = "dNdX";
string const detJ = "detJ";
}
}



class FiniteElementDiscretization : public dataRepository::Group
{
public:



  FiniteElementDiscretization() = delete;

  explicit FiniteElementDiscretization( string const & name, Group * const parent );

  ~FiniteElementDiscretization() override;

  /**
   * @name Static Factory Catalog Functions
   */
  ///@{
  static string catalogName() { return "FiniteElementSpace"; }

  ///@}

  template< typename SUBREGION_TYPE,
            typename FE_TYPE >
  void calculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                        SUBREGION_TYPE * const elementSubRegion,
                                        FE_TYPE & fe ) const;


  /**
   * @brief Factory method to instantiate a type of finite element formulation.
   * @param parentElementShape String key that indicates the type of
   *   element/basis/formulation that should be instantiated.
   * @return A unique_ptr< FinteElementBase > which contains the new
   *   instantiation.
   */
  std::unique_ptr< finiteElement::FiniteElementBase >
  factory( string const & parentElementShape ) const;

private:

  struct viewKeyStruct
  {
    static constexpr char const * orderString() { return "order"; }
    static constexpr char const * formulationString() { return "formulation"; }
  };

  /// The order of the finite element basis
  int m_order;

  /// Optional string indicating any specialized formulation type.
  string m_formulation;

  void postProcessInput() override final;

};

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
void
FiniteElementDiscretization::
  calculateShapeFunctionGradients( arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const & X,
                                   SUBREGION_TYPE * const elementSubRegion,
                                   FE_TYPE & finiteElement ) const
{
  GEOSX_MARK_FUNCTION;

  array4d< real64 > & dNdX = elementSubRegion->dNdX();
  array2d< real64 > & detJ = elementSubRegion->detJ();
  auto const & elemsToNodes = elementSubRegion->nodeList().toViewConst();

  constexpr localIndex numNodesPerElem = FE_TYPE::numNodes;
  constexpr localIndex numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;
  dNdX.resizeWithoutInitializationOrDestruction( elementSubRegion->size(), numQuadraturePointsPerElem, numNodesPerElem, 3 );
  detJ.resize( elementSubRegion->size(), numQuadraturePointsPerElem );

  finiteElement.setGradNView( dNdX.toViewConst() );
  finiteElement.setDetJView( detJ.toViewConst() );

  for( localIndex k = 0; k < elementSubRegion->size(); ++k )
  {
    real64 xLocal[numNodesPerElem][3];
    for( localIndex a=0; a< numNodesPerElem; ++a )
    {
      localIndex const nodeIndex = elemsToNodes[ k][ a ];
      for( int i=0; i<3; ++i )
      {
        xLocal[ a ][ i ] = X[ nodeIndex ][ i ];
      }
    }



    for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
    {
      real64 dNdXLocal[numNodesPerElem][3];
      detJ( k, q ) = finiteElement.calcGradN( q, xLocal, dNdXLocal );

      for( localIndex b = 0; b < numNodesPerElem; ++b )
      {
        LvArray::tensorOps::copy< 3 >( dNdX[ k ][ q ][ b ], dNdXLocal[b] );
      }
    }
  }

}


} /* namespace geosx */

#endif /* GEOSX_FINITEELEMENT_FINITEELEMENTDISCRETIZATION_HPP_ */
