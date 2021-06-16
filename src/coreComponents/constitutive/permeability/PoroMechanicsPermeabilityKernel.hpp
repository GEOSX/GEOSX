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
 * @file ConstantPermeability.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_PERMEABILITY_POROMECHANICSPERMEABILITYKERNEL_HPP_
#define GEOSX_CONSTITUTIVE_PERMEABILITY_POROMECHANICSPERMEABILITYKERNEL_HPP_

#include "common/DataTypes.hpp"
#include "common/GEOS_RAJA_Interface.hpp"
#include "constitutive/permeability/PermeabilityBase.hpp"


/******************************** PermeabilityKernel ********************************/

namespace geosx
{

namespace PermeabilityKernels
{

template< typename SUBREGION_TYPE,
          typename FE_TYPE >
class PoroMechanicsPermeabilityKernel
{
public:

  static constexpr int numNodesPerElem = FE_TYPE::numNodes;

  static constexpr int numQuadraturePointsPerElem = FE_TYPE::numQuadraturePoints;

  PoroMechanicsPermeabilityKernel( SUBREGION_TYPE const & elementSubRegion,
                                   FE_TYPE const & finiteElementSpace,
                                   NodeManager const & nodeManager ):
    m_elemsToNodes( elementSubRegion.nodeList().toViewConst() ),
    m_finiteElementSpace( finiteElementSpace ),
    m_nodeLocations( nodeManager.referencePosition()),
    m_disp( nodeManager.totalDisplacement() )
  {}

  template< typename POLICY,
            typename KERNEL_TYPE,
            typename PERM_WRAPPER >
  static
  void
  launch( localIndex const size,
          KERNEL_TYPE const & permKernel,
          PERM_WRAPPER const & permWrapper,
          arrayView1d< real64 const > const & pressure,
          arrayView1d< real64 const > const & deltaPressure,
          arrayView2d< real64 const > const & porosity,
          arrayView2d< real64 const > const & dPorosity_dVolStrain,
          arrayView3d< real64 > const & dPerm_dDisplacement )
  {

    forAll< POLICY >( size, [=] GEOSX_HOST_DEVICE ( localIndex const k )
    {
      permKernel.update( k,
                         permWrapper,
                         pressure,
                         deltaPressure,
                         porosity,
                         dPorosity_dVolStrain,
                         dPerm_dDisplacement );

    } );
  }

  template< typename PERM_WRAPPER >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void update( localIndex const k,
               PERM_WRAPPER const & permWrapper,
               arrayView1d< real64 const > const & pressure,
               arrayView1d< real64 const > const & deltaPressure,
               arrayView2d< real64 const > const & porosity,
               arrayView2d< real64 const > const & GEOSX_UNUSED_PARAM ( dPorosity_dVolStrain ),
               arrayView3d< real64 > const & GEOSX_UNUSED_PARAM ( dPerm_dDisplacement ) ) const
  {
    real64 displacementLocal[numNodesPerElem][3];
    real64 xLocal[numNodesPerElem][3];

    for( localIndex a=0; a<numNodesPerElem; ++a )
    {
      localIndex const localNodeIndex = m_elemsToNodes( k, a );

      for( int i=0; i<3; ++i )
      {
        xLocal[a][i] = m_nodeLocations[localNodeIndex][i];
        displacementLocal[a][i] = m_disp[localNodeIndex][i];
      }
    }

    for( localIndex q = 0; q < numQuadraturePointsPerElem; ++q )
    {
      real64 N[numNodesPerElem];
      real64 dNdX[numNodesPerElem][3];
      FE_TYPE::calcN( q, N );
      m_finiteElementSpace.template getGradN< FE_TYPE >( k, q, xLocal, dNdX );

      real64 strainIncrement[6] = {0};
      real64 dPerm_dVolStrain[3] = {0};

      FE_TYPE::symmetricGradient( dNdX, displacementLocal, strainIncrement );

      real64 const volStrain = LvArray::tensorOps::symTrace< 3 >( strainIncrement );

      permWrapper.updatePorosity( k, q, porosity( k, q ) );

      permWrapper.updatePressureStrain( k, q, pressure[k] + deltaPressure[k], volStrain, dPerm_dVolStrain );

    }
    // TODO: chain rule all dependencies to get to dPerm_dDisplacement
  }

protected:
  /// The element to nodes map.
  traits::ViewTypeConst< typename SUBREGION_TYPE::NodeMapType::base_type > const m_elemsToNodes;

  /// The finite element space/discretization object for the element type in
  /// the SUBREGION_TYPE.
  FE_TYPE const & m_finiteElementSpace;

  /// The array containing the nodal position array.
  arrayView2d< real64 const, nodes::REFERENCE_POSITION_USD > const m_nodeLocations;

  /// The rank-global displacement array.
  arrayView2d< real64 const, nodes::TOTAL_DISPLACEMENT_USD > const m_disp;
};

}/* namespace Peremability */

} /* namespace geosx */


#endif // GEOSX_CONSTITUTIVE_PERMEABILITY_POROMECHANICSPERMEABILITYKERNEL_HPP_
