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
 * @file SurfaceElementStencil.hpp
 */

#ifndef GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_
#define GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_

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
 * @struct SurfaceElementStencil_Traits
 * Struct to predeclare the types and constexpr values of SurfaceElementStencil so that they may be used in
 * StencilBase.
 */
struct SurfaceElementStencil_Traits
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


class SurfaceElementStencilWrapper : public StencilWrapperBase< SurfaceElementStencil_Traits >,
  public SurfaceElementStencil_Traits
{
public:

  template< typename VIEWTYPE >
  using ElementViewConst = ElementRegionManager::ElementViewConst< VIEWTYPE >;

  template< typename VIEWTYPE >
  using CoefficientAccessor = ElementRegionManager::MaterialViewAccessor< VIEWTYPE >;

  SurfaceElementStencilWrapper( IndexContainerType & elementRegionIndices,
                                IndexContainerType & elementSubRegionIndices,
                                IndexContainerType & elementIndices,
                                WeightContainerType & weights,
                                ArrayOfArrays< R1Tensor > & cellCenterToEdgeCenters )

    : StencilWrapperBase( elementRegionIndices, elementSubRegionIndices, elementIndices, weights ),
    m_cellCenterToEdgeCenters( cellCenterToEdgeCenters.toView() )
  {}

  /// Default copy constructor
  SurfaceElementStencilWrapper( SurfaceElementStencilWrapper const & ) = default;

  /// Default move constructor
  SurfaceElementStencilWrapper( SurfaceElementStencilWrapper && ) = default;

  /// Deleted copy assignment operator
  SurfaceElementStencilWrapper & operator=( SurfaceElementStencilWrapper const & ) = delete;

  /// Deleted move assignment operator
  SurfaceElementStencilWrapper & operator=( SurfaceElementStencilWrapper && ) = delete;

  /**
   * @brief Give the number of stencil entries.
   * @return The number of stencil entries
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }


  /**
   * @brief Give the number of points between which the flux is.
   * @param[in] index of the stencil entry for which to query the size
   * @return the number of points.
   */
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  localIndex numPointsInFlux( localIndex index ) const
  {
    return stencilSize( index );
  }

  template< typename PERMTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeTransmissibility( localIndex iconn,
                                PERMTYPE permeability,
                                PERMTYPE dPerm_dPressure,
                                real64 ( &transmissibility )[2],
                                real64 ( &dTrans_dPressure )[2] ) const;


  template< typename PERMTYPE >
  GEOSX_HOST_DEVICE
  GEOSX_FORCE_INLINE
  void computeTransmissibility( localIndex iconn,
                                PERMTYPE permeability,
                                PERMTYPE dPerm_dPressure,
                                PERMTYPE dPerm_dAperture,
                                real64 ( &transmissibility )[2],
                                real64 ( &dTrans_dPressure )[2],
                                real64 ( &dTrans_dAperture )[2] ) const;

  ArrayOfArraysView< R1Tensor const > getCellCenterToEdgeCenters() const
  { return m_cellCenterToEdgeCenters.toViewConst(); }

private:
  ArrayOfArraysView< R1Tensor > m_cellCenterToEdgeCenters;
};

/**
 * @class SurfaceElementStencil
 *
 * Provides management of the interior stencil points for a face elements when using Two-Point flux approximation.
 */
class SurfaceElementStencil : public StencilBase< SurfaceElementStencil_Traits, SurfaceElementStencil >,
  public SurfaceElementStencil_Traits
{
public:

  /**
   * @brief Default constructor.
   */
  SurfaceElementStencil();

  virtual void move( LvArray::MemorySpace const space ) override final;

  virtual void add( localIndex const numPts,
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
  void add( localIndex const numPts,
            R1Tensor const * const cellCenterToEdgeCenter,
            localIndex const connectorIndex );


  /// Type of kernel wrapper for in-kernel update
  using StencilWrapper = SurfaceElementStencilWrapper;

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
                           m_cellCenterToEdgeCenters );
  }


  /**
   * @brief Return the stencil size.
   * @return the stencil size
   */
  virtual localIndex size() const override final
  { return m_elementRegionIndices.size(); }

  /**
   * @brief Give the number of stencil entries for the provided index.
   * @param[in] index the index of which the stencil size is request
   * @return The number of stencil entries for the provided index
   */
  localIndex stencilSize( localIndex index ) const
  { return m_elementRegionIndices.sizeOfArray( index ); }

  /**
   * @brief Give the array of vectors pointing from the cell center to the edge center.
   * @return The array of vectors pointing from the cell center to the edge center
   */
  ArrayOfArraysView< R1Tensor const > getCellCenterToEdgeCenters() const
  { return m_cellCenterToEdgeCenters.toViewConst(); }

private:

  ArrayOfArrays< R1Tensor > m_cellCenterToEdgeCenters;

};

template< typename PERMTYPE >
GEOSX_HOST_DEVICE
GEOSX_FORCE_INLINE
void SurfaceElementStencilWrapper::computeTransmissibility( localIndex iconn,
                                                            PERMTYPE permeability,
                                                            PERMTYPE dPerm_dPressure,
                                                            real64 (& transmissibility)[2],
                                                            real64 (& dTrans_dPressure )[2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * permeability[er0][esr0][ei0][0][0]; // this is a bit insane to access perm
  real64 const t1 = m_weights[iconn][1] * permeability[er1][esr1][ei1][0][0];

  real64 const harmonicWeight   = t0*t1 / (t0+t1);
  real64 const arithmeticWeight = (t0+t1)/2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;

  transmissibility[0] = value;
  transmissibility[1] = -value;

  real64 const dt0 = m_weights[iconn][0] * dPerm_dPressure[er0][esr0][ei0][0][0];
  real64 const dt1 = m_weights[iconn][1] * dPerm_dPressure[er1][esr1][ei1][0][0];

  real64 dHarmonic[2];
  dHarmonic[0] = dt0 * t1 / (t0 + t1 ) - dt0 / ( (t0 + t1) * (t0 + t1) );
  dHarmonic[1] = t0 * dt1 / (t0 + t1 ) - dt1 / ( (t0 + t1) * (t0 + t1) );

  real64 dArithmetic[2];
  dArithmetic[0] = dt0 / 2;
  dArithmetic[1] = dt1 / 2;

  dTrans_dPressure[0] = meanPermCoeff * dHarmonic[0] + (1 - meanPermCoeff) * dArithmetic[0];
  dTrans_dPressure[1] = -( meanPermCoeff * dHarmonic[1] + (1 - meanPermCoeff) * dArithmetic[1] );
}


template< typename PERMTYPE >
GEOSX_HOST_DEVICE
void SurfaceElementStencilWrapper::computeTransmissibility( localIndex iconn,
                                                            PERMTYPE permeability,
                                                            PERMTYPE dPerm_dPressure,
                                                            PERMTYPE dPerm_dAperture,
                                                            real64 (& transmissibility)[2],
                                                            real64 (& dTrans_dPressure )[2],
                                                            real64 (& dTrans_dAperture )[2] ) const
{
  localIndex const er0  =  m_elementRegionIndices[iconn][0];
  localIndex const esr0 =  m_elementSubRegionIndices[iconn][0];
  localIndex const ei0  =  m_elementIndices[iconn][0];

  localIndex const er1  =  m_elementRegionIndices[iconn][1];
  localIndex const esr1 =  m_elementSubRegionIndices[iconn][1];
  localIndex const ei1  =  m_elementIndices[iconn][1];

  real64 const t0 = m_weights[iconn][0] * permeability[er0][esr0][ei0][0][0]; // this is a bit insane to access perm
  real64 const t1 = m_weights[iconn][1] * permeability[er1][esr1][ei1][0][0];

  real64 const harmonicWeight   = t0*t1 / (t0+t1);
  real64 const arithmeticWeight = (t0+t1) / 2;

  real64 const meanPermCoeff = 1.0; //TODO make it a member

  real64 const value = meanPermCoeff * harmonicWeight + (1 - meanPermCoeff) * arithmeticWeight;

  transmissibility[0] = value;
  transmissibility[1] = -value;

  real64 const dt0_dp = m_weights[iconn][0] * dPerm_dPressure[er0][esr0][ei0][0][0];
  real64 const dt1_dp = m_weights[iconn][1] * dPerm_dPressure[er1][esr1][ei1][0][0];

  real64 dHarmonic_dpres[2];
  dHarmonic_dpres[0] = dt0_dp * t1 / (t0 + t1 ) - dt0_dp / ( (t0 + t1) * (t0 + t1) );
  dHarmonic_dpres[1] = t0 * dt1_dp / (t0 + t1 ) - dt1_dp / ( (t0 + t1) * (t0 + t1) );

  real64 dArithmetic_dpres[2];
  dArithmetic_dpres[0] = dt0_dp / 2;
  dArithmetic_dpres[1] = dt1_dp / 2;

  dTrans_dPressure[0] = meanPermCoeff * dHarmonic_dpres[0] + (1 - meanPermCoeff) * dArithmetic_dpres[0];
  dTrans_dPressure[1] = -( meanPermCoeff * dHarmonic_dpres[1] + (1 - meanPermCoeff) * dArithmetic_dpres[1] );

  real64 const dt0_da = m_weights[iconn][0] * dPerm_dAperture[er0][esr0][ei0][0][0];
  real64 const dt1_da = m_weights[iconn][1] * dPerm_dAperture[er1][esr1][ei1][0][0];

  real64 dHarmonic_dAper[2];
  dHarmonic_dAper[0] = dt0_da * t1 / (t0 + t1 ) - dt0_da / ( (t0 + t1) * (t0 + t1) );
  dHarmonic_dAper[1] = t0 * dt1_da / (t0 + t1 ) - dt1_da / ( (t0 + t1) * (t0 + t1) );

  real64 dArithmetic_dAper[2];
  dArithmetic_dAper[0] = dt0_da / 2;
  dArithmetic_dAper[1] = dt1_da / 2;

  dTrans_dAperture[0] = 0.0* ( meanPermCoeff * dHarmonic_dAper[0] + (1 - meanPermCoeff) * dArithmetic_dAper[0] );
  dTrans_dAperture[1] = 0.0* ( -( meanPermCoeff * dHarmonic_dAper[1] + (1 - meanPermCoeff) * dArithmetic_dAper[1] ) );
}



} /* namespace geosx */

#endif /* GEOSX_FINITEVOLUME_SURFACEELEMENTSTENCIL_HPP_ */
