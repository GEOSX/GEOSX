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
 * @file PorousSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/solid/porosity/BiotPorosity.hpp"
#include "constitutive/solid/SolidBase.hpp"
#include "constitutive/permeability/ConstantPermeability.hpp"

namespace geosx
{
namespace constitutive
{

/**
 * @brief Provides kernel-callable constitutive update routines
 *
 *
 * @tparam SOLID_TYPE
 */
template< typename SOLID_TYPE >
class PorousSolidUpdates : public CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  using DiscretizationOps = typename SOLID_TYPE::KernelWrapper::DiscretizationOps;

  /**
   * @brief Constructor
   */
  PorousSolidUpdates( SOLID_TYPE const & solidModel,
                      BiotPorosity const & porosityModel,
                      ConstantPermeability const & permModel ):
    CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >( solidModel, porosityModel, permModel )
  {}

  GEOSX_HOST_DEVICE
  void smallStrainUpdate( localIndex const k,
                          localIndex const q,
                          real64 const & deltaPressure,
                          real64 const ( &strainIncrement )[6],
                          real64 ( & stress )[6],
                          real64 & dPorosity_dPressure,
                          real64 & dPorosity_dVolStrain,
                          real64 & dTotalStress_dPressure,
                          DiscretizationOps & stiffness ) const
  {
    m_solidUpdate.smallStrainUpdate( k, q, strainIncrement, stress, stiffness );

    updateBiotCoefficient( k );

    m_porosityUpdate.updatePorosity( k,
                                     q,
                                     deltaPressure,
                                     strainIncrement,
                                     dPorosity_dPressure,
                                     dPorosity_dVolStrain,
                                     dTotalStress_dPressure );

// TODO uncomment once we start using permeability model in flow.
//    m_permUpdate.updateFromPressureStrain( k,
//                                           q,
//                                           pressure,
//                                           volStrain );
  }

  GEOSX_HOST_DEVICE
  void updateBiotCoefficient( localIndex const k ) const
  {
    // This call is not general like this.
    real64 const bulkModulus = m_solidUpdate.getBulkModulus( k );

    m_porosityUpdate.updateBiotCoefficient( k, bulkModulus );
  }

private:

  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_solidUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_porosityUpdate;
  using CoupledSolidUpdates< SOLID_TYPE, BiotPorosity, ConstantPermeability >::m_permUpdate;

};

/**
 * @brief PorousSolidBase class used for dispatch of all Porous solids.
 */
class PorousSolidBase : public SolidBase
{};

/**
 * @brief Class to represent a coupled solid model
 */
template< typename SOLID_TYPE >
class PorousSolid : public CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >
{
public:

  /// Alias for ElasticIsotropicUpdates
  using KernelWrapper = PorousSolidUpdates< SOLID_TYPE >;

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  PorousSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~PorousSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "Porous" ) + SOLID_TYPE::m_catalogNameString; }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }

  /**
   * @brief Create a instantiation of the PorousSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of PorousSolidUpdates.
   */
  KernelWrapper createKernelUpdates() const
  {
    return KernelWrapper( getSolidModel(),
                          getPorosityModel(),
                          getPermModel() );
  }

  /**
   * @brief Const/non-mutable accessor for density
   * @return Accessor
   */
  arrayView2d< real64 const > const getDensity() const
  {
    return getSolidModel().getDensity();
  }

private:
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getSolidModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPorosityModel;
  using CoupledSolid< SOLID_TYPE, BiotPorosity, ConstantPermeability >::getPermModel;
};



}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROUSSOLID_HPP_ */
