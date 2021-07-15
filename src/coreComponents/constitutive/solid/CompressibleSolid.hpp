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
 * @file CompressibleSolid.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_
#define GEOSX_CONSTITUTIVE_SOLID_COMPRESSIBLESOLILD_HPP_

#include "constitutive/solid/CoupledSolid.hpp"
#include "constitutive/NullModel.hpp"

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
template< typename PORO_TYPE,
          typename PERM_TYPE >
class CompressibleSolidUpdates : public CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   */
  CompressibleSolidUpdates( NullModel const & solidModel,
                            PORO_TYPE const & porosityModel,
                            PERM_TYPE const & permModel ):
    CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >( solidModel, porosityModel, permModel )
  {}

  GEOSX_HOST_DEVICE
  void pressureUpdate( localIndex const k,
                       localIndex const q,
                       real64 const & pressure,
                       real64 const & deltaPressure ) const
  {
    m_porosityUpdate.updatePorosity( k, q, pressure, deltaPressure );
    real64 const porosity = m_porosityUpdate.getPorosity( k, q );
    m_permUpdate.updateFromPorosity( k, q, porosity );
  }

private:
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_solidUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_porosityUpdate;
  using CoupledSolidUpdates< NullModel, PORO_TYPE, PERM_TYPE >::m_permUpdate;

};

template< typename PORO_TYPE,
          typename PERM_TYPE >
class CompressibleSolid : public CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >
{
public:

  /**
   * @brief Constructor
   * @param name Object name
   * @param parent Object's parent group
   */
  CompressibleSolid( string const & name, dataRepository::Group * const parent );

  /// Destructor
  virtual ~CompressibleSolid() override;

  /**
   * @brief Catalog name
   * @return Static catalog string
   */
  static string catalogName() { return string( "CompressibleSolid" ) + PORO_TYPE::catalogName() + PERM_TYPE::catalogName(); }

  /**
   * @brief Get catalog name
   * @return Catalog name string
   */
  virtual string getCatalogName() const override { return catalogName(); }


  /**
   * @brief Create a instantiation of the CompressibleSolidUpdates class
   *        that refers to the data in this.
   * @return An instantiation of CompressibleSolidUpdates.
   */
  CompressibleSolidUpdates< PORO_TYPE, PERM_TYPE > createKernelUpdates() const
  {

    return CompressibleSolidUpdates< PORO_TYPE, PERM_TYPE >( getSolidModel(),
                                                             getPorosityModel(),
                                                             getPermModel() );
  }
private:
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getSolidModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPorosityModel;
  using CoupledSolid< NullModel, PORO_TYPE, PERM_TYPE >::getPermModel;

};

}
} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_SOLID_POROELASTIC_HPP_ */
