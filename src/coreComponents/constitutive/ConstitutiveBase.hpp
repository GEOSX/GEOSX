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
 * @file ConstitutiveBase.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_
#define GEOSX_CONSTITUTIVE_CONSTITUTIVEBASE_HPP_

#include "dataRepository/ObjectCatalog.hpp"
#include "common/DataTypes.hpp"
#include "dataRepository/Group.hpp"

namespace geosx
{
namespace constitutive
{
class ConstitutiveBase : public dataRepository::Group
{
public:
  /**
   * Single point of reference for generated constitutive field naming convention
   *
   * @param prefix name prefix (e.g. constitutive model name)
   * @param name actual field name
   * @return prefixed field name that is used to access data
   */
  inline static string
  makeFieldName( string const & prefix, string const & name )
  {
    return prefix + "_" + name;
  }

  ConstitutiveBase( string const & name, Group * const parent );

  virtual ~ConstitutiveBase() override;

  /**
   * @brief create a clone of this constitutive model
   * @param[in]  name   The name of the clone in the repository
   * @param[in]  parent A pointer to the group that contains the instance of the new clone
   * @param[out] clone  A reference to a unique_ptr  that will hold the clone.
   */
  virtual void
  DeliverClone( string const & name,
                Group * const parent,
                std::unique_ptr< ConstitutiveBase > & clone ) const = 0;

  virtual void
  StateUpdatePointPressure( real64 const & GEOSX_UNUSED_PARAM( pres ),
                            localIndex const GEOSX_UNUSED_PARAM( k ),
                            localIndex const GEOSX_UNUSED_PARAM( q ) )
  {}

  virtual void
  StateUpdateBatchPressure( arrayView1d< real64 const > const & pres,
                            arrayView1d< real64 const > const & dPres )
  {
    GEOSX_UNUSED_VAR( pres )
    GEOSX_UNUSED_VAR( dPres )
  }

  /**
   * @brief function to resize the fields in this constitutive model
   * @param[in] newSize the new size of the fields
   */
  virtual void
  resize( localIndex newSize ) override;

  /**
   * @name Static Factory Catalog members and functions
   */
  ///@{

  /// @typedef An alias for the ConstitutiveBase catalog
  using CatalogInterface =
    dataRepository::CatalogInterface< ConstitutiveBase, std::string const &, Group * const >;

  /**
   * @brief Singleton accessor for catalog
   * @return
   */
  static typename CatalogInterface::CatalogType &
  GetCatalog();

  /**
   * @brief function to return the catalog name of the derived class
   * @return a string that contains the catalog name of the derived class
   */
  virtual string
  GetCatalogName() = 0;

  ///@}

  /**
   * @brief Allocate constitutive data and make views to data on parent objects
   * @param[in] parent pointer to the group that holds the constitutive relation
   * @param[in] numConstitutivePointsPerParentIndex number of quadrature points
   *
   * This function does 2 things:
   *   1) Allocate data according to the size of parent and numConstitutivePointsPerParentIndex
   *   2) Create wrappers to the constitutive data in the parent for easier access
   */
  virtual void
  AllocateConstitutiveData(
    dataRepository::Group * const parent,
    localIndex const numConstitutivePointsPerParentIndex );

  struct viewKeyStruct
  {
    static constexpr auto poreVolumeMultiplierString = "poreVolumeMultiplier";
    static constexpr auto dPVMult_dPresString = "dPVMult_dDensity";
  };

  localIndex
  numQuadraturePoints() const
  {
    return m_numQuadraturePoints;
  }

protected:
private:
  localIndex m_numQuadraturePoints;
  Group * m_constitutiveDataGroup = nullptr;

  ConstitutiveBase( ConstitutiveBase const & ) = delete;
  ConstitutiveBase( ConstitutiveBase && ) = delete;
  ConstitutiveBase const &
  operator=( ConstitutiveBase const & ) = delete;
  ConstitutiveBase const &
  operator=( ConstitutiveBase && ) = delete;
};

}  // namespace constitutive
}  // namespace geosx
#endif
