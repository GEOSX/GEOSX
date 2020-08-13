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

#ifndef SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_
#define SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_

#include "ConstitutiveBase.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @class NullModel A null constitutive relation
 */
class NullModel : public constitutive::ConstitutiveBase
{
public:
  /// @copydoc geosx::dataRepository::Group::Group
  NullModel( string const & name, Group * const parent );

  /// Destrutor
  virtual ~NullModel();

  virtual void
  DeliverClone( string const & GEOSX_UNUSED_PARAM( name ),
                Group * const GEOSX_UNUSED_PARAM( parent ),
                std::unique_ptr< ConstitutiveBase > & GEOSX_UNUSED_PARAM(
                  clone ) ) const override final;

  /// string name to use for this class in the catalog
  static constexpr auto m_catalogNameString = "NullModel";

  /**
   * @return A string that is used to register/lookup this class in the registry
   */
  static std::string
  CatalogName()
  {
    return m_catalogNameString;
  }

  virtual string
  GetCatalogName() override
  {
    return CatalogName();
  }

  /**
   * Empty struct to serve as a KernelWrapper for the constitutive model.
   */
  struct KernelWrapper
  {};

  /**
   * @brief Create a kernel wrapper for this constitutive relation.
   * @param includeState Whether or not to include the state in the wrapper.
   * @return 0
   */
  KernelWrapper
  createKernelUpdates( bool const includeState = false )
  {
    GEOSX_UNUSED_VAR( includeState );
    return KernelWrapper();
  }
};

}  // namespace constitutive
} /* namespace geosx */

#endif /* SRC_CORECOMPONENTS_CONSTITUTIVE_NULLMODEL_HPP_ */
