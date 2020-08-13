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
 * @file ConstitutivePassThru.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_
#define GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_

#include "NullModel.hpp"
#include "solid/LinearElasticIsotropic.hpp"
#include "solid/LinearElasticAnisotropic.hpp"
#include "solid/LinearElasticTransverseIsotropic.hpp"
#include "solid/PoroElastic.hpp"

namespace geosx
{
namespace constitutive
{
/**
 * @struct ConstitutivePassThru
 * @brief Struct to facilitate launching of lambda functions with a compile
 *   time knowledge of what constitutive model is used.
 *
 * This struct works by implementing an if-else or switch-case block for a
 * specific constitutive base type, and executing the lambda passing it a
 * casted pointer to the constitutive relation.
 */
template< typename BASETYPE >
struct ConstitutivePassThru;

/**
 * Specialization for models that derive from SolidBase.
 */
template<>
struct ConstitutivePassThru< SolidBase >
{
  template< typename LAMBDA >
  static void
  Execute( ConstitutiveBase * const constitutiveRelation,
           LAMBDA && lambda )
  {
    GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                    "ConstitutiveBase* == nullptr" );

    if( dynamic_cast< LinearElasticIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticIsotropic * >( constitutiveRelation ) );
    }
    else if( dynamic_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticTransverseIsotropic * >( constitutiveRelation ) );
    }
    else if( dynamic_cast< LinearElasticAnisotropic * >( constitutiveRelation ) )
    {
      lambda( static_cast< LinearElasticAnisotropic * >( constitutiveRelation ) );
    }
    else
    {
      string name;
      if( constitutiveRelation != nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<SolidBase>::Execute( "
                   << constitutiveRelation << " ) failed. ( "
                   << constitutiveRelation << " ) is named " << name );
    }
  }
};

/**
 * Specialization for the NullModel.
 */
template<>
struct ConstitutivePassThru< NullModel >
{
  template< typename LAMBDA >
  static void
  Execute( ConstitutiveBase * const constitutiveRelation,
           LAMBDA && lambda )
  {
    GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                    "ConstitutiveBase* == nullptr" );

    if( dynamic_cast< NullModel * >( constitutiveRelation ) )
    {
      lambda( static_cast< NullModel * >( constitutiveRelation ) );
    }
    else
    {
      string name;
      if( constitutiveRelation != nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<NullModel>::Execute( "
                   << constitutiveRelation << " ) failed. ( "
                   << constitutiveRelation << " ) is named " << name );
    }
  }
};

/**
 * Specialization for the PoroElastic models.
 */
template<>
struct ConstitutivePassThru< PoroElasticBase >
{
  template< typename LAMBDA >
  static void
  Execute( ConstitutiveBase * const constitutiveRelation,
           LAMBDA && lambda )
  {
    GEOSX_ERROR_IF( constitutiveRelation == nullptr,
                    "ConstitutiveBase* == nullptr" );

    if( dynamic_cast< PoroElastic< LinearElasticIsotropic > * >( constitutiveRelation ) )
    {
      lambda(
        static_cast< PoroElastic< LinearElasticIsotropic > * >( constitutiveRelation ) );
    }
    else if( dynamic_cast< PoroElastic< LinearElasticTransverseIsotropic > * >(
               constitutiveRelation ) )
    {
      lambda( static_cast< PoroElastic< LinearElasticTransverseIsotropic > * >(
        constitutiveRelation ) );
    }
    else if( dynamic_cast< PoroElastic< LinearElasticAnisotropic > * >(
               constitutiveRelation ) )
    {
      lambda( static_cast< PoroElastic< LinearElasticAnisotropic > * >(
        constitutiveRelation ) );
    }
    else
    {
      string name;
      if( constitutiveRelation != nullptr )
      {
        name = constitutiveRelation->getName();
      }
      GEOSX_ERROR( "ConstitutivePassThru<SolidBase>::Execute( "
                   << constitutiveRelation << " ) failed. ( "
                   << constitutiveRelation << " ) is named " << name );
    }
  }
};

}  // namespace constitutive
}  // namespace geosx

#endif /* GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_ */
