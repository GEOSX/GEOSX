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
 * @file ConstitutivePassThru.hpp
 */

#ifndef GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_
#define GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_

#include "NullModel.hpp"
#include "solid/DamageVolDev.hpp"
#include "solid/DamageSpectral.hpp"
#include "solid/DruckerPrager.hpp"
#include "solid/DruckerPragerExtended.hpp"
#include "solid/ElasticIsotropic.hpp"
#include "solid/ElasticTransverseIsotropic.hpp"
#include "solid/ElasticOrthotropic.hpp"
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

  // NOTE: The switch order here can be fragile if a model derives from another
  //       model, as the dynamic_cast will also cast to a base version.
  //       Models should be ordered such that children come before parents.
  //       For example, DruckerPrager before ElasticIsotropic, DamageVolDev before
  //       Damage, etc.

  template< typename LAMBDA >
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr1 = dynamic_cast< DamageSpectral< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr1 );
    }
    else if( auto * const ptr2 = dynamic_cast< DamageVolDev< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr2 );
    }
    else if( auto * const ptr3 = dynamic_cast< Damage< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr3 );
    }
    else if( auto * const ptr4 = dynamic_cast< DruckerPragerExtended * >( &constitutiveRelation ) )
    {
      lambda( *ptr4 );
    }
    else if( auto * const ptr5 = dynamic_cast< DruckerPrager * >( &constitutiveRelation ) )
    {
      lambda( *ptr5 );
    }
    else if( auto * const ptr6 = dynamic_cast< ElasticIsotropic * >( &constitutiveRelation ) )
    {
      lambda( *ptr6 );
    }
    else if( auto * const ptr7 = dynamic_cast< ElasticTransverseIsotropic * >( &constitutiveRelation ) )
    {
      lambda( *ptr7 );
    }
    else if( auto * const ptr8 = dynamic_cast< ElasticOrthotropic * >( &constitutiveRelation ) )
    {
      lambda( *ptr8 );
    }
    else
    {
      GEOSX_ERROR( "ConstitutivePassThru< SolidBase >::execute failed. The constitutive relation is named "
                   << constitutiveRelation.getName() << " with type "
                   << LvArray::system::demangleType( constitutiveRelation ) );
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
  static
  void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr = dynamic_cast< NullModel * >( &constitutiveRelation ) )
    {
      lambda( *ptr );
    }
    else
    {
      GEOSX_ERROR( "ConstitutivePassThru< NullModel >::execute failed. The constitutive relation is named "
                   << constitutiveRelation.getName() << " with type "
                   << LvArray::system::demangleType( constitutiveRelation ) );
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
  static void execute( ConstitutiveBase & constitutiveRelation, LAMBDA && lambda )
  {
    if( auto * const ptr1 = dynamic_cast< PoroElastic< DruckerPragerExtended > * >( &constitutiveRelation ) )
    {
      lambda( *ptr1 );
    }
    else if( auto * const ptr2 = dynamic_cast< PoroElastic< DruckerPrager > * >( &constitutiveRelation ) )
    {
      lambda( *ptr2 );
    }
    else if( auto * const ptr3 = dynamic_cast< PoroElastic< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr3 );
    }
    else if( auto * const ptr4 = dynamic_cast< PoroElastic< ElasticTransverseIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr4 );
    }
    else if( auto * const ptr5 = dynamic_cast< PoroElastic< ElasticOrthotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr5 );
    }
    else
    {
      GEOSX_ERROR( "ConstitutivePassThru< PoroElasticBase >::execute failed. The constitutive relation is named "
                   << constitutiveRelation.getName() << " with type "
                   << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};

/**
 * Specialization for the Damage models.
 */
template<>
struct ConstitutivePassThru< DamageBase >
{
  template< typename LAMBDA >
  static void execute( ConstitutiveBase & constitutiveRelation,
                       LAMBDA && lambda )
  {
    if( auto * const ptr1 = dynamic_cast< DamageSpectral< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr1 );
    }
    else if( auto * const ptr2 = dynamic_cast< DamageVolDev< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr2 );
    }
    else if( auto * const ptr3 = dynamic_cast< Damage< ElasticIsotropic > * >( &constitutiveRelation ) )
    {
      lambda( *ptr3 );
    }
    else
    {
      GEOSX_ERROR( "ConstitutivePassThru< DamageBase >::execute failed. The constitutive relation is named "
                   << constitutiveRelation.getName() << " with type "
                   << LvArray::system::demangleType( constitutiveRelation ) );
    }
  }
};

} /* namespace constitutive */

} /* namespace geosx */

#endif /* GEOSX_CONSTITUTIVE_CONSTITUTIVEPASSTHRU_HPP_ */
